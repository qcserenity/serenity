/**
 * @file WavefunctionEmbeddingTask.cpp
 *
 * @author Moritz Bensberg
 * @date Jul 8, 2020
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
/* Include Serenity Internal Headers */
#include "tasks/WavefunctionEmbeddingTask.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/OrbitalPair.h"
#include "data/OrbitalTriple.h"
#include "energies/EnergyContributions.h" //Get/add specific energy contributions.
#include "geometry/Geometry.h"
#include "io/FormattedOutputStream.h"
#include "misc/SystemSplittingTools.h"
#include "postHF/CC/DLPNO_CCSD.h" //Local coupled cluster.
#include "postHF/CC/DLPNO_CCSD_T0.h"
#include "potentials/bundles/PotentialBundle.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/SystemAdditionTask.h"

namespace Serenity {

WavefunctionEmbeddingTask::WavefunctionEmbeddingTask(std::shared_ptr<SystemController> supersystem,
                                                     std::vector<std::shared_ptr<SystemController>> systems)
  : _supersystem(supersystem), _systems(systems) {
  for (auto env : systems) {
    LocalCorrelationSettings newSettings;
    newSettings.enforceHFFockian = true;
    settings.lcSettings.push_back(newSettings);
  }

  _fragmentEnergies = Eigen::VectorXd::Zero(_systems.size());
  _interactionEnergies = Eigen::VectorXd::Zero(_systems.size());
}

WavefunctionEmbeddingTask::~WavefunctionEmbeddingTask() = default;

const Eigen::VectorXd& WavefunctionEmbeddingTask::getFragmentEnergies() {
  return _fragmentEnergies;
}
const Eigen::VectorXd& WavefunctionEmbeddingTask::getFragmentWiseInteractionEnergy() {
  return _interactionEnergies;
}

Eigen::SparseMatrix<int> WavefunctionEmbeddingTask::setUpSubsystems() {
  LocalizationTask locTask(_supersystem);
  locTask.settings = settings.loc;
  locTask.run();
  SystemSplittingTask<RESTRICTED> splittingTask(_supersystem, _systems);
  splittingTask.settings = settings.split;
  splittingTask.run();
  Eigen::VectorXi finalAssignment = splittingTask.getFinalAssignment();
  unsigned int nOccSup = finalAssignment.size();
  Eigen::MatrixXi orbitalIndexMap = Eigen::MatrixXi::Zero(nOccSup, _systems.size());
  for (unsigned int iOcc = 0; iOcc < nOccSup; ++iOcc) {
    orbitalIndexMap(iOcc, finalAssignment(iOcc)) = 1;
  } // for iOcc
  return orbitalIndexMap.sparseView();
}

Eigen::SparseMatrix<int> WavefunctionEmbeddingTask::setUpSupersystem() {
  for (const auto& sys : _systems) {
    if (not sys->hasElectronicStructure<RESTRICTED>()) {
      throw SerenityError((std::string) "ERROR: The fragment " + sys->getSystemName() + " has no electronic structure!\n" +
                          "       The WavefunctionEmbeddingTask does not know the logic on how to construct it\n" +
                          "       if 'fromFragments true' is given. Please set it to false or supply all fragments\n" +
                          "       with (reasonalble) orbitals!");
    }
  }
  SystemAdditionTask<RESTRICTED> additionTask(_supersystem, _systems);
  additionTask.settings.addOccupiedOrbitals = true;
  additionTask.run();
  unsigned int nOccSup = _supersystem->getNOccupiedOrbitals<RESTRICTED>();
  Eigen::MatrixXi orbitalIndexMap = Eigen::MatrixXi::Zero(nOccSup, _systems.size());
  unsigned int row = 0;
  for (unsigned int iSys = 0; iSys < _systems.size(); ++iSys) {
    unsigned int nOccSub = _systems[iSys]->getNOccupiedOrbitals<RESTRICTED>();
    orbitalIndexMap.block(row, iSys, nOccSub, 1) = Eigen::VectorXi::Constant(nOccSub, 1);
    row += nOccSub;
  } // for iSys
  return orbitalIndexMap.sparseView();
}

std::vector<std::shared_ptr<OrbitalPair>> WavefunctionEmbeddingTask::buildOrbialPairs(unsigned int iSys, unsigned int jSys) {
  const Eigen::VectorXi& coreOrbitals = _supersystem->getActiveOrbitalController<RESTRICTED>()->getCoreOrbitals();
  const bool usesFrozenCore = settings.lcSettings[iSys].useFrozenCore || settings.lcSettings[jSys].useFrozenCore;
  unsigned int thresholdSystem = (settings.accurateInteraction) ? iSys : jSys;
  const auto& lcSettings = settings.lcSettings[thresholdSystem];
  const double pnoThreshold = lcSettings.pnoThreshold;
  const double coreScaling = lcSettings.pnoCoreScaling;
  const double pnoCoreThreshold = pnoThreshold * coreScaling;
  bool eligibleForTriples = lcSettings.method == Options::PNO_METHOD::DLPNO_CCSD_T0;
  std::vector<std::shared_ptr<OrbitalPair>> pairs;
  if (settings.lcSettings[iSys].method == Options::PNO_METHOD::NONE ||
      settings.lcSettings[jSys].method == Options::PNO_METHOD::NONE)
    return pairs;
  for (Eigen::SparseMatrix<int>::InnerIterator itI(_orbitalIndexMap, iSys); itI; ++itI) {
    unsigned int iOcc = itI.row();
    bool iIsCore = coreOrbitals(iOcc);
    if (usesFrozenCore && iIsCore)
      continue;
    for (Eigen::SparseMatrix<int>::InnerIterator itJ(_orbitalIndexMap, jSys); itJ; ++itJ) {
      unsigned int jOcc = itJ.row();
      bool jIsCore = coreOrbitals(jOcc);
      if ((usesFrozenCore && jIsCore) || (jOcc > iOcc && iSys == jSys))
        continue;
      unsigned int i = iOcc;
      unsigned int j = jOcc;
      if (iOcc < jOcc) {
        i = jOcc;
        j = iOcc;
      }
      auto newPair = std::make_shared<OrbitalPair>(i, j, (iIsCore || jIsCore) ? pnoCoreThreshold : pnoThreshold,
                                                   lcSettings.ccsdPairThreshold, lcSettings.collinearDipoleScaling,
                                                   lcSettings.extendedDomainScaling);
      newPair->eligibleForTriples = eligibleForTriples;
      pairs.push_back(newPair);
    } // for itJ/jOcc
  }   // for itI/iOcc
  return pairs;
}

Eigen::VectorXi WavefunctionEmbeddingTask::buildOrbitalAssignments() {
  Eigen::VectorXi assignments = Eigen::VectorXi::Zero(_orbitalIndexMap.rows());
  for (unsigned int iSys = 0; iSys < _systems.size(); ++iSys) {
    for (Eigen::SparseMatrix<int>::InnerIterator itI(_orbitalIndexMap, iSys); itI; ++itI) {
      assignments(itI.row()) = iSys;
    } // for iOcc
  }   // for iSys
  return assignments;
}

double WavefunctionEmbeddingTask::calculateTriplesCorrection(std::shared_ptr<LocalCorrelationController> localCorrelationController) {
  return DLPNO_CCSD_T0::calculateEnergyCorrection(localCorrelationController);
}

void WavefunctionEmbeddingTask::sortOrbitalTriples(std::shared_ptr<LocalCorrelationController> localCorrelationController) {
  const unsigned int nSystems = _systems.size();
  std::vector<std::shared_ptr<OrbitalTriple>> emptyTriplesVector;
  _crossTriples = std::make_shared<std::vector<std::vector<std::shared_ptr<OrbitalTriple>>>>(nSystems, emptyTriplesVector);
  _diagonalTriples = std::make_shared<std::vector<std::vector<std::shared_ptr<OrbitalTriple>>>>(nSystems, emptyTriplesVector);
  std::vector<std::shared_ptr<OrbitalTriple>> allTriples = localCorrelationController->getOrbitalTriples();
  auto& crossTriples = *_crossTriples;
  auto& diagonalTriples = *_diagonalTriples;
  for (auto triple : allTriples) {
    unsigned int iSys = _orbitalAssignments(triple->getI());
    unsigned int jSys = _orbitalAssignments(triple->getJ());
    unsigned int kSys = _orbitalAssignments(triple->getK());
    if (iSys == jSys && iSys == kSys) {
      // It a diagonal triple.
      diagonalTriples[iSys].push_back(triple);
    }
    else {
      // Non diagonal! Save it for i. Save it only for j and k if not already done.
      crossTriples[iSys].push_back(triple);
      if (iSys != jSys)
        crossTriples[jSys].push_back(triple);
      if (iSys != kSys && jSys != kSys)
        crossTriples[kSys].push_back(triple);
    }
  } // for triple
}

void WavefunctionEmbeddingTask::performFullDecomposition(unsigned int iSys) {
  bool oldIOOption = iOOptions.printFinalOrbitalEnergies;
  iOOptions.printFinalOrbitalEnergies = false;
  auto activeSystem = _systems[iSys];
  auto environmentSystems = _systems;
  environmentSystems.erase(std::remove(environmentSystems.begin(), environmentSystems.end(), activeSystem),
                           environmentSystems.end());
  FDETask<RESTRICTED> fde(activeSystem, environmentSystems);
  fde.settings.calculateEnvironmentEnergy = true;
  fde.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::HF;
  fde.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
  fde.settings.skipSCF = true;
  fde.run();
  iOOptions.printFinalOrbitalEnergies = oldIOOption;
}

void WavefunctionEmbeddingTask::buildIntegralThresholdVectors() {
  unsigned int nOccSuper = _orbitalAssignments.size();
  _orbitalWiseMullikenThresholds = Eigen::VectorXd::Zero(nOccSuper);
  _orbitalToShellThresholds = Eigen::VectorXd::Zero(nOccSuper);
  _orbitalWiseDOIPAOThresholds = Eigen::VectorXd::Zero(nOccSuper);
  for (unsigned int iOcc = 0; iOcc < nOccSuper; ++iOcc) {
    _orbitalWiseMullikenThresholds(iOcc) = settings.lcSettings[_orbitalAssignments(iOcc)].mullikenThreshold;
    _orbitalToShellThresholds(iOcc) = settings.lcSettings[_orbitalAssignments(iOcc)].orbitalToShellThreshold;
    _orbitalWiseDOIPAOThresholds(iOcc) = settings.lcSettings[_orbitalAssignments(iOcc)].doiPAOThreshold;
  } // for iOcc
}

void WavefunctionEmbeddingTask::run() {
  for (auto& lcSettings : settings.lcSettings)
    lcSettings.resolvePNOSettings();
  // Build supersystem
  if (settings.fromFragments) {
    _orbitalIndexMap = setUpSupersystem();
  }
  else {
    _orbitalIndexMap = setUpSubsystems();
  }
  bool usesTriples = false;
  for (unsigned int i = 0; i < settings.lcSettings.size(); ++i) {
    auto& lcSettings = settings.lcSettings[i];
    if (lcSettings.method == Options::PNO_METHOD::DLPNO_CCSD_T0 && _systems[i]->getNOccupiedOrbitals<RESTRICTED>() != 0)
      usesTriples = true;
  }
  _orbitalAssignments = buildOrbitalAssignments();
  buildIntegralThresholdVectors();
  // Build all pair combinations.
  std::vector<std::shared_ptr<OrbitalPair>> emptyPairVector;
  _pairMatrix = std::make_shared<Matrix<std::vector<std::shared_ptr<OrbitalPair>>>>(_systems.size(), _systems.size(),
                                                                                    emptyPairVector);
  auto& pairMatrix = *_pairMatrix;
  for (unsigned int iSys = 0; iSys < _systems.size(); ++iSys) {
    for (unsigned int jSys = iSys; jSys < _systems.size(); ++jSys) {
      auto pairSet = buildOrbialPairs(iSys, jSys);
      pairMatrix(iSys, jSys) = pairSet;
      pairMatrix(jSys, iSys) = pairSet;
      _allPairs.insert(_allPairs.end(), pairSet.begin(), pairSet.end());
    } // for jSys
  }   // for iSys
  std::shared_ptr<FockMatrix<RESTRICTED>> fockMatrix = nullptr;
  if (settings.fromFragments || _supersystem->getSettings().method != Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    auto eComp = _supersystem->getElectronicStructure<RESTRICTED>()->getEnergyComponentController();
    fockMatrix = std::make_shared<FockMatrix<RESTRICTED>>(
        _supersystem->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>()->getFockMatrix(
            _supersystem->getElectronicStructure<RESTRICTED>()->getDensityMatrix(), eComp));
  }
  double totalHFEnergy = _supersystem->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY);
  // Run CC calculation.
  std::vector<std::shared_ptr<SystemController>> dummy = {};
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(
      _supersystem, settings.lcSettings[0], dummy, fockMatrix, _allPairs, _orbitalWiseMullikenThresholds,
      _orbitalToShellThresholds, _orbitalWiseDOIPAOThresholds);
  DLPNO_CCSD dlpnoCCSD(localCorrelationController, settings.normThreshold, settings.maxCycles);
  Eigen::VectorXd energies = dlpnoCCSD.calculateElectronicEnergyCorrections();
  double totalCCSD = energies.sum();
  double triplesCorrection = 0.0;
  if (usesTriples) {
    sortOrbitalTriples(localCorrelationController);
    triplesCorrection = calculateTriplesCorrection(localCorrelationController);
  }

  for (unsigned int iSys = 0; iSys < _systems.size(); ++iSys) {
    auto sys = _systems[iSys];
    double totalACorrelationEnergy = 0.0;
    double totalCorrelationInteractionEnergy = 0.0;
    if (settings.fullDecomposition)
      performFullDecomposition(iSys);
    double aaPairEnergies = 0.0;
    for (auto& pair : pairMatrix(iSys, iSys)) {
      aaPairEnergies += pair->getCCSDPairEnergy();
    }
    double abPairEnergies = 0.0;
    unsigned int nABPairs = 0;
    for (unsigned int jSys = 0; jSys < _systems.size(); ++jSys) {
      if (jSys == iSys)
        continue;
      for (auto& pair : pairMatrix(iSys, jSys))
        abPairEnergies += pair->getCCSDPairEnergy();
      nABPairs += pairMatrix(iSys, jSys).size();
    }

    OutputControl::mOut << "Energy decomposition for system: " << sys->getSystemName() << std::endl;
    printf("%-80s %18.10f \n", "AA pair energy correction: ", aaPairEnergies);
    printf("%-80s %18.10f \n", "AB pair energy correction: ", abPairEnergies);
    totalACorrelationEnergy += aaPairEnergies;
    totalCorrelationInteractionEnergy += abPairEnergies;
    if (usesTriples) {
      auto& diagonalTriples = *_diagonalTriples;
      auto& crossTriples = *_crossTriples;
      double aaaTripleEnergy = 0.0;
      for (const auto& triple : diagonalTriples[iSys])
        aaaTripleEnergy += triple->getTripleEnergy();
      double abcTripleEnergy = 0.0;
      for (const auto& triple : crossTriples[iSys])
        abcTripleEnergy += triple->getTripleEnergy();
      printf("%-80s %18.10f \n", "AAA triples energy correction: ", aaaTripleEnergy);
      printf("%-80s %18.10f \n", "ABC triples energy correction: ", abcTripleEnergy);
      totalACorrelationEnergy += aaaTripleEnergy;
      totalCorrelationInteractionEnergy += abcTripleEnergy;
    }
    OutputControl::mOut
        << "----------------------------------------------------------------------------------------------------"
        << std::endl;

    _fragmentEnergies(iSys) = totalACorrelationEnergy;
    _interactionEnergies(iSys) = totalCorrelationInteractionEnergy;
  }
  _totalEnergy = std::make_unique<double>(triplesCorrection + totalCCSD + totalHFEnergy);
  OutputControl::mOut << "Total supersystem energies: " << _supersystem->getSystemName() << std::endl;
  OutputControl::mOut
      << "----------------------------------------------------------------------------------------------------" << std::endl;
  printf("%-80s %18.10f \n", "Total Hartree Fock Energy:", totalHFEnergy);
  printf("%-80s %18.10f \n", "Total LCCSD Energy:", totalCCSD);
  if (usesTriples)
    printf("%-80s %18.10f \n", "Total Semi-Can. Triples correction:", triplesCorrection);
  printf("%-80s %18.10f \n", "Total Correlation Energy:", triplesCorrection + totalCCSD);
  printf("%-80s %18.10f \n", "Total Energy (HF + corr):", *_totalEnergy);
  OutputControl::mOut
      << "----------------------------------------------------------------------------------------------------" << std::endl;
}

double WavefunctionEmbeddingTask::getTotalEnergy() {
  if (!_totalEnergy)
    this->run();
  return *_totalEnergy;
}

} /* namespace Serenity */
