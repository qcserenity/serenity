/**
 * @file DOSCCTask.cpp
 *
 * @author Moritz Bensberg
 * @date Jul 8, 2021
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
/* Include Class Header*/
#include "tasks/DOSCCTask.h"
/* Include Serenity Internal Headers */
#include "analysis/directOrbitalSelection/DirectOrbitalSelection.h" //Orbital group definition.
#include "data/ElectronicStructure.h"                               //Get energy component controller.
#include "data/OrbitalController.h"                                 //Core orbitals
#include "data/OrbitalPair.h"                                       //Orbital pair construction.
#include "energies/EnergyContributions.h"                           //Get/add specific energy contributions.
#include "geometry/Geometry.h"                                      //Dummy constructor for the geometry.
#include "io/FormattedOutput.h"                                     //Print nice section titles.
#include "io/FormattedOutputStream.h"                               //Filtered output.
#include "misc/SerenityError.h"                                     //Error messages.
#include "postHF/CC/DLPNO_CCSD.h"                                   //DLPNO-CCSD for pair selected CC.
#include "postHF/CC/DLPNO_CCSD_T0.h"                                //Triples correction for pair selected CC.
#include "potentials/bundles/PotentialBundle.h"                     //Get Fock matrix.
#include "settings/Settings.h"                                      //New system settings.
#include "system/SystemController.h"                                //Construct dummy systems.
#include "tasks/LocalizationTask.h"                                 //Run localization.
/* Include Std and External Headers */
#include <iomanip> //setw(...)

namespace Serenity {

DOSCCTask::DOSCCTask(std::vector<std::shared_ptr<SystemController>> supersystems) : _supersystems(supersystems) {
  settings.wfemb.lcSettings = std::vector<LocalCorrelationSettings>(_maxFragments);
  settings.lcPairSelection = std::vector<LocalCorrelationSettings>(_maxFragments);
}

DOSCCTask::~DOSCCTask() = default;

void DOSCCTask::run() {
  // Resolve DOS settings.
  settings.resolveDOSSettings();
  _nFragments = settings.gdos.similarityLocThreshold.size() + 1;
  // Resize the settings vector the number of desired fragments.
  settings.wfemb.lcSettings.resize(_nFragments);
  // Build subsystems
  setUpSubsystems();
  // Run the orbital localization for each supersystem.
  localizeOrbitals();
  // Run generalized DOS.
  const auto& orbitalGroups = runGDOS();
  // Run multi-level CC
  _totalEnergies = std::make_unique<Eigen::VectorXd>(runCC(_pairEnergyMatrices));
  for (auto& group : orbitalGroups) {
    group->updateOrbitalIndices(_sortingMatrices);
  } // for group
  // Print the final result / energy differences with respect to the 0-th energy.
  _relativeEnergies = std::make_unique<Eigen::VectorXd>(Eigen::VectorXd::Zero(_supersystems.size()));
  printResults(*_totalEnergies, *_relativeEnergies, orbitalGroups);
  if (settings.orbitalPairAnalysis) {
    this->resolvePNOSettingsForPairSelection();
    auto orbitalPairs = selectOrbitalPairsFromPairEnergies(orbitalGroups, _pairEnergyMatrices);
    _totalEnergiesFromPairSelected = std::make_unique<Eigen::VectorXd>(runPairSelectedCC(orbitalPairs));
    _relativeEnergiesFromPairSelected = std::make_unique<Eigen::VectorXd>(Eigen::VectorXd::Zero(_supersystems.size()));
    printResults(*_totalEnergiesFromPairSelected, *_relativeEnergiesFromPairSelected, {});
  }
}

Eigen::VectorXd DOSCCTask::runCC(std::vector<Eigen::MatrixXd>& pairEnergyMatrices) {
  Eigen::VectorXd energies(_supersystems.size());
  for (unsigned int i = 0; i < _supersystems.size(); ++i) {
    auto supersystem = _supersystems[i];
    auto subsystemSet = _subSystems[i];
    WavefunctionEmbeddingTask wfemb(supersystem, subsystemSet);
    wfemb.settings = settings.wfemb;
    wfemb.settings.fromFragments = true;
    wfemb.run();
    energies(i) = wfemb.getTotalEnergy();
    pairEnergyMatrices.push_back(wfemb.getLocalCorrelationController()->getPairEnergyMatrix());
  } // for i
  return energies;
}

const SpinPolarizedData<RESTRICTED, std::vector<std::shared_ptr<DOSOrbitalGroup>>> DOSCCTask::runGDOS() {
  std::vector<std::shared_ptr<SystemController>> allSubsystems;
  for (const auto& subsystemSet : _subSystems) {
    allSubsystems.insert(allSubsystems.end(), subsystemSet.begin(), subsystemSet.end());
  }
  GeneralizedDOSTask<RESTRICTED> gdos(_supersystems, allSubsystems);
  gdos.settings = settings.gdos;
  gdos.run();
  const auto& tmp = gdos.getSuperToSubsystemOccSortingMatrices();
  for (const auto& mat : tmp)
    _sortingMatrices.push_back(mat);
  return gdos.getOrbitalGroups();
}

void DOSCCTask::localizeOrbitals() {
  // Run localization for the first point
  LocalizationTask locOne(_supersystems[0]);
  locOne.settings = settings.wfemb.loc;
  locOne.run();
  // Run localization for all other points + optional alignment
  // with respect to the first point.
  for (unsigned int i = 1; i < _supersystems.size(); ++i) {
    auto& system = _supersystems[i];
    if (settings.alignOrbitals) {
      auto reference = _supersystems[0];
      LocalizationTask align(system, {reference});
      align.settings = settings.wfemb.loc;
      align.settings.maxSweeps = 10;
      align.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN;
      align.run();
    }
    LocalizationTask locI(system);
    locI.settings = settings.wfemb.loc;
    locI.run();
  } // for i
}

void DOSCCTask::setUpSubsystems() {
  for (const auto& sys : _supersystems) {
    std::vector<std::shared_ptr<SystemController>> newSubsystemSet;
    for (unsigned int i = 0; i < _nFragments; ++i) {
      Settings newSettings = sys->getSettings();
      newSettings.charge = 0;
      newSettings.spin = 0;
      newSettings.name = newSettings.name + std::to_string(i);
      newSettings.identifier = newSettings.identifier + std::to_string(i);
      newSubsystemSet.push_back(std::make_shared<SystemController>(std::make_shared<Geometry>(), newSettings));
    } // for i
    _subSystems.push_back(newSubsystemSet);
  } // for sys
}

void DOSCCTask::printResults(const Eigen::VectorXd& energies, Eigen::VectorXd& relativeEnergies,
                             const SpinPolarizedData<RESTRICTED, std::vector<std::shared_ptr<DOSOrbitalGroup>>>& orbitalGroups) {
  double zeroPoint = energies(0);
  printSubSectionTitle("DOS-Multi-Level DLPNO Results");
  OutputControl::mOut << "  Name                          Total Energy     Relative Energy" << std::endl;
  OutputControl::mOut << "----------------------------------------------------------------" << std::endl;
  for (unsigned int i = 0; i < _supersystems.size(); ++i) {
    double relativeEnergy = energies[i] - zeroPoint;
    relativeEnergies(i) = relativeEnergy;
    OutputControl::mOut << "  " << std::left << std::setw(22) << std::setprecision(10)
                        << _supersystems[i]->getSystemName() << std::right << std::setw(20) << energies[i]
                        << std::setw(20) << std::setprecision(6) << relativeEnergy << std::endl;
  }
  OutputControl::mOut << "----------------------------------------------------------------" << std::endl;
  if (settings.printGroupAnalysis && orbitalGroups.size() > 0) {
    OutputControl::mOut << "  Number of unique orbital groups           " << orbitalGroups.size() << std::endl;
    OutputControl::mOut << "----------------------------------------------------------------------" << std::endl;
    OutputControl::mOut << "        index      n-orbitals(indices of sys 0)       max|delta e_ij| " << std::endl;
    OutputControl::mOut << "----------------------------------------------------------------------" << std::endl;
    for (unsigned int iGroup = 0; iGroup < orbitalGroups.size(); ++iGroup) {
      OutputControl::mOut << "  " << std::setw(10) << iGroup << std::right << std::setw(10)
                          << orbitalGroups[iGroup]->getNOrbitals() << " (";
      const std::vector<unsigned int>& indices = orbitalGroups[iGroup]->getOrbitalIndicesForSystem(0);
      for (const auto& index : indices)
        OutputControl::mOut << index << " ";
      OutputControl::mOut << ") " << std::setw(42);
      OutputControl::mOut
          << std::right
          << orbitalGroups[iGroup]->getGroupEnergyDifferenceMatrix(_pairEnergyMatrices).array().abs().maxCoeff()
          << std::endl;
    }
    OutputControl::mOut << "----------------------------------------------------------------------" << std::endl;
  } // if settings.printGroupAnalysis
}

Eigen::MatrixXd DOSCCTask::calculateMaximumDifferencePairEnergyMatrix(std::vector<std::shared_ptr<DOSOrbitalGroup>> orbitalGroups,
                                                                      const std::vector<Eigen::MatrixXd>& pairEnergyMatrices) {
  const unsigned int nOcc = pairEnergyMatrices[0].cols();
  const unsigned int nSys = pairEnergyMatrices.size();
  const unsigned int nGroups = orbitalGroups.size();
  // occ(wfemb) x groups
  const std::vector<Eigen::SparseMatrix<int>> sortingMaps =
      DOSOrbitalGroup::buildOrbitalToOrbitalGroupSortingMatrix(orbitalGroups, nOcc);
  // Transform to group-wise representation.
  std::vector<Eigen::MatrixXd> groupWisePairEnergyMatrices;
  for (unsigned int iSys = 0; iSys < nSys; ++iSys) {
    const auto sortingMatrix = sortingMaps[iSys].cast<double>();
    groupWisePairEnergyMatrices.push_back(sortingMatrix.transpose() * pairEnergyMatrices[iSys] * sortingMatrix);
  } // for iSys
  // Calculate the maximum change in pair energies in group-wise representation.
  Eigen::MatrixXd maxDifferenceMatrix = Eigen::MatrixXd::Zero(nGroups, nGroups);
  for (unsigned int iSys = 0; iSys < nSys; ++iSys) {
    for (unsigned int jSys = 0; jSys < iSys; ++jSys) {
      const Eigen::MatrixXd difference =
          (groupWisePairEnergyMatrices[iSys] - groupWisePairEnergyMatrices[jSys]).array().abs();
      maxDifferenceMatrix.array() = maxDifferenceMatrix.array().max(difference.array());
    } // for jSys
  }   // for iSys
  return maxDifferenceMatrix;
}

Eigen::MatrixXi DOSCCTask::buildSettingsMatrix(const Eigen::MatrixXd& maxDifferenceMatrix) {
  const unsigned int nGroups = maxDifferenceMatrix.cols();
  const unsigned int nThresholds = settings.pairCutOff.size();
  Eigen::MatrixXi settingsMatrix = Eigen::MatrixXi::Constant(nGroups, nGroups, nThresholds);
  for (unsigned int iGroup = 0; iGroup < nGroups; ++iGroup) {
    for (unsigned int jGroup = 0; jGroup <= iGroup; ++jGroup) {
      const double maxDeltaGroupPairEnergy = maxDifferenceMatrix(iGroup, jGroup);
      for (unsigned int iThresh = 0; iThresh < nThresholds; ++iThresh) {
        const double threshold = settings.pairCutOff[iThresh];
        if (threshold < maxDeltaGroupPairEnergy) {
          settingsMatrix(iGroup, jGroup) = iThresh;
          settingsMatrix(jGroup, iGroup) = iThresh;
          break;
        } // if threshold < maxDeltaGroupPairEnergy
      }   // for iThresh
    }     // for jGroup
  }       // for iGroup
  return settingsMatrix;
}

std::vector<unsigned int> DOSCCTask::getSkipableOrbitalGroupIndices(const Eigen::MatrixXi& settingsMatrix) {
  const unsigned int nGroups = settingsMatrix.cols();
  std::vector<unsigned int> nonSkipableSettings;
  std::vector<unsigned int> skipableGroups;
  for (unsigned int iSetting = 0; iSetting < settings.pairCutOff.size(); ++iSetting) {
    if (settings.lcPairSelection[iSetting].method != Options::PNO_METHOD::NONE)
      nonSkipableSettings.push_back(iSetting);
  } // for iSetting
  for (unsigned int iGroup = 0; iGroup < nGroups; ++iGroup) {
    bool skipGroup = true;
    for (const auto& iSetting : nonSkipableSettings) {
      if ((settingsMatrix.col(iGroup).array() == iSetting).count() > 0) {
        skipGroup = false;
        break;
      }
    } // for iSetting
    if (skipGroup)
      skipableGroups.push_back(iGroup);
  } // for iGroup
  return skipableGroups;
}

std::vector<std::vector<std::shared_ptr<OrbitalPair>>>
DOSCCTask::selectOrbitalPairsFromPairEnergies(std::vector<std::shared_ptr<DOSOrbitalGroup>> orbitalGroups,
                                              std::vector<Eigen::MatrixXd> pairEnergyMatrices) {
  const unsigned int nOcc = pairEnergyMatrices[0].cols();
  const unsigned int nSys = pairEnergyMatrices.size();
  const unsigned int nGroups = orbitalGroups.size();
  // Calculate the maximum change in pair energies in group-wise representation.
  const Eigen::MatrixXd maxDifferenceMatrix = calculateMaximumDifferencePairEnergyMatrix(orbitalGroups, pairEnergyMatrices);
  // Check vs the pair selection thresholds and assign the settings to be used in the orbital pair construction.
  const Eigen::MatrixXi settingsMatrix = buildSettingsMatrix(maxDifferenceMatrix);
  OutputControl::dOut << "Orbital-group wise settings matrix\n" << settingsMatrix << std::endl;

  // Build the orbital pairs and keep track of orbital-wise thresholds as the minimum threshold
  // assigned to a pair containing the orbital in question.
  _orbitalWiseMullikenThresholds = std::vector<Eigen::VectorXd>(nSys, Eigen::VectorXd::Constant(nOcc, 1.0));
  _orbitalToShellThresholds = std::vector<Eigen::VectorXd>(nSys, Eigen::VectorXd::Constant(nOcc, 1.0));
  _orbitalWiseDOIPAOThresholds = std::vector<Eigen::VectorXd>(nSys, Eigen::VectorXd::Constant(nOcc, 1.0));
  std::vector<std::vector<std::shared_ptr<OrbitalPair>>> systemWiseOrbitalPairs(nSys);
  const auto skipableGroups = getSkipableOrbitalGroupIndices(settingsMatrix);
  LocalCorrelationSettings diagonalDummySettings;
  diagonalDummySettings.pnoSettings = Options::PNO_SETTINGS::LOOSE;
  diagonalDummySettings.method = Options::PNO_METHOD::DLPNO_CCSD;
  diagonalDummySettings.resolvePNOSettings();
  const unsigned int maxNumberOfSettings = settings.pairCutOff.size();
  std::vector<unsigned int> numberOfOrbitalPairs(maxNumberOfSettings + 1, 0);
  for (unsigned int iSys = 0; iSys < nSys; ++iSys) {
    Eigen::MatrixXi test = Eigen::MatrixXi::Zero(nOcc, nOcc);
    const Eigen::VectorXi& coreOrbitals = _supersystems[iSys]->getActiveOrbitalController<RESTRICTED>()->getOrbitalFlags();
    for (unsigned int iGroup = 0; iGroup < nGroups; ++iGroup) {
      if (std::find(skipableGroups.begin(), skipableGroups.end(), iGroup) != skipableGroups.end()) {
        continue;
      }
      const auto& orbitalIndicesI = orbitalGroups[iGroup]->getOrbitalIndicesForSystem(iSys);
      const unsigned int minSettingsIndexIGroup = settingsMatrix.col(iGroup).minCoeff();
      for (unsigned int jGroup = 0; jGroup <= iGroup; ++jGroup) {
        if (std::find(skipableGroups.begin(), skipableGroups.end(), jGroup) != skipableGroups.end()) {
          continue;
        }
        const auto& orbitalIndicesJ = orbitalGroups[jGroup]->getOrbitalIndicesForSystem(iSys);
        const unsigned int settingsIndex = settingsMatrix(iGroup, jGroup);
        auto lcSettings = settings.lcPairSelection[settingsIndex];
        double ccsdPairThreshold = lcSettings.ccsdPairThreshold;
        const unsigned int minSettingsIndexJGroup = settingsMatrix.col(jGroup).minCoeff();
        const unsigned int minSettingsIndex = std::min(minSettingsIndexIGroup, minSettingsIndexJGroup);
        if (this->settings.strictPairEnergyThresholds) {
          ccsdPairThreshold = settings.lcPairSelection[minSettingsIndex].ccsdPairThreshold;
        }
        bool eligibleForTriples = lcSettings.method == Options::PNO_METHOD::DLPNO_CCSD_T0;
        if (this->settings.strictTriples) {
          eligibleForTriples = settings.lcPairSelection[minSettingsIndex].method == Options::PNO_METHOD::DLPNO_CCSD_T0;
        }
        for (const auto& i : orbitalIndicesI) {
          bool iIsCore = coreOrbitals[i];
          if (lcSettings.useFrozenCore && iIsCore)
            continue;
          for (const auto& j : orbitalIndicesJ) {
            // Make sure that we do not construct the orbital pairs from the diagonal groups twice.
            if (test(i, j))
              continue;
            bool jIsCore = coreOrbitals[j];
            if (lcSettings.useFrozenCore && jIsCore)
              continue;
            std::shared_ptr<OrbitalPair> newPair;
            if (lcSettings.method == Options::PNO_METHOD::NONE) {
              if (i == j) {
                newPair = std::make_shared<OrbitalPair>(i, j,
                                                        (iIsCore || jIsCore) ? diagonalDummySettings.pnoCoreScaling *
                                                                                   diagonalDummySettings.pnoThreshold
                                                                             : diagonalDummySettings.pnoThreshold,
                                                        ccsdPairThreshold, diagonalDummySettings.collinearDipoleScaling,
                                                        diagonalDummySettings.extendedDomainScaling);
                if (iSys == 0)
                  numberOfOrbitalPairs[maxNumberOfSettings] += 1;
              }
              else {
                continue;
              }
            }
            else {
              newPair = std::make_shared<OrbitalPair>(
                  i, j, (iIsCore || jIsCore) ? lcSettings.pnoCoreScaling * lcSettings.pnoThreshold : lcSettings.pnoThreshold,
                  ccsdPairThreshold, lcSettings.collinearDipoleScaling, lcSettings.extendedDomainScaling);
              if (iSys == 0)
                numberOfOrbitalPairs[settingsIndex] += 1;
            }
            newPair->eligibleForTriples = eligibleForTriples;
            newPair->scMP2PairEnergy = pairEnergyMatrices[iSys](i, j);
            test(newPair->i, newPair->j) += 1;
            if (newPair->i != newPair->j)
              test(newPair->j, newPair->i) += 1;
            systemWiseOrbitalPairs[iSys].push_back(newPair);
            if (newPair->eligibleForTriples)
              _triplesInPairSelectedCalc = true;
            // Set the orbital-wise thresholds.
            setOrbitalWiseThreshold(i, j, iSys, lcSettings);
          } // for j
        }   // for i
      }     // for jGroup
    }       // for iGroup
    if ((test.array() > 1).count() != 0) {
      std::cout << "ERROR Test matrix\n" << test << std::endl;
      throw SerenityError("ERROR: Duplicated pair detected");
    }
  } // for iSys
  OutputControl::nOut << "  setting index         number of orbital group pairs" << std::endl;
  for (unsigned int i = 0; i < maxNumberOfSettings; ++i) {
    OutputControl::nOut << "    " << i << "                                    " << numberOfOrbitalPairs[i] << std::endl;
  }
  OutputControl::nOut << "    "
                      << "kept diagonal pairs (needed for singles) " << numberOfOrbitalPairs[maxNumberOfSettings]
                      << std::endl;
  OutputControl::nOut << "    number of orbital pair groups skipped " << skipableGroups.size() << std::endl;

  return systemWiseOrbitalPairs;
}

void DOSCCTask::setOrbitalWiseThreshold(const unsigned int i, const unsigned int j, const unsigned int iSys,
                                        const LocalCorrelationSettings& lcSettings) {
  _orbitalWiseMullikenThresholds[iSys](i) = std::min(_orbitalWiseMullikenThresholds[iSys](i), lcSettings.mullikenThreshold);
  _orbitalToShellThresholds[iSys](i) = std::min(_orbitalToShellThresholds[iSys](i), lcSettings.orbitalToShellThreshold);
  _orbitalWiseDOIPAOThresholds[iSys](i) = std::min(_orbitalWiseDOIPAOThresholds[iSys](i), lcSettings.doiPAOThreshold);
  _orbitalWiseMullikenThresholds[iSys](j) = std::min(_orbitalWiseMullikenThresholds[iSys](j), lcSettings.mullikenThreshold);
  _orbitalToShellThresholds[iSys](j) = std::min(_orbitalToShellThresholds[iSys](j), lcSettings.orbitalToShellThreshold);
  _orbitalWiseDOIPAOThresholds[iSys](j) = std::min(_orbitalWiseDOIPAOThresholds[iSys](j), lcSettings.doiPAOThreshold);
}

Eigen::VectorXd DOSCCTask::runPairSelectedCC(std::vector<std::vector<std::shared_ptr<OrbitalPair>>> orbitalPairs) {
  const unsigned int nSys = _supersystems.size();
  Eigen::VectorXd energies = Eigen::VectorXd::Zero(nSys);
  for (unsigned int iSys = 0; iSys < nSys; ++iSys) {
    auto supersystem = _supersystems[iSys];
    auto eComp = supersystem->getElectronicStructure<RESTRICTED>()->getEnergyComponentController();
    // Enforce a new Fock matrix build since the one one disk is probably incorrect.
    std::shared_ptr<FockMatrix<RESTRICTED>> fockMatrix = std::make_shared<FockMatrix<RESTRICTED>>(
        supersystem->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>()->getFockMatrix(
            supersystem->getElectronicStructure<RESTRICTED>()->getDensityMatrix(), eComp));
    double totalHFEnergy = supersystem->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY);
    std::vector<std::shared_ptr<SystemController>> dummy = {};
    auto localCorrelationController = std::make_shared<LocalCorrelationController>(
        supersystem, settings.lcPairSelection[0], dummy, fockMatrix, orbitalPairs[iSys],
        _orbitalWiseMullikenThresholds[iSys], _orbitalToShellThresholds[iSys], _orbitalWiseDOIPAOThresholds[iSys]);
    DLPNO_CCSD dlpnoCCSD(localCorrelationController, settings.normThreshold, settings.maxCycles);
    dlpnoCCSD.dlpnoCCSDSettings.skipCrudePrescreening = settings.skipCrudePresPairSelected;
    Eigen::VectorXd ccsdEnergies = dlpnoCCSD.calculateElectronicEnergyCorrections();
    double triplesCorrection = 0.0;
    if (_triplesInPairSelectedCalc) {
      triplesCorrection = DLPNO_CCSD_T0::calculateEnergyCorrection(localCorrelationController);
    } // if _triplesInPairSelectedCalc
    energies(iSys) = totalHFEnergy + ccsdEnergies.sum() + triplesCorrection;
  } // for iSys
  return energies;
}

const Eigen::VectorXd& DOSCCTask::getRelativeEnergies() {
  if (!_relativeEnergies)
    this->run();
  return *_relativeEnergies;
}

const Eigen::VectorXd& DOSCCTask::getRelativeEnergiesFromPairSelected() {
  if (!_relativeEnergiesFromPairSelected) {
    if (settings.orbitalPairAnalysis) {
      this->run();
    }
    else {
      throw SerenityError("Implementation Error! This task can not perform a pair selected CC calculation. The task "
                          "settings need to be adjusted");
    }
  }
  return *_relativeEnergiesFromPairSelected;
}

const Eigen::VectorXd& DOSCCTask::getTotalEnergies() {
  if (!_totalEnergies)
    this->run();
  return *_totalEnergies;
}

const Eigen::VectorXd& DOSCCTask::getTotalEnergiesFromPairSelected() {
  if (!_totalEnergiesFromPairSelected) {
    if (settings.orbitalPairAnalysis) {
      this->run();
    }
    else {
      throw SerenityError("Implementation Error! This task can not perform a pair selected CC calculation. The task "
                          "settings need to be adjusted");
    }
  }
  return *_totalEnergiesFromPairSelected;
}

void DOSCCTaskSettings::resolveDOSSettings() {
  const std::map<Options::DOS_SETTINGS, std::vector<double>> m = {
      {Options::DOS_SETTINGS::LOOSE, {1e-1, 1e-2, 1e-3}},   {Options::DOS_SETTINGS::NORMAL, {5e-2, 5e-3, 1e-3}},
      {Options::DOS_SETTINGS::TIGHT, {2e-2, 2e-3, 8e-4}},   {Options::DOS_SETTINGS::VERY_TIGHT, {8e-3, 1e-3, 1e-4}},
      {Options::DOS_SETTINGS::EXTREME, {5e-3, 5e-4, 1e-4}}, {Options::DOS_SETTINGS::SPREAD, {1e-1, 5e-3, 1e-4}}};

  if (gdos.similarityKinEnergyThreshold.size() == 0) {
    gdos.similarityKinEnergyThreshold = m.at(dosSettings);
  }
  if (gdos.similarityLocThreshold.size() == 0) {
    gdos.similarityLocThreshold = m.at(dosSettings);
  }
}

void DOSCCTask::resolvePNOSettingsForPairSelection() {
  for (unsigned int i = 0; i < settings.pairCutOff.size(); ++i)
    settings.lcPairSelection[i].resolvePNOSettings();
}

} /* namespace Serenity */
