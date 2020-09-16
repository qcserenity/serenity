/**
 * @file LocalCorrelationController.cpp
 *
 * @date May 13, 2019
 * @author Moritz Bensberg
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
/* Include Class Header*/
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/BoughtonPulayAlgorithm.h"             //BP-Algorithm
#include "analysis/PAOSelection/DOIBasedSelecter.h"                   //DOI-based PAO selection.
#include "analysis/PAOSelection/PNOConstructor.h"                     //PNO/QCPAO construction.
#include "analysis/PAOSelection/TNOConstructor.h"                     //TNO construction.
#include "analysis/populationAnalysis/IAOPopulationCalculator.h"      //Occ to core
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h" //BP-Algorithm
#include "basis/AtomCenteredBasisController.h"                        //Atom centered basis controller.
#include "data/ElectronicStructure.h"                                 //Access to density and coefficient matrix.
#include "data/OrbitalController.h"                                   //Internal coefficients.
#include "data/OrbitalTriple.h"                                       //Orbital-triplet definition.
#include "data/PAOController.h"                                       //PAO controller.
#include "data/SparseMapsController.h"                                //Map based prescreening.
#include "data/grid/BasisFunctionOnGridControllerFactory.h"           //Pair selection/DOIs
#include "data/matrices/DensityMatrixController.h"                    //Access to the density matrix.
#include "data/matrices/MatrixInBasis.h"                              //MatrixInBasis definition.
#include "geometry/Geometry.h"                                        //Geometry definition.
#include "integrals/MO3CenterIntegralController.h"                    //Linear scaling integral calculation.
#include "integrals/OneElectronIntegralController.h"                  //Overlap integrals.
#include "io/FormattedOutputStream.h"                                 //Filtered output.
#include "memory/MemoryManager.h"                                     //Memory handling
#include "misc/Timing.h"                                              //Timings.
#include "misc/WarningTracker.h"                                      //Warnings.
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"               //Definition of a k-set.
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h"    //Overlap matrices between domains.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"                     //KL-Set definition.
#include "postHF/MPn/OrbitalPairSelector.h"                           //Pair prescreening.
#include "potentials/bundles/FDEPotentialBundleFactory.h"             //Fock matrix construction.
#include "settings/Settings.h"                                        //Settings definition.
#include "tasks/SystemAdditionTask.h" //Supersystem construction for potential reconstruction.
/* Include Std and External Headers */
#include <iomanip> //setw(...) for ostream

namespace Serenity {

LocalCorrelationController::LocalCorrelationController(std::shared_ptr<SystemController> system,
                                                       LocalCorrelationSettings settings,
                                                       std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                                       std::shared_ptr<FockMatrix<Options::SCF_MODES::RESTRICTED>> fockMatrix,
                                                       std::vector<std::shared_ptr<OrbitalPair>> initialPairs)
  : _activeSystem(system),
    _settings(settings),
    _environmentSystems(environmentSystems),
    _fock(fockMatrix),
    _allOrbitalPairs(initialPairs) {
  resolvePNOSettings();
  printSettings();
  assert(_activeSystem->getLastSCFMode() == Options::SCF_MODES::RESTRICTED);
  if (_activeSystem->getSystemContinuumModelMode() && environmentSystems.size() > 0)
    throw SerenityError((std::string) "ERROR: Embedded local correlation calculations with an explicit environment and "
                                      "a PCM are not supported yet" +
                        "       due to a lack of testing.");

  if (!_fock)
    constructFockMatrix();
  _paoController = std::make_shared<PAOController>(
      std::make_shared<DensityMatrix<RESTRICTED>>(_activeSystem->getElectronicStructure<RESTRICTED>()->getDensityMatrix()),
      std::make_shared<MatrixInBasis<RESTRICTED>>(_activeSystem->getOneElectronIntegralController()->getOverlapIntegrals()),
      _settings.paoNormalizationThreshold);
  auto paoSelecter = producePAOSelecter();
  takeTime("PAO Selection");
  _occupiedToPAOOrbitalMap = paoSelecter->selectPAOs();
  printPAOInfo(_occupiedToPAOOrbitalMap);
  timeTaken(2, "PAO Selection");
  if (initialPairs.size() == 0)
    _allOrbitalPairs = buildInitialOrbitalPairs();
  OrbitalPairSelector orbitalPairSelector(_activeSystem, _paoController);
  auto screenedPairs = orbitalPairSelector.selectOrbitalPairs(
      _allOrbitalPairs, _settings.doiNetThreshold,
      std::max(_activeSystem->getSettings().scf.canOrthThreshold, _settings.paoOrthogonalizationThreshold), _fock,
      _occupiedToPAOOrbitalMap, _settings.doiPairThreshold, _settings.collinearDipoleScaling * _settings.ccsdPairThreshold);
  _closeOrbitalPairs = screenedPairs.first;
  _veryDistantOrbitalPairs = screenedPairs.second;
  buildPAOPairDomains(_closeOrbitalPairs, *_occupiedToPAOOrbitalMap);
}

std::shared_ptr<SparseMapsController> LocalCorrelationController::getSparseMapController(bool closeOnly) {
  if (!_sparseMapsController) {
    const auto closePairs = getOrbitalPairs(OrbitalPairTypes::CLOSE);
    auto distantPairs = getOrbitalPairs(OrbitalPairTypes::DISTANT_TRIPLES);
    if (closeOnly)
      distantPairs = {};
    _sparseMapsController = std::make_shared<SparseMapsController>(
        _activeSystem, _paoController, _occupiedToPAOOrbitalMap, closePairs, distantPairs, _settings.mullikenThreshold,
        _settings.orbitalToShellThreshold, _settings.crudeStrongTripFactor * _settings.mullikenThreshold,
        _settings.crudeWeakTripFactor * _settings.mullikenThreshold);
  }
  return _sparseMapsController;
}

void LocalCorrelationController::checkSetting(double& setting, double defaultValue) {
  if (setting == std::numeric_limits<double>::infinity())
    setting = defaultValue;
}

void LocalCorrelationController::resolvePNOSettingsSet(const Eigen::VectorXd values) {
  assert(values.size() == 5);
  checkSetting(_settings.ccsdPairThreshold, values(0));
  checkSetting(_settings.pnoThreshold, values(1));
  checkSetting(_settings.orbitalToShellThreshold, values(2));
  checkSetting(_settings.mullikenThreshold, values(3));
  checkSetting(_settings.doiPAOThreshold, values(4));
}

void LocalCorrelationController::resolvePNOSettings() {
  Eigen::VectorXd pnoSettings = Eigen::VectorXd::Zero(5);
  switch (_settings.method) {
    case PNO_METHOD::DLPNO_CCSD_T0:
    case PNO_METHOD::DLPNO_CCSD:
      switch (_settings.pnoSettings) {
        case Options::PNO_SETTINGS::LOOSE:
          pnoSettings << 1e-3, 1e-6, 1e-3, 1e-3, 1e-2;
          break;
        case Options::PNO_SETTINGS::NORMAL:
          pnoSettings << 1e-4, 3.33e-7, 1e-3, 1e-3, 1e-2;
          break;
        case Options::PNO_SETTINGS::TIGHT:
          pnoSettings << 1e-5, 1e-7, 1e-3, 1e-4, 1e-3;
          break;
      }
      break;
    case PNO_METHOD::DLPNO_MP2:
      switch (_settings.pnoSettings) {
        case Options::PNO_SETTINGS::LOOSE:
          pnoSettings << 1e-3, 3.33e-7, 1e-3, 1e-4, 1e-2;
          break;
        case Options::PNO_SETTINGS::NORMAL:
          pnoSettings << 1e-4, 1e-8, 1e-3, 1e-4, 1e-3;
          break;
        case Options::PNO_SETTINGS::TIGHT:
          pnoSettings << 1e-5, 1e-9, 1e-4, 1e-5, 1e-3;
          break;
      }
      break;
  }
  resolvePNOSettingsSet(pnoSettings);
}

void LocalCorrelationController::printSettings() {
  std::string macroPNOFlag = "";
  switch (_settings.pnoSettings) {
    case Options::PNO_SETTINGS::LOOSE:
      macroPNOFlag = "LOOSE-PNO";
      break;
    case Options::PNO_SETTINGS::NORMAL:
      macroPNOFlag = "NORMAL-PNO";
      break;
    case Options::PNO_SETTINGS::TIGHT:
      macroPNOFlag = "TIGHT-PNO";
      break;
  }
  std::string methodFlag = "";
  switch (_settings.method) {
    case PNO_METHOD::DLPNO_CCSD:
      methodFlag = "DLPNO-CCSD";
      break;
    case PNO_METHOD::DLPNO_CCSD_T0:
      methodFlag = "DLPNO-CCSD(T_0)";
      break;
    case PNO_METHOD::DLPNO_MP2:
      methodFlag = "DLPNO-MP2";
      break;
  }
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
  OutputControl::nOut << "            PNO-Settings: -- " << macroPNOFlag << " --" << std::endl;
  OutputControl::vOut << "Method: " << methodFlag << std::endl;
  OutputControl::vOut << "Specific Thresholds: " << std::endl;
  OutputControl::vOut << "  T_PNO     " << _settings.pnoThreshold << std::endl;
  OutputControl::vOut << "  T_PNOSing " << _settings.pnoThreshold * _settings.singlesPNOFactor << std::endl;
  if (_settings.method == PNO_METHOD::DLPNO_CCSD_T0)
    OutputControl::vOut << "  T_TNO     " << _settings.tnoThreshold << std::endl;
  OutputControl::vOut << "  S_Core    " << _settings.pnoCoreScaling << std::endl;
  OutputControl::vOut << "  T_CCPair  " << _settings.ccsdPairThreshold << std::endl;
  if (_settings.method == PNO_METHOD::DLPNO_CCSD_T0)
    OutputControl::vOut << "  T_MP2Pair " << _settings.ccsdPairThreshold * _settings.triplesSCMP2Scaling << std::endl;
  OutputControl::vOut << "  T_DOIPair " << _settings.doiPairThreshold << std::endl;
  OutputControl::vOut << "  T_DipPair " << _settings.ccsdPairThreshold * _settings.collinearDipoleScaling << std::endl;
  OutputControl::vOut << "  T_Fitt    " << _settings.mullikenThreshold << std::endl;
  OutputControl::vOut << "  T_Shell   " << _settings.orbitalToShellThreshold << std::endl;
  OutputControl::vOut << "  T_PAO     " << _settings.doiPAOThreshold << std::endl;
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
}

std::shared_ptr<MO3CenterIntegralController> LocalCorrelationController::getMO3CenterIntegralController(bool closeOnly) {
  if (!_accurateMo3CenterIntegralController) {
    auto denseMaps = getSparseMapController(closeOnly);
    auto nOcc = _activeSystem->getNOccupiedOrbitals<RESTRICTED>();
    auto coefficients = _activeSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
    auto occCoeff = std::make_shared<Eigen::MatrixXd>(coefficients.leftCols(nOcc).eval());
    _accurateMo3CenterIntegralController = std::make_shared<MO3CenterIntegralController>(
        _activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL), _activeSystem->getBasisController(),
        denseMaps, _paoController, occCoeff, _activeSystem->getHDF5BaseName(), _activeSystem->getSettings().identifier);
  }
  return _accurateMo3CenterIntegralController;
}

std::shared_ptr<MO3CenterIntegralController> LocalCorrelationController::getApproximateMO3CenterIntegralController() {
  if (_settings.crudeDomainFactor == 0.0)
    return getMO3CenterIntegralController();
  if (!_approximateMo3CenterIntegralController) {
    auto nOcc = _activeSystem->getNOccupiedOrbitals<RESTRICTED>();
    auto coefficients = _activeSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
    auto occCoeff = std::make_shared<Eigen::MatrixXd>(coefficients.leftCols(nOcc).eval());
    const auto closePairs = getOrbitalPairs(OrbitalPairTypes::CLOSE);
    const auto distantPairs = getOrbitalPairs(OrbitalPairTypes::DISTANT);
    auto curdeMaps = std::make_shared<SparseMapsController>(
        _activeSystem, _paoController, _occupiedToPAOOrbitalMap, closePairs, distantPairs,
        _settings.mullikenThreshold * _settings.crudeDomainFactor, _settings.orbitalToShellThreshold);
    _approximateMo3CenterIntegralController = std::make_shared<MO3CenterIntegralController>(
        _activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL), _activeSystem->getBasisController(),
        curdeMaps, _paoController, occCoeff, _activeSystem->getHDF5BaseName(), _activeSystem->getSettings().identifier);
  }
  return _approximateMo3CenterIntegralController;
}

void LocalCorrelationController::printPAOInfo(const std::shared_ptr<Eigen::SparseMatrix<int>> occupiedToPAOOrbitalMap) {
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
  OutputControl::nOut << " PAO Selection" << std::endl;
  if (_settings.useBPAlgorithm) {
    OutputControl::nOut << "  Employing the Boughton--Pulay algorithm for PAO selection" << std::endl;
  }
  else {
    OutputControl::nOut << "  Employing DOI prescreening for PAO selection" << std::endl;
  }
  OutputControl::nOut << "  Number of PAOs                       " << occupiedToPAOOrbitalMap->rows() << std::endl;
  OutputControl::nOut << "  Average number of PAOs per orbital   "
                      << occupiedToPAOOrbitalMap->sum() / occupiedToPAOOrbitalMap->cols() << std::endl;
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
}

std::shared_ptr<PAOSelecter> LocalCorrelationController::producePAOSelecter() {
  auto coefficients = _activeSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  auto nOcc = _activeSystem->getNOccupiedOrbitals<RESTRICTED>();
  std::shared_ptr<PAOSelecter> paoSelecter;
  if (_settings.useBPAlgorithm) {
    auto mullikenPops = std::make_shared<Eigen::MatrixXd>(
        MullikenPopulationCalculator<RESTRICTED>::calculateAtomwiseOrbitalPopulations(
            coefficients, _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals(),
            _activeSystem->getAtomCenteredBasisController()->getBasisIndices())
            .leftCols(nOcc)
            .eval());
    auto occCoeff = std::make_shared<Eigen::MatrixXd>(coefficients.leftCols(nOcc).eval());
    paoSelecter = std::shared_ptr<PAOSelecter>(new BoughtonPulayAlgorithm(
        _activeSystem->getOneElectronIntegralController(), _activeSystem->getAtomCenteredBasisController(),
        mullikenPops, occCoeff, _settings.completenessThreshold));
  }
  else {
    auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
        _activeSystem->getSettings(), _activeSystem->getBasisController(), _activeSystem->getGridController());
    auto occCoeff = std::make_shared<Eigen::MatrixXd>(coefficients.leftCols(nOcc).eval());
    paoSelecter = std::shared_ptr<PAOSelecter>(new DOIBasedSelecter(occCoeff, _paoController, basFuncOnGridController,
                                                                    _activeSystem->getAtomCenteredBasisController(),
                                                                    _settings.doiPAOThreshold, _settings.doiNetThreshold));
  }
  return paoSelecter;
}

void LocalCorrelationController::constructFockMatrix() {
  takeTime("Fock matrix construction");
  /*
   * Create Potentials
   */
  OutputControl::nOut << "  Fock matrix construction                        ..." << std::endl;
  ;
  OutputControl::nOut.flush();
  if (iOOptions.printGridInfo && _environmentSystems.size() > 0)
    OutputControl::nOut << std::endl;
  _fock = std::make_shared<FockMatrix<RESTRICTED>>(_activeSystem->getBasisController());
  auto& f = *_fock;
  // settings
  std::shared_ptr<PotentialBundle<RESTRICTED>> activeFockConstructor;
  if (_environmentSystems.size() > 0) {
    // Prepare the densities.
    // list of environment density matrices (their controllers)
    std::vector<std::shared_ptr<DensityMatrixController<RESTRICTED>>> envDensities;
    for (auto sys : _environmentSystems) {
      if (sys->getSettings().scfMode == RESTRICTED) {
        envDensities.push_back(sys->template getElectronicStructure<RESTRICTED>()->getDensityMatrixController());
      }
      else {
        assert(false && "An unrestricted environment is not supported in local correlation approaches!");
      }
    }

    std::shared_ptr<SystemController> supersystem = nullptr;
    if (_settings.embeddingSettings.embeddingMode == Options::KIN_EMBEDDING_MODES::RECONSTRUCTION &&
        _settings.topDownReconstruction) {
      Settings supersystemSettings = _environmentSystems[0]->getSettings();
      supersystemSettings.name = "TMP-Supersystem";
      supersystemSettings.charge = 0;
      supersystemSettings.spin = 0;
      supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersystemSettings);
      std::vector<std::shared_ptr<SystemController>> fragments = _environmentSystems;
      fragments.push_back(_activeSystem);
      SystemAdditionTask<RESTRICTED> additionTask(supersystem, fragments);
      additionTask.settings.addOccupiedOrbitals = true;
      additionTask.run();
    }

    activeFockConstructor = FDEPotentialBundleFactory<RESTRICTED>::produce(
        _activeSystem, _activeSystem->getElectronicStructure<RESTRICTED>()->getDensityMatrixController(),
        _environmentSystems, envDensities, std::make_shared<EmbeddingSettings>(_settings.embeddingSettings),
        _activeSystem->getGridController(), supersystem,
        _settings.embeddingSettings.embeddingMode == Options::KIN_EMBEDDING_MODES::RECONSTRUCTION &&
            _settings.topDownReconstruction);
  }
  else {
    const auto& actsettings = _activeSystem->getSettings();
    if (actsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
      activeFockConstructor = _activeSystem->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
    }
    else if (actsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      activeFockConstructor = _activeSystem->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
    }
    else {
      std::cout << "ERROR: None existing electronicStructureTheory requested." << std::endl;
      assert(false);
    }
  }
  f = activeFockConstructor->getFockMatrix(_activeSystem->getElectronicStructure<RESTRICTED>()->getDensityMatrix(),
                                           _activeSystem->getElectronicStructure<RESTRICTED>()->getEnergyComponentController());
  OutputControl::nOut << " done" << std::endl;
  timeTaken(1, "Fock matrix construction");
}

std::vector<std::shared_ptr<OrbitalPair>> LocalCorrelationController::buildInitialOrbitalPairs() {
  takeTime("Initial orbital pairs");
  unsigned int nOcc = _activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  const Eigen::VectorXd& eigenvalues = _activeSystem->getActiveOrbitalController<RESTRICTED>()->getEigenvalues();
  std::vector<std::shared_ptr<OrbitalPair>> initialPairs;
  if (_settings.useFrozenCore) {
    OutputControl::nOut << "  The frozen core approximation is used!\n";
    OutputControl::nOut << "  Please ensure that core orbitals have been localized independently!\n";
    OutputControl::nOut << "  This may lead to artifacts in the core-orbital selection\n";
    OutputControl::nOut << "  otherwise, since it is based on an energy cut-off.\n";
    OutputControl::nOut << "  Cut off used: " << _settings.energyCutOff << " E_h" << std::endl;
  }
  unsigned int nFrozen = 0;
  for (unsigned int iOcc = 0; iOcc < nOcc; ++iOcc) {
    if (_settings.useFrozenCore && eigenvalues(iOcc) <= _settings.energyCutOff) {
      nFrozen++;
      continue;
    }
    for (unsigned int jOcc = 0; jOcc <= iOcc; ++jOcc) {
      if (_settings.useFrozenCore && eigenvalues(jOcc) <= _settings.energyCutOff)
        continue;
      auto newPair = std::make_shared<OrbitalPair>(iOcc, jOcc);
      initialPairs.push_back(newPair);
    } // for jOcc
  }   // for iOcc
  timeTaken(2, "Initial orbital pairs");
  if (_settings.useFrozenCore)
    OutputControl::nOut << "  Omitted " << nFrozen << " orbital(s) from the correlation treatment.\n" << std::endl;
  return initialPairs;
}

void LocalCorrelationController::buildSingles() {
  unsigned int nOcc = _activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  std::vector<std::shared_ptr<OrbitalPair>> closePairs = getOrbitalPairs(OrbitalPairTypes::CLOSE);
  const Eigen::MatrixXi& indices = getOrbitalPairIndices();
  _singlesIndices = std::make_shared<Eigen::VectorXi>(Eigen::VectorXi::Constant(nOcc, -1));
  auto& singleIndices = *_singlesIndices;
  unsigned int counter = 0;
  const Eigen::VectorXd& eigenvalues = _activeSystem->getActiveOrbitalController<RESTRICTED>()->getEigenvalues();
  for (unsigned int iOcc = 0; iOcc < nOcc; ++iOcc) {
    if (indices(iOcc, iOcc) >= 0) {
      auto newSingles = std::make_shared<SingleSubstitution>(closePairs[indices(iOcc, iOcc)]);
      newSingles->coreLikeOrbital = (eigenvalues(iOcc) <= _settings.energyCutOff);
      if (newSingles->coreLikeOrbital)
        OutputControl::dOut << "Core like orbital: " << newSingles->i + 1 << std::endl;
      _singles.push_back(newSingles);
      singleIndices(iOcc) = counter;
      newSingles->getDiagonalPair()->singles_i = newSingles;
      newSingles->getDiagonalPair()->singles_j = newSingles;
      ++counter;
    } // if indices(iOcc,iOcc)
  }   // for iOcc
  // Idea: Move to pair construction?!?!
  for (auto& pair : _allOrbitalPairs) {
    if (pair->type == OrbitalPairTypes::VERY_DISTANT)
      continue;
    pair->singles_i = _singles[singleIndices(pair->i)];
    pair->singles_j = _singles[singleIndices(pair->j)];
  } // for pair
}

std::shared_ptr<OrbitalPair> LocalCorrelationController::getOrbitalPair(unsigned int i, unsigned int j,
                                                                        std::vector<OrbitalPairTypes> types) {
  std::shared_ptr<OrbitalPair> pair = nullptr;
  for (const auto type : types) {
    const auto& indices = getOrbitalPairIndices(type);
    int index = indices(i, j);
    if (index < 0)
      continue;
    pair = getOrbitalPairs(type)[index];
    break;
  } // for type
  return pair;
}

std::vector<std::shared_ptr<OrbitalTriple>> LocalCorrelationController::constructOrbitalTriplets() {
  OutputControl::vOut << "Constructing orbital triples ...";
  OutputControl::vOut.flush();
  std::vector<std::shared_ptr<OrbitalTriple>> triples;
  /*
   * Construct orbital triplets from three orbital pairs ij, ik and jk, where at most one pair
   * is part of the weak pair list and no pair i==j==k is constructed.
   */
  const std::vector<OrbitalPairTypes> types = {OrbitalPairTypes::CLOSE, OrbitalPairTypes::DISTANT_TRIPLES};
  const unsigned int nOcc = _activeSystem->getNOccupiedOrbitals<RESTRICTED>();
  for (unsigned int i = 0; i < nOcc; ++i) {
    std::vector<std::shared_ptr<OrbitalPair>> ilPairs;
    for (unsigned int l = 0; l < nOcc; ++l) {
      auto ilPair = getOrbitalPair(i, l, types);
      if (ilPair)
        ilPairs.push_back(ilPair);
    }
    for (unsigned int j = 0; j <= i; ++j) {
      std::shared_ptr<OrbitalPair> ijPair = getOrbitalPair(i, j, types);
      if (!ijPair)
        continue;
      std::vector<std::shared_ptr<OrbitalPair>> jlPairs;
      for (unsigned int l = 0; l < nOcc; ++l) {
        auto jlPair = getOrbitalPair(j, l, types);
        if (jlPair)
          jlPairs.push_back(jlPair);
      }

      for (unsigned int k = 0; k <= j; ++k) {
        if (i == j && j == k)
          continue;
        std::shared_ptr<OrbitalPair> jkPair = getOrbitalPair(j, k, types);
        std::shared_ptr<OrbitalPair> ikPair = getOrbitalPair(i, k, types);
        if (!jkPair || !ikPair)
          continue;
        unsigned int nDistant = 0;
        if (ijPair->type != OrbitalPairTypes::CLOSE)
          ++nDistant;
        if (jkPair->type != OrbitalPairTypes::CLOSE)
          ++nDistant;
        if (ikPair->type != OrbitalPairTypes::CLOSE)
          ++nDistant;
        if (nDistant > 1)
          continue;
        std::vector<std::shared_ptr<OrbitalPair>> klPairs;
        for (unsigned int l = 0; l < nOcc; ++l) {
          auto klPair = getOrbitalPair(k, l, types);
          if (klPair)
            klPairs.push_back(klPair);
        }
        triples.push_back(std::make_shared<OrbitalTriple>(ikPair, jkPair, ijPair, ilPairs, jlPairs, klPairs));
      } // for k
    }   // for i
  }     // for j
  OutputControl::vnOut << " done" << std::endl;
  return triples;
}

void LocalCorrelationController::initializeSingles() {
  if (_singles.size() == 0)
    buildSingles();
  std::vector<std::shared_ptr<OrbitalPair>> closePairs = getOrbitalPairs(OrbitalPairTypes::CLOSE);
  for (auto& pair : closePairs) {
    if (pair->type != OrbitalPairTypes::CLOSE)
      continue;
    pair->singles_i->orbitalPairs.push_back(pair);
    if (pair->singles_i != pair->singles_j && pair->type == OrbitalPairTypes::CLOSE)
      pair->singles_j->orbitalPairs.push_back(pair);
  } // for pair
}

void LocalCorrelationController::buildOrbitalPairCouplingMap() {
  takeTime("IK-Pair coupling map");
  const auto& orbitalPairIndices = getOrbitalPairIndices();
  for (unsigned int iPair = 0; iPair < _closeOrbitalPairs.size(); ++iPair) {
    std::shared_ptr<OrbitalPair> pair = _closeOrbitalPairs[iPair];
    unsigned int i = pair->i;
    unsigned int j = pair->j;
    for (unsigned int k = 0; k < orbitalPairIndices.rows(); ++k) {
      if (orbitalPairIndices(i, k) < 0 || orbitalPairIndices(k, j) < 0)
        continue;
      std::shared_ptr<OrbitalPair> ikPair = _closeOrbitalPairs[orbitalPairIndices(i, k)];
      std::shared_ptr<OrbitalPair> kjPair = _closeOrbitalPairs[orbitalPairIndices(k, j)];
      if (ikPair->type == OrbitalPairTypes::VERY_DISTANT || kjPair->type == OrbitalPairTypes::VERY_DISTANT ||
          ikPair->type == OrbitalPairTypes::DISTANT || kjPair->type == OrbitalPairTypes::DISTANT)
        continue;
      auto newKSet = std::make_shared<CouplingOrbitalSet>(pair, ikPair, kjPair, k);
      pair->coupledPairs.push_back(newKSet);
    } // for k
  }   // for pair
  timeTaken(3, "IK-Pair coupling map");
}

void LocalCorrelationController::buildKLOrbitalPairs() {
  const auto& orbitalPairIndices = getOrbitalPairIndices();
  unsigned int nOcc = _activeSystem->getNOccupiedOrbitals<RESTRICTED>();
  for (unsigned int iPair = 0; iPair < _closeOrbitalPairs.size(); ++iPair) {
    std::shared_ptr<OrbitalPair> pair = _closeOrbitalPairs[iPair];
    for (unsigned int k = 0; k < nOcc; ++k) {
      for (unsigned int l = 0; l <= k; ++l) {
        if (orbitalPairIndices(k, l) < 0)
          continue;
        std::shared_ptr<OrbitalPair> klPair = _closeOrbitalPairs[orbitalPairIndices(k, l)];
        if (klPair->type == OrbitalPairTypes::VERY_DISTANT || klPair->type == OrbitalPairTypes::DISTANT)
          assert(false);
        // Skip pairs that have no domain overlap.
        if ((pair->paoDomain.transpose() * klPair->paoDomain).sum() == 0)
          continue;
        auto newKLPairSet = std::make_shared<KLOrbitalSet>(pair, klPair);
        pair->klPairSets.push_back(newKLPairSet);
      } // for iLSet
    }   // for iKSet
  }     // for pair
}

Eigen::MatrixXi LocalCorrelationController::buildOrbitalPairIndices(OrbitalPairTypes orbitalPairType) {
  const auto& orbitalPairs = getOrbitalPairs(orbitalPairType);
  auto nOcc = _activeSystem->getNOccupiedOrbitals<RESTRICTED>();
  Eigen::MatrixXi map = Eigen::MatrixXi::Constant(nOcc, nOcc, -1);
  unsigned int pairIndex = 0;
  for (const auto& pair : orbitalPairs) {
    unsigned int i = pair->i;
    unsigned int j = pair->j;
    map(i, j) = pairIndex;
    map(j, i) = pairIndex;
    pairIndex++;
  } // for pair
  return map;
}

const Eigen::MatrixXi& LocalCorrelationController::getOrbitalPairIndices(OrbitalPairTypes orbitalPairType) {
  switch (orbitalPairType) {
    case OrbitalPairTypes::CLOSE:
      if (!_orbitalPairIndices)
        _orbitalPairIndices = std::make_shared<Eigen::MatrixXi>(buildOrbitalPairIndices(orbitalPairType));
      return *_orbitalPairIndices;
    case OrbitalPairTypes::DISTANT_TRIPLES:
      if (!_distantTriplesOrbitalPairIndices)
        _distantTriplesOrbitalPairIndices = std::make_shared<Eigen::MatrixXi>(buildOrbitalPairIndices(orbitalPairType));
      return *_distantTriplesOrbitalPairIndices;
    case OrbitalPairTypes::DISTANT:
      if (!_distantOrbitalPairIndices)
        _distantOrbitalPairIndices = std::make_shared<Eigen::MatrixXi>(buildOrbitalPairIndices(orbitalPairType));
      return *_distantOrbitalPairIndices;
    case OrbitalPairTypes::VERY_DISTANT:
      if (!_veryDistantOrbitalPairIndices)
        _veryDistantOrbitalPairIndices = std::make_shared<Eigen::MatrixXi>(buildOrbitalPairIndices(orbitalPairType));
      return *_veryDistantOrbitalPairIndices;
  }
  throw SerenityError("Orbital Pair type not handled in switch.");
}

const Eigen::MatrixXi& LocalCorrelationController::getOrbitalPairIndices() {
  if (!_orbitalPairIndices) {
    _orbitalPairIndices = std::make_shared<Eigen::MatrixXi>(buildOrbitalPairIndices(OrbitalPairTypes::CLOSE));
  }
  return *_orbitalPairIndices;
}

inline void LocalCorrelationController::buildPAOPairDomains(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs,
                                                            const SparseMap& occupiedToPAOOrbitalMap) {
  for (const auto& pair : orbitalPairs)
    pair->paoDomain = (occupiedToPAOOrbitalMap.col(pair->i) + occupiedToPAOOrbitalMap.col(pair->j)).pruned().eval();
}

std::shared_ptr<PNOConstructor> LocalCorrelationController::producePNOConstructor(double ssScaling, double osScaling) {
  auto pnoConstructor = std::make_shared<PNOConstructor>(
      _activeSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients(),
      std::make_shared<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>(
          _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals()),
      _fock, _paoController,
      std::max(_activeSystem->getSettings().scf.canOrthThreshold, _settings.paoOrthogonalizationThreshold),
      _settings.pnoThreshold, _settings.singlesPNOFactor, _settings.pnoCoreScaling, _environmentSystems,
      (_settings.projectedEnvironment) ? _settings.embeddingSettings.levelShiftParameter : 0.0, _settings.setFaiZero,
      ssScaling, osScaling);
  return pnoConstructor;
}

std::shared_ptr<QuasiCanonicalPAODomainConstructor> LocalCorrelationController::produceQCPAOConstructor(double ssScaling,
                                                                                                        double osScaling) {
  auto paoConstructor = std::make_shared<QuasiCanonicalPAODomainConstructor>(
      _activeSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients(),
      std::make_shared<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>(
          _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals()),
      _fock, _paoController,
      std::max(_activeSystem->getSettings().scf.canOrthThreshold, _settings.paoOrthogonalizationThreshold),
      _environmentSystems, _settings.embeddingSettings.levelShiftParameter, ssScaling, osScaling);
  return paoConstructor;
}

std::shared_ptr<TNOConstructor> LocalCorrelationController::produceTNOConstructor() {
  auto tnoConstructor = std::make_shared<TNOConstructor>(
      std::make_shared<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>(
          _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals()),
      _fock, _paoController,
      std::max(_activeSystem->getSettings().scf.canOrthThreshold, _settings.paoOrthogonalizationThreshold),
      _settings.tnoThreshold, (_settings.useTriplesCoreScaling) ? _settings.pnoCoreScaling : 1.0, _environmentSystems,
      (_settings.projectedEnvironment) ? _settings.embeddingSettings.levelShiftParameter : 0.0);
  return tnoConstructor;
}

void LocalCorrelationController::selectDistantOrbitalPairs() {
  _orbitalPairIndices = nullptr;
  _distantOrbitalPairIndices = nullptr;
  _orbitalTriples = {};
  _distantTriplesOrbitalPairIndices = nullptr;
  double veryDistantThreshold = _settings.collinearDipoleScaling * _settings.ccsdPairThreshold;
  double triplesDistantCutOff = _settings.ccsdPairThreshold * _settings.triplesSCMP2Scaling;
  std::vector<std::shared_ptr<OrbitalPair>> newClosePairs = {};
  OutputControl::vOut << "Selecting distant orbital pairs based on SC-MP2 pair energy." << std::endl;
  OutputControl::vOut << "Note that diagonal pairs will always be retained in the calculation!" << std::endl;
  for (auto& pair : _closeOrbitalPairs) {
    bool diagonalPair = pair->i == pair->j;
    if ((std::fabs(pair->scMP2PairEnergy) < veryDistantThreshold || pair->type == OrbitalPairTypes::VERY_DISTANT) &&
        not diagonalPair) {
      _veryDistantOrbitalPairs.push_back(pair);
      pair->t_ij = Eigen::MatrixXd(0, 0);
      pair->k_ij = Eigen::MatrixXd(0, 0);
      pair->type = OrbitalPairTypes::VERY_DISTANT;
      continue;
    }
    if (std::fabs(pair->scMP2PairEnergy) < _settings.ccsdPairThreshold && not diagonalPair) {
      if (pair->type != OrbitalPairTypes::VERY_DISTANT)
        pair->type = OrbitalPairTypes::DISTANT;
      _distantOrbitalPairs.push_back(pair);
      pair->type = OrbitalPairTypes::DISTANT;
      if (std::fabs(pair->scMP2PairEnergy) >= triplesDistantCutOff) {
        _distantOrbitalPairsTriples.push_back(pair);
        pair->type = OrbitalPairTypes::DISTANT_TRIPLES;
      }
    }
    else {
      newClosePairs.push_back(pair);
    } // if/else pair->scMP2PairEnergy < _settings.ccsdPairThreshold
  }   // for pair
  _closeOrbitalPairs = newClosePairs;
  if (_closeOrbitalPairs.size() < 1)
    throw SerenityError("Local-Correlation: All orbital pairs have been screened away! Adjust the SC-MP2 screening "
                        "threshold <ccsdPairThreshold>!");
}

std::vector<OrbitalPairSet> LocalCorrelationController::getCloseOrbitalPairSets() {
  if (_closeOrbitalPairSets.size() == 0) {
    double usedMemory = 0.0;
    double totalMemoty = 0.0;
    double maximumMemory = _settings.maximumMemoryRatio * MemoryManager::getInstance()->getAvailableSystemMemory();
    OutputControl::nOut << std::fixed;
    double overlapMatrixSize = getDomainOverlapMatrixController()->getOverlapMatrixSize();
    OutputControl::nOut << "  PNO-Overlap matrix size " << setw(26) << overlapMatrixSize * 1e-6 << " MB" << std::endl;
    OutputControl::nOut << "  Memory per block        " << setw(26) << maximumMemory * 1e-9 << " GB" << std::endl;

    OrbitalPairSet newOrbitalPairSet;
    for (auto pair : _closeOrbitalPairs) {
      newOrbitalPairSet.push_back(pair);
      usedMemory += pair->getMemoryRequirement(_settings.method == PNO_METHOD::DLPNO_MP2);
      if (usedMemory >= maximumMemory) {
        _closeOrbitalPairSets.push_back(newOrbitalPairSet);
        newOrbitalPairSet = {};
        totalMemoty += usedMemory;
        usedMemory = 0.0;
      } // if usedMemory >= maximumMemory
    }   // for _closeOrbitalPairs
    totalMemoty += usedMemory;
    _closeOrbitalPairSets.push_back(newOrbitalPairSet);
    OutputControl::nOut << "  Memory for integrals    " << setw(26) << totalMemoty * 1e-9 << " GB" << std::endl;
    OutputControl::nOut << "  Orbital pairs are divided into " << _closeOrbitalPairSets.size() << " set(s)." << std::endl;
    OutputControl::nOut << std::scientific;
    OutputControl::nOut << "-----------------------------------------------------" << std::endl;
    if (_closeOrbitalPairSets.size() > 1) {
      WarningTracker::printWarning(
          (std::string) "Warning: Writing and reading of integrals from file in DLPNO-CCSD is not fully supported "
                        "yet!" +
              " This may cause a crash!\n Solutions:\n" + "        * Increase the 'maximumMemoryRatio' to 1.\n" +
              "        * Use a computer with more main memory." + "        * Remember to localize the orbitals." +
              "        * Use the frozen core approximation 'useFrozenCore true'.",
          true);
    }
  }
  return _closeOrbitalPairSets;
}

std::shared_ptr<DomainOverlapMatrixController> LocalCorrelationController::getDomainOverlapMatrixController() {
  if (!_domainOverlapMatrixController) {
    _domainOverlapMatrixController = std::make_shared<DomainOverlapMatrixController>(
        _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals(), getPAOController(),
        getOrbitalPairs(OrbitalPairTypes::CLOSE), _singles, getOrbitalPairIndices(),
        _activeSystem->getNOccupiedOrbitals<RESTRICTED>());
  }
  return _domainOverlapMatrixController;
}

} /* namespace Serenity */
