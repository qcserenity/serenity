/**
 * @file ActiveSpaceSelectionTask.cpp
 *
 * @date Sep 11, 2018
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
#include "tasks/ActiveSpaceSelectionTask.h"
/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/OrbitalAligner.h"
#include "analysis/populationAnalysis/IAOPopulationCalculator.h"
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutputStream.h"
#include "misc/SystemSplittingTools.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"
#include "tasks/ScfTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ActiveSpaceSelectionTask<SCFMode>::ActiveSpaceSelectionTask(std::vector<std::shared_ptr<SystemController>> supersystems,
                                                            std::vector<std::shared_ptr<SystemController>> acitveSystems,
                                                            std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _supersystems(supersystems), _activeSystems(acitveSystems), _environmentSystems(environmentSystems) {
  // Sanity checks
  if (_supersystems.size() <= 1)
    throw SerenityError("Comparison between systems is only sensible if there are more than one!");
  // Check atom ordering
  auto firstGeom = _supersystems[0]->getGeometry();
  auto firstAtomList = firstGeom->getAtoms();
  const auto& firstOccupations = _supersystems[0]->getNOccupiedOrbitals<SCFMode>();
  for (unsigned int sys = 1; sys < _supersystems.size(); ++sys) {
    auto atomList = _supersystems[sys]->getGeometry()->getAtoms();
    if (atomList.size() != firstAtomList.size())
      throw SerenityError("The number of atoms have to be identical in all structures!");
    for (unsigned int iAtom = 0; iAtom < firstAtomList.size(); ++iAtom) {
      if (atomList[iAtom]->getAtomType()->getName() != firstAtomList[iAtom]->getAtomType()->getName())
        throw SerenityError("The atoms have to be ordered in the same way for all structures!");
    }
    // Check the number of electrons in the subsystems
    const auto& occupations = _supersystems[sys]->getNOccupiedOrbitals<SCFMode>();
    for_spin(firstOccupations, occupations) {
      if (firstOccupations_spin != occupations_spin)
        throw SerenityError("Inconsistent number of electrons along the reaction coordinate!");
    }; // for_spin
  }
  // Initialize the orbital map
  for (auto& sys : _supersystems) {
    auto nOcc = sys->getNOccupiedOrbitals<SCFMode>();
    SpinPolarizedData<SCFMode, Eigen::VectorXi> hasPartner;
    for_spin(nOcc, hasPartner) {
      hasPartner_spin = Eigen::VectorXi::Constant(nOcc_spin, 1);
    };
    _unpairedOrbitals.push_back(hasPartner);
    std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>> mapVector;
    for (unsigned int j = 0; j < _supersystems.size(); ++j) {
      SpinPolarizedData<SCFMode, Eigen::MatrixXi> newMap;
      auto nOccJ = _supersystems[j]->getNOccupiedOrbitals<SCFMode>();
      for_spin(newMap, nOccJ, nOcc) {
        newMap_spin = Eigen::MatrixXi::Constant(nOcc_spin, nOccJ_spin, 0);
      };
      mapVector.push_back(newMap);
    }
    _completeOrbitalMap.push_back(mapVector);
  }
  auto supersystemBasisLabel = _supersystems[0]->getSettings().basis.label;
  for (unsigned int iSys = 0; iSys < _supersystems.size(); ++iSys) {
    if (_supersystems[iSys]->getSettings().basis.label != supersystemBasisLabel)
      throw SerenityError("ERROR: The basis-set label for all supersystems has to be identical!");
    if (_activeSystems.size() > iSys) {
      if (_activeSystems[iSys]->getSettings().basis.label != supersystemBasisLabel)
        throw SerenityError(
            (std::string) "ERROR: The basis-set label for supersystem and subsystems has to be identical!" +
            "       Check the basis-set label for subsystem " + _activeSystems[iSys]->getSystemName());
    }
    if (_environmentSystems.size() > iSys) {
      if (_environmentSystems[iSys]->getSettings().basis.label != supersystemBasisLabel)
        throw SerenityError(
            (std::string) "ERROR: The basis-set label for supersystem and subsystems has to be identical!" +
            "       Check the basis-set label for subsystem " + _environmentSystems[iSys]->getSystemName());
    }
  }
}
template<Options::SCF_MODES SCFMode>
std::shared_ptr<SystemController>
ActiveSpaceSelectionTask<SCFMode>::initializeSubsystem(std::shared_ptr<SystemController> supersystemController,
                                                       std::shared_ptr<Geometry> geometry, bool active, int charge, int spin) {
  std::string namePostfix = "_Env";
  if (active)
    namePostfix = "_Act";
  Settings subSettings = supersystemController->getSettings();
  subSettings.charge = charge;
  subSettings.spin = spin;
  subSettings.name = subSettings.name + namePostfix;
  subSettings.path =
      subSettings.path.substr(0, subSettings.path.size() - (supersystemController->getSystemName().size() + 1));
  subSettings.identifier = subSettings.identifier + namePostfix;
  return std::make_shared<SystemController>(geometry, subSettings);
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::run() {
  // Generate comparison criteria between orbitals
  // and run SCF calculation for every system.
  if (!settings.skipLocalization)
    prepareOrbitals();
  produceComparisonCriteria();
  // Create mapping between the orbitals from similarity criteria
  compareOrbitals();

  for (unsigned int i = 0; i < _supersystems.size(); ++i) {
    const auto& unpairedVectorI = _unpairedOrbitals[i];
    OutputControl::vOut << "Active space selection for system " << _supersystems[i]->getSystemName() << std::endl;
    for_spin(unpairedVectorI) {
      OutputControl::vOut << unpairedVectorI_spin.transpose() << std::endl;
    };
  }

  /* Build Systems */
  for (unsigned int isys = 0; isys < _supersystems.size(); ++isys) {
    auto system = _supersystems[isys];
    auto& partnerVector = _unpairedOrbitals[isys];
    // Split the electronic structure
    SpinPolarizedData<SCFMode, std::vector<bool>> activeOrbitals;
    for_spin(partnerVector, activeOrbitals) {
      for (unsigned int iOrb = 0; iOrb < partnerVector_spin.size(); ++iOrb) {
        if (!partnerVector_spin(iOrb)) {
          activeOrbitals_spin.push_back(true);
        }
        else {
          activeOrbitals_spin.push_back(false);
        }
      }
    };
    auto electronicStructurePair = SystemSplittingTools<SCFMode>::splitElectronicStructure(system, activeOrbitals);
    // Build geometry from orbital localization
    auto geometryPair =
        SystemSplittingTools<SCFMode>::splitGeometry(system, electronicStructurePair.first, settings.localizationThreshold);
    // Occupations
    auto nOccPair = getOrbitalOccupations(partnerVector);

    // Calculate the charge of the subsystems
    int actNuclearCharge = 0;
    int envNuclearCharge = 0;
    for (const auto& atom : geometryPair.first->getAtoms())
      actNuclearCharge += atom->getEffectiveCharge();
    for (const auto& atom : geometryPair.second->getAtoms())
      envNuclearCharge += atom->getEffectiveCharge();
    double occupation = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
    int actCharge = actNuclearCharge - occupation * nOccPair.first.total();
    int envCharge = envNuclearCharge - occupation * nOccPair.second.total();
    // Sanity check. All active system must contain the same number of electrons!
    assertOccupationConsistency(nOccPair);
    std::shared_ptr<SystemController> environmentSystemController;
    std::shared_ptr<SystemController> activeSystemController;
    if (GLOBAL_PRINT_LEVEL != Options::GLOBAL_PRINT_LEVELS::MINIMUM)
      printSectionTitle((std::string) "DOS for " + system->getSystemName());
    if (isys + 1 > _activeSystems.size()) {
      activeSystemController = initializeSubsystem(system, geometryPair.first, true, actCharge, getSpin(nOccPair.first));
    }
    else {
      activeSystemController = _activeSystems[isys];
      while (activeSystemController->getGeometry()->getNAtoms() > 0)
        activeSystemController->getGeometry()->deleteAtom(0);
      *activeSystemController->getGeometry() += *geometryPair.first;
      activeSystemController->setCharge(actCharge);
      activeSystemController->setSpin(getSpin(nOccPair.first));
      activeSystemController->getGeometry()->printToFile(activeSystemController->getHDF5BaseName(),
                                                         activeSystemController->getSettings().identifier);
      activeSystemController->print();
    }
    printSystemInformations(activeSystemController);
    if (isys + 1 > _environmentSystems.size()) {
      environmentSystemController =
          initializeSubsystem(system, geometryPair.second, false, envCharge, getSpin(nOccPair.second));
    }
    else {
      environmentSystemController = _environmentSystems[isys];
      while (environmentSystemController->getGeometry()->getNAtoms() > 0)
        environmentSystemController->getGeometry()->deleteAtom(0);
      *environmentSystemController->getGeometry() += *geometryPair.second;
      environmentSystemController->setCharge(envCharge);
      environmentSystemController->setSpin(getSpin(nOccPair.second));
      environmentSystemController->getGeometry()->printToFile(environmentSystemController->getHDF5BaseName(),
                                                              environmentSystemController->getSettings().identifier);
      environmentSystemController->print();
    }
    printSystemInformations(environmentSystemController);
    // Set electronic structures
    activeSystemController->setBasisController(system->getAtomCenteredBasisController());
    activeSystemController->setElectronicStructure(electronicStructurePair.first);
    environmentSystemController->setBasisController(system->getAtomCenteredBasisController());
    environmentSystemController->setElectronicStructure(electronicStructurePair.second);
    assert(activeSystemController->getHDF5BaseName() != system->getHDF5BaseName());
    // Save to file
    std::cout << "Act-base name: " << activeSystemController->getHDF5BaseName() << std::endl;
    std::cout << "Env-base name: " << environmentSystemController->getHDF5BaseName() << std::endl;
    electronicStructurePair.first->getDensityMatrixController()->updateDensityMatrix();
    electronicStructurePair.second->getDensityMatrixController()->updateDensityMatrix();
    electronicStructurePair.first->toHDF5(activeSystemController->getHDF5BaseName(),
                                          activeSystemController->getSettings().identifier);
    electronicStructurePair.second->toHDF5(environmentSystemController->getHDF5BaseName(),
                                           environmentSystemController->getSettings().identifier);
    Settings actSettings = activeSystemController->getSettings();
    Settings envSettings = environmentSystemController->getSettings();
    actSettings.printSettings();
    envSettings.printSettings();
    if (keepSystemPairs) {
      _systemPairs.push_back(std::pair<std::shared_ptr<SystemController>, std::shared_ptr<SystemController>>(
          activeSystemController, environmentSystemController));
    }
  }
}

/* ************************** */
/*      Helper Functions      */
/* ************************** */
template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::printSystemInformations(std::shared_ptr<SystemController> system) {
  OutputControl::nOut << "Number of occupied orbitals: ";
  const auto& occupations = system->getNOccupiedOrbitals<SCFMode>();
  for_spin(occupations) {
    OutputControl::nOut << occupations_spin << " ";
  };
  OutputControl::nOut << std::endl;
  OutputControl::nOut << "Spin:                        " << system->getSpin() << std::endl;
  OutputControl::nOut << "Charge:                      " << system->getCharge() << std::endl;
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::assertOccupationConsistency(
    std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>> nOccPair) {
  if (!_finalOccupation)
    _finalOccupation =
        std::make_shared<std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>>(nOccPair);
  const auto& finalOccAct = _finalOccupation->first;
  const auto& finalOccEnv = _finalOccupation->second;
  const auto& currentOccAct = nOccPair.first;
  const auto& currentOccEnv = nOccPair.second;
  for_spin(finalOccAct, finalOccEnv, currentOccAct, currentOccEnv) {
    if (currentOccAct_spin != finalOccAct_spin || currentOccEnv_spin != finalOccEnv_spin)
      throw SerenityError("The occupation numbers between the systems are not consistent!");
  };
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::prepareOrbitals() {
  for (auto& sys : _supersystems) {
    if (!settings.load) {
      ScfTask<SCFMode> scfTask(sys);
      scfTask.run();
    }
    else {
      if (!sys->hasElectronicStructure<SCFMode>())
        throw SerenityError((std::string) "No electronic structure available. However load=true was set! System " +
                            sys->getSystemName());
    } // else if !settings.load
  }   // for sys
  {
    LocalizationTask locTask(_supersystems[0]);
    locTask.settings.locType = settings.locType;
    locTask.settings.splitValenceAndCore = settings.splitValenceAndCore;
    locTask.run();
  }

  if (settings.alignPiOrbitals)
    alignPiOrbitals();

  for (unsigned int i = 1; i < _supersystems.size(); ++i) {
    LocalizationTask locTask(_supersystems[i]);
    locTask.settings.locType = settings.locType;
    locTask.settings.splitValenceAndCore = settings.splitValenceAndCore;
    locTask.run();
  }
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::buildOrbitalPopulations() {
  _orbitalPopulations.resize(0);
  for (auto& sys : _supersystems) {
    SPMatrix<SCFMode> orbitalPopulations;
    if (settings.populationAlgorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN) {
      orbitalPopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
          sys->getActiveOrbitalController<SCFMode>()->getCoefficients(),
          sys->getOneElectronIntegralController()->getOverlapIntegrals(),
          sys->getAtomCenteredBasisController()->getBasisIndices());
    }
    else if (settings.populationAlgorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::IAO) {
      orbitalPopulations = IAOPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(sys);
    }
    else if (settings.populationAlgorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell) {
      // This is basically the square of the CIAO-coefficient summed over a shell.
      orbitalPopulations = IAOPopulationCalculator<SCFMode>::calculateShellwiseOrbitalPopulations(sys);
    }
    else {
      throw SerenityError("The algorithm used for evaluating orbital-wise populations is not supported.");
    }
    SPMatrix<SCFMode> occPopulations;
    const auto nOcc = sys->getNOccupiedOrbitals<SCFMode>();
    for_spin(occPopulations, orbitalPopulations, nOcc) {
      occPopulations_spin = orbitalPopulations_spin.block(0, 0, orbitalPopulations_spin.rows(), nOcc_spin);
    };
    _orbitalPopulations.push_back(occPopulations);
  } /* for sys : _supersystems */
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::reduceOrbitalMap() {
  // Reduce the orbital maps to a single vector for each system
  for (unsigned int i = 0; i < _supersystems.size(); ++i) {
    auto& unpairedVectorI = _unpairedOrbitals[i];
    for (unsigned int j = 0; j < _supersystems.size(); ++j) {
      if (i == j)
        continue; // i==j is always the identity
      // rowwise sum for ij-map reduction
      auto& mapIJ = _completeOrbitalMap[i][j];
      for_spin(mapIJ, unpairedVectorI) {
        Eigen::VectorXi ijImportant = mapIJ_spin.rowwise().sum();
        unpairedVectorI_spin = unpairedVectorI_spin.cwiseProduct(ijImportant);
      };
    } // for j
  }   // for i
  // Back mapping to ensure consistency in the orbital selection.
  backMapping();
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::orbitalClustering() {
  for (unsigned int i = 0; i < _supersystems.size(); ++i) {
    auto sys = _supersystems[i];
    auto& unpairedVectorI = _unpairedOrbitals[i];
    // Remove isolated MOs
    // Build orbital-wise IAO charges (if not already present)
    auto orbitalWiseCharges = _orbitalPopulations[i];
    if (settings.populationAlgorithm != Options::POPULATION_ANALYSIS_ALGORITHMS::IAO) {
      orbitalWiseCharges = IAOPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(sys);
    }
    for_spin(orbitalWiseCharges, unpairedVectorI) {
      // Build atom to orbital map
      Eigen::MatrixXi atomToOrbitalMap =
          Eigen::MatrixXi::Zero(orbitalWiseCharges_spin.rows(), orbitalWiseCharges_spin.cols());
      for (unsigned int iOrb = 0; iOrb < orbitalWiseCharges_spin.cols(); ++iOrb) {
        for (unsigned int iAtom = 0; iAtom < orbitalWiseCharges_spin.rows(); ++iAtom) {
          if (std::fabs(orbitalWiseCharges_spin(iAtom, iOrb)) > settings.clusterThreshold) {
            atomToOrbitalMap(iAtom, iOrb) = 1;
          } // if clusterThreshold
        }   // iAtom
      }     // iOrb
      // Loop over environment orbitals
      for (unsigned int iOrb = 0; iOrb < unpairedVectorI_spin.size(); ++iOrb) {
        if (unpairedVectorI_spin(iOrb)) {
          Eigen::VectorXi atomToOrbitalVector = atomToOrbitalMap.col(iOrb);
          bool connectedOrbital = false;
          for (unsigned int jOrb = 0; jOrb < unpairedVectorI_spin.size(); ++jOrb) {
            if (unpairedVectorI_spin(jOrb)) {
              int test = atomToOrbitalVector.cwiseProduct(atomToOrbitalMap.col(jOrb)).sum();
              if (test != 0)
                connectedOrbital = true;
            } // if jOrb-->env
          }   // for jOrb
          // The orbital connection is strictly binary! Thus, we are fine with setting this
          // to zero without checking the connection of the other orbitals again.
          if (!connectedOrbital)
            unpairedVectorI_spin(iOrb) = 0;
        } // if iOrb-->env
      }   // for iOrb
    };    // for_spin
  }       // for i
  backMapping();
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::backMapping() {
  // Back mapping. Every orbital that is not found at every point along
  // the reaction coordinate is set to 0
  while (true) {
    bool occChanged = false;
    for (unsigned int i = 0; i < _supersystems.size(); ++i) {
      auto& unpairedVectorI = _unpairedOrbitals[i];
      for (unsigned int j = 0; j < _supersystems.size(); ++j) {
        if (i == j && !settings.checkDegeneracies)
          continue; // i==j is always the identity
        auto& unpairedVectorJ = _unpairedOrbitals[j];
        auto& mapIJ = _completeOrbitalMap[i][j];
        for_spin(mapIJ, unpairedVectorI, unpairedVectorJ) {
          for (unsigned int iOrb = 0; iOrb < unpairedVectorI_spin.size(); ++iOrb) {
            if (unpairedVectorI_spin(iOrb)) {
              // get partner/partners in j
              auto ijPartnerMap = mapIJ_spin.row(iOrb);
              for (unsigned int jOrb = 0; jOrb < ijPartnerMap.size(); ++jOrb) {
                if (ijPartnerMap(jOrb)) {
                  // Check whether jOrb is not active
                  if (!unpairedVectorJ_spin(jOrb)) {
                    // If this is reached the orbital is considered active in one and passive in an other system.
                    // --> Set the orbital to active!
                    unpairedVectorI_spin(iOrb) = 0;
                    occChanged = true;
                    if (i == j) {
                      OutputControl::vOut << "Occupation was changed due to degeneracy of orbitals! "
                                          << _supersystems[i]->getSystemName() << " " << iOrb << " <-> " << jOrb
                                          << std::endl;
                    }
                  } // if !unpairedVectorJ_spin(jOrb)
                  /*
                   * Check for consistent handling of degeneracy!
                   * Consider two systems. In on system a set of orbitals is degenerate (Orb1, Orb2, Orb3) and
                   * identically localized (e.g. 2p^3 in the core region of a heavy atom). This degeneracy is broken in
                   * the reaction. However, at least one orbital (OrbU) is not changed in the process. The resulting
                   * partner map would look like this: Orb1, Orb2, Orb3 -> 1; OrbU -> 3 Thus no orbital would be
                   * considered active even though at least two of them are changed! Solution: If the orbital map is
                   * inconsistent, meaning that a 3 maps to three 1s or a 2 to two 1s. All of these orbitals have to be
                   * considered as active!
                   * --> The number of mapped orbitals has to be identical on both sides of the map!
                   */
                  if (unpairedVectorI_spin(iOrb) != unpairedVectorJ_spin(jOrb)) {
                    unpairedVectorI_spin(iOrb) = 0;
                    occChanged = true;
                  }
                } // if ijPartnerMap(jOrb)
              }   // for jOrb
            }     // if unpairedVectorI_spin(iOrb)
          }       // for iOrb
        };        // for_spin
      }           // for j
    }             // for i
    if (!occChanged)
      break;
  } // while true
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::buildOrbitalMap() {
  unsigned int nPopulations = _orbitalPopulations.size();
  for (unsigned int i = 0; i < nPopulations; ++i) {
    // Diagonal map is always fulfilled
    auto& iiMap = _completeOrbitalMap[i][i];
    for_spin(iiMap) {
      auto nOcci = iiMap_spin.cols();
      iiMap_spin.diagonal() = Eigen::VectorXi::Constant(nOcci, 1);
    };
    unsigned int jStart = (settings.checkDegeneracies) ? i : i + 1;
    for (unsigned int j = jStart; j < nPopulations; ++j) {
      auto& mapIJ = _completeOrbitalMap[i][j];
      auto& mapJI = _completeOrbitalMap[j][i];
      const auto& popsI = _orbitalPopulations[i];
      const auto& popsJ = _orbitalPopulations[j];
      const auto& kinEI = _kineticEnergies[i];
      const auto& kinEJ = _kineticEnergies[j];

      for_spin(popsI, popsJ, mapIJ, mapJI) {
        auto nOccI = mapIJ_spin.rows();
        auto nOccJ = mapIJ_spin.cols();
        OutputControl::dOut << "Best match loc " << _supersystems[i]->getSystemName() << "<->"
                            << _supersystems[j]->getSystemName() << std::endl;
        for (unsigned int k = 0; k < nOccI; ++k) {
          double bestMatch = 1.0;
          for (unsigned int l = 0; l < nOccJ; ++l) {
            double populationMeasure = (popsI_spin.col(k) - popsJ_spin.col(l)).array().abs().sum();
            double threshold = (i == j) ? settings.degeneracyFactor * settings.similarityLocThreshold
                                        : settings.similarityLocThreshold;
            if (settings.usePiBias && i != j) {
              int nMajorCont_k = (popsI_spin.col(k).array().abs() >= settings.biasThreshold).count();
              int nMajorCont_l = (popsJ_spin.col(l).array().abs() >= settings.biasThreshold).count();
              double bias = nMajorCont_k * nMajorCont_l / settings.biasAverage;
              OutputControl::vOut << "Bias: " << nMajorCont_k << " " << nMajorCont_l << " " << bias << std::endl;
              threshold *= bias;
            }
            if (populationMeasure < bestMatch)
              bestMatch = populationMeasure;
            if (populationMeasure < threshold) {
              mapIJ_spin(k, l) = 1;
              mapJI_spin(l, k) = 1;
            }
          }
          OutputControl::dOut << k << "  " << bestMatch << std::endl;
        }
        OutputControl::dOut << "Best match loc END" << std::endl << std::endl;
        OutputControl::dOut << "Population map for system " << _supersystems[i]->getSystemName() << " (rows) to system "
                            << _supersystems[j]->getSystemName() << " (columns)" << std::endl;
        OutputControl::dOut << mapIJ_spin << std::endl;
      }; /* for_spin */
      OutputControl::dOut << "Kinetic energy map " << _supersystems[i]->getSystemName() << " (rows) to system "
                          << _supersystems[j]->getSystemName() << " (columns)" << std::endl;
      for_spin(kinEI, kinEJ, mapIJ, mapJI) {
        auto nOccI = mapIJ_spin.rows();
        auto nOccJ = mapIJ_spin.cols();
        for (unsigned int k = 0; k < nOccI; ++k) {
          double bestMatch = 10;
          for (unsigned int l = 0; l < nOccJ; ++l) {
            // inverse check to remove wrong matches from population analysis
            double kineticEnergyMeasure = std::fabs(kinEI_spin[k] - kinEJ_spin[l]);
            double threshold = (i == j) ? settings.degeneracyFactor * settings.similarityKinEnergyThreshold
                                        : settings.similarityKinEnergyThreshold;
            OutputControl::dOut << (kineticEnergyMeasure) << " ";
            if (kineticEnergyMeasure > threshold) {
              mapIJ_spin(k, l) = 0;
              mapJI_spin(l, k) = 0;
            } // if
            else {
              if (mapIJ_spin(k, l) && mapJI_spin(l, k)) {
                if (kineticEnergyMeasure < bestMatch)
                  bestMatch = kineticEnergyMeasure;
              }
            }
          } // for l
          OutputControl::dOut << std::endl;
        } // for k
      };  /* for_spin */
    }     /* system j */
  }       /* system i */
}

template<Options::SCF_MODES SCFMode>
std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>
ActiveSpaceSelectionTask<SCFMode>::getOrbitalOccupations(SpinPolarizedData<SCFMode, Eigen::VectorXi>& partnerVector) {
  SpinPolarizedData<SCFMode, unsigned int> nOccAct(0);
  SpinPolarizedData<SCFMode, unsigned int> nOccEnv(0);
  for_spin(nOccAct, nOccEnv, partnerVector) {
    for (unsigned int i = 0; i < partnerVector_spin.size(); ++i) {
      if (partnerVector_spin(i))
        ++nOccEnv_spin;
    }
    nOccAct_spin = partnerVector_spin.size() - nOccEnv_spin;
  };
  return std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>(nOccAct, nOccEnv);
}
template<>
int ActiveSpaceSelectionTask<Options::SCF_MODES::RESTRICTED>::getSpin(
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int> nOcc) {
  (void)nOcc;
  return 0;
}
template<>
int ActiveSpaceSelectionTask<Options::SCF_MODES::UNRESTRICTED>::getSpin(
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int> nOcc) {
  return nOcc.alpha - nOcc.beta;
}
template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::calculateKineticEnergy() {
  _kineticEnergies = {};
  auto libint = Libint::getSharedPtr();
  for (unsigned int i = 0; i < _supersystems.size(); ++i) {
    const auto& system = _supersystems[i];
    // Calculate the kinetic energies of the occupied orbitals
    const auto& coeff = system->getActiveOrbitalController<SCFMode>()->getCoefficients();
    auto kinIntegrals = libint->compute1eInts(libint2::Operator::kinetic, system->getBasisController());
    auto nOcc = system->getNOccupiedOrbitals<SCFMode>();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> kineticEnergies;
    for_spin(nOcc, coeff, kineticEnergies) {
      kineticEnergies_spin.resize(nOcc_spin);
      for (unsigned int occI = 0; occI < nOcc_spin; ++occI) {
        const auto& coeffI = coeff_spin.col(occI);
        kineticEnergies_spin[occI] = coeffI.transpose() * kinIntegrals * coeffI;
      } // for occI
    };
    // save the kinetic energies.
    _kineticEnergies.push_back(kineticEnergies);
  } // for i
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::produceComparisonCriteria() {
  buildOrbitalPopulations();
  calculateKineticEnergy();
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::compareOrbitals() {
  buildOrbitalMap();
  // Reduce the orbital map to a single vector for each system
  // which contains a list of the important orbitals.
  reduceOrbitalMap();
  // Clustering
  if (settings.noIsolatedMOs)
    orbitalClustering();
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::alignPiOrbitals() {
  for (unsigned int i = 1; i < _supersystems.size(); ++i) {
    LocalizationTask alignTask(_supersystems[i], {_supersystems[0]});
    alignTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN;
    alignTask.settings.useKineticAlign = settings.kineticAlign;
    alignTask.settings.alignExponent = settings.alignExponent;
    alignTask.settings.splitValenceAndCore = settings.splitValenceAndCore;
    alignTask.run();
  }
}

template class ActiveSpaceSelectionTask<Options::SCF_MODES::RESTRICTED>;
template class ActiveSpaceSelectionTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
