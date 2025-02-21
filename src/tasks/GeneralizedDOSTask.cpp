/**
 * @file GeneralizedDOSTask.cpp
 *
 * @author Moritz Bensberg
 * @date Sep 21, 2020
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
#include "tasks/GeneralizedDOSTask.h"
/* Include Serenity Internal Headers */
#include "analysis/directOrbitalSelection/DirectOrbitalSelection.h"
#include "analysis/populationAnalysis/IAOPopulationCalculator.h"
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h" //Atom to basis mapping
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h"
#include "misc/SerenityError.h"
#include "misc/SystemSplittingTools.h"
#include "settings/Settings.h" //Basis label
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <cmath>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GeneralizedDOSTask<SCFMode>::GeneralizedDOSTask(std::vector<std::shared_ptr<SystemController>> activeSystems,
                                                std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _supersystems(activeSystems),
    _allFragments(environmentSystems),
    _fragments(Matrix<std::shared_ptr<SystemController>>(0, 0, nullptr)) {
}

template<Options::SCF_MODES SCFMode>
void GeneralizedDOSTask<SCFMode>::run() {
  for (auto& sys : _supersystems) {
    if (sys->getSCFMode() != SCFMode)
      throw SerenityError("ERROR: The active system '" + sys->getSystemName() + "' needs the same SCFMode as the other systems.");
  }
  for (auto& sys : _allFragments) {
    if (sys->getSCFMode() != SCFMode)
      throw SerenityError("ERROR: The environment system '" + sys->getSystemName() +
                          "' needs the same SCFMode as the other systems.");
  }
  printSubSectionTitle("Generalized DOS");
  checkInput();
  _fragments = assignFragmentsToSupersystem(_allFragments, _supersystems.size(), _nFragments);

  OutputControl::vOut << "calculating orbital properties          ...";
  OutputControl::vOut.flush();
  // Prepare the DOS
  const auto populations = calculateOrbitalPopulations(_supersystems, settings.populationAlgorithm);
  const auto kinEnergies = calculateKineticEnergies(_supersystems, this->getMaxNumberOfOrbitalsFromPopulations(populations));
  OutputControl::vnOut << " done" << std::endl;
  OutputControl::vOut << "comparing orbitals                      ...";
  OutputControl::vOut.flush();
  const auto nOcc = _supersystems[0]->template getNOccupiedOrbitals<SCFMode>();
  DirectOrbitalSelection<SCFMode> dos(populations, kinEnergies, nOcc, settings.usePiBias, settings.biasThreshold,
                                      settings.biasAverage, settings.checkDegeneracies, settings.degeneracyFactor);
  OutputControl::vnOut << " done" << std::endl;

  if (this->settings.writeScores) {
    OutputControl::vOut << "Writing orbital scores to file   ...\n";
    writeScoresToFile(dos);
    OutputControl::vnOut << "                                        ... done" << std::endl;
    return;
  }

  // Run the DOS.
  _finalAssignments = std::make_unique<std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>>();
  *_finalAssignments =
      this->getIterativeAssignments(settings.similarityLocThreshold, settings.similarityKinEnergyThreshold, dos);
  OutputControl::vOut << "Splitting geometries into fragments     ...";
  OutputControl::vOut.flush();
  splitGeometryIntoFragments(*_finalAssignments);
  OutputControl::vnOut << " done" << std::endl;
  // Split the final electronic structure based on the final assignments
  splitElectronicStructureIntoFragments(*_finalAssignments);
  // Check if the result makes sense.
  if (not this->checkFinalOccupations())
    throw SerenityError("ERROR: GDOS detected an inconsistent orbital occupation for the given molecules.");
  if (settings.writeGroupsToFile)
    this->writeGroupsToFile();
}

template<Options::SCF_MODES SCFMode>
void GeneralizedDOSTask<SCFMode>::splitGeometryIntoFragments(const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments) {
  for (unsigned int iSuper = 0; iSuper < _supersystems.size(); ++iSuper) {
    auto geoms = SystemSplittingTools<SCFMode>::splitGeometry(_supersystems[iSuper], assignments[iSuper], settings.prioFirst,
                                                              settings.localizationThreshold, _fragments.cols());
    for (unsigned int iFrag = 0; iFrag < geoms.size(); ++iFrag) {
      auto fragment = _fragments(iSuper, iFrag);
      auto fragGeom = fragment->getGeometry();
      while (fragGeom->getNAtoms() > 0)
        fragGeom->deleteAtom(0);
      *fragGeom += *geoms[iFrag];
    } // for iFrag
  }   // for iSuper
}

template<Options::SCF_MODES SCFMode>
void GeneralizedDOSTask<SCFMode>::splitElectronicStructureIntoFragments(
    const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments) {
  for (unsigned int iSuper = 0; iSuper < _supersystems.size(); ++iSuper) {
    SystemSplittingTools<SCFMode>::splitSupersystemBasedOnAssignment(_supersystems[iSuper], _fragments.row(iSuper),
                                                                     assignments[iSuper]);
  } // for iSuper
}

template<Options::SCF_MODES SCFMode>
void GeneralizedDOSTask<SCFMode>::checkInput() {
  // Check atom ordering
  auto firstGeom = _supersystems[0]->getGeometry();
  auto firstAtomList = firstGeom->getAtoms();
  const auto& firstOccupations = _supersystems[0]->template getNOccupiedOrbitals<SCFMode>();
  for (unsigned int sys = 1; sys < _supersystems.size(); ++sys) {
    auto atomList = _supersystems[sys]->getGeometry()->getAtoms();
    if (atomList.size() != firstAtomList.size())
      throw SerenityError("The number of atoms have to be identical in all structures!");
    for (unsigned int iAtom = 0; iAtom < firstAtomList.size(); ++iAtom) {
      if (atomList[iAtom]->getAtomType()->getName() != firstAtomList[iAtom]->getAtomType()->getName())
        throw SerenityError("The atoms have to be ordered in the same way for all structures!");
    }
    // Check the number of electrons in the subsystems
    const auto& occupations = _supersystems[sys]->template getNOccupiedOrbitals<SCFMode>();
    for_spin(firstOccupations, occupations) {
      if (firstOccupations_spin != occupations_spin)
        throw SerenityError("Inconsistent number of electrons along the reaction coordinate!");
    }; // for_spin
  }
  auto supersystemBasisLabel = _supersystems[0]->getSettings().basis.label;
  for (unsigned int iSys = 1; iSys < _supersystems.size(); ++iSys) {
    if (_supersystems[iSys]->getSettings().basis.label != supersystemBasisLabel)
      throw SerenityError("ERROR: The basis-set label for all supersystems has to be identical!");
  }
  if (settings.similarityLocThreshold.size() != settings.similarityKinEnergyThreshold.size())
    throw SerenityError(
        (std::string) "ERROR: The number of provided similarityLocThreshold and similarityKinEnergyThreshold\n" +
        "       has to be identical.");
  if (settings.similarityLocThreshold.size() < 1)
    throw SerenityError(
        (std::string) "ERROR: At least one threshold has to be provided for any system partitioning!\n" +
        "       Please provide at least one similarityLocThreshold and similarityKinEnergyThreshold\n" + "       threshold.");
  for (unsigned int iThr = 1; iThr < settings.similarityLocThreshold.size(); ++iThr) {
    if ((settings.similarityLocThreshold[iThr] >= settings.similarityLocThreshold[iThr - 1]) ||
        (settings.similarityKinEnergyThreshold[iThr] >= settings.similarityKinEnergyThreshold[iThr - 1]))
      throw SerenityError((std::string) "ERROR: DOS thresholds have to be tightend/lowered in each step!\n" +
                          "       Please adjust the order/values of similarityKinEnergyThreshold\n" +
                          "       and/or similarityLocThreshold!");
  }
  _nFragments = settings.similarityLocThreshold.size() + 1;
  if (_nFragments * _supersystems.size() > _allFragments.size())
    throw SerenityError((std::string) "ERROR: GDOS expects " + std::to_string(_nFragments * _supersystems.size()) +
                        " environments systems, but only " + std::to_string(_allFragments.size()) + " were provided!");
}

template<Options::SCF_MODES SCFMode>
Matrix<std::shared_ptr<SystemController>>
GeneralizedDOSTask<SCFMode>::assignFragmentsToSupersystem(std::vector<std::shared_ptr<SystemController>> allFragments,
                                                          unsigned int nSupersystems, unsigned int nFragmentsEach) {
  Matrix<std::shared_ptr<SystemController>> fragments(nSupersystems, nFragmentsEach, nullptr);
  unsigned int counter = 0;
  OutputControl::nOut << "--------------------------------------------------------" << std::endl;
  OutputControl::nOut << "  The supersystem to fragment partitioning:" << std::endl;
  for (unsigned int iSuper = 0; iSuper < nSupersystems; ++iSuper) {
    std::string name = _supersystems[iSuper]->getSystemName();
    OutputControl::nOut << _supersystems[iSuper]->getSystemName() << "       ";
    for (unsigned int iFrag = 0; iFrag < nFragmentsEach; ++iFrag) {
      fragments(iSuper, iFrag) = allFragments[counter];
      OutputControl::nOut << allFragments[counter]->getSystemName() << "  ";
      ++counter;
    } // for iFrag
    OutputControl::nOut << std::endl;
  } // for iSuper
  OutputControl::nOut << "--------------------------------------------------------" << std::endl;
  return fragments;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, unsigned int>
GeneralizedDOSTask<SCFMode>::getMaxNumberOfOrbitalsFromPopulations(const std::vector<SPMatrix<SCFMode>>& populations) {
  SpinPolarizedData<SCFMode, unsigned int> nOrbitals(0);
  const auto& populationZero = populations[0];
  for_spin(populationZero, nOrbitals) {
    nOrbitals_spin = populationZero_spin.cols();
  };
  return nOrbitals;
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>>
GeneralizedDOSTask<SCFMode>::calculateKineticEnergies(std::vector<std::shared_ptr<SystemController>> supersystems,
                                                      SpinPolarizedData<SCFMode, unsigned int> nOrbitals) {
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>> kineticEnergies = {};
  for (const auto& system : supersystems) {
    // Calculate the kinetic energies of the occupied orbitals
    const auto& coeff = system->getActiveOrbitalController<SCFMode>()->getCoefficients();
    auto kinIntegrals = system->getOneElectronIntegralController()->getKinIntegrals();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> kin;
    for_spin(nOrbitals, coeff, kin) {
      kin_spin.resize(nOrbitals_spin);
      for (unsigned int occI = 0; occI < nOrbitals_spin; ++occI) {
        const auto& coeffI = coeff_spin.col(occI);
        kin_spin[occI] = coeffI.transpose() * kinIntegrals * coeffI;
      } // for occI
    };
    kineticEnergies.push_back(kin);
  } // for i
  return kineticEnergies;
}

template<Options::SCF_MODES SCFMode>
std::vector<SPMatrix<SCFMode>>
GeneralizedDOSTask<SCFMode>::calculateOrbitalPopulations(std::vector<std::shared_ptr<SystemController>> supersystems,
                                                         Options::POPULATION_ANALYSIS_ALGORITHMS alg) {
  std::vector<SPMatrix<SCFMode>> populations = {};
  for (auto& system : supersystems) {
    SPMatrix<SCFMode> orbitalPopulations;
    if (alg == Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN) {
      orbitalPopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
          system->getActiveOrbitalController<SCFMode>()->getCoefficients(),
          system->getOneElectronIntegralController()->getOverlapIntegrals(),
          system->getAtomCenteredBasisController()->getBasisIndices());
    }
    else if (alg == Options::POPULATION_ANALYSIS_ALGORITHMS::IAO) {
      orbitalPopulations = IAOPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(system, settings.mapVirtuals);
    }
    else if (alg == Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell) {
      // This is basically the square of the CIAO-coefficient summed over a shell.
      orbitalPopulations =
          IAOPopulationCalculator<SCFMode>::calculateShellwiseOrbitalPopulations(system, settings.mapVirtuals);
    }
    else {
      throw SerenityError("The algorithm used for evaluating orbital-wise populations is not supported.");
    }
    SPMatrix<SCFMode> occPopulations;
    SpinPolarizedData<SCFMode, unsigned int> nOrbitals(system->getBasisController()->getNBasisFunctions());
    if (!settings.mapVirtuals) {
      nOrbitals = system->getNOccupiedOrbitals<SCFMode>();
    }
    for_spin(occPopulations, orbitalPopulations, nOrbitals) {
      const unsigned int nMaxOrbitals = std::min(nOrbitals_spin, (unsigned int)orbitalPopulations_spin.cols());
      occPopulations_spin = orbitalPopulations_spin.leftCols(nMaxOrbitals);
    };
    populations.push_back(occPopulations);
  } /* for sys : _supersystems */
  return populations;
}

template<Options::SCF_MODES SCFMode>
void GeneralizedDOSTask<SCFMode>::setAssignIndices(std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& finalAssignment,
                                                   const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& newAssignment,
                                                   int assignIndex) {
  unsigned int nSystems = finalAssignment.size();
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    const auto& newAssign = newAssignment[iSys];
    auto& finalAssign = finalAssignment[iSys];
    for_spin(newAssign, finalAssign) {
      for (unsigned int iOcc = 0; iOcc < newAssign_spin.size(); ++iOcc) {
        bool wasAssignedBefore = finalAssign_spin[iOcc] < assignIndex;
        bool isAssignedNow = newAssign_spin[iOcc] == 0;
        // If the orbital was already assigned, skip it.
        if (wasAssignedBefore)
          continue;
        // If it was considered as changed in this iteration (and no iteration before), select it.
        if (isAssignedNow) {
          finalAssign_spin[iOcc] = assignIndex;
        }
        else {
          // Else: Set the orbital-index to the next higher fragment-index.
          finalAssign_spin[iOcc] = assignIndex + 1;
        }
      } // for iOcc
    };
  } // for iSys
}
template<Options::SCF_MODES SCFMode>
bool GeneralizedDOSTask<SCFMode>::checkFinalOccupations() {
  for (unsigned int iSlec = 0; iSlec < _nFragments; ++iSlec) {
    auto fragmentSet = _fragments.col(iSlec);
    auto refOcc = fragmentSet[0]->template getNOccupiedOrbitals<SCFMode>();
    for (const auto& frag : fragmentSet) {
      auto fragOcc = frag->template getNOccupiedOrbitals<SCFMode>();
      // The for_spin below is a lambda function, thus I will avoid calling the
      // return false for inconsistent occupations from within.
      bool correct = true;
      for_spin(refOcc, fragOcc) {
        if (refOcc_spin != fragOcc_spin)
          correct = false;
      };
      if (not correct)
        return false;
    } // for frag
  }   // for iSlec
  return true;
}

template<Options::SCF_MODES SCFMode>
GeneralizedDOSTask<SCFMode>::~GeneralizedDOSTask() = default;

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
GeneralizedDOSTask<SCFMode>::getIterativeAssignments(std::vector<double> locs, std::vector<double> kins,
                                                     DirectOrbitalSelection<SCFMode>& dos) {
  OutputControl::vOut << "running iterative selection             ...";
  OutputControl::vOut.flush();
  // Run the initial DOS. We will use this object to collect all future indices.
  // In the first iteration changed orbitals will have the index 0 and all other the index 1.
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> lastMap;
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> finalAssignments;
  double threshold = locs[0];
  double virtThreshold = locs[0];
  if (settings.bestMatchMapping && !settings.writeScores) {
    finalAssignments = dos.getBestMatchAssignments(threshold, virtThreshold, lastMap,
                                                   _supersystems[0]->template getNOccupiedOrbitals<SCFMode>(),
                                                   settings.scoreEnd, settings.scoreStart);
  }
  else {
    finalAssignments = dos.getOrbitalAssignments(locs[0], kins[0], lastMap);
  }

  // Keep track of the orbital groups during the assignment. The DOS algorithm will try to map
  // orbitals between the systems, and we want to bring this mapping into a better representation.
  // The DOS algorithm will not be able to map all orbitals between the systems. These initially
  // selected orbitals constitute their own, initial orbital group.
  const auto& nOcc = _supersystems[0]->template getNOccupiedOrbitals<SCFMode>();
  auto initialGroups = DOSOrbitalGroup::orbitalGroupFromAssignment<SCFMode>(finalAssignments, 0);
  for_spin(initialGroups, _orbitalGroups, nOcc) {
    auto virtOrbitals = initialGroups_spin.splitOffVirtualOrbitals(nOcc_spin);
    _orbitalGroups_spin.push_back(std::make_shared<DOSOrbitalGroup>(initialGroups_spin));
    if (virtOrbitals.getNOrbitals())
      _orbitalGroups_spin.push_back(std::make_shared<DOSOrbitalGroup>(virtOrbitals));
  };
  _unmappableOrbitalGroups = _orbitalGroups;
  // Rerun the DOS with increasingly tighter thresholds.
  for (unsigned int iAssign = 1; iAssign < locs.size(); ++iAssign) {
    // Run the next assignment with the new thresholds.
    std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> currentMap;
    const auto newAssignments = dos.getOrbitalAssignments(locs[iAssign], kins[iAssign], currentMap);
    // Set the fragment indices for the orbitals assigned in this iteration to iAssign.
    setAssignIndices(finalAssignments, newAssignments, iAssign);
    // Keep track of the orbital groups assigned just now.
    auto newGroups = DOSOrbitalGroup::orbitalGroupsFromMapWithAssignment<SCFMode>(lastMap, finalAssignments[0], iAssign);
    for_spin(newGroups, _orbitalGroups) {
      _orbitalGroups_spin.insert(_orbitalGroups_spin.end(), newGroups_spin.begin(), newGroups_spin.end());
    };
    lastMap = currentMap;
  } // for iAssign
  auto lastGroups = DOSOrbitalGroup::orbitalGroupsFromMapWithAssignment<SCFMode>(lastMap, finalAssignments[0], locs.size());
  for_spin(lastGroups, _orbitalGroups) {
    _orbitalGroups_spin.insert(_orbitalGroups_spin.end(), lastGroups_spin.begin(), lastGroups_spin.end());
  };
  OutputControl::vnOut << " done" << std::endl;
  if (settings.bestMatchMapping && !settings.writeScores)
    OutputControl::nOut << "  Optimized (occupied) threshold " << threshold << "  (virtual) threshold " << virtThreshold
                        << std::endl;
  return finalAssignments;
}

template<Options::SCF_MODES SCFMode>
void GeneralizedDOSTask<SCFMode>::writeScoresToFile(DirectOrbitalSelection<SCFMode>& dos) {
  if (this->settings.scoreStart < this->settings.scoreEnd)
    throw SerenityError("ERROR: The setting <scoreStart> has to be larger than <scoreEnd>.");
  double factor = pow(settings.scoreEnd / settings.scoreStart, 1 / (double)settings.nTest);
  std::vector<double> thresholds = {};
  double threshold = this->settings.scoreStart;
  for (unsigned int i = 0; i < this->settings.nTest; ++i) {
    thresholds.push_back(threshold);
    threshold *= factor;
  }
  OutputControl::nOut << " inversion start:       " << settings.scoreStart << std::endl;
  OutputControl::nOut << " inversion end:         " << settings.scoreEnd << std::endl;
  OutputControl::nOut << " inversion step-size:   " << factor << std::endl;
  OutputControl::nOut << " n-Steps:               " << settings.nTest << std::endl;
  unsigned int counter = 0;
  OutputControl::dOut << " Selection thresholds:  " << std::endl;
  for (const auto& tr : thresholds) {
    OutputControl::dOut << tr << " ";
    ++counter;
    if (counter % 10 == 0)
      OutputControl::dOut << std::endl;
  }
  OutputControl::dOut << std::endl;
  // Run selection
  auto assignments = getIterativeAssignments(thresholds, thresholds, dos);
  // We will have nTest+1 different assignments! The score for the last system will simply be the scoreEnd.
  thresholds.push_back(settings.scoreEnd);
  // Initialize
  SPMatrix<SCFMode> systemWiseScores;
  auto zeroAssign = assignments[0];
  for_spin(zeroAssign, systemWiseScores) {
    unsigned int nOcc = zeroAssign_spin.size();
    systemWiseScores_spin = Eigen::MatrixXd::Zero(nOcc, assignments.size());
  };
  // Calculate scores
  for (unsigned int iSys = 0; iSys < assignments.size(); ++iSys) {
    auto assign = assignments[iSys];
    for_spin(assign, systemWiseScores) {
      Eigen::VectorXd scores = Eigen::VectorXd::Zero(assign_spin.size());
      for (unsigned int iOrb = 0; iOrb < assign_spin.size(); ++iOrb) {
        double score = thresholds[assign_spin[iOrb]];
        scores[iOrb] = score;
      }
      systemWiseScores_spin.col(iSys) = scores;
    };
  } // for iSys
  // Write to file
  writeFile(systemWiseScores);
}

template<>
void GeneralizedDOSTask<RESTRICTED>::writeFile(SPMatrix<RESTRICTED> scores) {
  this->writeMatrix(scores, "dosScores.dat");
}

template<>
void GeneralizedDOSTask<UNRESTRICTED>::writeFile(SPMatrix<UNRESTRICTED> scores) {
  this->writeMatrix(scores.alpha, "dosScores.alpha.dat");
  this->writeMatrix(scores.beta, "dosScores.beta.dat");
}

template<Options::SCF_MODES SCFMode>
void GeneralizedDOSTask<SCFMode>::writeMatrix(Eigen::MatrixXd scores, std::string fileName) {
  std::ofstream ofs;
  ofs.open(fileName, std::ofstream::out | std::ofstream::trunc);
  for (const auto& sys : _supersystems)
    ofs << sys->getSystemName() << " ";
  ofs << std::endl;
  ofs << scores << std::endl;
  ofs.close();
}

template<Options::SCF_MODES SCFMode>
const SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>>&
GeneralizedDOSTask<SCFMode>::getOrbitalGroups() {
  bool isInitialized = true;
  for_spin(_orbitalGroups) {
    isInitialized = (_orbitalGroups_spin.size() > 0 && isInitialized) ? true : false;
  };
  if (not isInitialized)
    this->run();
  return _orbitalGroups;
}

template<Options::SCF_MODES SCFMode>
const SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>>&
GeneralizedDOSTask<SCFMode>::getUnmappableOrbitalGroups() {
  this->getOrbitalGroups();
  return _unmappableOrbitalGroups;
}

template<Options::SCF_MODES SCFMode>
void GeneralizedDOSTask<SCFMode>::writeGroupSetToFile(std::string fileName,
                                                      std::vector<std::shared_ptr<DOSOrbitalGroup>> groups) {
  std::ofstream ofs;
  ofs.open(fileName, std::ofstream::out | std::ofstream::trunc);
  ofs << "Columns: System index, rows: orbitals in set as unordered list." << std::endl;
  for (const auto& group : groups) {
    if (group->getNOrbitals())
      ofs << *group << std::endl;
  }
  ofs.close();
}

template<>
void GeneralizedDOSTask<RESTRICTED>::writeGroupsToFile() {
  const auto& orbitalGroups = this->getOrbitalGroups();
  this->writeGroupSetToFile("orbital_groups.dat", orbitalGroups);
}

template<>
void GeneralizedDOSTask<UNRESTRICTED>::writeGroupsToFile() {
  const auto& orbitalGroups = this->getOrbitalGroups();
  this->writeGroupSetToFile("orbital_groups.alpha.dat", orbitalGroups.alpha);
  this->writeGroupSetToFile("orbital_groups.beta.dat", orbitalGroups.beta);
}

template<Options::SCF_MODES SCFMode>
const std::vector<SpinPolarizedData<SCFMode, Eigen::SparseMatrix<int>>>&
GeneralizedDOSTask<SCFMode>::getSuperToSubsystemOccSortingMatrices() {
  if (!_superToSubsystemOccSortingMatrices) {
    if (!_finalAssignments)
      this->run();
    _superToSubsystemOccSortingMatrices =
        std::make_unique<std::vector<SpinPolarizedData<SCFMode, Eigen::SparseMatrix<int>>>>();
    const unsigned int nSuper = _finalAssignments->size();

    for (unsigned int iSuper = 0; iSuper < nSuper; ++iSuper) {
      _superToSubsystemOccSortingMatrices->push_back(SpinPolarizedData<SCFMode, Eigen::SparseMatrix<int>>(0, 0));
      const auto& assignment = (*_finalAssignments)[iSuper];
      auto& sortingMatrix = (*_superToSubsystemOccSortingMatrices)[iSuper];
      for_spin(assignment, sortingMatrix) {
        std::vector<unsigned int> subsystemWiseCounter = {0};
        for (unsigned int iSub = 1; iSub < _nFragments; ++iSub) {
          subsystemWiseCounter.push_back((assignment_spin.array() < iSub).count());
        }
        const unsigned int nOcc = assignment_spin.size();
        std::vector<Eigen::Triplet<int>> triplets;
        for (unsigned int iOcc = 0; iOcc < nOcc; ++iOcc) {
          auto& counter = subsystemWiseCounter[assignment_spin(iOcc)];
          triplets.push_back(Eigen::Triplet<int>(counter, iOcc, 1));
          counter++;
        } // for iOcc
        Eigen::SparseMatrix<int> newSortingMatrix(nOcc, nOcc);
        newSortingMatrix.setFromTriplets(triplets.begin(), triplets.end());
        sortingMatrix_spin = newSortingMatrix;
      };
    } // for iSuper
  }
  return *_superToSubsystemOccSortingMatrices;
}

template class GeneralizedDOSTask<Options::SCF_MODES::RESTRICTED>;
template class GeneralizedDOSTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
