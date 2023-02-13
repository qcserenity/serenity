/**
 * @file DirectOrbitalSelection.cpp
 *
 * @author Moritz Bensberg
 * @date Sep 18, 2020
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
#include "analysis/directOrbitalSelection/DirectOrbitalSelection.h"
/* Include Serenity Internal Headers */
#include "io/FormattedOutputStream.h"
#include "misc/SerenityError.h" //Errors.

namespace Serenity {

DOSOrbitalGroup::DOSOrbitalGroup(std::vector<std::vector<unsigned int>> systemWiseOrbitalIndices)
  : _systemWiseOrbitalIndices(systemWiseOrbitalIndices) {
  for (const auto& indices : _systemWiseOrbitalIndices) {
    if (indices.size() != _systemWiseOrbitalIndices[0].size())
      throw SerenityError(
          "Implementation error! The orbital group has to contain the same number of orbitals for each system!");
  }
}

const std::vector<unsigned int>& DOSOrbitalGroup::getOrbitalIndicesForSystem(unsigned int systemIndex) {
  if (systemIndex >= _systemWiseOrbitalIndices.size())
    throw SerenityError("Implementation error! Orbital-wise indices for the given DOS group are either not initialized "
                        "or the system index is incorrect/to large.");
  return _systemWiseOrbitalIndices[systemIndex];
}

void DOSOrbitalGroup::addSystemWiseIndices(std::vector<unsigned int> systemIndices) {
  if (_systemWiseOrbitalIndices.size() > 0) {
    if (_systemWiseOrbitalIndices[0].size() != systemIndices.size())
      throw SerenityError(
          "Implementation error! The orbital group has to contain the same number of orbitals for each system!");
  }
  _systemWiseOrbitalIndices.push_back(systemIndices);
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, DOSOrbitalGroup>
DOSOrbitalGroup::orbitalGroupFromAssignment(const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignment,
                                            unsigned int targetIndex) {
  SpinPolarizedData<SCFMode, DOSOrbitalGroup> group;
  for (const auto& assignmentSet : assignment) {
    for_spin(assignmentSet, group) {
      std::vector<unsigned int> systemWiseOrbitalIndices;
      for (unsigned int iOrb = 0; iOrb < assignmentSet_spin.size(); ++iOrb) {
        if (assignmentSet_spin(iOrb) == (int)targetIndex)
          systemWiseOrbitalIndices.push_back(iOrb);
      }
      group_spin.addSystemWiseIndices(systemWiseOrbitalIndices);
    };
  }
  return group;
}

DOSOrbitalGroup DOSOrbitalGroup::splitOffVirtualOrbitals(const unsigned int nOcc) {
  std::vector<std::vector<unsigned int>> systemWiseVirtualOrbitalIndices(this->getNSystems());
  std::vector<std::vector<unsigned int>> systemWiseOccupiedOrbitalIndices(this->getNSystems());
  for (unsigned int iSys = 0; iSys < this->getNSystems(); ++iSys) {
    const auto& systemOrbitalIndices = _systemWiseOrbitalIndices[iSys];
    for (const auto& iOrb : systemOrbitalIndices) {
      if (iOrb < nOcc) {
        systemWiseOccupiedOrbitalIndices[iSys].push_back(iOrb);
      }
      else {
        systemWiseVirtualOrbitalIndices[iSys].push_back(iOrb);
      }
    } // for iOrb
  }   // for systemOrbitalIndices
  _systemWiseOrbitalIndices = systemWiseOccupiedOrbitalIndices;
  return DOSOrbitalGroup(systemWiseVirtualOrbitalIndices);
}

void DOSOrbitalGroup::updateOrbitalIndices(std::vector<Eigen::SparseMatrix<int>> newXorigMap) {
  for (unsigned int iSys = 0; iSys < this->getNSystems(); ++iSys) {
    auto& systemIndexSet = _systemWiseOrbitalIndices[iSys];
    for (auto& iOrb : systemIndexSet) {
      unsigned int counter = 0;
      // This loop will only ever have one iteration!.
      for (Eigen::SparseMatrix<int>::InnerIterator itOrb(newXorigMap[iSys], iOrb); itOrb; ++itOrb) {
        iOrb = itOrb.row();
        ++counter;
      } // for itOrb
      if (counter > 1)
        throw SerenityError(
            "Implementation ERROR! The new to original orbital map has more than one assignment per column!");
    } // for iOrb
  }   // for iSys
}

std::vector<Eigen::SparseMatrix<int>>
DOSOrbitalGroup::buildOrbitalToOrbitalGroupSortingMatrix(std::vector<std::shared_ptr<DOSOrbitalGroup>> orbitalGroups,
                                                         const unsigned int nOcc) {
  if (orbitalGroups.size() < 1)
    throw SerenityError(
        "Implementation Error! Orbital groups were not initialized correctly before sorting matrix construction.");
  const unsigned int nSystems = orbitalGroups[0]->getNSystems();
  const unsigned int nGroups = orbitalGroups.size();
  std::vector<Eigen::SparseMatrix<int>> maps(nSystems, Eigen::SparseMatrix<int>(nOcc, nGroups));
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    std::vector<Eigen::Triplet<int>> triplets;
    for (unsigned int iGroup = 0; iGroup < nGroups; ++iGroup) {
      const auto& group = orbitalGroups[iGroup];
      const std::vector<unsigned int> orbitalIndices = group->getOrbitalIndicesForSystem(iSys);
      for (const auto& index : orbitalIndices) {
        triplets.push_back(Eigen::Triplet<int>(index, iGroup, 1));
      } // for index
    }   // for iGroup
    maps[iSys].setFromTriplets(triplets.begin(), triplets.end());
  } // for iSys
  return maps;
}

template SpinPolarizedData<RESTRICTED, DOSOrbitalGroup> DOSOrbitalGroup::orbitalGroupFromAssignment<RESTRICTED>(
    const std::vector<SpinPolarizedData<RESTRICTED, Eigen::VectorXi>>& assignment, unsigned int targetIndex);
template SpinPolarizedData<UNRESTRICTED, DOSOrbitalGroup> DOSOrbitalGroup::orbitalGroupFromAssignment<UNRESTRICTED>(
    const std::vector<SpinPolarizedData<UNRESTRICTED, Eigen::VectorXi>>& assignment, unsigned int targetIndex);

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>> DOSOrbitalGroup::orbitalGroupsFromMapWithAssignment(
    const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& map,
    const SpinPolarizedData<SCFMode, Eigen::VectorXi>& assignment, unsigned int targetIndex) {
  const unsigned int nSystems = map.size();
  SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>> groups;
  // The order of std::vector and SpinPolarized data in map is not practical here. Hence, I will reorder the map.
  SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXi>> reorderedMap;
  for_spin(reorderedMap) {
    reorderedMap_spin = std::vector<Eigen::MatrixXi>(nSystems, Eigen::MatrixXi(0, 0));
  };
  // We can not directly use the map 0/0 because of the way it is handled/largely ignored in the orignial DOS.
  // Thus, we use the chain for the maps 0->1->0 instead of the 0/0 map.
  const auto& map10 = map[1][0];
  for_spin(reorderedMap, map10) {
    reorderedMap_spin[0] = (map10_spin.transpose() * map10_spin).eval();
  };

  for (unsigned int iSys = 1; iSys < nSystems; ++iSys) {
    const auto& mapI0 = map[iSys][0];
    for_spin(mapI0, reorderedMap) {
      reorderedMap_spin[iSys] = mapI0_spin;
    };
  } // for iSys

  for_spin(assignment, groups, reorderedMap) {
    const unsigned int nOrbs = assignment_spin.size();
    std::vector<bool> alreadyAssigned(nOrbs, false);
    for (unsigned int iOrb = 0; iOrb < nOrbs; ++iOrb) {
      if (assignment_spin(iOrb) == (int)targetIndex && not alreadyAssigned[iOrb]) {
        std::vector<std::vector<unsigned int>> systemWiseOrbitalIndices;
        for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
          const auto& mapI0 = reorderedMap_spin[iSys];
          std::vector<unsigned int> orbitalIndices_iSys;
          // Get the orbital it is mapped to in each system.
          for (unsigned int jOrb = 0; jOrb < nOrbs; ++jOrb) {
            if (mapI0(jOrb, iOrb) > 0) {
              orbitalIndices_iSys.push_back(jOrb);
              if (iSys == 0) {
                alreadyAssigned[jOrb] = true;
              } // if iSys == 0
            }   // if mapped
          }     // for jOrb
          systemWiseOrbitalIndices.push_back(orbitalIndices_iSys);
        } // for iSys
        groups_spin.push_back(std::make_shared<DOSOrbitalGroup>(systemWiseOrbitalIndices));
      } // if assignment_spin(iOrb) == targetIndex
    }   // for iOrb
  };
  return groups;
}
template SpinPolarizedData<RESTRICTED, std::vector<std::shared_ptr<DOSOrbitalGroup>>>
DOSOrbitalGroup::orbitalGroupsFromMapWithAssignment<RESTRICTED>(
    const std::vector<std::vector<SpinPolarizedData<RESTRICTED, Eigen::MatrixXi>>>& map,
    const SpinPolarizedData<RESTRICTED, Eigen::VectorXi>& assignment, unsigned int targetIndex);

template SpinPolarizedData<UNRESTRICTED, std::vector<std::shared_ptr<DOSOrbitalGroup>>>
DOSOrbitalGroup::orbitalGroupsFromMapWithAssignment<UNRESTRICTED>(
    const std::vector<std::vector<SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXi>>>& map,
    const SpinPolarizedData<UNRESTRICTED, Eigen::VectorXi>& assignment, unsigned int targetIndex);

Eigen::VectorXd DOSOrbitalGroup::getGroupEnergiesFromPairEnergies(const std::vector<Eigen::MatrixXd>& pairEnergies) const {
  if (pairEnergies.size() != _systemWiseOrbitalIndices.size())
    throw SerenityError("Implementation error! The number of systems for which pair energies are available and the "
                        "number of systems for which orbital indices were saved have to identical.");
  const unsigned int nSystems = _systemWiseOrbitalIndices.size();
  Eigen::VectorXd groupEnergies = Eigen::VectorXd::Zero(nSystems);
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    for (const auto& iOrb : _systemWiseOrbitalIndices[iSys]) {
      groupEnergies(iSys) += pairEnergies[iSys].col(iOrb).sum();
    } // for iOrb
  }   // for iSys
  return groupEnergies;
}

Eigen::MatrixXd DOSOrbitalGroup::getGroupEnergyDifferenceMatrix(const std::vector<Eigen::MatrixXd>& pairEnergies) const {
  const Eigen::VectorXd groupEnergies = getGroupEnergiesFromPairEnergies(pairEnergies);
  const unsigned int nSystems = _systemWiseOrbitalIndices.size();
  Eigen::MatrixXd differenceMatrix = Eigen::MatrixXd::Zero(nSystems, nSystems);
  differenceMatrix.rowwise() += groupEnergies.transpose();
  differenceMatrix.colwise() -= groupEnergies;
  return differenceMatrix;
}

unsigned int DOSOrbitalGroup::getNOrbitals() const {
  return _systemWiseOrbitalIndices[0].size();
}

unsigned int DOSOrbitalGroup::getNSystems() const {
  return _systemWiseOrbitalIndices.size();
}

const std::vector<std::vector<unsigned int>>& DOSOrbitalGroup::getSystemWiseIndices() const {
  return _systemWiseOrbitalIndices;
}

std::vector<std::vector<std::vector<unsigned int>>>
DOSOrbitalGroup::getIndicesFromGroups(std::vector<std::shared_ptr<DOSOrbitalGroup>> groups) {
  std::vector<std::vector<std::vector<unsigned int>>> allIndices;
  for (const auto& group : groups) {
    allIndices.push_back(group->getSystemWiseIndices());
  }
  return allIndices;
}

template<Options::SCF_MODES SCFMode>
DirectOrbitalSelection<SCFMode>::DirectOrbitalSelection(const std::vector<SPMatrix<SCFMode>>& orbitalPopulations,
                                                        const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>> kineticEnergies,
                                                        const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals,
                                                        bool usePiBias, double biasThreshold, double biasAverage,
                                                        bool checkDegeneracies, double degeneracyFactor,
                                                        bool excludeOccVirtMappings)
  : _orbitalPopulations(orbitalPopulations),
    _kineticEnergies(kineticEnergies),
    _usePiBias(usePiBias),
    _biasThreshold(biasThreshold),
    _biasAverage(biasAverage),
    _checkDegeneracies(checkDegeneracies),
    _degeneracyFactor(degeneracyFactor),
    _excludeOccVirtMappings(excludeOccVirtMappings),
    _nOrbitals(initializeNOrbitals()) {
  if (orbitalPopulations.size() != kineticEnergies.size() || orbitalPopulations.size() < 1)
    throw SerenityError("ERROR: Insufficient information provided for the direct orbital selection. Not enough orbital "
                        "populations or kinetic energies.");
  // Prepare comparison by precalculating all differences.
  calculateLocScore(nOccupiedOrbitals);
  calculateKinScore(nOccupiedOrbitals);
}
template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, unsigned int> DirectOrbitalSelection<SCFMode>::initializeNOrbitals() {
  SpinPolarizedData<SCFMode, unsigned int> nOrbitals(0);
  const auto& orbitalPopulationsZero = _orbitalPopulations[0];
  for_spin(orbitalPopulationsZero, nOrbitals) {
    nOrbitals_spin = orbitalPopulationsZero_spin.cols();
  };
  return nOrbitals;
}

template<Options::SCF_MODES SCFMode>
DirectOrbitalSelection<SCFMode>::~DirectOrbitalSelection() = default;

template<Options::SCF_MODES SCFMode>
void DirectOrbitalSelection<SCFMode>::calculateLocScore(const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals) {
  initializeScores(_locScore, _nOrbitals);
  unsigned int nSystems = _orbitalPopulations.size();
  for (unsigned int i = 0; i < nSystems; ++i) {
    const auto& popsI = _orbitalPopulations[i];
    for (unsigned int j = i; j < nSystems; ++j) {
      const auto& popsJ = _orbitalPopulations[j];
      auto& scoreIJ = _locScore[i][j];
      auto& scoreJI = _locScore[j][i];
      for_spin(popsI, popsJ, scoreIJ, scoreJI, nOccupiedOrbitals) {
        const unsigned int nOrbsI = scoreIJ_spin.rows();
        const unsigned int nOrbsJ = scoreIJ_spin.cols();
        // occ--occ block
        for (unsigned int k = 0; k < nOccupiedOrbitals_spin; ++k) {
          for (unsigned int l = 0; l < nOccupiedOrbitals_spin; ++l) {
            double populationMeasure = (popsI_spin.col(k) - popsJ_spin.col(l)).array().abs().sum();
            scoreIJ_spin(k, l) = populationMeasure;
            scoreJI_spin(l, k) = populationMeasure;
          } // for l
        }   // for k
        // virt--virt block
        for (unsigned int k = nOccupiedOrbitals_spin; k < nOrbsI; ++k) {
          for (unsigned int l = nOccupiedOrbitals_spin; l < nOrbsJ; ++l) {
            double populationMeasure = (popsI_spin.col(k) - popsJ_spin.col(l)).array().abs().sum();
            scoreIJ_spin(k, l) = populationMeasure;
            scoreJI_spin(l, k) = populationMeasure;
          } // for l
        }   // for k
        if (!_excludeOccVirtMappings) {
          // occ--virt blocks
          for (unsigned int k = nOccupiedOrbitals_spin; k < nOrbsI; ++k) {
            for (unsigned int l = 0; l < nOccupiedOrbitals_spin; ++l) {
              double populationMeasure = (popsI_spin.col(k) - popsJ_spin.col(l)).array().abs().sum();
              scoreIJ_spin(k, l) = populationMeasure;
              scoreJI_spin(l, k) = populationMeasure;
            } // for l
          }   // for k
          for (unsigned int k = 0; k < nOccupiedOrbitals_spin; ++k) {
            for (unsigned int l = nOccupiedOrbitals_spin; l < nOrbsJ; ++l) {
              double populationMeasure = (popsI_spin.col(k) - popsJ_spin.col(l)).array().abs().sum();
              scoreIJ_spin(k, l) = populationMeasure;
              scoreJI_spin(l, k) = populationMeasure;
            } // for l
          }   // for k
        }
      };
    } // for j
  }   // for i
}

template<Options::SCF_MODES SCFMode>
void DirectOrbitalSelection<SCFMode>::calculateKinScore(const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals) {
  initializeScores(_kinScore, _nOrbitals);
  unsigned int nSystems = _orbitalPopulations.size();
  for (unsigned int i = 0; i < nSystems; ++i) {
    const auto& kinsI = _kineticEnergies[i];
    for (unsigned int j = i; j < nSystems; ++j) {
      const auto& kinsJ = _kineticEnergies[j];
      auto& scoreIJ = _kinScore[i][j];
      auto& scoreJI = _kinScore[j][i];
      for_spin(kinsI, kinsJ, scoreIJ, scoreJI, nOccupiedOrbitals) {
        const unsigned int nOrbsI = scoreIJ_spin.rows();
        const unsigned int nOrbsJ = scoreIJ_spin.cols();
        // occ--occ block
        for (unsigned int k = 0; k < nOccupiedOrbitals_spin; ++k) {
          for (unsigned int l = 0; l < nOccupiedOrbitals_spin; ++l) {
            double kineticEnergyMeasure = std::fabs(kinsI_spin[k] - kinsJ_spin[l]);
            scoreIJ_spin(k, l) = kineticEnergyMeasure;
            scoreJI_spin(l, k) = kineticEnergyMeasure;
          } // for l
        }   // for k
        // virt--virt block
        for (unsigned int k = nOccupiedOrbitals_spin; k < nOrbsI; ++k) {
          for (unsigned int l = nOccupiedOrbitals_spin; l < nOrbsJ; ++l) {
            double kineticEnergyMeasure = std::fabs(kinsI_spin[k] - kinsJ_spin[l]);
            scoreIJ_spin(k, l) = kineticEnergyMeasure;
            scoreJI_spin(l, k) = kineticEnergyMeasure;
          } // for l
        }   // for k
        if (!this->_excludeOccVirtMappings) {
          // virt--occ blocks
          for (unsigned int k = nOccupiedOrbitals_spin; k < nOrbsI; ++k) {
            for (unsigned int l = 0; l < nOccupiedOrbitals_spin; ++l) {
              double kineticEnergyMeasure = std::fabs(kinsI_spin[k] - kinsJ_spin[l]);
              scoreIJ_spin(k, l) = kineticEnergyMeasure;
              scoreJI_spin(l, k) = kineticEnergyMeasure;
            } // for l
          }   // for k
          for (unsigned int k = 0; k < nOccupiedOrbitals_spin; ++k) {
            for (unsigned int l = nOccupiedOrbitals_spin; l < nOrbsJ; ++l) {
              double kineticEnergyMeasure = std::fabs(kinsI_spin[k] - kinsJ_spin[l]);
              scoreIJ_spin(k, l) = kineticEnergyMeasure;
              scoreJI_spin(l, k) = kineticEnergyMeasure;
            } // for l
          }   // for k
        }
      };
    } // for j
  }   // for i
}

template<Options::SCF_MODES SCFMode>
void DirectOrbitalSelection<SCFMode>::initializeScores(std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>>& scores,
                                                       const SpinPolarizedData<SCFMode, unsigned int>& nOrbitals) {
  unsigned int nSystems = _orbitalPopulations.size();
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> scoreVector;
    for (unsigned int jSys = 0; jSys < nSystems; ++jSys) {
      SpinPolarizedData<SCFMode, Eigen::MatrixXd> newScores;
      for_spin(newScores, nOrbitals) {
        newScores_spin = Eigen::MatrixXd::Constant(nOrbitals_spin, nOrbitals_spin, std::numeric_limits<double>::infinity());
      };
      scoreVector.push_back(newScores);
    } // for jSys
    scores.push_back(scoreVector);
  } // for iSys
}

template<Options::SCF_MODES SCFMode>
std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> DirectOrbitalSelection<SCFMode>::initializeOrbitalMap() {
  unsigned int nSystems = _orbitalPopulations.size();
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> maps;
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>> mapVector;
    for (unsigned int jSys = 0; jSys < nSystems; ++jSys) {
      SpinPolarizedData<SCFMode, Eigen::MatrixXi> newMap;
      const auto& locScoreIJ = _locScore[iSys][jSys];
      for_spin(newMap, locScoreIJ) {
        unsigned int nOcc = locScoreIJ_spin.cols();
        newMap_spin = Eigen::MatrixXi::Constant(nOcc, nOcc, 0);
      };
      mapVector.push_back(newMap);
    } // for jSys
    maps.push_back(mapVector);
  } // for iSys
  return maps;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>
DirectOrbitalSelection<SCFMode>::buildOrbitalMap(double locThreshold, double kinThreshold) {
  unsigned int nSystems = _orbitalPopulations.size();
  auto maps = initializeOrbitalMap();
  for (unsigned int i = 0; i < nSystems; ++i) {
    for (unsigned int j = i; j < nSystems; ++j) {
      auto& mapIJ = maps[i][j];
      auto& mapJI = maps[j][i];
      const auto& locScoreIJ = _locScore[i][j];
      const auto& kinScoreIJ = _kinScore[i][j];
      const auto& popsI = _orbitalPopulations[i];
      const auto& popsJ = _orbitalPopulations[j];
      for_spin(locScoreIJ, kinScoreIJ, mapIJ, mapJI, popsI, popsJ) {
        unsigned int nOcc = locScoreIJ_spin.cols();
        for (unsigned int k = 0; k < nOcc; ++k) {
          for (unsigned int l = 0; l < nOcc; ++l) {
            double locTr = (i == j && _checkDegeneracies) ? _degeneracyFactor * locThreshold : locThreshold;
            if (_usePiBias && i != j) {
              int nMajorCont_k = (popsI_spin.col(k).array().abs() >= _biasThreshold).count();
              int nMajorCont_l = (popsJ_spin.col(l).array().abs() >= _biasThreshold).count();
              double bias = nMajorCont_k * nMajorCont_l / _biasAverage;
              locTr *= bias;
            } // _usePiBias && i != j
            if (locScoreIJ_spin(k, l) < locTr) {
              mapIJ_spin(k, l) = 1;
              mapJI_spin(l, k) = 1;
            } // if locScoreIJ < locTr
            double kinTr = (i == j && _checkDegeneracies) ? _degeneracyFactor * kinThreshold : kinThreshold;
            if (kinScoreIJ_spin(k, l) > kinTr) {
              mapIJ_spin(k, l) = 0;
              mapJI_spin(l, k) = 0;
            }
          } // for l
        }   // for k
      };
    } // for j
  }   // for i
  return maps;
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> DirectOrbitalSelection<SCFMode>::initializeAssignments() {
  unsigned int nSystems = _orbitalPopulations.size();
  const auto& locScore00 = _locScore[0][0];
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> assignments;
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    SpinPolarizedData<SCFMode, Eigen::VectorXi> hasPartner;
    for_spin(locScore00, hasPartner) {
      unsigned int nOcc = locScore00_spin.cols();
      hasPartner_spin = Eigen::VectorXi::Constant(nOcc, 1);
    };
    assignments.push_back(hasPartner);
  } // for i
  return assignments;
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
DirectOrbitalSelection<SCFMode>::reduceOrbitalMap(std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> map) {
  auto assignments = initializeAssignments();
  unsigned int nSystems = _orbitalPopulations.size();
  for (unsigned int i = 0; i < nSystems; ++i) {
    auto& unpairedVectorI = assignments[i];
    for (unsigned int j = 0; j < nSystems; ++j) {
      auto& mapIJ = map[i][j];
      // row-wise sum for ij-map reduction
      for_spin(mapIJ, unpairedVectorI) {
        Eigen::VectorXi ijImportant = mapIJ_spin.rowwise().sum();
        if (j == 0) {
          unpairedVectorI_spin = ijImportant;
        }
        else {
          // Ensure that the number of orbitals is always the same!
          for (unsigned int iOrb = 0; iOrb < ijImportant.size(); ++iOrb) {
            if (unpairedVectorI_spin(iOrb) != ijImportant(iOrb))
              unpairedVectorI_spin(iOrb) = 0;
          } // for iOrb
        }
        unpairedVectorI_spin = unpairedVectorI_spin.cwiseProduct(ijImportant);
      };
    } // for j
  }   // for i
  return assignments;
}

template<Options::SCF_MODES SCFMode>
void DirectOrbitalSelection<SCFMode>::backMapping(std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments,
                                                  const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> map) {
  unsigned int nSystems = _orbitalPopulations.size();
  // Back mapping. Every orbital that is not found at every point along
  // the reaction coordinate is set to 0
  while (true) {
    bool occChanged = false;
    for (unsigned int i = 0; i < nSystems; ++i) {
      auto& unpairedVectorI = assignments[i];
      for (unsigned int j = 0; j < nSystems; ++j) {
        if (i == j && !_checkDegeneracies)
          continue; // i==j is always the identity
        auto& unpairedVectorJ = assignments[j];
        const auto& mapIJ = map[i][j];
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
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> DirectOrbitalSelection<SCFMode>::getOrbitalAssignments(
    double locThreshold, double kinThreshold, std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& map) {
  // Build initial orbital map and determine non-mappable orbitals.
  map = buildOrbitalMap(locThreshold, kinThreshold);
  auto assignments = reduceOrbitalMap(map);
  backMapping(assignments, map);
  // Build orbital set maps and determine inconsistent orbital mappings.
  auto orbitalSets = buildSystemWiseOrbitalSets(map, assignments);
  auto orbitalSetMap = buildOrbitalSetMap(map, orbitalSets);
  auto orbitalSetAssignments = enforceConsistentOrbitalSetMap(orbitalSetMap);
  // Translate the orbital set mapping to the orbital mapping.
  map = setMapsToOrbitalMap(orbitalSetMap, orbitalSetAssignments, orbitalSets);
  updateOrbitalAssignmentWithSetAssignment(orbitalSets, orbitalSetAssignments, assignments);
  // Ensure that we only have 1 and 0 in the vector.
  for (auto& assign : assignments) {
    for_spin(assign) {
      assign_spin = assign_spin.template cast<bool>().template cast<int>();
    };
  }
  return assignments;
}
template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> DirectOrbitalSelection<SCFMode>::getBestMatchAssignments(
    double& threshold, double& virtThreshold, std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& map,
    const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals, double minThreshold, double maxThreshold) {
  double step = minThreshold / 4.0;
  double current = maxThreshold;
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> occAssignment;
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> virtAssignment;
  if (minThreshold >= maxThreshold)
    throw SerenityError("Minimum threshold must be smaller than the maximum threshold!");

  threshold = maxThreshold;
  virtThreshold = maxThreshold;
  unsigned int bestOccSize = std::numeric_limits<unsigned int>::infinity() - 10;
  unsigned int bestVirtSize = std::numeric_limits<unsigned int>::infinity() - 10;
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> occMap;
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> virtMap;
  bool mapVirtuals = false;
  for_spin(_nOrbitals, nOccupiedOrbitals) {
    if (_nOrbitals_spin > nOccupiedOrbitals_spin)
      mapVirtuals = true;
  };
  while (current >= minThreshold) {
    std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> currentMap;
    auto currentAssignment = this->getOrbitalAssignments(current, current, currentMap);
    auto assignmentZero = currentAssignment[0];
    unsigned int newOccSize = 0;
    unsigned int newVirtSize = 0;
    for_spin(assignmentZero, nOccupiedOrbitals) {
      const auto occAssignment = assignmentZero_spin.head(nOccupiedOrbitals_spin);
      newOccSize += nOccupiedOrbitals_spin - (occAssignment.array() > 0).count();
      if (mapVirtuals) {
        const unsigned int nVirt = assignmentZero_spin.size() - nOccupiedOrbitals_spin;
        const auto virtAssignment = assignmentZero_spin.segment(nOccupiedOrbitals_spin, nVirt);
        newVirtSize += nVirt - (virtAssignment.array() > 0).count();
      }
    };
    if (newOccSize < bestOccSize || !occAssignment.size()) {
      bestOccSize = newOccSize;
      threshold = current;
      occMap = currentMap;
      occAssignment = currentAssignment;
    }
    if (mapVirtuals && (newVirtSize < bestVirtSize || !virtAssignment.size())) {
      bestVirtSize = newVirtSize;
      virtThreshold = current;
      virtMap = currentMap;
      virtAssignment = currentAssignment;
    }
    OutputControl::vOut << current << "  " << newOccSize << " " << newVirtSize << std::endl;
    current -= step;
  }
  map = occMap;
  if (!mapVirtuals) {
    return occAssignment;
  }
  // Combine assignment and maps.
  const unsigned int nSystems = map.size();
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    for (unsigned int jSys = 0; jSys < nSystems; ++jSys) {
      auto& ijMap = map[iSys][jSys];
      const auto& ijVirtMap = virtMap[iSys][jSys];
      for_spin(ijMap, ijVirtMap, nOccupiedOrbitals) {
        const unsigned int nVirt = ijMap_spin.cols() - nOccupiedOrbitals_spin;
        ijMap_spin.block(nOccupiedOrbitals_spin, nOccupiedOrbitals_spin, nVirt, nVirt) =
            ijVirtMap_spin.block(nOccupiedOrbitals_spin, nOccupiedOrbitals_spin, nVirt, nVirt);
      };
    }
    auto& iOccAssignement = occAssignment[iSys];
    const auto& iVirtAssignement = virtAssignment[iSys];
    for_spin(iOccAssignement, iVirtAssignement, nOccupiedOrbitals) {
      const unsigned int nVirt = iOccAssignement_spin.size() - nOccupiedOrbitals_spin;
      iOccAssignement_spin.segment(nOccupiedOrbitals_spin, nVirt) =
          iVirtAssignement_spin.segment(nOccupiedOrbitals_spin, nVirt);
    };
  }
  return occAssignment;
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
DirectOrbitalSelection<SCFMode>::getOrbitalAssignments(double locThreshold, double kinThreshold) {
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> map;
  return getOrbitalAssignments(locThreshold, kinThreshold, map);
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>>
DirectOrbitalSelection<SCFMode>::buildSystemWiseOrbitalSets(
    const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& map,
    const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments) {
  const unsigned int nSystems = map.size();
  std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>> allSets;
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    const auto& zeroMap = map[iSys][iSys];
    SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>> newListOfSets;
    const auto& iAssignments = assignments[iSys];
    for_spin(zeroMap, newListOfSets, iAssignments) {
      const unsigned int nOrbs = zeroMap_spin.cols();
      std::vector<bool> alreadyAssigned(nOrbs, false);
      for (unsigned int iOrb = 0; iOrb < nOrbs; ++iOrb) {
        std::vector<unsigned int> newSet;
        // Do not consider the orbital for orbital sets if it is considered unmappable.
        if (alreadyAssigned[iOrb] || !iAssignments_spin(iOrb))
          continue;
        for (unsigned int jOrb = 0; jOrb < nOrbs; ++jOrb) {
          if (zeroMap_spin(iOrb, jOrb)) {
            newSet.push_back(jOrb);
            alreadyAssigned[jOrb] = true;
          }
        }
        if ((int)newSet.size() != iAssignments_spin(iOrb))
          throw SerenityError("Implementation error! The orbital map construction was probably incorrect!");
        if (newSet.size())
          newListOfSets_spin.push_back(newSet);
      }
    };
    allSets.push_back(newListOfSets);
  }
  return allSets;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> DirectOrbitalSelection<SCFMode>::buildOrbitalSetMap(
    std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& orbitalMap,
    const std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>>& orbitalSets) {
  const unsigned int nSystems = orbitalSets.size();
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> orbitalSetMaps(
      nSystems, std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>(nSystems));
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    const auto& iSets = orbitalSets[iSys];
    for (unsigned int jSys = 0; jSys < iSys; ++jSys) {
      const auto& jSets = orbitalSets[jSys];
      const auto& ijOrbitalMap = orbitalMap[iSys][jSys];
      auto& ijOrbitalSetMap = orbitalSetMaps[iSys][jSys];
      auto& jiOrbitalSetMap = orbitalSetMaps[jSys][iSys];
      for_spin(iSets, jSets, ijOrbitalMap, ijOrbitalSetMap, jiOrbitalSetMap) {
        ijOrbitalSetMap_spin = Eigen::MatrixXi::Zero(iSets_spin.size(), jSets_spin.size());
        jiOrbitalSetMap_spin = Eigen::MatrixXi::Zero(jSets_spin.size(), iSets_spin.size());
        for (unsigned int iSetIndex = 0; iSetIndex < iSets_spin.size(); ++iSetIndex) {
          const auto& iSet = iSets_spin[iSetIndex];
          for (unsigned int jSetIndex = 0; jSetIndex < jSets_spin.size(); ++jSetIndex) {
            const auto& jSet = jSets_spin[jSetIndex];
            if (this->mapBetweenSets(iSet, jSet, ijOrbitalMap_spin)) {
              ijOrbitalSetMap_spin(iSetIndex, jSetIndex) = 1;
              jiOrbitalSetMap_spin(jSetIndex, iSetIndex) = 1;
              break;
            }
          }
        }
      };
    }
  }
  return orbitalSetMaps;
}

template<Options::SCF_MODES SCFMode>
bool DirectOrbitalSelection<SCFMode>::mapBetweenSets(const std::vector<unsigned int>& iSet,
                                                     const std::vector<unsigned int>& jSet, const Eigen::MatrixXi& map) {
  /*
   * Note that we do not check whether there any sets not explicitly between the sets i and j.
   * We already ensure that this is not the case in the function reduceOrbitalMap(...).
   */
  if (iSet.size() != jSet.size())
    return false;
  for (const auto& iOrb : iSet) {
    for (const auto& jOrb : jSet) {
      if (!map(iOrb, jOrb))
        return false;
    }
  }
  return true;
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> DirectOrbitalSelection<SCFMode>::enforceConsistentOrbitalSetMap(
    const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& orbitalSetMap) {
  unsigned int nSystems = orbitalSetMap.size();
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> setAssignments(nSystems);
  for (unsigned iSys = 0; iSys < nSystems; ++iSys) {
    auto& iSetAssignment = setAssignments[iSys];
    const unsigned int jSizeIndex = (iSys == 0) ? 1 : 0;
    const auto& iiSetMap = orbitalSetMap[jSizeIndex][iSys];
    for_spin(iSetAssignment, iiSetMap) {
      iSetAssignment_spin = Eigen::VectorXi::Ones(iiSetMap_spin.cols());
    };
    for (unsigned int jSys = 0; jSys < iSys; ++jSys) {
      const auto& ijSetMap = orbitalSetMap[iSys][jSys];
      for_spin(iSetAssignment, ijSetMap) {
        const Eigen::VectorXi ijMapExists = ijSetMap_spin.colwise().sum();
        iSetAssignment_spin.array() *= ijMapExists.array();
      };
    }
  }
  return setAssignments;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> DirectOrbitalSelection<SCFMode>::setMapsToOrbitalMap(
    const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& orbitalSetMap,
    const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& orbitalSetAssignment,
    const std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>>& orbitalSets) {
  auto orbitalMap = initializeOrbitalMap();
  const unsigned int nSystems = orbitalSetMap.size();
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    const auto& iSets = orbitalSets[iSys];
    const auto& iAssign = orbitalSetAssignment[iSys];
    for (unsigned int jSys = 0; jSys <= iSys; ++jSys) {
      const auto& jSets = orbitalSets[jSys];
      const auto& jAssign = orbitalSetAssignment[jSys];
      const auto& ijSetMap = orbitalSetMap[iSys][jSys];
      auto& ijOrbitalMap = orbitalMap[iSys][jSys];
      for_spin(ijSetMap, ijOrbitalMap, iSets, jSets, iAssign, jAssign) {
        for (unsigned int iSetIndex = 0; iSetIndex < iSets_spin.size(); ++iSetIndex) {
          if (!iAssign_spin(iSetIndex))
            continue;
          for (unsigned int jSetIndex = 0; jSetIndex < jSets_spin.size(); ++jSetIndex) {
            const bool validMapping = (iSys == jSys) ? true : ijSetMap_spin(iSetIndex, jSetIndex);
            if (!jAssign_spin(jSetIndex) || !validMapping)
              continue;
            for (auto& i : iSets_spin[iSetIndex]) {
              for (auto& j : jSets_spin[jSetIndex]) {
                ijOrbitalMap_spin(i, j) = 1;
              }
            }
          }
        }
      };
    }
  }
  return orbitalMap;
}

template<Options::SCF_MODES SCFMode>
void DirectOrbitalSelection<SCFMode>::updateOrbitalAssignmentWithSetAssignment(
    const std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>>& orbitalSets,
    const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& orbitalSetAssignment,
    std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments) {
  const auto nSystems = assignments.size();
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    const auto& iSets = orbitalSets[iSys];
    const auto& iSetAssignments = orbitalSetAssignment[iSys];
    auto& orbitalAssignments = assignments[iSys];
    for_spin(iSets, iSetAssignments, orbitalAssignments) {
      for (unsigned int iSetIndex = 0; iSetIndex < iSets_spin.size(); ++iSetIndex) {
        if (!iSetAssignments_spin[iSetIndex]) {
          for (const auto& i : iSets_spin[iSetIndex]) {
            orbitalAssignments_spin(i) = 0;
          }
        }
      }
    };
  }
}

template class DirectOrbitalSelection<Options::SCF_MODES::RESTRICTED>;
template class DirectOrbitalSelection<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
