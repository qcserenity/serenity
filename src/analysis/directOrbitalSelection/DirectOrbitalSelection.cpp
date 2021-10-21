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
#include "misc/SerenityError.h" //Errors.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
DirectOrbitalSelection<SCFMode>::DirectOrbitalSelection(const std::vector<SPMatrix<SCFMode>>& orbitalPopulations,
                                                        const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>> kineticEnergies,
                                                        const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals,
                                                        bool usePiBias, double biasThreshold, double biasAverage,
                                                        bool checkDegeneracies, double degeneracyFactor)
  : _orbitalPopulations(orbitalPopulations),
    _kineticEnergies(kineticEnergies),
    _usePiBias(usePiBias),
    _biasThreshold(biasThreshold),
    _biasAverage(biasAverage),
    _checkDegeneracies(checkDegeneracies),
    _degeneracyFactor(degeneracyFactor) {
  if (orbitalPopulations.size() != kineticEnergies.size() || orbitalPopulations.size() < 1)
    throw SerenityError("ERROR: Insufficient information provided for the direct orbital selection. Not enough orbital "
                        "populations or kinetic energies.");
  // Prepare comparison by precalculating all differences.
  calculateLocScore(nOccupiedOrbitals);
  calculateKinScore(nOccupiedOrbitals);
}

template<Options::SCF_MODES SCFMode>
DirectOrbitalSelection<SCFMode>::~DirectOrbitalSelection() = default;

template<Options::SCF_MODES SCFMode>
void DirectOrbitalSelection<SCFMode>::calculateLocScore(const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals) {
  initializeScores(_locScore, nOccupiedOrbitals);
  unsigned int nSystems = _orbitalPopulations.size();
  for (unsigned int i = 0; i < nSystems; ++i) {
    const auto& popsI = _orbitalPopulations[i];
    for (unsigned int j = i; j < nSystems; ++j) {
      const auto& popsJ = _orbitalPopulations[j];
      auto& scoreIJ = _locScore[i][j];
      auto& scoreJI = _locScore[j][i];
      for_spin(popsI, popsJ, scoreIJ, scoreJI) {
        const unsigned int nOccI = scoreIJ_spin.rows();
        const unsigned int nOccJ = scoreIJ_spin.cols();
        for (unsigned int k = 0; k < nOccI; ++k) {
          for (unsigned int l = 0; l < nOccJ; ++l) {
            double populationMeasure = (popsI_spin.col(k) - popsJ_spin.col(l)).array().abs().sum();
            scoreIJ_spin(k, l) = populationMeasure;
            scoreJI_spin(l, k) = populationMeasure;
          } // for l
        }   // for k
      };
    } // for j
  }   // for i
}

template<Options::SCF_MODES SCFMode>
void DirectOrbitalSelection<SCFMode>::calculateKinScore(const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals) {
  initializeScores(_kinScore, nOccupiedOrbitals);
  unsigned int nSystems = _orbitalPopulations.size();
  for (unsigned int i = 0; i < nSystems; ++i) {
    const auto& kinsI = _kineticEnergies[i];
    for (unsigned int j = i; j < nSystems; ++j) {
      const auto& kinsJ = _kineticEnergies[j];
      auto& scoreIJ = _kinScore[i][j];
      auto& scoreJI = _kinScore[j][i];
      for_spin(kinsI, kinsJ, scoreIJ, scoreJI) {
        const unsigned int nOccI = scoreIJ_spin.rows();
        const unsigned int nOccJ = scoreIJ_spin.cols();
        for (unsigned int k = 0; k < nOccI; ++k) {
          for (unsigned int l = 0; l < nOccJ; ++l) {
            double kineticEnergyMeasure = std::fabs(kinsI_spin[k] - kinsJ_spin[l]);
            scoreIJ_spin(k, l) = kineticEnergyMeasure;
            scoreJI_spin(l, k) = kineticEnergyMeasure;
          } // for l
        }   // for k
      };
    } // for j
  }   // for i
}

template<Options::SCF_MODES SCFMode>
void DirectOrbitalSelection<SCFMode>::initializeScores(std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>>& scores,
                                                       const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals) {
  unsigned int nSystems = _orbitalPopulations.size();
  for (unsigned int iSys = 0; iSys < nSystems; ++iSys) {
    std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> scoreVector;
    for (unsigned int jSys = 0; jSys < nSystems; ++jSys) {
      SpinPolarizedData<SCFMode, Eigen::MatrixXd> newScores;
      for_spin(newScores, nOccupiedOrbitals) {
        newScores_spin = Eigen::MatrixXd::Constant(nOccupiedOrbitals_spin, nOccupiedOrbitals_spin, 0);
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
            double locTr = (i == j) ? _degeneracyFactor * locThreshold : locThreshold;
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
            double kinTr = (i == j) ? _degeneracyFactor * kinThreshold : kinThreshold;
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
      if (i == j)
        continue; // i==j is always the identity
      auto& mapIJ = map[i][j];
      // row-wise sum for ij-map reduction
      for_spin(mapIJ, unpairedVectorI) {
        Eigen::VectorXi ijImportant = mapIJ_spin.rowwise().sum();
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
        if (i == j && _checkDegeneracies)
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
  map = buildOrbitalMap(locThreshold, kinThreshold);
  auto assignments = reduceOrbitalMap(map);
  backMapping(assignments, map);
  // Ensure that we only have 1 and 0 in the vector.
  for (auto& assign : assignments) {
    for_spin(assign) {
      assign_spin = assign_spin.template cast<bool>().template cast<int>();
    };
  }
  return assignments;
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
DirectOrbitalSelection<SCFMode>::getOrbitalAssignments(double locThreshold, double kinThreshold) {
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> map;
  return getOrbitalAssignments(locThreshold, kinThreshold, map);
}

template class DirectOrbitalSelection<Options::SCF_MODES::RESTRICTED>;
template class DirectOrbitalSelection<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
