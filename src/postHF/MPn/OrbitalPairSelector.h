/**
 * @file OrbitalPairSelector.h
 *
 * @date Apr 2, 2019
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

#ifndef POSTHF_MPN_ORBITALPAIRSELECTOR_H_
#define POSTHF_MPN_ORBITALPAIRSELECTOR_H_
/* Include Serenity Internal Headers */
#include "data/matrices/FockMatrix.h" //Fock matrix definition
/* Include Std and External Headers */
#include <memory> //smart ptr.
#include <vector> //std::vector

namespace Serenity {
/* Forward Declarations */
class OrbitalPair;
class SystemController;
class PAOController;

/**
 * @class OrbitalPairSelector OrbitalPairSelector.h
 * @brief A class that selects significant orbital pairs from a given orbital pair list.
 *        DOI and dipole approximation are used as criteria.
 */
class OrbitalPairSelector {
 public:
  /**
   * @brief Constructor.
   * @param system The system controller associated to the orbital pairs
   *               (Needed for Settings, overlap matrix, basis, occupied orbitals).
   * @param paoController The PAOController for the given orbital pairs.
   */
  OrbitalPairSelector(std::shared_ptr<SystemController> system, std::shared_ptr<PAOController> paoController)
    : _systemController(system), _paoController(paoController){};
  /**
   * @brief Default destructor.
   */
  virtual ~OrbitalPairSelector() = default;

  /**
   * @brief Worker function. Separates the initial pairs into significant and insignificant pairs.
   * @param initialPairs The initial orbital pair list.
   * @param mnpPreThreshold Shell prescreening threshold for the DOIs.
   * @param paoOrthoThreshold The renormalization threshold for removing linear dependencies from
   *                          the PAO domains (domains used in the dipole approximation).
   * @param f The fock matrix of the system.
   * @param paoToOccupiedOrbitalMap Map between PAOs and occupied orbitals.
   * @param doiThreshold Significants threshold for the DOI.
   * @param collinearDipolePairThreshold Significants threshold for the dipole approximation.
   * @return A pair of orbital pair lists. .first--> significant pairs; .second-->insignificant pairs.
   */
  std::pair<std::vector<std::shared_ptr<OrbitalPair>>, std::vector<std::shared_ptr<OrbitalPair>>>
  selectOrbitalPairs(std::vector<std::shared_ptr<OrbitalPair>> initialPairs, double mnpPreThreshold,
                     double paoOrthoThreshold, std::shared_ptr<FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                     std::shared_ptr<Eigen::SparseMatrix<int>> paoToOccupiedOrbitalMap, double doiThreshold);

 private:
  ///@brief The system controller.
  std::shared_ptr<SystemController> _systemController;
  ///@brief The PAO controller.
  std::shared_ptr<PAOController> _paoController;

  /* Helper Functions */
  /**
   * @brief Calculates the DOIs for the given initialPairs and writes them into the actual pair.
   * @param initialPairs The initial pairs.
   * @param mnpPreThreshold The shell prescreening threshold.
   */
  void calculateDOIs(std::vector<std::shared_ptr<OrbitalPair>> initialPairs, double mnpPreThreshold);

  /**
   * @brief Calculates the dipole and collinear dipole approximation for the given initialPairs.
   * @param initialPairs The initial pairs.
   * @param paoOrthoThreshold The PAO renormalization threshold after removal of linear dependencies.
   * @param f The fock matrix.
   * @param paoToOccupiedOrbitalMap Map between PAOs and occupied orbitals.
   */
  void calculateDipoleApproximation(std::vector<std::shared_ptr<OrbitalPair>> initialPairs, double paoOrthoThreshold,
                                    std::shared_ptr<FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                                    std::shared_ptr<Eigen::SparseMatrix<int>> paoToOccupiedOrbitalMap);
};

} /* namespace Serenity */

#endif /* POSTHF_MPN_ORBITALPAIRSELECTOR_H_ */
