/**
 * @file DipoleApproximationToPairEnergies.h
 *
 * @date Feb 1, 2019
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

#ifndef POSTHF_MPN_DIPOLEAPPROXIMATIONTOPAIRENERGIES_H_
#define POSTHF_MPN_DIPOLEAPPROXIMATIONTOPAIRENERGIES_H_
/* Include Serenity Internal Headers */
#include "data/matrices/FockMatrix.h"
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dense matrices
#include <Eigen/SparseCore> //Sparse maps
#include <memory>           //smart ptr
#include <vector>           //std::vector

namespace Serenity {
class OrbitalPair;
class BasisController;
class PAOController;

/**
 * @class
 * @brief Calculates a dipole approximations to the pair energies of orbital pairs ij.
 *
 * Reference:\n
 *       J.Chem.Phys. 143, 034108 (2015).\n
 */
class DipoleApproximationToPairEnergies {
 public:
  /**
   * @brief Constructor.
   * @param basisController The basis controller.
   * @param overlapMatrix The overlap matrix.
   * @param occCoefficients The coefficients of the occupied orbitals. Sorted by columns.
   * @param paoController The controller for the PAO coefficients.
   * @param fockMatrix The fock matrix.
   * @param paoToOccupiedOrbitalMap The map between PAOs and occupied orbitals.
   * @param paoOrthogonalizationThreshold The orthogonalization threshold for the construction
   *        of the quasi canonical, linear independent PAO domains.
   */
  DipoleApproximationToPairEnergies(std::shared_ptr<BasisController> basisController,
                                    std::shared_ptr<MatrixInBasis<RESTRICTED>> overlapMatrix,
                                    std::shared_ptr<Eigen::MatrixXd> occCoefficients,
                                    std::shared_ptr<PAOController> paoController,
                                    std::shared_ptr<FockMatrix<RESTRICTED>> fockMatrix,
                                    std::shared_ptr<Eigen::SparseMatrix<int>> paoToOccupiedOrbitalMap,
                                    double paoOrthogonalizationThreshold);
  /**
   * @brief Default Destructor.
   */
  virtual ~DipoleApproximationToPairEnergies() = default;
  /**
   * @brief Calculates the dipole approximation.
   * @param orbitalPairs The orbital pairs.
   */
  void calculateDipoleApproximation(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs);
  /**
   * @brief Calculates the dipole approximation by assuming collinear orientation between
   *        the dipoles (i|r|i) and (j|r|j).
   * @param orbitalPairs The orbital pairs.
   */
  void calculateDipoleApproximationCollinear(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs);

 private:
  // The basis controller.
  std::shared_ptr<BasisController> _basisController;
  // The overlap matrix.
  std::shared_ptr<MatrixInBasis<RESTRICTED>> _overlapMatrix;
  // The coefficients of the occupied orbitals.
  std::shared_ptr<Eigen::MatrixXd> _occCoefficients;
  // The PAO controller.
  std::shared_ptr<PAOController> _paoController;
  // The fock matrix.
  std::shared_ptr<FockMatrix<RESTRICTED>> _fockMatrix;
  // The map between PAOs and occupied orbitals.
  std::shared_ptr<Eigen::SparseMatrix<int>> _paoToOccupiedOrbitalMap;
  // The orthogonalization threshold for the construction
  // of the quasi canonical, linear independent PAO domains.
  double _paoOrthogonalizationThreshold;
  // The dipole integrals for every orbital. Columns transformed to occupied basis.
  std::vector<Eigen::MatrixXd> _dipoleInts;
  /**
   * @brief Calculates the dipole integrals (mu|x|nu),(mu|y|nu) and (mu|z|nu).
   * @return A vector containing the dipole integrals (x-->0, y-->1, z-->2).
   */
  inline std::vector<Eigen::MatrixXd> calculateDipoleIntegrals();
  ///@brief QC-eigenvalues for the diagonal pairs.
  std::vector<std::shared_ptr<Eigen::VectorXd>> _eigenvalues_ii;
  ///@brief Integrals (i|r|i)
  Eigen::Matrix3Xd _iri;
  ///@brief Integrals (mu|r|i), mu in [ii] PAO domain.
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _pao_r_i;

  /**
   * @brief Calculate non-redundant, quasi-canonical PAOs for the given orbital i.
   * @param i The orbital index.
   */
  inline void calculateTransformationAndEigenvalues(unsigned int i);

  /**
   * @brief Precalculate the transformations to the non-redundant, quasi-canonical PAOs for
   *        all orbitals contained in at least one orbital pair.
   * @param orbitalPairs The orbital pairs.
   */
  inline void preCalculateTransformations(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs);
};

} /* namespace Serenity */

#endif /* POSTHF_MPN_DIPOLEAPPROXIMATIONTOPAIRENERGIES_H_ */
