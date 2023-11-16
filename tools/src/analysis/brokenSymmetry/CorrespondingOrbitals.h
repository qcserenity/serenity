/**
 * @file   CorrespondingOrbitals.h
 *
 * @date   Dec 14, 2017
 * @author Anja Massolle
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

#ifndef ANALYSIS_CORRESPONDINGORBITALS_H_
#define ANALYSIS_CORRESPONDINGORBITALS_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

/**
 * @class CorrespondingOrbitals CorrespondingOrbitals.h
 * @brief Performs a corresponding orbital transformation of two given orbital sets. The orbitals are transformed
 * in such a way that their MO overlap matrix is diagonal. The orbital set resulting from this transformation
 * are the corresponding orbitals. \n
 * According to: \n
 * - Proc. R. Soc. Lond. A. 263, 483-493
 * - J. Chem. Phys. 47, 1936-1941
 * - Int. J. Quantum Chem. 45, 587-590
 * - J. Phys. Chem. Solids 65, 781-785
 */
template<Options::SCF_MODES SCFMode>
class CorrespondingOrbitals {
 public:
  /**
   * @brief Constructor
   * @param sys1 First system for the corresponding orbital transformation.
   * @param sys2 Second system for the corresponding orbital transformation.
   */
  CorrespondingOrbitals(std::shared_ptr<SystemController> sys1, std::shared_ptr<SystemController> sys2);

  /**
   * @brief Constructor
   * @param coeff1 First coefficient matrix for the corresponding orbital transformation.
   * @param coeff2 Second coefficient matrix for the corresponding orbital transformation.
   * @param nOcc1 Number of occupied orbitals in the first coefficient matrix
   * @param nOcc2 Number of occupied orbitals in the second coefficient matrix
   */
  CorrespondingOrbitals(CoefficientMatrix<SCFMode> coeff1, CoefficientMatrix<SCFMode> coeff2,
                        SpinPolarizedData<SCFMode, unsigned int> nOcc1, SpinPolarizedData<SCFMode, unsigned int> nOcc2);

  /**
   * @brief Destructor
   */
  virtual ~CorrespondingOrbitals() = default;

  /**
   * @return the corresponding orbitals.
   */
  std::pair<SpinPolarizedData<SCFMode, Eigen::MatrixXd>, SpinPolarizedData<SCFMode, Eigen::MatrixXd>> getCorrespondingOrbitals();

  /**
   * @return the overlap matrix between the corresponding orbitals.
   */
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> getSPrimePrime();

  /**
   * @brief Prints the overlap of the corresponding orbitals.
   * @param n Number of the orbital of interest
   */
  SpinPolarizedData<SCFMode, double> getOverlap(SpinPolarizedData<SCFMode, unsigned int> n);

 private:
  // occupied orbitals of the original set 1
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _coeff1;
  // occupied orbitals of the original set 2
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _coeff2;
  // amount of occupied orbitals of the original set 1
  SpinPolarizedData<SCFMode, unsigned int> _nOcc1;
  // amount of occupied orbitals of the original set 2
  SpinPolarizedData<SCFMode, unsigned int> _nOcc2;
  // U matrix of the SVD decomposition of _SMO
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _U;
  // V matrix of the SVD decomposition of _SMO
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _V;
  // The MO overlap matrix between _coeff1 and _coeff2
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _SMO;
  // The AO overlap matrix
  Eigen::MatrixXd _SAO;
};
} /* namespace Serenity */

#endif /* ANALYSIS_CORRESPONDINGORBITALS_H_ */
