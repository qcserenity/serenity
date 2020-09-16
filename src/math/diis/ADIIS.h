/**
 * @file ADIIS.h
 *
 * @date Apr 24, 2017
 * @author Jan Unsleber
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

#ifndef MATH_DIIS_ADIIS_H_
#define MATH_DIIS_ADIIS_H_

/* Include Serenity Internal Headers */
#include "math/optimizer/LBFGS.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <vector>

namespace Serenity {
/**
 * @class ADIIS ADIIS.h
 * @brief An implementation of the A-DIIS.
 *
 * For the main reference please see:
 *  Accelerating self-consistent field convergence with the
 *    augmented Roothaan–Hall energy function
 *  Xiangqian Hu and Weitao Yang,
 *    J. Chem. Phys. (2010) 132(5):054109.
 *    https://dx.doi.org/10.1063%2F1.3304922
 */
class ADIIS {
 public:
  /// @brief Constructor.
  ADIIS();
  /// @brief Default destructor.
  virtual ~ADIIS() = default;
  /**
   * @brief Function to optimize the Fock matrix.
   * @param energy The current (electronic) energy (\f$ E_n \f$).
   * @param F The current Fock matrix (\f$ F_n \f$).
   * @param P The current density matrix (\f$ F_n \f$).
   * @return The new Fock Matrix.
   */
  Eigen::MatrixXd optimize(const Eigen::MatrixXd& F, const Eigen::MatrixXd& P);

 private:
  /// @brief Number of matrices stored;
  unsigned int _n;
  /// @brief Maximum number of matrices stored;
  static constexpr unsigned int _maxn = 6;
  /// @brief Stored density matrices.
  std::vector<Eigen::MatrixXd> _dmat;
  /// @brief Stored fock matrices.
  std::vector<Eigen::MatrixXd> _fmat;
};

} /* namespace Serenity */

#endif /* MATH_DIIS_ADIIS_H_ */
