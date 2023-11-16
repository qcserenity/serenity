/**
 * @file   PadeApproximation.h
 *
 * @date   09.10, 2020
 * @author Johannes Toelle
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

#ifndef PADEAPPROXIMATION_H_
#define PADEAPPROXIMATION_H_

/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {
/**
 * @class PadeApproximation PadeApproximation.h
 * @brief A class, which calculates determines the Pade coefficients for numerical \n
 *
 * Reference: H. J. Vidberg and J. W. Serene, J. Low Temp. Phys. 29, 179 (1977)\n
 *
 */
class PadeApproximation {
 public:
  /**
   * @brief Constructor
   * @param points The Points for which the padecoefficients are fitted
   * @param functionValues The corresponding functionvalues to the points
   */
  PadeApproximation(Eigen::VectorXcd points, Eigen::VectorXcd functionValues);
  /// @brief Default destructor.
  virtual ~PadeApproximation() = default;
  /**
   * @brief Calculates the complex pade coefficients
   * @return The complex pade coefficients
   */
  const Eigen::VectorXcd calculatePadeCoeffs();
  /**
   * @brief Performs the pade approximation for a sepcific value
   * @param value The value for which the pade approximation is performed
   * @return The corresponding function value via the pade approximation
   */
  const std::complex<double> padeApproximation(std::complex<double> value);

 private:
  /// @brief The points used for the pade approximation
  Eigen::VectorXcd _points;
  /// @brief The corresponding functionvalues to the points
  Eigen::VectorXcd _functionValues;
  /// @brief The number of points
  unsigned int _numberOfPoints;
  /// @brief The fitted pade coefficients
  Eigen::VectorXcd _padeCoeffs;
};

} // namespace Serenity

#endif /* PADEAPPROXIMATION_H_ */