/**
 * @file SphericalHarmonicsRotations.h
 *
 * @author Moritz Bensberg
 * @date Feb 11, 2020
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

#ifndef BASIS_SPHERICALHARMONICSROTATIONS_H_
#define BASIS_SPHERICALHARMONICSROTATIONS_H_

/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices
#include <map>         //std::map (storing of transformation matrices)
#include <memory>      //smrt_ptr

namespace Serenity {
/**
 * @class SphericalHarmonicsRotations SphericalHarmonicsRotations.h
 * @brief Calculate the transformation matrices for a rotation of a set of real, spherical harmonics.
 *        The resulting matrix contains the linear coefficients which are used to expand the rotated
 *        function in the basis of the (fixed axis) initial, real spherical harmonics.\n\n
 *
 *        The transformation \f$ \mathbf{R}^l(\alpha , \beta , \gamma) \f$ has the property\n
 *        \f$ \mathbf{C}^l_\mathrm{rot} = \mathbf{R}^l(\alpha , \beta , \gamma) \mathbf{C}^l_\mathrm{initial}\f$ .\n
 *
 *        The transformation matrices are documentet in the manual and taken from:\n
 *           Pinchon, D.; Hoggan, P. E. J. Phys. A: Math. Theor. 2007, 40, 1597–1610.
 *
 *
 */
class SphericalHarmonicsRotations {
 private:
  // Purely static, never instantiated.
  SphericalHarmonicsRotations();
  virtual ~SphericalHarmonicsRotations();
  ///@brief A map that stores the transformation matrices J_l.
  static std::map<unsigned int, std::shared_ptr<Eigen::MatrixXd>> _jMatrices;
  ///@brief Calculate the helper matrix hat(G)_z.
  static Eigen::MatrixXd calculateHatGZMatrix(unsigned int l);
  ///@brief Calculate the helper matrix G_y.
  static Eigen::MatrixXd calculateGYMatrix(unsigned int l);
  ///@brief Calculate the rotation matrix X(angle) for a given l.
  static Eigen::MatrixXd calculateXMatrix(unsigned int l, double alpha);

 public:
  /**
   * @brief Getter for the transformation matrix J_l. Note that this function is only
   *        public so that it can be tested easily.
   * @param l The angular momentum.
   * @return J_l.
   */
  static const Eigen::MatrixXd& getJMatrix(unsigned int l);
  /**
   * @brief Getter for the rotation matrix in the basis of the solid spherical harmonics.
   * @param l     The angular momentum.
   * @param alpha The Euler angle alpha-->rotation around z-axis.
   * @param beta  The Euler angle beta -->rotation around y-axis.
   * @param gamma The Euler angle gamm -->rotation around z-axis.
   * @return The rotation matrix.
   */
  static Eigen::MatrixXd getTransformationMatrix(unsigned int l, double alpha, double beta, double gamma);
};

} /* namespace Serenity */

#endif /* BASIS_SPHERICALHARMONICSROTATIONS_H_ */
