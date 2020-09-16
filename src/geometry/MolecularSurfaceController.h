/**
 * @file MolecularSurfaceController.h
 *
 * @author Moritz Bensberg
 * @date May 25, 2020
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

#ifndef GEOMETRY_MOLECULARSURFACECONTROLLER_H_
#define GEOMETRY_MOLECULARSURFACECONTROLLER_H_

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <string>

namespace Serenity {

class GridController;

enum class MOLECULAR_SURFACE_TYPES { FDE, ACTIVE };

/**
 * @class MolecularSurfaceController MolecularSurfaceController.h
 * @brief A class that holds all information of a molecular surface required to perform
 *        CPCM and IEFPCM calculations. It needs to be rebuild if the geometry changes!\n\n
 *
 *   The implementation is based on: [1] Chemical Reviews, 2005, Vol. 105, No. 83013\n\n
 *   All notations are taken from the reference.
 *   Notations:\n
 *     \f$ a_i \f$ : Area of tessarae i.\n
 *     \f$ \pmb{s}_i \f$ : Coordinates of the representative point of tessarae i.\n
 *     \f$ \pmb{n}_i \f$ : Normal vector on the surface at point \f$ \pmb{s}_i \f$.\n
 *     \f$ k \f$ : Constant factor k. Values of \f$ k = 1.0694 \f$ or \f$ k = 1.07 \f$
 *                 are commonly used.
 *     \f$ R_I\f$ : Radius of the sphere associated to the point \f$ \pmb{s}_i \f$.\n.
 */
class MolecularSurfaceController {
 public:
  /**
   * @brief Constructor.
   * @param gridController   The grid controller of the surface grid.
   * @param normalVectors    The normal vectors for each surface point, pointing
   *                         outwards from the molecular surface away.
   * @param surfaceLabel     The label of the controlled surface.
   */
  MolecularSurfaceController(std::shared_ptr<GridController> gridController,
                             std::unique_ptr<Eigen::Matrix3Xd> normalVectors, std::string surfaceLabel);
  /**
   * @brief Constructor.
   * @param coordinates      The coordinates of the surface points.
   * @param weights          The areas/integration weights associated to each point.
   * @param normalVectors    The normal vectors for each surface point, pointing
   *                         outwards from the molecular surface away.
   */
  MolecularSurfaceController(std::unique_ptr<Eigen::Matrix3Xd> coordinates, std::unique_ptr<Eigen::VectorXd> weights,
                             std::unique_ptr<Eigen::Matrix3Xd> normalVectors, std::string surfaceLabel);
  virtual ~MolecularSurfaceController() = default;
  /**
   * @brief Getter for the underlying grid controller.
   * @return The grid controller.
   */
  std::shared_ptr<GridController> getGridController();
  /**
   * @brief Forwarded getter for the surface point coordinates.
   * @return The surface coordinates (\f$ \pmb{s}_i \f$).
   */
  const Eigen::Matrix3Xd& getGridPoints();
  /**
   * @brief Forwarded getter for the weights/areas associated to each surface point.
   * @return The surface weights/areas \f$ a_i \f$.
   */
  const Eigen::VectorXd& getWeights();
  /**
   * @brief Getter for the normal vectors.
   * @return The normal vectors \f$ \pmb{n}_i \f$.
   */
  const Eigen::Matrix3Xd& getNormalVectors();
  /**
   * @brief Getter for the inverse distance matrix \f$ \pmb{S} \f$.\n
   *           \f$ S_{ij} = \frac{1}{|\pmb{s}_i - \pmb{s}_j|}, i \neq j \f$\n
   *           \f$ S_{ii} = k \sqrt(\frac{4\pi}{a_i}) \f$\n
   *           see [1] tab. 1 and eq. 35.
   * @return The matrix S.
   */
  const Eigen::MatrixXd& getMatrixS();
  /**
   * @brief Getter for the inverse of the inverse distance matrix \f$ \pmb{S} \f$.\n
   *           \f$ S_{ij} = \frac{1}{|\pmb{s}_i - \pmb{s}_j|}, i \neq j \f$\n
   *           \f$ S_{ii} = k \sqrt(\frac{4\pi}{a_i}) \f$\n
   *           see [1] tab. 1 and eq. 35.
   * @return The matrix S^-1.
   */
  const Eigen::MatrixXd& getMatrixSinv();
  /**
   * @brief Getter for the matrix \f$ \pmb{D} \f$
   *           \f$ D_ij = \frac{(\pmb{s}_i-\pmb{s}_j)\pmb{n}_j{|\pmb{s}_i-\pmb{s}_j|^3}, i \neq j \f$.\n
   *           For the diagonal entries two approximations are used:
   *             \f$ D_ii = k \frac{sqrt(4\pi a_i)}{2 R_I} \f$\n or \n
   *             \f$ D_ii = -(2\pi+\sum_{i\neq j}D_{ij}a_j)/a_i \f$.\n
   *           We use the latter one, since it does not require any (meta) information about
   *           the spheres which may have been used in the construction of the original surface.
   *           However, there is no reason to constrain us to using a sphere based construction
   *           algorithm.\n
   *           See [1] tab. 1 and eqs. 36,37.
   * @return The matrix D.
   */
  const Eigen::MatrixXd& getMatrixD();
  /**
   * @brief Getter for the diagonal matrix A. The diagonal entries are the weights \f$ a_i \f$.
   * @return The matrix A.
   */
  const Eigen::MatrixXd& getMatrixA();
  /**
   * @brief Getter for the inverse of the diagonal matrix A. The diagonal entries are \f$ 1/a_i \f$.
   * @return The matrix A^-1.
   */
  const Eigen::MatrixXd& getMatrixAinv();

 private:
  // The grid controller. Holds areas and coordinates.
  std::shared_ptr<GridController> _gridController;
  // The normal vectors.
  std::unique_ptr<Eigen::Matrix3Xd> _normalVectors;
  // The matrix S.
  std::unique_ptr<Eigen::MatrixXd> _S;
  // The inverse of S.
  std::unique_ptr<Eigen::MatrixXd> _invS;
  // The matrix D.
  std::unique_ptr<Eigen::MatrixXd> _D;
  // The matrix A.
  std::unique_ptr<Eigen::MatrixXd> _A;
  // The inverse of A.
  std::unique_ptr<Eigen::MatrixXd> _Ainv;
  // The constant factor used for the diagonal entries of S and D.
  const double _k = 1.0694;
  // The surface label.
  std::string _surfaceLabel;
  // Print some information about the surface.
  void printInfo();
};

} /* namespace Serenity */

#endif /* GEOMETRY_MOLECULARSURFACECONTROLLER_H_ */
