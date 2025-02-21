/**
 * @file DelleySurfaceConstructor.h
 *
 * @date   May 27, 2020
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

#ifndef GEOMETRY_DELLEYSURFACECONSTRUCTOR_H_
#define GEOMETRY_DELLEYSURFACECONSTRUCTOR_H_
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward Declarations */
class Sphere;
class GridController;
class Ellipse;
class Plane;
class MolecularSurface;

/**
 * @struct
 * @brief Container for the coordinates and weights associated with a unit sphere.
 */
struct UnitSphere {
 public:
  /**
   * @brief Default constructor.
   */
  UnitSphere() = default;
  /**
   * @brief The coordinates on the unit sphere.
   */
  Eigen::Matrix3Xd coordinates;
  /**
   * @brief The weights of the points.
   */
  Eigen::VectorXd weights;
};

/**
 * @class
 * @brief A class that handles the construction of a molecular surface according to
 *        B. Delley:\n
 *          Molecular Simulation, Vol. 32, No. 2, 15 February 2006, 117–123
 *          https://doi.org/10.1080/08927020600589684
 *
 *        The molecular surface is defined as a continuous, analytic function that
 *        varies with the atom positions. The resulting surface has no discontinuities
 *        at the atom-sphere intersections. The Grid points are obtained by atom-wise
 *        projection of spherical grid points used for DFT-grids to the molecular surface.
 *        See the SerenityUserManual.tex/pdf for a detailed description.
 */
class DelleySurfaceConstructor {
 public:
  /**
   * @brief Constructor. Constructs auxiliary parameters.
   * @param spheres           The spheres used for the surface construction.
   * @param r_s               Radius of the solvent probe.
   * @param alpha             The sharpness parameter.
   * @param projectionCutOff  Upper limit for the binary search for the 0 surface.
   */
  DelleySurfaceConstructor(std::vector<Sphere> spheres, double r_s, double alpha, double projectionCutOff,
                           double minimalDistance, bool oneCavity, double connectivityFactor);
  /**
   * @brief Default destructor.
   */
  ~DelleySurfaceConstructor();
  /**
   * @brief Getter for the molecular surface..
   * @return The grid controller.
   */
  std::unique_ptr<MolecularSurface> getMolecularSurface();

 private:
  /**
   * @brief Construct the surface.
   */
  void buildSurface();
  /**
   * @brief Get the spherical grid points on a unit sphere for a given angular momentum.
   * @param l  The angular momentum.
   * @return The unit sphere.
   */
  const UnitSphere& getUnitSphere(unsigned int l);
  /**
   * @brief Construct the unit sphere.
   * @param The angular momentum for the spherical points.
   */
  std::shared_ptr<UnitSphere> buildUnitSphere(unsigned int l);
  /**
   * @brief Project a point to the molecular surface. Performs a binary search.
   * @param center        The sphere center.
   * @param unitPoint     The initial point on the unit sphere.
   * @param initialGuess  The initial guess (VdW radius).
   * @param r             The coordinates up the point. Set within the function.
   * @param grad          The function gradient at the position of the point. Set within the function.
   * @return If True: The point can not be projected to the surface. If false: The projection was successful.
   */
  bool projectCenterOntoSurface(const Eigen::Vector3d& center, const Eigen::Vector3d& unitPoint, double initialGuess,
                                Eigen::Vector3d& r, Eigen::Vector3d& grad);
  /**
   * @brief Project a circle from the plane touching the unit surface to a plane on the surface with the given
   *        normal vector. This projection will result in an ellipse if both planes are not parallel.
   * @param centerOnUnit      The center of the circle on the unit sphere.
   * @param normalOnUnit      The normal vector defining the plane on the unit sphere.
   * @param centerOnSurface   The center of the ellipse.
   * @param normalOnSurface   The normal vector on the surface at the center of the ellipse.
   * @param weightOnUnit      The initial weight on the unit surface.
   * @return The ellipse resulting from the projection.
   */
  std::shared_ptr<Ellipse> projectCircleOntoSurface(const Eigen::Vector3d& centerOnUnit, const Eigen::Vector3d& normalOnUnit,
                                                    const Eigen::Vector3d& centerOnSurface,
                                                    const Eigen::Vector3d& normalOnSurface, double weightOnUnit);
  /**
   * @brief Calculate the function value and its gradient at position r.
   * @param r      The point.
   * @param value  The value. Set in this function.
   * @return The gradient at point r.
   */
  Eigen::Vector3d calculateGradientAndValue(const Eigen::Vector3d& r, double& value);
  /**
   * @brief Calculate the function value  at position r.
   * @param r      The point.
   * @return The function value F(r).
   */
  double calculateFunctionValue(const Eigen::Vector3d& r);
  /**
   * @brief Get the projection of r onto the bond between r_i and r_j.
   *        The point r is assumed to be associated to the center r_i.
   * @param r_i   The first center.
   * @param r_j   The second center.
   * @param r     The point.
   * @return The projection of r onto the bond. If r does not fall on the bond,
   *         the position of the bond end in question is returned.
   */
  Eigen::Vector3d getBondProjection(const Eigen::Vector3d& r_i, const Eigen::Vector3d& r_j, const Eigen::Vector3d& r);
  /**
   * @brief Get the plane orthogonal to the bond between two spheres.
   * @param sphere_i   The first sphere.
   * @param sphere_j   The second sphere.
   * @return The normal plane.
   */
  Plane getSphereSphereTouchingPlane(const Sphere& sphere_i, const Sphere& sphere_j);
  /**
   * @brief Calculate the radii and parameters of all bond-cylinders within the molecule.
   */
  void calculateCylinderRadiiAndParameters();
  /// The bond-cylinder radii.
  Eigen::SparseMatrix<double> _cylinderRadii;
  /// The cylinder parameters (2.0 * _r_s * std::max(R_z, 1.0))
  /// for each pair of spheres.
  /// Here, _r_s is the solvent radius and R_z is the cylinder radius.
  /// If no bond exists between spheres, the value is set to -1.
  Eigen::MatrixXd _cylinderParameters;
  /// The spheres.
  std::vector<Sphere> _spheres;
  /// The solvent radius.
  const double _r_s;
  /// The sharpness parameter.
  const double _alpha;
  /// Upper limit for the binary search for the 0 surface.
  const double _projectionCutOff;
  /// Minimal distance between points belonging to one sphere.
  double _minimalDistance;
  /// Assume that all cavity points are "connected" to each other.
  /// Remove all other points.
  bool _oneCavity;
  /// The connection between point is given for all points closer
  /// _connectivityFactor * r_s.
  double _connectivityFactor;
  /// The auxiliary parameter a.
  double _a;
  /// The auxiliary parameter b.
  double _b;
  /// The auxiliary function cut off (4.0/_alpha)
  const double _cutOff;
  /// Map between unit spheres and angular momentum.
  static std::map<unsigned int, std::shared_ptr<UnitSphere>> _unitSpheres;
  /// The auxiliary function f(x).
  inline double f(double x) const {
    return -exp(-_alpha * x) + _a * x + _b * x * x;
  }
  /// The derivative of the auxiliary function f(x) with respect to x.
  inline double df(double x) const {
    return _alpha * exp(-_alpha * x) + _a + 2.0 * _b * x;
  }
  /// @brief The radius of a cylinder touching spheres with radius R_i, R_j and _r_s
  ///        with the distance d_ij between the sphere i and j.
  /// @param R_i   Radius of sphere i.
  /// @param R_j   Radius of sphere j.
  /// @param d_ij  Distance between spheres i and j.
  /// @return The cyllinder radius.
  inline double R_z_ij(double R_i, double R_j, double d_ij) const {
    const double S_i = R_i + _r_s;
    const double S_j = R_j + _r_s;
    const double S_i2 = S_i * S_i;
    const double S_j2 = S_j * S_j;
    const double tmp = (S_i2 + d_ij * d_ij - S_j2) / (2.0 * d_ij);
    const double R_z_ij = sqrt(S_i2 - tmp * tmp) - _r_s;
    return R_z_ij;
  }
  /**
   * @brief Remove all points that are not connected to the point with the most
   *        extreme x-coordinates.
   * @param coords    The coordinates. They will be changed.
   * @param norms     The normal vectors. They will be changed.
   * @param weights   The weights. They will be changed.
   *
   * The diameter(or the _connectivityFactor times _r_s) of the probe-sphere is used
   * to construct a connectivity matrix between points. All points are kept that can
   * be reached via this matrix from the point with the most extreme x-coordinate.
   */
  void groupPoints(std::vector<Eigen::Vector3d>& coords, std::vector<Eigen::Vector3d>& norms,
                   std::vector<double>& weights, std::vector<unsigned int>& sphereIndices);
  /**
   * @brief Get the indices of all points connected arcording to connections to the given seed.
   * @param connections   The connections. 0 for non-connected.
   * @param seed          The seed.
   * @return  The connected indices.
   */
  std::vector<unsigned int> buildPointCloud(Eigen::SparseMatrix<int>& connections, unsigned int seed);
  /**
   * @biref Construct the mapping between spheres and grid points.
   * @param sphereIndices The sphere indices.
   * @param nSpheres      The total number of spheres.
   * @return The map from spheres to cavity points as pairs. All cavity points are ordered according to
   *         their sphere index. Each pair holds the first and last index of the cavity point of the sphere
   *         that corresponds to its vector-index.
   */
  std::vector<std::pair<unsigned int, unsigned int>> collectSphereIndices(std::vector<unsigned int> sphereIndices,
                                                                          unsigned int nSpheres);
  /**
   * @brief Construct from the thread-wise lists of cavity points, the final list.
   * @param allCoordinates   Thread-wise coordinates.
   * @param allNorms         Thread-wise norms.
   * @param allWeights       Thread-wise weights.
   * @param allSphereIndices Thread-wise sphere indices.
   * @param coords           Final list of coordinates.
   * @param norms            Final list of norms.
   * @param weights          Final list of weights.
   * @param sphereIndices    Final list of sphere indices.
   * @param spheres          List of spheres.
   */
  void mergeGridPoints(const std::vector<std::vector<Eigen::Vector3d>>& allCoordinates,
                       const std::vector<std::vector<Eigen::Vector3d>>& allNorms,
                       const std::vector<std::vector<double>>& allWeights,
                       const std::vector<std::vector<unsigned int>>& allSphereIndices, std::vector<Eigen::Vector3d>& coords,
                       std::vector<Eigen::Vector3d>& norms, std::vector<double>& weights,
                       std::vector<unsigned int>& sphereIndices, const std::vector<Sphere>& spheres);
};

} /* namespace Serenity */

#endif /* GEOMETRY_DELLEYSURFACECONSTRUCTOR_H_ */
