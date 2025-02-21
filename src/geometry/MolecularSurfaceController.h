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

/* Include Serenity Internal Headers */
#include "grid/GridController.h"               //Inherits from.
#include "notification/ObjectSensitiveClass.h" //Notify.
#include "settings/PCMSettings.h"              //PCMSettings as member
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <string>

namespace Serenity {

class Geometry;
class MolecularSurface;
/**
 * @brief A class for the possible usage scenarios of the surface.
 *   FDE:        Use for an embedding calculations. The surface may cover multiple systems.
 *   ACTIVE:     Use for an isolated system. The surface will cover only one system.
 *   ACTIVE_VDW: Use in the case of calculating the cavity formation energy. The Van-der-Waals surface for one system.
 */
enum class MOLECULAR_SURFACE_TYPES { FDE, ACTIVE, ACTIVE_VDW };

/**
 * @class MolecularSurfaceController MolecularSurfaceController.h
 * @brief A class that holds all information of a molecular surface required to perform
 *        CPCM and IEFPCM calculations. It can construct the molecular surface by itself.
 *        This, it will reconstruct it if changes to the geometry occur.
 *
 *
 *
 *   The implementation is based on: [1] Chemical Reviews, 2005, Vol. 105, No. 83013\n\n
 *   All notations are taken from the reference.
 *   Notations:\n
 *     \f$ a_i \f$ : Area of tesserae i.\n
 *     \f$ \pmb{s}_i \f$ : Coordinates of the representative point of tesserae i.\n
 *     \f$ \pmb{n}_i \f$ : Normal vector on the surface at point \f$ \pmb{s}_i \f$.\n
 *     \f$ k \f$ : Constant factor k. Values of \f$ k = 1.0694 \f$ or \f$ k = 1.07 \f$
 *                 are commonly used.
 *     \f$ R_I \f$ : The radius of sphere \f$ I \f$.
 *
 *    This class inherits from GridController. The associated grid is the molecular surface.
 *      Any object that may depend on the surface, can use the NotifyingClass<Grid> functionality
 *      of GridController.
 *
 *    This class inherits from ObjectSensitiveClass<Geometry>. It will reconstruct the grid
 *      if any changed to the underlying geometry occur.
 */
class MolecularSurfaceController : public GridController, ObjectSensitiveClass<Geometry> {
 public:
  /**
   * @brief Constructor.
   * @param geometry     The geometry.
   * @param pcmSettings  The PCM settings.
   */
  MolecularSurfaceController(std::shared_ptr<Geometry> geometry, const PCMSettings& pcmSettings);

  virtual ~MolecularSurfaceController() = default;
  /**
   * @brief Getter for the normal vectors.
   * @return The normal vectors \f$ \pmb{n}_i \f$.
   */
  const Eigen::Matrix3Xd& getNormalVectors();
  /**
   * @brief Getter for the normal vectors as a pointer.
   * @return The normal vectors \f$ \pmb{n}_i \f$.
   */
  std::shared_ptr<Eigen::Matrix3Xd> getNormalVectors_ptr();
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
  /**
   * @brief The grid and any associated matrices will be deleted.
   */
  void notify() override;
  /**
   * @brief Getter for the underlying grid as a molecular surface.
   */
  const MolecularSurface& getMolecularSurface();
  /**
   * @brief Forwarded getter for the grid point to sphere/atom mapping of the surface.
   */
  const std::vector<unsigned int>& getPointToSphereMapping();
  /**
   * @brief Forwarded getter for the sphere to surface point mapping.
   */
  const std::vector<std::pair<unsigned int, unsigned int>>& getSphereToPointMapping();
  /**
   * @brief Getter for the grid point coordinates. Overrides implementation in GridController.h
   * @return The grid point coordinates.
   */
  const Eigen::Matrix3Xd& getGridPoints() override;
  /**
   * @brief Getter for the grid point weight. Overrides implementation in GridController.h
   * @return The grid point weight.
   */
  const Eigen::VectorXd& getWeights() override;
  /**
   * @brief Getter for the number of grid points. Overrides implementation in GridController.h
   * @return The number of grid points..
   */
  unsigned int getNGridPoints() override;
  /**
   * @brief Getter for the underlying geometry.
   * @return The geometry.
   */
  std::shared_ptr<Geometry> getGeometry();
  /**
   * @brief Getter for the cavity formation energy of the cavity.
   *        The Pierotti--Claverie  expression derived from scaled particle theory is used.
   *        J. Phys. Chem., 1988, 92, 1617-1631
   * @return The cavity formation energy.
   */
  double getCavityEnergy();

  /**
   * @brief Setter for the underlying MolecularSurface.
   * @param surface The molecular surface.
   */
  void setSurface(std::unique_ptr<MolecularSurface>&& surface);

  /**
   * @brief Boolean value to save the information if the molecular surface has been loaded from a file.
   */
  bool isLoaded();

  /**
   * @brief Getter for the path of the charges belonging to a molecular surface.
   */
  std::string getChargesPath();

  /**
   * @brief Print some information about the surface.
   */
  void printInfo();

 private:
  // The geometry.
  std::shared_ptr<Geometry> _geometry;
  // The PCMSettings.
  PCMSettings _pcmSettings;
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

  std::shared_ptr<double> _cavityEnergy;
  // The constant factor used for the diagonal entries of S and D.
  const double _k = 1.0694;
  // (Re-)Build the molecular surface.
  void buildSurface();
  // Calculates the cavity energy.
  void calculateCavityEnergy();
  // Calculates the cavity energy for a spherical cavity with radius r_m.
  double calculateG_cav_sphere(double r_s, double rho, double r_m, double t);
  // Calculates the part of the surface of the sphere exposed to the solvent.
  double calculateExposedSurface_sphere(const std::pair<unsigned int, unsigned int> indices, const Eigen::VectorXd& weights);
  // Getter for the number density of the solvent. Takes care of the 'EXPLICIT' solvent.
  double getNumberDensity();

  double getCavityFormationProbeRadius();
};

} /* namespace Serenity */

#endif /* GEOMETRY_MOLECULARSURFACECONTROLLER_H_ */
