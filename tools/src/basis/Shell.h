/**
 * @file Shell.h
 *
 * @date Nov 7, 2016
 * @author Jan Unsleber, Thomas Dresselhaus
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

#ifndef SHELL_H_
#define SHELL_H_

/* Include Serenity Internal Headers */
#include "notification/NotifyingClass.h"
#include "parameters/Constants.h"

/* Include Std and External Headers */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <libint2/shell.h>
#pragma GCC diagnostic pop
#include <Eigen/Dense>
#include <vector>

namespace Serenity {
/**
 * @class Shell Shell.h
 * @brief decorated version of a Libint2 shell.
 *
 * Each Shell can consist of several primitive functions. With the angular momentum
 * more than one basis function is easily implied. So for an angular momentum of 1, i.e. a p-type
 * function implicitly p_x, p_y and p_z are defined.
 * This results in the once in a while encountered 'reduced' and 'extended' counting of the basis
 * functions. 'Reduced' then means that shells of functions are counted, while 'extended' counts
 * the actually used (contracted) basis functions.\n
 * Assume e.g. a basis with 2 s-type, 2 p-type and 1 d-type function. Reduced counting then results
 * in 2+2+1=5, and extended counting in 2\*1 + 2\*3 + 1\*6 = 14 (Cartesian functions are assumed
 * here).
 *
 * For details of what a basis function is consider the book 'Modern Quantum Chemistry' by Szabo
 * and Ostlund.
 */
class Shell : public libint2::Shell, public NotifyingClass<Shell> {
 public:
  /**
   * @brief Constructor.
   *
   * Shells will be normalize automatically by libint.
   *
   * @param exponents        The exponents for the primitives.
   * @param contractions     The contractions of the primitives (only one set).
   * @param angularMomentum  The angular momentum.
   * @param spherical        Boolean, true if the shell is given in spherical functions.
   * @param coords           The reference coordinates of the shell.
   */
  Shell(libint2::svector<double> exponents, libint2::svector<double> contractions, unsigned int angularMomentum,
        bool spherical, std::array<double, 3> coords, std::string element = "");
  Shell(libint2::svector<double> exponents, libint2::svector<double> exponents1, libint2::svector<double> contractions,
        unsigned int angularMomentum, bool spherical, std::array<double, 3> coords, std::string element = "");
  /**
   * @brief Copy constructor.
   *
   * @param other            The original shell.
   */
  Shell(const Shell& other);

  /**
   * @brief Default destructor.
   */
  virtual ~Shell() = default;

  /// @returns Returns the number of primitives forming the basis function.
  inline unsigned int getNPrimitives() const {
    return this->alpha.size();
  }

  /// @returns Returns the number of contracted functions defined by this shell.
  inline unsigned int getNContracted() const {
    return !(this->contr[0].pure) ? N_SHELL_CART[this->contr[0].l] : N_SHELL_SPH[this->contr[0].l];
  }

  /// @returns the angular momentum
  inline unsigned int getAngularMomentum() const {
    return (unsigned int)(this->contr[0].l);
  }
  /**
   * @returns Returns true iff the basis function shell is to be interpreted with cartesian gaussians, i.e.
   *          false if it shall be interpreted with pure spherical harmonics.
   */
  inline bool isCartesian() const {
    return !(this->contr[0].pure);
  }
  /// @returns the opposite of isCartesian
  inline bool isSpherical() const {
    return this->contr[0].pure;
  }
  /// @returns the x-Coordinate at which the basis function is centered.
  inline const double& getX() const {
    return this->O[0];
  }
  /// @returns the y-Coordinate at which the basis function is centered.
  inline const double& getY() const {
    return this->O[1];
  }
  /// @returns the z-Coordinate at which the basis function is centered.
  inline const double& getZ() const {
    return this->O[2];
  }

  /// @param x new x coordinate
  inline void setX(double x) {
    this->O[0] = x;
    this->notifyObjects();
  }

  /// @param y new y coordinate
  inline void setY(double y) {
    this->O[1] = y;
    this->notifyObjects();
  }

  /// @param z new z coordinate
  inline void setZ(double z) {
    this->O[2] = z;
    this->notifyObjects();
  }
  /// @param add_to_x The component to be added to the x coordinate.
  inline void addToX(double add_to_x) {
    setX(this->O[0] + add_to_x);
  }
  /// @param add_to_y The component to be added to the y coordinate.
  inline void addToY(double add_to_y) {
    setY(this->O[1] + add_to_y);
  }
  /// @param add_to_z The component to be added to the z coordinate.
  inline void addToZ(double add_to_z) {
    setZ(this->O[2] + add_to_z);
  }

  inline const std::string getElement() const {
    return _element;
  }

  /**
   * @return Returns the normalization factors needed to calculate the values of
   *         the contracted functions on the grid.
   *         The data is stored in a vector with the contracted functions ordered
   *         the same way as in the libint2:shells.
   *         The normalization factors are given for each basis function e.g. px, py ... .
   */
  const Eigen::VectorXd& getNormFactors() const {
    return *_normFactors;
  }

  /**
   * @return Returns the contractions of the basis functions.
   */
  const libint2::svector<double> getContractions() const {
    return _contractions;
  }
  /**
   * @return Returns the contractions of the normalized basis functions.
   */
  const libint2::svector<double> getNormContractions() const {
    return this->contr[0].coeff;
  }
  /**
   * @return Returns the exponents of the basis functions.
   */
  const libint2::svector<double> getExponents() const {
    return _exponents;
  }
  /**
   * @brief Equal operator. Comparison is done by:\n
   *          origin,
   *          angular momentum,
   *          exponents,
   *          contractions,
   *          spherical/cartesian
   * @param other The other shell.
   * @return True if considered to be equal, else false.
   */
  bool operator==(const Shell& other) const;

 private:
  ///@brief Normalization factors for the basis functions.
  std::unique_ptr<Eigen::VectorXd> _normFactors;
  ///@brief Contractions of the primitive Gaussians.
  libint2::svector<double> _contractions;
  ///@brief Exponents of the primitive Gaussians.
  libint2::svector<double> _exponents;
  ///@brief Element identifier as a string.
  std::string _element;
};

} /* namespace Serenity */

#endif /* SHELL_H_ */
