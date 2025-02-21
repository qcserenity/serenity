/**
 * @file   Geometry.h
 *
 * @date   Mar 19, 2013
 * @author Thomas Dresselhaus
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
#ifndef GEOMETRY_H_
#define GEOMETRY_H_
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <cmath>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
template<typename>
class Matrix;

/**
 * @class Geometry Geometry.h
 *
 * @brief Holds a set of atoms.
 *
 */
class Geometry : public ObjectSensitiveClass<Atom>, public NotifyingClass<Geometry> {
 public:
  /**
   * @brief Constructor.
   * @param atoms
   */
  explicit Geometry(std::vector<std::shared_ptr<Atom>> atoms);
  /**
   * @brief Constructor.
   * @param atomSymbols    The atomic symbols.
   * @param atomPositions  The atom positions, dimensions: (nAtoms,3)
   */
  explicit Geometry(std::vector<std::string> atomSymbols, Matrix<double> atomPositions);
  /**
   * @brief Constructor for Geometries that are filled via addition.
   */
  Geometry();
  /**
   * @brief Destructor.
   */
  virtual ~Geometry() = default;
  /**
   * @returns the atoms forming this Geometry
   */
  inline const std::vector<std::shared_ptr<Atom>>& getAtoms() const {
    return _atoms;
  }
  /**
   * @returns number of atoms
   */
  inline unsigned int getNAtoms() const {
    return _atoms.size();
  }
  /**
   * @returns The energy due to the Coulomb repulsion of the atoms from each other
   */
  double getCoreCoreRepulsion() const;
  /**
   * @returns The center of mass
   */
  Point getCenterOfMass() const;
  /**
   * @brief Calculate and return the sum of all effective charges.
   * @return The total effective charge.
   */
  int getTotalEffectiveCharge();

  /**
   * @brief A shortcut to acess atoms directly.
   *
   * @param   i
   * @returns the ith atom of the Geometry.
   * TODO see getAtoms(), this method would need a const-overload then as well. (Here it's actually easy.)
   */
  inline std::shared_ptr<Atom> operator[](unsigned int i) const {
    return _atoms[i];
  }

  /**
   * @returns the largest x coordinate of the underlying atoms
   */
  inline double getMaxX() const {
    return _maxX;
  }
  /**
   * @returns the largest y coordinate of the underlying atoms
   */
  inline double getMaxY() const {
    return _maxY;
  }
  /**
   * @returns the largest z coordinate of the underlying atoms
   */
  inline double getMaxZ() const {
    return _maxZ;
  }
  /**
   * @returns the smallest (or most negative) x coordinate of the underlying atoms
   */
  inline double getMinX() const {
    return _minX;
  }
  /**
   * @returns the smallest (or most negative) y coordinate of the underlying atoms
   */
  inline double getMinY() const {
    return _minY;
  }
  /**
   * @returns the smallest (or most negative) z coordinate of the underlying atoms
   */
  inline double getMinZ() const {
    return _minZ;
  }

  /**
   * @return A list of all atoms symbols.
   */
  std::vector<std::string> getAtomSymbols() const;
  /**
   * @brief Deletes an atom from the geometry.
   * @param i The number of the atom to be deleted within the array
   */
  void deleteAtom(unsigned int i);

  /**
   * @brief Prints the current geometry to the screen.
   */
  void print() const;

  /**
   * @brief Prints the current geometry to file.
   */
  void printToFile(std::string baseName, std::string id) const;

  /**
   * @brief Adds the current geometry to the trajectory file
   * @param baseName The basename consisting of path + systemname
   * @param energy The energy associated to this geometry
   * @param gradNorm The gradient associated to this geometry
   */
  void updateTrajFile(std::string baseName, double energy, double gradNorm = -1.0) const;

  /**
   * @brief Prints the current geometry gradients to the screen.
   */
  void printGradients() const;

  /**
   * @brief Get the coordinates of all atoms in a Matrix
   * @return The coordinates, dimensions: (nAtoms,3)
   */
  Matrix<double> getCoordinates() const;

  /**
   * @brief Get the aligned coordinates of all atoms in a Matrix
   *
   * These coordinates are way of removing major RMSDs between similar structures.
   * For any system the the heaviest atom is taken, as a tie breaker the closest
   * heaviest atom to the center of mass is used.
   * This atom is placed in the origin, the second heaviest atom is then placed
   * on the z-axis (-z) and the third heaviest atom is placed on the zy-plane
   * (-y,-z).
   *
   * @return The coordinates, dimensions: (nAtoms,3)
   */
  Matrix<double> getAlignedCoordinates() const;

  /**
   * @brief Set all coordinates at once
   * @param newCoordinates The new coordinates, dimensions: (nAtoms,3)
   */
  void setCoordinates(const Eigen::MatrixXd& newCoordinates) const;
  /**
   * @brief Get the gradients if they are all up to date
   * @return The geometry gradient, dimensions: (nAtoms,3)
   */
  Matrix<double> getGradients() const;
  /**
   * @brief Set the Gradients.
   * @param newGradients The new geometry gradient, dimensions: (nAtoms,3)
   */
  void setGradients(const Eigen::MatrixXd& newGradients) const;

  /**
   * @brief Deletes coinciding atoms. Deletion of dummy atoms prioritized.
   */
  void deleteIdenticalAtoms();
  /**
   * @brief Construct and get translational modes
   * @return the 3 translational modes
   */
  Eigen::MatrixXd getTransModes();

  /**
   * @brief Construct and get rotational modes
   * @return the 3 rotational modes (one of them is a zero vector for linear molecules)
   */
  Eigen::MatrixXd getRotModes();

  /**
   * @brief Make gradients translationally invariant
   */
  void makeGradientsTranslationallyInvariant();

  /**
   * @brief Linearity check
   * @return True if the molecules is linear.
   */
  bool isLinear();

  /**
   * @brief Check for atoms with the same coordinates.
   * @return The result of the check.
   */
  bool hasIdenticalAtoms() const;

  /**
   * @param rhs is added to this instance.
   */
  void operator+=(const Geometry& rhs);

  /**
   * @brief checks if two geometries are equal.
   * @return True if geometries are equal.
   */
  bool operator==(const Geometry& rhs);
  /**
   * @brief checks if two geometries are not equal.
   * @return True if geometries are not equal.
   */
  bool operator!=(const Geometry& rhs) {
    return (!(*this == rhs));
  }

  /**
   * @brief Add geometry as dummy atoms.
   * @param add The other geometry.
   * @param toFront Add at the front of the atom list.
   */
  void addAsDummy(const Geometry& add, bool toFront = false);

  /**
   * @brief Add the geometry dummy atoms.
   * @param add The other geometry.
   * @param toFront Add at the front of the atom list.
   */
  void addDummy(const Geometry& add, bool toFront = false);

  /**
   * @brief Checks for atoms with ECPs.
   * @return The result of the check.
   */
  bool hasAtomsWithECPs() const;
  /**
   * @brief Delete all ghost atoms.
   */
  void deleteGhostAtoms();
  /**
   * @brief Enforce an update of the core--core repulsion.
   */
  void updateCoreCoreRepulsion();
  /**
   * @brief Getter for the number of non-valence electrons associated to this geometry.
   *        This already omit the electrons contained in ECPs.
   * @return The number of non-valence/core electroncs.
   */
  unsigned int getNumberOfCoreElectrons();
  /**
   * @brief Getter for the number of minimal-basis functions associated to this geometry.
   * @param excludeDummyAtoms If true, the basis functions on dummy atoms are not counted. By default false.
   * @return The number of minimal basis functions.
   */
  unsigned int getNMinimalBasisFunctions(bool excludeDummyAtoms = false) const;

  /**
   * @brief Calculates the distance between the two closest atoms of two geometries.
   * @return The number of non-valence/core electroncs.
   */
  double getMinimumDistance(const Geometry& rhs);

 private:
  /**
   * @brief Calculates the core core repulsion of the current geometry
   */
  void calcCoreCoreRepulsion() const;
  /**
   * @brief Calculates the center of mass of the current geometry
   */
  void calcCenterOfMass() const;
  /**
   * @brief Initializes updates.
   */
  void notify();

  std::vector<std::shared_ptr<Atom>> _atoms;

  double _minX, _minY, _minZ;
  double _maxX, _maxY, _maxZ;
  mutable double _coreCoreRepulsion;
  mutable Point _centerOfMass;
  mutable bool _coreCoreRepulsionUpToDate;
  mutable bool _centerOfMassUpToDate;
};

} /* namespace Serenity */
#endif /* GEOMETRY_H_ */
