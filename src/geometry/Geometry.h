/**
 * @file   Geometry.h
 *
 * @date   Mar 19, 2013
 * @author Thomas Dresselhaus
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
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


/* Forward declarations */
namespace HDF5 {
class Filepath;
}
namespace Serenity {
template<typename> class Matrix;

/**
 * @class Geometry Geometry.h
 *
 * @brief Holds a set of atoms.
 *
 */
class Geometry : public ObjectSensitiveClass<Atom> {
  friend class GradientCalculator;
public:
  /**
   * @brief Constructor.
   * @param atoms
   */
  explicit Geometry(std::vector<std::shared_ptr<Atom> > atoms);
  /**
   * @brief Constructor.
   * @param atomSymbols    The atoic symbols.
   * @param atomPositions  The atom positions, dimensions: (nAtoms,3)
   */
  explicit Geometry(std::vector<std::string> atomSymbols,
                    Matrix<double> atomPositions);
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
   * TODO see SDGeomOptimizer.cpp, maybe we should const-overload this function and in the const
   * case return a vector<shared_ptr<const Atom> >.
   */
  inline const std::vector<std::shared_ptr<Atom> >& getAtoms() const {
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
  inline double getCoreCoreRepulsion() const {
    if (!_coreCoreRepulsionUpToDate) calcCoreCoreRepulsion();
    return _coreCoreRepulsion;
  }
  /**
   * @returns The center of mass
   */
  inline Point getCenterOfMass() const {
    if (!_centerOfMassUpToDate) calcCenterOfMass();
    return _centerOfMass;
  }

  void calcMomentOfInertiaX();
  void calcMomentOfInertiaY();
  void calcMomentOfInertiaZ();

  double getMomentOfInertia(){
    if (!_momentOfInertia){
      _momentOfInertia = getMomentOfInertiaX() + getMomentOfInertiaY() + getMomentOfInertiaZ();
    }
    return _momentOfInertia;
  }
  double getMomentOfInertiaX(){
    if (!_momentOfInertiaX) calcMomentOfInertiaX();
    return _momentOfInertiaX;
  }
  double getMomentOfInertiaY(){
    if (!_momentOfInertiaY) calcMomentOfInertiaY();
    return _momentOfInertiaY;
  }
  double getMomentOfInertiaZ(){
    if (!_momentOfInertiaZ) calcMomentOfInertiaZ();
    return _momentOfInertiaZ;
  }

  /**
   * @brief A schortcut to acess atoms directly.
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
  std::vector<std::string> getAtomSymbols() const{
    std::vector<std::string> symbols;
    for (auto& atom: this->_atoms){
      symbols.push_back(atom->getAtomType()->getElementSymbol());
    }
    return symbols;
  }

  /**
   * @brief Deletes an atom from the geometry.
   * @param i The number of the atom to be deleted within the array
   */
  void deleteAtom(unsigned int i){
    std::vector<std::shared_ptr<Atom> >::iterator it=_atoms.begin();
    std::advance(it, i);
    _atoms.erase(it);
    notify();
  };

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
  void setCoordinates(Eigen::MatrixXd& newCoordinates) const;

//  /**
//   * @brief Set all coordinates at once
//   * @param newCoordinates The new coordinates, dimensions: (nAtoms,3)
//   */
//  void setCoordinates(const Eigen::VectorXd& newCoordinates) const;

  /**
   * @brief Get the gradients if they are all up to date
   * @return The geometry gradient, dimensions: (nAtoms,3)
   */
  Matrix<double> getGradients() const;

  /**
   * @brief Set the Gradients.
   * @param newGradients The new geometry gradient, dimensions: (nAtoms,3)
   */
  void setGradients(Matrix<double>& newGradients) const;

//  /**
//   * @brief Set the Gradients.
//   * @param newGradients The new geometry gradient, dimensions: (nAtoms*3)
//   */
//  void setGradients(Eigen::VectorXd& newGradients) const;

  /**
   * @brief Check if there are gradients that are up-to-date
   */
  bool checkGradients(){
    bool check = true;
    for (unsigned int i=0;i!=this->getNAtoms();++i){
      check *= this->_atoms[i]->gradientsUpToDate();
    }
    return check;
  }

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
   * @param rhs is added to this instance.
   */
  Geometry& operator+=(const Geometry& rhs) {
    for (auto atom : rhs.getAtoms()) {
      this->_atoms.push_back(atom);
      atom->addSensitiveObject(this->_self);
      if (atom->getX() < this->_minX) this->_minX = atom->getX();
      if (atom->getY() < this->_minY) this->_minY = atom->getY();
      if (atom->getZ() < this->_minZ) this->_minZ = atom->getZ();
      if (atom->getX() > this->_maxX) this->_maxX = atom->getX();
      if (atom->getY() > this->_maxY) this->_maxY = atom->getY();
      if (atom->getZ() > this->_maxZ) this->_maxZ = atom->getZ();
    }
    return *this;
  }

  void addAsDummy(const Geometry& add, bool toFront = false);

  bool hasAtomsWithECPs() const {
    bool atomsWithECPs =false;
    for (auto& atom : _atoms) {
      if (atom->usesECP()) atomsWithECPs = true;
    }
    return atomsWithECPs;
  }

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

  std::vector<std::shared_ptr<Atom> > _atoms;

  double _minX, _minY, _minZ;
  double _maxX, _maxY, _maxZ;
  mutable double _coreCoreRepulsion;
  mutable Point _centerOfMass;
  mutable bool _coreCoreRepulsionUpToDate;
  mutable bool _centerOfMassUpToDate;
  double _momentOfInertia = 0.0;
  double _momentOfInertiaX = 0.0;
  double _momentOfInertiaY = 0.0;
  double _momentOfInertiaZ = 0.0;
};

} /* namespace Serenity */
#endif /* GEOMETRY_H_ */
