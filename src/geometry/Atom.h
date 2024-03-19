/**
 * @file   Atom.h
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
#ifndef ATOM_H_
#define ATOM_H_
/* Include Serenity Internal Headers */
#include "geometry/AtomType.h"
#include "geometry/Point.h"
#include "math/FloatMaths.h"
#include "misc/SerenityError.h"
#include "notification/NotifyingClass.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <vector>

/* External forward declaration */
namespace libecpint {
struct ECP;
}
namespace Serenity {
/* Forward declarations */
class AtomGrid;
class Shell;
/**
 * @class Atom Atom.h
 * @brief Model class for an atom
 */
class Atom : public Point, public NotifyingClass<Atom> {
  friend class Geometry;

 public:
  /**
   * @param atomType is typically resolved to an element symbol, but will maybe be generalized
   *                 to an atom type as it is used e.g. in force fields.
   * @param x,y,z    coordinates of the atom (atomic units)
   */
  Atom(std::shared_ptr<const AtomType> atomType, double x, double y, double z);
  /**
   * @param atomType is typically resolved to an element symbol, but will maybe be generalized
   *                 to an atom type as it is used e.g. in force fields.
   * @param x,y,z    coordinates of the atom (atomic units)
   * @param basisFunctions can be directly associated to the atom upon construction.
   */
  Atom(std::shared_ptr<const AtomType> atomType, double x, double y, double z,
       std::pair<std::string, std::vector<std::shared_ptr<Shell>>> basisFunctions);
  /**
   * @brief checks if two atoms are equal.
   * @return True if atoms are equal.
   * @param rhs other atom.
   */
  bool operator==(Atom rhs);
  /**
   * @param symbol   The atom symbol in a string.
   * @param x,y,z    Coordinates of the atom (atomic units).
   */
  Atom(const std::string symbol, double x, double y, double z);

  virtual ~Atom() = default;

  /// @returns the atom type (NOT necessarily an element symbol!!)
  std::shared_ptr<const AtomType> getAtomType() const {
    return _atomType;
  }

  /**
   * @brief Dummy check.
   * @return Returns true if the atom type is a dummy atom.
   */
  bool isDummy() {
    return _atomType->isDummy();
  }
  /**
   * @returns the number of basis functions for this atom of the active basis
   * CAUTION! REDUCED COUNTING
   */
  unsigned int getNBasisFunctions() const;
  /**
   * @param   label a string identifying the basis type
   * @returns the number of basis functions for this atom of the basis with the label label
   * CAUTION! REDUCED COUNTING
   */
  unsigned int getNBasisFunctions(std::string label) const;
  /// @returns the basis functions for this atom of the active basis
  std::vector<std::shared_ptr<Shell>>& getBasisFunctions();
  /// @returns the basis functions for this atom of the active basis
  const std::vector<std::shared_ptr<Shell>>& getBasisFunctions() const;
  /**
   * @param   label a string identifying the basis type
   * @returns the basis functions for this atom of the basis with the label label
   */
  std::vector<std::shared_ptr<Shell>>& getBasisFunctions(std::string label);
  /**
   * @param   label a string identifying the basis type
   * @returns Returns if a basis with the given label is present
   */
  bool basisFunctionsExist(std::string label);

  /**
   * @param   label a string identifying the basis type
   * @returns the basis functions for this atom of the basis with the label label
   */
  const std::vector<std::shared_ptr<Shell>>& getBasisFunctions(std::string label) const;
  /**
   * @returns the nuclear charge which determines the element symbol. Caution: the effective
   *          charge of the nucleus may be different if effective core potentials are used.
   */
  inline int getNuclearCharge() const {
    return _atomType->getNuclearCharge();
  }
  /**
   * @returns the effective charge of the nucleus which is the same as getNuclearCharge() unless
   *          an effective core potential is used.
   */
  inline int getEffectiveCharge() const {
    return getNuclearCharge() - _nECPElectrons;
  }
  /**
   * @brief Supplies this atom with another set of basis functions.
   *
   * @param newBasis  the additional basis with an identifying label
   * @param isPrimary if true the new basis is marked as being the active one from now on.
   * TODO extract these functions into the cpp file!!!
   */
  void addBasis(std::pair<std::string, std::vector<std::shared_ptr<Shell>>> newBasis, bool isPrimary) {
    /*
     * There is actually no problem if the same basis is tried to be added twice, but it seems
     * quite weird. Thus nothing happens here.
     */
    if (_associatedBasis.find(newBasis.first) == _associatedBasis.end()) {
      _associatedBasis.insert(newBasis);
      if (isPrimary) {
        _primaryBasisLabel = newBasis.first;
      }
    }
  }
  /**
   * @brief Supplies this atom with a set of basis functions and associated ECP functions.
   *
   * Note that adding effective core potential (ECP) functions forces the basis to become
   * the primary basis and no other primary basis may have been specified before! Otherwise
   * one would be able to specify two sets of ECP functions which renders the effective charge
   * of the atom (and other properties) ambiguous. A similar problem arises if the ecp is
   * only added to the atom after its effective charge has been used somewhere.
   *
   * @param newBasis  the primary basis with an identifying label
   * @param ecpSet    A set of effective core potential functions which is used along with
   *                  the newBasis.
   * @param nECPElectrons The number of electrons contained in the ECP.
   */
  void addBasis(std::pair<std::string, std::vector<std::shared_ptr<Shell>>> newBasis,
                std::shared_ptr<libecpint::ECP> ecp, unsigned int nECPElectrons);
  /**
   * @brief Deletes a previously attached set of basis functions.
   *
   * @param label the identifying label used at construction or the addBasis() call.
   */
  inline void eraseBasis(std::string label) {
    _associatedBasis.erase(label);
  }
  /// @returns the name given to the currently active set of basis functions
  inline const std::string& getPrimaryBasisLabel() const {
    return _primaryBasisLabel;
  }
  /**
   * @param newGrid which will be attached to this atom
   * @param isPrimary indicates whether newGrid will be the new 'default' grid for this atom
   */
  void addGrid(std::pair<std::string, AtomGrid*> newGrid, bool isPrimary);
  /**
   * @returns the primary grid attached to this atom
   */
  AtomGrid* getGrid();
  /**
   * @param label a (hopefully) unique identifier
   * @returns the grid identified by the label
   */
  AtomGrid* getGrid(std::string label);

  /**
   * @brief sets the new x coordinate, also moves the basis functions accordingly
   * @param x  new x coordinate
   */
  void setX(double x) override final;

  /**
   * @brief sets the new y coordinate, also moves the basis functions accordingly
   * @param y new y coordinate
   */
  void setY(double y) override final;

  /**
   * @brief sets the new z coordinate, also moves the basis functions accordingly
   * @param z new z coordinate
   */
  void setZ(double z) override final;

  void addToX(double add_to_x) override final;
  void addToY(double add_to_y) override final;
  void addToZ(double add_to_z) override final;

  /**
   * @brief Returns the Cartesian gradients for this atom.
   * @returns Cartesian gradients
   */
  const Eigen::Vector3d& getGradient() const;

  inline bool gradientsUpToDate() {
    return _gradientsUpToDate;
  }
  /**
   * @returns true if the atom is defined with an effective core potential
   * and thus getEffectiveCharge() != getNuclearCharge(), false otherwise
   */
  inline bool usesECP() const {
    assert(_primaryBasisLabel != "-");
    return _nECPElectrons != 0;
  }
  /**
   * @returns the number of core electrons. These core electrons are purely described by an
   *          effective potential / pseudopotential (ECP) and thus appear as being nonexistent
   *          in most other contexts (e.g. they do not enter the electron density)
   */
  inline unsigned int getNECPElectrons() const {
    assert(_primaryBasisLabel != "-");
    return _nECPElectrons;
  }
  /**
   * @brief Getter for the number of non-valence electrons. Respects ECPs.
   * @return The number of non-valence/core orbitals.
   */
  inline unsigned int getNCoreElectrons() const {
    int nCore = _atomType->getNCoreElectrons() - _nECPElectrons;
    return (nCore < 0) ? 0 : nCore;
  }

  /**
   * @brief Getter for the chemical hardness.
   * @return Chemical hardness.
   */
  inline double getChemicalHardness() const {
    return _atomType->getChemicalHardness();
  }

  std::shared_ptr<libecpint::ECP> getCorePotential() const;

 private:
  const std::shared_ptr<const AtomType> _atomType;

  bool _gradientsUpToDate;
  Eigen::Vector3d _gradient;
  /*
   * Most often basis functions are centered on atoms.
   * For Mulliken-type analyses an actual relationship
   * between an atom and basis functions which are centered
   * there is formed. Scince multiple basis sets can in
   * principle be present at the same time a map of vectors
   * is used here (the map connects a label to the basis
   * set for this atom).
   */
  std::map<std::string, std::vector<std::shared_ptr<Shell>>> _associatedBasis;
  /*
   * The set of functions representing the core potential for this atom
   */
  std::shared_ptr<libecpint::ECP> _corePotential;
  /*
   * The label of the active basis.
   */
  std::string _primaryBasisLabel;
  // TODO think about ownership and use the corresponding smart pointer. Probably unique.
  std::map<std::string, AtomGrid*> _grids;

  std::string _primaryGridLabel;
  /*
   * The number of core electrons in case the atom is defined with an effective core potential
   */
  unsigned int _nECPElectrons = 0;
};

} /* namespace Serenity */
#endif /* ATOM_H_ */
