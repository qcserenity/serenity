/**
 * @file   AtomType.h
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
#ifndef ATOMTYPE_H_
#define ATOMTYPE_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <cassert>
#include <map>
#include <string>
#include <vector>

namespace Serenity {
/* Forward Declarations */
enum class ANGULAR_QUANTUM_NUMBER;
namespace Options {
enum class SCF_MODES;
}

/**
 * @class AtomType AtomType.h
 *
 * @brief Model class which holds data about an atom type, i.e. element symbol, weight, ...
 *
 * More exactly typically this is only information related to the nucleus.
 */
class AtomType {
 public:
  /**
   * @param name a unique identifier
   * @param nuclearCharge Dummys also need their nuclear charge/PSE number given
   *                      (the resulting nuclear charge will be zero though)
   * @param mass
   * @param braggSlaterRadius is needed for grid setups
   * @param vanDerWaalsRadius The van der Waals radius.
   * @param uffRadius The UFF radius.
   * @param nCoreElectrons The number of non-valence electrons.
   * @param occupations Isolated atom occupations.
   * @param chemicalHardness Chemical hardness.
   * @param isDummy Flag for dummy atom types.
   */
  AtomType(std::string name, int nuclearCharge, double mass, double braggSlaterRadius, double vanDerWaalsRadius, double uffRadius,
           unsigned int nCoreElectrons, const std::vector<std::map<ANGULAR_QUANTUM_NUMBER, unsigned int>>& occupations,
           double chemicalHardness, const bool isDummy = false);
  virtual ~AtomType() = default;
  /// @returns the atom's mass, i.e. the mass of the selected isotope (no average), in u (g/mol)
  inline double getMass() const {
    return _mass;
  }
  /// @returns a unique identifier
  inline const std::string& getName() const {
    return _name;
  }
  /// @returns the charge of the nucleus, NOT the charge of an atom/ion!
  inline int getNuclearCharge() const {
    return _nuclearCharge;
  }
  /// @returns The position in the PSE (a version of getNuclearCharge() that also works for dummy atoms).
  inline int getPSEPosition() const {
    return _psePosition;
  }
  /// @returns the element symbol
  std::string getElementSymbol() const;
  /// @returns this atom type's Bragg-Slater radius, e.g. used for setting up integration grids
  double getBraggSlaterRadius() const;
  /// @returns this atom type's Van der Waals radius
  double getVanDerWaalsRadius() const;
  /// @returns this atom type's UFF radius
  double getUFFRadius() const;
  /// @returns The number of core electrons.
  unsigned int getNCoreElectrons() const;
  /// @returns The chemical hardness of this atom type.
  double getChemicalHardness() const;

  /**
   * @brief Dummy check.
   * @return Returns true if the atom type is a dummy atom.
   */
  inline bool isDummy() const {
    return _isDummy;
  }

  /**
   * @returns A vector resembling the expected occupations of an atom, i.e. how many s, p, d...-
   *          electrons it has as an isolated atom in the ground state. Each entry in the vector
   *          holds the occupations within one period.
   */
  inline const std::vector<std::map<ANGULAR_QUANTUM_NUMBER, unsigned int>>& getOccupations() const {
    return _occupations;
  }

  /** @returns in which period the element of the atom type is found */
  inline unsigned int getPeriod() const {
    return _occupations.size();
  }
  inline unsigned int getRow() const {
    assert(_psePosition && "Even dummy atoms need to have a position in the PSE, that of the actual atom they mimic.");
    if (_psePosition < 3) {
      return 1;
    }
    else if (_psePosition < 11) {
      return 2;
    }
    else if (_psePosition < 19) {
      return 3;
    }
    else if (_psePosition < 37) {
      return 4;
    }
    else if (_psePosition < 55) {
      return 5;
    }
    else if (_psePosition < 88) {
      return 6;
    }
    else {
      return 7;
    }
  }
  /**
   * @brief Getter for the minimal basis size of the given element.
   * @return The minimal basis size.
   */
  unsigned int getMinimalBasisSize() const;

 private:
  const std::string _name;
  const int _nuclearCharge;
  const int _psePosition;
  const double _mass;
  const double _braggSlaterRadius;
  const double _vanDerWaalsRadius;
  const double _uffRadius;
  const unsigned int _nCoreElectrons;
  const std::vector<std::map<ANGULAR_QUANTUM_NUMBER, unsigned int>> _occupations;
  const double _chemicalHardness;
  const bool _isDummy;
};

/**
 * @param   atomType
 * @returns for each spin a vector indicating how orbitals are filled in the neutral atom
 *          resulting in a spherical distribution of electrons. Spins are according to Hund's rule.
 */
template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, std::vector<double>> getOccupationFactors(const AtomType& atomType);
/**
 * @param   atomType
 * @returns the spin of an atom of atomType in its ground state electronic configuration, based on
 *          Hund's rule (always positive or zero, i.e. excess of alpha electrons).
 */
int getAtomSpin(const AtomType& atomType);

} /* namespace Serenity */
#endif /* ATOMTYPE_H_ */
