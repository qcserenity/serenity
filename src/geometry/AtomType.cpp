/**
 * @file   AtomType.cpp
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
/* Include Class Header*/
#include "geometry/AtomType.h"
/* Include Serenity Internal Headers */
#include "misc/WarningTracker.h"
#include "parameters/AtomicParameters.h"
#include "parameters/Constants.h"
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <array>
#include <type_traits>
#include <utility>

namespace Serenity {

AtomType::AtomType(const std::string name, const int nuclearCharge, const double mass, const double braggSlaterRadius,
                   const double vanDerWaalsRadius, double uffRadius, unsigned int nCoreElectrons,
                   const std::vector<std::map<ANGULAR_QUANTUM_NUMBER, unsigned int>>& occupations,
                   double chemicalHardness, const bool isDummy)
  : _name(name),
    _nuclearCharge(isDummy ? 0 : nuclearCharge),
    _psePosition(nuclearCharge),
    _mass(mass),
    _braggSlaterRadius(braggSlaterRadius),
    _vanDerWaalsRadius(vanDerWaalsRadius),
    _uffRadius(uffRadius),
    _nCoreElectrons(isDummy ? 0 : nCoreElectrons),
    _occupations(occupations),
    _chemicalHardness(chemicalHardness),
    _isDummy(isDummy) {
  assert(_psePosition && "Even dummy atoms need to have a position in the PSE, that of the actual atom they mimic.");
}

std::string AtomType::getElementSymbol() const {
  std::string copy = _name;
  if (copy.substr(copy.size() - 1) == ":") {
    copy.pop_back();
  }
  return copy;
}

double AtomType::getBraggSlaterRadius() const {
  if (_braggSlaterRadius < 0.0) {
    WarningTracker::printWarning("Warning: No tabulated Bragg-Slater radius available. Simply guessing 1.5 Angstrom.", true);
    return 1.5 * ANGSTROM_TO_BOHR;
  }
  else {
    return _braggSlaterRadius;
  }
}
/// @returns this atom type's Van der Waals radius
double AtomType::getVanDerWaalsRadius() const {
  if (_vanDerWaalsRadius < 0.0) {
    WarningTracker::printWarning("Warning: No tabulated van der Waals radius available. Simply guessing 2.0 Angstrom.", true);
    return 2.0 * ANGSTROM_TO_BOHR;
  }
  else {
    return _vanDerWaalsRadius;
  }
}

double AtomType::getUFFRadius() const {
  if (_uffRadius < 0.0) {
    WarningTracker::printWarning("Warning: No tabulated UFF radius available. Simply guessing 2.0 Angstrom.", true);
    return 2.0 * ANGSTROM_TO_BOHR;
  }
  else {
    return _uffRadius;
  }
}

unsigned int AtomType::getNCoreElectrons() const {
  return _nCoreElectrons;
}

template<>
SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::vector<double>> getOccupationFactors(const AtomType& atomType) {
  unsigned int nElectronsInIncompleteShells = 0;
  unsigned int nFunctionsInIncompleteShells = 0;
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::vector<double>> occVec;
  for (const auto& period : atomType.getOccupations()) {
    for (const auto& elem : period) {
      const unsigned int nElectronsInShell = elem.second;
      const unsigned int nFunctionsInShell =
          N_SHELL_SPH[static_cast<typename std::underlying_type<ANGULAR_QUANTUM_NUMBER>::type>(elem.first)];
      if (nElectronsInShell < 2 * nFunctionsInShell) {
        nElectronsInIncompleteShells += nElectronsInShell;
        nFunctionsInIncompleteShells += nFunctionsInShell;
      }
      else {
        for (unsigned int i = 0; i < nFunctionsInShell; ++i)
          occVec.push_back(2.0);
      }
    }
  }
  for (unsigned int i = 0; i < nFunctionsInIncompleteShells; ++i)
    occVec.push_back((double)nElectronsInIncompleteShells / nFunctionsInIncompleteShells);
  return occVec;
}

/*
 * The unrestricted version is built to respect Hund's rule. All electrons of incomplete shells
 * are first distributed into the alpha orbitals and only after that (if electrons are still left)
 * into beta orbitals
 */
template<>
SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<double>> getOccupationFactors(const AtomType& atomType) {
  unsigned int nElectronsInIncompleteShells = 0;
  unsigned int nFunctionsInIncompleteShells = 0;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<double>> occVec;
  for (const auto& period : atomType.getOccupations()) {
    for (const auto& elem : period) {
      const unsigned int nElectronsInShell = elem.second;
      const unsigned int nFunctionsInShell =
          N_SHELL_SPH[static_cast<typename std::underlying_type<ANGULAR_QUANTUM_NUMBER>::type>(elem.first)];
      if (nElectronsInShell < 2 * nFunctionsInShell) {
        nElectronsInIncompleteShells += nElectronsInShell;
        nFunctionsInIncompleteShells += nFunctionsInShell;
      }
      else {
        for (unsigned int i = 0; i < nFunctionsInShell; ++i) {
          occVec.alpha.push_back(1.0);
          occVec.beta.push_back(1.0);
        }
      }
    }
  }
  double alphaElectronsPerFunc;
  double betaElectronsPerFunc;
  if (nElectronsInIncompleteShells < nFunctionsInIncompleteShells) {
    alphaElectronsPerFunc = (double)nElectronsInIncompleteShells / nFunctionsInIncompleteShells;
    betaElectronsPerFunc = 0;
  }
  else {
    alphaElectronsPerFunc = 1.0;
    betaElectronsPerFunc = (double)(nElectronsInIncompleteShells - nFunctionsInIncompleteShells) / nFunctionsInIncompleteShells;
  }
  for (unsigned int i = 0; i < nFunctionsInIncompleteShells; ++i) {
    occVec.alpha.push_back(alphaElectronsPerFunc);
    occVec.beta.push_back(betaElectronsPerFunc);
  }
  return occVec;
}

double AtomType::getChemicalHardness() const {
  if (_chemicalHardness < 0.0) {
    WarningTracker::printWarning("Warning: No tabulated chemical hardness available. Simply guessing 0.2 a.u.", true);
    return 0.2;
  }
  else {
    return _chemicalHardness;
  }
}

int getAtomSpin(const AtomType& atomType) {
  unsigned int nElectronsInIncompleteShells = 0;
  unsigned int nFunctionsInIncompleteShells = 0;
  for (const auto& period : atomType.getOccupations()) {
    for (const auto& elem : period) {
      const unsigned int nElectronsInShell = elem.second;
      const unsigned int nFunctionsInShell =
          N_SHELL_SPH[static_cast<typename std::underlying_type<ANGULAR_QUANTUM_NUMBER>::type>(elem.first)];
      assert(nElectronsInShell <= 2 * nFunctionsInShell);
      if (nElectronsInShell < 2 * nFunctionsInShell) {
        nElectronsInIncompleteShells += nElectronsInShell;
        nFunctionsInIncompleteShells += nFunctionsInShell;
      }
    }
  }
  if (nElectronsInIncompleteShells < nFunctionsInIncompleteShells) {
    return nElectronsInIncompleteShells;
  }
  else {
    return 2 * nFunctionsInIncompleteShells - nElectronsInIncompleteShells;
  }
}

unsigned int AtomType::getMinimalBasisSize() const {
  if (_psePosition < 3) {
    return 1;
  }
  else if (_psePosition < 11) {
    return 5;
  }
  else if (_psePosition < 19) {
    return 9;
  }
  else if (_psePosition < 21) {
    return 13;
  }
  else if (_psePosition < 37) {
    return 18;
  }
  else if (_psePosition < 39) {
    return 22;
  }
  else if (_psePosition < 55) {
    return 27;
  }
  else if (_psePosition < 57) {
    return 31;
  }
  else if (_psePosition < 87) {
    return 43;
  }
  else if (_psePosition < 89) {
    return 47;
  }
  else {
    return 59;
  }
}

} /* namespace Serenity */
