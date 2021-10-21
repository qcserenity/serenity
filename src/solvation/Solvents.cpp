/**
 * @file Solvents.cpp
 *
 * @author Moritz Bensberg
 * @date May 19, 2020
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
#include "solvation/Solvents.h"
/* Include Serenity Internal Headers */
#include "geometry/AtomTypeFactory.h"
#include "misc/SerenityError.h"
#include "parameters/Constants.h"
#include "settings/PCMSettings.h"
/* Include Std and External Headers */
#include <string>

namespace Serenity {

void Solvents::printSolventInfo(const PCMSettings& settings) {
  std::string model;
  auto mod = settings.solverType;
  Options::resolve<Options::PCM_SOLVER_TYPES>(model, mod);
  printf("%4s Solvation Model:       %15s\n", "", model.c_str());
  std::string solvent;
  auto solv = settings.solvent;
  Options::resolve<Options::PCM_SOLVENTS>(solvent, solv);
  printf("%4s PCM Solvent:           %15s\n", "", solvent.c_str());
  double eps = settings.eps;
  if (settings.solvent != Options::PCM_SOLVENTS::EXPLICIT)
    eps = Solvents::getStaticPermittivity(settings.solvent);
  printf("%4s Static Permittivity:             %.2f\n", "", eps);
  double probeRadius = settings.probeRadius;
  if (settings.solvent != Options::PCM_SOLVENTS::EXPLICIT)
    probeRadius = Solvents::getProbeRadius(settings.solvent);
  printf("%4s Solvent Radius:                  %.3f\n", "", probeRadius);
}

double Solvents::getProbeRadius(Options::PCM_SOLVENTS solvent) {
  double radius = 0;
  switch (solvent) {
    case Options::PCM_SOLVENTS::WATER:
      radius = 1.385;
      break;
    case Options::PCM_SOLVENTS::PROPYLENE_CARBONATE:
      radius = 1.385;
      break;
    case Options::PCM_SOLVENTS::DMSO:
      radius = 2.455;
      break;
    case Options::PCM_SOLVENTS::NITROMETHANE:
      radius = 2.155;
      break;
    case Options::PCM_SOLVENTS::ACETONITRILE:
      return 2.155;
      break;
    case Options::PCM_SOLVENTS::METHANOL:
      radius = 1.855;
      break;
    case Options::PCM_SOLVENTS::ETHANOL:
      radius = 2.180;
      break;
    case Options::PCM_SOLVENTS::ACETONE:
      radius = 2.38;
      break;
    case Options::PCM_SOLVENTS::DICHLORETHANE:
      radius = 2.505;
      break;
    case Options::PCM_SOLVENTS::METHYLENECHLORIDE:
      radius = 2.27;
      break;
    case Options::PCM_SOLVENTS::THF:
      radius = 2.9;
      break;
    case Options::PCM_SOLVENTS::ANILINE:
      radius = 2.80;
      break;
    case Options::PCM_SOLVENTS::CHLOROBENZENE:
      radius = 2.805;
      break;
    case Options::PCM_SOLVENTS::CHLOROFORM:
      radius = 2.48;
      break;
    case Options::PCM_SOLVENTS::TOLUENE:
      radius = 2.82;
      break;
    case Options::PCM_SOLVENTS::DIOXANE:
      radius = 2.630;
      break;
    case Options::PCM_SOLVENTS::BENZENE:
      radius = 2.630;
      break;
    case Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE:
      radius = 2.685;
      break;
    case Options::PCM_SOLVENTS::CYCLOHEXANE:
      radius = 2.815;
      break;
    case Options::PCM_SOLVENTS::N_HEPTANE:
      radius = 3.125;
      break;
    default:
      throw SerenityError("Solvent data not tabulated.");
  }
  return radius * ANGSTROM_TO_BOHR;
}
double Solvents::getStaticPermittivity(Options::PCM_SOLVENTS solvent) {
  switch (solvent) {
    case Options::PCM_SOLVENTS::WATER:
      return 78.39;
    case Options::PCM_SOLVENTS::PROPYLENE_CARBONATE:
      return 64.96;
    case Options::PCM_SOLVENTS::DMSO:
      return 46.7;
    case Options::PCM_SOLVENTS::NITROMETHANE:
      return 38.20;
    case Options::PCM_SOLVENTS::ACETONITRILE:
      return 36.64;
    case Options::PCM_SOLVENTS::METHANOL:
      return 32.63;
    case Options::PCM_SOLVENTS::ETHANOL:
      return 24.55;
    case Options::PCM_SOLVENTS::ACETONE:
      return 20.7;
    case Options::PCM_SOLVENTS::DICHLORETHANE:
      return 10.36;
    case Options::PCM_SOLVENTS::METHYLENECHLORIDE:
      return 8.93;
    case Options::PCM_SOLVENTS::THF:
      return 7.58;
    case Options::PCM_SOLVENTS::ANILINE:
      return 6.89;
    case Options::PCM_SOLVENTS::CHLOROBENZENE:
      return 5.621;
    case Options::PCM_SOLVENTS::CHLOROFORM:
      return 4.90;
    case Options::PCM_SOLVENTS::TOLUENE:
      return 2.379;
    case Options::PCM_SOLVENTS::DIOXANE:
      return 2.250;
    case Options::PCM_SOLVENTS::BENZENE:
      return 2.247;
    case Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE:
      return 2.228;
    case Options::PCM_SOLVENTS::CYCLOHEXANE:
      return 2.023;
    case Options::PCM_SOLVENTS::N_HEPTANE:
      return 1.92;
    default:
      throw SerenityError("Solvent data not tabulated.");
  }
  return 78.39;
}

double Solvents::getCavityFormationProbeRadius(Options::PCM_SOLVENTS solvent) {
  switch (solvent) {
    case Options::PCM_SOLVENTS::BENZENE:
      return 0.5 * 5.24 * ANGSTROM_TO_BOHR; // Pierotti, R. A.J. Phys. Chem.1963,67, 1840.
    case Options::PCM_SOLVENTS::CYCLOHEXANE:
      return 0.5 * 5.59 * ANGSTROM_TO_BOHR; // Pierotti, R. A.J. Phys. Chem.1963,67, 1840.
    case Options::PCM_SOLVENTS::WATER:
      return 0.5 * 2.75 * ANGSTROM_TO_BOHR; // Pierotti, R. A.J. Phys. Chem.1965,69, 281.
    case Options::PCM_SOLVENTS::PROPYLENE_CARBONATE:
    case Options::PCM_SOLVENTS::DMSO:
    case Options::PCM_SOLVENTS::NITROMETHANE:
    case Options::PCM_SOLVENTS::ACETONITRILE:
    case Options::PCM_SOLVENTS::METHANOL:
    case Options::PCM_SOLVENTS::ETHANOL:
    case Options::PCM_SOLVENTS::ACETONE:
    case Options::PCM_SOLVENTS::DICHLORETHANE:
    case Options::PCM_SOLVENTS::METHYLENECHLORIDE:
    case Options::PCM_SOLVENTS::THF:
    case Options::PCM_SOLVENTS::ANILINE:
    case Options::PCM_SOLVENTS::CHLOROBENZENE:
    case Options::PCM_SOLVENTS::CHLOROFORM:
    case Options::PCM_SOLVENTS::TOLUENE:
    case Options::PCM_SOLVENTS::DIOXANE:
    case Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE:
    case Options::PCM_SOLVENTS::N_HEPTANE:
      throw SerenityError("Solvent data not tabulated.");
    default:
      throw SerenityError("Solvent data not tabulated.");
  }
  return 1.385;
}

double Solvents::getMolarVolume(Options::PCM_SOLVENTS solvent) {
  return getMolecularWeight(solvent) / getDensity(solvent);
}

double Solvents::getNumberDensity(Options::PCM_SOLVENTS solvent) {
  return 1.0 / getMolarVolume(solvent);
}

double Solvents::getMolecularWeight(Options::PCM_SOLVENTS solvent) {
  return resolveMass(getAtomTypes(solvent)) * U_TO_ELEC_MASS;
}

double Solvents::resolveMass(std::vector<std::pair<unsigned int, std::shared_ptr<const AtomType>>> types) {
  double mass = 0.0;
  for (const auto& type : types) {
    mass += (double)type.first * type.second->getMass();
  }
  return mass;
}

double Solvents::getDensity(Options::PCM_SOLVENTS solvent) {
  // TODO: Get same more reliable references for the data here.
  const double kgPerm3_to_mePerBohr3 = BOHR * BOHR * BOHR / ELEC_MASS;
  switch (solvent) {
    case Options::PCM_SOLVENTS::WATER: {
      return 0.997295 * 1000 * kgPerm3_to_mePerBohr3; // Wikipedia, 25 celsius
    }
    case Options::PCM_SOLVENTS::CYCLOHEXANE:
      return 0.78 * 1000 * kgPerm3_to_mePerBohr3; // Wikipedia
    case Options::PCM_SOLVENTS::PROPYLENE_CARBONATE:
    case Options::PCM_SOLVENTS::DMSO:
    case Options::PCM_SOLVENTS::NITROMETHANE:
    case Options::PCM_SOLVENTS::ACETONITRILE:
    case Options::PCM_SOLVENTS::METHANOL:
    case Options::PCM_SOLVENTS::ETHANOL:
    case Options::PCM_SOLVENTS::ACETONE:
    case Options::PCM_SOLVENTS::DICHLORETHANE:
    case Options::PCM_SOLVENTS::METHYLENECHLORIDE:
    case Options::PCM_SOLVENTS::THF:
    case Options::PCM_SOLVENTS::ANILINE:
    case Options::PCM_SOLVENTS::CHLOROBENZENE:
    case Options::PCM_SOLVENTS::CHLOROFORM:
    case Options::PCM_SOLVENTS::TOLUENE:
    case Options::PCM_SOLVENTS::DIOXANE:
    case Options::PCM_SOLVENTS::BENZENE:
    case Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE:
    case Options::PCM_SOLVENTS::N_HEPTANE:
      throw SerenityError("Solvent data not tabulated. Density");
    default:
      throw SerenityError("Solvent data not tabulated. Density");
  }
}

std::vector<std::pair<unsigned int, std::shared_ptr<const AtomType>>> Solvents::getAtomTypes(Options::PCM_SOLVENTS solvent) {
  auto H = AtomTypeFactory::getAtomType("H");
  auto C = AtomTypeFactory::getAtomType("C");
  auto O = AtomTypeFactory::getAtomType("O");
  auto S = AtomTypeFactory::getAtomType("S");
  auto N = AtomTypeFactory::getAtomType("N");
  auto Cl = AtomTypeFactory::getAtomType("Cl");
  std::vector<std::pair<unsigned int, std::shared_ptr<const AtomType>>> types;
  switch (solvent) {
    case Options::PCM_SOLVENTS::WATER: {
      types.push_back(std::make_pair(2, H));
      types.push_back(std::make_pair(1, O));
      break;
    }
    case Options::PCM_SOLVENTS::PROPYLENE_CARBONATE:
      types.push_back(std::make_pair(4, C));
      types.push_back(std::make_pair(6, H));
      types.push_back(std::make_pair(3, O));
      break;
    case Options::PCM_SOLVENTS::DMSO:
      types.push_back(std::make_pair(2, C));
      types.push_back(std::make_pair(6, H));
      types.push_back(std::make_pair(1, O));
      types.push_back(std::make_pair(1, S));
      break;
    case Options::PCM_SOLVENTS::NITROMETHANE:
      types.push_back(std::make_pair(1, C));
      types.push_back(std::make_pair(3, H));
      types.push_back(std::make_pair(1, N));
      types.push_back(std::make_pair(2, O));
      break;
    case Options::PCM_SOLVENTS::ACETONITRILE:
      types.push_back(std::make_pair(2, C));
      types.push_back(std::make_pair(3, H));
      types.push_back(std::make_pair(1, N));
      break;
    case Options::PCM_SOLVENTS::METHANOL:
      types.push_back(std::make_pair(1, C));
      types.push_back(std::make_pair(4, H));
      types.push_back(std::make_pair(1, O));
      break;
    case Options::PCM_SOLVENTS::ETHANOL:
      types.push_back(std::make_pair(2, C));
      types.push_back(std::make_pair(6, H));
      types.push_back(std::make_pair(1, O));
      break;
    case Options::PCM_SOLVENTS::ACETONE:
      types.push_back(std::make_pair(3, C));
      types.push_back(std::make_pair(6, H));
      types.push_back(std::make_pair(1, O));
      break;
    case Options::PCM_SOLVENTS::DICHLORETHANE:
      types.push_back(std::make_pair(2, C));
      types.push_back(std::make_pair(4, H));
      types.push_back(std::make_pair(2, Cl));
      break;
    case Options::PCM_SOLVENTS::METHYLENECHLORIDE:
      types.push_back(std::make_pair(1, C));
      types.push_back(std::make_pair(2, H));
      types.push_back(std::make_pair(2, Cl));
      break;
    case Options::PCM_SOLVENTS::THF:
      types.push_back(std::make_pair(4, C));
      types.push_back(std::make_pair(8, H));
      types.push_back(std::make_pair(1, O));
      break;
    case Options::PCM_SOLVENTS::ANILINE:
      types.push_back(std::make_pair(6, C));
      types.push_back(std::make_pair(7, H));
      types.push_back(std::make_pair(1, N));
      break;
    case Options::PCM_SOLVENTS::CHLOROBENZENE:
      types.push_back(std::make_pair(6, C));
      types.push_back(std::make_pair(5, H));
      types.push_back(std::make_pair(1, Cl));
      break;
    case Options::PCM_SOLVENTS::CHLOROFORM:
      types.push_back(std::make_pair(1, C));
      types.push_back(std::make_pair(1, H));
      types.push_back(std::make_pair(3, Cl));
      break;
    case Options::PCM_SOLVENTS::TOLUENE:
      types.push_back(std::make_pair(7, C));
      types.push_back(std::make_pair(8, H));
      break;
    case Options::PCM_SOLVENTS::DIOXANE:
      types.push_back(std::make_pair(4, C));
      types.push_back(std::make_pair(8, H));
      types.push_back(std::make_pair(2, O));
      break;
    case Options::PCM_SOLVENTS::BENZENE:
      types.push_back(std::make_pair(6, C));
      types.push_back(std::make_pair(6, H));
      break;
    case Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE:
      types.push_back(std::make_pair(1, C));
      types.push_back(std::make_pair(4, Cl));
      break;
    case Options::PCM_SOLVENTS::CYCLOHEXANE:
      types.push_back(std::make_pair(6, C));
      types.push_back(std::make_pair(12, H));
      break;
    case Options::PCM_SOLVENTS::N_HEPTANE:
      types.push_back(std::make_pair(7, C));
      types.push_back(std::make_pair(16, H));
      break;
      throw SerenityError("Solvent data not tabulated. AtomTypes.");
    default:
      throw SerenityError("Solvent data not tabulated. AtomTypes.");
  }
  return types;
}

} /* namespace Serenity */
