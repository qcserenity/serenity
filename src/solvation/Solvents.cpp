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
  const std::map<Options::PCM_SOLVENTS, double> radii = {{Options::PCM_SOLVENTS::WATER, 1.385},
                                                         {Options::PCM_SOLVENTS::PROPYLENE_CARBONATE, 1.385},
                                                         {Options::PCM_SOLVENTS::DMSO, 2.455},
                                                         {Options::PCM_SOLVENTS::NITROMETHANE, 2.155},
                                                         {Options::PCM_SOLVENTS::ACETONITRILE, 2.155},
                                                         {Options::PCM_SOLVENTS::METHANOL, 1.855},
                                                         {Options::PCM_SOLVENTS::ETHANOL, 2.180},
                                                         {Options::PCM_SOLVENTS::ACETONE, 2.38},
                                                         {Options::PCM_SOLVENTS::DICHLORETHANE, 2.505},
                                                         {Options::PCM_SOLVENTS::METHYLENECHLORIDE, 2.27},
                                                         {Options::PCM_SOLVENTS::THF, 2.9},
                                                         {Options::PCM_SOLVENTS::ANILINE, 2.80},
                                                         {Options::PCM_SOLVENTS::CHLOROBENZENE, 2.805},
                                                         {Options::PCM_SOLVENTS::CHLOROFORM, 2.48},
                                                         {Options::PCM_SOLVENTS::TOLUENE, 2.82},
                                                         {Options::PCM_SOLVENTS::DIOXANE, 2.630},
                                                         {Options::PCM_SOLVENTS::BENZENE, 2.630},
                                                         {Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE, 2.685},
                                                         {Options::PCM_SOLVENTS::CYCLOHEXANE, 2.815},
                                                         {Options::PCM_SOLVENTS::N_HEPTANE, 3.125}};
  if (radii.find(solvent) == radii.end())
    throw SerenityError("Solvent data not tabulated.");
  return radii.at(solvent) * ANGSTROM_TO_BOHR;
}

double Solvents::getStaticPermittivity(Options::PCM_SOLVENTS solvent) {
  const std::map<Options::PCM_SOLVENTS, double> eps = {{Options::PCM_SOLVENTS::WATER, 78.39},
                                                       {Options::PCM_SOLVENTS::PROPYLENE_CARBONATE, 64.96},
                                                       {Options::PCM_SOLVENTS::DMSO, 46.7},
                                                       {Options::PCM_SOLVENTS::NITROMETHANE, 38.20},
                                                       {Options::PCM_SOLVENTS::ACETONITRILE, 36.64},
                                                       {Options::PCM_SOLVENTS::METHANOL, 32.63},
                                                       {Options::PCM_SOLVENTS::ETHANOL, 24.55},
                                                       {Options::PCM_SOLVENTS::ACETONE, 20.7},
                                                       {Options::PCM_SOLVENTS::DICHLORETHANE, 10.36},
                                                       {Options::PCM_SOLVENTS::METHYLENECHLORIDE, 8.93},
                                                       {Options::PCM_SOLVENTS::THF, 7.58},
                                                       {Options::PCM_SOLVENTS::ANILINE, 6.89},
                                                       {Options::PCM_SOLVENTS::CHLOROBENZENE, 5.621},
                                                       {Options::PCM_SOLVENTS::CHLOROFORM, 4.90},
                                                       {Options::PCM_SOLVENTS::TOLUENE, 2.379},
                                                       {Options::PCM_SOLVENTS::DIOXANE, 2.250},
                                                       {Options::PCM_SOLVENTS::BENZENE, 2.247},
                                                       {Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE, 2.228},
                                                       {Options::PCM_SOLVENTS::CYCLOHEXANE, 2.023},
                                                       {Options::PCM_SOLVENTS::N_HEPTANE, 1.92}};
  if (eps.find(solvent) == eps.end())
    throw SerenityError("Solvent data not tabulated.");
  return eps.at(solvent);
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
    default:
      throw SerenityError("Solvent data not tabulated. AtomTypes.");
  }
  return types;
}

} /* namespace Serenity */
