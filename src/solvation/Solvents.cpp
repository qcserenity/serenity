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
#include "misc/SerenityError.h"
#include "settings/Options.h"
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
  switch (solvent) {
    case Options::PCM_SOLVENTS::WATER:
      return 1.385;
    case Options::PCM_SOLVENTS::PROPYLENE_CARBONATE:
      return 1.385;
    case Options::PCM_SOLVENTS::DMSO:
      return 2.455;
    case Options::PCM_SOLVENTS::NITROMETHANE:
      return 2.155;
    case Options::PCM_SOLVENTS::ACETONITRILE:
      return 2.155;
    case Options::PCM_SOLVENTS::METHANOL:
      return 1.855;
    case Options::PCM_SOLVENTS::ETHANOL:
      return 2.180;
    case Options::PCM_SOLVENTS::ACETONE:
      return 2.38;
    case Options::PCM_SOLVENTS::DICHLORETHANE:
      return 2.505;
    case Options::PCM_SOLVENTS::METHYLENECHLORIDE:
      return 2.27;
    case Options::PCM_SOLVENTS::THF:
      return 2.9;
    case Options::PCM_SOLVENTS::ANILINE:
      return 2.80;
    case Options::PCM_SOLVENTS::CHLOROBENZENE:
      return 2.805;
    case Options::PCM_SOLVENTS::CHLOROFORM:
      return 2.48;
    case Options::PCM_SOLVENTS::TOLUENE:
      return 2.82;
    case Options::PCM_SOLVENTS::DIOXANE:
      return 2.630;
    case Options::PCM_SOLVENTS::BENZENE:
      return 2.630;
    case Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE:
      return 2.685;
    case Options::PCM_SOLVENTS::CYCLOHEXANE:
      return 2.815;
    case Options::PCM_SOLVENTS::N_HEPTANE:
      return 3.125;
    default:
      throw SerenityError("Solvent data not tabulated.");
  }
  return 1.385;
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

} /* namespace Serenity */
