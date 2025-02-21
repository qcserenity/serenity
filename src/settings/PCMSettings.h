/**
 * @file PCMSettings.h
 *
 * @date Apr 30, 2020
 * @author Moritz Bensberg
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

#ifndef SETTINGS_PCMSETTINGS_H_
#define SETTINGS_PCMSETTINGS_H_
/* Include Serenity Internal Headers */
#include "settings/PCMOptions.h"
#include "settings/Reflection.h"
/* Include Std and External Headers */
#include <string>

namespace Serenity {
using namespace Serenity::Reflection;

/**
 * @brief Settings used for specifying the PCM model.
 *   --use                If true, the PCM model is employed.
 *   --cavity             The cavity type.
 *   --scaling            If true, the radii for the spheres will be scaled by 1.2.
 *   --radiiType          Atomic radii used for cavity construction.
 *                        Note that not all radii are available for all atoms.
 *   --minRadius          Minimal radius for additional spheres not centered on atoms (GEPOL only).
 *   --solverType         Type of solver to be used.
 *   --solvent            The solvent string the specifies the solvent used.
 *   --correction         Correction k for the apparent surface charge scaling factor in the CPCM solver.
 *   --probeRadius        Radius of the spherical probe approximating a solvent molecule (Overridden by the built-in
 * value for the chosen solvent).
 *   --eps                Static dielectric permittivity of the medium.
 *   --patchLevel         Wavelet cavity mesh patch level (GEPOL only).
 *   --overlapFactor      Maximum ratio of a new sphere to be allowed to be covered within the already present ones.
 *   --cacheSize          Maximum number of two center integrals to be stored in memory.
 *   --minDistance        Minimal distance between sampling points.
 *   --lLarge             Angular momentum used for the spherical grid construction for non-hydrogen atoms (DELLEY
 * only).
 *   --lSmall             Angular momentum used for the spherical grid construction for hydrogen atoms (DELLEY only).
 *   --alpha              Sharpness parameter for molecular surface model function (DELLEY only).
 *   --projectionCutOff   Cut off for the projection to the molecular surface (DELLEY only).
 *   --oneCavity          All surface points have to be connected.
 *   --connectivityFactor Connectivity between surface points is determined by a distance cut off. The distance is
 *                        given by connectivityFactor * probeRadius
 *   --numberDensity      The number density of the solvent (number of particles per volume)
 *   --temperature        The temperature of the ensemble.
 *   --cavityFormation    Calculate the cavity formation energy using the scaled particle theory.
 *   --cavityProbeRadius  The solvent probe radius to be used in the calculation of the cavity formation.
 *   --saveCharges        Switch to determine whether we want to save the PCM charges.
 *
 */
struct PCMSettings {
  PCMSettings()
    : use(false),
      cavity(Options::PCM_CAVITY_TYPES::DELLEY),
      scaling(true),
      radiiType(Options::PCM_ATOMIC_RADII_TYPES::BONDI),
      minRadius(0.377), // 0.2 Angstrom
      solverType(Options::PCM_SOLVER_TYPES::CPCM),
      solvent(Options::PCM_SOLVENTS::WATER),
      correction(0.0),
      probeRadius(1.0),
      eps(1.0),
      patchLevel(0),
      minDistance(0.2),
      overlapFactor(0.7),
      cacheSize(128),
      lLarge(7),
      lSmall(4),
      alpha(50.0),
      projectionCutOff(5.0),
      oneCavity(false),
      connectivityFactor(2.0),
      numberDensity(1.0),
      temperature(298.15),
      cavityFormation(false),
      cavityProbeRadius(5.0),
      saveCharges(false) {
  }
  REFLECTABLE((bool)use, (Options::PCM_CAVITY_TYPES)cavity, (bool)scaling, (Options::PCM_ATOMIC_RADII_TYPES)radiiType,
              (double)minRadius, (Options::PCM_SOLVER_TYPES)solverType, (Options::PCM_SOLVENTS)solvent,
              (double)correction, (double)probeRadius, (double)eps, (int)patchLevel, (double)minDistance,
              (double)overlapFactor, (unsigned int)cacheSize, (unsigned int)lLarge, (unsigned int)lSmall, (double)alpha,
              (double)projectionCutOff, (bool)oneCavity, (double)connectivityFactor, (double)numberDensity,
              (double)temperature, (bool)cavityFormation, (double)cavityProbeRadius, (bool)saveCharges);

 public:
  /**
   * @brief Parse the settings from the visitor to this object.
   * @param v The visitor.
   * @param blockname The block name.
   * @return True if the block name corresponds to this block. Otherwise false.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("PCM")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
  // Switch to determine whether we are dealing with a PCM loaded from another calculation.
  bool loadedPCM = false;
  // The file path of the PCM to be loaded from another calculation.
  std::string cavityPath = "";
};

} /* namespace Serenity */

#endif /* SETTINGS_PCMSETTINGS_H_ */
