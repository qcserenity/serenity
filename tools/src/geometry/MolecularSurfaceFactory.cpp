/**
 * @file MolecularSurfaceFactory.cpp
 *
 * @author Moritz Bensberg
 * @date May 25, 2020
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
#include "geometry/MolecularSurfaceFactory.h"
/* Include Serenity Internal Headers */
#include "geometry/DelleySurfaceConstructor.h"
#include "geometry/GEPOLSurfaceConstructor.h"
#include "geometry/Geometry.h"
#include "geometry/MolecularSurface.h"
#include "geometry/Sphere.h"
#include "misc/Timing.h" //Timings.
#include "settings/PCMSettings.h"
#include "solvation/Solvents.h"

namespace Serenity {

std::unique_ptr<MolecularSurface>
MolecularSurfaceFactory::produce(std::shared_ptr<const Geometry> geometry, Options::PCM_CAVITY_TYPES cavityType,
                                 Options::PCM_ATOMIC_RADII_TYPES radiiType, unsigned int calveLevel, double minDistance,
                                 double minRadius, double solvRadius, double overlapFactor, bool scaling,
                                 unsigned int sphericalAngularMomentum, unsigned int smallSphericalAngularMomentum,
                                 double alpha, double projectionCutOff, bool oneCavity, double connectivityFactor) {
  Timings::takeTime(" Tech. -    Cavity Construction");
  auto atoms = geometry->getAtoms();
  std::vector<Sphere> initialSpheres;
  std::unique_ptr<MolecularSurface> molecularSurface;
  for (auto atom : atoms) {
    unsigned int l = sphericalAngularMomentum;
    if (atom->getAtomType()->getElementSymbol() == "H")
      l = smallSphericalAngularMomentum;
    initialSpheres.push_back(Sphere(*atom, getAtomRadius(atom, radiiType, scaling), SphereType::VanDerWaals, l));
  } // for atom
  switch (cavityType) {
    case Options::PCM_CAVITY_TYPES::GEPOL_SAS: {
      GEPOLSurfaceConstructor gepol(initialSpheres, true, calveLevel, minDistance, minRadius, solvRadius, overlapFactor);
      molecularSurface = gepol.getMolecularSurface();
      break;
    }
    case Options::PCM_CAVITY_TYPES::GEPOL_SES: {
      GEPOLSurfaceConstructor gepol(initialSpheres, false, calveLevel, minDistance, minRadius, solvRadius, overlapFactor);
      molecularSurface = gepol.getMolecularSurface();
      break;
    }
    case Options::PCM_CAVITY_TYPES::DELLEY: {
      DelleySurfaceConstructor delley(initialSpheres, solvRadius, alpha, projectionCutOff, minDistance, oneCavity,
                                      connectivityFactor);
      molecularSurface = delley.getMolecularSurface();
      break;
    }
  }
  Timings::timeTaken(" Tech. -    Cavity Construction");
  return molecularSurface;
}

std::unique_ptr<MolecularSurface> MolecularSurfaceFactory::produce(std::shared_ptr<const Geometry> geometry,
                                                                   const PCMSettings& pcmSettings) {
  double solvRad = pcmSettings.probeRadius;
  if (pcmSettings.solvent != Options::PCM_SOLVENTS::EXPLICIT)
    solvRad = Solvents::getProbeRadius(pcmSettings.solvent);
  return produce(geometry, pcmSettings.cavity, pcmSettings.radiiType, pcmSettings.patchLevel, pcmSettings.minDistance,
                 pcmSettings.minRadius, solvRad, pcmSettings.overlapFactor, pcmSettings.scaling, pcmSettings.lLarge,
                 pcmSettings.lSmall, pcmSettings.alpha, pcmSettings.projectionCutOff, pcmSettings.oneCavity,
                 pcmSettings.connectivityFactor);
}

double MolecularSurfaceFactory::getAtomRadius(std::shared_ptr<Atom> atom, Options::PCM_ATOMIC_RADII_TYPES radiiType,
                                              bool scaling) {
  double radius = -1.0;
  switch (radiiType) {
    case Options::PCM_ATOMIC_RADII_TYPES::BONDI:
      radius = atom->getAtomType()->getVanDerWaalsRadius();
      break;
    case Options::PCM_ATOMIC_RADII_TYPES::UFF:
      radius = atom->getAtomType()->getUFFRadius();
      break;
  } // switch radiiType
  if (scaling)
    radius *= 1.2;
  return radius;
}

} /* namespace Serenity */
