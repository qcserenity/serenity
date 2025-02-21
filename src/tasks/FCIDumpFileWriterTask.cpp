/**
 * @file FCIDumpFileWriterTask.cpp
 *
 * @author Moritz Bensberg
 * @date Feb. 12, 2024
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

/* Include Class Header*/
#include "tasks/FCIDumpFileWriterTask.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"   //Minimal basis access to set the ranges for the valence orbitals.
#include "data/OrbitalController.h"  //Valence orbital ranges.
#include "io/FCIDumpFileWriter.h"    //Write FCIDump file
#include "settings/Settings.h"       //Check the electronic structure method.
#include "system/SystemController.h" //System properties.
/* Include Std and External Headers */
#include <iostream> //Output file

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FCIDumpFileWriterTask<SCFMode>::FCIDumpFileWriterTask(std::shared_ptr<SystemController> activeSystem,
                                                      std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _activeSystem(activeSystem), _environmentSystems(environmentSystems) {
  if (activeSystem->getSettings().method != Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    throw SerenityError("The active system must be described with HF!");
  }
}

template<Options::SCF_MODES SCFMode>
void FCIDumpFileWriterTask<SCFMode>::run() {
  FCIDumpFileWriter<SCFMode> writer(_activeSystem, _environmentSystems, settings);

  auto orbitalRanges = this->getOrbitalRanges();
  std::ofstream outputFile(settings.outputFilePath);
  writer.writeFCIDumpFile(outputFile, orbitalRanges);
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, std::vector<unsigned int>> FCIDumpFileWriterTask<SCFMode>::getOrbitalRanges() {
  const auto nOcc = _activeSystem->getNOccupiedOrbitals<SCFMode>();
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> ranges;
  if (!settings.orbitalRangeAlpha.empty() or !settings.orbitalRangeBeta.empty()) {
    std::vector<std::vector<unsigned int>> bothRanges = {settings.orbitalRangeAlpha, settings.orbitalRangeBeta};
    unsigned int scfModeCounter = 0;
    for_spin(ranges) {
      ranges_spin = bothRanges[scfModeCounter];
      scfModeCounter++;
    };
    return ranges;
  }
  if (settings.onlyValenceOrbitals) {
    if (settings.valenceOrbitalsFromEnergyCutOff) {
      this->_activeSystem->template getActiveOrbitalController<SCFMode>()->setCoreOrbitalsByEnergyCutOff(settings.energyCutOff);
      this->_activeSystem->template getActiveOrbitalController<SCFMode>()->setRydbergOrbitalsByEnergyCutOff(
          settings.virtualEnergyCutOff);
    }
    else {
      const int nMinimalBasisFunctions =
          this->_activeSystem->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION)->getNBasisFunctions();
      const int nBasisFunctions = this->_activeSystem->getBasisController()->getNBasisFunctions();
      const unsigned int nRydbergOrbitals = std::max(0, nBasisFunctions - nMinimalBasisFunctions);
      this->_activeSystem->template getActiveOrbitalController<SCFMode>()->setRydbergOrbitalsByNumber(nRydbergOrbitals);
    }

    auto occRange = this->_activeSystem->template getActiveOrbitalController<SCFMode>()->getValenceOrbitalIndices(nOcc).first;
    auto virtRange =
        this->_activeSystem->template getActiveOrbitalController<SCFMode>()->getVirtualValenceOrbitalIndices(nOcc).first;
    for_spin(occRange, virtRange, ranges) {
      ranges_spin = occRange_spin;
      ranges_spin.insert(ranges_spin.end(), virtRange_spin.begin(), virtRange_spin.end());
    };
    return ranges;
  }
  const auto coefficients = _activeSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  for_spin(ranges, coefficients) {
    std::vector<unsigned int> allOrbitals;
    for (unsigned int iOrb = 0; iOrb < coefficients_spin.cols(); ++iOrb) {
      if (coefficients_spin.col(iOrb).norm() > 1e-1) {
        allOrbitals.push_back(iOrb);
      }
    }
    ranges_spin = allOrbitals;
  };
  return ranges;
}

template class FCIDumpFileWriterTask<Options::SCF_MODES::RESTRICTED>;
template class FCIDumpFileWriterTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
