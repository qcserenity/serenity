/**
 * @file   LocalizationTask.cpp
 *
 * @date   Apr 22, 2014
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
#include "tasks/LocalizationTask.h"
/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/EdmistonRuedenbergLocalization.h"
#include "analysis/orbitalLocalization/FosterBoysLocalization.h"
#include "analysis/orbitalLocalization/IBOLocalization.h"
#include "analysis/orbitalLocalization/NonOrthogonalLocalization.h"
#include "analysis/orbitalLocalization/OrbitalAligner.h" //Orbital alignment.
#include "analysis/orbitalLocalization/PipekMezeyLocalization.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/DensityMatrixController.h"
#include "geometry/Geometry.h" //getNMinimalBasisFunctions()
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h" //Filtered output streams.
#include "math/Matrix.h"
#include "misc/WarningTracker.h" //Warnings.
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <Eigen/StdVector>
#include <algorithm> // std::max
#include <iostream>

namespace Serenity {

LocalizationTask::LocalizationTask(std::shared_ptr<SystemController> systemController,
                                   std::vector<std::shared_ptr<SystemController>> templateSystem)
  : _systemController(systemController), _templateSystem(templateSystem) {
  assert(_systemController);
}

void LocalizationTask::run() {
  if (_systemController->getLastSCFMode() == Options::SCF_MODES::UNRESTRICTED or
      _systemController->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
    runByLastSCFMode<Options::SCF_MODES::UNRESTRICTED>();
  }
  else {
    runByLastSCFMode<Options::SCF_MODES::RESTRICTED>();
  }
}
template<>
std::pair<SpinPolarizedData<RESTRICTED, std::vector<unsigned int>>, SpinPolarizedData<RESTRICTED, std::vector<unsigned int>>>
LocalizationTask::separateSOMOs(const SpinPolarizedData<RESTRICTED, std::vector<unsigned int>>& valenceOrbitalRange) {
  return std::make_pair(valenceOrbitalRange, SpinPolarizedData<RESTRICTED, std::vector<unsigned int>>());
}

template<>
std::pair<SpinPolarizedData<UNRESTRICTED, std::vector<unsigned int>>, SpinPolarizedData<UNRESTRICTED, std::vector<unsigned int>>>
LocalizationTask::separateSOMOs(const SpinPolarizedData<UNRESTRICTED, std::vector<unsigned int>>& valenceOrbitalRange) {
  OutputControl::nOut
      << "  Separating SOMOs from normal orbitals. Note that this assumes that the orbitals are (quasi) restricted."
      << std::endl;
  std::vector<unsigned int> nonSOMOs = (valenceOrbitalRange.alpha.size() > valenceOrbitalRange.beta.size())
                                           ? valenceOrbitalRange.beta
                                           : valenceOrbitalRange.alpha;
  SpinPolarizedData<UNRESTRICTED, std::vector<unsigned int>> doMORange(nonSOMOs);
  SpinPolarizedData<UNRESTRICTED, std::vector<unsigned int>> soMORange;
  for_spin(doMORange, valenceOrbitalRange, soMORange) {
    for (const auto& iOrb : valenceOrbitalRange_spin) {
      if (std::find(doMORange_spin.begin(), doMORange_spin.end(), iOrb) == doMORange_spin.end()) {
        soMORange_spin.push_back(iOrb);
      }
    }
  };
  return std::make_pair(doMORange, soMORange);
}
template<>
SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::string> LocalizationTask::getOutputPraefix() {
  return "  ";
}
template<>
SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::string> LocalizationTask::getOutputPraefix() {
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::string> toReturn;
  toReturn.alpha = "  alpha orbitals";
  toReturn.beta = "  beta orbitals";
  return toReturn;
}
template<Options::SCF_MODES SCFMode>
std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>>
LocalizationTask::getValenceOrbitalIndices() {
  auto nOcc = _systemController->getNOccupiedOrbitals<SCFMode>();
  std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>> ranges;
  SpinPolarizedData<SCFMode, std::vector<unsigned int>>& valenceRange = ranges.first;
  SpinPolarizedData<SCFMode, std::vector<unsigned int>>& coreRange = ranges.second;
  if (!settings.splitValenceAndCore) {
    for_spin(nOcc, valenceRange) {
      for (unsigned int iOcc = 0; iOcc < nOcc_spin; ++iOcc)
        valenceRange_spin.push_back(iOcc);
    };
    return ranges;
  }

  if (settings.useEnergyCutOff) {
    this->_systemController->template getActiveOrbitalController<SCFMode>()->setCoreOrbitalsByEnergyCutOff(settings.energyCutOff);
  }
  else if (settings.nCoreOrbitals != std::numeric_limits<unsigned int>::infinity()) {
    this->_systemController->template getActiveOrbitalController<SCFMode>()->setCoreOrbitalsByNumber(settings.nCoreOrbitals);
  }
  ranges = this->_systemController->template getActiveOrbitalController<SCFMode>()->getValenceOrbitalIndices(nOcc);

  OutputControl::mOut << "  Note: Valence and core orbitals are localized separately." << std::endl;
  OutputControl::mOut << "        Please ensure that the initial orbitals were canonical orbitals!\n" << std::endl;
  OutputControl::vOut << "The following orbitals have been identified as valence orbitals:" << std::endl;
  SpinPolarizedData<SCFMode, std::string> praefix = getOutputPraefix<SCFMode>();
  for_spin(valenceRange, praefix, nOcc) {
    OutputControl::vOut << praefix_spin << std::endl << "  ";
    for (const auto& orb : valenceRange_spin)
      OutputControl::vOut << orb << " ";
    OutputControl::vOut << std::endl;
    OutputControl::nOut << "  Number of core-like orbitals: " << nOcc_spin - valenceRange_spin.size() << "  "
                        << praefix_spin << std::endl;
  };
  OutputControl::nOut << "  Core-like orbitals will be localized separately." << std::endl;
  OutputControl::nOut << "  Running orbital localization for valence orbitals." << std::endl;
  return std::make_pair(valenceRange, coreRange);
}
template<Options::SCF_MODES SCFMode>
std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>>
LocalizationTask::getVirtualValenceOrbitalIndices() {
  const auto nOcc = _systemController->getNOccupiedOrbitals<SCFMode>();
  const auto nBasisFunctions = _systemController->getBasisController()->getNBasisFunctions();
  std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>> ranges;
  SpinPolarizedData<SCFMode, std::vector<unsigned int>>& valenceRange = ranges.first;
  SpinPolarizedData<SCFMode, std::vector<unsigned int>>& rydbergRange = ranges.second;
  if (!settings.splitVirtuals and settings.locType != Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO) {
    for_spin(nOcc, valenceRange) {
      for (unsigned int iVirt = nOcc_spin; iVirt < nBasisFunctions; ++iVirt)
        valenceRange_spin.push_back(iVirt);
    };
    return ranges;
  }
  if (settings.locType == Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO ||
      settings.locType == Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN ||
      settings.locType == Options::ORBITAL_LOCALIZATION_ALGORITHMS::IAO) {
    OutputControl::mOut << "  IBO localization is used. Ignoring any specifications for Rydberg orbitals and"
                           " referencing\n  the IBO minimal basis set instead."
                        << std::endl;
    const int nMinimalBasisFunctions =
        this->_systemController->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION)->getNBasisFunctions();
    const unsigned int nRydbergOrbitals = std::max(0, int(nBasisFunctions) - nMinimalBasisFunctions);
    this->_systemController->template getActiveOrbitalController<SCFMode>()->setRydbergOrbitalsByNumber(nRydbergOrbitals);
  }
  else if (settings.useEnergyCutOff) {
    this->_systemController->template getActiveOrbitalController<SCFMode>()->setRydbergOrbitalsByEnergyCutOff(
        settings.virtualEnergyCutOff);
  }
  else if (settings.nRydbergOrbitals != std::numeric_limits<unsigned int>::infinity()) {
    this->_systemController->template getActiveOrbitalController<SCFMode>()->setRydbergOrbitalsByNumber(settings.nRydbergOrbitals);
  }
  else {
    this->_systemController->template getActiveOrbitalController<SCFMode>()->setRydbergOrbitalsByNumber(
        this->_systemController->getGeometry()->getNMinimalBasisFunctions(true));
  }
  ranges = this->_systemController->template getActiveOrbitalController<SCFMode>()->getVirtualValenceOrbitalIndices(nOcc);
  OutputControl::vOut << "The following orbitals have been identified as valence orbitals:" << std::endl;
  SpinPolarizedData<SCFMode, std::string> praefix = getOutputPraefix<SCFMode>();
  for_spin(valenceRange, praefix, rydbergRange) {
    OutputControl::vOut << praefix_spin << std::endl << "  ";
    for (const auto& orb : valenceRange_spin)
      OutputControl::vOut << orb << " ";
    OutputControl::vOut << std::endl;
    OutputControl::nOut << "  Number of Rydberg-like orbitals: " << rydbergRange_spin.size() << "  " << praefix_spin
                        << std::endl;
    OutputControl::nOut << "  Number of virtual valence orbitals: " << valenceRange_spin.size() << "  " << praefix_spin
                        << std::endl;
  };
  OutputControl::nOut << "  Rydberg-like orbitals will be localized separately or not at all." << std::endl;
  return std::make_pair(valenceRange, rydbergRange);
}
template<Options::SCF_MODES SCFMode>
void LocalizationTask::runByLastSCFMode() {
  std::shared_ptr<Localization<SCFMode>> localizationRoutine;

  std::shared_ptr<FockMatrix<SCFMode>> fockMatrixPtr = nullptr;
  if (_systemController->template getElectronicStructure<SCFMode>()->checkFock()) {
    fockMatrixPtr = std::make_shared<FockMatrix<SCFMode>>(
        _systemController->template getElectronicStructure<SCFMode>()->getFockMatrix());
  }

  switch (settings.locType) {
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS:
      localizationRoutine = std::make_shared<FosterBoysLocalization<SCFMode>>(_systemController);
      printSubSectionTitle("Running Foster-Boys Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY:
      localizationRoutine = std::make_shared<PipekMezeyLocalization<SCFMode>>(_systemController);
      printSubSectionTitle("Running Pipek-Mezey Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO:
      localizationRoutine = std::make_shared<IBOLocalization<SCFMode>>(_systemController, false,
                                                                       settings.replaceVirtuals, settings.separateSOMOs);
      printSubSectionTitle("Running IBO Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::IAO:
      // TODO: Does this actually do anything? Projecting on the IAOs and back should not change the orbitals.
      localizationRoutine = std::make_shared<IBOLocalization<SCFMode>>(_systemController, true, settings.replaceVirtuals);
      printSubSectionTitle("Running IAO Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::EDMISTON_RUEDENBERG:
      localizationRoutine = std::make_shared<EdmistonRuedenbergLocalization<SCFMode>>(_systemController);
      printSubSectionTitle("Running Edmiston-Ruedenberg Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::NON_ORTHOGONAL:
      localizationRoutine = std::make_shared<NonOrthogonalLocalization<SCFMode>>(_systemController);
      printSubSectionTitle("Running Non-Orthogonal Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN:
      if (_templateSystem.size() == 0) {
        throw SerenityError("Orbital alignment was requested without any template system.\n"
                            "Please specify a template system in the task input via the env keyword.");
      }
      localizationRoutine =
          std::make_shared<OrbitalAligner<SCFMode>>(_systemController, _templateSystem[0], settings.alignExponent,
                                                    settings.useKineticAlign, settings.replaceVirtuals);
      printSubSectionTitle("Running IBO-Like Orbital Alignment");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::NONE:
      // Do nothing
      printSubSectionTitle("Skip localization...");
      return;
  }

  auto densMatrix(_systemController->getElectronicStructure<SCFMode>()->getDensityMatrix());
  auto orbitals = _systemController->getActiveOrbitalController<SCFMode>();

  const auto orbitalRangesToLocalize = separateOrbitalRanges<SCFMode>();
  for (const auto& range : orbitalRangesToLocalize) {
    localizationRoutine->localizeOrbitals(*orbitals, settings.maxSweeps, range);
  }

  DensityMatrixController<SCFMode> densMatContr(orbitals, _systemController->getNOccupiedOrbitals<SCFMode>());
  densMatContr.updateDensityMatrix();
  auto densMatNew = densMatContr.getDensityMatrix();
  double densMatrixDiff = 0.0;
  for_spin(densMatrix, densMatNew) {
    Matrix<double> densMatrixTmp = densMatrix_spin;
    densMatrixDiff += densMatrixTmp.maxCoeff();
    densMatrixTmp = densMatNew_spin;
    densMatrixDiff -= densMatrixTmp.maxCoeff();
  };
  if (_systemController->getElectronicStructure<SCFMode>()->getDiskMode() != true) {
    _systemController->getElectronicStructure<SCFMode>()->toHDF5(_systemController->getHDF5BaseName(),
                                                                 _systemController->getSystemIdentifier());
  }
  if (fabs(densMatrixDiff) >= 1.0e-10) {
    WarningTracker::printWarning(
        (std::string) "Warning: Density Matrix Changed! Largest absolute change: " + fabs(densMatrixDiff), true);
  }
  else if (fockMatrixPtr) {
    _systemController->getElectronicStructure<SCFMode>()->setFockMatrix(*fockMatrixPtr);
  }
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, std::vector<unsigned int>>> LocalizationTask::separateOrbitalRanges() {
  std::vector<SpinPolarizedData<SCFMode, std::vector<unsigned int>>> allOrbitalRanges;
  const auto valenceAndCoreOrbitalIndices = this->template getValenceOrbitalIndices<SCFMode>();
  if (settings.separateSOMOs) {
    OutputControl::nOut << "  SOMO orbitals will be localized separately." << std::endl;
    auto doMOAndSOMORanges = this->template separateSOMOs<SCFMode>(valenceAndCoreOrbitalIndices.first);
    allOrbitalRanges.push_back(doMOAndSOMORanges.first);
    allOrbitalRanges.push_back((doMOAndSOMORanges.second));
  }
  else {
    allOrbitalRanges.push_back(valenceAndCoreOrbitalIndices.first);
  }

  if (settings.splitValenceAndCore) {
    OutputControl::nOut << "  Core orbitals will be localized separately." << std::endl;
    allOrbitalRanges.push_back(valenceAndCoreOrbitalIndices.second);
  }
  if (settings.localizeVirtuals) {
    if (settings.locType == Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO ||
        settings.locType == Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN ||
        settings.locType == Options::ORBITAL_LOCALIZATION_ALGORITHMS::IAO) {
      auto virtualOrbitalRanges = this->template getVirtualValenceOrbitalIndices<SCFMode>();
      OutputControl::nOut << "  Valence virtual orbitals will be localized." << std::endl;
      if (settings.separateSOMOs) {
        auto doMOAndSOMORanges = this->template separateSOMOs<SCFMode>(virtualOrbitalRanges.first);
        allOrbitalRanges.push_back(doMOAndSOMORanges.first);
        allOrbitalRanges.push_back((doMOAndSOMORanges.second));
      }
      else {
        allOrbitalRanges.push_back(virtualOrbitalRanges.first);
      }
    }
    else {
      throw SerenityError("Virtual orbital localization is only implemented for the IBO localization scheme.");
    }
  }
  return allOrbitalRanges;
}

} /* namespace Serenity */
