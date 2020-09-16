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
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h" //Filtered output streams.
#include "math/Matrix.h"
#include "misc/WarningTracker.h" //Warnings.
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <Eigen/StdVector>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>

namespace Serenity {
using namespace std;

LocalizationTask::LocalizationTask(std::shared_ptr<SystemController> systemController,
                                   std::vector<std::shared_ptr<SystemController>> templateSystem)
  : _systemController(systemController), _templateSystem(templateSystem) {
  assert(_systemController);
}

void LocalizationTask::run() {
  if (_systemController->getLastSCFMode() == Options::SCF_MODES::UNRESTRICTED or
      _systemController->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
    runByLastSCFMode<Options::SCF_MODES::UNRESTRICTED>();
  }
  else {
    runByLastSCFMode<Options::SCF_MODES::RESTRICTED>();
  }
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
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> valenceRange;
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> coreRange;
  if (!settings.splitValenceAndCore) {
    for_spin(nOcc, valenceRange) {
      for (unsigned int iOcc = 0; iOcc < nOcc_spin; ++iOcc)
        valenceRange_spin.push_back(iOcc);
    };
    return std::make_pair(valenceRange, coreRange);
  }
  const auto eigenvalues = _systemController->getActiveOrbitalController<SCFMode>()->getEigenvalues();
  for_spin(nOcc, valenceRange, eigenvalues, coreRange) {
    for (unsigned int iOcc = 0; iOcc < nOcc_spin; ++iOcc) {
      if (eigenvalues_spin[iOcc] > settings.energyCutOff) {
        valenceRange_spin.push_back(iOcc);
      }
      else {
        coreRange_spin.push_back(iOcc);
      }
    } // for iOcc
  };  // for spin
  OutputControl::mOut << "  Note: Valence and core orbitals are localized separately." << std::endl;
  OutputControl::mOut << "        Please ensure that the initial orbitals were canonical orbitals!\n" << std::endl;
  OutputControl::vOut << "The following orbitals have been identified as valence orbitals:" << std::endl;
  SpinPolarizedData<SCFMode, std::string> praefix = getOutputPraefix<SCFMode>();
  for_spin(valenceRange, praefix, nOcc) {
    OutputControl::vOut << praefix_spin << std::endl << "  ";
    for (const auto orb : valenceRange_spin)
      OutputControl::vOut << orb + 1 << " ";
    OutputControl::vOut << std::endl;
    OutputControl::nOut << "  Number of core-like orbitals: " << nOcc_spin - valenceRange_spin.size() << "  "
                        << praefix_spin << std::endl;
  };
  OutputControl::nOut << "  Core-like orbitals will be localized separately." << std::endl;
  OutputControl::nOut << "  Running orbital localization for valence orbitals." << std::endl;
  return std::make_pair(valenceRange, coreRange);
}
template<Options::SCF_MODES T>
void LocalizationTask::runByLastSCFMode() {
  std::shared_ptr<Localization<T>> localizationRoutine;
  switch (settings.locType) {
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS:
      localizationRoutine = std::make_shared<FosterBoysLocalization<T>>(_systemController);
      printSubSectionTitle("Running Foster-Boys Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY:
      localizationRoutine = std::make_shared<PipekMezeyLocalization<T>>(_systemController);
      printSubSectionTitle("Running Pipek-Mezey Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO:
      localizationRoutine = std::make_shared<IBOLocalization<T>>(_systemController, false);
      printSubSectionTitle("Running IBO Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::IAO:
      // TODO: Does this actually do anything? Projecting on the IAOs and back should not change the orbitals.
      localizationRoutine = std::make_shared<IBOLocalization<T>>(_systemController, true);
      printSubSectionTitle("Running IAO Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::EDMISTON_RUEDENBERG:
      localizationRoutine = std::make_shared<EdmistonRuedenbergLocalization<T>>(_systemController);
      printSubSectionTitle("Running Edmiston-Ruedenberg Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::NON_ORTHOGONAL:
      localizationRoutine = std::make_shared<NonOrthogonalLocalization<T>>(_systemController);
      printSubSectionTitle("Running Non-Orthogonal Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN:
      if (_templateSystem.size() == 0) {
        throw SerenityError("Orbital alignment was requested without any template system.\
                           Please specify a template system in the task input via the env keyword.");
      }
      localizationRoutine = std::make_shared<OrbitalAligner<T>>(_systemController, _templateSystem[0],
                                                                settings.alignExponent, settings.useKineticAlign);
      printSubSectionTitle("Running IBO-Like Orbital Alignment");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::NONE:
      // Do nothing
      printSubSectionTitle("Skip localization...");
      return;
  }

  auto densMatrix(_systemController->getElectronicStructure<T>()->getDensityMatrix());
  auto orbitals = _systemController->getActiveOrbitalController<T>();
  auto valenceAndCoreOrbitalIndices = getValenceOrbitalIndices<T>();
  localizationRoutine->localizeOrbitals(*orbitals, settings.maxSweeps, valenceAndCoreOrbitalIndices.first);
  if (settings.splitValenceAndCore) {
    OutputControl::nOut << "  Running orbital localization for core orbitals." << std::endl;
    localizationRoutine->localizeOrbitals(*orbitals, settings.maxSweeps, valenceAndCoreOrbitalIndices.second);
  }

  DensityMatrixController<T> densMatContr(orbitals, _systemController->getNOccupiedOrbitals<T>());
  densMatContr.updateDensityMatrix();
  auto densMatNew = densMatContr.getDensityMatrix();
  double densMatrixDiff = 0.0;
  for_spin(densMatrix, densMatNew) {
    Matrix<double> densMatrixTmp = densMatrix_spin;
    densMatrixDiff += densMatrixTmp.maxCoeff();
    densMatrixTmp = densMatNew_spin;
    densMatrixDiff -= densMatrixTmp.maxCoeff();
  };
  if (fabs(densMatrixDiff) >= 1.0e-10) {
    WarningTracker::printWarning(
        (std::string) "Warning: Density Matrix Changed! Largest absolute change: " + fabs(densMatrixDiff), true);
  };
  const auto systemSettings = _systemController->getSettings();
  _systemController->getElectronicStructure<T>()->toHDF5(systemSettings.path + systemSettings.name, systemSettings.identifier);
}

} /* namespace Serenity */
