/**
 * @file   LocalizationTask.cpp
 *
 * @date   Apr 22, 2014
 * @author Thomas Dresselhaus
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
#include "tasks/LocalizationTask.h"
/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "analysis/orbitalLocalization/EdmistonRuedenbergLocalization.h"
#include "data/ElectronicStructure.h"
#include "io/FormattedOutput.h"
#include "analysis/orbitalLocalization/FosterBoysLocalization.h"
#include "analysis/orbitalLocalization/IBOLocalization.h"
#include "math/Matrix.h"
#include "analysis/orbitalLocalization/NonOrthogonalLocalization.h"
#include "integrals/OneElectronIntegralController.h"
#include "data/OrbitalController.h"
#include "analysis/orbitalLocalization/PipekMezeyLocalization.h"
#include "system/SystemController.h"
#include "misc/WarningTracker.h"
/* Include Std and External Headers */
#include <Eigen/StdVector>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>


namespace Serenity {
using namespace std;

LocalizationTask::LocalizationTask(shared_ptr<SystemController> systemController) :
    _systemController(systemController){
  assert(_systemController);
}

void LocalizationTask::run() {
  if (_systemController->getLastSCFMode() == Options::SCF_MODES::UNRESTRICTED
      or _systemController->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
    runByLastSCFMode<Options::SCF_MODES::UNRESTRICTED>();
  } else {
    runByLastSCFMode<Options::SCF_MODES::RESTRICTED>();
  }
}

template<Options::SCF_MODES T> void LocalizationTask::runByLastSCFMode() {


  std::shared_ptr<Localization<T> > localizationRoutine;
  switch(settings.locType){
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS:
      localizationRoutine = std::make_shared<FosterBoysLocalization<T> >(_systemController);
      printSubSectionTitle("Running Foster-Boys Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY:
      localizationRoutine = std::make_shared<PipekMezeyLocalization<T> >(_systemController);
      printSubSectionTitle("Running Pipek-Mezey Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO:
      localizationRoutine = std::make_shared<IBOLocalization<T> >(_systemController, false);
      printSubSectionTitle("Running IBO Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::IAO:
      localizationRoutine = std::make_shared<IBOLocalization<T> >(_systemController, true);
      printSubSectionTitle("Running IAO Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::EDMISTON_RUEDENBERG:
      localizationRoutine = std::make_shared<EdmistonRuedenbergLocalization<T> >(_systemController);
      printSubSectionTitle("Running Edmiston-Ruedenberg Localization");
      break;
    case Options::ORBITAL_LOCALIZATION_ALGORITHMS::NON_ORTHOGONAL:
      localizationRoutine = std::make_shared<NonOrthogonalLocalization<T> >(_systemController);
      printSubSectionTitle("Running Non-Orthogonal Localization");
      break;
  }

  auto densMatrix(
      _systemController->getElectronicStructure<T>()->getDensityMatrix());
  auto orbitals = _systemController->getActiveOrbitalController<T>();
  localizationRoutine->localizeOrbitals(*orbitals,settings.maxSweeps);

  DensityMatrixController<T> densMatContr(orbitals,_systemController->getNOccupiedOrbitals<T>());
  densMatContr.updateDensityMatrix();
  auto densMatNew = densMatContr.getDensityMatrix();
  double densMatrixDiff = 0.0;
  for_spin(densMatrix,densMatNew){
    Matrix<double> densMatrixTmp = densMatrix_spin;
    densMatrixDiff += densMatrixTmp.maxCoeff();
    densMatrixTmp = densMatNew_spin;
    densMatrixDiff -= densMatrixTmp.maxCoeff();
  };
    if (fabs(densMatrixDiff) >= 1.0e-10) {
      WarningTracker::printWarning("Warning: Density Matrix Changed!",true);
    };
}

} /* namespace Serenity */
