/**
 * @file NumericalDipoleMomentCalculator.cpp
 *
 * @date Oct 30, 2015
 * @author David Schnieders
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
#include "analysis/multipoles/NumericalDipoleMomentCalculator.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/DensityOnGridController.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "geometry/Point.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
Eigen::Vector3d NumericalDipoleMomentCalculator::calculateDipoleMoment(std::shared_ptr<SystemController> system) {
  /*
   * Gather Data
   */

  // geometry
  auto geometry = system->getGeometry();

  // density on grid
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      system->getSettings(), system->getBasisController(), system->getGridController());
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
      basisFunctionOnGridController, system->getSettings().grid.blockAveThreshold);
  auto densityOnGrid =
      densityOnGridCalculator->calcDensityOnGrid(system->getElectronicStructure<SCFMode>()->getDensityMatrix());
  auto grid = densityOnGrid.getGridController();
  /*
   * Calculate Dipole Moment: mu=SUM(Qi*Ri)-SUM(Wj*RHOj*Rj)
   */

  // Cores
  Eigen::Vector3d dipoleMoment(3);
  dipoleMoment.setZero();
  for (const auto& atom : geometry->getAtoms()) {
    const double charge = atom->getEffectiveCharge();
    dipoleMoment[0] += charge * atom->getX();
    dipoleMoment[1] += charge * atom->getY();
    dipoleMoment[2] += charge * atom->getZ();
  };

  // Electrons
  const auto& weights = grid->getWeights();
  const auto& points = grid->getGridPoints();
  for (unsigned int i = 0; i < grid->getNGridPoints(); i++) {
    for_spin(densityOnGrid) {
      dipoleMoment.array() -= weights[i] * points.col(i).array() * densityOnGrid_spin[i];
    };
  };

  return dipoleMoment;
}

template Eigen::Vector3d
NumericalDipoleMomentCalculator::calculateDipoleMoment<Options::SCF_MODES::RESTRICTED>(std::shared_ptr<SystemController> system);
template Eigen::Vector3d
NumericalDipoleMomentCalculator::calculateDipoleMoment<Options::SCF_MODES::UNRESTRICTED>(std::shared_ptr<SystemController> system);

} // namespace Serenity
