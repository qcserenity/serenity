/**
 * @file ExportGridTask.cpp
 *
 * @date: Mar 24, 2016
 * @author: Jan Unsleber
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
#include "tasks/ExportGridTask.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "grid/AtomCenteredGrid.h"
#include "grid/AtomCenteredGridController.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <stdio.h>

namespace Serenity {

ExportGridTask::ExportGridTask(const std::shared_ptr<SystemController> system) : _systemController(system) {
}

void ExportGridTask::run() {
  auto systemSettings = _systemController->getSettings();
  if (settings.withAtomInfo and systemSettings.grid.gridPointSorting == true)
    WarningTracker::printWarning("WARNING: The exported grid Points do not necessarily correspond to the atomindex displayed!\n \
      Set the gridPointSorting option to false in the grid block!",
                                 iOOptions.printGridInfo);
  auto gridController = dynamic_cast<AtomCenteredGridController*>(_systemController->getGridController().get());
  auto atomGrid = gridController->getAtomGrid();
  auto points = gridController->getGridPoints();
  auto weights = gridController->getWeights();
  auto indices = atomGrid->getGridIndicesOfAtoms();
  unsigned int nAtoms = _systemController->getGeometry()->getNAtoms();
  auto systemName = _systemController->getSystemName();
  auto outptFile = fopen((_systemController->getSystemPath() + systemName + ".grid").c_str(), "w");
  if (settings.withAtomInfo) {
    fprintf(outptFile, "# %4s  %10s  %16s  %16s  %16s\n", "Atom", "x", "y", "z", "w");
    for (unsigned int iAtom = 0; iAtom != nAtoms; ++iAtom) {
      for (unsigned int iPoint = indices[iAtom].first; iPoint != indices[iAtom].second; ++iPoint) {
        fprintf(outptFile, "%4d  %+14.10E  %+14.10E  %+14.10E  %+14.10E\n", iAtom, points(0, iPoint), points(1, iPoint),
                points(2, iPoint), weights[iPoint]);
      }
    }
  }
  else {
    fprintf(outptFile, "# %10s  %16s  %16s  %16s\n", "x", "y", "z", "w");
    for (unsigned int iAtom = 0; iAtom != nAtoms; ++iAtom) {
      for (unsigned int iPoint = indices[iAtom].first; iPoint != indices[iAtom].second; ++iPoint) {
        fprintf(outptFile, " %+14.10E  %+14.10E  %+14.10E  %+14.10E\n", points(0, iPoint), points(1, iPoint),
                points(2, iPoint), weights[iPoint]);
      }
    }
  }
  fclose(outptFile);
}

} /* namespace Serenity */
