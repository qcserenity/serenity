/**
 * @file CHELPGPopulationCalculator.cpp
 *
 * @date October 7, 2024
 * @author: Thorben Wiegmann
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
#include "CHELPGPopulationCalculator.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/grid/ElectrostaticPotentialOnGridController.h"
#include "geometry/Geometry.h"
#include "grid/Grid.h"
#include "grid/GridController.h"
#include "io/DataOnGridWriter.h"
#include "io/GeneralGridFileWriter.h"
#include "parameters/AtomicParameters.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <algorithm>
#include <iostream>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CHELPGPopulationCalculator<SCFMode>::CHELPGPopulationCalculator(std::shared_ptr<SystemController> system)
  : _system(system), _atomPopulations(nullptr) {
  _atoms = _system->getGeometry()->getAtoms();
  _nAtoms = _atoms.size();
  _headspace = 2.8 * ANGSTROM_TO_BOHR;
  _pointDistance = 0.3 * ANGSTROM_TO_BOHR;
}

template<Options::SCF_MODES SCFMode>
void CHELPGPopulationCalculator<SCFMode>::calculateCHELPGPopulations() {
  Eigen::VectorXd charges = this->calculateCharges();
  std::shared_ptr<Eigen::VectorXd> atomPopulations = std::make_shared<Eigen::VectorXd>(Eigen::VectorXd::Zero(_nAtoms));
  for (unsigned int iAtom = 0; iAtom < _nAtoms; iAtom++) {
    (*atomPopulations)[iAtom] = charges[iAtom];
  }
  _atomPopulations = atomPopulations;
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd CHELPGPopulationCalculator<SCFMode>::calculateCharges() {
  /*
   * The grid will be constructed according to the specifications given in [1]. In the following a brief summary will be
   * given: 1) Equidistant points (0.3 angstrom) 2) The grid must include all atoms 3) The grid has 2.8 angstrom
   * headspace on all sides 4) All points that lie in a VDW radius of any atom will be removed 5) All points that are
   * further than 2.8 angstrom away from any atom will be removed
   */

  // Determine the space the grid will take up
  double xMax = _system->getGeometry()->getMaxX() + _headspace;
  double yMax = _system->getGeometry()->getMaxY() + _headspace;
  double zMax = _system->getGeometry()->getMaxZ() + _headspace;

  double xMin = _system->getGeometry()->getMinX() - _headspace;
  double yMin = _system->getGeometry()->getMinY() - _headspace;
  double zMin = _system->getGeometry()->getMinZ() - _headspace;

  double gridX = xMax - xMin;
  double gridY = yMax - yMin;
  double gridZ = zMax - zMin;
  // Determine number of grid points
  unsigned int nPoints =
      std::ceil(gridX / _pointDistance) * std::ceil(gridY / _pointDistance) * std::ceil(gridZ / _pointDistance);
  // Calculate the coordinates of the grid points - the grid will be built from the minimum corner to the maximum corner
  Eigen::MatrixXd points = Eigen::MatrixXd(3, nPoints).setZero();
  unsigned int iCol = 0;
  double x = xMin;
  while (x <= xMax) {
    double y = yMin;
    while (y <= yMax) {
      double z = zMin;
      while (z <= zMax) {
        points.col(iCol) << x, y, z;
        iCol++;
        z += _pointDistance;
      } // z loop
      y += _pointDistance;
    } // y loop
    x += _pointDistance;
  } // x loop
  // Remove points according to 4) and 5)
  // Find points that need to be removed
  std::vector<unsigned int> removePoints;
  for (unsigned int iPoint = 0; iPoint < points.cols(); iPoint++) {
    unsigned int atomsExceeded = 0;
    for (unsigned int iAtom = 0; iAtom < _nAtoms; iAtom++) {
      double vdwRadius = ELEMENTAL_VAN_DER_WAALS_RADII[_atoms[iAtom]->getNuclearCharge()] * ANGSTROM_TO_BOHR;
      double dist = Serenity::distance(*_atoms[iAtom], Point(points(0, iPoint), points(1, iPoint), points(2, iPoint)));
      if (dist <= vdwRadius) {
        removePoints.push_back(iPoint);
      }
      if (dist > _headspace) {
        atomsExceeded++;
      }
    } // atom loop
    if (atomsExceeded == _nAtoms) {
      removePoints.push_back(iPoint);
    }
  } // point loop
  // Make sure the remove vector has no duplicates
  std::sort(removePoints.begin(), removePoints.end());
  auto newEnd = std::unique(removePoints.begin(), removePoints.end());
  removePoints.erase(newEnd, removePoints.end());
  // Actually remove points
  unsigned int nCols = points.cols() - removePoints.size();
  Eigen::MatrixXd reducedPoints = Eigen::MatrixXd(3, nCols).setZero();
  unsigned int nRemoved = 0;
  for (unsigned int iCol = 0; iCol < points.cols(); iCol++) {
    // Check if the current col is in the remove list
    if (std::find(removePoints.begin(), removePoints.end(), iCol) != removePoints.end()) {
      nRemoved++;
      continue;
    }
    reducedPoints.col(iCol - nRemoved) = points.col(iCol);
  }
  // Construct grid
  std::unique_ptr<Eigen::Matrix3Xd> gridPoints =
      std::make_unique<Eigen::Matrix3Xd>(Eigen::Matrix3Xd::Zero(3, reducedPoints.cols()));
  for (unsigned int i = 0; i < reducedPoints.cols(); i++) {
    gridPoints->col(i) = reducedPoints.col(i);
  }
  std::unique_ptr<Eigen::VectorXd> weights = std::make_unique<Eigen::VectorXd>(Eigen::VectorXd::Ones(gridPoints->cols()));
  auto grid = std::make_unique<Grid>(std::move(gridPoints), std::move(weights));
  std::shared_ptr<GridController> gridController = std::make_shared<GridController>(std::move(grid));
  // Calculate Charge
  auto densMatrixController = _system->getElectronicStructure<SCFMode>()->getDensityMatrixController();
  ElectrostaticPotentialOnGridController<SCFMode> electrostaticPotentialController(
      gridController, densMatrixController, _system->getGeometry(), _system->getSettings().name);
  Eigen::VectorXd potential = electrostaticPotentialController.getPotential();
  Eigen::Matrix3Xd finalPoints = gridController->getGridPoints();
  Eigen::MatrixXd reciprocalDistances = Eigen::MatrixXd(finalPoints.cols(), _nAtoms).setZero();
  for (unsigned int iAtom = 0; iAtom < _nAtoms; iAtom++) {
    for (unsigned int iPoint = 0; iPoint < finalPoints.cols(); iPoint++) {
      Point p(finalPoints(0, iPoint), finalPoints(1, iPoint), finalPoints(2, iPoint));
      reciprocalDistances(iPoint, iAtom) = 1.0 / Serenity::distance(*_atoms[iAtom], p);
    } // atom loop
  }   // point loop
  Eigen::VectorXd charges = reciprocalDistances.colPivHouseholderQr().solve(potential);

  return charges;
}

template class CHELPGPopulationCalculator<Options::SCF_MODES::RESTRICTED>;
template class CHELPGPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
