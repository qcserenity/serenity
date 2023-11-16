/**
 * @file   PlaneFileWriter.cpp
 *
 * @date   Apr 7, 2019
 * @author Anja Massolle
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
#include "io/PlaneFileWriter.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"        //Geometry
#include "grid/GridController.h"      //GridController
#include "io/FormattedOutputStream.h" //OutputControl
#include "math/Matrix.h"              //Matrix
#include "parameters/Constants.h"     //ANGSTROM_TO_BOHR

namespace Serenity {

PlaneFileWriter::PlaneFileWriter(const Settings& settings, const PlotTaskSettings& planeGridSettings)
  : DataOnGridWriter(settings),
    _planeGridSettings(planeGridSettings),
    _xUnitVector(_planeGridSettings.xUnitVector.data()),
    _yUnitVector(_planeGridSettings.yUnitVector.data()),
    _gridSpacing(_planeGridSettings.gridSpacing.data()),
    _borderWidth(_planeGridSettings.borderWidth * ANGSTROM_TO_BOHR),
    _projectCutOffRadius(_planeGridSettings.projectCutOffRadius * ANGSTROM_TO_BOHR) {
}

std::shared_ptr<GridController> PlaneFileWriter::writeHeaderAndCreateGrid(std::string filename,
                                                                          std::shared_ptr<const Geometry> geometry) {
  return this->writeHeaderAndCreateGrid(std::vector<std::string>{filename}, geometry);
} /* writeHeaderAndCreateGrid */

std::shared_ptr<GridController> PlaneFileWriter::writeHeaderAndCreateGrid(std::vector<std::string> filenames,
                                                                          std::shared_ptr<const Geometry> geometry) {
  _geometry = geometry;
  auto coordinates = _geometry->getCoordinates();
  /*
   * Save the points which define the plane (either from atoms or points)
   */
  if (_planeGridSettings.atom1 != std::numeric_limits<int>::infinity()) {
    _pointsOnPlane.push_back(coordinates.row(_planeGridSettings.atom1 - 1));
  }
  else if (std::find(_planeGridSettings.p1.begin(), _planeGridSettings.p1.end(), std::numeric_limits<double>::infinity()) ==
           _planeGridSettings.p1.end()) {
    _pointsOnPlane.push_back(convertPointToBohr(_planeGridSettings.p1));
  }
  else {
    throw SerenityError("The point 1 of the plane is not defined!");
  }
  if (_planeGridSettings.atom2 != std::numeric_limits<int>::infinity()) {
    _pointsOnPlane.push_back(coordinates.row(_planeGridSettings.atom2 - 1));
  }
  else if (std::find(_planeGridSettings.p2.begin(), _planeGridSettings.p2.end(), std::numeric_limits<double>::infinity()) ==
           _planeGridSettings.p2.end()) {
    _pointsOnPlane.push_back(convertPointToBohr(_planeGridSettings.p2));
  }
  else {
    throw SerenityError("The point 2 of the plane is not defined!");
  }
  if (_planeGridSettings.atom3 != std::numeric_limits<int>::infinity()) {
    _pointsOnPlane.push_back(coordinates.row(_planeGridSettings.atom3 - 1));
  }
  else if (std::find(_planeGridSettings.p3.begin(), _planeGridSettings.p3.end(), std::numeric_limits<double>::infinity()) ==
           _planeGridSettings.p3.end()) {
    _pointsOnPlane.push_back(convertPointToBohr(_planeGridSettings.p3));
  }
  else {
    throw SerenityError("The point 3 of the plane is not defined!");
  }

  Eigen::Vector3d e1 = (_pointsOnPlane[1] - _pointsOnPlane[0]);
  Eigen::Vector3d e2 = (_pointsOnPlane[2] - _pointsOnPlane[0]);

  // check if the three points are non-collinear
  if ((e1.cross(e2)).norm() < 1.0e-6) {
    throw SerenityError("The given points do not define a plane, they are collinear!");
  }

  // calculate unitvectors and normal vector of the plane
  // _pointsOnPlane[0] will be the origin of the plane
  _normalVector = Eigen::Vector3d((e1.cross(e2)).normalized());
  _xUnitVector = Eigen::Vector3d(e1.normalized());
  _yUnitVector = Eigen::Vector3d((_normalVector.cross(e1)).normalized());

  /*
   * Find the rotation matrix which rotates the unit vectors of the plane onto (1.0, 0.0, 0.0) and (0.0, 1.0, 0.0).
   * The Triad method is used to calculate the rotation matrix, see also https://en.wikipedia.org/wiki/Triad_method
   */

  e1 = Eigen::Vector3d(1.0, 0.0, 0.0);
  e2 = Eigen::Vector3d(0.0, 1.0, 0.0);

  Eigen::Vector3d M = e1.cross(e2) / ((e1.cross(e2)).norm());
  Eigen::Vector3d m = _xUnitVector.cross(_yUnitVector) / ((_xUnitVector.cross(_yUnitVector)).norm());
  Eigen::Vector3d SxM = e1.cross(M);
  Eigen::Vector3d sxm = _xUnitVector.cross(m);

  Eigen::Matrix3d SMSxM(3, 3);
  SMSxM.col(0) = e1;
  SMSxM.col(1) = M;
  SMSxM.col(2) = SxM;

  Eigen::Matrix3d smsxm(3, 3);
  smsxm.col(0) = _xUnitVector;
  smsxm.col(1) = m;
  smsxm.col(2) = sxm;

  _rotationMatrix = SMSxM * smsxm.transpose();

  /*
   * Project the atoms on the plane and save them in _projectedPointsOnPlane
   */
  unsigned int nProjectedAtoms = 0;
  _projectedPointsOnPlane = Eigen::Matrix3Xd(0, 0);
  for (unsigned int i = 0; i < coordinates.rows(); i++) {
    Eigen::Vector3d vecToPoint = (coordinates.row(i)).transpose() - _pointsOnPlane[0];
    double distance = _normalVector.dot(vecToPoint);
    if (std::fabs(distance) <= _projectCutOffRadius) {
      _projectedPointsOnPlane.conservativeResize(_projectedPointsOnPlane.rows(), _projectedPointsOnPlane.cols() + 1);
      _projectedPointsOnPlane.col(nProjectedAtoms) = (coordinates.row(i)).transpose() - distance * _normalVector;
      nProjectedAtoms += 1;
    } /* check if the distance of an atom is within the cut off radius */
  }   /* Loop over all coordinates */

  /*
   * Find min / max values on the plane
   */
  Eigen::VectorXd xValues(_projectedPointsOnPlane.cols());
  Eigen::VectorXd yValues(_projectedPointsOnPlane.cols());

  for (unsigned int i = 0; i < _projectedPointsOnPlane.cols(); i++) {
    xValues[i] = (_projectedPointsOnPlane.col(i) - _pointsOnPlane[0]).dot(_xUnitVector);
    yValues[i] = (_projectedPointsOnPlane.col(i) - _pointsOnPlane[0]).dot(_yUnitVector);
  }

  Eigen::Vector3d min = (xValues.minCoeff() - _borderWidth) * _xUnitVector + (yValues.minCoeff() - _borderWidth) * _yUnitVector;

  /*
   * Calculate number of points in both directions of the plane
   */
  const double stepX = _gridSpacing[0] * ANGSTROM_TO_BOHR;
  const double stepY = _gridSpacing[1] * ANGSTROM_TO_BOHR;
  _nPointsX = (abs(xValues.maxCoeff() - xValues.minCoeff()) + 2 * _borderWidth) / stepX;
  _nPointsY = (abs(yValues.maxCoeff() - yValues.minCoeff()) + 2 * _borderWidth) / stepY;
  const unsigned int nPointsTotal = _nPointsX * _nPointsY;

  // Create the matrix for storing the gridpoints on the plane and the weight
  // vector (since the grid is equidistant every point has the same weight)
  std::unique_ptr<Eigen::Matrix3Xd> points(new Eigen::Matrix3Xd(3, nPointsTotal));
  std::unique_ptr<Eigen::VectorXd> weights(new Eigen::VectorXd(nPointsTotal));
  (*weights) = Eigen::VectorXd::Constant(nPointsTotal, 1.0 / (double)nPointsTotal);

  /*
   * Create plane grid
   */
  Eigen::Vector3d gridPointX;
  gridPointX.setZero();
  unsigned long long npt = 0;
  for (unsigned int ix = 0; ix < _nPointsX; ++ix, gridPointX += stepX * _xUnitVector) {
    Eigen::Vector3d gridPointY;
    gridPointY.setZero();
    for (unsigned int iy = 0; iy < _nPointsY; ++iy, gridPointY += stepY * _yUnitVector, npt++) {
      (*points).col(iy + ix * _nPointsY) = min + _pointsOnPlane[0] + gridPointX + gridPointY;
    } /* Loop over y direction */
  }   /* Loop over x direction */

  /*
   * Write header
   */
  for (const auto& filename : filenames) {
    std::string fullFileName = filename + ".dat";
    FILE* file = fopen(fullFileName.data(), "w");
    fprintf(file, "%s %12s %12s %16s \n", "X", "Y", "Z", "value");
    fclose(file);

    if (_planeGridSettings.xyHeatmap) {
      std::string fullXYfileName = filename + "_XYPLANE.dat";
      FILE* xyFile = fopen(fullXYfileName.data(), "w");
      fprintf(xyFile, "%s %12s %16s \n", "X", "Y", "value");
      fclose(xyFile);
    }
  }
  _planeGridController =
      std::make_shared<GridController>(std::unique_ptr<Grid>(new Grid(std::move(points), std::move(weights))));
  return _planeGridController;
}; /* writeHeaderAndCreateGrid */

void PlaneFileWriter::writeData(std::string filename, const Eigen::VectorXd& data) {
  /*
   * Print data (x, y, z, value)
   */
  auto gridPoints = _planeGridController->getGridPoints();
  std::string fullFileName = filename + ".dat";
  FILE* file = fopen(fullFileName.data(), "a");
  for (unsigned int i = 0; i < data.size(); i++) {
    fprintf(file, "%E %E %E %E \n", gridPoints(0, i) * BOHR_TO_ANGSTROM, gridPoints(1, i) * BOHR_TO_ANGSTROM,
            gridPoints(2, i) * BOHR_TO_ANGSTROM, data[i]);
  }
  fclose(file);

  /*
   * Print grid data, rotated to x y plane
   */
  if (_planeGridSettings.xyHeatmap) {
    Eigen::Vector3d point = _rotationMatrix * gridPoints.col(0);
    double xDistance = point[0] * BOHR_TO_ANGSTROM;
    double yDistance = point[1] * BOHR_TO_ANGSTROM;
    double zDistance = point[2] * BOHR_TO_ANGSTROM;
    std::string fullXYfileName = filename + "_XYPLANE.dat";
    FILE* xyFile = fopen(fullXYfileName.data(), "a");
    for (unsigned int i = 0; i < data.size(); i++) {
      Eigen::Vector3d xyPoint = _rotationMatrix * gridPoints.col(i);
      fprintf(xyFile, "%E %E %E \n", xyPoint[0] * BOHR_TO_ANGSTROM - xDistance,
              xyPoint[1] * BOHR_TO_ANGSTROM - yDistance, data[i]);
    }
    fclose(xyFile);

    /*
     * Print molecule rotated the same way, as the plane
     */
    std::string fullGeometryXYfileName = filename + "_MOLECULE_ROTATED_TO_XYPLANE.xyz";
    FILE* xyGeometryFile = fopen(fullGeometryXYfileName.data(), "w");
    auto atoms = _geometry->getAtoms();
    auto coordinates = _geometry->getCoordinates();
    auto elements = _geometry->getAtomSymbols();
    bool containDummyAtom = false;
    fprintf(xyGeometryFile, "%ld \n\n", coordinates.rows());
    for (unsigned int i = 0; i < coordinates.rows(); i++) {
      Eigen::Vector3d xyGeometryPoint = _rotationMatrix * (coordinates.row(i)).transpose();
      const char* elementSymbol = (elements[i]).c_str();
      /*
       * ASE can not interpret elements with :. Dummy atoms are
       * marked with a X in ASE
       */
      if (atoms[i]->isDummy()) {
        elementSymbol = "X";
        containDummyAtom = true;
      }
      fprintf(xyGeometryFile, "%s %g %g %g\n", elementSymbol, xyGeometryPoint[0] * BOHR_TO_ANGSTROM - xDistance,
              xyGeometryPoint[1] * BOHR_TO_ANGSTROM - yDistance, xyGeometryPoint[2] * BOHR_TO_ANGSTROM - zDistance);
    }
    if (containDummyAtom) {
      OutputControl::vOut << "Element of the dummy atoms is changed to X in " << fullGeometryXYfileName << std::endl;
    }

    fclose(xyGeometryFile);
  } /*writeHeatmap*/

} /* writeData*/

} /* namespace Serenity */
