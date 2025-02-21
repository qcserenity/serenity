/**
 * @file   CubeFileWriter.cpp
 *
 * @date   Oct 7, 2014
 * @author Jan Unsleber
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
#include "io/CubeFileWriter.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"   //Geometry
#include "grid/GridController.h" //GridController
#include "parameters/Constants.h"

namespace Serenity {

CubeFileWriter::CubeFileWriter(const Settings& settings, const PlotTaskSettings& cubeGridSettings)
  : DataOnGridWriter(settings),
    _cubeGridSettings(cubeGridSettings),
    _xUnitVector(_cubeGridSettings.xUnitVector.data()),
    _yUnitVector(_cubeGridSettings.yUnitVector.data()),
    _zUnitVector(_cubeGridSettings.zUnitVector.data()),
    _gridSpacing(_cubeGridSettings.gridSpacing.data()),
    _borderWidth(_cubeGridSettings.borderWidth * ANGSTROM_TO_BOHR){};

std::shared_ptr<GridController> CubeFileWriter::writeHeaderAndCreateGrid(std::string filename,
                                                                         std::shared_ptr<const Geometry> geometry) {
  return writeHeaderAndCreateGrid(std::vector<std::string>{filename}, geometry);
} /* writeHeaderAndCreateGrid */

std::shared_ptr<GridController> CubeFileWriter::writeHeaderAndCreateGrid(std::vector<std::string> filenames,
                                                                         std::shared_ptr<const Geometry> geometry) {
  /*
   * Create grid
   */
  const double minX = geometry->getMinX() - _borderWidth;
  const double minY = geometry->getMinY() - _borderWidth;
  const double minZ = geometry->getMinZ() - _borderWidth;
  const double maxX = geometry->getMaxX() + _borderWidth;
  const double maxY = geometry->getMaxY() + _borderWidth;
  const double maxZ = geometry->getMaxZ() + _borderWidth;
  const double stepX = _gridSpacing[0] * ANGSTROM_TO_BOHR;
  const double stepY = _gridSpacing[1] * ANGSTROM_TO_BOHR;
  const double stepZ = _gridSpacing[2] * ANGSTROM_TO_BOHR;
  const unsigned int nPointsX = (maxX - minX) / stepX;
  const unsigned int nPointsY = (maxY - minY) / stepY;
  const unsigned int nPointsZ = (maxZ - minZ) / stepZ;
  const unsigned int nPointsTotal = nPointsX * nPointsY * nPointsZ;
  if (_cubeGridSettings.maxGridPoints < nPointsTotal) {
    throw SerenityError("You are about to print files of at least " +
                        std::to_string(_cubeGridSettings.maxGridPoints * sizeof(double)) + " bytes.\n\
        Cubic grid construction is stopped to prevent a possible memory bottleneck. If you want to print these files, you can adjust the maximum number of grid points in the input.");
  }
  std::unique_ptr<Eigen::Matrix3Xd> points(new Eigen::Matrix3Xd(3, nPointsTotal));
  std::unique_ptr<Eigen::VectorXd> weights(new Eigen::VectorXd(nPointsTotal));
  (*weights) = Eigen::VectorXd::Constant(nPointsTotal, 1.0 / (double)nPointsTotal);
  double x = minX;
  unsigned long long npt = 0;
  for (unsigned int ix = 0; ix < nPointsX;
       ++ix, x += stepX * _xUnitVector[0] + stepY * _yUnitVector[0] + stepZ * _zUnitVector[0]) {
    double y = minY;
    for (unsigned int iy = 0; iy < nPointsY;
         ++iy, y += stepX * _xUnitVector[1] + stepY * _yUnitVector[1] + stepZ * _zUnitVector[1]) {
      double z = minZ;
      for (unsigned int iz = 0; iz < nPointsZ;
           npt++, ++iz, z += stepX * _xUnitVector[2] + stepY * _yUnitVector[2] + stepZ * _zUnitVector[2]) {
        (*points)(0, npt) = x;
        (*points)(1, npt) = y;
        (*points)(2, npt) = z;
      } /* Loop over z direction */
    }   /* Loop over y direction */
  }     /* Loop over x direction */
  /*
   * Print headers
   */
  for (const auto& filename : filenames) {
    std::string fullFileName = filename + ".cube";
    FILE* file = fopen(fullFileName.data(), "w");
    fprintf(file, "%s\n\n", "CUBE FILE");
    const unsigned int nAtoms = geometry->getNAtoms();
    fprintf(file, "%u %f %f %f\n", nAtoms, minX, minY, minZ);
    fprintf(file, "%u %f %f %f\n", nPointsX, stepX * _xUnitVector[0], stepX * _xUnitVector[1], stepX * _xUnitVector[2]);
    fprintf(file, "%u %f %f %f\n", nPointsY, stepY * _yUnitVector[0], stepY * _yUnitVector[1], stepY * _yUnitVector[2]);
    fprintf(file, "%u %f %f %f\n", nPointsZ, stepZ * _zUnitVector[0], stepZ * _zUnitVector[1], stepZ * _zUnitVector[2]);
    const auto& atoms = geometry->getAtoms();
    for (unsigned int i = 0; i < nAtoms; ++i) {
      fprintf(file, "%u %f %f %f %f\n", atoms[i]->getNuclearCharge(), 0.0, atoms[i]->getX(), atoms[i]->getY(),
              atoms[i]->getZ());
    }
    fclose(file);
  }
  return std::make_shared<GridController>(std::unique_ptr<Grid>(new Grid(std::move(points), std::move(weights))));
} /* writeHeaderAndCreateGrid */

void CubeFileWriter::writeData(std::string filename, const Eigen::VectorXd& data) {
  /*
   * Print data
   */
  std::string fullFileName = filename + ".cube";
  FILE* file = fopen(fullFileName.data(), "a");
  for (unsigned int i = 0; i < data.size(); i++) {
    fprintf(file, "%g ", data[i]);
    if (i % 6 == 5)
      fprintf(file, "\n");
  }
  fclose(file);
} /* writeData */

} /*namespace Serenity*/
