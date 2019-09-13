/**
 * @file   CubeFileWriter.cpp
 *
 * @date   Oct 7, 2014
 * @author Jan Unsleber
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
#include "io/CubeFileWriter.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "math/FloatMaths.h"
#include "geometry/Geometry.h"
#include "grid/Grid.h"
#include "grid/GridController.h"
#include "data/grid/MOCalculator.h"
#include "data/matrices/MatrixInBasis.h"
#include "data/grid/MatrixOperatorToGridTransformer.h"
#include "settings/Options.h"
#include "geometry/Point.h"
#include "settings/Settings.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <cassert>


namespace Serenity {
using namespace std;

CubeFileWriter::CubeFileWriter(const Settings& settings) :
    _xUnitVector(1.0, 0.0, 0.0),
    _yUnitVector(0.0, 1.0, 0.0),
    _zUnitVector(0.0, 0.0, 1.0),
    _stepSize(0.1, 0.1, 0.1),
    _settings(settings){
}

void CubeFileWriter::writeMatrixToCube(
    const string& filename,
    shared_ptr<const Geometry> geometry,
    const MatrixInBasis<RESTRICTED>& matrix) {
  auto basisController = matrix.getBasisController();
  this->writeMatrixToCube(filename,geometry,basisController,matrix);
}

void CubeFileWriter::writeMatrixToCube(
      const string& filename,
      shared_ptr<const Geometry> geometry,
      std::shared_ptr<BasisController> basisController,
      const Eigen::MatrixXd& matrix) {
  /*
   * define lambda function
   */
  auto calcValuesOnGrid =
    [&basisController,&matrix,this](shared_ptr<GridController> cubeGridController) -> Eigen::VectorXd {
      auto basFuncOnGridController =
        BasisFunctionOnGridControllerFactory::produce(
            _settings, basisController, cubeGridController);
      Eigen::VectorXd result;
      MatrixInBasis<RESTRICTED> tmp(basisController);
      tmp = matrix;
      MatrixOperatorToGridTransformer::transform(
        tmp,
        result,
        *basFuncOnGridController);
      return result;
    };
  // write out the data
  writeCube(filename, geometry, calcValuesOnGrid);
}

void CubeFileWriter::writeMatrixToCube(
    const vector<string>& filenames,
    shared_ptr<const Geometry> geometry,
    const vector<std::reference_wrapper<const MatrixInBasis<RESTRICTED> > >& matrices) {
  assert (filenames.size() == matrices.size());
  // Make sure all matrices are defined in the same basis
  for (auto& matrix : matrices){
    if (!isDefinedInSameBasis(matrix.get(), matrices[0].get()))
      throw SerenityError("CubeFileWriter: Components are not defined on the same grid.");
  }
  auto basisController = matrices[0].get().getBasisController();
  /*
   * define lambda function
   */
  auto calcValuesOnGrid =
    [&basisController,&matrices,this](shared_ptr<GridController> cubeGridController)
            -> vector<Eigen::VectorXd> {
      auto basFuncOnGridController =
        BasisFunctionOnGridControllerFactory::produce(
            _settings, basisController, cubeGridController);
      vector<Eigen::VectorXd> results(matrices.size());
      vector<reference_wrapper<Eigen::VectorXd> > tmp;
      for (auto& result : results) tmp.push_back(result);
      MatrixOperatorToGridTransformer::transform(
        matrices,
        tmp,
        *basFuncOnGridController);
      return results;
    };
  // write out the data
  writeCube(filenames, geometry, calcValuesOnGrid);
}

void CubeFileWriter::writeVectorSetToCube(
    const std::vector<string>& filenames,
    shared_ptr<const Geometry> geometry,
    shared_ptr<BasisController> basisController,
    Eigen::MatrixXd& inVectorSet) {
  assert(filenames.size()==inVectorSet.cols());
  /*
   * define lambda function
   */
  auto calcValuesOnGrid =
      [&basisController,&inVectorSet,this](shared_ptr<GridController> cubeGridController) -> Eigen::MatrixXd{

    //get useful data
    auto basFuncOnGridController =
            BasisFunctionOnGridControllerFactory::produce(
                _settings, basisController, cubeGridController);

    //use MOCalculator to calculate MO on grid
    MOCalculator moCalc(basFuncOnGridController);
    Eigen::MatrixXd valuesOnGrid=moCalc.calcMOValuesOnGrid(inVectorSet);


    return valuesOnGrid;

    };
  // write out the data
  writeCube(filenames, geometry, calcValuesOnGrid);
}

void CubeFileWriter::writeCube(
    const vector<string>& filenames,
    shared_ptr<const Geometry> geometry,
    std::function<Eigen::MatrixXd(shared_ptr<GridController>)> calcProperties) {
  /*
   * write headers
   */
  shared_ptr<GridController> cubeGridController = this->writeHeaderAndCreateCubeGrid(
      filenames, geometry);
  /*
   * calc data
   */
  auto data = calcProperties(cubeGridController);
  assert(data.size() == filenames.size());
  /*
   * write data
   */
  for (unsigned int i=0; i<filenames.size(); ++i) {
    this->writeData(filenames[i], data.col(i));
  }
}

void CubeFileWriter::writeCube(
    const string& filename,
    shared_ptr<const Geometry> geometry,
    std::function<Eigen::VectorXd(shared_ptr<GridController>)> calcProperty) {
  /*
   * write header
   */
  shared_ptr<GridController> cubeGridController = this->writeHeaderAndCreateCubeGrid(
      filename, geometry);
  /*
   * calc data
   */

  Eigen::VectorXd data = calcProperty(cubeGridController);
  /*
   * write data
   */
  this->writeData(filename, data);
}

void CubeFileWriter::writeCube(
    const vector<string>& filenames,
    shared_ptr<const Geometry> geometry,
    std::function<vector<Eigen::VectorXd>(shared_ptr<GridController>)> calcProperties) {
  /*
   * write headers
   */
  shared_ptr<GridController> cubeGridController = this->writeHeaderAndCreateCubeGrid(
      filenames, geometry);
  /*
   * calc data
   */
  auto data = calcProperties(cubeGridController);
  assert(data.size() == filenames.size());
  /*
   * write data
   */
  for (unsigned int i=0; i<filenames.size(); ++i) {
    this->writeData(filenames[i], data[i]);
  }
}

shared_ptr<GridController> CubeFileWriter::writeHeaderAndCreateCubeGrid(
    string filename,
    shared_ptr<const Geometry> geometry) {
  return writeHeaderAndCreateCubeGrid( vector<string>{filename}, geometry);
}

shared_ptr<GridController> CubeFileWriter::writeHeaderAndCreateCubeGrid(
    vector<string> filenames,
    shared_ptr<const Geometry> geometry) {
  /*
   * Create grid
   */
  const double minX = geometry->getMinX() - _borderWidth;
  const double minY = geometry->getMinY() - _borderWidth;
  const double minZ = geometry->getMinZ() - _borderWidth;
  const double maxX = geometry->getMaxX() + _borderWidth;
  const double maxY = geometry->getMaxY() + _borderWidth;
  const double maxZ = geometry->getMaxZ() + _borderWidth;
  const double stepX = _stepSize[0];
  const double stepY = _stepSize[1];
  const double stepZ = _stepSize[2];
  const unsigned int nPointsX = (maxX - minX) / stepX;
  const unsigned int nPointsY = (maxY - minY) / stepY;
  const unsigned int nPointsZ = (maxZ - minZ) / stepZ;
  const unsigned int nPointsTotal = nPointsX * nPointsY * nPointsZ;
  std::unique_ptr<Eigen::Matrix3Xd > points(new Eigen::Matrix3Xd(3,nPointsTotal));
  std::unique_ptr<Eigen::VectorXd> weights(new Eigen::VectorXd(nPointsTotal));
  (*weights) = Eigen::VectorXd::Constant(nPointsTotal,1.0 / (double)nPointsTotal);
  double x = minX;
  unsigned long long npt=0;
  for (unsigned int ix = 0; ix < nPointsX;
      ++ix, x += stepX * _xUnitVector[0] + stepY * _yUnitVector[0]
                 + stepZ * _zUnitVector[0]) {
    double y = minY;
    for (unsigned int iy = 0; iy < nPointsY;
        ++iy, y += stepX * _xUnitVector[1] + stepY * _yUnitVector[1]
                   + stepZ * _zUnitVector[1]) {
      double z = minZ;
      for (unsigned int iz = 0; iz < nPointsZ; npt++,
          ++iz, z += stepX * _xUnitVector[2] + stepY * _yUnitVector[2]
                     + stepZ * _zUnitVector[2]) {
        (*points)(0,npt) = x;
        (*points)(1,npt) = y;
        (*points)(2,npt) = z;
      }
    }
  }
  /*
   * Print headers
   */
  for (const auto& filename : filenames) {
    string fullFileName = filename + ".cube";
    FILE* file = fopen(fullFileName.data(), "w");
    fprintf(file, "%s\n\n", "CUBE FILE");
    const unsigned int nAtoms = geometry->getNAtoms();
    fprintf(file, "%u %f %f %f\n", nAtoms, minX, minY, minZ);
    fprintf(
        file, "%u %f %f %f\n", nPointsX, stepX * _xUnitVector[0], stepX * _xUnitVector[1],
        stepX * _xUnitVector[2]);
    fprintf(
        file, "%u %f %f %f\n", nPointsY, stepY * _yUnitVector[0], stepY * _yUnitVector[1],
        stepY * _yUnitVector[2]);
    fprintf(
        file, "%u %f %f %f\n", nPointsZ, stepZ * _zUnitVector[0], stepZ * _zUnitVector[1],
        stepZ * _zUnitVector[2]);
    const auto& atoms = geometry->getAtoms();
    for (unsigned int i = 0; i < nAtoms; ++i) {
      fprintf(
          file, "%u %f %f %f %f\n", atoms[i]->getNuclearCharge(), 0.0, atoms[i]->getX(),
          atoms[i]->getY(), atoms[i]->getZ());
    }
    fclose(file);
  }
  return make_shared<GridController>(unique_ptr<Grid>(new Grid(move(points), move(weights))));
}

void CubeFileWriter::writeData(string filename, const Eigen::VectorXd& data) {
  /*
   * Print data
   */
  string fullFileName = filename + ".cube";
  FILE* file = fopen(fullFileName.data(), "a");
  for (unsigned int i = 0; i < data.size(); i++) {
    fprintf(file, "%g ", data[i]);
    if (i % 6 == 5) fprintf(file, "\n");
  }
  fclose(file);
}

}
