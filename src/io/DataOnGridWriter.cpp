/**
 * @file   DataOnGridWriter.cpp
 *
 * @date   Apr 29, 2019
 * @author Jan Unsleber, Anja Massolle
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
#include "io/DataOnGridWriter.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"                          //Basis controller definition
#include "data/grid/BasisFunctionOnGridControllerFactory.h" //Basis function values.
#include "data/grid/MOCalculator.h"                         //MO on grid evaluation.
#include "data/grid/MatrixOperatorToGridTransformer.h"      //Density on grid evaluation.
#include "data/matrices/MatrixInBasis.h"                    //MatrixInBasis definition.
#include "geometry/Geometry.h"                              //Geometry definition.
#include "grid/GridController.h"                            //Access to the grid for the evaluations.
#include "parameters/Constants.h"

namespace Serenity {

DataOnGridWriter::DataOnGridWriter(const Settings& settings) : _settings(settings) {
}

void DataOnGridWriter::writeMatrixToGrid(const std::string& filename, std::shared_ptr<const Geometry> geometry,
                                         const MatrixInBasis<RESTRICTED>& matrix) {
  auto basisController = matrix.getBasisController();
  this->writeMatrixToGrid(filename, geometry, basisController, matrix);
}

void DataOnGridWriter::writeMatrixToGrid(const std::string& filename, std::shared_ptr<const Geometry> geometry,
                                         std::shared_ptr<BasisController> basisController, const Eigen::MatrixXd& matrix) {
  /*
   * define lambda function
   */
  auto calcValuesOnGrid = [&basisController, &matrix, this](std::shared_ptr<GridController> gridController) -> Eigen::VectorXd {
    auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(_settings, basisController, gridController);
    Eigen::VectorXd result;
    MatrixInBasis<RESTRICTED> tmp(basisController);
    tmp = matrix;
    MatrixOperatorToGridTransformer::transform(tmp, result, *basFuncOnGridController);
    return result;
  };
  // write out the data
  writeFile(filename, geometry, calcValuesOnGrid);
}

void DataOnGridWriter::writeVectorSetToGrid(const std::vector<std::string>& filenames, std::shared_ptr<const Geometry> geometry,
                                            std::shared_ptr<BasisController> basisController, Eigen::MatrixXd& inVectorSet) {
  assert((long int)(filenames.size()) == inVectorSet.cols());
  /*
   * define lambda function
   */
  auto calcValuesOnGrid = [&basisController, &inVectorSet, this](std::shared_ptr<GridController> gridController) -> Eigen::MatrixXd {
    // get useful data
    auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(_settings, basisController, gridController);

    // use MOCalculator to calculate MOs on grid
    MOCalculator moCalc(basFuncOnGridController);
    Eigen::MatrixXd valuesOnGrid = moCalc.calcMOValuesOnGrid(inVectorSet);

    return valuesOnGrid;
  };
  // write out the data
  writeFile(filenames, geometry, calcValuesOnGrid);
}

void DataOnGridWriter::writeFile(const std::string& filename, std::shared_ptr<const Geometry> geometry,
                                 std::function<Eigen::VectorXd(std::shared_ptr<GridController>)> calcProperty) {
  /*
   * write header and generate Grid
   */
  std::shared_ptr<GridController> plotGridController = this->writeHeaderAndCreateGrid(filename, geometry);
  /*
   * calc data
   */
  Eigen::VectorXd data = calcProperty(plotGridController);
  /*
   * write data
   */
  this->writeData(filename, data);
} /* writeFile */

void DataOnGridWriter::writeFile(const std::vector<std::string>& filenames, std::shared_ptr<const Geometry> geometry,
                                 std::function<Eigen::MatrixXd(std::shared_ptr<GridController>)> calcProperties) {
  /*
   * write header and generate Grid
   */
  std::shared_ptr<GridController> plotGridController = this->writeHeaderAndCreateGrid(filenames, geometry);
  /*
   * calc data
   */
  auto data = calcProperties(plotGridController);
  /*
   * write data
   */
  for (unsigned int i = 0; i < filenames.size(); ++i) {
    this->writeData(filenames[i], data.col(i));
  }
} /* writeFile */

Eigen::Vector3d DataOnGridWriter::convertPointToBohr(const std::vector<double> point) {
  Eigen::Vector3d pointBohr;
  for (unsigned int i = 0; i < 3; i++) {
    pointBohr[i] = point[i] * ANGSTROM_TO_BOHR;
  }
  return pointBohr;
}; /* convertPointToBohr */

} /*namespace Serenity*/
