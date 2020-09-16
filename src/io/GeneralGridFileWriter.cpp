/**
 * @file GeneralGridFileWriter.cpp
 *
 * @author Moritz Bensberg
 * @date May 19, 2020
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
#include "io/GeneralGridFileWriter.h"
/* Include Serenity Internal Headers */
#include "grid/GridController.h" //GridController
#include "parameters/Constants.h"

namespace Serenity {

GeneralGridFileWriter::GeneralGridFileWriter(const Settings& settings, std::shared_ptr<GridController> gridController)
  : DataOnGridWriter(settings), _gridController(gridController) {
}

std::shared_ptr<GridController> GeneralGridFileWriter::writeHeaderAndCreateGrid(std::string filename,
                                                                                std::shared_ptr<const Geometry> geometry) {
  return writeHeaderAndCreateGrid(std::vector<std::string>{filename}, geometry);
}

std::shared_ptr<GridController> GeneralGridFileWriter::writeHeaderAndCreateGrid(std::vector<std::string> filenames,
                                                                                std::shared_ptr<const Geometry> geometry) {
  (void)geometry;
  for (const auto& filename : filenames) {
    std::string fullFileName = filename + ".grid.data.xyz";
    std::FILE* file = std::fopen(fullFileName.data(), "w");
    fprintf(file, "%10d\n", _gridController->getNGridPoints());
    fprintf(file, "%10s %10s %10s %10s %10s\n", "Label ", "x ", "y ", "z ", "value");
    fclose(file);
  }
  return _gridController;
}

void GeneralGridFileWriter::writeData(std::string filename, const Eigen::VectorXd& data) {
  std::string fullFileName = filename + ".grid.data.xyz";
  std::FILE* file = std::fopen(fullFileName.data(), "a");
  const Eigen::MatrixXd coordinates = BOHR_TO_ANGSTROM * _gridController->getGridPoints();
  for (unsigned int i = 0; i < data.size(); i++) {
    fprintf(file, "%3s %10g %10g %10g %10g\n", "P ", coordinates(0, i), coordinates(1, i), coordinates(2, i), data[i]);
  }
  fclose(file);
}

} /* namespace Serenity */
