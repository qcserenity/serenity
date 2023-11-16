/**
 * @file GeneralGridFileWriter.h
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

#ifndef IO_GENERALGRIDFILEWRITER_H_
#define IO_GENERALGRIDFILEWRITER_H_
/* Include Serenity Internal Headers */
#include "io/DataOnGridWriter.h" //DataOnGridWriter (base class)
/* Include Std and External Headers */
#include <memory>
#include <string>
#include <vector>

namespace Serenity {

class GridController;
struct Settings;
/**
 * @class
 * @brief A class that writes values on any grid to a text file.
 */
class GeneralGridFileWriter : public DataOnGridWriter {
 public:
  /**
   * @brief Constructor.
   * @param settings       The settings.
   * @param gridController The grid controller.
   */
  GeneralGridFileWriter(const Settings& settings, std::shared_ptr<GridController> gridController);
  /**
   * @brief Default destructor.
   */
  virtual ~GeneralGridFileWriter() = default;

 private:
  /**
   * @brief Writes the header of the cube file and generates the cube grid.
   *
   * @param filename The base name of the cube file.
   * @param geometry The geometry of the system printed to the cube file.
   *
   * @return The grid controller of the cube grid.
   */
  std::shared_ptr<GridController> writeHeaderAndCreateGrid(std::string filename,
                                                           std::shared_ptr<const Geometry> geometry) override final;
  /**
   * @brief Writes the headers of multiple cube files and generates the cube
   * grid.
   *
   * @param filenames The base names of the cube files.
   * @param geometry The geometry of the system printed to the cube file.
   *
   * @return The grid controller of the cubic grid.
   */
  std::shared_ptr<GridController> writeHeaderAndCreateGrid(std::vector<std::string> filenames,
                                                           std::shared_ptr<const Geometry> geometry) override final;
  /**
   * @brief Writes the data to the cube file
   *
   * @param filename The base name of the cube file.
   * @param data The data which shall be printed to the cube file.
   */
  void writeData(std::string filename, const Eigen::VectorXd& data) override;
  ///@brief The grid controller.
  std::shared_ptr<GridController> _gridController;
};

} /* namespace Serenity */

#endif /* IO_GENERALGRIDFILEWRITER_H_ */
