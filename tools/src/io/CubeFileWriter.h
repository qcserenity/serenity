/**
 * @file   CubeFileWriter.h
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
#ifndef CUBEFILEWRITER_H_
#define CUBEFILEWRITER_H_
/* Include Serenity Internal Headers */
#include "io/DataOnGridWriter.h" //DataOnGridWriter (base class)
#include "tasks/PlotTask.h"      //PlotTaskSettings

namespace Serenity {

/* Forward declarations */
struct Settings;

/**
 * @class CubeFileWriter CubeFileWriter.h
 *
 * @brief Prints data into a cube file.
 *
 * This class provides the methods to write properties into a cube file.
 * The properties can be expanded in the AO basis or they can be on a
 * different grid.
 */
class CubeFileWriter : public DataOnGridWriter {
 public:
  /**
   * @brief Constructor.
   *
   * @param settings The settings of the system.
   * @param cubeGridSettings The settings of the PlotTask.
   */
  CubeFileWriter(const Settings& settings, const PlotTaskSettings& cubeGridSettings);
  /**
   * @brief default destructor.
   */
  virtual ~CubeFileWriter() = default;

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
  void writeData(std::string filename, const Eigen::VectorXd& data) override final;
  /*
   * Variables
   */
  ///@brief settings of the PlotTask
  const PlotTaskSettings& _cubeGridSettings;
  ///@brief x unit vector of the cube grid
  Eigen::Vector3d _xUnitVector;
  ///@brief y unit vector of the cube grid
  Eigen::Vector3d _yUnitVector;
  ///@brief z unit vector of the cube grid
  Eigen::Vector3d _zUnitVector;
  ///@brief step size in x, y and z direction of the cube grid
  Eigen::Vector3d _gridSpacing;
  ///@brief additional width around the molecule
  double _borderWidth;
};

} // namespace Serenity
#endif /* CUBEFILEWRITER_H_ */
