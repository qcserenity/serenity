/**
 * @file   PlaneFileWriter.h
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
#ifndef PLANEFILEWRITER_H_
#define PLANEFILEWRITER_H_

/* Include Serenity Internal Headers */
#include "io/DataOnGridWriter.h" //DataOnGridWriter (base class)
#include "tasks/PlotTask.h"      //PlotTaskSettings

namespace Serenity {

/* Forward declarations */
struct Settings;

/**
 * @class PlaneFileWriter PlaneFileWriter.h
 *
 * @brief Print things into a dat file
 *
 * This class provides the methods to write properties on a plane grid into a\n
 * dat file. The properties can be expressed in coefficients of basis functions\n
 * or they can be on a different grid.
 *
 */
class PlaneFileWriter : public DataOnGridWriter {
 public:
  /**
   * @brief creates an instance of the PlaneFileWriter
   *
   * @param settings The settings of the system.
   * @param planeGridSettings The settings of the PlotTask.
   */
  PlaneFileWriter(const Settings& settings, const PlotTaskSettings& planeGridSettings);
  /**
   * @brief Destructor
   */
  virtual ~PlaneFileWriter() = default;

 private:
  /**
   * @brief Creates the plane grid from the three points which defines the plane
   * and writes the header for the dat file.
   *
   * @param filename basename of the dat file
   * @param geometry the geometry of the system printed to the dat file
   * @return returns the grid controller from the plane grid
   */
  std::shared_ptr<GridController> writeHeaderAndCreateGrid(std::string filename,
                                                           std::shared_ptr<const Geometry> geometry) override final;

  /**
   * @brief Creates the plane grid from the three points which defines the plane
   * and writes the header for the dat files.
   *
   * @param filenames basenames of the dat files
   * @param geometry the geometry of the system printed to the dat files
   * @return returns the grid controller from the plane grid
   */
  std::shared_ptr<GridController> writeHeaderAndCreateGrid(std::vector<std::string> filenames,
                                                           std::shared_ptr<const Geometry> geometry) override final;
  /**
   * @brief Writes the data calculated on the plane grid to a dat file (file
   * format: x, y, z, value)
   *
   * @param filename the basename of the file
   * @param data the data calculated on these gridpoints
   */
  void writeData(std::string filename, const Eigen::VectorXd& data) override final;

  ///@brief the settings of the PlotTask
  const PlotTaskSettings& _planeGridSettings;
  ///@brief the first unitvector of the plane
  Eigen::Vector3d _xUnitVector;
  ///@brief the second unitvector of the plane
  Eigen::Vector3d _yUnitVector;
  ///@brief the normal vector of the plane
  Eigen::Vector3d _normalVector;
  ///@brief the rotation matrix to rotate the plane in the xy plane
  Eigen::Matrix3d _rotationMatrix;
  ///@brief the step size to generate the plane grid along the two unit vectors
  Eigen::Vector2d _gridSpacing;
  ///@brief number of grid points along the first unit vector
  unsigned int _nPointsX;
  ///@brief number of grid points along the second unit vector
  unsigned int _nPointsY;
  ///@brief the three points which defines the plane
  std::vector<Eigen::Vector3d> _pointsOnPlane;
  ///@brief the grid Controller of the plane grid
  std::shared_ptr<GridController> _planeGridController;
  ///@brief the atoms of the molecule which were projected on the plane
  Eigen::Matrix3Xd _projectedPointsOnPlane;
  ///@brief the border with around the geometry
  double _borderWidth;
  ///@brief the distance from the plane in which an atom has to be, to be
  /// projected on the plane
  double _projectCutOffRadius;
  ///@brief the geometry which shall be printed
  std::shared_ptr<const Geometry> _geometry;

}; /* class PlaneFileWriter */

} /*namespace Serenity*/
#endif /* PLANEFILEWRITER_H_ */
