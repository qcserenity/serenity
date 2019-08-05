/**
 * @file   CubeFileWriter.h
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
#ifndef CUBEFILEWRITER_H_
#define CUBEFILEWRITER_H_
/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "geometry/Point.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <string>
#include <vector>


namespace Serenity {
/* Forward declarations */
class BasisController;
class Geometry;
struct Settings;
class Grid;
class GridController;
template<Options::SCF_MODES>class MatrixInBasis;
/**
 * @class CubeFileWriter CubeFileWriter.h
 *
 * @brief Print things into a cube file
 *
 * This class provides the methods to write properties into a cube file.
 * The properties can be expressed in coefficients of basis functions
 * or they can be on a different grid.
 *
 */
class CubeFileWriter {

public:
  /**
   * @brief creates an instance of the CubeFileWriter
   *
   * @param settings The settings.
   */
  CubeFileWriter(const Settings& settings);

  virtual ~CubeFileWriter() = default;
  /**
   * @brief Prints sum of sandwiched matrix element into a cube file (e.g. electron density)
   *
   * This method prints the sum of all products of two basis functions and their
   * corresponding matrix element onto a grid in a cube file
   * \f$\sum_{a,b} <a(r)|M_{a,b}|{b(r)}>\f$
   * An example for this is the electron density. Note that this is *not* a reasonable thing to do
   * for the matrix representation of *any* operator.
   *
   * @param filename the file name for the generated cube file
   * @param geometry the geometry that shall be printed into the header
   * @param matrix the matrix holding the coefficients of the property that shall be printed
   */
  void writeMatrixToCube(
      const std::string& filename,
      std::shared_ptr<const Geometry> geometry,
      const MatrixInBasis<RESTRICTED>& matrix);
  /**
   * @brief Prints sum of sandwiched matrix element into a cube file (e.g. electron density)
   *
   * This method prints the sum of all products of two basis functions and their
   * corresponding matrix element onto a grid in a cube file
   * \f$\sum_{a,b} <a(r)|M_{a,b}|{b(r)}>\f$
   * An example for this is the electron density. Note that this is *not* a reasonable thing to do
   * for the matrix representation of *any* operator.
   *
   * @param filename the file name for the generated cube file
   * @param geometry the geometry that shall be printed into the header
   * @param basis the basis for the matrix
   * @param matrix the matrix holding the coefficients of the property that shall be printed
   */
  void writeMatrixToCube(
      const std::string& filename,
      std::shared_ptr<const Geometry> geometry,
      std::shared_ptr<BasisController> basis,
      const Eigen::MatrixXd& matrix);
  /**
   * @brief See the other method, extension for multiple matrices.
   */
  void writeMatrixToCube(
      const std::vector<std::string>& filenames,
      std::shared_ptr<const Geometry> geometry,
      const std::vector<std::reference_wrapper<const MatrixInBasis<RESTRICTED> > >& matrices);
  /**
   * @brief Prints dot product of inVector and basis into a cube file (e.g. MOs)
   *
   * This method prints the sum of all products of a basis functions its
   * corresponding coefficient in the vector onto a grid in a cube file
   * \f$\sum_{a} Vec_{a}|a(r)>\f$
   * An example for this are the MOs.
   *
   * @param filename the file name for the generated cube file
   * @param geometry the geometry that shall be printed into the header
   * @param basis the corresponding basis for the coefficients in the vector
   * @param inVector the vector holding the coefficients of the property that shall be printed
   */
  void writeVectorToCube(
      const std::string& filename,
      std::shared_ptr<const Geometry> geometry,
      std::shared_ptr<BasisController> basisController,
      Eigen::VectorXd& inVector);

  /*
   * Getters and setters
   */

  /**
   * @brief Returns the border width around the geometry for which the cube grid is created.
   *
   * Returns the border width around the geometry for which the cube grid is created.
   * The default is 5.0 bohr
   *
   * @return grid border width
   */
  inline double getBorderWidth() const {
    return _borderWidth;
  }
  /**
   * @brief Sets the border width around the geometry for which the cube grid is created.
   *
   * Sets the border width around the geometry for which the cube grid is created.
   * The default is 5.0 bohr. 
   *
   * @param borderWidth grid border width
   */
  inline void setBorderWidth(double borderWidth) {
    _borderWidth = borderWidth;
  }
  /**
   * @brief Returns the origin for the cube grid.
   *
   * Returns the origin for the cube grid.
   * The default is (0.0,0.0,0.0)
   *
   * @return origin for the cube grid.
   */
  inline const Point& getOrigin() const {
    return _origin;
  }
  /**
   * @brief Sets the origin for the cube grid.
   *
   * Sets the origin for the cube grid.
   * The default is (0.0,0.0,0.0)
   *
   * @param origin origin for the cube grid
   */
  inline void setOrigin(const Point& origin ) {
    _origin = origin;
  }
  /**
   * @brief Returns the step sizes for the three vectors spanning the cube mesh.
   *
   * Returns the step sizes for the three vectors spanning the cube mesh.
   * The default is (0.25,0.25,0.25)
   *
   * @return step sizes
   */
  inline Eigen::Vector3d getStepSize() const {
    return _stepSize;
  }
  /**
   * @brief Sets the step sizes for the three vectors spanning the cube mesh.
   *
   * Sets the step sizes for the three vectors spanning the cube mesh.
   * The default is (0.25,0.25,0.25)
   *
   * @param stepSize
   * TODO see above
   */
  inline void setStepSize(const Eigen::Vector3d& stepSize) {
    _stepSize = stepSize;
  }

  /**
   * @brief Returns the the first (x) vector used to define the cube grid.
   *
   * Returns the the first (x) vector used to define the cube grid.
   * The default is (1.0,0.0,0.0)
   *
   * @return x unit vector
   */
  inline const Eigen::Vector3d& getXUnitVector() const {
    return _xUnitVector;
  }
  /**
   * @brief Sets the the first (x) vector used to define the cube grid.
   *
   * Sets the the first (x) vector used to define the cube grid.
   * The default is (1.0,0.0,0.0)
   *
   * @return unitVector
   */
  inline void setXUnitVector(const Eigen::Vector3d& unitVector) {
    _xUnitVector = unitVector;
  }
  /**
   * @brief Returns the the second (y) vector used to define the cube grid.
   *
   * Returns the the second (y) vector used to define the cube grid.
   * The default is (0.0,1.0,0.0)
   *
   * @return y unit vector
   */
  inline const Eigen::Vector3d& getYUnitVector() const {
    return _yUnitVector;
  }
  /**
   * @brief Sets the the second (y) vector used to define the cube grid.
   *
   * Sets the the second (y) vector used to define the cube grid.
   * The default is (0.0,1.0,0.0)
   *
   * @return unitVector
   */
  inline void setYUnitVector(const Eigen::Vector3d& unitVector) {
    _yUnitVector = unitVector;
  }
  /**
   * @brief Returns the the third (z) vector used to define the cube grid.
   *
   * Returns the the third (z) vector used to define the cube grid.
   * The default is (0.0,0.0,1.0)
   *
   * @return z unit vector
   */
  inline const Eigen::Vector3d& getZUnitVector() const {
    return _zUnitVector;
  }
  /**
   * @brief Sets the the third (z) vector used to define the cube grid.
   *
   * Sets the the third (z) vector used to define the cube grid.
   * The default is (0.0,0.0,1.0)
   *
   * @return unitVector
   */
  inline void setZUnitVector(const Eigen::Vector3d& unitVector) {
    _zUnitVector = unitVector;
  }
  /*
   * Functions
   */
  void writeCube(
      const std::string& filename,
      std::shared_ptr<const Geometry> geometry,
      std::function<Eigen::VectorXd(std::shared_ptr<GridController>)> calcPropetry);
  void writeCube(
      const std::vector<std::string>& filenames,
      std::shared_ptr<const Geometry> geometry,
      std::function<std::vector<Eigen::VectorXd>(std::shared_ptr<GridController>)>
        calcPropetries);

private:

  std::shared_ptr<GridController> writeHeaderAndCreateCubeGrid(
      std::string filename,
      std::shared_ptr<const Geometry> geometry);

  std::shared_ptr<GridController> writeHeaderAndCreateCubeGrid(
      std::vector<std::string> filenames,
      std::shared_ptr<const Geometry> geometry);
  
  void writeData(std::string filename, const Eigen::VectorXd& data);
  /*
   * Variables
   */
  Point _origin = Point(0.0, 0.0, 0.0);
  Eigen::Vector3d _xUnitVector;
  Eigen::Vector3d _yUnitVector;
  Eigen::Vector3d _zUnitVector;
  Eigen::Vector3d _stepSize;
  double _borderWidth = 5.0;
  const Settings& _settings;
};

}
#endif /* CUBEFILEWRITER_H_ */
