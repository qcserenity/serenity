/**
 * @file   DataOnGridWriter.h
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
#ifndef DATAONGRIDWRITER_H_
#define DATAONGRIDWRITER_H_
/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h" //RESTRICTED/UNRESTRICTED
/* Include Std and External Headers */
#include <Eigen/Dense> // MatrixXd
#include <memory>      // std::shared_ptr
#include <vector>      //std::vector

namespace Serenity {
/* Forward declarations */
class BasisController;
class Geometry;
struct Settings;
class GridController;
template<Options::SCF_MODES>
class MatrixInBasis;
/**
 * @class DataOnGridWriter DataOnGridWriter.h
 *
 * @brief Base class for printing data on grids
 *
 * This class provides methods to write properties onto a grid for plotting.
 * The properties can be expressed in coefficients of basis functions
 * or they can be defined on a different grid.
 *
 */
class DataOnGridWriter {
 public:
  /**
   * @brief creates an instance of the DataOnGridWriter
   *
   * @param settings The settings of the system.
   */
  DataOnGridWriter(const Settings& settings);

  /**
   * @brief default destructor
   */
  virtual ~DataOnGridWriter() = default;
  /**
   * @brief Prints sum of sandwiched matrix element on the grid into a file
   * (e.g. electron density)
   *
   * This method prints the sum of all products of two basis functions and their
   * corresponding matrix element onto a grid in a file
   * \f$\sum_{a,b} M_{a,b} \cdot a(r) \cdot b(r) \f$
   * An example for this is the electron density. Note that this is *not* a
   * reasonable thing to do for the matrix representation of *any* operator.
   *
   * @param filename the file name for the generated file
   * @param geometry the geometry that shall be printed
   * @param matrix the matrix holding the coefficients of the property that
   * shall be printed
   */
  void writeMatrixToGrid(const std::string& filename, std::shared_ptr<const Geometry> geometry,
                         const MatrixInBasis<RESTRICTED>& matrix);
  /**
   * @brief Prints sum of sandwiched matrix element on a grid into a file (e.g.
   * electron density)
   *
   * This method prints the sum of all products of two basis functions and their
   * corresponding matrix element onto a grid in a file
   * \f$\sum_{a,b} M_{a,b} \cdot a(r) \cdot b(r) \f$
   * An example for this is the electron density. Note that this is *not* a
   * reasonable thing to do for the matrix representation of *any* operator.
   *
   * @param filename the file name for the generated file
   * @param geometry the geometry that shall be printed
   * @param basis the basis for the matrix
   * @param matrix the matrix holding the coefficients of the property that
   * shall be printed
   */
  void writeMatrixToGrid(const std::string& filename, std::shared_ptr<const Geometry> geometry,
                         std::shared_ptr<BasisController> basis, const Eigen::MatrixXd& matrix);
  /**
   * @brief Prints dot product of inVectorSet and basis into a cube file (e.g. MOs)
   *
   * This method prints the sum of all products of a basis functions its
   * corresponding coefficient in the vectors onto a grid in a cube file
   * \f$\sum_{a} Vec_{a} \cdot a(r) \f$
   * An example for this are the MOs.
   *
   * @param filenames the file names for the generated cube files
   * @param geometry the geometry that shall be printed into the header
   * @param basis the corresponding basis for the coefficients in the vector
   * @param inVectorSet the vectors holding the coefficients of the property that shall be printed
   */
  void writeVectorSetToGrid(const std::vector<std::string>& filenames, std::shared_ptr<const Geometry> geometry,
                            std::shared_ptr<BasisController> basisController, Eigen::MatrixXd& inVectorSet);
  /**
   * @brief Generates the file (.cube or .dat) with the plotted data.
   *
   * @param filename the base name of the file
   * @param geometry the geometry that shall be printed
   * @param calcPropetry lambda function of the property that
   * shall be printed
   */
  void writeFile(const std::string& filename, std::shared_ptr<const Geometry> geometry,
                 std::function<Eigen::VectorXd(std::shared_ptr<GridController>)> calcProperty);

  /**
   * @brief Generates the file (.cube or .dat) with the plotted data.
   *
   * @param filenames the base names of the files
   * @param geometry the geometry that shall be printed
   * @param calcPropetry lambda function of the property that
   * shall be printed
   */
  void writeFile(const std::vector<std::string>& filenames, std::shared_ptr<const Geometry> geometry,
                 std::function<Eigen::MatrixXd(std::shared_ptr<GridController>)> calcProperties);

 protected:
  /**
   * @brief Converts a point in Angstrom into a point in Bohr, the data
   * type is also changed from a std::vector to an Eigen::Vector3d.
   *
   * @param point The point in Angstrom
   *
   * @return the point in Bohr
   */
  Eigen::Vector3d convertPointToBohr(const std::vector<double> point);

  /**
   * @brief writes the headers of a cube or dat file and generates the
   * corresponding grid.
   *
   * @param filename base name of the file
   * @param geometry geometry of the system printed to the file
   *
   * @return the grid controller of the grid
   */
  virtual std::shared_ptr<GridController> writeHeaderAndCreateGrid(std::string filename,
                                                                   std::shared_ptr<const Geometry> geometry) = 0;
  /**
   * @brief writes the headers of multiple cube or dat files and generates the
   * corresponding grid
   *
   * @param filenames base names of the files
   * @param geometry geometry of the system printed to the files
   *
   * @return the grid controller of the grid
   */
  virtual std::shared_ptr<GridController> writeHeaderAndCreateGrid(std::vector<std::string> filenames,
                                                                   std::shared_ptr<const Geometry> geometry) = 0;

  /**
   * @brief writes the data to the file
   *
   * @param filename base name of the file
   * @param data data which shall be printed to the file
   */
  virtual void writeData(std::string filename, const Eigen::VectorXd& data) = 0;

  ///@brief Settings of the system
  const Settings& _settings;
};

} /*namespace Serenity*/
#endif /* DATAONGRIDWRITER_H_ */
