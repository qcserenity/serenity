/**
 * @file   GeometryFileReader.h
 *
 * @date   Mar 23, 2013
 * @author Thomas Dresselhaus
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
#ifndef GEOMETRYFILEREADER_H_
#define GEOMETRYFILEREADER_H_
/* Include Std and External Headers */
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace Serenity {
/* Forward declarations */
class Geometry;
/**
 * @class GeometryFileReader GeometryFileReader.h
 * @brief Abstract file reading.
 */
class GeometryFileReader {
 public:
  GeometryFileReader() = default;
  /**
   * @param filePath inside the file with this path the geometry resides in some form.
   */
  GeometryFileReader(std::string filePath);
  virtual ~GeometryFileReader();
  /**
   * @brief Loads a new file into this object (as an ifstream).
   *
   * This method may be overloaded if the file is treated in a non-standard way.
   *
   * @param filePath The full (relative) path to the file.
   */
  virtual void loadFile(std::string filePath);
  /**
   * @returns A new instance of a Geometry based on the contents of the loaded file. If more
   *          than one geometry is available in the file, the first one will be used.
   */
  std::unique_ptr<Geometry> readGeometry();
  /**
   * @param  index The number indicating which geometry out of a file should be read. Starts
   *               at zero.
   * @returns A new instance of a Geometry based on the contents of the loaded file iff the
   *          index is valid, NULL otherwise.
   *
   * This method may be overloaded iff the file is treated in a non-standard way.
   */
  virtual std::unique_ptr<Geometry> readGeometry(unsigned int index);
  /**
   * @param filePath The full (relative) path to the file from which the Geometry will be read.
   */
  std::unique_ptr<Geometry> readGeometry(std::string filePath);
  /**
   * @param filePath The full (relative) path to the file from which the Geometry will be read.
   * @param index    The number indicating which geometry out of the file should be read. Starts
   *                 at zero.
   * @returns A new instance of a Geometry based on the contents of the loaded file iff the
   *          index is valid, NULL otherwise.
   */
  std::unique_ptr<Geometry> readGeometry(std::string filePath, unsigned int index);
  /**
   * @returns New instances for all geometries present in the currently loaded file.
   */
  std::vector<std::unique_ptr<Geometry>> readAllGeometries();
  /**
   * @param   filePath The full (relative) path to the file from which the geometries will be read.
   * @returns New instances for all geometries present in the specified file.
   */
  std::vector<std::unique_ptr<Geometry>> readAllGeometries(std::string filePath);

 protected:
  /**
   * @brief The actually working method.
   *
   * This method must be overloaded in an actual implementation.
   *
   * @param  index The number indicating which geometry out of a file should be read. Starts
   *               at zero.
   * @returns A new instance of a Geometry based on the contents of the loaded file iff the
   *          index is valid, NULL otherwise.
   */
  virtual std::unique_ptr<Geometry> readGeometryFromLoadedFile(unsigned int index) = 0;
  /**
   * The actual content of the file.
   */
  std::ifstream _file;
  /**
   * The path to the file.
   */
  std::string _filePath;
};

} /* namespace Serenity */
#endif /* GEOMETRYFILEREADER_H_ */
