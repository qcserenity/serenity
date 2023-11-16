/**
 * @file   XyzFileToGeometryConverter.h
 *
 * @date   Mar 19, 2013
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
#ifndef XYZFILETOGEOMETRYCONVERTER_H_
#define XYZFILETOGEOMETRYCONVERTER_H_
/* Include Serenity Internal Headers */
#include "geometry/GeometryFileReader.h"

namespace Serenity {
/* Forward declarations */
/**
 * @class XyzFileToGeometryConverter XyzFileToGeometryConverter.h
 * @brief Generates a Geometry object out of an xyz file.
 *
 * Reads xyz-files and converts the coordinates from Angstrom to internally used atomic units.
 */
class XyzFileToGeometryConverter : public GeometryFileReader {
 public:
  XyzFileToGeometryConverter() = default;
  /**
   * @param filePath full (relative) path to an xyz file
   */
  XyzFileToGeometryConverter(std::string filePath);
  virtual ~XyzFileToGeometryConverter() = default;

  std::unique_ptr<Geometry> readGeometryFromLoadedFile(unsigned int index) override final;
};

} /* namespace Serenity */
#endif /* XYZFILETOGEOMETRYCONVERTER_H_ */
