/**
 * @file   GeometryFactory.h
 *
 * @date   24.03.2013
 * @author Thomas Dresselhaus
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
#ifndef GEOMETRYFACTORY_H_
#define GEOMETRYFACTORY_H_

namespace Serenity {
/* Forward declarations */
class Geometry;
/**
 * @class GeometryFactory GeometryFactory.h
 * @brief Interface to create geometries. Only a marker interface!
 *
 * No method is specified here, because different types of arguments are expected.
 */
class GeometryFactory {
public:
  GeometryFactory() = default;
  virtual ~GeometryFactory() = default;
};

} /* namespace Serenity */
#endif /* GEOMETRYFACTORY_H_ */
