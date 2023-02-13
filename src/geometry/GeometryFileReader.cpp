/**
 * @file   GeometryFileReader.cpp
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
/* Include Class Header*/
#include "geometry/GeometryFileReader.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
/* Include Std and External Headers */
#include <cassert>
#include <stdexcept>

namespace Serenity {

GeometryFileReader::GeometryFileReader(std::string filePath) : _filePath(filePath) {
}

GeometryFileReader::~GeometryFileReader() {
  assert(!_file.is_open());
}

void GeometryFileReader::loadFile(std::string filePath) {
  assert(!_file.is_open());
  _file.open(filePath.c_str(), std::ifstream::in);
  if (!_file.good()) {
    throw SerenityError((std::string) "Error while loading file: " + filePath + ". Does it exist?");
  }
}

std::unique_ptr<Geometry> GeometryFileReader::readGeometry() {
  loadFile(_filePath);
  auto result = readGeometry(0);
  _file.close();
  return result;
}

std::unique_ptr<Geometry> GeometryFileReader::readGeometry(std::string filePath) {
  loadFile(filePath);
  auto result = readGeometry(0);
  _file.close();
  return result;
}

std::unique_ptr<Geometry> GeometryFileReader::readGeometry(std::string filePath, unsigned int index) {
  loadFile(filePath);
  auto result = readGeometry(index);
  _file.close();
  return result;
}

std::vector<std::unique_ptr<Geometry>> GeometryFileReader::readAllGeometries() {
  const bool fileWasOpen = _file.is_open();
  if (!fileWasOpen)
    loadFile(_filePath);
  std::vector<std::unique_ptr<Geometry>> geometries;
  unsigned int index = 0;
  while (true) {
    auto newGeometry = readGeometry(index++);
    if (!newGeometry) {
      break;
    }
    else {
      geometries.push_back(std::move(newGeometry));
    }
  }
  if (!fileWasOpen)
    _file.close();
  return geometries;
}

std::unique_ptr<Geometry> GeometryFileReader::readGeometry(unsigned int index) {
  const bool fileWasOpen = _file.is_open();
  if (!fileWasOpen)
    loadFile(_filePath);
  auto result = readGeometryFromLoadedFile(index);
  if (!fileWasOpen)
    _file.close();
  return result;
}

std::vector<std::unique_ptr<Geometry>> GeometryFileReader::readAllGeometries(std::string filePath) {
  loadFile(filePath);
  auto result = readAllGeometries();
  _file.close();
  return result;
}

} /* namespace Serenity */
