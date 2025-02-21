/**
 * @file ExternalChargeReader.h
 *
 * @date Apr. 29, 2024
 * @author Moritz Bensberg
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

#ifndef SERENITY_EXTERNALCHARGEREADER_H
#define SERENITY_EXTERNALCHARGEREADER_H

/* Include Std and External Headers */
#include <memory>
#include <string>
#include <vector>

namespace Serenity {
class Point;

/**
 * @class ExternalChargeReader ExternalChargeReader.h
 * @brief Static class that provides read functions for external charge files.
 */
class ExternalChargeReader {
 private:
  ExternalChargeReader() = delete;

 public:
  /**
   * @brief Read an external charge file.
   *
   * The file format is expected to be of the following kind:
   *
   * <Number of charges>
   * <charge> <x-coord> <y-coord> <z-coord>
   * <charge> <x-coord> <y-coord> <z-coord>
   * ...
   *
   * For instance:
   *
   * 16808
   * 0.2943 1.79625 33.8873 1.77406
   * 0.1642 1.05325 33.8222 1.09206
   * 0.1642 2.69725 34.0243 1.33906
   * 0.1642 1.97625 33.0262 2.26906
   *
   * @param filePath The path to the file.
   * @return The external charges.
   */
  static std::vector<std::pair<double, Point>> readExternalChargeFile(const std::string& filePath);
};

} // namespace Serenity

#endif // SERENITY_EXTERNALCHARGEREADER_H
