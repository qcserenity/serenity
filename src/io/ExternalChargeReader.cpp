/**
 * @file ExternalChargeReader.cpp
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

/* Include Class Header*/
#include "io/ExternalChargeReader.h"
/* Include Serenity Internal Headers */
#include "geometry/Point.h"
#include "io/FormattedOutputStream.h"
#include "misc/SerenityError.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <fstream>

namespace Serenity {

std::vector<std::pair<double, Point>> ExternalChargeReader::readExternalChargeFile(const std::string& filePath) {
  OutputControl::nOut << "Reading external charges from file " << filePath << std::endl;
  std::ifstream input(filePath);
  if (!input.is_open()) {
    throw SerenityError("Unable to read point charge file " + filePath);
  }
  std::string line;
  std::string xCoordString;
  std::string yCoordString;
  std::string zCoordString;
  std::string chargeString;
  std::vector<std::pair<double, Point>> charges;
  std::getline(input, line); // Skip first line.
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream iss(line);
    try {
      iss >> chargeString;
      iss >> xCoordString;
      iss >> yCoordString;
      iss >> zCoordString;
      const double charge = std::stod(chargeString);
      const double xCoord = std::stod(xCoordString) * ANGSTROM_TO_BOHR;
      const double yCoord = std::stod(yCoordString) * ANGSTROM_TO_BOHR;
      const double zCoord = std::stod(zCoordString) * ANGSTROM_TO_BOHR;
      charges.emplace_back(charge, Point(xCoord, yCoord, zCoord));
    }
    catch (...) {
      SerenityError("Error: External charge file not formatted as expected. The format must be:\n"
                    "charge-value x-coord y-coord z-coord\n"
                    "All coordinates must be provided in Angstrom and all charges in atomic units.");
    };
  }
  return charges;
}

} // namespace Serenity