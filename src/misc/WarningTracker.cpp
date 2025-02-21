/**
 * @file WarningTracker.cpp
 *
 * @date Mar 9, 2017
 * @author David Schnieders
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
#include "misc/WarningTracker.h"
/* Include Std and External Headers */
#include <fstream>
#include <iostream>

namespace Serenity {

std::unique_ptr<WarningTracker> WarningTracker::_instance;
unsigned int WarningTracker::_taskNmbr = 0;
unsigned int WarningTracker::_nWarnings = 0;

WarningTracker& WarningTracker::getInstance() {
  if (!_instance) {
    _instance.reset(new WarningTracker());
  };
  return *_instance;
}

void WarningTracker::printWarning(std::string text, bool printToOutput) {
  // Print warning to screen
  if (printToOutput)
    std::cout << text << std::endl;

  // Print warning to file
  std::ofstream file;
  file.open("WARNING", std::ofstream::out | std::ofstream::app);
  file << std::endl;
  file << "In Task: " << std::to_string(_taskNmbr) << std::endl;
  file << std::endl;
  file << text << std::endl;
  file << std::endl;
  file << std::endl;
  file.close();

  _nWarnings++;
}

} /* namespace Serenity */
