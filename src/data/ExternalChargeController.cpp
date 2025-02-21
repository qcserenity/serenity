/**
 * @file ExternalChargeController.cpp
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
#include "ExternalChargeController.h"
/* Include Serenity Internal Headers */
#include "geometry/Point.h"
#include "geometry/efields/EFieldPlates.h"
#include "io/ExternalChargeReader.h"
#include "settings/Settings.h"

namespace Serenity {
ExternalChargeController::ExternalChargeController(const Settings& settings)
  : _settings(std::make_shared<Settings>(settings)) {
}

const std::vector<std::pair<double, Point>>& ExternalChargeController::getExternalCharges() {
  if (!this->_externalCharges) {
    this->produceExternalCharges();
  }
  return *_externalCharges;
}

void ExternalChargeController::produceExternalCharges() {
  std::vector<std::pair<double, Point>> allExtCharges;
  auto ef = this->_settings->efield;
  if (ef.pos1.size() != 3 || ef.pos2.size() != 3)
    throw SerenityError("Error: The electric field direction/position vector must have three coordinates.");
  if (ef.use && !ef.analytical) {
    EFieldPlates plates(Eigen::Map<Eigen::Vector3d>(&ef.pos1[0]), Eigen::Map<Eigen::Vector3d>(&ef.pos2[0]), ef.distance,
                        ef.nRings, ef.radius, ef.fieldStrength, ef.nameOutput);
    allExtCharges.insert(allExtCharges.end(), plates.getPairList().begin(), plates.getPairList().end());
  }
  if (!this->_settings->extCharges.externalChargesFile.empty()) {
    const auto externalCharges = ExternalChargeReader::readExternalChargeFile(this->_settings->extCharges.externalChargesFile);
    allExtCharges.insert(allExtCharges.end(), externalCharges.begin(), externalCharges.end());
  }
  this->_externalCharges = std::make_unique<std::vector<std::pair<double, Point>>>(allExtCharges);
}
void ExternalChargeController::clear() {
  _externalCharges = nullptr;
}
} // namespace Serenity