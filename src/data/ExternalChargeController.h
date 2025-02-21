/**
 * @file ExternalChargeController.h
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

#ifndef SERENITY_EXTERNALCHARGECONTROLLER_H
#define SERENITY_EXTERNALCHARGECONTROLLER_H

/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
class Point;
struct Settings;

/**
 * @class ExternalChargeController ExternalChargeController.h
 * @brief This controller manages all external charges, e.g., point charges in QM/MM or charges modeling
 *        electric fields.
 */
class ExternalChargeController {
 public:
  /**
   * @brief Constructor.
   * @param settings The system settings object which contains all paths/charge settings.
   */
  explicit ExternalChargeController(const Settings& settings);
  /**
   * @brief Getter for all external charges.
   * @return The external charges.
   */
  const std::vector<std::pair<double, Point>>& getExternalCharges();
  /**
   * @brief Clear any cached external charges to free memory.
   */
  void clear();

 private:
  const std::shared_ptr<Settings> _settings;
  std::unique_ptr<std::vector<std::pair<double, Point>>> _externalCharges;
  void produceExternalCharges();
};

} // namespace Serenity

#endif // SERENITY_EXTERNALCHARGECONTROLLER_H
