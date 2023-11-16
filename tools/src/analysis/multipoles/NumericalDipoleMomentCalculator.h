/**
 * @file NumericalDipoleMomentCalculator.h
 *
 * @date Oct 30, 2015
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
#ifndef NUMERICALDIPOLEMOMENTCALCULATOR_H_
#define NUMERICALDIPOLEMOMENTCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <array>
#include <memory>
#include <vector>

namespace Serenity {

namespace Options {
enum class SCF_MODES;
}

class SystemController;
/**
 * @class NumericalDipoleMomentCalculator NumericalDipoleMomentCalculator.h
 *
 * @brief Calculates Dipole Moment numerically on grid
 *
 */
class NumericalDipoleMomentCalculator {
 public:
  /**
   * @brief Constructor
   */
  NumericalDipoleMomentCalculator() = default;
  /**
   * @brief Destructor
   */
  virtual ~NumericalDipoleMomentCalculator() = default;
  /**
   * @brief Calculates the Dipole Moment from data given from SystemController
   *
   * @param system The system of which the dipole moment should be calculated
   * @return dipoleMoment array{x,y,z}
   */
  template<Options::SCF_MODES SCFMode>
  static Eigen::Vector3d calculateDipoleMoment(std::shared_ptr<SystemController> system);
};

} // namespace Serenity

#endif /* NUMERICALDIPOLEMOMENTCALCULATOR_H_ */
