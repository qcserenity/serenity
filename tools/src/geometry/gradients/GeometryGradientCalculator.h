/**
 * @file   GeometryGradientCalculator.h
 *
 * @date   May 18, 2015
 * @author Jan Unsleber
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
#ifndef GEOMETRY_GRADIENTS_GEOMETRYGRADIENTCALCULATOR_H_
#define GEOMETRY_GRADIENTS_GEOMETRYGRADIENTCALCULATOR_H_
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;
/**
 * @class GeometryGradientCalculator GeometryGradientCalculator.h
 * @brief Interface for all geometry gradient calculators.
 */
class GeometryGradientCalculator {
 public:
  GeometryGradientCalculator() = default;
  virtual ~GeometryGradientCalculator() = default;
  /**
   * @brief Calculates the gradients for the given geometry.
   *
   * This function calculates the gradients of all atoms in all
   * three dimensions of space.
   * Calculated geometry gradients are stored in each atom.
   *
   * @param system   The system holding the geometry of interest.
   */
  virtual void calcGradients(std::shared_ptr<SystemController> systemController) = 0;
};

} /* namespace Serenity */

#endif /* GEOMETRY_GRADIENTS_GEOMETRYGRADIENTCALCULATOR_H_ */
