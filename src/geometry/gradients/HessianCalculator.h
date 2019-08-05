/**
 * @file HessianCalculator.h
 *
 * @date Feb 23, 2017
 * @author Kevin Klahr
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
#ifndef GEOMETRY_GRADIENTS_HESSIANCALCULATOR_H_
#define GEOMETRY_GRADIENTS_HESSIANCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "math/Matrix.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;
/**
 * @class HessianCalculator HessianCalculator.h
 * @brief Interface for all hessian calculators.
 */
class HessianCalculator {
public:
	HessianCalculator() = default;
  virtual ~HessianCalculator() = default;
  /**
   * @brief Calculates the hessian for the given geometry.
   *
   * This function calculates the hessian of the given geometry
   *
   * @param system   The system holding the geometry of interest.
   *
   */
  virtual void calcHessian(std::shared_ptr<SystemController> systemController) = 0;

  virtual void calcFaTHessian( std::vector<std::shared_ptr<SystemController> > activeSystems,
  										  std::vector<std::shared_ptr<SystemController> > passiveSystems,
  										  Options::KINFUNCTIONALS FaTnaddKinFunc,
  										  Options::XCFUNCTIONALS FaTnaddXCFunc,
  								      int FatmaxCycles,
  								      double FaTenergyConvThresh,
  										  double FaTgridCutOff,
  								      Options::DFT_DISPERSION_CORRECTIONS dispersion) = 0;



  virtual void frequencyCalculation(Matrix<double> hessian,
                                    std::shared_ptr<Geometry> geometry,
                                    const std::vector<Settings> settings) = 0;
};

} /* namespace Serenity */



#endif /* GEOMETRY_GRADIENTS_HESSIANCALCULATOR_H_ */
