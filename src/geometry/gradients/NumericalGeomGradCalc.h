/**
 * @file   NumericalGeomGradCalc.h
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
#ifndef GEOMETRY_GRADIENTS_NUMERICALGEOMGRADCALC_H_
#define GEOMETRY_GRADIENTS_NUMERICALGEOMGRADCALC_H_
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "geometry/gradients/GeometryGradientCalculator.h"
#include "settings/DFTOptions.h"

namespace Serenity {
/**
 * @class NumericalGeomGradCalc NumericalGeomGradCalc.h
 * @brief A GeometryGradientCalculator that calculates the gradients numerically.
 */
template<Options::SCF_MODES SCFMode>
class NumericalGeomGradCalc : public GeometryGradientCalculator {
 public:
  /**
   * @brief Constructor.
   * @param stepsize The step size for each step in the gradient calculations (in bohr).
   */
  NumericalGeomGradCalc(double stepsize);
  virtual ~NumericalGeomGradCalc() = default;
  /**
   * @brief Calculates the gradients for the given geometry numerically.
   *
   * This function calculates the gradients of all atoms in all
   * three dimensions of space numerically.
   * Calculated geometry gradients are stored in each atom.
   *
   * @param system   The system holding the geometry of interest.
   */
  void calcGradients(std::shared_ptr<SystemController> systemController) override final;

  void calcFDEGradients(std::vector<std::shared_ptr<SystemController>> activeSystems,
                        std::vector<std::shared_ptr<SystemController>> environmentSystems,
                        CompositeFunctionals::KINFUNCTIONALS naddKinFunc,
                        CompositeFunctionals::XCFUNCTIONALS naddXCFunc, double FDEgridCutOff, int FaTmaxCycles,
                        double FaTenergyConvThresh, Options::DFT_DISPERSION_CORRECTIONS dispersion);

 private:
  const double _delta;
};

} /* namespace Serenity */

#endif /* GEOMETRY_GRADIENTS_NUMERICALGEOMGRADCALC_H_ */
