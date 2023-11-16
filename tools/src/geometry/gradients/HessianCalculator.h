/**
 * @file HessianCalculator.h
 *
 * @date Feb 23, 2017
 * @author Kevin Klahr
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
#ifndef GEOMETRY_GRADIENTS_HESSIANCALCULATOR_H_
#define GEOMETRY_GRADIENTS_HESSIANCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "settings/DFTOptions.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;
struct Settings;
struct EmbeddingSettings;
/**
 * @class
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
   * @param system The system holding the geometry of interest.
   * @return Eigen::MatrixXd The Hessian.
   */
  virtual Eigen::MatrixXd calcHessian(std::shared_ptr<SystemController> systemController) = 0;
  /**
   * @brief Calculates a subsystem Hessian.
   *
   * TODO: this function should be rebuild in order to allow arbitrary embedding schemes.
   *
   * @param activeSystems       The list of active systems
   * @param passiveSystems      The list of passive systems (not part of the hessian).
   * @param FaTnaddKinFunc      The non additive kinetic energy functional.
   * @param FaTnaddXCFunc       The non additive exchange correlation energy functional.
   * @param FatmaxCycles        The maximum number of freeze and thaw cycles to be used.
   * @param FaTenergyConvThresh The freeze and thaw convergence threshold.
   * @param FaTgridCutOff       The grid cut-off for the freeze and thaw procedure.
   * @param dispersion          The type of dispersion correction.
   * @return Eigen::MatrixXd    The Hessian.
   */
  virtual Eigen::MatrixXd calcFaTHessian(std::vector<std::shared_ptr<SystemController>> activeSystems,
                                         std::vector<std::shared_ptr<SystemController>> passiveSystems,
                                         EmbeddingSettings embedding, int FatmaxCycles, double FaTenergyConvThresh,
                                         double FaTgridCutOff) = 0;
  /**
   * @brief Calculates the frequencies and normal modes, dumps them to a file in the systems' directories.
   *
   * @param hessian  The hessian.
   * @param geometry The corresponding geometry.
   * @param settings The settings of all systems involved
   */
  virtual void frequencyCalculation(Eigen::MatrixXd hessian, std::shared_ptr<Geometry> geometry,
                                    const std::vector<Settings> settings) = 0;
};

} /* namespace Serenity */

#endif /* GEOMETRY_GRADIENTS_HESSIANCALCULATOR_H_ */
