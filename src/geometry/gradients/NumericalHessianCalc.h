/**
 * @file NumericalHessianCalc.h
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
#ifndef GEOMETRY_GRADIENTS_NUMERICALHESSIANCALC_H_
#define GEOMETRY_GRADIENTS_NUMERICALHESSIANCALC_H_

/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "geometry/gradients/HessianCalculator.h"
#include "math/Matrix.h"
#include "settings/DFTOptions.h"
#include "settings/Options.h"

namespace Serenity {

struct Settings;
struct EmbeddingSettings;

/**
 * @class
 * @brief A calculator for fully numerical and semi-numerical Hessians.
 */
template<Options::SCF_MODES SCFMode>
class NumericalHessianCalc : public HessianCalculator {
 public:
  /**
   * @brief Construct a new Numerical Hessian Calc object
   *
   * All gradients are calculated numerically or analytically depending on the settings given here.
   * If stepsizeGrad is greater than 0.0 numerical gradients are used.
   *
   * @param stepsizeGrad The step size for each step in the gradient calculations (in bohr).
   * @param stepsizeHess The step size for each step in the hessian calculations (in bohr).
   * @param printToFile Print the Hessian to a file.
   */
  NumericalHessianCalc(double stepsizeGrad, double stepsizeHess, bool printToFile);
  virtual ~NumericalHessianCalc() = default;
  /**
   * @brief Calculates the hessian for the given geometry.
   *
   * This function calculates the hessian of the given geometry
   *
   * @param system The system holding the geometry of interest.
   * @return Eigen::MatrixXd The Hessian.
   */
  virtual Eigen::MatrixXd calcHessian(std::shared_ptr<SystemController> systemController) override;
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
                                         double FaTgridCutOff) override;
  /**
   * @brief Calculates the frequencies and normal modes, dumps them to a file in the systems' directories.
   *
   * @param hessian  The hessian.
   * @param geometry The corresponding geometry.
   * @param settings The settings of all systems involved
   */
  void frequencyCalculation(Eigen::MatrixXd hessian, std::shared_ptr<Geometry> geometry,
                            const std::vector<Settings> settings) override;

 private:
  const double _deltaGrad;
  const double _deltaHess;
  const bool _printToFile;
};

} /* namespace Serenity */

#endif /* GEOMETRY_GRADIENTS_NUMERICALHESSIANCALC_H_ */
