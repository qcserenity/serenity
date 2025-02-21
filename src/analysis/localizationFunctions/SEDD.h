/**
 * @file SEDD.h
 *
 * @date Apr 4, 2016, reworked on Jul 12, 2017
 * @author Philipp Lenz, reworked by Jan Unsleber
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
#ifndef POSTSCF_LOCALIZATIONFUNCTIONS_SEDD_H_
#define POSTSCF_LOCALIZATIONFUNCTIONS_SEDD_H_
/* Include Serenity Internal Headers */
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declaration */
class GridController;
class SystemController;
namespace Options {
enum class SCF_MODES;
}
/**
 * @class SEDD SEDD.h
 * @brief Prints the calculated SEDD and/or DORI on a cubic grid.
 *
 * The Single Exponential Decay Detector is defined as:
 * \f$ \mathrm{SEDD}(r) = \mathrm{ln} \left( \frac{\nabla \left( \frac{\nabla \rho(r)}{\rho(r)} \right)^2}
 *                              {\rho(r)} \right)^2\f$
 *
 * Ref.: P. de Silva, J. Korchowiec, T. A. Wesolowski, ChemPhysChem 2012, 13, 3462
 *
 *
 * The Density Overlap Regions Indicator is defined as:
 * \f$ \mathrm{DORI}(r) = \frac{\Theta(r)}{1+\Theta(r)}\f$
 * with
 *
 * \f$ \Theta(r) = \frac{\left(\nabla \left( \frac{\nabla \rho(r)}{\rho(r)} \right)^2 \right)^2}
 *                              {\left( \frac{\nabla \rho(r)}{\rho(r)} \right)^6}\f$
 *
 * Ref.: P. de Silva, C. Corminboeuf,  J. Chem. Theory Comput. 2014, 10, 3745-3756
 */
template<Options::SCF_MODES SCFMode>
class SEDD {
 public:
  /**
   * @brief Default Constructor.
   */
  SEDD() = default;
  /**
   * @brief Default destructor.
   */
  virtual ~SEDD() = default;
  /**
   * @brief Returns a labmda function that evaluates and returns the SEDD on a given grid.
   * @param systemController The system controller.
   */
  std::function<Eigen::VectorXd(std::shared_ptr<GridController>)>
  getSEDDLambda(const std::shared_ptr<SystemController> systemController);
  /**
   * @brief Returns a labmda function that evaluates and returns the DORI on a given grid.
   * @param systemController The system controller.
   */
  std::function<Eigen::VectorXd(std::shared_ptr<GridController>)>
  getDORILambda(const std::shared_ptr<SystemController> systemController);
  /**
   * @brief Returns a labmda function that evaluates and returns the signed density on a given grid.
   * @param systemController The system controller.
   *
   * The signed density is the density multiplied by the sign of the density Laplacian.
   */
  std::function<Eigen::VectorXd(std::shared_ptr<GridController>)>
  getSignedDensityLambda(const std::shared_ptr<SystemController> systemController);
};
} /* namespace Serenity */
#endif /* POSTSCF_LOCALIZATIONFUNCTIONS_SEDD_H_ */
