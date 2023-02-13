/**
 * @file LTMP2.h
 *
 * @date Dec 30, 2022
 * @author Niklas Niemeyer
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

#ifndef POSTHF_MPN_LTMP2_H_
#define POSTHF_MPN_LTMP2_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

/* Forward declarations */
class SystemController;

/**
 * @class  LTMP2 LTMP2.h
 * @brief  Calculates the Spin-opposite scaled (SOS) MP2 correction of the
 * electronic energy using RI and a Laplace Transformation (LT) of the energy denominator.
 *
 * Largely based on the RICC2 code.
 */
template<Options::SCF_MODES SCFMode>
class LTMP2 {
 public:
  /**
   * @brief Constructor
   * @param systemController The system.
   * @param oss Opposite spin scaling.
   * @param laplaceConv Threshold for numerical integration of Laplace transformation.
   */
  LTMP2(std::shared_ptr<SystemController> systemController, const double oss = 1.3, const double laplaceConv = 1e-5);
  /**
   * @brief Default destructor;
   */
  virtual ~LTMP2() = default;

  /**
   * @brief Calculates the Spin-opposite scaled (SOS) MP2 correction of the
   * electronic energy using RI and a Laplace Transformation (LT) of the energy denominator.
   * @return Returns the (SOS) MP2 correction to the electronic energy.
   */
  double calculateCorrection();

 private:
  /**
   * @brief Calculates the correction
   * This function is the restricted/unrestricted part of calculateCorrection() that could not
   * be combined via for_spin macro.
   * @return Returns the MP2 correction to the electronic energy.
   */
  double calculateEnergy();

  std::shared_ptr<SystemController> _systemController;
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jia;

  const double _oss;
  const double _ltConv;
};

} /* namespace Serenity */

#endif /* POSTHF_MPN_LTMP2_H_ */
