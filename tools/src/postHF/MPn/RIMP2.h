/**
 * @file RIMP2.h
 *
 * @date Aug 14, 2017
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

#ifndef POSTHF_MPN_RIMP2_H_
#define POSTHF_MPN_RIMP2_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/matrices/DensityMatrix.h"
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

/* Forward declarations */
class SystemController;

/**
 * @class  RIMP2 RIMP2.h
 * @brief  Calculates the MP2 correction of the electronic energy using RI.
 *
 * Largely based on the RICC2 code.
 */
template<Options::SCF_MODES SCFMode>
class RIMP2 {
 public:
  /**
   * @brief Constructor
   * @param systemController
   */
  RIMP2(std::shared_ptr<SystemController> systemController, const double sss = 1.0, const double oss = 1.0);
  /**
   * @brief Default destructor;
   */
  virtual ~RIMP2() = default;

  /**
   * @brief Calculates the MP2 correction to the electronic energy using RI.
   * @return Returns the MP2 correction to the electronic energy.
   */
  double calculateCorrection();
  /**
   * @brief Calculates the unrelaxed density correction.
   * @return Returns the RI-MP2 correction of the AO density matrix.
   */
  DensityMatrix<SCFMode> calculateDensityCorrection();

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

  const double _sss;
  const double _oss;
};

} /* namespace Serenity */

#endif /* POSTHF_MPN_RIMP2_H_ */
