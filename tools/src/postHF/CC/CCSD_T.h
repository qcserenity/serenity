/**
 * @file CCSD_T.h
 *
 * @date Apr 14, 2016
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

#ifndef ELECTRONICSTRUCTURECALCULATIONS_POSTHF_CC_CCSD_T_H_
#define ELECTRONICSTRUCTURECALCULATIONS_POSTHF_CC_CCSD_T_H_

/* Include Serenity Internal Headers */
#include "postHF/CC/CCSD.h"

namespace Serenity {
/* Forward declarations */
class SystemController;
template<class T>
class Matrix;
/**
 * @class CCSD_T CCSD_T.h
 * @brief A class that can calculate the CCSD energy with an additional function for the triples correction.
 */
class CCSD_T : public CCSD {
 public:
  /**
   * @brief Constructor.
   * @param systemController    The system of interest.
   */
  CCSD_T(std::shared_ptr<SystemController> systemController, double normThreshold, unsigned int maxCycles)
    : CCSD(systemController, normThreshold, maxCycles){};
  /**
   * @brief Default Destructor.
   */
  virtual ~CCSD_T() = default;

  /**
   * @brief Calculates the triples correction from converged CCSD amplitudes.
   * @return The triples correction.
   */
  double calculateTripplesCorrection();

 private:
  void p6(Matrix<Matrix<Matrix<double>>>& mat);
  void r6(Matrix<Matrix<Matrix<double>>>& mat);
};

} /* namespace Serenity */

#endif /* ELECTRONICSTRUCTURECALCULATIONS_POSTHF_CC_CCSD_T_H_ */
