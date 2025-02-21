/**
 * @file LocalAOSigmaVectorWrapper.h
 *
 * @author Moritz Bensberg
 * @date Oct 21, 2019
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

#ifndef POSTHF_LOCALCORRELATION_LOCALAOSIGMAVECTORWRAPPER_H_
#define POSTHF_LOCALCORRELATION_LOCALAOSIGMAVECTORWRAPPER_H_
/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h" //MatrixInBasis definition.
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices.
#include <memory>      //smrt_ptr
#include <vector>      //std::vector.

namespace Serenity {

/* Forward Declarations */
class SingleSubstitution;
class PAOController;
class SystemController;

/**
 * @class LocalAOSigmaVectorWrapper LocalAOSigmaVectorWrapper.h
 * @brief Static class for the sigma vector that is used in local-coupled cluster.
 */
class LocalAOSigmaVectorWrapper {
 private:
  /*
   * Purely static. Constructed is never called.
   */
  LocalAOSigmaVectorWrapper() = default;
  virtual ~LocalAOSigmaVectorWrapper() = default;
  /**
   * @brief Calculated the perturbed density matrix.
   * @param system The system controller.
   * @param singles The singles.
   * @param paoController The PAO controller.
   * @param testRun Flag for a test run. This is only true for some specific tests.
   * @return The perturbed density matrix.
   */
  static MatrixInBasis<RESTRICTED> calculatePerturbedDensityMatrix(std::shared_ptr<SystemController> system,
                                                                   std::vector<std::shared_ptr<SingleSubstitution>> singles,
                                                                   std::shared_ptr<PAOController> paoController,
                                                                   bool testRun = false);

 public:
  /**
   * @brief Getter for the sigma vector in AO basis.
   * @param system The system controller.
   * @param singles The singles, which are needed for their amplitudes.
   * @param paoController The PAO controller.
   * @param testRun Flag for a test run. This is only true for some specific tests.
   * @return The sigma vector in AO basis.
   */
  static MatrixInBasis<RESTRICTED> getSigmaVector_AO_AO(std::shared_ptr<SystemController> system,
                                                        std::vector<std::shared_ptr<SingleSubstitution>> singles,
                                                        std::shared_ptr<PAOController> paoController, bool testRun = false);
};

} /* namespace Serenity */

#endif /* POSTHF_LOCALCORRELATION_LOCALAOSIGMAVECTORWRAPPER_H_ */
