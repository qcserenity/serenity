/**
 * @file   TwoElecFourCenterIntDecomposer.h
 *
 * @date   Jun 6, 2019
 * @author Lars Hellmann
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

#ifndef INTEGRALS_TWOELECFOURCENTERINTDECOMPOSER_H_
#define INTEGRALS_TWOELECFOURCENTERINTDECOMPOSER_H_

/* Include Serenity Internal Headers */
#include "integrals/CDIntegralController.h"
#include "integrals/wrappers/Libint.h"

namespace Serenity {

/* Forward declarations */
class CholeskyDecomposer;
class BasisController;
class CDStorageController;

/**
 * @class TwoElecFourCenterIntDecomposer TwoElecFourCenterIntDecomposer.h
 *
 *  @brief A class that wraps around the Cholesky Decomposer to calculate the Cholesky Vectors
 *  corresponding to the two-electron four-center repulsion integral matrix (defined by the
 *  system corresponding to the CDIntegralController and the BasisController)
 */
class TwoElecFourCenterIntDecomposer {
 public:
  /**
   * @brief Constructor
   * @param settings The system settings
   * @param basController The basisController holding the basis the matrix is defined in
   * @param cdIntController Cholesky Integral Controller for the system
   * @param label The label used to identify the Cholesky vectors
   */
  TwoElecFourCenterIntDecomposer(const Settings& settings, std::shared_ptr<BasisController> basController,
                                 std::shared_ptr<CDIntegralController> cdIntController, std::string label,
                                 LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb, double mu = 0.0);

  /**
   * @brief Run the decomposition
   */
  void run();

  /**
   * @brief Setter function for a non-default decomposition threshold
   * @param cdThresh The new decomposition threshold
   */
  void setThreshold(double cdThresh);

  /**
   * @brief Get the Cholesky basis. Performs the Cholesky decomposition if needed.
   * @return A vector of the indices contributin to the Cholesky Vectors.
   */
  std::vector<unsigned int> getCholeskyBasis();

 private:
  // The settings of the system
  const Settings& _settings;
  // The basis the system is described in
  std::shared_ptr<BasisController> _basisController;
  // The Controller handling all integral decomposition for a system
  std::shared_ptr<CDIntegralController> _cdIntController;
  // Label for the present matrix
  std::string _label;
  // Cholesky Decomposer performing the actual decomposition
  std::shared_ptr<CholeskyDecomposer> _decomposer;
  // Cholesky decomposition threshold
  double _cdThresh;
  // Operator to calculate the integrals for.
  LIBINT_OPERATOR _op;
  // Scaling for the long-range exchange contribution
  double _mu;
};

} /* namespace Serenity */

#endif /* INTEGRALS_TWOELECFOURCENTERINTDECOMPOSER_H_ */
