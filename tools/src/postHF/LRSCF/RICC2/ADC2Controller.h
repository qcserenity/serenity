/**
 * @file ADC2Controller.h
 *
 * @date Apr 20, 2020
 * @author Niklas Niemeyer
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

#ifndef LRSCF_ADC2CONTROLLER
#define LRSCF_ADC2CONTROLLER

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/RICC2/CC2Controller.h"

namespace Serenity {

/**
 * @class ADC2Controller ADC2Controller.h
 * @brief Implementation for ADC(2) Jacobian transformation and transition density matrices.
 *
 * For reference, please have a look into CC2Controller.h.
 */
template<Options::SCF_MODES SCFMode>
class ADC2Controller : public CC2Controller<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param lrscf LRSCFController holding all the necessary information.
   */
  ADC2Controller(std::shared_ptr<LRSCFController<SCFMode>> lrscf);

  /**
   * @brief Destructor.
   */
  virtual ~ADC2Controller() = default;

  /**
   * @brief Return the right Jacobian transformation of the incoming guess vector.\n
   *         Must be overriden by derived class.
   * @param guessVector The guess vector.
   * @param eigenvalue The current eigenvalue (non-linearity of the eigenvalue problem).
   * @return The sigma vector.
   */
  Eigen::VectorXd getRightCC2Sigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) override;

  /**
   * @brief Return the left Jacobian transformation of the incoming guess vector.\n
   *        Must be overriden by derived class.
   * @param guessVector The guess vector.
   * @param eigenvalue The current eigenvalue (non-linearity of the eigenvalue problem).
   * @return The sigma vector.
   */
  Eigen::VectorXd getLeftCC2Sigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) override {
    return this->getRightCC2Sigma(guessVector, eigenvalue);
  }

  /**
   * @brief Calculates the E-intermediate.
   */
  void calculateE() override;

  /**
   * @brief Calculates the one-particle density matrices.
   * @param eigenvectors The eigenvectors to be normalized.
   * @param eigenvalues The associated eigenvalues.
   * @param densityMatrices The densityMatrices to be calcualted.
   */
  void calculateTransitionDensities(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues,
                                    std::vector<Eigen::MatrixXd>& densityMatrices) override;
}; // class ADC2Controller
} // namespace Serenity

#endif /* LRSCF_ADC2CONTROLLER */
