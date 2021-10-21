/**
 * @file CC2Sigmavector.h
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

#ifndef LRSCF_CC2SIGMAVECTOR
#define LRSCF_CC2SIGMAVECTOR

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Sigmavectors/RICC2/XWFController.h"

/* Include Std and External Headers */

namespace Serenity {

/**
 * @class CC2Sigmavector CC2Sigmavector.h
 * @brief Implementation for CC2 (left and right) Jacobian transformations and CC2 transition density matrices.\n
 *
 * For reference, please have a look into XWFController.h.
 */
template<Options::SCF_MODES SCFMode>
class CC2Sigmavector : public XWFController<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param lrscf LRSCFController holding all the necessary information.
   */
  CC2Sigmavector(std::shared_ptr<LRSCFController<SCFMode>> lrscf);

  /**
   * @brief Destructor.
   */
  virtual ~CC2Sigmavector() = default;

  /**
   * @brief Return the right Jacobian transformation of the incoming guess vector.\n
   *         Must be overriden by derived class.
   * @param guessVector The guess vector.
   * @param eigenvalue The current eigenvalue (non-linearity of the eigenvalue problem).
   * @return The sigma vector.
   */
  Eigen::VectorXd getRightXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) override final;

  /**
   * @brief Return the left Jacobian transformation of the incoming guess vector.\n
   *        Must be overriden by derived class.
   * @param guessVector The guess vector.
   * @param eigenvalue The current eigenvalue (non-linearity of the eigenvalue problem).
   * @return The sigma vector.
   */
  Eigen::VectorXd getLeftXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) override final;

  /**
   * @brief Calculates the CC2 residual function.
   */
  void calculateResidual() override final;

  /**
   * @brief Calculates the E-intermediate.
   */
  void calculateE() override final;

  /**
   * @brief Calculates the ground-state Lagrange multiplier.
   */
  void calcGroundStateLagrangeMultiplier() override final;

  /**
   * @brief Calculates the transition-moment Lagrange multiplier.
   * @param eigenvectors The eigenvectors to be normalized.
   * @param eigenvectors The associated eigenvalues.
   */
  void calcTransitionMomentLagrangeMultiplier(std::vector<Eigen::MatrixXd>& eigenvectors,
                                              Eigen::VectorXd eigenvalues) override final;

  /**
   * @brief Calculates the one-particle density matrices.
   * @param eigenvectors The eigenvectors to be normalized.
   * @param eigenvectors The associated eigenvalues.
   * @param densityMatrices The densityMatrices to be calcualted.
   */
  void calcDensityMatrices(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues,
                           std::vector<Eigen::MatrixXd>& densityMatrices) override final;

  /**
   * @brief Solves a linear equation system A^T x = b for x, where A is the CC2 Jacobian.
   * @param rightHandSides The right-hand sides.
   * @param eigenvalues The eigenvalue of the Jacobian.
   * @return The solution vectors.
   */
  Eigen::MatrixXd lagrangianSubspaceSolve(Eigen::Ref<Eigen::MatrixXd> rightHandSides, Eigen::Ref<Eigen::VectorXd> eigenvalues);

  /**
   * @brief Returns the virt-virt part of the similarity-transformed Fock matrix.
   * @param trafoVector The singles part of the transformation.
   * @param XQ An intermediate.
   * @param Yia Another intermediate.
   * @param iSpin Spin identifier.
   * @return The virt-virt part of the similarity-transformed Fock matrix.
   */
  Eigen::MatrixXd transformedABFock(Eigen::Ref<Eigen::VectorXd> trafovector, Eigen::Ref<Eigen::VectorXd> XQ,
                                    Eigen::MatrixXd& Yia, int iSpin = 1);

}; // class CC2Sigmavector
} // namespace Serenity

#endif /* LRSCF_CC2SIGMAVECTOR */
