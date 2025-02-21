/**
 * @file NonOrthogonalLocalization.cpp
 *
 * @date Dec 7, 2015
 * @author David Schnieders
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

/* Include Class Header*/
#include "analysis/orbitalLocalization/NonOrthogonalLocalization.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"                                    //Loop shells.
#include "data/OrbitalController.h"                         // OrbitalController
#include "data/grid/BasisFunctionOnGridControllerFactory.h" // BasisFunctionOnGridControllerFactory
#include "data/matrices/CoefficientMatrix.h"                // CoefficientMatrix
#include "integrals/wrappers/Libint.h"                      // Libint
#include "math/Matrix.h"                                    // Matrix
#include "math/optimizer/BFGS.h"                            // BFGS
#include "system/SystemController.h"                        // SystemController
#include "tasks/LocalizationTask.h"                         // LocalizationTask

namespace Serenity {

template<Options::SCF_MODES SCFMode>
NonOrthogonalLocalization<SCFMode>::NonOrthogonalLocalization(std::shared_ptr<SystemController> systemController)
  : _system(systemController) {
}

template<Options::SCF_MODES SCFMode>
void NonOrthogonalLocalization<SCFMode>::localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                                          SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) {
  /*
   * Do a Foster-Boys localization first
   */
  LocalizationTask superSystemOrbLocalization(_system);
  superSystemOrbLocalization.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS;
  superSystemOrbLocalization.run();

  CoefficientMatrix<SCFMode> coefficients = orbitals.getCoefficients();
  const auto& nOccOrbs = _system->getNOccupiedOrbitals<SCFMode>();

  auto& basis = _system->getBasisController()->getBasis();
  auto basisController = _system->getBasisController();
  auto gridController = _system->getGridController();
  auto basisFuncOnGridController(
      BasisFunctionOnGridControllerFactory::produce(_system->getSettings(), basisController, gridController));

  /*
   * Get dipole and quadrupole integrals and oder them
   */

  Matrix<double> IntMatrix(basisController->getNBasisFunctions(), basisController->getNBasisFunctions());
  IntMatrix.setZero();
  std::vector<Matrix<double>> dipoles(3, IntMatrix);
  std::vector<Matrix<double>> quadrupoles(3, IntMatrix);

  auto& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::emultipole2, 0, 2);

  for (unsigned int i = 0; i < basis.size(); i++) {
    for (unsigned int j = 0; j < basis.size(); j++) {
      Eigen::MatrixXd quadrupoleInts;
      libint.compute(LIBINT_OPERATOR::emultipole2, 0, *basis[i], *basis[j], quadrupoleInts);
      for (unsigned int k = 0; k < basis[i]->getNContracted(); k++) {
        auto mu = basisController->extendedIndex(i) + k;
        for (unsigned int l = 0; l < basis[j]->getNContracted(); l++) {
          auto nu = basisController->extendedIndex(j) + l;
          const unsigned int nJ = basis[j]->getNContracted();
          dipoles[0](mu, nu) = quadrupoleInts((nJ * k + l), 1);
          dipoles[1](mu, nu) = quadrupoleInts((nJ * k + l), 2);
          dipoles[2](mu, nu) = quadrupoleInts((nJ * k + l), 3);
          quadrupoles[0](mu, nu) = quadrupoleInts((nJ * k + l), 4);
          quadrupoles[1](mu, nu) = quadrupoleInts((nJ * k + l), 7);
          quadrupoles[2](mu, nu) = quadrupoleInts((nJ * k + l), 9);
        }
      }
    }
  };
  libint.finalize(LIBINT_OPERATOR::emultipole2, 0, 2);
  /*
   * AO dipole Integrals -> MO dipole Integrals
   */

  for_spin(coefficients, nOccOrbs, orbitalRange) {
    std::vector<Matrix<double>> moDipoleInt(3, IntMatrix);
    std::vector<Matrix<double>> moQuadrupoleInt(3, IntMatrix);

    moDipoleInt[0] = coefficients_spin.transpose() * dipoles[0] * coefficients_spin;
    moDipoleInt[1] = coefficients_spin.transpose() * dipoles[1] * coefficients_spin;
    moDipoleInt[2] = coefficients_spin.transpose() * dipoles[2] * coefficients_spin;
    moQuadrupoleInt[0] = coefficients_spin.transpose() * quadrupoles[0] * coefficients_spin;
    moQuadrupoleInt[1] = coefficients_spin.transpose() * quadrupoles[1] * coefficients_spin;
    moQuadrupoleInt[2] = coefficients_spin.transpose() * quadrupoles[2] * coefficients_spin;

    /*
     * Build up an approximate transformation matrix as initial guess (see FB localization)
     */

    Matrix<double> transMatrix(nOccOrbs_spin, nOccOrbs_spin);
    transMatrix.setZero();
    for (unsigned int iOrb = 0; iOrb < orbitalRange_spin.size(); ++iOrb) {
      unsigned int i = orbitalRange_spin[iOrb];
      for (unsigned int jOrb = 0; jOrb < orbitalRange_spin.size(); ++jOrb) {
        unsigned int j = orbitalRange_spin[jOrb];
        double factorA = 0.0;
        double factorB = 0.0;
        for (unsigned int m = 0; m < 3; m++) {
          factorA += moDipoleInt[m](i, j) * moDipoleInt[m](i, j) -
                     ((moDipoleInt[m](i, i) - moDipoleInt[m](j, j)) * (moDipoleInt[m](i, i) - moDipoleInt[m](j, j))) / 4;
          factorB += ((moDipoleInt[m](i, i) - moDipoleInt[m](j, j)) * moDipoleInt[m](i, j));
        }
        if (i == j) {
          transMatrix(i, j) = 1;
        }
        else if (i > j) {
          transMatrix(i, j) = 0.01 * (-factorB / factorA) / 4;
        }
        else {
          transMatrix(i, j) = 0.01 * (-factorB / factorA) / 4;
        };
      }
    };

    Matrix<double> normMatrix(nOccOrbs_spin, nOccOrbs_spin);
    normMatrix.setZero();

    for (unsigned int iOrb = 0; iOrb < orbitalRange_spin.size(); ++iOrb) {
      unsigned int i = orbitalRange_spin[iOrb];
      for (unsigned int jOrb = 0; jOrb < orbitalRange_spin.size(); ++jOrb) {
        unsigned int j = orbitalRange_spin[jOrb];
        normMatrix(i, i) += transMatrix(j, i) * transMatrix(j, i);
      }
    };

    for (unsigned int iOrb = 0; iOrb < orbitalRange_spin.size(); ++iOrb) {
      unsigned int i = orbitalRange_spin[iOrb];
      normMatrix(i, i) = 1 / sqrt(normMatrix(i, i));
    };

    transMatrix = normMatrix * transMatrix;

    /*
     * Resize matrices
     */

    for (unsigned int i = 0; i < 3; i++) {
      moDipoleInt[i].conservativeResize(nOccOrbs_spin, nOccOrbs_spin);
      moQuadrupoleInt[i].conservativeResize(nOccOrbs_spin, nOccOrbs_spin);
    }

    /*
     * Optimization of the transformation matrix as described in
     *
     * H. Feng, J. Bian, L. Li, W. Yang: An efficient method for constructing
     * nonorthogonal localized molecular orbitals; J. Chem. Phys. 120, 9458 (2004)
     */

    Matrix<double> identity(nOccOrbs_spin, nOccOrbs_spin);
    identity.setIdentity();
    Eigen::VectorXd omegaGradient(nOccOrbs_spin);
    omegaGradient.setZero();
    Eigen::VectorXd wGradient(nOccOrbs_spin);
    wGradient.setZero();
    Eigen::VectorXd lambdaGradient(4);
    lambdaGradient.setZero();
    Eigen::VectorXd lambda(4);
    lambda.setZero();
    double omegaValue = 0.0;
    double wValue = 0.0;
    double wValueOld = 99999999.9;
    double dipoleScalar = 0.0;
    int cycle = 0;
    unsigned int cycle2 = 0;
    double wValueOldOuter = 99999999.9;

    for (unsigned int i = 0; i < nOccOrbs_spin; i++) {
      wValueOld = 99999999.9;

      cycle2 = 0;

      lambda.setZero();

      BFGS optimizer2(lambda);
      auto const updateFunction2 = [&](const Eigen::VectorXd& parameters2, double& value2, Eigen::VectorXd& gradients2,
                                       std::shared_ptr<Eigen::MatrixXd> hessian2, bool print2) {
        (void)print2;
        (void)hessian2;
        cycle = 0;

        lambda = parameters2;

        Eigen::VectorXd col = transMatrix.col(i);
        BFGS optimizer(col);
        auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                        std::shared_ptr<Eigen::MatrixXd> hessian, bool print) {
          (void)print;
          (void)hessian;

          transMatrix.col(i) = parameters;

          omegaValue = 0.0;
          wValue = 0.0;
          wGradient.setZero();
          omegaGradient.setZero();

          /*
           * Calculate value, gradient and Hessian
           */
          for (unsigned int j = 0; j < 3; j++) {
            dipoleScalar = (transMatrix.col(i)).transpose() * moDipoleInt[j] * transMatrix.col(i);
            omegaValue += (transMatrix.col(i)).transpose() * moQuadrupoleInt[j] * transMatrix.col(i);
            omegaValue -= dipoleScalar * dipoleScalar;
            omegaGradient += 2.0 * moQuadrupoleInt[j] * transMatrix.col(i);
            omegaGradient -= 4.0 * (dipoleScalar) * (moDipoleInt[j] * transMatrix.col(i));
            wValue -= lambda(j, 0) * (dipoleScalar - moDipoleInt[j](i, i));
            wValue += 20.0 * (dipoleScalar - moDipoleInt[j](i, i)) * (dipoleScalar - moDipoleInt[j](i, i));
            wGradient -= 2.0 * (lambda(j, 0) * moDipoleInt[j] * transMatrix.col(i));
            wGradient += 80.0 * (((dipoleScalar)-moDipoleInt[j](i, i)) * moDipoleInt[j] * transMatrix.col(i));
            lambdaGradient(j) = dipoleScalar - moDipoleInt[j](i, i);
          }

          double scalar = (transMatrix.col(i)).transpose() * transMatrix.col(i);
          wValue += omegaValue;
          wValue -= lambda[3] * (scalar - 1.0);
          wValue += 20.0 * (scalar - 1.0) * (scalar - 1.0);
          wGradient += omegaGradient;
          wGradient -= 2.0 * (lambda[3] * transMatrix.col(i));
          wGradient += 80.0 * ((scalar - 1.0) * transMatrix.col(i));
          lambdaGradient[3] = scalar - 1.0;

          gradients = wGradient;
          value = wValue;

          cycle += 1;
          bool converged = false;
          if ((fabs(wValueOld - wValue) < 1e-8) or cycle >= 100) {
            wValueOld = 99999999.9;
            converged = true;
          };

          wValueOld = wValue;

          return converged;
        };

        /*
         * Optimize inner function: Transformation matrix column
         */
        optimizer.optimize(updateFunction);

        omegaValue = 0.0;
        wValue = 0.0;
        wGradient.setZero();
        omegaGradient.setZero();

        /*
         * Calculate value, gradient, Hessian
         */
        for (unsigned int j = 0; j < 3; j++) {
          dipoleScalar = (transMatrix.col(i)).transpose() * moDipoleInt[j] * transMatrix.col(i);
          omegaValue += (transMatrix.col(i)).transpose() * moQuadrupoleInt[j] * transMatrix.col(i);
          omegaValue -= dipoleScalar * dipoleScalar;
          wValue -= lambda[j] * (dipoleScalar - moDipoleInt[j](i, i));
          wValue += 20.0 * (dipoleScalar - moDipoleInt[j](i, i)) * (dipoleScalar - moDipoleInt[j](i, i));
          lambdaGradient[j] = dipoleScalar - moDipoleInt[j](i, i);
        }

        double scalar = (transMatrix.col(i)).transpose() * transMatrix.col(i);
        wValue += omegaValue;
        wValue -= lambda[3] * (scalar - 1.0);
        wValue += 20.0 * (scalar - 1.0) * (scalar - 1.0);
        lambdaGradient[3] = scalar - 1.0;

        gradients2 = lambdaGradient;

        value2 = wValue;
        cycle2 += 1;
        bool converged2 = false;
        if ((fabs(wValueOldOuter - wValue) < 1e-10) or cycle2 >= maxSweeps or cycle < 2) {
          wValueOld = 99999999.9;
          converged2 = true;
        };
        wValueOldOuter = wValue;

        transMatrix.col(i) = col;

        return converged2;
      };

      /*
       * Optimize second function: Constraints
       */
      optimizer2.optimize(updateFunction2);
    };

    /*
     * Apply coefficientmatrix
     */
    Matrix<double> coeffMatrix(nOccOrbs_spin, basisController->getNBasisFunctions());

    /*
     * transform the coefficient matrix
     */
    for (unsigned int i = 0; i < basisController->getNBasisFunctions(); i++) {
      for (unsigned int jOrb = 0; jOrb < orbitalRange_spin.size(); ++jOrb) {
        unsigned int j = orbitalRange_spin[jOrb];
        coeffMatrix(j, i) = coefficients_spin(i, j);
      }
    };

    coeffMatrix = transMatrix.transpose() * coeffMatrix;

    for (unsigned int i = 0; i < basisController->getNBasisFunctions(); i++) {
      for (unsigned int jOrb = 0; jOrb < orbitalRange_spin.size(); ++jOrb) {
        unsigned int j = orbitalRange_spin[jOrb];
        coefficients_spin(i, j) = coeffMatrix(j, i);
      }
    };
  };

  orbitals.updateOrbitals(coefficients, orbitals.getEigenvalues());
}

template class NonOrthogonalLocalization<Options::SCF_MODES::RESTRICTED>;
template class NonOrthogonalLocalization<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
