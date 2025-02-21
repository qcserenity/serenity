/**
 * @file OrthogonalizationTask.cpp
 *
 * @date Aug 06, 2020
 * @author Anja Massolle
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
#include "tasks/OrthogonalizationTask.h"
/* Include Serenity Internal Headers */
#include "analysis/brokenSymmetry/CorrespondingOrbitals.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutputStream.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/SystemAdditionTask.h"
/* Include Std and External Headers */
#include <boost/math/special_functions/binomial.hpp>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
OrthogonalizationTask<SCFMode>::OrthogonalizationTask(std::vector<std::shared_ptr<SystemController>> systemController,
                                                      std::shared_ptr<SystemController> superSystem)
  : _systemController(systemController), _superSystem(superSystem) {
  Settings settingsSuper;
  if (_superSystem == nullptr) {
    settingsSuper = _systemController[0]->getSettings();
    settingsSuper.name = _systemController[0]->getSettings().name + "Ortho";
    settingsSuper.charge = 0;
    settingsSuper.spin = 0;
    _superSystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settingsSuper);
  }

  if (_systemController.size() > 1) {
    SystemAdditionTask<SCFMode> additionTask(_superSystem, _systemController);
    additionTask.settings.addOccupiedOrbitals = true;
    additionTask.run();
  }
  else {
    _superSystem = _systemController[0];
  }
}

template<Options::SCF_MODES SCFMode>
void OrthogonalizationTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _systemController, {_superSystem});
  printSectionTitle("Orthogonalize Orbitals");
  auto coeff = _superSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  auto nOcc = _superSystem->getNOccupiedOrbitals<SCFMode>();
  auto nBas = _superSystem->getBasisController()->getNBasisFunctions();
  auto& libint = Libint::getInstance();
  _sAO = libint.compute1eInts(LIBINT_OPERATOR::overlap, _superSystem->getBasisController());

  switch (settings.orthogonalizationScheme) {
    case Options::ORTHOGONALIZATION_ALGORITHMS::NONE: {
      OutputControl::nOut << "No orthogonalization scheme was selected" << std::endl;
      break;
    }
    /*
     * According to: Löwdin, Per‐Olov, J. Chem. Phys. 18, 1950, 365-375 and
     *               Löwdin, Per-Olov, Advances in quantum chemistry. Vol. 5. Academic Press, 1970. 185-199.
     */
    case Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN: {
      OutputControl::nOut << "Use Löwdin's orthogonalization scheme" << std::endl;
      for_spin(coeff, nOcc) {
        if (nOcc_spin > 0) {
          Eigen::MatrixXd tmp = (coeff_spin).leftCols(nOcc_spin);
          auto sMO = (tmp).transpose() * _sAO * tmp;
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> sInv(sMO);
          auto tmptmp = sInv.operatorSqrt();
          Eigen::MatrixXd Cortho = tmp * tmptmp.completeOrthogonalDecomposition().pseudoInverse();
          (coeff_spin).leftCols(nOcc_spin) = Cortho;
        }
      };
      break;
    } /*Loewdin*/
    /*
     * According to: Pipek, J. Int. J. Quantum Chem. 1985, 27, 527–546.
     */
    case Options::ORTHOGONALIZATION_ALGORITHMS::PIPEK: {
      OutputControl::nOut << "Use Pipek's orthogonalization scheme" << std::endl;
      auto coeffDum = coeff;
      unsigned int counter = 0;
      while (!checkOrtho(coeffDum, nOcc, true)) {
        OutputControl::nOut << "Iteration " << counter << std::endl;
        for_spin(coeffDum, nOcc) {
          if (nOcc_spin > 0) {
            auto T = this->calcT((coeffDum_spin).leftCols(nOcc_spin), nOcc_spin);
            Eigen::MatrixXd coeffNew = (coeffDum_spin).leftCols(nOcc_spin) * T;
            (coeffDum_spin).leftCols(nOcc_spin) = coeffNew;
          }
        };
        coeffDum = this->normalize(coeffDum);
        counter += 1;
        if (counter > settings.maxIterations) {
          break;
        }
      } /*while*/
      OutputControl::nOut << "Converged after " << counter << " iterations" << std::endl;
      for_spin(coeff, coeffDum, nOcc) {
        if (nOcc_spin > 0) {
          (coeff_spin).leftCols(nOcc_spin) = (coeffDum_spin).leftCols(nOcc_spin);
        }
      };
      break;
    } /*Pipek*/
    /*
     * According to: Broer, R. Int. J. Quantum Chem. 1993, 45, 587-590.
     */
    case Options::ORTHOGONALIZATION_ALGORITHMS::BROER: {
      OutputControl::nOut << "Use Broer's orthogonalization scheme" << std::endl;
      if (_systemController.size() == 2) {
        auto nOccAct = _systemController[0]->template getNOccupiedOrbitals<SCFMode>();
        auto nOccEnv = _systemController[1]->template getNOccupiedOrbitals<SCFMode>();
        auto nBasAct = _systemController[0]->getBasisController()->getNBasisFunctions();
        auto nBasEnv = _systemController[1]->getBasisController()->getNBasisFunctions();

        auto actSuperSysCoeff = _systemController[0]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
        auto envSuperSysCoeff = _systemController[1]->template getActiveOrbitalController<SCFMode>()->getCoefficients();

        if (nBasAct != nBas and nBasEnv != nBas) {
          actSuperSysCoeff = coeff;
          envSuperSysCoeff = coeff;

          auto actSysCoeff = _systemController[0]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
          auto envSysCoeff = _systemController[1]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
          for_spin(actSysCoeff, envSysCoeff, actSuperSysCoeff, envSuperSysCoeff, nOccAct, nOccEnv) {
            actSuperSysCoeff_spin.setZero();
            actSuperSysCoeff_spin.block(0, 0, nBasAct, nOccAct_spin) = (actSysCoeff_spin).leftCols(nOccAct_spin);
            envSuperSysCoeff_spin.setZero();
            envSuperSysCoeff_spin.block(nBasAct, 0, nBasEnv, nOccEnv_spin) = (envSysCoeff_spin).leftCols(nOccEnv_spin);
          };
        }
        CorrespondingOrbitals<SCFMode> corresOrb(actSuperSysCoeff, envSuperSysCoeff, nOccAct, nOccEnv);
        auto corrOrbs = corresOrb.getCorrespondingOrbitals();
        auto corrOrbs1 = corrOrbs.first;
        auto corrOrbs2 = corrOrbs.second;

        auto sMOs = corresOrb.getSPrimePrime();

        for_spin(corrOrbs1, corrOrbs2, sMOs, nOccAct, nOccEnv, coeff) {
          Eigen::MatrixXd orthogonalOrbitals = corrOrbs1_spin;
          nOccAct_spin = corrOrbs1_spin.cols();
          nOccEnv_spin = corrOrbs2_spin.cols();
          unsigned int max = (corrOrbs1_spin).cols();
          if ((corrOrbs1_spin).cols() > (corrOrbs2_spin).cols()) {
            max = (corrOrbs2_spin).cols();
          }
          for (unsigned int i = 0; i < max; i++) {
            orthogonalOrbitals.col(i) = ((corrOrbs1_spin).col(i) - (sMOs_spin)(i, i) * (corrOrbs2_spin).col(i));
          }
          (coeff_spin).setZero();
          (coeff_spin).leftCols(nOccAct_spin) = (orthogonalOrbitals).leftCols(nOccAct_spin);
          (coeff_spin).block(0, nOccAct_spin, nBas, nOccEnv_spin) = (corrOrbs2_spin).leftCols(nOccEnv_spin);
          for (unsigned int i = 0; i < (coeff_spin).cols(); i++) {
            double length = ((coeff_spin).col(i)).transpose() * _sAO * (coeff_spin).col(i);
            if (length >= 1e-10) {
              (coeff_spin).col(i) = (coeff_spin).col(i) / sqrt(length);
            }
          }
        };
        coeff = this->normalize(coeff);
      }
      else {
        throw SerenityError((std::string) "The orthogonalization of orbitals using corresponding orbitals is only "
                                          "implemented for two subsystems.");
      }
      break;
    } /*Broer*/
  }

  /*
   * Check
   */
  _ortho = checkOrtho(coeff, nOcc, true);
  if (!_ortho) {
    WarningTracker::printWarning((std::string) "Warning: The deviation of the MO overlap matrix from the unit matrix "
                                               "is larger than 1e-9 so the orbitals are not fully orthonormal.",
                                 true);
  }

  _superSystem->setDiskMode(true);
  _superSystem->getActiveOrbitalController<SCFMode>()->updateOrbitals(
      coeff, _superSystem->getActiveOrbitalController<SCFMode>()->getEigenvalues());
} /*run*/

template<Options::SCF_MODES SCFMode>
bool OrthogonalizationTask<SCFMode>::checkOrtho(CoefficientMatrix<SCFMode> coeff,
                                                SpinPolarizedData<SCFMode, unsigned int> nOcc, bool print) {
  auto error = Eigen::VectorXd(2);
  unsigned int i = 0;
  for_spin(coeff, nOcc) {
    error.setZero();
    if (nOcc_spin > 0) {
      auto sMO = ((coeff_spin).leftCols(nOcc_spin)).transpose() * _sAO * (coeff_spin).leftCols(nOcc_spin);
      Eigen::MatrixXd unity = Eigen::MatrixXd::Identity(sMO.rows(), sMO.cols());
      if (abs((sMO - unity).maxCoeff()) < abs((sMO - unity).minCoeff())) {
        error[i] = abs((sMO - unity).minCoeff());
      }
      else {
        error[i] = abs((sMO - unity).maxCoeff());
      }
      i += 1;
    }
  };
  if (print) {
    OutputControl::nOut << "Largest MO overlap: " << error.maxCoeff() << std::endl;
  }
  return ((error.maxCoeff() < 1e-9) ? true : false);
} /*_checkOrtho*/

template<Options::SCF_MODES SCFMode>
CoefficientMatrix<SCFMode> OrthogonalizationTask<SCFMode>::normalize(CoefficientMatrix<SCFMode> coeff) {
  auto normCoeff = coeff;
  for_spin(normCoeff, coeff) {
    for (unsigned int i = 0; i < (coeff_spin).cols(); i++) {
      double length = ((coeff_spin).col(i)).transpose() * _sAO * (coeff_spin).col(i);
      if (length >= 1e-10) {
        (normCoeff_spin).col(i) = (coeff_spin).col(i) / sqrt(length);
      }
    }
  };
  return normCoeff;
} /*normalize*/

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd OrthogonalizationTask<SCFMode>::calcT(Eigen::MatrixXd coeff, unsigned int nOcc) {
  Eigen::MatrixXd T(nOcc, nOcc);
  T = T.setZero();
  Eigen::MatrixXd unity = T;
  unity = unity.setIdentity();
  Eigen::MatrixXd sMO = coeff.transpose() * _sAO * coeff;
  auto epsilon = this->calcEpsilon(nOcc, coeff);
  T += this->calcGamma(0, nOcc, epsilon) * unity;
  T += this->calcGamma(1, nOcc, epsilon) * sMO;
  return T;
}; /*calcT*/

template<Options::SCF_MODES SCFMode>
double OrthogonalizationTask<SCFMode>::calcGamma(unsigned int l, unsigned int n, double epsilon) {
  double gammaLin = 0;
  double gammaOrt = 0;
  /*
   * interpolate between near-orthogonal case and near linearly-dependent case
   * Eq 20 and 21 in the Pipek paper (Pipek, J. Int. J. Quantum Chem. 1985, 27, 527–546).
   */
  switch (l) {
    case 0:
      gammaOrt = 1.5;
      gammaLin = 1.0 * n / (n - 1.0);
      break;
    case 1:
      gammaOrt = -0.5;
      gammaLin = -1.0 / (n - 1.0);
      break;
  }
  double gamma = epsilon * gammaLin + (1 - epsilon) * gammaOrt;
  return gamma;
}; /*calcGamma*/

template<Options::SCF_MODES SCFMode>
double OrthogonalizationTask<SCFMode>::calcEpsilon(unsigned int nOcc, Eigen::MatrixXd coeff) {
  Eigen::MatrixXd sMO = coeff.transpose() * _sAO * coeff;
  auto binom = boost::math::binomial_coefficient<double>(nOcc, 2);
  double epsilon = 0.0;
  for (unsigned int i = 0; i < nOcc; i++) {
    for (unsigned int j = 0; j < i; j++) {
      epsilon += abs(sMO(i, j));
    }
  }
  epsilon *= 1 / binom;
  return epsilon;
}; /*calcEpsilon*/

template class OrthogonalizationTask<Options::SCF_MODES::RESTRICTED>;
template class OrthogonalizationTask<Options::SCF_MODES::UNRESTRICTED>;

} /*namespace Serenity*/
