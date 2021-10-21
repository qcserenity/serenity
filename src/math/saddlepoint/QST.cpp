/**
 * @file QST.cpp
 *
 * @date Apr 24, 2017
 * @author Marabel Riesmeier, cleaning and merge Jan Unsleber
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
#include "math/saddlepoint/QST.h"
/* Include Serenity Internal Headers */
#include "integrals/wrappers/Libint.h"
#include "math/optimizer/BFGS.h"
#include "misc/SerenityError.h"

namespace Serenity {

QST::QST(const Eigen::VectorXd& min1, const Eigen::VectorXd& min2, bool preopt, std::unique_ptr<Eigen::VectorXd> guess)
  : _min1(min1),
    _min2(min2),
    _preopt(preopt),
    _guess(std::move(guess)),
    _ts(min2.rows(), min2.cols()),
    _gradient(min2.rows(), min2.cols()),
    _tangent(min2.rows(), min2.cols()),
    _energy(0.0){};

void QST::optimize(std::function<void(const Eigen::VectorXd&, double&, Eigen::VectorXd&, bool)> updateFunction) {
  /*==========================================
   *  Initializations
   *==========================================*/

  // sanity checks
  assert(_min2.size() == _min1.size());
  if (_guess)
    assert(_guess->size() == _min1.size());
  _tangent.setZero();

  // Initial saddle point guess
  unsigned int preoptbreak = 0;
  if (_guess) {
    _ts = *_guess;
    preoptbreak = 1;
    _guess.reset(nullptr);
  }
  else {
    _ts = (_min2 - _min1) * 0.5 + _min1;
    preoptbreak = 2;
  }
  const double distanceMinima = (_min2 - _min1).norm();

  // minimization variables
  double hesmax = 5.0;
  double tStep = 0.0;
  double gradt = 0.0;

  // variables tracking previous iteration
  Eigen::VectorXd oldX;
  Eigen::VectorXd oldGradient;
  Eigen::VectorXd oldTangent;
  double oldEnergy;

  // initial energy and gradients
  updateFunction(_ts, _energy, _gradient, true);

  /*==========================================
   *  Maximization - Minimization Supercycle
   *==========================================*/
  unsigned int nSuperCycles = 0;
  for (unsigned int nTotCycles = 0; nTotCycles < 100; nTotCycles++) {
    // calculate t and throw error if t is not between minima
    double t = (_ts - _min1).dot(_min2 - _min1) / (distanceMinima * distanceMinima);
    if (t > 1.0 or t < 0.0) {
      auto& libint = Libint::getInstance();
      libint.clearAllEngines();
      throw SerenityError((std::string) "The TS guess is not on a QST path between the minima (t=" + std::to_string(t) + ")!");
    }

    // calculate c (parameter for quadratic path)
    const Eigen::VectorXd c((_ts - (1 - t) * _min1 - t * _min2) / (t * (t - 1)));

    /*===========================
     *  Maximization Along Path
     *===========================*/
    for (unsigned int nMaxCycles = 0; nMaxCycles < 15; nMaxCycles++, nTotCycles++) {
      // calculate the tangent to the quadratic path at a given t
      oldTangent = _tangent;
      _tangent = (_min2 - _min1 - c) + 2 * c * t;
      _tangent *= 1.0 / _tangent.norm();

      // calculate gradient of the tangent, set current gradient to old gradient for later use in calculation of hesmax
      const double oldGradt = gradt;
      gradt = _gradient.cwiseProduct(_tangent).sum() / _tangent.norm();

      // Carry over hesmax but double it just in case.
      if (nMaxCycles == 0) {
        hesmax *= 2;
      }
      else {
        hesmax = fabs(gradt - oldGradt) / fabs(tStep);
      }

      // calculate step along t using gradt and hesmax, calculate new t
      tStep = gradt / (hesmax * distanceMinima);
      t += tStep;
      oldX = _ts;
      _ts = (1 - t) * _min1 + t * _min2 + c * t * (t - 1);
      oldGradient = _gradient;
      oldEnergy = _energy;
      updateFunction(_ts, _energy, _gradient, true);
      std::cout << "Maximization Cycle " << nMaxCycles + 1 << "(" << nTotCycles + 1 << ")" << std::endl;

      // backtracking
      if (_energy < oldEnergy and tStep <= 0.1 and nMaxCycles > 0) {
        updateFunction(_ts, _energy, _gradient, true);
        t -= tStep;
        _tangent = (_min2 - _min1 - c) + 2 * c * t;
        _tangent *= 1.0 / _tangent.norm();
        std::cout << "    Backtracking: maximization stopped after " << nMaxCycles + 1 << " cycles." << std::endl;
        break;
      }

      // Convergence check for maximization
      const double convCheck1(_gradient.cwiseProduct(_tangent).sum());
      const double convCheck2(oldGradient.cwiseProduct(oldTangent).sum());
      if (fabs(convCheck1 / convCheck2) <= 0.5) {
        std::cout << "    Maximization converged after " << nMaxCycles + 1 << " cycles." << std::endl;
        nTotCycles++;
        break;
      };
    }

    if (nSuperCycles == preoptbreak and _preopt) {
      std::cout << "    Finished LST/QST pre-optimization." << std::endl;
      break;
    }

    // Orthogonalize gradient against tangent
    Eigen::VectorXd orthGrad(_gradient - _gradient.dot(_tangent) / _tangent.dot(_tangent) * _tangent);

    /*======================
     *  Convergence Checks
     *======================*/
    unsigned int convCriteriaMet = 0;
    double deltaE = fabs(_energy - oldEnergy);
    double rmsGrad = (orthGrad).norm() / sqrt(_gradient.size());
    double maxGrad = (orthGrad).array().abs().maxCoeff();
    std::cout << "Transition State Convergence" << std::endl;
    printf("%4s Energy Change  %15.10f\n", "", deltaE);
    printf("%4s RMS O-Gradient %15.10f\n", "", rmsGrad);
    printf("%4s Max O-Gradient %15.10f\n", "", maxGrad);
    if (deltaE < 1e-06)
      convCriteriaMet++;
    if (rmsGrad < 1e-02)
      convCriteriaMet++;
    if (maxGrad < 2e-02)
      convCriteriaMet++;
    bool converged = (convCriteriaMet > 2);

    /*===========================
     *  Orthogonal Minimization
     *===========================*/
    BFGS minimization(_ts);

    oldX = _ts;
    unsigned int nMinCycles = 0;
    auto const updateFunction2 = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                     std::shared_ptr<Eigen::MatrixXd> hessian, bool print) {
      (void)hessian;
      oldEnergy = value;
      oldGradient = gradients;
      updateFunction(parameters, value, gradients, print);
      orthGrad = gradients.dot(orthGrad) / orthGrad.dot(orthGrad) * orthGrad;
      gradients = 10.0 * orthGrad;

      if (print)
        std::cout << "Minimization Cycle " << nMinCycles + 1 << "(" << nTotCycles + 1 << ")" << std::endl;
      bool minConverged = false;
      // convergence tests
      if (nMinCycles > 0 and orthGrad.norm() < 1e-05) {
        if (print)
          std::cout << "minimization stopped bc gradnorm" << std::endl;
        nTotCycles++;
        minConverged = true;
        _energy = value;
        _gradient = gradients;
      }
      oldX = parameters;
      nMinCycles++;
      return minConverged;
    };

    minimization.optimize(updateFunction2);

    deltaE = fabs(_energy - oldEnergy);
    rmsGrad = (orthGrad).norm() / sqrt(_gradient.size());
    maxGrad = (orthGrad).array().abs().maxCoeff();
    double rmsStep = (_ts - oldX).norm() / sqrt(_ts.size());
    double maxStep = (_ts - oldX).array().abs().maxCoeff();

    std::cout << "Transition State Convergence" << std::endl;
    printf("%4s Energy Change  %15.10f\n", "", deltaE);
    printf("%4s RMS O-Gradient %15.10f\n", "", rmsGrad);
    printf("%4s Max O-Gradient %15.10f\n", "", maxGrad);
    printf("%4s RMS Step       %15.10f\n", "", rmsStep);
    printf("%4s Max Step       %15.10f\n\n", "", maxStep);

    if (converged) {
      std::cout << "Convergence reached. Exiting..." << std::endl;
      break;
    }
    nSuperCycles++;
  }
}

} /* namespace Serenity */
