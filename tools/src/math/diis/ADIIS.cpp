/**
 * @file ADIIS.cpp
 *
 * @date Apr 24, 2017
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

/* Include Class Header*/
#include "math/diis/ADIIS.h"
/* Include Serenity Internal Headers */
#include "misc/WarningTracker.h"
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

ADIIS::ADIIS() : _n(0) {
}

Eigen::MatrixXd ADIIS::optimize(const Eigen::MatrixXd& F, const Eigen::MatrixXd& P) {
  // ADIIS storage and sanity checks
  assert(P.cols() == F.cols());
  assert(P.rows() == F.rows());
  _dmat.push_back(P);
  _fmat.push_back(F);
  if (_dmat.size() > _maxn) {
    _dmat.erase(_dmat.begin());
    _fmat.erase(_fmat.begin());
  }
  assert(_dmat.size() == _fmat.size());
  _n = _dmat.size();

  if (_n == 1) {
    return F;
  }

  // initial coefficients
  Eigen::VectorXd coeff(_n);
  for (unsigned int i = 0; i < _n; ++i) {
    coeff[i] = 1.0 / (double)_n;
  }
  Eigen::VectorXd t(coeff);
  double oldE = 0.0;

  // pre-calculate traces
  Eigen::MatrixXd traceij = Eigen::MatrixXd::Zero(_n, _n);
  Eigen::VectorXd tracei = Eigen::VectorXd::Zero(_n);
  for (unsigned int i = 0; i < _n - 1; ++i) {
    tracei[i] = ((_dmat[i] - _dmat[_n - 1]).transpose() * _fmat[_n - 1]).trace();
    for (unsigned int j = 0; j < _n - 1; ++j) {
      traceij(i, j) = ((_dmat[i] - _dmat[_n - 1]).transpose() * (_fmat[j] - _fmat[_n - 1])).trace();
    }
  }

  /*
   * Get the LBFGS optimizer
   */
  std::shared_ptr<Optimizer> optimizer = std::make_shared<LBFGS>(t);

  unsigned int lbfgscounter = 0;
  auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                  std::shared_ptr<Eigen::MatrixXd> hessian, bool printInfo) {
    (void)hessian;
    (void)printInfo;
    lbfgscounter++;
    // calculate ADIIS coefficent gradient
    //  gradients are calculated for substituted t
    //  instead of the actual coefficients c
    //  substitution: c_i = t_i*t_i / sum(t_i*t_i)
    gradients.setZero();
    // derivative tr[ (D_i-D_n)^T*F ] part
    const double tnorm = parameters.squaredNorm();
    for (unsigned int i = 0; i < _n - 1; ++i) {
      const double pref = 4.0 * parameters[i] * (tnorm - parameters[i] * parameters[i]) / (tnorm * tnorm);
      gradients[i] += pref * tracei[i];
      for (unsigned int j = 0; j < _n; ++j) {
        if (i == j)
          continue;
        gradients[j] -= 4.0 * parameters[i] * parameters[i] * parameters[j] / (tnorm * tnorm) * tracei[i];
      }
    }
    // derivative tr[ (D_i-D_n)^T*(F_j-F_n) ] part
    for (unsigned int i = 0; i < _n - 1; ++i) {
      for (unsigned int j = 0; j < _n - 1; ++j) {
        const double tnorm3 = tnorm * tnorm * tnorm / traceij(i, j);
        for (unsigned int k = 0; k < _n; ++k) {
          if (i == k)
            continue;
          if (j == k)
            continue;
          gradients[k] -= 4.0 * parameters[i] * parameters[i] * parameters[j] * parameters[j] * parameters[k] / tnorm3;
        }
        if (i == j) {
          gradients[i] += 4.0 * parameters[i] * parameters[i] * parameters[i] * (tnorm - parameters[i] * parameters[i]) / tnorm3;
        }
        else {
          const double prefi = 2.0 * (tnorm - 2.0 * parameters[i] * parameters[i]) / tnorm3;
          const double prefj = 2.0 * (tnorm - 2.0 * parameters[j] * parameters[j]) / tnorm3;
          gradients[i] += parameters[i] * parameters[j] * parameters[j] * prefi;
          gradients[j] += parameters[j] * parameters[i] * parameters[i] * prefj;
        }
      }
    }

    coeff = t.cwiseProduct(t) / t.squaredNorm();

    // convergence check
    bool converged = false;
    value = 2.0 * tracei.dot(coeff) + coeff.transpose() * traceij * coeff;
    if (gradients.norm() < 1e-8 or lbfgscounter > 1000) {
      converged = true;
      if (lbfgscounter > 1000)
        WarningTracker::printWarning((std::string) "WARNING: LBFGS reached 1000 cycles.", true);
    }
    oldE = value;

    return converged;
  };

  /* Optimize. */
  optimizer->optimize(updateFunction);
  // compute final F
  Eigen::MatrixXd newF = Eigen::MatrixXd::Zero(_fmat[0].rows(), _fmat[0].cols());
  for (unsigned int i = 0; i < _n; ++i) {
    newF += coeff[i] * _fmat[i];
  }
  return newF;
}

} /* namespace Serenity */
