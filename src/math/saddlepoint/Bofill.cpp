/**
 * @file Bofill.cpp
 *
 * @date Jul 13, 2017
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
#include "math/saddlepoint/Bofill.h"
/* Include Std and External Headers */
#include <iostream>
namespace Serenity {

Bofill::Bofill(const Eigen::VectorXd& params, double trustRadius, std::unique_ptr<Eigen::VectorXd> searchDirection)
  : _x(params), _ev(std::move(searchDirection)), _trustradius(trustRadius) {
}

void Bofill::optimize(std::function<void(const Eigen::VectorXd&, double&, Eigen::VectorXd&, bool)> updateFunction) {
  /*===========================================
   *  Implemented, as described in:
   *  Phys. Chem. Chem. Phys., 2002, 4, 11–15
   *  The original paper by Bofill is the following:
   *  J. Comput. Chem., 1994, 15, 1
   *===========================================*/

  // One SD step to start with, note that a semi-empirical hessian
  //  as starting guess would eliminate this step and would also improve
  //  the algorithm significantly.
  const unsigned int nparams = _x.size();
  double value;
  Eigen::VectorXd g(nparams);
  updateFunction(_x, value, g, true);
  Eigen::MatrixXd B = Eigen::MatrixXd::Identity(nparams, nparams);
  Eigen::VectorXd xOld(_x);
  Eigen::VectorXd gOld(g);
  if (_ev != nullptr) {
    //????????????????????????????????????????????????????????
    //    _x -= 0.01*(g - 2.0 * g.dot(*_ev)/((*_ev).dot(*_ev)) *(*_ev) ) ;
    //????????????????????????????????????????????????????????
    _x -= 0.01 * (*_ev);
  }
  else {
    _x -= 0.01 * g;
  }
  updateFunction(_x, value, g, true);
  // DFP Iterations
  unsigned int counter = 0;
  while (true) {
    counter++;
    Eigen::VectorXd dx = _x - xOld;
    Eigen::VectorXd dg = g - gOld;

    /*==================
     *  update Hessian
     *==================*/
    const double tmp1 = dx.transpose() * dx;
    Eigen::VectorXd tmp2(dg + B * dx);
    const double tmpdotdx = tmp2.dot(dx);

    // the Bofill weight factor, improved by a sqrt as described
    //   in the reference paper
    double bofillFacor = (tmpdotdx * tmpdotdx) / (tmp2.dot(tmp2) * dx.dot(dx));
    bofillFacor = sqrt(bofillFacor);

    // Powell  symmetric Broyden (PSB) part of the Bofill algorithm
    B -= (1 - bofillFacor) * (tmp2 * dx.transpose() + dx * tmp2.transpose()) * (1.0 / tmp1);
    B += (1 - bofillFacor) * (tmpdotdx / (tmp1 * tmp1)) * dx * dx.transpose();

    // SR1 part of the Bofill algorithm
    B -= bofillFacor * tmp2 * tmp2.transpose() / tmpdotdx;

    // Hessian decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
    Eigen::MatrixXd ew = Eigen::MatrixXd::Zero(nparams, nparams);
    ew.diagonal().array() = 1.0 / es.eigenvalues().array();

    // Generate inverse Hessian
    Eigen::MatrixXd hInvmax = es.eigenvectors() * ew * es.eigenvectors().transpose();

    // update parameters
    xOld = _x;
    Eigen::VectorXd step = hInvmax * g;

    // check trustradius
    double norm = step.norm();
    if (norm > _trustradius)
      step *= _trustradius / norm;

    // take a step
    _x += step;

    // update gradient
    gOld = g;
    double oldVal = value;
    updateFunction(_x, value, g, true);

    // convergence
    if (fabs(oldVal - value) < 1e-8)
      break;
    if (counter > 100) {
      break;
    }
  }
}
} /* namespace Serenity */
