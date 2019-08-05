/**
 * @file   SteepestDescent.cpp
 *
 * @date   Oct 12, 2015
 * @author Linus Scholz
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
/* Include Class Header*/
#include "math/optimizer/SteepestDescent.h"


namespace Serenity {
using namespace std;

SteepestDescent::SteepestDescent(Eigen::VectorXd& parameters, double stepWidth) :
    Optimizer(parameters),
    _stepWidth(stepWidth) {
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter" // value is ignored -> not needed
void SteepestDescent::optimize(std::function<bool(
    const Eigen::VectorXd&,
    double&,
    Eigen::VectorXd&,
    std::shared_ptr<Eigen::MatrixXd> hessian,
    bool print)>updateFunction) {
  double value=999999.9;
  Eigen::VectorXd gradients(this->_parameters.size());
  gradients.setZero();
  bool converged=false;
  converged=updateFunction(this->_parameters,value,gradients,nullptr,false);

  while(true){
    if(converged) break;
    this->_parameters -= _stepWidth*gradients;
    converged=updateFunction(this->_parameters,value,gradients,nullptr,false);
  }
}
#pragma GCC diagnostic pop // "-Wunused-parameter"

} /* namespace Serenity */
