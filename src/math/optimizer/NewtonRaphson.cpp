/**
 * @file   NewtonRaphson.cpp
 *
 * @date   Jul 22, 2015
 * @author Jan Unsleber
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
#include "math/optimizer/NewtonRaphson.h"
/* Include Serenity Internal Headers */
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <Eigen/Core>
#include <Eigen/SVD>
#include <cassert>


namespace Serenity {
using namespace std;

NewtonRaphson::NewtonRaphson(Eigen::VectorXd& parameters, double threshold) :
  Optimizer(parameters),
  _threshold(threshold){
}


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter" // value is ignored
void NewtonRaphson::optimize(std::function<bool(
    const Eigen::VectorXd&,
    double&,
    Eigen::VectorXd&,
    std::shared_ptr<Eigen::MatrixXd> hessian,
    bool print)>updateFunction)  {

  double value=999999.9;
  bool converged=true;
  unsigned int nParams=this->_parameters.size();
  Eigen::VectorXd gradients(nParams);
  gradients.setZero();
  auto hessian=std::make_shared<Eigen::MatrixXd>(nParams,nParams);
  hessian->setZero();
  converged=updateFunction(this->_parameters,value,gradients,hessian,false);

  while(true){
    if(converged) break;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(*hessian, Eigen::ComputeFullV | Eigen::ComputeFullU);
    svd.setThreshold(_threshold);
    Eigen::VectorXd steps = svd.solve(gradients);
    this->_parameters -= steps;
    converged=updateFunction(this->_parameters,value,gradients,hessian,false);
  }

}
#pragma GCC diagnostic pop

} /* namespace Serenity */
