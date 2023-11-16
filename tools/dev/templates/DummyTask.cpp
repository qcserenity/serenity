/**
 * @file   DummyTask.cpp
 *
 * @date   Aug 12, 2014
 * @author Thomas Dresselhaus
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
#include "tasks/DummyTask.h"
/* Include Serenity Internal Headers */
#include "math/RegularRankFourTensor.h"
#include "system/SystemController.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
/* Include Std and External Headers */
#include <Eigen/Eigen>

namespace Serenity {

DummyTask::DummyTask( const std::vector<std::shared_ptr<SystemController> > systems) :
    _systemController(systems) {
}

void DummyTask::run() {
  //Example code fragment to calculate and print two-electron four-center integrals. May be replaced at will.

  //Get the basisController for one system.
  const auto& basisController = (_systemController[0])->getBasisController();
  //Obtain the number of basis functions in that system.
  const unsigned int nBasisFunc = basisController->getNBasisFunctions();
  //Initialize a symmetric rank four tensor to store the integrals.
  RegularRankFourTensor<double> eriTensor(nBasisFunc, 0.0);
  //Initialize a pre-build looper to calculate the integrals in an efficient way.
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb,0,basisController,1E-10);

  //Lambda-function used by the looper to distribute the calculated integrals.
  auto const storeERIS = [&eriTensor]
                          (const unsigned int&  a,
                              const unsigned int&  b,
                              const unsigned int&  i,
                              const unsigned int&  j,
                              const Eigen::VectorXd&  integral,
                              const unsigned int threadId) {
    //avoid warnings from unused variables
    (void)threadId;

    eriTensor(b,a,i,j) = integral(0);
    eriTensor(b,a,j,i) = integral(0);
    eriTensor(a,b,j,i) = integral(0);
    eriTensor(a,b,i,j) = integral(0);
    eriTensor(i,j,b,a) = integral(0);
    eriTensor(i,j,a,b) = integral(0);
    eriTensor(j,i,b,a) = integral(0);
    eriTensor(j,i,a,b) = integral(0);
  };

  //Calculate the integrals.
  looper.loop(storeERIS);

  //Initialize a matrix using the Eigen3 library
  Eigen::MatrixXd eriMatrix(nBasisFunc*nBasisFunc,nBasisFunc*nBasisFunc);
  eriMatrix.setZero();

  //Copy the integrals from the Tensor to the Matrix
  for (unsigned int a = 0, ab = 0; a < nBasisFunc; a++){
    for (unsigned int b = 0; b < nBasisFunc; b++, ab++){
      for (unsigned int i = 0, ij = 0; i < nBasisFunc; i++){
        for (unsigned int j = 0; j < nBasisFunc; j++, ij++){
          eriMatrix(ab,ij) = eriTensor(a,b,i,j);
        }
      }
    }
  }

  //Print the matrix if it is not to big.
  if (nBasisFunc <= 4){
    printSmallCaption("Two-Electron Four-Center Integral Matrix");
    std::cout << eriMatrix << std::endl << std::endl;
  }

}

} /* namespace QCpack */
