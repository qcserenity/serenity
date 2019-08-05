/**
 * @file RI_J_IntegralController.cpp
 * @author: Kevin Klahr
 *
 * @date 8. March 2016
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
#include "integrals/RI_J_IntegralController.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/wrappers/Libint.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <cassert>


namespace Serenity {
using namespace std;

RI_J_IntegralController::RI_J_IntegralController(
    std::shared_ptr<BasisController> basisControllerA,
    std::shared_ptr<BasisController> auxBasisController,
    std::shared_ptr<BasisController> basisControllerB) :
        _basisControllerA(basisControllerA),
        _basisControllerB(basisControllerB),
        _auxBasisController(auxBasisController),
        _nBasisFunctions(basisControllerA->getNBasisFunctions()),
        _nAuxFunctions(auxBasisController->getNBasisFunctions()),
        _nAuxFunctionsRed(auxBasisController->getReducedNBasisFunctions()),
        _inverseM(0,0),
        _memManager(MemoryManager::getInstance()){
  assert(_basisControllerA);
  assert(_auxBasisController);
  initialize();
  _basisControllerA->addSensitiveObject(this->_self);
  _auxBasisController->addSensitiveObject(this->_self);
  if (_basisControllerB) _basisControllerB->addSensitiveObject(this->_self);
}


const Matrix<double>& RI_J_IntegralController::getInverseM(){

  if(!_2cIntsAvailable){
    calculate2CenterIntegrals();
  }

  return _inverseM;
}

const std::shared_ptr<BasisController> RI_J_IntegralController::getBasisController(){
  return _basisControllerA;
}

const std::shared_ptr<BasisController> RI_J_IntegralController::getBasisControllerB(){
  return _basisControllerB;
}

const std::shared_ptr<BasisController> RI_J_IntegralController::getAuxBasisController(){
  return _auxBasisController;
}

void RI_J_IntegralController::calculate2CenterIntegrals(){
  assert(!_basisControllerB && "Two center integrals are only available in single basis mode!");
	  takeTime("Calc 2-center ints");
	  /*
	  * Calculate Coulomb metric of aux _basis (P|1/r|Q) (-> eq. [1].(3) )
	  */
	  const auto& auxBasis = _auxBasisController->getBasis();
	  /*
	   * With the current setup this matrix is only needed temporarily. It can be quite large for large
	   * systems, because the auxiliary basis sets are typically quite a bit larger than the normal ones.
	   */
	  MatrixInBasis<RESTRICTED> _M(_auxBasisController);
	  _M.setZero();
	  _libint->initialize(libint2::Operator::coulomb,0,2);
#pragma omp parallel for schedule(static)
	  for (unsigned int i=0; i<_nAuxFunctionsRed; ++i) {
	    const unsigned int pStart = _auxBasisController->extendedIndex(i);
	    const unsigned int nI = auxBasis[i]->getNContracted();
	    const auto& shellA = *auxBasis[i];
	    for (unsigned int j=0; j<=i; ++j) {

	      const unsigned int qStart = _auxBasisController->extendedIndex(j);
	      const unsigned int nJ = auxBasis[j]->getNContracted();
	      const auto& shellB = *auxBasis[j];


	      Eigen::MatrixXd ints;

	      if(_libint->compute(libint2::Operator::coulomb,0,shellA,shellB,ints)){

	        Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(),nJ,nI);

	        _M.block(qStart,pStart,nJ,nI) = tmp;

	        _M.block(pStart,qStart,nI,nJ) = tmp.transpose();

	      }
	    }
	  }

	  _libint->finalize(libint2::Operator::coulomb,0,2);

	  timeTaken(3,"Calc 2-center ints");

	  takeTime("Inversion");
	  {
	  /*
	   * Changed to a pseudo-inverse here. The direct inversion is numerically unstable.
	   */
	  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(_M,
	      Eigen::DecompositionOptions::ComputeEigenvectors);
	  auto eigenvals = eigensolver.eigenvalues().eval();
	  auto eigenvectors = eigensolver.eigenvectors().eval();
	  for (unsigned int i=0; i< _nAuxFunctions; ++i) {
	    if (fabs(eigenvals(i)) < 1e-6) {
	      eigenvals(i) = 0.0;
	    } else {
	      eigenvals(i) = 1.0/eigenvals(i);
	    }
	  }
	  _inverseM.resize(_nAuxFunctions,_nAuxFunctions);
	  _inverseM = (eigenvectors * eigenvals.asDiagonal() * eigenvectors.transpose()).eval();
	  }
	  timeTaken(3,"Inversion");
	  _2cIntsAvailable = true;
}

void RI_J_IntegralController::initialize(){
  _inverseM.resize(0,0);
  _2cIntsAvailable = false;
}

} /* namespace Serenity */
