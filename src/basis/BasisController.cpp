/**
 * @file   BasisController.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 * 
 * @date   30. Juli 2015, 17:19
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
#include "basis/BasisController.h"
/* Include Serenity Internal Headers */
#include "parameters/Constants.h"
#include "integrals/wrappers/Libint.h"
#include "basis/Shell.h"
/* Include Std and External Headers */
#include <cassert>


namespace Serenity {

void BasisController::produceBasis() {
  assert(!_basis);
  _basis = produceBasisFunctionVector();
  assert(_basis);
  for (auto& basFunc : *_basis) {
    basFunc->addSensitiveObject(this->_self);
  }
  /*
   * For a nice handling of the basis some further information is very helpful which is
   * precalculated and stored here. In particular, some fiddling around with indices is necessary,
   * because a BasisFunction is in fact a complete shell of basis functions, i.e. the Basis contains
   * less entries than the actual number of basis functions. The quantities calculated here make
   * it possible to always convert an index in the vector of shells to an actual index for the
   * basis functions (which e.g. determines the size of matrices).
   */
  const unsigned int nBasRed = _basis->size();
  _extendedIndex.resize(nBasRed);
  _extendedIndexCart.resize(nBasRed);
  _extendedIndexSpherical.resize(nBasRed);
  _nBasisFunctions = 0;
  _nBasisFunctionsCart = 0;
  _nBasisFunctionsSpherical = 0;
  _maxAngularMomentum = 0;
  /*
   * Determine _nBasisFunctions and store the relationship
   * between extended->reduced indices.
   */
  for (unsigned int i = 0; i < nBasRed; ++i) {
    const unsigned int angMom = (*_basis)[i]->getAngularMomentum();
    _extendedIndex[i] = _nBasisFunctions;
    _extendedIndexCart[i] = _nBasisFunctionsCart;
    _extendedIndexSpherical[i] = _nBasisFunctionsSpherical;
    _nBasisFunctions += (*_basis)[i]->getNContracted();
    _nBasisFunctionsCart += N_SHELL_CART[angMom];
    _nBasisFunctionsSpherical += N_SHELL_SPH[angMom];
    /*
     * Check for maxAM
     */
    if (angMom > _maxAngularMomentum) {
      _maxAngularMomentum = angMom;
    }
  }
  /*
   * Precalculate the backrelation reduced->extended indices
   */
  _reducedIndex.reserve(_nBasisFunctionsCart);
  for (unsigned int i = 0, iRed = 0; iRed < nBasRed; ++iRed) {
    for (unsigned int j = 0; j < (*_basis)[iRed]->getNContracted(); ++j) {
      _reducedIndex[i + j] = iRed;
    }
    i += (*_basis)[iRed]->getNContracted();
  }
  _pureCartesian = (_nBasisFunctions == _nBasisFunctionsCart);
  _pureSpherical = (_nBasisFunctions == _nBasisFunctionsSpherical);
  postConstruction();
}

void BasisController::calculateRIPrescreeningFactors(){
  _RIPrescreeningFactors = std::make_shared<std::vector<ShellPairData> >();
  // intialize libint
  auto& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::coulomb,0,2);
  Eigen::MatrixXd integrals;
  const unsigned int nShells =  (*_basis).size();
      for (unsigned int i=0;i<nShells;++i){
        const auto& shellI =  *(*_basis)[i];
          // calculate integrals
          if (libint.compute(libint2::Operator::coulomb,0,shellI,shellI,integrals)){
            (*_RIPrescreeningFactors).push_back(ShellPairData(i, i,sqrt(integrals.maxCoeff())));
          } /* if (prescreen) */
      } /* i/shellI */
      // finalize libint
      libint.finalize(libint2::Operator::coulomb,0,2);
      // sort the list
//      std::sort((*_RIPrescreeningFactors).begin(), (*_RIPrescreeningFactors).end());
//      std::reverse((*_RIPrescreeningFactors).begin(),(*_RIPrescreeningFactors).end());

}

void BasisController::createShellPairData(){
  _shellPairList = std::make_shared<std::vector<ShellPairData> >();
  // intialize libint
  auto& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::coulomb,0,4);
  // loops over shells
  Eigen::MatrixXd integrals;
  const unsigned int nShells =  (*_basis).size();
      for (unsigned int i=0;i<nShells;++i){
        const auto& shellI =  *(*_basis)[i];
        for (unsigned int j=0;j<=i;++j){
          const auto& shellJ = *(*_basis)[j];
          // calculate integrals
          if (libint.compute(libint2::Operator::coulomb,0,shellI,shellJ,shellI,shellJ,integrals)){
            (*_shellPairList).push_back(ShellPairData(i, j,sqrt(integrals.maxCoeff())));
          } /* if (prescreen) */
        } /* j/shellJ */
      } /* i/shellI */
  // finalize libint
  libint.finalize(libint2::Operator::coulomb,0,4);
  // sort the list
  std::sort((*_shellPairList).begin(), (*_shellPairList).end());
  std::reverse((*_shellPairList).begin(),(*_shellPairList).end());

  // delete everything that is definitely too small in its contribution
//  for(unsigned int i=0;i<(*_shellPairList).size();++i){
//    if (((*_shellPairList)[0].factor*(*_shellPairList)[i].factor)<1E-10){
//      (*_shellPairList).erase((*_shellPairList).begin()+i,(*_shellPairList).end());
//      break;
//    }
//  }
}

} /* namespace Serenity */
