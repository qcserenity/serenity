/**
 * @file CustomBasisController.cpp
 *
 * @date Sep 11, 2015
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
#include "basis/CustomBasisController.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h" //Constructor.

namespace Serenity {

CustomBasisController::CustomBasisController(std::vector<std::shared_ptr<const Shell>>& basisFunctions,
                                             const std::string& basisString)
  : BasisController(basisString) {
  _basis = std::unique_ptr<Basis>(new Basis());
  _basis->insert(_basis->end(), basisFunctions.begin(), basisFunctions.end());

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
}

std::unique_ptr<Basis> CustomBasisController::produceBasisFunctionVector() {
  return std::move(_basis);
}

void CustomBasisController::postConstruction() {
}

} // namespace Serenity
