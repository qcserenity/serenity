/**
 * @file ABERIPotential.cpp
 *
 * @date Jun 23, 2018
 * @author Moritz Bensberg
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
#include "potentials/ABFockMatrixConstruction/ABERIPotential.h"
/* Include Serenity Internal Headers */
#include "basis/ABShellPairCalculator.h"
#include "data/ElectronicStructure.h"
#include "potentials/ABFockMatrixConstruction/ABCoulombInteractionPotential.h"
#include "potentials/ABFockMatrixConstruction/ABExchangePotential.h"
#include "potentials/ABFockMatrixConstruction/ABLRExchangePotential.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ABERIPotential<SCFMode>::ABERIPotential(std::shared_ptr<SystemController> system,
                                        std::shared_ptr<BasisController> basisA, std::shared_ptr<BasisController> basisB,
                                        std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> dMats,
                                        double exchangeRatio, double LRexchangeRatio, double mu, bool topDown,
                                        Options::DENS_FITS densFitJ, std::shared_ptr<BasisController> auxBasisAB,
                                        std::vector<std::shared_ptr<BasisController>> envAuxBasisController)
  : ABPotential<SCFMode>(basisA, basisB), _exchangeRatio(exchangeRatio), _lrExchangeRatio(LRexchangeRatio), _mu(mu) {
  // Basis
  this->_basisA->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_basisB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  // Densities
  for (const auto& densityController : dMats) {
    densityController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  }
  // Build Coulomb and exchange parts.
  if (_exchangeRatio != 0.0) {
    _abExchange = std::make_shared<ABExchangePotential<SCFMode>>(system, this->_basisA, this->_basisB, dMats, _exchangeRatio);
  }
  _abCoulomb = std::make_shared<ABCoulombInteractionPotential<SCFMode>>(system, this->_basisA, this->_basisB, dMats, topDown,
                                                                        densFitJ, auxBasisAB, envAuxBasisController);
  if (_lrExchangeRatio != 0.0) {
    _abLRExchange = std::make_shared<ABLRExchangePotential<SCFMode>>(system, this->_basisA, this->_basisB, dMats,
                                                                     _lrExchangeRatio, _mu);
  }
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode>& ABERIPotential<SCFMode>::getMatrix() {
  if (!_abPotential) {
    const unsigned int nBasisA = this->_basisA->getNBasisFunctions();
    const unsigned int nBasisB = this->_basisB->getNBasisFunctions();
    _abPotential.reset(new SPMatrix<SCFMode>(nBasisA, nBasisB));
    SPMatrix<SCFMode>& f_AB = *_abPotential;
    if (_exchangeRatio != 0.0)
      f_AB += _abExchange->getMatrix();
    f_AB += _abCoulomb->getMatrix();
    if (_abLRExchange)
      f_AB += _abLRExchange->getMatrix();
  } /* if !_abPotential */
  return *_abPotential;
}

template class ABERIPotential<Options::SCF_MODES::RESTRICTED>;
template class ABERIPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
