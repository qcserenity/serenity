/**
 * @file ABEmbeddedBundle.cpp
 *
 * @date 18 Aug 2019
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
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundle.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ABEmbeddedBundle<SCFMode>::ABEmbeddedBundle(std::shared_ptr<ABPotential<SCFMode>> hcore,
                                            std::shared_ptr<ABPotential<SCFMode>> activeCoulomb,
                                            std::shared_ptr<ABPotential<SCFMode>> environmentCoulomb,
                                            std::shared_ptr<ABPotential<SCFMode>> activeExchangeCorrelation,
                                            std::shared_ptr<ABPotential<SCFMode>> naddExchangeCorrelation,
                                            std::shared_ptr<ABPotential<SCFMode>> naddKinetic)
  : _hcore(hcore),
    _activeCoulomb(activeCoulomb),
    _environmentCoulomb(environmentCoulomb),
    _activeExchangeCorrelation(activeExchangeCorrelation),
    _naddExchangeCorrelation(naddExchangeCorrelation),
    _naddKinetic(naddKinetic) {
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode> ABEmbeddedBundle<SCFMode>::getABMatrix() {
  SPMatrix<SCFMode> f_AB = _hcore->getMatrix();
  f_AB += _activeCoulomb->getMatrix();
  f_AB += _environmentCoulomb->getMatrix();
  f_AB += _activeExchangeCorrelation->getMatrix();
  f_AB += _naddExchangeCorrelation->getMatrix();
  if (_naddKinetic)
    f_AB += _naddKinetic->getMatrix();
  return f_AB;
}

template class ABEmbeddedBundle<Options::SCF_MODES::RESTRICTED>;
template class ABEmbeddedBundle<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
