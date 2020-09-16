/**
 * @file ABEmbeddedBundle.h
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

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABBUNDLES_ABEMBEDDEDBUNDLE_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABBUNDLES_ABEMBEDDEDBUNDLE_H_
/* Include Serenity Internal Headers */
#include "data/matrices/SPMatrix.h"                          //Definition of an spin polarized matrix.
#include "potentials/ABFockMatrixConstruction/ABPotential.h" //Definition of an ABPotential.
#include "settings/Options.h"                                //Spin polarization.
/* Include Std and External Headers */
#include <memory> //smrt_ptr

namespace Serenity {

/**
 * @class ABEmbeddedBundle ABEmbeddedBundle.h
 * @brief An ABEmbeddedBundle that sorts the various contributions of the off-diagonal fock matrix block.
 */
template<Options::SCF_MODES SCFMode>
class ABEmbeddedBundle {
 public:
  /**
   * @brief Constructor.
   * @param hcore                      The HCore part.
   * @param activeCoulomb              The coulomb interaction of the active system.
   * @param environmentCoulomb         The coulomb interaction of the environment systems.
   * @param activeExchangeCorrelation  The exchange--correlation part of the active system.
   * @param naddExchangeCorrelation    The exchange--correlation part of the environment systems.
   * @param naddKinetic                The non--additive kinetic part.
   */
  ABEmbeddedBundle(std::shared_ptr<ABPotential<SCFMode>> hcore, std::shared_ptr<ABPotential<SCFMode>> activeCoulomb,
                   std::shared_ptr<ABPotential<SCFMode>> environmentCoulomb,
                   std::shared_ptr<ABPotential<SCFMode>> activeExchangeCorrelation,
                   std::shared_ptr<ABPotential<SCFMode>> naddExchangeCorrelation,
                   std::shared_ptr<ABPotential<SCFMode>> naddKinetic);
  /**
   * @brief Getter for the total off-diagonal fock matrix block.
   * @return The off-diagonal fock matrix block.
   */
  SPMatrix<SCFMode> getABMatrix();

 private:
  // The different underlying ABPotentials.
  std::shared_ptr<ABPotential<SCFMode>> _hcore;
  std::shared_ptr<ABPotential<SCFMode>> _activeCoulomb;
  std::shared_ptr<ABPotential<SCFMode>> _environmentCoulomb;
  std::shared_ptr<ABPotential<SCFMode>> _activeExchangeCorrelation;
  std::shared_ptr<ABPotential<SCFMode>> _naddExchangeCorrelation;
  std::shared_ptr<ABPotential<SCFMode>> _naddKinetic;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABBUNDLES_ABEMBEDDEDBUNDLE_H_ */
