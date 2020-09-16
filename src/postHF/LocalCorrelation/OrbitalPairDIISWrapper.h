/**
 * @file OrbitalPairDIISWrapper.h
 *
 * @date Jul 30, 2019
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

#ifndef POSTHF_LOCALCORRELATION_ORBITALPAIRDIISWRAPPER_H_
#define POSTHF_LOCALCORRELATION_ORBITALPAIRDIISWRAPPER_H_

/* Include Serenity Internal Headers */
#include "math/diis/DIIS.h" //The actual DIIS.
/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {

/* Forward Declarations */
class OrbitalPair;
class SingleSubstitution;

/**
 * @class OrbitalPairDIISWrapper OrbitalPairDIISWrapper.h
 * @brief A wrapper for the convergence acceleration of the amplitude optimization of
 *        correlation methods based directly on orbital pairs or singles.
 */
class OrbitalPairDIISWrapper {
 public:
  /**
   * @brief Constructor.
   * @param maxStore Maximum number of DIIS vectors stored.
   */
  OrbitalPairDIISWrapper(unsigned int maxStore);
  /**
   * @brief Performs a DIIS step.
   * @param energy The current energy/correlation energy.
   * @param orbitalPairs The orbital pairs.
   * @param singles The singles.
   *
   * Note: One or both vectors (singles/orbitalPairs) can be
   *       empty without causing problems.
   */
  void optimize(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs, std::vector<std::shared_ptr<SingleSubstitution>> singles);

 private:
  ///@brief The actual wrapped DIIS.
  DIIS _diis;
};

} /* namespace Serenity */

#endif /* POSTHF_LOCALCORRELATION_ORBITALPAIRDIISWRAPPER_H_ */
