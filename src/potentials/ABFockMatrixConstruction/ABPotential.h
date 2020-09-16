/**
 * @file ABPotential.h
 *
 * @date May 8, 2018
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

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABPOTENTIAL_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/matrices/SPMatrix.h"

namespace Serenity {
/**
 * @class ABPotential ABPotential.h
 * @brief An interface for all outer diagonal fock matrix blocks.
 */
template<Options::SCF_MODES SCFMode>
class ABPotential {
 public:
  /**
   * @brief Constructor.
   * @param basisA The first basis.
   * @param basisB The second basis.
   */
  ABPotential(std::shared_ptr<BasisController> basisA, std::shared_ptr<BasisController> basisB)
    : _basisA(basisA), _basisB(basisB) {
  }
  /**
   * @brief Default destructor.
   */
  virtual ~ABPotential() = default;

  virtual SPMatrix<SCFMode>& getMatrix() = 0;

  /**
   * @param Small function for sanity checks.
   * @param basisA The first basis to compare with.
   * @param basisB The second basis to compare with.
   * @return Boolean, true if the basis sets match.
   */
  virtual bool compareBasis(std::shared_ptr<BasisController> basisA, std::shared_ptr<BasisController> basisB) {
    bool AEqA(basisA == _basisA);
    bool BEqB(basisB == _basisB);
    return AEqA && BEqB;
  };

 protected:
  ///@brief One of the basis sets this potential is defined in.
  std::shared_ptr<BasisController> _basisA;
  ///@brief The second of the basis sets this potential is defined in.
  std::shared_ptr<BasisController> _basisB;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABPOTENTIAL_H_ */
