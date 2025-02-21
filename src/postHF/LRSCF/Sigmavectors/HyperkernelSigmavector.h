/**
 * @file HyperkernelSigmavector.h
 *
 * @date Apr 30, 2024
 * @author Anton Rikus
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

#ifndef LRSCF_HYPERKERNELSIGMAVECTOR
#define LRSCF_HYPERKERNELSIGMAVECTOR

/* Include Serenity Internal Headers */
#include "data/grid/GridData.h"
#include "data/matrices/MatrixInBasis.h"
#include "math/Derivatives.h" // Gradient class
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

class BasisFunctionOnGridController;

template<Options::SCF_MODES SCFMode>
class FunctionalData;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

template<Options::SCF_MODES SCFMode>
class HyperkernelSigmavector {
 public:
  HyperkernelSigmavector(std::shared_ptr<LRSCFController<SCFMode>> lrscf,
                         std::shared_ptr<FunctionalData<SCFMode>> funcData, double screeningThreshold);

  /**
   * @brief Destructor.
   */
  virtual ~HyperkernelSigmavector() = default;

  ///@brief Function to calculate and return Fock-like matrix F_IJ.
  std::unique_ptr<MatrixInBasis<SCFMode>> calcF(std::shared_ptr<MatrixInBasis<SCFMode>> P);

  /**
   * @brief Calculates (weighted)
   *          pbb = \sum_{kl} P_{k l} \phi_k \phi_l
   * and
   *          pnbb = \sum_{kl} P_{k l} \nabla(\phi_k \phi_l)
   * for each density matrix and contract with kernel.
   */
  void contractKernel(GridData<SCFMode>& scalar, Gradient<GridData<SCFMode>>& gradient,
                      const GridData<SCFMode>& density, const Gradient<GridData<SCFMode>>& densitygrad);

 private:
  /// @brief Sets third derivatives to zero at points where the density falls below the _screeningThreshold.
  void screen();

  std::shared_ptr<LRSCFController<SCFMode>> _lrscf;

  ///@brief Underlying kernel object.
  std::shared_ptr<FunctionalData<SCFMode>> _funcData;

  const double _screeningThreshold;

  ///@brief Density thresholds for prescreening.
  const double _blockAveThreshold;

  const bool _isGGA = false;

  std::shared_ptr<GridController> _gridController;

  GridData<SCFMode> _scalar;
  ///@brief Number of threads used for Fock contractions.
  unsigned int _nThreads = omp_get_max_threads();

}; // class HyperkernelSigmavector
} // namespace Serenity

#endif /* LRSCF_KERNELSIGMAVECTOR */
