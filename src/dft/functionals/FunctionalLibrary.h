/**
 * @file   FunctionalLibrary.h
 *
 * @date   Sep 3, 2020
 * @author Jan Unsleber, Thomas Weymuth
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
#ifndef FUNCTIONALLIBRARY_H_
#define FUNCTIONALLIBRARY_H_

/* Include Serenity Internal Headers */
#include "dft/functionals/wrappers/LibXC.h"
#include "dft/functionals/wrappers/XCFun.h"
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/* Forward decalartions */
template<Options::SCF_MODES SCFMode>
class FunctionalData;
template<Options::SCF_MODES SCFMode>
class DensityOnGridController;
class Functional;
enum class FUNCTIONAL_DATA_TYPE;

/**
 * @class FunctionalLibrary FunctionalLibrary.h
 * @brief An interface for LibXC and XCFun wrappers.
 */
template<Options::SCF_MODES SCFMode>
class FunctionalLibrary {
 public:
  /**
   * @brief Constructor.
   * @param maxBlockSize Maximum number of grid points per block for the blockwise evaluation
   */
  FunctionalLibrary(unsigned int maxBlockSize);
  virtual ~FunctionalLibrary() = default;

  /**
   * @brief Evaluates functional energy expression and derivatives according to specified type up to chosen order and
   *        stores data in FunctionalData object.
   * @param type The type of density derivatives to use in the evaluation.
   * @param functional A density functional.
   * @param densityOnGridController A controller for the density to be processed.
   * @param order Highest derivative to be evaluated. Currently, only derivatives up to 2nd order are available.
   */
  FunctionalData<SCFMode> calcData(FUNCTIONAL_DATA_TYPE type, const Functional functional,
                                   const std::shared_ptr<DensityOnGridController<SCFMode>> densityOnGridController,
                                   unsigned int order = 1);

 private:
#ifdef SERENITY_USE_LIBXC
  std::unique_ptr<LibXC<SCFMode>> _libxc;
#endif /* SERENITY_USE_LIBXC */
#ifdef SERENITY_USE_XCFUN
  std::unique_ptr<XCFun<SCFMode>> _xcfun;
#endif /* SERENITY_USE_XCFUN */
};

} /* namespace Serenity */

#endif /* FUNCTIONALLIBRARY_H_ */
