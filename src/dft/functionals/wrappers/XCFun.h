/**
 * @file   XCFun.h
 *
 * @date   Feb 3, 2017
 * @author M. Boeckers
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
#ifdef SERENITY_USE_XCFUN
#ifndef XCFUN_H_
#define XCFUN_H_

/* Include Serenity Internal Headers */
#include "data/DoublySpinPolarizedData.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridController.h"
#include "data/grid/GridPotential.h"
#include "dft/Functional.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "dft/functionals/wrappers/PartialDerivatives.h"
#include "math/Derivatives.h"
/* GTest friend macro without include */
#define FRIEND_TEST(test_case_name, test_name) friend class test_case_name##_##test_name##_Test
/* Include Std and External Headers */
#include <XCFun/xcfun.h>

namespace Serenity {

/**
 * @class XCFun XCFun.h
 * @brief A wrapper for the XCFun library
 */
template<Options::SCF_MODES SCFMode>
class XCFun {
 public:
  /**
   * @brief Constructor for XCFun wrapper
   * @param maxBlockSize Maximum number of grid points per block for the blockwise evaluation
   */
  XCFun(unsigned int maxBlockSize);
  virtual ~XCFun() = default;

  /**
   * @brief Evaluates functional energy expression and derivatives according to specified type up to chosen order and
   * stores data in FunctionalData object.
   * @param type XCFun can either evaluate derivatives w.r.t. gradients,  gradient invariants or directly evaluate the
   * LDA/GGA potential.
   * @param functional A density functional.
   * @param densityOnGridController A controller for the density to be processed.
   * @param order Highest derivative to be evaluated. Currently, only derivatives up to 2nd order are available.
   */
  FunctionalData<SCFMode> calcData(FUNCTIONAL_DATA_TYPE type, const Functional functional,
                                   const std::shared_ptr<DensityOnGridController<SCFMode>> densityOnGridController,
                                   unsigned int order = 1);

 private:
  unsigned int _maxBlockSize;

  FRIEND_TEST(XCFunTest, FunctionalAliases);

  // Produces XCFun functional from basic functionals
  xcfun_t* getFunctional(Functional functional);

  // Returns size of block with blockIndex given the number of gridpoints and the total number of blocks
  unsigned int determineBlockSize(unsigned int blockIndex, unsigned int nPoints, unsigned int nBlocks);

  // Calculates gradient invariants
  Eigen::MatrixXd calculateSigma(const Gradient<DensityOnGrid<SCFMode>>& gradient, const unsigned int& iGridStart,
                                 const unsigned int& blocksize);

  // Prepare XCFun input from Serenity objects
  void prepareInput(const unsigned int iGridStart, const xcfun_vars xcVars, const DensityOnGrid<SCFMode>& density,
                    const std::shared_ptr<Gradient<DensityOnGrid<SCFMode>>> gradient,
                    const std::shared_ptr<Hessian<DensityOnGrid<SCFMode>>> hessian, Eigen::MatrixXd& input);

  // Stores XCFun output in FunctionalData object
  // from XCFun website: output is given in graded reverse lexicographical order
  void parseOutput(const unsigned int order, const xcfun_vars xcVars, const unsigned int iGridStart,
                   Eigen::MatrixXd& output, FunctionalData<SCFMode>& funcData);

  // Evaluates and returns energy from energy density
  double calcEnergy(std::shared_ptr<GridData<RESTRICTED>> epuv, const Eigen::VectorXd& weights);
};

} /* namespace Serenity */

#endif /* XCFUN_H_ */
#endif /* SERENITY_USE_XCFUN */
