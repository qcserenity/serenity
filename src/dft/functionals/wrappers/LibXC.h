/**
 * @file   LibXC.h
 *
 * @date   Feb 24, 2020
 * @author Jan P. Unsleber
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
#ifdef SERENITY_USE_LIBXC
#ifndef LIBXC_H_
#define LIBXC_H_

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

struct xc_func_type;

namespace Serenity {

/**
 * @class LibXC LibXC.h
 * @brief A wrapper for the LibXC library
 */
template<Options::SCF_MODES SCFMode>
class LibXC {
 public:
  /**
   * @brief Constructor for LibXC wrapper
   * @param maxBlockSize Maximum number of grid points per block for the blockwise evaluation
   */
  LibXC(unsigned int maxBlockSize);
  virtual ~LibXC() = default;

  /**
   * @brief Evaluates functional energy expression and derivatives according to specified type up to chosen order and
   * stores data in FunctionalData object.
   * @param type LibXC can either evaluate derivatives w.r.t. gradients,  gradient invariants or directly evaluate the
   * LDA/GGA potential.
   * @param functional A density functional.
   * @param densityOnGridController A controller for the density to be processed.
   * @param order Highest derivative to be evaluated. Currently, only derivatives up to 2nd order are available.
   */
  FunctionalData<SCFMode> calcData(FUNCTIONAL_DATA_TYPE type, const Functional functional,
                                   const std::shared_ptr<DensityOnGridController<SCFMode>> densityOnGridController,
                                   unsigned int order = 1);

  /**
   * @brief Returns the version number of Libxc.
   */
  static std::string version();

  /**
   * @brief Returns the publication describing Libxc.
   */
  static std::string reference();

  /**
   * @brief Returns the DOI of the Libxc publication.
   */
  static std::string referenceDOI();

 private:
  unsigned int _maxBlockSize;

  FRIEND_TEST(LibXCTest, Comparison_PBE_Unrestricted_Cartesian_Cross);
  FRIEND_TEST(LibXCTest, Comparison_PBE_Restricted_Cartesian_Cross);

  // Produces LibXC functional from basic functionals
  void getFunctional(Functional functional);

  // Returns size of block with blockIndex given the number of gridpoints and the total number of blocks
  unsigned int determineBlockSize(unsigned int blockIndex, unsigned int nPoints, unsigned int nBlocks);

  // Calculates gradient invariants
  Eigen::MatrixXd calculateSigma(const Gradient<DensityOnGrid<SCFMode>>& gradient, const unsigned int& iGridStart,
                                 const unsigned int& blocksize);

  // calculates quantities in terms of gradients from gradient invariants
  void complete(const FunctionalData<SCFMode>& f, const Gradient<DensityOnGrid<SCFMode>>& gradient,
                const unsigned int& firstIndex, const unsigned int& blockSize);

  // evaluates density functionals and derivatives, i.e. passes the necessary ingredients to Libxc and stores the output
  void eval(const FunctionalData<SCFMode>& funcData, const DensityOnGrid<SCFMode>& density,
            std::shared_ptr<Gradient<DensityOnGrid<SCFMode>>> gradient, const double& prefactor,
            const CompositeFunctionals::CLASSES ftype, const xc_func_type& func, unsigned int order);

  // Evaluates and returns energy from energy density
  double calcEnergy(std::shared_ptr<GridData<RESTRICTED>> epuv, const Eigen::VectorXd& weights);
};

} /* namespace Serenity */

#endif /* LIBXC_H_ */
#endif /* SERENITY_USE_LIBXC */
