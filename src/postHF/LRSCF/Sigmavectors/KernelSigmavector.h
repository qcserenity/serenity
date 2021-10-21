/**
 * @file KernelSigmavector.h
 *
 * @date Dec 07, 2018
 * @author Johannes Tölle
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

#ifndef LRSCF_KERNELSIGMAVECTOR
#define LRSCF_KERNELSIGMAVECTOR

/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "postHF/LRSCF/Kernel/Kernel.h"
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/**
 * @class KernelSigmavector KernelSigmavector.h
 * @brief Calculates the Kernel Sigma Vector. Performs a numerical integration of the XC-kernel with the basisfunctions
 *        of the subsystems.
 *
 * i.e.: \f$ \sigma_{ia} = \sum_{jb} M_{ia,jb} b_{jb} \f$
 *
 */
template<Options::SCF_MODES SCFMode>
class KernelSigmavector : public Sigmavector<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param lrscf A controller for the lrscf calculation
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of
   * guess vectors\n for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1],
   * respectively.\n For TDA-like problems, the guesses for X are stored in b[0].\n Note that sigma vectors can, albeit
   * not used in the present implementation, also be calculated for more than two sets of test vectors.
   * @param kernel The XC-kernel Object; holds the information of the kernel contributions stored on each grid point\n
   *        (calculated at the beginning of the TDDFT calculation).
   */
  KernelSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b,
                    std::shared_ptr<Kernel<SCFMode>> kernel);
  virtual ~KernelSigmavector() = default;

 private:
  // Function to calculate and return Fock-like matrix F_IJ
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned int I, unsigned int J, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_J) final;
  // Calculates
  //          pbb = \sum_{kl} P_{k l} \phi_k \phi_l
  // and
  //          pnbb = \sum_{kl} P_{k l} \nabla(\phi_k \phi_l)
  // Add weights!
  // for each density matrix and contract with kernel
  void contractKernel(std::vector<std::vector<MatrixInBasis<SCFMode>>>& dens,
                      std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                      std::vector<std::vector<GridPotential<SCFMode>>>& scalarPart,
                      std::vector<std::vector<Gradient<GridPotential<SCFMode>>>>& gradientPart, unsigned int I,
                      unsigned int J);
  // Performs the contraction of pbb/pnbb with the kernel stored on the grid
  void contractBlock(const unsigned int iGridStart, const unsigned int blockSize,
                     SpinPolarizedData<SCFMode, Eigen::VectorXd>& pbb,
                     Gradient<SpinPolarizedData<SCFMode, Eigen::VectorXd>>& pnbb, GridPotential<SCFMode>& scalarPart,
                     Gradient<GridPotential<SCFMode>>& gradientPart, unsigned int I, unsigned int J);
  // Performs numerical integration for GGAs
  void numIntSigma(std::vector<std::vector<std::vector<MatrixInBasis<SCFMode>>>>& focklikeMatrix,
                   std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                   std::vector<std::vector<GridPotential<SCFMode>>>& scalarPart,
                   std::vector<std::vector<Gradient<GridPotential<SCFMode>>>>& gradientPart, unsigned int I);
  // Performs numerical integration for LDAs
  void numIntSigma(std::vector<std::vector<std::vector<MatrixInBasis<SCFMode>>>>& focklikeMatrix,
                   std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                   std::vector<std::vector<GridPotential<SCFMode>>>& scalarPart, unsigned int I);

  // Kernel object
  std::shared_ptr<Kernel<SCFMode>> _kernel;

}; // class KernelSigmavector
} // namespace Serenity

#endif /* LRSCF_KERNELSIGMAVECTOR */
