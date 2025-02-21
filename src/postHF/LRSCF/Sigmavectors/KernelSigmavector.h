/**
 * @file KernelSigmavector.h
 *
 * @date Dec 07, 2018
 * @author Niklas Niemeyer, Johannes Toelle
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
#include "data/SpinPolarizedData.h"
#include "data/grid/GridPotential.h"
#include "math/Derivatives.h"
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

class BasisFunctionOnGridController;

template<Options::SCF_MODES SCFMode>
class Kernel;

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
   * @param ukernel An unrestricted kernel. Used to calculate triplet sigmavector in the restricted case.
   */
  KernelSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b,
                    std::shared_ptr<Kernel<SCFMode>> kernel, std::shared_ptr<Kernel<UNRESTRICTED>> ukernel = nullptr);

  /**
   * @brief Constructor if only the AO representation is required.
   * @param lrscf A set of LRSCF Controllers.
   * @param kernel The XC-Kernel holding the second functional derivatives of the XC energy at each gridpoint.
   */
  KernelSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::shared_ptr<Kernel<SCFMode>> kernel);

  /**
   * @brief Destructor.
   */
  virtual ~KernelSigmavector() = default;

  ///@brief Function to calculate and return Fock-like matrix F_IJ.
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned I, unsigned J, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_J) final;

 private:
  /**
   * @brief Calculates (weighted)
   *          pbb = \sum_{kl} P_{k l} \phi_k \phi_l
   * and
   *          pnbb = \sum_{kl} P_{k l} \nabla(\phi_k \phi_l)
   * for each density matrix and contract with kernel.
   */
  void contractKernel(std::vector<std::vector<MatrixInBasis<SCFMode>>>& dens,
                      std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                      std::vector<std::vector<GridPotential<SCFMode>>>& scalarPart,
                      std::vector<std::vector<Gradient<GridPotential<SCFMode>>>>& gradientPart, unsigned I, unsigned int J);

  ///@brief Performs the contraction of pbb/pnbb with the kernel stored on the grid.
  void contractBlock(const unsigned iGridStart, const unsigned blockSize, SpinPolarizedData<SCFMode, Eigen::VectorXd>& pbb,
                     Gradient<SpinPolarizedData<SCFMode, Eigen::VectorXd>>& pnbb, GridPotential<SCFMode>& scalarPart,
                     Gradient<GridPotential<SCFMode>>& gradientPart, unsigned I, unsigned J);

  ///@brief Performs numerical integration.
  void numericalIntegration(std::vector<std::vector<std::vector<MatrixInBasis<SCFMode>>>>& Fxc,
                            std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                            std::vector<std::vector<GridPotential<SCFMode>>>& scalarPart,
                            std::vector<std::vector<Gradient<GridPotential<SCFMode>>>>& gradientPart);

  ///@brief Contracts supersystem density (which needs to be done only once).
  void contractSupersystemDensity(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf);

  ///@brief Underlying kernel object.
  std::shared_ptr<Kernel<SCFMode>> _kernel;
  std::shared_ptr<Kernel<UNRESTRICTED>> _ukernel;

  ///@brief Density thresholds for prescreening.
  const double _blockAveThreshold;

  const bool _isGGA;

  std::vector<std::vector<GridPotential<SCFMode>>> _supersystem_scalar;
  std::vector<std::vector<Gradient<GridPotential<SCFMode>>>> _supersystem_gradient;

}; // class KernelSigmavector
} // namespace Serenity

#endif /* LRSCF_KERNELSIGMAVECTOR */
