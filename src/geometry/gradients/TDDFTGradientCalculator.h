/**
 * @file   TDDFTGradientCalculator.h
 *
 * @date   Mar 11, 2024
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
#ifndef GEOMETRY_GRADIENTS_TDDFTGRADIENTCALCULATOR_H_
#define GEOMETRY_GRADIENTS_TDDFTGRADIENTCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "data/DoublySpinPolarizedData.h"
#include "data/matrices/DensityMatrix.h"
#include "dft/functionals/wrappers/PartialDerivatives.h"
#include "math/Derivatives.h"
#include "postHF/LRSCF/Tools/SigmaCalculator.h"

namespace Serenity {
/* Forward Declarations */
template<Options::SCF_MODES>
class HyperkernelSigmavector;
template<Options::SCF_MODES>
class LRSCFController;
template<Options::SCF_MODES>
class Kernel;
template<Options::SCF_MODES>
class ScalarOperatorToMatrixAdder;
class BasisFunctionOnGridController;
template<Options::SCF_MODES>
class MatrixInBasis;

/**
 * @class TDDFTGradientCalculator TDDFTGradientCalculator.h
 * @brief Calculates analytic nuclear Cartesian gradients for excited states obtained from a prior LRSCFTask using
 * TDDFT, TDA, TDHF or CIS.
 * Reference: [1] F. Furche, R. Ahlrichs, J. Chem. Phys. 117, 16 (2002)
 */
template<Options::SCF_MODES SCFMode>
class TDDFTGradientCalculator {
 public:
  /**
   * @brief Constructor.
   * @param lrscf LRSCFController with access to the excitation energies and excitation vectors.
   * @param hypThresh A screening threshold for third derivatives of the exchange--correlation kernel. To increase the
   * numerical stability, the third derivatives are set to zero in places where the density falls below this threshold.
   */
  TDDFTGradientCalculator(std::shared_ptr<LRSCFController<SCFMode>> lrscf, double hypThresh = 1e-12);
  /**
   * @brief Calculates the nuclear gradients for all the excited states that are specified in the LRSCFTask's
   * excGradList keyword.
   * @returns The excited-state gradients as (nAtoms, 3) matrices that already include the groundstate contribution.
   */
  std::vector<Eigen::MatrixXd> calculateGradients();

 private:
  // After the Z-vector equation has been solved, this function takes all necessary matrices and contracts them with the
  // integral derivatives to return the gradient.
  const Eigen::MatrixXd evaluateGradients(const MatrixInBasis<SCFMode>& PAO, const MatrixInBasis<SCFMode>& WAO,
                                          const MatrixInBasis<SCFMode>& D, const MatrixInBasis<SCFMode>& XpYAO,
                                          const MatrixInBasis<SCFMode>& XmYAO, double omega);

  // evaluates Coulomb and exchange gradients using full 4-index integrals given the density matrices as arguments
  const Eigen::MatrixXd evaluateFullTwoCenterGradient(const MatrixInBasis<SCFMode>& PAO, const MatrixInBasis<SCFMode>& D,
                                                      const MatrixInBasis<SCFMode>& xpyAO,
                                                      const MatrixInBasis<SCFMode>& xmyAO);

  // evaluates exchange gradients using full 4-index integrals given the density matrices as arguments
  const Eigen::MatrixXd evaluateFullExchangeGradient(const MatrixInBasis<SCFMode>& PAO, const MatrixInBasis<SCFMode>& D,
                                                     const MatrixInBasis<SCFMode>& xpyAO, const MatrixInBasis<SCFMode>& xmyAO);

  // evaluates long-range exchange gradients using full 4-index integrals given the density matrices as arguments
  const Eigen::MatrixXd evaluateFullLRExchangeGradient(const MatrixInBasis<SCFMode>& PAO, const MatrixInBasis<SCFMode>& D,
                                                       const MatrixInBasis<SCFMode>& xpyAO,
                                                       const MatrixInBasis<SCFMode>& xmyAO);

  // evaluates gradient contributions of the first functional derivative of the exchange--correlation functional
  const Eigen::MatrixXd evaluateVxcGradientContribution(const DensityMatrix<SCFMode>& densitymat);

  // contract d2f/drho2 with the density and store the result in scalar - this is for LDA functionals
  void prepareFxcGradient(GridData<SCFMode>& scalar, const std::shared_ptr<d2F_dRho2<SCFMode>> pp,
                          const GridData<SCFMode>& dens);

  // contract the density and the density gradient with the mixed second Cartesian functional derivatives of the xc
  // functional - this is for GGA functionals
  void prepareFxcGradient(GridData<SCFMode>& scalar, Gradient<GridData<SCFMode>>& gradient,
                          const std::shared_ptr<d2F_dRho2<SCFMode>> pp,
                          const std::shared_ptr<Gradient<DoublySpinPolarizedData<SCFMode, GridData<RESTRICTED>>>> pg,
                          const std::shared_ptr<Hessian<DoublySpinPolarizedData<SCFMode, GridData<RESTRICTED>>>> gg,
                          const GridData<SCFMode>& dens, const Gradient<GridData<SCFMode>>& densg);

  // evaluates gradient contributions of the second functional derivative of the exchange--correlation functional
  const Eigen::MatrixXd evaluateFxcGradientContribution(const DensityMatrix<SCFMode>& densitymat,
                                                        const GridData<SCFMode>& densityOnGrid,
                                                        const Gradient<GridData<SCFMode>>& gradient);

  // evaluates gradient contributions of the third functional derivative of the exchange--correlation functional (the
  // "preparation" is done in the HyperkernelSigmavector)
  const Eigen::MatrixXd evaluateKxcGradientContribution(const DensityMatrix<SCFMode>& densitymat,
                                                        const GridData<SCFMode>& densityOnGrid,
                                                        const Gradient<GridData<SCFMode>>& gradient);
  // linear transformations (see [1] eq. 20a and eq. 20b)
  MatrixInBasis<SCFMode> hPlus(MatrixInBasis<SCFMode> V);
  MatrixInBasis<SCFMode> hMinus(MatrixInBasis<SCFMode> V);
  // the LRSFController used to construct this
  std::shared_ptr<LRSCFController<SCFMode>> _lrscf;
  // the density threshold for the hyperkernel
  double _hypThresh;
  // lambda function to calculate matrix-vector products in the iterative solution of the Z-Vector equation
  NonlinearSigmaCalculator _nlSolver;
  SigmaCalculator _lSolver;
  bool _usesXC = false;
  std::shared_ptr<Kernel<SCFMode>> _kernel = nullptr;
  std::shared_ptr<HyperkernelSigmavector<SCFMode>> _hyperSig = nullptr;
  // the excitations specified in the LRSCFTask
  std::vector<unsigned int> _excGradList;
  // the density functional used here
  Functional _func;
  std::shared_ptr<FunctionalData<SCFMode>> _funcData = nullptr;
  std::shared_ptr<BasisFunctionOnGridController> _basisFunctionOnGridController = nullptr;
  // A ScalarOperatorToMatrixAdder which performs numerical integration on the grid for XC integrals and integral
  // derivatives.
  std::shared_ptr<ScalarOperatorToMatrixAdder<SCFMode>> _gridToMatrix = nullptr;
  Eigen::MatrixXd _groundstateGradient;
};

} /* namespace Serenity */

#endif /* GEOMETRY_GRADIENTS_TDDFTGRADIENTCALCULATOR_H_ */
