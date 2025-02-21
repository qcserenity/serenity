/**
 * @file HyperkernelSigmavector.cpp
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

/* Include Class Header*/
#include "postHF/LRSCF/Sigmavectors/HyperkernelSigmavector.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/DensityMatrixController.h"
#include "dft/functionals/wrappers/PartialDerivatives.h" // FunctionalData
#include "geometry/GeometryAdderFactory.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h" // LRSCFTaskSettings
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HyperkernelSigmavector<SCFMode>::HyperkernelSigmavector(std::shared_ptr<LRSCFController<SCFMode>> lrscf,
                                                        std::shared_ptr<FunctionalData<SCFMode>> funcData,
                                                        double screeningThreshold)
  : _lrscf(lrscf),
    _funcData(funcData),
    _screeningThreshold(screeningThreshold),
    _blockAveThreshold(_lrscf->getLRSCFSettings().grid.blockAveThreshold),
    _isGGA(_funcData->dFdGradRho),
    _gridController(AtomCenteredGridControllerFactory::produce(
        GeometryAdderFactory::produce({_lrscf->getSys()}, _lrscf->getLRSCFSettings().subsystemgrid),
        _lrscf->getLRSCFSettings().grid, Options::GRID_PURPOSES::SMALL)),
    _scalar(_gridController) {
  this->screen();
}

template<>
void HyperkernelSigmavector<RESTRICTED>::screen() {
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _lrscf->getSys()->getSettings(), _lrscf->getSys()->getBasisController(), _gridController);
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<RESTRICTED>>(
      basisFunctionOnGridController, _lrscf->getSys()->getSettings().grid.blockAveThreshold);
  std::shared_ptr<DensityMatrixController<RESTRICTED>> densityMatrixController =
      std::make_shared<DensityMatrixController<RESTRICTED>>(_lrscf->getSys()->getActiveOrbitalController<RESTRICTED>(),
                                                            _lrscf->getSys()->getNOccupiedOrbitals<RESTRICTED>());
  auto densityOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(densityOnGridCalculator, densityMatrixController);
  const auto& density = densityOnGridController->getDensityOnGrid();

  for (unsigned int iPoint = 0; iPoint < _gridController->getNGridPoints(); iPoint++) {
    if (density[iPoint] < _screeningThreshold) {
      (*_funcData->d3FdRho3)(iPoint) = 0.0;
      if (_isGGA) {
        _funcData->d3FdRho2dGradRho->setZero(iPoint);
        _funcData->d3FdRhodGradRho2->setZero(iPoint);
        _funcData->d3FdGradRho3->setZero(iPoint);
      }
    }
  }
}

template<>
void HyperkernelSigmavector<UNRESTRICTED>::screen() {
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _lrscf->getSys()->getSettings(), _lrscf->getSys()->getBasisController(), _gridController);
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<UNRESTRICTED>>(
      basisFunctionOnGridController, _lrscf->getSys()->getSettings().grid.blockAveThreshold);
  std::shared_ptr<DensityMatrixController<UNRESTRICTED>> densityMatrixController =
      std::make_shared<DensityMatrixController<UNRESTRICTED>>(_lrscf->getSys()->getActiveOrbitalController<UNRESTRICTED>(),
                                                              _lrscf->getSys()->getNOccupiedOrbitals<UNRESTRICTED>());
  auto densityOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<UNRESTRICTED>>(densityOnGridCalculator, densityMatrixController);
  const auto& density = densityOnGridController->getDensityOnGrid();

  for (unsigned int iPoint = 0; iPoint < _gridController->getNGridPoints(); iPoint++) {
    if (density.alpha[iPoint] < _screeningThreshold) {
      _funcData->d3FdRho3->setAlphaZero(iPoint);
      if (_isGGA) {
        _funcData->d3FdRho2dGradRho->setAlphaZero(iPoint);
        _funcData->d3FdRhodGradRho2->setAlphaZero(iPoint);
        _funcData->d3FdGradRho3->setAlphaZero(iPoint);
      }
    }
    if (density.beta[iPoint] < _screeningThreshold) {
      _funcData->d3FdRho3->setBetaZero(iPoint);
      if (_isGGA) {
        _funcData->d3FdRho2dGradRho->setBetaZero(iPoint);
        _funcData->d3FdRhodGradRho2->setBetaZero(iPoint);
        _funcData->d3FdGradRho3->setBetaZero(iPoint);
      }
    }
  }
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<MatrixInBasis<SCFMode>>
HyperkernelSigmavector<SCFMode>::calcF(std::shared_ptr<MatrixInBasis<SCFMode>> densityMatrix) {
  // Set dimensions for Fock-like matrices.
  auto fock = std::make_unique<MatrixInBasis<SCFMode>>(this->_lrscf->getBasisController());

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      this->_lrscf->getLRSCFSettings().grid.blocksize, this->_lrscf->getLRSCFSettings().grid.basFuncRadialThreshold,
      _isGGA ? 1 : 0, this->_lrscf->getBasisController(), this->_gridController);

  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
      basisFunctionOnGridController, this->_lrscf->getLRSCFSettings().grid.blockAveThreshold);

  const DensityMatrix<SCFMode>& D = (*densityMatrix);

  std::shared_ptr<DensityMatrixController<SCFMode>> dmatcontroller(std::make_shared<DensityMatrixController<SCFMode>>(D));
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densityOnGridCalculator, dmatcontroller);

  // calculates  \rho(r) = \sum_{\mu\nu} \chi_\mu(r) \chi_\nu(r) D_{\mu\nu} for each gridpoint
  const auto& density = densOnGridController->getDensityOnGrid();

  GridData<SCFMode> scalar(this->_gridController);
  Gradient<GridData<SCFMode>> gradient(makeGradient<GridData<SCFMode>>(this->_gridController)),
      densitygradient(makeGradient<GridData<SCFMode>>(this->_gridController));
  if (_isGGA) {
    densitygradient = densOnGridController->getDensityGradientOnGrid();
  }
  // calculates scalar(r) = \sum_{\sigma'\sigma''} \rho^X_\sigma'(r) \rho^X_\sigma''(r) w_i(r)
  // g^xc_{\sigma\sigma'\sigma''}(r)
  this->contractKernel(scalar, gradient, density, densitygradient);

  std::shared_ptr<ScalarOperatorToMatrixAdder<SCFMode>> gridToMatrix =
      std::make_shared<ScalarOperatorToMatrixAdder<SCFMode>>(basisFunctionOnGridController,
                                                             this->_lrscf->getLRSCFSettings().grid.blockAveThreshold);
  // calculates F_{\mu\nu} = \int dr \chi_\mu(r) \chi_\nu(r) _scalar(r)
  if (_isGGA) {
    gridToMatrix->addScalarOperatorToMatrix(*fock, scalar, gradient);
  }
  else {
    gridToMatrix->addScalarOperatorToMatrix(*fock, scalar);
  }

  return fock;
}

template<>
void HyperkernelSigmavector<RESTRICTED>::contractKernel(GridData<RESTRICTED>& scalar, Gradient<GridData<RESTRICTED>>& gradient,
                                                        const GridData<RESTRICTED>& d,
                                                        const Gradient<GridData<RESTRICTED>>& dg) {
  // first term
  scalar = 0.5 * d.cwiseProduct(d.cwiseProduct(*_funcData->d3FdRho3));
  if (_isGGA) {
    const auto& ppg = *_funcData->d3FdRho2dGradRho;
    const auto& pgg = *_funcData->d3FdRhodGradRho2;
    const auto& ggg = *_funcData->d3FdGradRho3;

    // second term
    scalar += d.array() * (dg.x.array() * ppg.x.array() + dg.y.array() * ppg.y.array() + dg.z.array() * ppg.z.array());

    gradient.x += Eigen::VectorXd(0.5 * d.cwiseProduct(d.cwiseProduct(ppg.x)));
    gradient.y += Eigen::VectorXd(0.5 * d.cwiseProduct(d.cwiseProduct(ppg.y)));
    gradient.z += Eigen::VectorXd(0.5 * d.cwiseProduct(d.cwiseProduct(ppg.z)));

    // third term
    scalar += Eigen::VectorXd(0.5 * dg.x.cwiseProduct(dg.x.cwiseProduct(pgg.xx)));
    scalar += Eigen::VectorXd(0.5 * dg.y.cwiseProduct(dg.y.cwiseProduct(pgg.yy)));
    scalar += Eigen::VectorXd(0.5 * dg.z.cwiseProduct(dg.z.cwiseProduct(pgg.zz)));
    scalar += Eigen::VectorXd(dg.x.cwiseProduct(dg.y.cwiseProduct(pgg.xy)));
    scalar += Eigen::VectorXd(dg.x.cwiseProduct(dg.z.cwiseProduct(pgg.xz)));
    scalar += Eigen::VectorXd(dg.y.cwiseProduct(dg.z.cwiseProduct(pgg.yz)));

    // fourth Term
    gradient.x += Eigen::VectorXd(dg.x.cwiseProduct(d.cwiseProduct(pgg.xx)));
    gradient.x += Eigen::VectorXd(dg.y.cwiseProduct(d.cwiseProduct(pgg.xy)));
    gradient.x += Eigen::VectorXd(dg.z.cwiseProduct(d.cwiseProduct(pgg.xz)));
    gradient.y += Eigen::VectorXd(dg.x.cwiseProduct(d.cwiseProduct(pgg.xy)));
    gradient.y += Eigen::VectorXd(dg.y.cwiseProduct(d.cwiseProduct(pgg.yy)));
    gradient.y += Eigen::VectorXd(dg.z.cwiseProduct(d.cwiseProduct(pgg.yz)));
    gradient.z += Eigen::VectorXd(dg.x.cwiseProduct(d.cwiseProduct(pgg.xz)));
    gradient.z += Eigen::VectorXd(dg.y.cwiseProduct(d.cwiseProduct(pgg.yz)));
    gradient.z += Eigen::VectorXd(dg.z.cwiseProduct(d.cwiseProduct(pgg.zz)));

    // fifth Term
    gradient.x += Eigen::VectorXd(0.5 * dg.x.cwiseProduct(dg.x.cwiseProduct(ggg.xxx)));
    gradient.x += Eigen::VectorXd(0.5 * dg.y.cwiseProduct(dg.y.cwiseProduct(ggg.xyy)));
    gradient.x += Eigen::VectorXd(0.5 * dg.z.cwiseProduct(dg.z.cwiseProduct(ggg.xzz)));
    gradient.x += Eigen::VectorXd(dg.x.cwiseProduct(dg.y.cwiseProduct(ggg.xxy)));
    gradient.x += Eigen::VectorXd(dg.x.cwiseProduct(dg.z.cwiseProduct(ggg.xxz)));
    gradient.x += Eigen::VectorXd(dg.y.cwiseProduct(dg.z.cwiseProduct(ggg.xyz)));

    gradient.y += Eigen::VectorXd(0.5 * dg.x.cwiseProduct(dg.x.cwiseProduct(ggg.xxy)));
    gradient.y += Eigen::VectorXd(0.5 * dg.y.cwiseProduct(dg.y.cwiseProduct(ggg.yyy)));
    gradient.y += Eigen::VectorXd(0.5 * dg.z.cwiseProduct(dg.z.cwiseProduct(ggg.yzz)));
    gradient.y += Eigen::VectorXd(dg.x.cwiseProduct(dg.y.cwiseProduct(ggg.xyy)));
    gradient.y += Eigen::VectorXd(dg.x.cwiseProduct(dg.z.cwiseProduct(ggg.xyz)));
    gradient.y += Eigen::VectorXd(dg.y.cwiseProduct(dg.z.cwiseProduct(ggg.yyz)));

    gradient.z += Eigen::VectorXd(0.5 * dg.x.cwiseProduct(dg.x.cwiseProduct(ggg.xxz)));
    gradient.z += Eigen::VectorXd(0.5 * dg.y.cwiseProduct(dg.y.cwiseProduct(ggg.yyz)));
    gradient.z += Eigen::VectorXd(0.5 * dg.z.cwiseProduct(dg.z.cwiseProduct(ggg.zzz)));
    gradient.z += Eigen::VectorXd(dg.x.cwiseProduct(dg.y.cwiseProduct(ggg.xyz)));
    gradient.z += Eigen::VectorXd(dg.x.cwiseProduct(dg.z.cwiseProduct(ggg.xzz)));
    gradient.z += Eigen::VectorXd(dg.y.cwiseProduct(dg.z.cwiseProduct(ggg.yzz)));
  }
}

template<>
void HyperkernelSigmavector<UNRESTRICTED>::contractKernel(GridData<UNRESTRICTED>& scalar,
                                                          Gradient<GridData<UNRESTRICTED>>& gradient,
                                                          const GridData<UNRESTRICTED>& d,
                                                          const Gradient<GridData<UNRESTRICTED>>& dg) {
  const auto& hyperkernel = *_funcData->d3FdRho3;
  scalar.alpha = 0.5 * (d.alpha.cwiseProduct(d.alpha.cwiseProduct(hyperkernel.aaa)) +
                        2 * d.alpha.cwiseProduct(d.beta.cwiseProduct(hyperkernel.aab)) +
                        d.beta.cwiseProduct(d.beta.cwiseProduct(hyperkernel.abb)));
  scalar.beta = 0.5 * (d.alpha.cwiseProduct(d.alpha.cwiseProduct(hyperkernel.aab)) +
                       2 * d.alpha.cwiseProduct(d.beta.cwiseProduct(hyperkernel.abb)) +
                       d.beta.cwiseProduct(d.beta.cwiseProduct(hyperkernel.bbb)));
  if (_isGGA) {
    const auto& ppg = *_funcData->d3FdRho2dGradRho;
    const auto& pgg = *_funcData->d3FdRhodGradRho2;
    const auto& ggg = *_funcData->d3FdGradRho3;

    // second term
    scalar.alpha.array() +=
        d.alpha.array() * (dg.x.alpha.array() * ppg.x.aaa.array() + dg.y.alpha.array() * ppg.y.aaa.array() +
                           dg.z.alpha.array() * ppg.z.aaa.array() + dg.x.beta.array() * ppg.x.aab.array() +
                           dg.y.beta.array() * ppg.y.aab.array() + dg.z.beta.array() * ppg.z.aab.array()) +
        d.beta.array() * (dg.x.alpha.array() * ppg.x.aba.array() + dg.y.alpha.array() * ppg.y.aba.array() +
                          dg.z.alpha.array() * ppg.z.aba.array() + dg.x.beta.array() * ppg.x.abb.array() +
                          dg.y.beta.array() * ppg.y.abb.array() + dg.z.beta.array() * ppg.z.abb.array());

    scalar.beta.array() +=
        d.alpha.array() * (dg.x.alpha.array() * ppg.x.aba.array() + dg.y.alpha.array() * ppg.y.aba.array() +
                           dg.z.alpha.array() * ppg.z.aba.array() + dg.x.beta.array() * ppg.x.abb.array() +
                           dg.y.beta.array() * ppg.y.abb.array() + dg.z.beta.array() * ppg.z.abb.array()) +
        d.beta.array() * (dg.x.alpha.array() * ppg.x.bba.array() + dg.y.alpha.array() * ppg.y.bba.array() +
                          dg.z.alpha.array() * ppg.z.bba.array() + dg.x.beta.array() * ppg.x.bbb.array() +
                          dg.y.beta.array() * ppg.y.bbb.array() + dg.z.beta.array() * ppg.z.bbb.array());

    gradient.x.alpha.array() += 0.5 * (d.alpha.array() * d.alpha.array() * ppg.x.aaa.array() +
                                       2 * d.alpha.array() * d.beta.array() * ppg.x.aba.array() +
                                       d.beta.array() * d.beta.array() * ppg.x.bba.array());
    gradient.x.beta.array() += 0.5 * (d.alpha.array() * d.alpha.array() * ppg.x.aab.array() +
                                      2 * d.alpha.array() * d.beta.array() * ppg.x.abb.array() +
                                      d.beta.array() * d.beta.array() * ppg.x.bbb.array());
    gradient.y.alpha.array() += 0.5 * (d.alpha.array() * d.alpha.array() * ppg.y.aaa.array() +
                                       2 * d.alpha.array() * d.beta.array() * ppg.y.aba.array() +
                                       d.beta.array() * d.beta.array() * ppg.y.bba.array());
    gradient.y.beta.array() += 0.5 * (d.alpha.array() * d.alpha.array() * ppg.y.aab.array() +
                                      2 * d.alpha.array() * d.beta.array() * ppg.y.abb.array() +
                                      d.beta.array() * d.beta.array() * ppg.y.bbb.array());
    gradient.z.alpha.array() += 0.5 * (d.alpha.array() * d.alpha.array() * ppg.z.aaa.array() +
                                       2 * d.alpha.array() * d.beta.array() * ppg.z.aba.array() +
                                       d.beta.array() * d.beta.array() * ppg.z.bba.array());
    gradient.z.beta.array() += 0.5 * (d.alpha.array() * d.alpha.array() * ppg.z.aab.array() +
                                      2 * d.alpha.array() * d.beta.array() * ppg.z.abb.array() +
                                      d.beta.array() * d.beta.array() * ppg.z.bbb.array());

    // third term
    scalar.alpha.array() += 0.5 * (dg.x.alpha.array() * dg.x.alpha.array() * pgg.xx.aaa.array() +
                                   2 * dg.x.alpha.array() * dg.x.beta.array() * pgg.xx.aab.array() +
                                   dg.x.beta.array() * dg.x.beta.array() * pgg.xx.abb.array());
    scalar.beta.array() += 0.5 * (dg.x.alpha.array() * dg.x.alpha.array() * pgg.xx.baa.array() +
                                  2 * dg.x.alpha.array() * dg.x.beta.array() * pgg.xx.bab.array() +
                                  dg.x.beta.array() * dg.x.beta.array() * pgg.xx.bbb.array());
    scalar.alpha.array() += 0.5 * (dg.y.alpha.array() * dg.y.alpha.array() * pgg.yy.aaa.array() +
                                   2 * dg.y.alpha.array() * dg.y.beta.array() * pgg.yy.aab.array() +
                                   dg.y.beta.array() * dg.y.beta.array() * pgg.yy.abb.array());
    scalar.beta.array() += 0.5 * (dg.y.alpha.array() * dg.y.alpha.array() * pgg.yy.baa.array() +
                                  2 * dg.y.alpha.array() * dg.y.beta.array() * pgg.yy.bab.array() +
                                  dg.y.beta.array() * dg.y.beta.array() * pgg.yy.bbb.array());
    scalar.alpha.array() += 0.5 * (dg.z.alpha.array() * dg.z.alpha.array() * pgg.zz.aaa.array() +
                                   2 * dg.z.alpha.array() * dg.z.beta.array() * pgg.zz.aab.array() +
                                   dg.z.beta.array() * dg.z.beta.array() * pgg.zz.abb.array());
    scalar.beta.array() += 0.5 * (dg.z.alpha.array() * dg.z.alpha.array() * pgg.zz.baa.array() +
                                  2 * dg.z.alpha.array() * dg.z.beta.array() * pgg.zz.bab.array() +
                                  dg.z.beta.array() * dg.z.beta.array() * pgg.zz.bbb.array());
    scalar.alpha.array() += dg.x.alpha.array() * dg.y.alpha.array() * pgg.xy.aaa.array() +
                            dg.x.alpha.array() * dg.y.beta.array() * pgg.xy.aab.array() +
                            dg.x.beta.array() * dg.y.alpha.array() * pgg.xy.aba.array() +
                            dg.x.beta.array() * dg.y.beta.array() * pgg.xy.abb.array();
    scalar.beta.array() += (dg.x.alpha.array() * dg.y.alpha.array() * pgg.xy.baa.array() +
                            dg.x.alpha.array() * dg.y.beta.array() * pgg.xy.bab.array() +
                            dg.x.beta.array() * dg.y.alpha.array() * pgg.xy.bba.array() +
                            dg.x.beta.array() * dg.y.beta.array() * pgg.xy.bbb.array());
    scalar.alpha.array() += (dg.x.alpha.array() * dg.z.alpha.array() * pgg.xz.aaa.array() +
                             dg.x.alpha.array() * dg.z.beta.array() * pgg.xz.aab.array() +
                             dg.x.beta.array() * dg.z.alpha.array() * pgg.xz.aba.array() +
                             dg.x.beta.array() * dg.z.beta.array() * pgg.xz.abb.array());
    scalar.beta.array() += (dg.x.alpha.array() * dg.z.alpha.array() * pgg.xz.baa.array() +
                            dg.x.alpha.array() * dg.z.beta.array() * pgg.xz.bab.array() +
                            dg.x.beta.array() * dg.z.alpha.array() * pgg.xz.bba.array() +
                            dg.x.beta.array() * dg.z.beta.array() * pgg.xz.bbb.array());
    scalar.alpha.array() += (dg.y.alpha.array() * dg.z.alpha.array() * pgg.yz.aaa.array() +
                             dg.y.alpha.array() * dg.z.beta.array() * pgg.yz.aab.array() +
                             dg.y.beta.array() * dg.z.alpha.array() * pgg.yz.aba.array() +
                             dg.y.beta.array() * dg.z.beta.array() * pgg.yz.abb.array());
    scalar.beta.array() += (dg.y.alpha.array() * dg.z.alpha.array() * pgg.yz.baa.array() +
                            dg.y.alpha.array() * dg.z.beta.array() * pgg.yz.bab.array() +
                            dg.y.beta.array() * dg.z.alpha.array() * pgg.yz.bba.array() +
                            dg.y.beta.array() * dg.z.beta.array() * pgg.yz.bbb.array());

    // fourth Term
    gradient.x.alpha.array() +=
        (dg.x.alpha.array() * d.alpha.array() * pgg.xx.aaa.array() + dg.x.alpha.array() * d.beta.array() * pgg.xx.baa.array() +
         dg.x.beta.array() * d.alpha.array() * pgg.xx.aab.array() + dg.x.beta.array() * d.beta.array() * pgg.xx.bab.array());
    gradient.x.beta.array() +=
        (dg.x.alpha.array() * d.alpha.array() * pgg.xx.aab.array() + dg.x.alpha.array() * d.beta.array() * pgg.xx.bab.array() +
         dg.x.beta.array() * d.alpha.array() * pgg.xx.abb.array() + dg.x.beta.array() * d.beta.array() * pgg.xx.bbb.array());
    gradient.x.alpha.array() +=
        (dg.y.alpha.array() * d.alpha.array() * pgg.xy.aaa.array() + dg.y.alpha.array() * d.beta.array() * pgg.xy.baa.array() +
         dg.y.beta.array() * d.alpha.array() * pgg.xy.aab.array() + dg.y.beta.array() * d.beta.array() * pgg.xy.bab.array());
    gradient.x.beta.array() +=
        (dg.y.alpha.array() * d.alpha.array() * pgg.xy.aba.array() + dg.y.alpha.array() * d.beta.array() * pgg.xy.bba.array() +
         dg.y.beta.array() * d.alpha.array() * pgg.xy.abb.array() + dg.y.beta.array() * d.beta.array() * pgg.xy.bbb.array());
    gradient.x.alpha.array() +=
        (dg.z.alpha.array() * d.alpha.array() * pgg.xz.aaa.array() + dg.z.alpha.array() * d.beta.array() * pgg.xz.baa.array() +
         dg.z.beta.array() * d.alpha.array() * pgg.xz.aab.array() + dg.z.beta.array() * d.beta.array() * pgg.xz.bab.array());
    gradient.x.beta.array() +=
        (dg.z.alpha.array() * d.alpha.array() * pgg.xz.aba.array() + dg.z.alpha.array() * d.beta.array() * pgg.xz.bba.array() +
         dg.z.beta.array() * d.alpha.array() * pgg.xz.abb.array() + dg.z.beta.array() * d.beta.array() * pgg.xz.bbb.array());
    gradient.y.alpha.array() +=
        (dg.x.alpha.array() * d.alpha.array() * pgg.xy.aaa.array() + dg.x.alpha.array() * d.beta.array() * pgg.xy.baa.array() +
         dg.x.beta.array() * d.alpha.array() * pgg.xy.aba.array() + dg.x.beta.array() * d.beta.array() * pgg.xy.bba.array());
    gradient.y.beta.array() +=
        (dg.x.alpha.array() * d.alpha.array() * pgg.xy.aab.array() + dg.x.alpha.array() * d.beta.array() * pgg.xy.bab.array() +
         dg.x.beta.array() * d.alpha.array() * pgg.xy.abb.array() + dg.x.beta.array() * d.beta.array() * pgg.xy.bbb.array());
    gradient.y.alpha.array() +=
        (dg.y.alpha.array() * d.alpha.array() * pgg.yy.aaa.array() + dg.y.alpha.array() * d.beta.array() * pgg.yy.baa.array() +
         dg.y.beta.array() * d.alpha.array() * pgg.yy.aab.array() + dg.y.beta.array() * d.beta.array() * pgg.yy.bab.array());
    gradient.y.beta.array() +=
        (dg.y.alpha.array() * d.alpha.array() * pgg.yy.aab.array() + dg.y.alpha.array() * d.beta.array() * pgg.yy.bab.array() +
         dg.y.beta.array() * d.alpha.array() * pgg.yy.abb.array() + dg.y.beta.array() * d.beta.array() * pgg.yy.bbb.array());
    gradient.y.alpha.array() +=
        (dg.z.alpha.array() * d.alpha.array() * pgg.yz.aaa.array() + dg.z.alpha.array() * d.beta.array() * pgg.yz.baa.array() +
         dg.z.beta.array() * d.alpha.array() * pgg.yz.aab.array() + dg.z.beta.array() * d.beta.array() * pgg.yz.bab.array());
    gradient.y.beta.array() +=
        (dg.z.alpha.array() * d.alpha.array() * pgg.yz.aba.array() + dg.z.alpha.array() * d.beta.array() * pgg.yz.bba.array() +
         dg.z.beta.array() * d.alpha.array() * pgg.yz.abb.array() + dg.z.beta.array() * d.beta.array() * pgg.yz.bbb.array());
    gradient.z.alpha.array() +=
        (dg.x.alpha.array() * d.alpha.array() * pgg.xz.aaa.array() + dg.x.alpha.array() * d.beta.array() * pgg.xz.baa.array() +
         dg.x.beta.array() * d.alpha.array() * pgg.xz.aab.array() + dg.x.beta.array() * d.beta.array() * pgg.xz.bab.array());
    gradient.z.beta.array() +=
        (dg.x.alpha.array() * d.alpha.array() * pgg.xz.aba.array() + dg.x.alpha.array() * d.beta.array() * pgg.xz.bba.array() +
         dg.x.beta.array() * d.alpha.array() * pgg.xz.abb.array() + dg.x.beta.array() * d.beta.array() * pgg.xz.bbb.array());
    gradient.z.alpha.array() +=
        (dg.y.alpha.array() * d.alpha.array() * pgg.yz.aaa.array() + dg.y.alpha.array() * d.beta.array() * pgg.yz.baa.array() +
         dg.y.beta.array() * d.alpha.array() * pgg.yz.aab.array() + dg.y.beta.array() * d.beta.array() * pgg.yz.bab.array());
    gradient.z.beta.array() +=
        (dg.y.alpha.array() * d.alpha.array() * pgg.yz.aba.array() + dg.y.alpha.array() * d.beta.array() * pgg.yz.bba.array() +
         dg.y.beta.array() * d.alpha.array() * pgg.yz.abb.array() + dg.y.beta.array() * d.beta.array() * pgg.yz.bbb.array());
    gradient.z.alpha.array() +=
        (dg.z.alpha.array() * d.alpha.array() * pgg.zz.aaa.array() + dg.z.alpha.array() * d.beta.array() * pgg.zz.baa.array() +
         dg.z.beta.array() * d.alpha.array() * pgg.zz.aab.array() + dg.z.beta.array() * d.beta.array() * pgg.zz.bab.array());
    gradient.z.beta.array() +=
        (dg.z.alpha.array() * d.alpha.array() * pgg.zz.aab.array() + dg.z.alpha.array() * d.beta.array() * pgg.zz.bab.array() +
         dg.z.beta.array() * d.alpha.array() * pgg.zz.abb.array() + dg.z.beta.array() * d.beta.array() * pgg.zz.bbb.array());

    // fifth Term
    gradient.x.alpha.array() += 0.5 * (dg.x.alpha.array() * dg.x.alpha.array() * ggg.xxx.aaa.array() +
                                       2 * dg.x.alpha.array() * dg.x.beta.array() * ggg.xxx.aab.array() +
                                       dg.x.beta.array() * dg.x.beta.array() * ggg.xxx.abb.array());
    gradient.x.beta.array() += 0.5 * (dg.x.alpha.array() * dg.x.alpha.array() * ggg.xxx.aab.array() +
                                      2 * dg.x.alpha.array() * dg.x.beta.array() * ggg.xxx.abb.array() +
                                      dg.x.beta.array() * dg.x.beta.array() * ggg.xxx.bbb.array());
    gradient.x.alpha.array() += 0.5 * (dg.y.alpha.array() * dg.y.alpha.array() * ggg.xyy.aaa.array() +
                                       2 * dg.y.alpha.array() * dg.y.beta.array() * ggg.xyy.aab.array() +
                                       dg.y.beta.array() * dg.y.beta.array() * ggg.xyy.abb.array());
    gradient.x.beta.array() += 0.5 * (dg.y.alpha.array() * dg.y.alpha.array() * ggg.xyy.baa.array() +
                                      2 * dg.y.alpha.array() * dg.y.beta.array() * ggg.xyy.bab.array() +
                                      dg.y.beta.array() * dg.y.beta.array() * ggg.xyy.bbb.array());
    gradient.x.alpha.array() += 0.5 * (dg.z.alpha.array() * dg.z.alpha.array() * ggg.xzz.aaa.array() +
                                       2 * dg.z.alpha.array() * dg.z.beta.array() * ggg.xzz.aab.array() +
                                       dg.z.beta.array() * dg.z.beta.array() * ggg.xzz.abb.array());
    gradient.x.beta.array() += 0.5 * (dg.z.alpha.array() * dg.z.alpha.array() * ggg.xzz.baa.array() +
                                      2 * dg.z.alpha.array() * dg.z.beta.array() * ggg.xzz.bab.array() +
                                      dg.z.beta.array() * dg.z.beta.array() * ggg.xzz.bbb.array());
    gradient.x.alpha.array() += (dg.x.alpha.array() * dg.y.alpha.array() * ggg.xxy.aaa.array() +
                                 dg.x.alpha.array() * dg.y.beta.array() * ggg.xxy.aab.array() +
                                 dg.x.beta.array() * dg.y.alpha.array() * ggg.xxy.aba.array() +
                                 dg.x.beta.array() * dg.y.beta.array() * ggg.xxy.abb.array());
    gradient.x.beta.array() += (dg.x.alpha.array() * dg.y.alpha.array() * ggg.xxy.aba.array() +
                                dg.x.alpha.array() * dg.y.beta.array() * ggg.xxy.abb.array() +
                                dg.x.beta.array() * dg.y.alpha.array() * ggg.xxy.bba.array() +
                                dg.x.beta.array() * dg.y.beta.array() * ggg.xxy.bbb.array());
    gradient.x.alpha.array() += (dg.x.alpha.array() * dg.z.alpha.array() * ggg.xxz.aaa.array() +
                                 dg.x.alpha.array() * dg.z.beta.array() * ggg.xxz.aab.array() +
                                 dg.x.beta.array() * dg.z.alpha.array() * ggg.xxz.aba.array() +
                                 dg.x.beta.array() * dg.z.beta.array() * ggg.xxz.abb.array());
    gradient.x.beta.array() += (dg.x.alpha.array() * dg.z.alpha.array() * ggg.xxz.aba.array() +
                                dg.x.alpha.array() * dg.z.beta.array() * ggg.xxz.abb.array() +
                                dg.x.beta.array() * dg.z.alpha.array() * ggg.xxz.bba.array() +
                                dg.x.beta.array() * dg.z.beta.array() * ggg.xxz.bbb.array());
    gradient.x.alpha.array() += (dg.y.alpha.array() * dg.z.alpha.array() * ggg.xyz.aaa.array() +
                                 dg.y.alpha.array() * dg.z.beta.array() * ggg.xyz.aab.array() +
                                 dg.y.beta.array() * dg.z.alpha.array() * ggg.xyz.aba.array() +
                                 dg.y.beta.array() * dg.z.beta.array() * ggg.xyz.abb.array());
    gradient.x.beta.array() += (dg.y.alpha.array() * dg.z.alpha.array() * ggg.xyz.aab.array() +
                                dg.y.alpha.array() * dg.z.beta.array() * ggg.xyz.bab.array() +
                                dg.y.beta.array() * dg.z.alpha.array() * ggg.xyz.bba.array() +
                                dg.y.beta.array() * dg.z.beta.array() * ggg.xyz.bbb.array());
    gradient.y.alpha.array() += 0.5 * (dg.x.alpha.array() * dg.x.alpha.array() * ggg.xxy.aaa.array() +
                                       2 * dg.x.alpha.array() * dg.x.beta.array() * ggg.xxy.aba.array() +
                                       dg.x.beta.array() * dg.x.beta.array() * ggg.xxy.bba.array());
    gradient.y.beta.array() += 0.5 * (dg.x.alpha.array() * dg.x.alpha.array() * ggg.xxy.aab.array() +
                                      2 * dg.x.alpha.array() * dg.x.beta.array() * ggg.xxy.abb.array() +
                                      dg.x.beta.array() * dg.x.beta.array() * ggg.xxy.bbb.array());
    gradient.y.alpha.array() += 0.5 * (dg.y.alpha.array() * dg.y.alpha.array() * ggg.yyy.aaa.array() +
                                       2 * dg.y.alpha.array() * dg.y.beta.array() * ggg.yyy.aab.array() +
                                       dg.y.beta.array() * dg.y.beta.array() * ggg.yyy.abb.array());
    gradient.y.beta.array() += 0.5 * (dg.y.alpha.array() * dg.y.alpha.array() * ggg.yyy.aab.array() +
                                      2 * dg.y.alpha.array() * dg.y.beta.array() * ggg.yyy.abb.array() +
                                      dg.y.beta.array() * dg.y.beta.array() * ggg.yyy.bbb.array());
    gradient.y.alpha.array() += 0.5 * (dg.z.alpha.array() * dg.z.alpha.array() * ggg.yzz.aaa.array() +
                                       2 * dg.z.alpha.array() * dg.z.beta.array() * ggg.yzz.aab.array() +
                                       dg.z.beta.array() * dg.z.beta.array() * ggg.yzz.abb.array());
    gradient.y.beta.array() += 0.5 * (dg.z.alpha.array() * dg.z.alpha.array() * ggg.yzz.baa.array() +
                                      2 * dg.z.alpha.array() * dg.z.beta.array() * ggg.yzz.bab.array() +
                                      dg.z.beta.array() * dg.z.beta.array() * ggg.yzz.bbb.array());
    gradient.y.alpha.array() += (dg.x.alpha.array() * dg.y.alpha.array() * ggg.xyy.aaa.array() +
                                 dg.x.alpha.array() * dg.y.beta.array() * ggg.xyy.aab.array() +
                                 dg.x.beta.array() * dg.y.alpha.array() * ggg.xyy.baa.array() +
                                 dg.x.beta.array() * dg.y.beta.array() * ggg.xyy.bab.array());
    gradient.y.beta.array() += (dg.x.alpha.array() * dg.y.alpha.array() * ggg.xyy.aab.array() +
                                dg.x.alpha.array() * dg.y.beta.array() * ggg.xyy.abb.array() +
                                dg.x.beta.array() * dg.y.alpha.array() * ggg.xyy.bab.array() +
                                dg.x.beta.array() * dg.y.beta.array() * ggg.xyy.bbb.array());
    gradient.y.alpha.array() += (dg.x.alpha.array() * dg.z.alpha.array() * ggg.xyz.aaa.array() +
                                 dg.x.alpha.array() * dg.z.beta.array() * ggg.xyz.aab.array() +
                                 dg.x.beta.array() * dg.z.alpha.array() * ggg.xyz.baa.array() +
                                 dg.x.beta.array() * dg.z.beta.array() * ggg.xyz.bab.array());
    gradient.y.beta.array() += (dg.x.alpha.array() * dg.z.alpha.array() * ggg.xyz.aba.array() +
                                dg.x.alpha.array() * dg.z.beta.array() * ggg.xyz.abb.array() +
                                dg.x.beta.array() * dg.z.alpha.array() * ggg.xyz.bba.array() +
                                dg.x.beta.array() * dg.z.beta.array() * ggg.xyz.bbb.array());
    gradient.y.alpha.array() += (dg.y.alpha.array() * dg.z.alpha.array() * ggg.yyz.aaa.array() +
                                 dg.y.alpha.array() * dg.z.beta.array() * ggg.yyz.aab.array() +
                                 dg.y.beta.array() * dg.z.alpha.array() * ggg.yyz.aba.array() +
                                 dg.y.beta.array() * dg.z.beta.array() * ggg.yyz.abb.array());
    gradient.y.beta.array() += (dg.y.alpha.array() * dg.z.alpha.array() * ggg.yyz.aba.array() +
                                dg.y.alpha.array() * dg.z.beta.array() * ggg.yyz.abb.array() +
                                dg.y.beta.array() * dg.z.alpha.array() * ggg.yyz.bba.array() +
                                dg.y.beta.array() * dg.z.beta.array() * ggg.yyz.bbb.array());
    gradient.z.alpha.array() += 0.5 * (dg.x.alpha.array() * dg.x.alpha.array() * ggg.xxz.aaa.array() +
                                       2 * dg.x.alpha.array() * dg.x.beta.array() * ggg.xxz.aba.array() +
                                       dg.x.beta.array() * dg.x.beta.array() * ggg.xxz.bba.array());
    gradient.z.beta.array() += 0.5 * (dg.x.alpha.array() * dg.x.alpha.array() * ggg.xxz.aab.array() +
                                      2 * dg.x.alpha.array() * dg.x.beta.array() * ggg.xxz.abb.array() +
                                      dg.x.beta.array() * dg.x.beta.array() * ggg.xxz.bbb.array());
    gradient.z.alpha.array() += 0.5 * (dg.y.alpha.array() * dg.y.alpha.array() * ggg.yyz.aaa.array() +
                                       2 * dg.y.alpha.array() * dg.y.beta.array() * ggg.yyz.aba.array() +
                                       dg.y.beta.array() * dg.y.beta.array() * ggg.yyz.bba.array());
    gradient.z.beta.array() += 0.5 * (dg.y.alpha.array() * dg.y.alpha.array() * ggg.yyz.aab.array() +
                                      2 * dg.y.alpha.array() * dg.y.beta.array() * ggg.yyz.abb.array() +
                                      dg.y.beta.array() * dg.y.beta.array() * ggg.yyz.bbb.array());
    gradient.z.alpha.array() += 0.5 * (dg.z.alpha.array() * dg.z.alpha.array() * ggg.zzz.aaa.array() +
                                       2 * dg.z.alpha.array() * dg.z.beta.array() * ggg.zzz.aab.array() +
                                       dg.z.beta.array() * dg.z.beta.array() * ggg.zzz.abb.array());
    gradient.z.beta.array() += 0.5 * (dg.z.alpha.array() * dg.z.alpha.array() * ggg.zzz.aab.array() +
                                      2 * dg.z.alpha.array() * dg.z.beta.array() * ggg.zzz.abb.array() +
                                      dg.z.beta.array() * dg.z.beta.array() * ggg.zzz.bbb.array());
    gradient.z.alpha.array() += (dg.x.alpha.array() * dg.y.alpha.array() * ggg.xyz.aaa.array() +
                                 dg.x.alpha.array() * dg.y.beta.array() * ggg.xyz.aba.array() +
                                 dg.x.beta.array() * dg.y.alpha.array() * ggg.xyz.baa.array() +
                                 dg.x.beta.array() * dg.y.beta.array() * ggg.xyz.bba.array());
    gradient.z.beta.array() += (dg.x.alpha.array() * dg.y.alpha.array() * ggg.xyz.aab.array() +
                                dg.x.alpha.array() * dg.y.beta.array() * ggg.xyz.abb.array() +
                                dg.x.beta.array() * dg.y.alpha.array() * ggg.xyz.bab.array() +
                                dg.x.beta.array() * dg.y.beta.array() * ggg.xyz.bbb.array());
    gradient.z.alpha.array() += (dg.x.alpha.array() * dg.z.alpha.array() * ggg.xzz.aaa.array() +
                                 dg.x.alpha.array() * dg.z.beta.array() * ggg.xzz.aab.array() +
                                 dg.x.beta.array() * dg.z.alpha.array() * ggg.xzz.baa.array() +
                                 dg.x.beta.array() * dg.z.beta.array() * ggg.xzz.bab.array());
    gradient.z.beta.array() += (dg.x.alpha.array() * dg.z.alpha.array() * ggg.xzz.aab.array() +
                                dg.x.alpha.array() * dg.z.beta.array() * ggg.xzz.abb.array() +
                                dg.x.beta.array() * dg.z.alpha.array() * ggg.xzz.bab.array() +
                                dg.x.beta.array() * dg.z.beta.array() * ggg.xzz.bbb.array());
    gradient.z.alpha.array() += (dg.y.alpha.array() * dg.z.alpha.array() * ggg.yzz.aaa.array() +
                                 dg.y.alpha.array() * dg.z.beta.array() * ggg.yzz.aab.array() +
                                 dg.y.beta.array() * dg.z.alpha.array() * ggg.yzz.baa.array() +
                                 dg.y.beta.array() * dg.z.beta.array() * ggg.yzz.bab.array());
    gradient.z.beta.array() += (dg.y.alpha.array() * dg.z.alpha.array() * ggg.yzz.aab.array() +
                                dg.y.alpha.array() * dg.z.beta.array() * ggg.yzz.abb.array() +
                                dg.y.beta.array() * dg.z.alpha.array() * ggg.yzz.bab.array() +
                                dg.y.beta.array() * dg.z.beta.array() * ggg.yzz.bbb.array());
  }
}

template class HyperkernelSigmavector<Options::SCF_MODES::RESTRICTED>;
template class HyperkernelSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
