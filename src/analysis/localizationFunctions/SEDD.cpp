/**
 * @file SEDD.cpp
 *
 * @date Apr 4, 2016, reworked on Jul 12, 2017
 * @author Philipp Lenz, reworked by Jan Unsleber
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
#include "analysis/localizationFunctions/SEDD.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
std::function<Eigen::VectorXd(std::shared_ptr<GridController>)>
SEDD<SCFMode>::getSEDDLambda(const std::shared_ptr<SystemController> systemController) {
  auto calculateSEDD = [systemController](std::shared_ptr<GridController> gridController) {
    auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
        systemController->getSettings(), systemController->getAtomCenteredBasisController(), gridController);
    auto oldmaxderiv = basisFunctionOnGridController->getHighestDerivative();
    basisFunctionOnGridController->setHighestDerivative(2);
    auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
        basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
    const auto& densMatrix(systemController->getElectronicStructure<SCFMode>()->getDensityMatrix());
    DensityOnGrid<SCFMode> dens(gridController);
    // creates an object of type Gradient<densOnGrid<SCFMode>> and Hessian<densOnGrid<SCFMode>>
    auto densGradOnGrid(makeGradient<DensityOnGrid<SCFMode>>(gridController));
    auto densHessOnGrid(makeHessian<DensityOnGrid<SCFMode>>(gridController));

    Eigen::VectorXd result(gridController->getNGridPoints());
    densOnGridCalculator->calcDensityAndDerivativesOnGrid(densMatrix, dens, densGradOnGrid, densHessOnGrid);
    Eigen::VectorXd scalarGradient(gridController->getNGridPoints());
    Eigen::VectorXd tmp(gridController->getNGridPoints());
    scalarGradient.array() = densGradOnGrid.x.total().array() * densGradOnGrid.x.total().array() +
                             densGradOnGrid.y.total().array() * densGradOnGrid.y.total().array() +
                             densGradOnGrid.z.total().array() * densGradOnGrid.z.total().array();
    tmp.array() = (-1.0 * densGradOnGrid.x.total().array() * scalarGradient.array() +
                   dens.total().array() * (densGradOnGrid.x.total().array() * densHessOnGrid.xx.total().array() +
                                           densGradOnGrid.y.total().array() * densHessOnGrid.xy.total().array() +
                                           densGradOnGrid.z.total().array() * densHessOnGrid.xz.total().array()));
    result.array() = tmp.array() * tmp.array();
    tmp.array() = (-1.0 * densGradOnGrid.y.total().array() * scalarGradient.array() +
                   dens.total().array() * (densGradOnGrid.x.total().array() * densHessOnGrid.xy.total().array() +
                                           densGradOnGrid.y.total().array() * densHessOnGrid.yy.total().array() +
                                           densGradOnGrid.z.total().array() * densHessOnGrid.yz.total().array()));
    result.array() += tmp.array() * tmp.array();
    tmp.array() = (-1.0 * densGradOnGrid.z.total().array() * scalarGradient.array() +
                   dens.total().array() * (densGradOnGrid.x.total().array() * densHessOnGrid.xz.total().array() +
                                           densGradOnGrid.y.total().array() * densHessOnGrid.yz.total().array() +
                                           densGradOnGrid.z.total().array() * densHessOnGrid.zz.total().array()));
    result.array() += tmp.array() * tmp.array();
    result.array() *= 4.0 / dens.total().array().pow(8);
    result.array() = result.array().log();
    result.array() = result.array().unaryExpr([](double c) { return std::isfinite(c) ? c : 0.0; });
    basisFunctionOnGridController->setHighestDerivative(oldmaxderiv);
    return result;
  };
  return calculateSEDD;
}

template<Options::SCF_MODES SCFMode>
std::function<Eigen::VectorXd(std::shared_ptr<GridController>)>
SEDD<SCFMode>::getDORILambda(const std::shared_ptr<SystemController> systemController) {
  auto calculateDORI = [systemController](std::shared_ptr<GridController> gridController) {
    auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
        systemController->getSettings(), systemController->getAtomCenteredBasisController(), gridController);
    auto oldmaxderiv = basisFunctionOnGridController->getHighestDerivative();
    basisFunctionOnGridController->setHighestDerivative(2);
    auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
        basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
    const auto& densMatrix(systemController->getElectronicStructure<SCFMode>()->getDensityMatrix());
    DensityOnGrid<SCFMode> dens(gridController);
    // creates an object of type Gradient<densOnGrid<SCFMode>> and Hessian<densOnGrid<SCFMode>>
    auto densGradOnGrid(makeGradient<DensityOnGrid<SCFMode>>(gridController));
    auto densHessOnGrid(makeHessian<DensityOnGrid<SCFMode>>(gridController));

    Eigen::VectorXd result(gridController->getNGridPoints());
    densOnGridCalculator->calcDensityAndDerivativesOnGrid(densMatrix, dens, densGradOnGrid, densHessOnGrid);
    Eigen::VectorXd scalarGradient(gridController->getNGridPoints());
    Eigen::VectorXd tmp(gridController->getNGridPoints());
    scalarGradient.array() = densGradOnGrid.x.total().array() * densGradOnGrid.x.total().array() +
                             densGradOnGrid.y.total().array() * densGradOnGrid.y.total().array() +
                             densGradOnGrid.z.total().array() * densGradOnGrid.z.total().array();
    tmp.array() = (-1.0 * densGradOnGrid.x.total().array() * scalarGradient.array() +
                   dens.total().array() * (densGradOnGrid.x.total().array() * densHessOnGrid.xx.total().array() +
                                           densGradOnGrid.y.total().array() * densHessOnGrid.xy.total().array() +
                                           densGradOnGrid.z.total().array() * densHessOnGrid.xz.total().array()));
    result.array() = tmp.array() * tmp.array();
    tmp.array() = (-1.0 * densGradOnGrid.y.total().array() * scalarGradient.array() +
                   dens.total().array() * (densGradOnGrid.x.total().array() * densHessOnGrid.xy.total().array() +
                                           densGradOnGrid.y.total().array() * densHessOnGrid.yy.total().array() +
                                           densGradOnGrid.z.total().array() * densHessOnGrid.yz.total().array()));
    result.array() += tmp.array() * tmp.array();
    tmp.array() = (-1.0 * densGradOnGrid.z.total().array() * scalarGradient.array() +
                   dens.total().array() * (densGradOnGrid.x.total().array() * densHessOnGrid.xz.total().array() +
                                           densGradOnGrid.y.total().array() * densHessOnGrid.yz.total().array() +
                                           densGradOnGrid.z.total().array() * densHessOnGrid.zz.total().array()));
    result.array() += tmp.array() * tmp.array();
    result.array() *= 4.0 / dens.total().array().pow(6);
    result.array() *= 1.0 / ((densGradOnGrid.x.total().array() / dens.total().array()).pow(6) +
                             (densGradOnGrid.y.total().array() / dens.total().array()).pow(6) +
                             (densGradOnGrid.z.total().array() / dens.total().array()).pow(6));
    result.array() *= 1.0 / (1.0 + result.array());
    result.array() = result.array().unaryExpr([](double c) { return std::isfinite(c) ? c : 0.0; });
    basisFunctionOnGridController->setHighestDerivative(oldmaxderiv);
    return result;
  };
  return calculateDORI;
}

template<Options::SCF_MODES SCFMode>
std::function<Eigen::VectorXd(std::shared_ptr<GridController>)>
SEDD<SCFMode>::getSignedDensityLambda(const std::shared_ptr<SystemController> systemController) {
  auto calculateSignedDensity = [systemController](std::shared_ptr<GridController> gridController) {
    auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
        systemController->getSettings(), systemController->getAtomCenteredBasisController(), gridController);
    auto oldmaxderiv = basisFunctionOnGridController->getHighestDerivative();
    basisFunctionOnGridController->setHighestDerivative(2);
    auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
        basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
    const auto& densMatrix(systemController->getElectronicStructure<SCFMode>()->getDensityMatrix());
    DensityOnGrid<SCFMode> dens(gridController);
    // creates an object of type Gradient<densOnGrid<SCFMode>> and Hessian<densOnGrid<SCFMode>>
    auto densGradOnGrid(makeGradient<DensityOnGrid<SCFMode>>(gridController));
    auto densHessOnGrid(makeHessian<DensityOnGrid<SCFMode>>(gridController));

    densOnGridCalculator->calcDensityAndDerivativesOnGrid(densMatrix, dens, densGradOnGrid, densHessOnGrid);

    auto result = dens.total();
    auto result_ptr = result.data();
    unsigned int nPoints = result.size();
    unsigned int nBlocks = omp_get_max_threads();
    auto xx = densHessOnGrid.xx.total();
    auto yy = densHessOnGrid.yy.total();
    auto zz = densHessOnGrid.zz.total();
    auto xx_ptr = xx.data();
    auto yy_ptr = yy.data();
    auto zz_ptr = zz.data();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      for (unsigned int i = start; i < start + n; i++) {
        double val = xx_ptr[i] + yy_ptr[i] + zz_ptr[i];
        result_ptr[i] *= (double(0) < val) - (val < double(0));
      }
    }
    basisFunctionOnGridController->setHighestDerivative(oldmaxderiv);
    return result;
  };
  return calculateSignedDensity;
}

template class SEDD<Options::SCF_MODES::RESTRICTED>;
template class SEDD<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
