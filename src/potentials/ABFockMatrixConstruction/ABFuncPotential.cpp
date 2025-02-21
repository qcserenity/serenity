/**
 * @file ABFuncPotential.cpp
 *
 * @date May 16, 2018
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

/* Include Class Header*/
#include "potentials/ABFockMatrixConstruction/ABFuncPotential.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ExternalDensityOnGridController.h"
#include "dft/functionals/FunctionalLibrary.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ABFuncPotential<SCFMode>::ABFuncPotential(std::shared_ptr<SystemController> activeSystem,
                                          std::shared_ptr<BasisController> basisA,
                                          std::shared_ptr<BasisController> basisB, std::shared_ptr<GridController> grid,
                                          std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> dMats,
                                          Functional functional)
  : ABPotential<SCFMode>(basisA, basisB), _actSystem(activeSystem), _densityMatrices(dMats), _grid(grid), _functional(functional) {
  // Setting the notifying system up.
  // Basis
  this->_basisA->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_basisB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  // density matrices
  for (const auto& dMat : _densityMatrices) {
    dMat->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  }
  // grid
  _grid->addSensitiveObject(ObjectSensitiveClass<Grid>::_self);

  auto basisFunctionsOnGridA = BasisFunctionOnGridControllerFactory::produce(activeSystem->getSettings(), this->_basisA, grid);
  auto basisFunctionsOnGridB = BasisFunctionOnGridControllerFactory::produce(activeSystem->getSettings(), this->_basisB, grid);

  _gridToMatrix_AB = std::make_shared<ScalarOperatorToMatrixAdder<SCFMode>>(
      basisFunctionsOnGridA, basisFunctionsOnGridB, activeSystem->getSettings().grid.blockAveThreshold);
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode>& ABFuncPotential<SCFMode>::getMatrix() {
  if (!_abPotential) {
    // build the density on grid
    std::shared_ptr<ExternalDensityOnGridController<SCFMode>> densOnGridController;
    bool isGGA = _functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA;
    auto activeSystem = _actSystem.lock();
    const auto& actSettings = activeSystem->getSettings();
    if (!isGGA) {
      auto densityOnGrid = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_grid));
      // Loop over density matrices and add their density on the grid.
      for (const auto& dMat : _densityMatrices) {
        auto basisFunctionsOnGridC =
            BasisFunctionOnGridControllerFactory::produce(actSettings, dMat->getDensityMatrix().getBasisController(), _grid);
        DensityOnGridCalculator<SCFMode> densOnGridCalc(basisFunctionsOnGridC, actSettings.grid.blockAveThreshold);
        *densityOnGrid += densOnGridCalc.calcDensityOnGrid(dMat->getDensityMatrix());
      }
      densOnGridController = std::make_shared<ExternalDensityOnGridController<SCFMode>>(densityOnGrid);
    }
    else {
      auto densityOnGrid = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_grid));
      auto totDensityGradOnGrid = makeGradientPtr<DensityOnGrid<SCFMode>>(_grid);
      // Loop over density matrices and add their density and gradient on the grid.
      for (const auto& dMat : _densityMatrices) {
        auto basisFunctionsOnGridC =
            BasisFunctionOnGridControllerFactory::produce(actSettings, dMat->getDensityMatrix().getBasisController(), _grid);
        DensityOnGridCalculator<SCFMode> densOnGridCalc(basisFunctionsOnGridC, actSettings.grid.blockAveThreshold);
        auto densityGradOnGridC = makeGradientPtr<DensityOnGrid<SCFMode>>(_grid);
        *densityOnGrid += densOnGridCalc.calcDensityAndGradientOnGrid(dMat->getDensityMatrix(), *densityGradOnGridC);
        *totDensityGradOnGrid += *densityGradOnGridC;
      }
      densOnGridController = std::make_shared<ExternalDensityOnGridController<SCFMode>>(densityOnGrid, totDensityGradOnGrid);
    }

    // initialize functional library
    FunctionalLibrary<SCFMode> flib(128);
    // Calculate data
    auto funcData = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, _functional, densOnGridController);
    // initialize fock matrix
    unsigned int nBasisA = this->_basisA->getNBasisFunctions();
    unsigned int nBasisB = this->_basisB->getNBasisFunctions();
    _abPotential.reset(new SPMatrix<SCFMode>(nBasisA, nBasisB));
    auto& f_AB = *_abPotential;
    // convert the scalar potential into a matrix
    if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::LDA) {
      _gridToMatrix_AB->addScalarOperatorToMatrix(f_AB, *funcData.dFdRho);
    }
    else if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
      _gridToMatrix_AB->addScalarOperatorToMatrix(f_AB, *funcData.dFdRho, *funcData.dFdGradRho);
    }
    else {
      assert(false && "Unsupported functional type!");
    }
  } /* if !_abPotential */
  return *_abPotential;
}

template class ABFuncPotential<Options::SCF_MODES::RESTRICTED>;
template class ABFuncPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
