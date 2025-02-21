/**
 * @file ABNAddFuncPotential.cpp
 *
 * @date May 17, 2018
 * @author Moritz Bensberg
 *
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
#include "potentials/ABFockMatrixConstruction/ABNAddFuncPotential.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ExternalDensityOnGridController.h"
#include "dft/functionals/FunctionalLibrary.h"
#include "dft/functionals/wrappers/PartialDerivatives.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ABNAddFuncPotential<SCFMode>::ABNAddFuncPotential(std::shared_ptr<SystemController> activeSystem,
                                                  std::shared_ptr<BasisController> basisA,
                                                  std::shared_ptr<BasisController> basisB,
                                                  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> dMats,
                                                  std::shared_ptr<GridController> grid, Functional functional)
  : ABPotential<SCFMode>(basisA, basisB), _actSystem(activeSystem), _densityMatrices(dMats), _grid(grid), _functional(functional) {
  // Basis
  this->_basisA->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_basisB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  // density matrices
  // active density
  activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->addSensitiveObject(
      ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  // Environment densities
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
SPMatrix<SCFMode>& ABNAddFuncPotential<SCFMode>::getMatrix() {
  auto activeSystem = _actSystem.lock();
  if (!_abPotential) {
    // Only recalculate the environment density if it is not available yet.
    // TODO: The environment density is not recalculated if it changes during the
    // calculation for the active system ... this may lead to errors.
    if (!_envDensityOnGrid) {
      // build env density on grid.
      bool isGGA = _functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA;
      const auto& actSettings = activeSystem->getSettings();
      if (!isGGA) {
        auto envDensityOnGrid = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_grid));
        // Loop over density matrices and add their density on the grid.
        for (const auto& dMat : _densityMatrices) {
          auto basisFunctionsOnGridC =
              BasisFunctionOnGridControllerFactory::produce(actSettings, dMat->getDensityMatrix().getBasisController(), _grid);
          DensityOnGridCalculator<SCFMode> densOnGridCalc(basisFunctionsOnGridC, actSettings.grid.blockAveThreshold);
          *envDensityOnGrid += densOnGridCalc.calcDensityOnGrid(dMat->getDensityMatrix());
        }
        _envDensityOnGrid = std::make_shared<ExternalDensityOnGridController<SCFMode>>(envDensityOnGrid);
      }
      else {
        auto envDensityOnGrid = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_grid));
        auto envDensityGradOnGrid = makeGradientPtr<DensityOnGrid<SCFMode>>(_grid);
        // Loop over density matrices and add their density and gradient on the grid.
        for (const auto& dMat : _densityMatrices) {
          auto basisFunctionsOnGridC =
              BasisFunctionOnGridControllerFactory::produce(actSettings, dMat->getDensityMatrix().getBasisController(), _grid);
          DensityOnGridCalculator<SCFMode> densOnGridCalc(basisFunctionsOnGridC, actSettings.grid.blockAveThreshold);
          auto densityGradOnGridC = makeGradientPtr<DensityOnGrid<SCFMode>>(_grid);
          *envDensityOnGrid += densOnGridCalc.calcDensityAndGradientOnGrid(dMat->getDensityMatrix(), *densityGradOnGridC);
          *envDensityGradOnGrid += *densityGradOnGridC;
        }
        _envDensityOnGrid = std::make_shared<ExternalDensityOnGridController<SCFMode>>(envDensityOnGrid, envDensityGradOnGrid);
      } /* else !isGGA */
    }   /* if !_envDensityOnGrid */

    // build active and supersystem density on grid
    std::shared_ptr<ExternalDensityOnGridController<SCFMode>> actDensOnGridController;
    std::shared_ptr<ExternalDensityOnGridController<SCFMode>> superDensOnGridController;
    bool isGGA = _functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA;
    const auto& actSettings = activeSystem->getSettings();
    if (!isGGA) {
      // initializing grids
      auto actDensityOnGrid = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_grid));
      auto superDensityOnGrid = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_grid));

      auto basisFunctionsOnGridAct = BasisFunctionOnGridControllerFactory::produce(
          actSettings, activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix().getBasisController(), _grid);
      DensityOnGridCalculator<SCFMode> densOnGridCalc(basisFunctionsOnGridAct, actSettings.grid.blockAveThreshold);
      // Calculating and adding density.
      *actDensityOnGrid =
          densOnGridCalc.calcDensityOnGrid(activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix());
      *superDensityOnGrid = _envDensityOnGrid->getDensityOnGrid();
      *superDensityOnGrid += *actDensityOnGrid;
      // Building controllers.
      actDensOnGridController = std::make_shared<ExternalDensityOnGridController<SCFMode>>(actDensityOnGrid);
      superDensOnGridController = std::make_shared<ExternalDensityOnGridController<SCFMode>>(superDensityOnGrid);
    }
    else {
      // initializing grids
      auto actDensityOnGrid = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_grid));
      auto actDensityGradOnGrid = makeGradientPtr<DensityOnGrid<SCFMode>>(_grid);
      auto superDensityOnGrid = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_grid));
      auto superDensityGradOnGrid = makeGradientPtr<DensityOnGrid<SCFMode>>(_grid);

      auto basisFunctionsOnGridC = BasisFunctionOnGridControllerFactory::produce(
          actSettings, activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix().getBasisController(), _grid);
      DensityOnGridCalculator<SCFMode> densOnGridCalc(basisFunctionsOnGridC, actSettings.grid.blockAveThreshold);
      // Calculating and adding density.
      *actDensityOnGrid = densOnGridCalc.calcDensityAndGradientOnGrid(
          activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix(), *actDensityGradOnGrid);
      *superDensityOnGrid = _envDensityOnGrid->getDensityOnGrid();
      *superDensityOnGrid += *actDensityOnGrid;
      *superDensityGradOnGrid = _envDensityOnGrid->getDensityGradientOnGrid();
      *superDensityGradOnGrid += *actDensityGradOnGrid;
      // Building controllers.
      actDensOnGridController =
          std::make_shared<ExternalDensityOnGridController<SCFMode>>(actDensityOnGrid, actDensityGradOnGrid);
      superDensOnGridController =
          std::make_shared<ExternalDensityOnGridController<SCFMode>>(superDensityOnGrid, superDensityGradOnGrid);
    } /* else !isGGA */

    // initializing the functional library
    FunctionalLibrary<SCFMode> flib(128);
    // Calculating supersystem part
    auto superFuncDat = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, _functional, superDensOnGridController, 1);
    // Calculating active system part
    auto activeFuncDat = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, _functional, actDensOnGridController, 1);
    // initialize fock matrix
    unsigned int nBasisA = this->_basisA->getNBasisFunctions();
    unsigned int nBasisB = this->_basisB->getNBasisFunctions();
    _abPotential.reset(new SPMatrix<SCFMode>(nBasisA, nBasisB));
    auto& f_AB = *_abPotential;
    // Calculating non-additive part.
    if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::NONE) {
      // Nothing to be done here
    }
    else if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::LDA) {
      *superFuncDat.dFdRho -= *activeFuncDat.dFdRho;
      _gridToMatrix_AB->addScalarOperatorToMatrix(f_AB, *superFuncDat.dFdRho);
    }
    else if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
      *superFuncDat.dFdRho -= *activeFuncDat.dFdRho;
      superFuncDat.dFdGradRho->x -= activeFuncDat.dFdGradRho->x;
      superFuncDat.dFdGradRho->y -= activeFuncDat.dFdGradRho->y;
      superFuncDat.dFdGradRho->z -= activeFuncDat.dFdGradRho->z;
      _gridToMatrix_AB->addScalarOperatorToMatrix(f_AB, *superFuncDat.dFdRho, *superFuncDat.dFdGradRho);
    }
    else {
      assert(false && "Unsupported functional type as nadd. functional!");
    }
  } /* if !_abPotential */
  return *_abPotential;
}

template class ABNAddFuncPotential<Options::SCF_MODES::RESTRICTED>;
template class ABNAddFuncPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
