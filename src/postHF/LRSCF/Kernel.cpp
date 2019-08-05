/**
 * @file   Kernel.cpp
 *
 * @date   Mar 19, 2017
 * @author M. Boeckers
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "postHF/LRSCF/Kernel.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/grid/DensityAdder.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/ElectronicStructure.h"
#include "data/grid/ExternalDensityOnGridController.h"
#include "tasks/FDETask.h"
#include "grid/GridControllerFactory.h"
#include "data/OrbitalController.h"
#include "geometry/Point.h"
#include "dft/functionals/wrappers/XCFun.h"


namespace Serenity {

template<Options::SCF_MODES T> Kernel<T>::Kernel(
    std::shared_ptr<SystemController> activeSystem,
    std::vector<std::shared_ptr<SystemController> > environmentSystems,
    bool superSystemGrid,
    bool noNaddKernel,
    Options::FUNCTIONALS func,
    Options::FUNCTIONALS naddKinFunc,
    Options::FUNCTIONALS naddXCFunc):
      _activeSystem(activeSystem),
      _environmentSystems(environmentSystems),
      _superSystemGrid(superSystemGrid),
      _noNaddKernel(noNaddKernel),
      _func(func),
      _naddKinFunc(naddKinFunc),
      _naddXCFunc(naddXCFunc),
      _fde((environmentSystems.size() > 0) ? true : false),
      _d2FdRho2(nullptr),
      _dFdSigma(nullptr),
      _d2FdSigma2(nullptr),
      _d2FdRhodSigma(nullptr),
      _totalNaddD2FdSigma2(nullptr),
      _totalNaddD2FdRhodSigma(nullptr),
      _activeDensityGradient(nullptr),
      _totalDensityGradient(nullptr),
      _gga(false),
      _naddGGA(false),
      _hasBeenCalculated(false){
  assert(_activeSystem);
  if (_fde) {
    for (unsigned int iSys = 0; iSys < _environmentSystems.size(); ++iSys) {
      assert(_environmentSystems[iSys]);
    }
  }


  //Get gridController
  if (_superSystemGrid && _fde) {
    //Use super-system grid (copied and adapted from FDETask.cpp)

    FDETask<T> fdeTask(_activeSystem,_environmentSystems);

    // atoms of all subsystems
    std::vector<std::shared_ptr<Atom> > superSystemAtoms;
    superSystemAtoms.insert(superSystemAtoms.end(), _activeSystem->getAtoms().begin(), _activeSystem->getAtoms().end());

    double cutoff = fdeTask.settings.gridCutOff;
    if (cutoff<0.0) cutoff = std::numeric_limits<double>::infinity();

    for (auto sys : _environmentSystems){
      for (auto atom : sys->getAtoms()){
        for (auto check : _activeSystem->getAtoms()){
          if (distance(*atom,*check)< cutoff) {
            superSystemAtoms.push_back(atom);
            break;
          }
        }
      }
    }
    //geometry of the entire system
    auto superSystemGeometry = std::make_shared<Geometry>(superSystemAtoms);

    //supersystem grid
    _gridController = GridControllerFactory::produce(
            superSystemGeometry, _activeSystem->getSettings(), Options::GRID_PURPOSES::DEFAULT);
  } else {
    //Use grid of active system (default)
    _gridController = _activeSystem->getGridController();
  }
  assert(_gridController);

  //initialize data objects
  _d2FdRho2 = std::make_shared<d2F_dRho2<T> >(_gridController);

  if(FunctionalClassResolver::resolveFunctional(_func).getFunctionalClass() == FUNCTIONAL_CLASSES::GGA) {
    _gga = true;
    _dFdSigma = std::make_shared<dF_dSigma<T> >(_gridController);
    _d2FdSigma2 = std::make_shared<d2F_dSigma2<T> >(_gridController);
    _d2FdRhodSigma = std::make_shared<d2F_dRhodSigma<T> >(_gridController);
  }
  if (_fde && !_noNaddKernel) {
    if(FunctionalClassResolver::resolveFunctional(_naddKinFunc).getFunctionalClass() == FUNCTIONAL_CLASSES::GGA ||
        FunctionalClassResolver::resolveFunctional(_naddXCFunc).getFunctionalClass() == FUNCTIONAL_CLASSES::GGA) {
      _naddGGA = true;
      _totalNaddD2FdSigma2 = std::make_shared<d2F_dSigma2<T> >(_gridController);
      _totalNaddD2FdRhodSigma = std::make_shared<d2F_dRhodSigma<T> >(_gridController);
    }
  }
}

template<Options::SCF_MODES T>
std::shared_ptr<GridController> Kernel<T>::getGridController(){
  assert(_gridController);
  return _gridController;
}

template<Options::SCF_MODES T>
std::shared_ptr<d2F_dRho2<T> > Kernel<T>::getD2F_dRho2() {
  if (!_hasBeenCalculated) calculateKernelOnGrid();
  return _d2FdRho2;
}

template<Options::SCF_MODES T>
std::shared_ptr<dF_dSigma<T> > Kernel<T>::getDF_dSigma() {
  //Returns nullptr if not calculated
  if (!_hasBeenCalculated) calculateKernelOnGrid();
  return _dFdSigma;
}

template<Options::SCF_MODES T>
std::shared_ptr<d2F_dSigma2<T> > Kernel<T>::getD2F_dSigma2() {
  //Returns nullptr if not calculated
  if (!_hasBeenCalculated) calculateKernelOnGrid();
  return _d2FdSigma2;
}

template<Options::SCF_MODES T>
std::shared_ptr<d2F_dRhodSigma<T> > Kernel<T>::getD2F_dRhodSigma() {
  //Returns nullptr if not calculated
  if (!_hasBeenCalculated) calculateKernelOnGrid();
  return _d2FdRhodSigma;
}

template<Options::SCF_MODES T>
std::shared_ptr<d2F_dSigma2<T> > Kernel<T>::getTotalNaddD2F_dSigma2() {
  //Returns nullptr if not calculated
  if (!_hasBeenCalculated) calculateKernelOnGrid();
  return _totalNaddD2FdSigma2;
}

template<Options::SCF_MODES T>
std::shared_ptr<d2F_dRhodSigma<T> > Kernel<T>::getTotalNaddD2F_dRhodSigma() {
  //Returns nullptr if not calculated
  if (!_hasBeenCalculated) calculateKernelOnGrid();
  return _totalNaddD2FdRhodSigma;
}

template<Options::SCF_MODES T>
std::shared_ptr<Gradient<DensityOnGrid<T> > > Kernel<T>::getActiveDensityGradient() {
  //Object is initialized in calculateKernelOnGridData() when density on grid of choice is known
  if (!_hasBeenCalculated) calculateKernelOnGrid();
  return _activeDensityGradient;
}

template<Options::SCF_MODES T>
std::shared_ptr<Gradient<DensityOnGrid<T> > > Kernel<T>::getTotalDensityGradient() {
  //Object is initialized in calculateKernelOnGridData() when density on grid of choice and total density are known
  if (!_hasBeenCalculated) calculateKernelOnGrid();
  return _totalDensityGradient;
}


template<>
void Kernel<Options::SCF_MODES::RESTRICTED>::calculateKernelOnGrid() {
  //Alias for better readability
  const auto R = Options::SCF_MODES::RESTRICTED;

  //calculate densities of active system on grid of choice and build density controller
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _activeSystem->getSettings(), _activeSystem->getBasisController(), _gridController);
  auto densityOnGridCalculator =
      std::make_shared<DensityOnGridCalculator<R> >(
          basisFunctionOnGridController, _activeSystem->getSettings().grid.blockAveThreshold);
  auto densityMatrixController = std::make_shared<DensityMatrixController<R> > (
      _activeSystem->getActiveOrbitalController<R>(),
      _activeSystem->getNOccupiedOrbitals<R>());
  auto activeDensityOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<R> >(
          densityOnGridCalculator,
          densityMatrixController);
  auto activeDensity = std::unique_ptr<DensityOnGrid<R> > (
      new DensityOnGrid<R>(_gridController));
  (*activeDensity) = activeDensityOnGridController->getDensityOnGrid();

  //Calculate gradient of active density on grid of choice if needed
  if (_gga || (_naddGGA && _fde)) {
    _activeDensityGradient = makeGradientPtr<DensityOnGrid<R> >(_gridController);
    (*_activeDensityGradient) = activeDensityOnGridController->getDensityGradientOnGrid();
  }

  XCFun<Options::SCF_MODES::RESTRICTED> xcfun(_activeSystem->getSettings().grid.blocksize);
  auto activeFuncData = xcfun.calcData(
      FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
      FunctionalClassResolver::resolveFunctional(_func),
      activeDensityOnGridController,
      2);

    (*_d2FdRho2) = (*activeFuncData.d2FdRho2);
    if (_gga) {
      (*_dFdSigma) = (*activeFuncData.dFdSigma);
      (*_d2FdSigma2) = (*activeFuncData.d2FdSigma2);
      (*_d2FdRhodSigma) = (*activeFuncData.d2FdRhodSigma);
    }

  //initialize FDE objects and add to kernel objects if kernel is requested
  if (_fde && !_noNaddKernel) {
    //calculate total density and gradient
    std::unique_ptr<DensityOnGrid<R> > totalDensity = std::unique_ptr<DensityOnGrid<R> > (
        new DensityOnGrid<R>(_gridController));
    (*totalDensity) = (*activeDensity);
    if (_naddGGA) {
      _totalDensityGradient = makeGradientPtr<DensityOnGrid<R> >(_gridController);
      (*_totalDensityGradient) = (*_activeDensityGradient);
    }
    for (auto sys : _environmentSystems) {
      if (sys->getLastSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
        assert(false && "UNRESTRICTED ENVIRONMENT FOR RESTRICTED FDE-TDDFT CALCULATIONS NOT IMPLEMENTED");
      } else {
        auto sysBasisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
            sys->getSettings(), sys->getBasisController(), _gridController);
        auto sysDensityOnGridCalculator =
            std::make_shared<DensityOnGridCalculator<R> >(
                sysBasisFunctionOnGridController, sys->getSettings().grid.blockAveThreshold);
        auto sysDensityMatrixController = std::make_shared<DensityMatrixController<R> > (
            sys->getActiveOrbitalController<R>(),
            sys->getNOccupiedOrbitals<R>());
        auto sysDensityOnGridController =
            std::make_shared<DensityMatrixDensityOnGridController<R> >(
                sysDensityOnGridCalculator,
                sysDensityMatrixController);
        (*totalDensity) += sysDensityOnGridController->getDensityOnGrid();
        if (_naddGGA) (*_totalDensityGradient) += sysDensityOnGridController->getDensityGradientOnGrid();
      }
    }
    //build DensityOnGridController for total density
    std::shared_ptr<ExternalDensityOnGridController<R> > totalDensityOnGridController;
    if (!_naddGGA) {
      totalDensityOnGridController = std::make_shared<ExternalDensityOnGridController<R> > (totalDensity);
    } else {
      auto totDens = makeGradientPtr<DensityOnGrid<R> >(_gridController);
      (*totDens) = (*_totalDensityGradient);
      totalDensityOnGridController = std::make_shared<ExternalDensityOnGridController<R> > (totalDensity,totDens);
    }

    auto activeNaddKinFuncData = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        FunctionalClassResolver::resolveFunctional(_naddKinFunc),
        activeDensityOnGridController,
        2);
    auto activeNaddXCFuncData = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        FunctionalClassResolver::resolveFunctional(_naddXCFunc),
        activeDensityOnGridController,
        2);
    auto totalNaddKinFuncData = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        FunctionalClassResolver::resolveFunctional(_naddKinFunc),
        totalDensityOnGridController,
        2);
    auto totalNaddXCFuncData = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        FunctionalClassResolver::resolveFunctional(_naddXCFunc),
        totalDensityOnGridController,
        2);

    (*_d2FdRho2) += (*totalNaddKinFuncData.d2FdRho2);
    (*_d2FdRho2) -= (*activeNaddKinFuncData.d2FdRho2);
    if (_naddGGA) {
      (*_dFdSigma) -= (*activeNaddKinFuncData.dFdSigma);
      (*_d2FdSigma2) -= (*activeNaddKinFuncData.d2FdSigma2);
      (*_d2FdRhodSigma) -= (*activeNaddKinFuncData.d2FdRhodSigma);
      (*_dFdSigma) += (*totalNaddKinFuncData.dFdSigma);
      (*_totalNaddD2FdSigma2) = (*totalNaddKinFuncData.d2FdSigma2);
      (*_totalNaddD2FdRhodSigma) = (*totalNaddKinFuncData.d2FdRhodSigma);
    }
    (*_d2FdRho2) += (*totalNaddXCFuncData.d2FdRho2);
    (*_d2FdRho2) -= (*activeNaddXCFuncData.d2FdRho2);
    if (_naddGGA) {
      (*_dFdSigma) -= (*activeNaddXCFuncData.dFdSigma);
      (*_d2FdSigma2) -= (*activeNaddXCFuncData.d2FdSigma2);
      (*_d2FdRhodSigma) -= (*activeNaddXCFuncData.d2FdRhodSigma);
      (*_dFdSigma) += (*totalNaddXCFuncData.dFdSigma);
      (*_totalNaddD2FdSigma2) += (*totalNaddXCFuncData.d2FdSigma2);
      (*_totalNaddD2FdRhodSigma) += (*totalNaddXCFuncData.d2FdRhodSigma);
    }
  }
  //Screen Kernel:
  //If the density is close to zero, p^-x goes to infinity and will cause numerical
  //problems. However, if the density at this grid points is low, also the molecular
  //orbitals at this point will be small so that the contribution to the XC kernel
  //integrals can be neglected. Since the kernel at this points goes to infinity,
  //while the molecular orbitals go to zero, it is not guaranteed that the integrals
  //will numerically be zero. Thus, we set the kernel to zero at positions
  //where the density falls below a threshold.
  //ToDo: Make threshold input parameter
  double screeningThreshold = 1.0e-8;
  for (unsigned int iPoint = 0; iPoint < _gridController->getNGridPoints(); ++iPoint) {
    if ((*activeDensity)[iPoint] < screeningThreshold) {
      (*_d2FdRho2)[iPoint] = 0.0;
      if (_gga) {
        (*_dFdSigma)[iPoint] = 0.0;
        (*_d2FdSigma2)[iPoint] = 0.0;
        (*_d2FdRhodSigma)[iPoint] = 0.0;
      }
      if (_fde && _naddGGA && !_noNaddKernel) {
        (*_totalNaddD2FdSigma2)[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma)[iPoint] = 0.0;
      }
    }
  }

  _hasBeenCalculated = true;
}



template<>
void Kernel<Options::SCF_MODES::UNRESTRICTED>::calculateKernelOnGrid() {
  //Alias for better readability
  const auto U = Options::SCF_MODES::UNRESTRICTED;

  //calculate densities of active system on grid of choice and build density controller
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _activeSystem->getSettings(), _activeSystem->getBasisController(), _gridController);
  auto densityOnGridCalculator =
      std::make_shared<DensityOnGridCalculator<U> >(
          basisFunctionOnGridController, _activeSystem->getSettings().grid.blockAveThreshold);
  auto densityMatrixController = std::make_shared<DensityMatrixController<U> > (
      _activeSystem->getActiveOrbitalController<U>(),
      _activeSystem->getNOccupiedOrbitals<U>());
  auto activeDensityOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<U> >(
          densityOnGridCalculator,
          densityMatrixController);
  auto activeDensity = std::unique_ptr<DensityOnGrid<U> > (
      new DensityOnGrid<U>(_gridController));
  (*activeDensity) = activeDensityOnGridController->getDensityOnGrid();

  //Calculate gradient of active density on grid of choice if needed
  if (_gga || (_naddGGA && _fde)) {
    _activeDensityGradient = makeGradientPtr<DensityOnGrid<U> >(_gridController);
    (*_activeDensityGradient) = activeDensityOnGridController->getDensityGradientOnGrid();
  }

  XCFun<Options::SCF_MODES::UNRESTRICTED> xcfun(128);
  auto activeFuncData = xcfun.calcData(
      FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
      FunctionalClassResolver::resolveFunctional(_func),
      activeDensityOnGridController,
      2);

    (*_d2FdRho2) = (*activeFuncData.d2FdRho2);
    if (_gga) {
      (*_dFdSigma) = (*activeFuncData.dFdSigma);
      (*_d2FdSigma2) = (*activeFuncData.d2FdSigma2);
      (*_d2FdRhodSigma) = (*activeFuncData.d2FdRhodSigma);
    }

  //initialize FDE objects and add to kernel objects if kernel is requested
  if (_fde && !_noNaddKernel) {
    //calculate total density and gradient
    std::unique_ptr<DensityOnGrid<U> > totalDensity = std::unique_ptr<DensityOnGrid<U> > (
        new DensityOnGrid<U>(_gridController));
    (*totalDensity) = (*activeDensity);
    if (_naddGGA) {
      _totalDensityGradient = makeGradientPtr<DensityOnGrid<U> >(_gridController);
      (*_totalDensityGradient) = (*_activeDensityGradient);
    }
    for (auto sys : _environmentSystems) {
      auto sysBasisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
          sys->getSettings(), sys->getBasisController(), _gridController);
      if (sys->getLastSCFMode() == Options::SCF_MODES::RESTRICTED) {
        auto rsysDensityOnGridCalculator =
            std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED> >(
                sysBasisFunctionOnGridController, sys->getSettings().grid.blockAveThreshold);
        auto rsysDensityMatrixController = std::make_shared<DensityMatrixController<Options::SCF_MODES::RESTRICTED> > (
            sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>(),
            sys->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>());
        auto rsysDensityOnGridController =
            std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED> >(
                rsysDensityOnGridCalculator,
                rsysDensityMatrixController);
        (*totalDensity).alpha += 0.5 * rsysDensityOnGridController->getDensityOnGrid();
        (*totalDensity).beta += 0.5 * rsysDensityOnGridController->getDensityOnGrid();
        if (_naddGGA) {
          (*_totalDensityGradient).x.alpha += 0.5 * rsysDensityOnGridController->getDensityGradientOnGrid().x;
          (*_totalDensityGradient).x.beta += 0.5 * rsysDensityOnGridController->getDensityGradientOnGrid().x;
          (*_totalDensityGradient).y.alpha += 0.5 * rsysDensityOnGridController->getDensityGradientOnGrid().y;
          (*_totalDensityGradient).y.beta += 0.5 * rsysDensityOnGridController->getDensityGradientOnGrid().y;
          (*_totalDensityGradient).z.alpha += 0.5 * rsysDensityOnGridController->getDensityGradientOnGrid().z;
          (*_totalDensityGradient).z.beta += 0.5 * rsysDensityOnGridController->getDensityGradientOnGrid().z;
        }
      } else {
        auto sysDensityOnGridCalculator =
            std::make_shared<DensityOnGridCalculator<U> >(
                sysBasisFunctionOnGridController, sys->getSettings().grid.blockAveThreshold);
        auto sysDensityMatrixController = std::make_shared<DensityMatrixController<U> > (
            sys->getActiveOrbitalController<U>(),
            sys->getNOccupiedOrbitals<U>());
        auto sysDensityOnGridController =
            std::make_shared<DensityMatrixDensityOnGridController<U> >(
                sysDensityOnGridCalculator,
                sysDensityMatrixController);
        (*totalDensity) += sysDensityOnGridController->getDensityOnGrid();
        if (_naddGGA) (*_totalDensityGradient) += sysDensityOnGridController->getDensityGradientOnGrid();
      }
    }

    //build DensityOnGridController for total density
    std::shared_ptr<ExternalDensityOnGridController<U> > totalDensityOnGridController;
    if (!_naddGGA) {
      totalDensityOnGridController = std::make_shared<ExternalDensityOnGridController<U> > (totalDensity);
    } else {
      auto totDens = makeGradientPtr<DensityOnGrid<U> >(_gridController);
      (*totDens) = (*_totalDensityGradient);
      totalDensityOnGridController = std::make_shared<ExternalDensityOnGridController<U> > (totalDensity,totDens);
    }

    auto activeNaddKinFuncData = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        FunctionalClassResolver::resolveFunctional(_naddKinFunc),
        activeDensityOnGridController,
        2);
    auto activeNaddXCFuncData = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        FunctionalClassResolver::resolveFunctional(_naddXCFunc),
        activeDensityOnGridController,
        2);
    auto totalNaddKinFuncData = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        FunctionalClassResolver::resolveFunctional(_naddKinFunc),
        totalDensityOnGridController,
        2);
    auto totalNaddXCFuncData = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        FunctionalClassResolver::resolveFunctional(_naddXCFunc),
        totalDensityOnGridController,
        2);

      (*_d2FdRho2).aa += (*totalNaddKinFuncData.d2FdRho2).aa;
      (*_d2FdRho2).ab += (*totalNaddKinFuncData.d2FdRho2).ab;
      (*_d2FdRho2).bb += (*totalNaddKinFuncData.d2FdRho2).bb;
      (*_d2FdRho2).aa -= (*activeNaddKinFuncData.d2FdRho2).aa;
      (*_d2FdRho2).ab -= (*activeNaddKinFuncData.d2FdRho2).ab;
      (*_d2FdRho2).bb -= (*activeNaddKinFuncData.d2FdRho2).bb;
      if (_naddGGA) {
        (*_dFdSigma).aa -= (*activeNaddKinFuncData.dFdSigma).aa;
        (*_dFdSigma).ab -= (*activeNaddKinFuncData.dFdSigma).ab;
        (*_dFdSigma).bb -= (*activeNaddKinFuncData.dFdSigma).bb;
        (*_d2FdSigma2).aaaa -= (*activeNaddKinFuncData.d2FdSigma2).aaaa;
        (*_d2FdSigma2).aaab -= (*activeNaddKinFuncData.d2FdSigma2).aaab;
        (*_d2FdSigma2).aabb -= (*activeNaddKinFuncData.d2FdSigma2).aabb;
        (*_d2FdSigma2).abab -= (*activeNaddKinFuncData.d2FdSigma2).abab;
        (*_d2FdSigma2).abbb -= (*activeNaddKinFuncData.d2FdSigma2).abbb;
        (*_d2FdSigma2).bbbb -= (*activeNaddKinFuncData.d2FdSigma2).bbbb;
        (*_d2FdRhodSigma).aaa -= (*activeNaddKinFuncData.d2FdRhodSigma).aaa;
        (*_d2FdRhodSigma).aab -= (*activeNaddKinFuncData.d2FdRhodSigma).aab;
        (*_d2FdRhodSigma).abb -= (*activeNaddKinFuncData.d2FdRhodSigma).abb;
        (*_d2FdRhodSigma).baa -= (*activeNaddKinFuncData.d2FdRhodSigma).baa;
        (*_d2FdRhodSigma).bab -= (*activeNaddKinFuncData.d2FdRhodSigma).bab;
        (*_d2FdRhodSigma).bbb -= (*activeNaddKinFuncData.d2FdRhodSigma).bbb;
        (*_dFdSigma).aa += (*totalNaddKinFuncData.dFdSigma).aa;
        (*_dFdSigma).ab += (*totalNaddKinFuncData.dFdSigma).ab;
        (*_dFdSigma).bb += (*totalNaddKinFuncData.dFdSigma).bb;
        (*_totalNaddD2FdSigma2).aaaa = (*totalNaddKinFuncData.d2FdSigma2).aaaa;
        (*_totalNaddD2FdSigma2).aaab = (*totalNaddKinFuncData.d2FdSigma2).aaab;
        (*_totalNaddD2FdSigma2).aabb = (*totalNaddKinFuncData.d2FdSigma2).aabb;
        (*_totalNaddD2FdSigma2).abab = (*totalNaddKinFuncData.d2FdSigma2).abab;
        (*_totalNaddD2FdSigma2).abbb = (*totalNaddKinFuncData.d2FdSigma2).abbb;
        (*_totalNaddD2FdSigma2).bbbb = (*totalNaddKinFuncData.d2FdSigma2).bbbb;
        (*_totalNaddD2FdRhodSigma).aaa = (*totalNaddKinFuncData.d2FdRhodSigma).aaa;
        (*_totalNaddD2FdRhodSigma).aab = (*totalNaddKinFuncData.d2FdRhodSigma).aab;
        (*_totalNaddD2FdRhodSigma).abb = (*totalNaddKinFuncData.d2FdRhodSigma).abb;
        (*_totalNaddD2FdRhodSigma).baa = (*totalNaddKinFuncData.d2FdRhodSigma).baa;
        (*_totalNaddD2FdRhodSigma).bab = (*totalNaddKinFuncData.d2FdRhodSigma).bab;
        (*_totalNaddD2FdRhodSigma).bbb = (*totalNaddKinFuncData.d2FdRhodSigma).bbb;
      }

      (*_d2FdRho2).aa += (*totalNaddXCFuncData.d2FdRho2).aa;
      (*_d2FdRho2).ab += (*totalNaddXCFuncData.d2FdRho2).ab;
      (*_d2FdRho2).bb += (*totalNaddXCFuncData.d2FdRho2).bb;
      (*_d2FdRho2).aa -= (*activeNaddXCFuncData.d2FdRho2).aa;
      (*_d2FdRho2).ab -= (*activeNaddXCFuncData.d2FdRho2).ab;
      (*_d2FdRho2).bb -= (*activeNaddXCFuncData.d2FdRho2).bb;
      if (_naddGGA) {
        (*_dFdSigma).aa -= (*activeNaddXCFuncData.dFdSigma).aa;
        (*_dFdSigma).ab -= (*activeNaddXCFuncData.dFdSigma).ab;
        (*_dFdSigma).bb -= (*activeNaddXCFuncData.dFdSigma).bb;
        (*_d2FdSigma2).aaaa -= (*activeNaddXCFuncData.d2FdSigma2).aaaa;
        (*_d2FdSigma2).aaab -= (*activeNaddXCFuncData.d2FdSigma2).aaab;
        (*_d2FdSigma2).aabb -= (*activeNaddXCFuncData.d2FdSigma2).aabb;
        (*_d2FdSigma2).abab -= (*activeNaddXCFuncData.d2FdSigma2).abab;
        (*_d2FdSigma2).abbb -= (*activeNaddXCFuncData.d2FdSigma2).abbb;
        (*_d2FdSigma2).bbbb -= (*activeNaddXCFuncData.d2FdSigma2).bbbb;
        (*_d2FdRhodSigma).aaa -= (*activeNaddXCFuncData.d2FdRhodSigma).aaa;
        (*_d2FdRhodSigma).aab -= (*activeNaddXCFuncData.d2FdRhodSigma).aab;
        (*_d2FdRhodSigma).abb -= (*activeNaddXCFuncData.d2FdRhodSigma).abb;
        (*_d2FdRhodSigma).baa -= (*activeNaddXCFuncData.d2FdRhodSigma).baa;
        (*_d2FdRhodSigma).bab -= (*activeNaddXCFuncData.d2FdRhodSigma).bab;
        (*_d2FdRhodSigma).bbb -= (*activeNaddXCFuncData.d2FdRhodSigma).bbb;
        (*_dFdSigma).aa += (*totalNaddXCFuncData.dFdSigma).aa;
        (*_dFdSigma).ab += (*totalNaddXCFuncData.dFdSigma).ab;
        (*_dFdSigma).bb += (*totalNaddXCFuncData.dFdSigma).bb;
        (*_totalNaddD2FdSigma2).aaaa += (*totalNaddXCFuncData.d2FdSigma2).aaaa;
        (*_totalNaddD2FdSigma2).aaab += (*totalNaddXCFuncData.d2FdSigma2).aaab;
        (*_totalNaddD2FdSigma2).aabb += (*totalNaddXCFuncData.d2FdSigma2).aabb;
        (*_totalNaddD2FdSigma2).abab += (*totalNaddXCFuncData.d2FdSigma2).abab;
        (*_totalNaddD2FdSigma2).abbb += (*totalNaddXCFuncData.d2FdSigma2).abbb;
        (*_totalNaddD2FdSigma2).bbbb += (*totalNaddXCFuncData.d2FdSigma2).bbbb;
        (*_totalNaddD2FdRhodSigma).aaa += (*totalNaddXCFuncData.d2FdRhodSigma).aaa;
        (*_totalNaddD2FdRhodSigma).aab += (*totalNaddXCFuncData.d2FdRhodSigma).aab;
        (*_totalNaddD2FdRhodSigma).abb += (*totalNaddXCFuncData.d2FdRhodSigma).abb;
        (*_totalNaddD2FdRhodSigma).baa += (*totalNaddXCFuncData.d2FdRhodSigma).baa;
        (*_totalNaddD2FdRhodSigma).bab += (*totalNaddXCFuncData.d2FdRhodSigma).bab;
        (*_totalNaddD2FdRhodSigma).bbb += (*totalNaddXCFuncData.d2FdRhodSigma).bbb;
      }
  }
  //Screen Kernel
  double screeningThreshold = 1.0e-8;
  for (unsigned int iPoint = 0; iPoint < _gridController->getNGridPoints(); ++iPoint) {
    if ((*activeDensity).alpha[iPoint] < screeningThreshold) {
      (*_d2FdRho2).aa[iPoint] = 0.0;
      (*_d2FdRho2).ab[iPoint] = 0.0;
      if (_gga) {
        (*_dFdSigma).aa[iPoint] = 0.0;
        (*_dFdSigma).ab[iPoint] = 0.0;
        (*_d2FdSigma2).aaaa[iPoint] = 0.0;
        (*_d2FdSigma2).aaab[iPoint] = 0.0;
        (*_d2FdSigma2).aabb[iPoint] = 0.0;
        (*_d2FdSigma2).abab[iPoint] = 0.0;
        (*_d2FdSigma2).abbb[iPoint] = 0.0;
        (*_d2FdRhodSigma).aaa[iPoint] = 0.0;
        (*_d2FdRhodSigma).aab[iPoint] = 0.0;
        (*_d2FdRhodSigma).abb[iPoint] = 0.0;
        (*_d2FdRhodSigma).baa[iPoint] = 0.0;
        (*_d2FdRhodSigma).bab[iPoint] = 0.0;
      }
      if (_fde && _naddGGA && !_noNaddKernel) {
        (*_totalNaddD2FdSigma2).aaaa[iPoint] = 0.0;
        (*_totalNaddD2FdSigma2).aaab[iPoint] = 0.0;
        (*_totalNaddD2FdSigma2).aabb[iPoint] = 0.0;
        (*_totalNaddD2FdSigma2).abab[iPoint] = 0.0;
        (*_totalNaddD2FdSigma2).abbb[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).aaa[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).aab[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).abb[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).baa[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).bab[iPoint] = 0.0;
      }
    }
    if ((*activeDensity).beta[iPoint] < screeningThreshold) {
      (*_d2FdRho2).ab[iPoint] = 0.0;
      (*_d2FdRho2).bb[iPoint] = 0.0;
      if (_gga) {
        (*_dFdSigma).ab[iPoint] = 0.0;
        (*_dFdSigma).bb[iPoint] = 0.0;
        (*_d2FdSigma2).aabb[iPoint] = 0.0;
        (*_d2FdSigma2).aaab[iPoint] = 0.0;
        (*_d2FdSigma2).abab[iPoint] = 0.0;
        (*_d2FdSigma2).abbb[iPoint] = 0.0;
        (*_d2FdSigma2).bbbb[iPoint] = 0.0;
        (*_d2FdRhodSigma).aab[iPoint] = 0.0;
        (*_d2FdRhodSigma).abb[iPoint] = 0.0;
        (*_d2FdRhodSigma).baa[iPoint] = 0.0;
        (*_d2FdRhodSigma).bab[iPoint] = 0.0;
        (*_d2FdRhodSigma).bbb[iPoint] = 0.0;
      }
      if (_fde && _naddGGA && !_noNaddKernel) {
        (*_totalNaddD2FdSigma2).aaab[iPoint] = 0.0;
        (*_totalNaddD2FdSigma2).aabb[iPoint] = 0.0;
        (*_totalNaddD2FdSigma2).abab[iPoint] = 0.0;
        (*_totalNaddD2FdSigma2).abbb[iPoint] = 0.0;
        (*_totalNaddD2FdSigma2).bbbb[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).aab[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).abb[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).baa[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).bab[iPoint] = 0.0;
        (*_totalNaddD2FdRhodSigma).bbb[iPoint] = 0.0;
      }
    }
  }
  _hasBeenCalculated = true;
}


template<Options::SCF_MODES T>
std::unique_ptr<DensityOnGrid<T> >
Kernel<T>::getDensity(
    std::shared_ptr<SystemController> systemController,
    std::shared_ptr<DensityOnGridController<T> >& densityOnGridController){
  assert(systemController);
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
        systemController->getSettings(), systemController->getBasisController(), _gridController);
  auto densityOnGridCalculator =
        std::make_shared<DensityOnGridCalculator<T> >(
            basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
  auto densityMatrixController = std::make_shared<DensityMatrixController<T> > (
      systemController->getActiveOrbitalController<T>(),
      systemController->getNOccupiedOrbitals<T>());
   densityOnGridController =
          std::make_shared<DensityMatrixDensityOnGridController<T> >(
              densityOnGridCalculator,
              densityMatrixController);
   auto density = std::unique_ptr<DensityOnGrid<T> > (
           new DensityOnGrid<T>(_gridController));
   (*density) = densityOnGridController->getDensityOnGrid();
   return density;
}


template class Kernel<Options::SCF_MODES::RESTRICTED> ;
template class Kernel<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
