/**
 * @file   SupersystemDensityOnGridController.cpp
 * 
 * @date Dezember 30, 2014
 * @author Thomas Dresselhaus
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
#include "data/grid/SupersystemDensityOnGridController.h"
/* Include Std and External Headers */
#include <algorithm>
#include <cassert>


namespace Serenity {

template<Options::SCF_MODES T>
SupersystemDensityOnGridController<T>::SupersystemDensityOnGridController(
      const std::vector<std::shared_ptr<DensityOnGridController<T> > >& subsystemDensOnGridControllers) :
    DensityOnGridController<T>(
      subsystemDensOnGridControllers[0]->getGridController(),
      min_element(
        subsystemDensOnGridControllers.begin(),
        subsystemDensOnGridControllers.end(),
        [](std::shared_ptr<DensityOnGridController<T> > i,
            std::shared_ptr<DensityOnGridController<T> > j) {
          return (i->getHighestDerivative() < j->getHighestDerivative());
        })->get()->getHighestDerivative()),
    _subsystemDensOnGridControllers(subsystemDensOnGridControllers),
    _upToDate(false) {
  this->_densityOnGrid.reset(new DensityOnGrid<T>(this->getGridController()));
  if (this->_highestDerivative >= 1)
    this->_densityGradientOnGrid = makeGradientPtr<DensityOnGrid<T> >(this->getGridController());
  if (this->_highestDerivative >= 2)
    this->_densityHessianOnGrid = makeHessianPtr<DensityOnGrid<T> >(this->getGridController());
  for (const auto& subsystemController : subsystemDensOnGridControllers) {
    subsystemController->addSensitiveObject(this->ObjectSensitiveClass<DensityOnGrid<T> >::_self);
    assert(isDefinedOnSameGrid(*subsystemController, *this));
  }
}

template<Options::SCF_MODES T>
SupersystemDensityOnGridController<T>::~SupersystemDensityOnGridController() {
}

template<Options::SCF_MODES T>
const DensityOnGrid<T>& SupersystemDensityOnGridController<T>::getDensityOnGrid() {
  if (!_upToDate) updateData();
  return *this->_densityOnGrid;
}

template<Options::SCF_MODES T>
const Gradient<DensityOnGrid<T> >&
    SupersystemDensityOnGridController<T>::getDensityGradientOnGrid() {
  if (!_upToDate) updateData();
  return *this->_densityGradientOnGrid;
}

template<Options::SCF_MODES T>
const Hessian<DensityOnGrid<T> >&
    SupersystemDensityOnGridController<T>::getDensityHessianOnGrid() {
  if (!_upToDate) updateData();
  return *this->_densityHessianOnGrid;
}

template<Options::SCF_MODES T>
void SupersystemDensityOnGridController<T>::notify() {
  _upToDate = false;
  this->notifyObjects();
}

template<Options::SCF_MODES T>
void SupersystemDensityOnGridController<T>::setHighestDerivative(
      unsigned int newHighestDerivative) {
  for (auto& subsystemController : _subsystemDensOnGridControllers) {
    if (subsystemController->getHighestDerivative() < newHighestDerivative)
          subsystemController->setHighestDerivative(newHighestDerivative);
  }
  this->_highestDerivative=newHighestDerivative;
  _upToDate = false;
  this->notifyObjects();
}

template<>
void SupersystemDensityOnGridController<RESTRICTED>::updateData() {
  // Re-init data
  for (const auto& subsystemController : _subsystemDensOnGridControllers) {
    (void)subsystemController->getDensityOnGrid();
    if (this->_highestDerivative >= 1) (void)subsystemController->getDensityGradientOnGrid();
    if (this->_highestDerivative >= 2) (void)subsystemController->getDensityHessianOnGrid();
  }

  unsigned int nPoints =  this->_densityOnGrid->size();
  unsigned int nBlocks = omp_get_max_threads();

#pragma omp parallel for schedule (dynamic)
  for (unsigned int iBlock = 0; iBlock<nBlocks;iBlock++){
    unsigned int n = (unsigned int)(nPoints/nBlocks);
    const unsigned int start = iBlock*n;
    if(iBlock==nBlocks-1) n += nPoints%nBlocks;
    // Set first data instead of setting to zero
      (*this->_densityOnGrid).segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityOnGrid().segment(start,n);
      if (this->_highestDerivative >= 1) {
        this->_densityGradientOnGrid->x.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().x.segment(start,n);
        this->_densityGradientOnGrid->y.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().y.segment(start,n);
        this->_densityGradientOnGrid->z.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().z.segment(start,n);
      }
      if (this->_highestDerivative >= 2) {
        this->_densityHessianOnGrid->xx.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xx.segment(start,n);
        this->_densityHessianOnGrid->xy.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xy.segment(start,n);
        this->_densityHessianOnGrid->xz.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xz.segment(start,n);
        this->_densityHessianOnGrid->yy.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().yy.segment(start,n);
        this->_densityHessianOnGrid->yz.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().yz.segment(start,n);
        this->_densityHessianOnGrid->zz.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().zz.segment(start,n);
      }

    // add the rest
    for (unsigned int i=1;i<_subsystemDensOnGridControllers.size();i++) {
      (*this->_densityOnGrid).segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityOnGrid().segment(start,n);
      if (this->_highestDerivative >= 1) {
        this->_densityGradientOnGrid->x.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().x.segment(start,n);
        this->_densityGradientOnGrid->y.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().y.segment(start,n);
        this->_densityGradientOnGrid->z.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().z.segment(start,n);
      }
      if (this->_highestDerivative >= 2) {
        this->_densityHessianOnGrid->xx.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xx.segment(start,n);
        this->_densityHessianOnGrid->xy.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xy.segment(start,n);
        this->_densityHessianOnGrid->xz.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xz.segment(start,n);
        this->_densityHessianOnGrid->yy.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().yy.segment(start,n);
        this->_densityHessianOnGrid->yz.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().yz.segment(start,n);
        this->_densityHessianOnGrid->zz.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().zz.segment(start,n);
      }
    }
  }
}

template<>
void SupersystemDensityOnGridController<UNRESTRICTED>::updateData() {
  // Re-init data
  for (const auto& subsystemController : _subsystemDensOnGridControllers) {
    (void)subsystemController->getDensityOnGrid();
    if (this->_highestDerivative >= 1) (void)subsystemController->getDensityGradientOnGrid();
    if (this->_highestDerivative >= 2) (void)subsystemController->getDensityHessianOnGrid();
  }

  unsigned int nPoints =  this->_densityOnGrid->alpha.size();
  unsigned int nBlocks = omp_get_max_threads();

#pragma omp parallel for schedule (dynamic)
  for (unsigned int iBlock = 0; iBlock<nBlocks;iBlock++){
    unsigned int n = (unsigned int)(nPoints/nBlocks);
    const unsigned int start = iBlock*n;
    if(iBlock==nBlocks-1) n += nPoints%nBlocks;

    // Set first data instead of setting to zero
      (*this->_densityOnGrid).alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityOnGrid().alpha.segment(start,n);
      (*this->_densityOnGrid).beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityOnGrid().beta.segment(start,n);
      if (this->_highestDerivative >= 1) {
        this->_densityGradientOnGrid->x.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().x.alpha.segment(start,n);
        this->_densityGradientOnGrid->y.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().y.alpha.segment(start,n);
        this->_densityGradientOnGrid->z.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().z.alpha.segment(start,n);
        this->_densityGradientOnGrid->x.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().x.beta.segment(start,n);
        this->_densityGradientOnGrid->y.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().y.beta.segment(start,n);
        this->_densityGradientOnGrid->z.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityGradientOnGrid().z.beta.segment(start,n);
      }
      if (this->_highestDerivative >= 2) {
        this->_densityHessianOnGrid->xx.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xx.alpha.segment(start,n);
        this->_densityHessianOnGrid->xy.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xy.alpha.segment(start,n);
        this->_densityHessianOnGrid->xz.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xz.alpha.segment(start,n);
        this->_densityHessianOnGrid->yy.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().yy.alpha.segment(start,n);
        this->_densityHessianOnGrid->yz.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().yz.alpha.segment(start,n);
        this->_densityHessianOnGrid->zz.alpha.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().zz.alpha.segment(start,n);
        this->_densityHessianOnGrid->xx.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xx.beta.segment(start,n);
        this->_densityHessianOnGrid->xy.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xy.beta.segment(start,n);
        this->_densityHessianOnGrid->xz.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().xz.beta.segment(start,n);
        this->_densityHessianOnGrid->yy.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().yy.beta.segment(start,n);
        this->_densityHessianOnGrid->yz.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().yz.beta.segment(start,n);
        this->_densityHessianOnGrid->zz.beta.segment(start,n) = _subsystemDensOnGridControllers[0]->getDensityHessianOnGrid().zz.beta.segment(start,n);
      }

    // add the rest
    for (unsigned int i=1;i<_subsystemDensOnGridControllers.size();i++) {
      (*this->_densityOnGrid).alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityOnGrid().alpha.segment(start,n);
      (*this->_densityOnGrid).beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityOnGrid().beta.segment(start,n);
      if (this->_highestDerivative >= 1) {
        this->_densityGradientOnGrid->x.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().x.alpha.segment(start,n);
        this->_densityGradientOnGrid->y.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().y.alpha.segment(start,n);
        this->_densityGradientOnGrid->z.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().z.alpha.segment(start,n);
        this->_densityGradientOnGrid->x.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().x.beta.segment(start,n);
        this->_densityGradientOnGrid->y.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().y.beta.segment(start,n);
        this->_densityGradientOnGrid->z.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityGradientOnGrid().z.beta.segment(start,n);
      }
      if (this->_highestDerivative >= 2) {
        this->_densityHessianOnGrid->xx.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xx.alpha.segment(start,n);
        this->_densityHessianOnGrid->xy.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xy.alpha.segment(start,n);
        this->_densityHessianOnGrid->xz.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xz.alpha.segment(start,n);
        this->_densityHessianOnGrid->yy.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().yy.alpha.segment(start,n);
        this->_densityHessianOnGrid->yz.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().yz.alpha.segment(start,n);
        this->_densityHessianOnGrid->zz.alpha.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().zz.alpha.segment(start,n);
        this->_densityHessianOnGrid->xx.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xx.beta.segment(start,n);
        this->_densityHessianOnGrid->xy.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xy.beta.segment(start,n);
        this->_densityHessianOnGrid->xz.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().xz.beta.segment(start,n);
        this->_densityHessianOnGrid->yy.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().yy.beta.segment(start,n);
        this->_densityHessianOnGrid->yz.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().yz.beta.segment(start,n);
        this->_densityHessianOnGrid->zz.beta.segment(start,n) += _subsystemDensOnGridControllers[i]->getDensityHessianOnGrid().zz.beta.segment(start,n);

      }
    }
  }
}

template class SupersystemDensityOnGridController<Options::SCF_MODES::RESTRICTED>;
template class SupersystemDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
