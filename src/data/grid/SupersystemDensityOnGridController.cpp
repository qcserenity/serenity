/**
 * @file   SupersystemDensityOnGridController.cpp
 *
 * @date Dezember 30, 2014
 * @author Thomas Dresselhaus
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
#include "data/grid/SupersystemDensityOnGridController.h"
/* Include Std and External Headers */
#include <Eigen/SparseCore>
#include <algorithm>
#include <cassert>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SupersystemDensityOnGridController<SCFMode>::SupersystemDensityOnGridController(
    const std::vector<std::shared_ptr<DensityOnGridController<SCFMode>>>& subsystemDensOnGridControllers)
  : DensityOnGridController<SCFMode>(subsystemDensOnGridControllers[0]->getGridController(),
                                     min_element(subsystemDensOnGridControllers.begin(), subsystemDensOnGridControllers.end(),
                                                 [](std::shared_ptr<DensityOnGridController<SCFMode>> i,
                                                    std::shared_ptr<DensityOnGridController<SCFMode>> j) {
                                                   return (i->getHighestDerivative() < j->getHighestDerivative());
                                                 })
                                         ->get()
                                         ->getHighestDerivative()),
    _subsystemDensOnGridControllers(subsystemDensOnGridControllers),
    _upToDate(false) {
  this->_densityOnGrid.reset(new DensityOnGrid<SCFMode>(this->getGridController()));
  if (this->_highestDerivative >= 1)
    this->_densityGradientOnGrid = makeGradientPtr<DensityOnGrid<SCFMode>>(this->getGridController());
  if (this->_highestDerivative >= 2)
    this->_densityHessianOnGrid = makeHessianPtr<DensityOnGrid<SCFMode>>(this->getGridController());
  for (const auto& subsystemController : subsystemDensOnGridControllers) {
    subsystemController->addSensitiveObject(this->ObjectSensitiveClass<DensityOnGrid<SCFMode>>::_self);
    assert(isDefinedOnSameGrid(*subsystemController, *this));
  }
}

template<Options::SCF_MODES SCFMode>
SupersystemDensityOnGridController<SCFMode>::~SupersystemDensityOnGridController() = default;

template<Options::SCF_MODES SCFMode>
const DensityOnGrid<SCFMode>& SupersystemDensityOnGridController<SCFMode>::getDensityOnGrid() {
  if (!_upToDate)
    updateData();
  return *this->_densityOnGrid;
}

template<Options::SCF_MODES SCFMode>
const Gradient<DensityOnGrid<SCFMode>>& SupersystemDensityOnGridController<SCFMode>::getDensityGradientOnGrid() {
  if (!_upToDate)
    updateData();
  return *this->_densityGradientOnGrid;
}

template<Options::SCF_MODES SCFMode>
const Hessian<DensityOnGrid<SCFMode>>& SupersystemDensityOnGridController<SCFMode>::getDensityHessianOnGrid() {
  if (!_upToDate)
    updateData();
  return *this->_densityHessianOnGrid;
}

template<Options::SCF_MODES SCFMode>
void SupersystemDensityOnGridController<SCFMode>::notify() {
  _upToDate = false;
  this->notifyObjects();
}

template<Options::SCF_MODES SCFMode>
void SupersystemDensityOnGridController<SCFMode>::setHighestDerivative(unsigned int newHighestDerivative) {
  for (auto& subsystemController : _subsystemDensOnGridControllers) {
    if (subsystemController->getHighestDerivative() < newHighestDerivative)
      subsystemController->setHighestDerivative(newHighestDerivative);
  }
  this->_highestDerivative = newHighestDerivative;
  _upToDate = false;
  this->notifyObjects();
}

template<>
void SupersystemDensityOnGridController<RESTRICTED>::updateData() {
  // Re-init data
  /*
   * MB: I have changed the following part of the code quite drastically.
   *     Its tasks are: 1. Calculate the density, gradient etc. for each subsystem
   *                       on the grid.
   *                    2. Add them up.
   *     For subsystem grids and a very large number of subsystems (~1000) step
   *     two became very expensive, since the addition was performed for every
   *     subsystem regardless if it actually had any values different from zero
   *     on the grid. Now, only grid blocks with significant entries are treated
   *     in the addition. The prescreening is performed during the density evaluation
   *     based on the basis function values.
   *     Prescreening information is saved and handled via sparse maps similarly
   *     to the grid-based basis function evaluation.
   */
  std::vector<std::shared_ptr<DensityOnGridController<RESTRICTED>>> screenedSubDensOnGridCont;
  for (const auto& subsystemController : _subsystemDensOnGridControllers) {
    (void)subsystemController->getDensityOnGrid();
    if (this->_highestDerivative >= 1)
      (void)subsystemController->getDensityGradientOnGrid();
    if (this->_highestDerivative >= 2)
      (void)subsystemController->getDensityHessianOnGrid();
    if (subsystemController->getNonNegligibleBlocks().nonZeros() > 0)
      screenedSubDensOnGridCont.push_back(subsystemController);
  }
  // make sure that there is at least one element in the vector in order to avoid 0-dimensional Eigen objects.
  if (screenedSubDensOnGridCont.size() == 0)
    screenedSubDensOnGridCont.push_back(_subsystemDensOnGridControllers[0]);
  const unsigned int nBlocks = screenedSubDensOnGridCont[0]->getNonNegligibleBlocks().size();
  const unsigned int maxBlockSize = screenedSubDensOnGridCont[0]->getMaxBlockSize();
  // Construct the sparse map subsystems-->blocks
  Eigen::SparseMatrix<int> sparseMapSubToBlock(nBlocks, screenedSubDensOnGridCont.size());
  for (unsigned int iSub = 0; iSub < screenedSubDensOnGridCont.size(); ++iSub) {
    const Eigen::SparseVector<int> tmp = screenedSubDensOnGridCont[iSub]->getNonNegligibleBlocks();
    sparseMapSubToBlock.col(iSub) += screenedSubDensOnGridCont[iSub]->getNonNegligibleBlocks();
  }
  // Invert it: Map blocks-->subsystems
  const Eigen::SparseMatrix<int> sparseMapBlockToSub = sparseMapSubToBlock.transpose().eval();
  // Set zero.
  this->_densityOnGrid->setZero();
  if (this->_highestDerivative >= 1) {
    this->_densityGradientOnGrid->x.setZero();
    this->_densityGradientOnGrid->y.setZero();
    this->_densityGradientOnGrid->z.setZero();
  }
  if (this->_highestDerivative >= 2) {
    this->_densityHessianOnGrid->xx.setZero();
    this->_densityHessianOnGrid->xy.setZero();
    this->_densityHessianOnGrid->xz.setZero();
    this->_densityHessianOnGrid->yy.setZero();
    this->_densityHessianOnGrid->yz.setZero();
    this->_densityHessianOnGrid->zz.setZero();
  }

  /*
   * - Loop blocks.
   * - Loop subsystems with significant values on the given block based on the
   *   sparse maps above.
   */
  unsigned int nPoints = this->_densityOnGrid->size();
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned int n = maxBlockSize;
    const unsigned int start = iBlock * n;
    if (iBlock == nBlocks - 1) {
      n = nPoints % maxBlockSize;
      if (n == 0)
        n = maxBlockSize;
    }
    for (Eigen::SparseMatrix<int>::InnerIterator itSub(sparseMapBlockToSub, iBlock); itSub; ++itSub) {
      unsigned int iSub = itSub.row();
      const auto& subsystemController = screenedSubDensOnGridCont[iSub];
      this->_densityOnGrid->segment(start, n) += subsystemController->getDensityOnGrid().segment(start, n);
      if (this->_highestDerivative >= 1) {
        this->_densityGradientOnGrid->x.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().x.segment(start, n);
        this->_densityGradientOnGrid->y.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().y.segment(start, n);
        this->_densityGradientOnGrid->z.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().z.segment(start, n);
      }
      if (this->_highestDerivative >= 2) {
        this->_densityHessianOnGrid->xx.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xx.segment(start, n);
        this->_densityHessianOnGrid->xy.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xy.segment(start, n);
        this->_densityHessianOnGrid->xz.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xz.segment(start, n);
        this->_densityHessianOnGrid->yy.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().yy.segment(start, n);
        this->_densityHessianOnGrid->yz.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().yz.segment(start, n);
        this->_densityHessianOnGrid->zz.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().zz.segment(start, n);
      }
    }
  } // for iBlock
}

template<>
void SupersystemDensityOnGridController<UNRESTRICTED>::updateData() {
  // Re-init data
  std::vector<std::shared_ptr<DensityOnGridController<UNRESTRICTED>>> screenedSubDensOnGridCont;
  for (const auto& subsystemController : _subsystemDensOnGridControllers) {
    (void)subsystemController->getDensityOnGrid();
    if (this->_highestDerivative >= 1)
      (void)subsystemController->getDensityGradientOnGrid();
    if (this->_highestDerivative >= 2)
      (void)subsystemController->getDensityHessianOnGrid();
    if (subsystemController->getNonNegligibleBlocks().size() > 0)
      screenedSubDensOnGridCont.push_back(subsystemController);
  }
  if (screenedSubDensOnGridCont.size() == 0)
    screenedSubDensOnGridCont.push_back(_subsystemDensOnGridControllers[0]);
  const unsigned int nBlocks = screenedSubDensOnGridCont[0]->getNonNegligibleBlocks().size();
  const unsigned int maxBlockSize = screenedSubDensOnGridCont[0]->getMaxBlockSize();
  Eigen::SparseMatrix<int> sparseMapSubToBlock(nBlocks, screenedSubDensOnGridCont.size());
  for (unsigned int iSub = 0; iSub < screenedSubDensOnGridCont.size(); ++iSub) {
    const Eigen::SparseVector<int> tmp = screenedSubDensOnGridCont[iSub]->getNonNegligibleBlocks();
    sparseMapSubToBlock.col(iSub) += tmp.eval();
  }
  const Eigen::SparseMatrix<int> sparseMapBlockToSub = sparseMapSubToBlock.transpose().eval();
  this->_densityOnGrid->alpha.setZero();
  this->_densityOnGrid->beta.setZero();
  if (this->_highestDerivative >= 1) {
    this->_densityGradientOnGrid->x.alpha.setZero();
    this->_densityGradientOnGrid->y.alpha.setZero();
    this->_densityGradientOnGrid->z.alpha.setZero();

    this->_densityGradientOnGrid->x.beta.setZero();
    this->_densityGradientOnGrid->y.beta.setZero();
    this->_densityGradientOnGrid->z.beta.setZero();
  }
  if (this->_highestDerivative >= 2) {
    this->_densityHessianOnGrid->xx.alpha.setZero();
    this->_densityHessianOnGrid->xy.alpha.setZero();
    this->_densityHessianOnGrid->xz.alpha.setZero();
    this->_densityHessianOnGrid->yy.alpha.setZero();
    this->_densityHessianOnGrid->yz.alpha.setZero();
    this->_densityHessianOnGrid->zz.alpha.setZero();

    this->_densityHessianOnGrid->xx.beta.setZero();
    this->_densityHessianOnGrid->xy.beta.setZero();
    this->_densityHessianOnGrid->xz.beta.setZero();
    this->_densityHessianOnGrid->yy.beta.setZero();
    this->_densityHessianOnGrid->yz.beta.setZero();
    this->_densityHessianOnGrid->zz.beta.setZero();
  }

  unsigned int nPoints = this->_densityOnGrid->alpha.size();
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned int n = maxBlockSize;
    const unsigned int start = iBlock * n;
    if (iBlock == nBlocks - 1) {
      n = nPoints % maxBlockSize;
      if (n == 0)
        n = maxBlockSize;
    }
    for (Eigen::SparseMatrix<int>::InnerIterator itSub(sparseMapBlockToSub, iBlock); itSub; ++itSub) {
      unsigned int iSub = itSub.row();
      const auto& subsystemController = screenedSubDensOnGridCont[iSub];
      (*this->_densityOnGrid).alpha.segment(start, n) += subsystemController->getDensityOnGrid().alpha.segment(start, n);
      (*this->_densityOnGrid).beta.segment(start, n) += subsystemController->getDensityOnGrid().beta.segment(start, n);
      if (this->_highestDerivative >= 1) {
        this->_densityGradientOnGrid->x.alpha.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().x.alpha.segment(start, n);
        this->_densityGradientOnGrid->y.alpha.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().y.alpha.segment(start, n);
        this->_densityGradientOnGrid->z.alpha.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().z.alpha.segment(start, n);
        this->_densityGradientOnGrid->x.beta.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().x.beta.segment(start, n);
        this->_densityGradientOnGrid->y.beta.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().y.beta.segment(start, n);
        this->_densityGradientOnGrid->z.beta.segment(start, n) +=
            subsystemController->getDensityGradientOnGrid().z.beta.segment(start, n);
      }
      if (this->_highestDerivative >= 2) {
        this->_densityHessianOnGrid->xx.alpha.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xx.alpha.segment(start, n);
        this->_densityHessianOnGrid->xy.alpha.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xy.alpha.segment(start, n);
        this->_densityHessianOnGrid->xz.alpha.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xz.alpha.segment(start, n);
        this->_densityHessianOnGrid->yy.alpha.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().yy.alpha.segment(start, n);
        this->_densityHessianOnGrid->yz.alpha.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().yz.alpha.segment(start, n);
        this->_densityHessianOnGrid->zz.alpha.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().zz.alpha.segment(start, n);
        this->_densityHessianOnGrid->xx.beta.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xx.beta.segment(start, n);
        this->_densityHessianOnGrid->xy.beta.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xy.beta.segment(start, n);
        this->_densityHessianOnGrid->xz.beta.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().xz.beta.segment(start, n);
        this->_densityHessianOnGrid->yy.beta.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().yy.beta.segment(start, n);
        this->_densityHessianOnGrid->yz.beta.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().yz.beta.segment(start, n);
        this->_densityHessianOnGrid->zz.beta.segment(start, n) +=
            subsystemController->getDensityHessianOnGrid().zz.beta.segment(start, n);
      }
    }
  } // for iBlock
}

template class SupersystemDensityOnGridController<Options::SCF_MODES::RESTRICTED>;
template class SupersystemDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
