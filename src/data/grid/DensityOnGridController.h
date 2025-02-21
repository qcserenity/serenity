/**
 * @file   DensityOnGridController.h
 * @author Thomas Dresselhaus
 *
 * @date   22. November 2014, 10:55
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
#ifndef DENSITYONGRIDCONTROLLER_H
#define DENSITYONGRIDCONTROLLER_H
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/grid/DensityOnGrid.h"
#include "math/Derivatives.h"
#include "misc/SerenityError.h"
#include "notification/NotifyingClass.h"
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <Eigen/SparseCore>
#include <iostream>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward Declarations */
class GridController;
class Grid;
/**
 * @class DensityOnGridController DensityOnGridController.h
 * @brief Interface to manage the density and its derivatives represented on a grid.
 *
 * Caution: Make sure in every implementation that the notifyObjects() method is called whenever
 * the managed data changes!
 */
template<Options::SCF_MODES SCFMode>
class DensityOnGridController : public NotifyingClass<DensityOnGrid<SCFMode>>, public ObjectSensitiveClass<Grid> {
 public:
  /**
   * @param gridController
   * @param highestDerivative see the getter function
   */
  DensityOnGridController(std::shared_ptr<GridController> gridController, const unsigned int highestDerivative);

  virtual ~DensityOnGridController() = default;
  /**
   * @returns the current (lazily evaluated) density on the grid
   */
  virtual const DensityOnGrid<SCFMode>& getDensityOnGrid() = 0;
  /**
   * @returns the current (lazily evaluated) gradient of the density on the grid
   */
  virtual const Gradient<DensityOnGrid<SCFMode>>& getDensityGradientOnGrid() = 0;
  /**
   * @returns the current (lazily evaluated) Hessian of the density on the grid
   */
  virtual const Hessian<DensityOnGrid<SCFMode>>& getDensityHessianOnGrid() = 0;
  /**
   * @returns the number of grid points of the Grid on which the density will be represented.
   */
  unsigned int getNGridPoints() const;
  /*
   * @returns the highest derivative allowed by this controller, i.e. 0 => density only, 1=> density
   *          and gradient, ...
   */
  unsigned int getHighestDerivative() const {
    return _highestDerivative;
  }
  /**
   * @brief sets a new highest derivative.
   *
   * This indicates up to which order derivatives are calculated by this object.
   * TODO this solution does not look elegant.
   */
  virtual void setHighestDerivative(unsigned int newHighestDerivative) = 0;
  /**
   * @returns the Grid(controller) the object is working on
   */
  std::shared_ptr<GridController> getGridController() const {
    return _gridController;
  }

  void notify() override {
    _upToDate = false;
    this->notifyObjects();
  }
  /**
   * @brief Getter for the blocks with entries different from zero.
   *        Note that the grid may become extremely sparsely occupied during
   *        embedding calculations.
   *
   *        This default function should never be used!
   */
  virtual Eigen::SparseVector<int> getNonNegligibleBlocks() {
    throw SerenityError("Overwrite this function in the specific realizations, please!");
    Eigen::VectorXi tmp = Eigen::VectorXi::Constant(1, 1);
    return tmp.sparseView();
  }
  /**
   * @brief Getter for the number of blocks.
   *        This default function should never be used!
   */
  virtual unsigned int getNBlocks() {
    throw SerenityError("Overwrite this function in the specific realizations, please!");
    return 1;
  }
  /**
   * @brief Getter for the maximum block size.
   *        This default function should never be used!
   */
  virtual unsigned int getMaxBlockSize() {
    throw SerenityError("Overwrite this function in the specific realizations, please!");
    return _nGridPoints;
  }

 protected:
  const std::shared_ptr<GridController> _gridController;

  unsigned int _highestDerivative;

  unsigned int _nGridPoints;

  bool _upToDate;

  std::unique_ptr<DensityOnGrid<SCFMode>> _densityOnGrid;
  std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>> _densityGradientOnGrid;
  std::unique_ptr<Hessian<DensityOnGrid<SCFMode>>> _densityHessianOnGrid;
};

} // namespace Serenity
#endif /* DENSITYONGRIDCONTROLLER_H */
