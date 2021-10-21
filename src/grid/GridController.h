/**
 * @file   GridController.h
 *
 * @date   Mar 15, 2014
 * @author Jan Unsleber
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
#ifndef GRIDCONTROLLER_H_
#define GRIDCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "grid/Grid.h"
#include "notification/NotifyingClass.h"
/* Include Std and External Headers */
#include <memory>
#include <utility>
#include <vector>

namespace Serenity {
/* Forward declarations */
class Point;
/**
 * @class GridController GridController.h
 * @brief Manages a complete integration grid, e.g. of a System
 *
 * This base class can be used to work with pre-existing grids. Most often probably the derived
 * class AtomCenteredGridController is used in practise.
 *
 * If you want to inherit from this, override the produceGrid function. Also make sure you use
 * the notification system properly.
 */
class GridController : public NotifyingClass<Grid> {
 public:
  /**
   * @brief Directly controls a supplied grid instead of producing one itself.
   * @param grid
   */
  GridController(std::unique_ptr<Grid> grid);

 public:
  virtual ~GridController() = default;
  /**
   * @brief   gathers grid points from all atoms of the underlying system
   *
   * @returns all grid points for the underlying system
   */
  virtual const Eigen::Matrix3Xd& getGridPoints();
  /**
   * @brief   gathers all the weights of the grid points from all atoms of the underlying system
   *
   * @returns the weights of all grid points of the underlying system
   */
  virtual const Eigen::VectorXd& getWeights();
  /**
   * @returns the number of grid points of the underlying grid
   *          (if not yet present a grid is created upon call).
   */
  virtual unsigned int getNGridPoints();

 protected:
  /**
   * @brief Construct the object without a grid. Only makes sense to use in derived classes.
   * This constructor makes a lazy grid construction easily possible.
   */
  GridController() : _grid(nullptr) {
  }
  /**
   * @brief Creates a new integration grid.
   * This base class cannot create an integration grid by itself, so this function must be
   * overridden in derived classes. However, if a pre-defined grid is already available, a
   * GridController may already be constructed. In that case, this method must not be called.
   */
  virtual void produceGrid() {
    assert(false);
  }
  /**
   * @brief the currently controlled grid
   */
  std::unique_ptr<Grid> _grid;
};

} // namespace Serenity
#endif /* GRIDCONTROLLER_H_ */
