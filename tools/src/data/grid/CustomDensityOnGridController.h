/**
 * @file   CustomDensityOnGridController.h
 * @author Jan Unsleber
 *
 * @date   Oct 24, 2017
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

#ifndef DATA_GRID_CUSTOMDENSITYONGRIDCONTROLLER_H_
#define DATA_GRID_CUSTOMDENSITYONGRIDCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGridController.h"
#include "math/Derivatives.h"

namespace Serenity {

/**
 * @class CustomDensityOnGridController CustomDensityOnGridController.h
 *
 * @brief A small class to implement a DensityOnGridController for self generated densities on a grid.
 */
template<Options::SCF_MODES SCFMode>
class CustomDensityOnGridController : public DensityOnGridController<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param dens The density as std::unique_ptr (will be moved).
   * @param grad The density gradient as std::unique_ptr (will be moved).
   */
  CustomDensityOnGridController(std::unique_ptr<DensityOnGrid<SCFMode>> dens,
                                std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>> grad)
    : DensityOnGridController<SCFMode>(dens->getGridController(), 1) {
    this->_densityOnGrid = std::move(dens);
    this->_densityGradientOnGrid = std::move(grad);
    this->_densityHessianOnGrid.reset(nullptr);
    assert(this->_densityOnGrid->getGridController() == this->_gridController);
    if (this->_densityGradientOnGrid)
      assert(this->_densityGradientOnGrid->x.getGridController() == this->_gridController);
  }
  /// @brief Default destructor.
  virtual ~CustomDensityOnGridController() = default;

  /// @brief @see  DensityOnGridController
  virtual const DensityOnGrid<SCFMode>& getDensityOnGrid() override final {
    return *(this->_densityOnGrid);
  }
  /// @brief @see  DensityOnGridController
  virtual const Gradient<DensityOnGrid<SCFMode>>& getDensityGradientOnGrid() override final {
    return *(this->_densityGradientOnGrid);
  }
  /// @brief @see  DensityOnGridController
  virtual const Hessian<DensityOnGrid<SCFMode>>& getDensityHessianOnGrid() override final {
    assert(false && "ERROR: The CustomDensityOnGridController wasn't build for other derivatives");
    return *(this->_densityHessianOnGrid);
  }
  /// @brief @see  DensityOnGridController
  virtual void setHighestDerivative(unsigned int newHighestDerivative) override final {
    (void)newHighestDerivative;
    assert(false && "ERROR: The CustomDensityOnGridController wasn't build for other derivatives");
  }
};

} /* namespace Serenity */

#endif /* DATA_GRID_CUSTOMDENSITYONGRIDCONTROLLER_H_ */
