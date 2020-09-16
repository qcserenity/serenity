/**
 * @file ExternalDensityOnGridController.h
 *
 * @date Dec 30, 2014
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
#ifndef EXTERNALDENSITYONGRIDCONTROLLER_H
#define EXTERNALDENSITYONGRIDCONTROLLER_H
/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGridController.h"
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <functional>

namespace Serenity {
/* Forward Declarations */

/**
 * @class ExternalDensityOnGridController ExternalDensityOnGridController.h
 * @brief Wraps a DensityOnGrid and its derivatives when provided externally or constant.
 *
 * It is an ObjectSensitiveClass to give the possibility to notify it when the externally provided
 * data changes. Because it is completely unknown how data are provided the template argument is
 * void (i.e. may be anything).
 */
template<Options::SCF_MODES T>
class ExternalDensityOnGridController : public DensityOnGridController<T>, public ObjectSensitiveClass<void> {
 public:
  /**
   * @brief Constructor.
   * @param densityOnGrid The density on a grid.
   */
  ExternalDensityOnGridController(std::unique_ptr<DensityOnGrid<T>>& densityOnGrid);

  /**
   * @brief Constructor.
   * @param densityOnGrid The density on a grid.
   * @param densityGradientOnGrid The gradient of the density on a grid.
   */
  ExternalDensityOnGridController(std::unique_ptr<DensityOnGrid<T>>& densityOnGrid,
                                  std::unique_ptr<Gradient<DensityOnGrid<T>>>& densityGradientOnGrid);
  /**
   * @brief Constructor
   * @param densityOnGrid The density on a grid.
   * @param densityGradientOnGrid The gradient of the density on a grid.
   * @param densityHessianOnGrid The hessian of the density on a grid.
   */
  ExternalDensityOnGridController(std::unique_ptr<DensityOnGrid<T>>& densityOnGrid,
                                  std::unique_ptr<Gradient<DensityOnGrid<T>>>& densityGradientOnGrid,
                                  std::unique_ptr<Hessian<DensityOnGrid<T>>>& densityHessianOnGrid);
  /**
   * @brief Constructor.
   * @param densityOnGrid The density on a grid.
   * @param densityGradientOnGrid The gradient of the density on a grid.
   * @param densityHessianOnGrid The hessian of the density on a grid.
   * @param externalUpdateFunction This function is called whenever an object of this class is notified (->notify()).
   */
  ExternalDensityOnGridController(std::unique_ptr<DensityOnGrid<T>>& densityOnGrid,
                                  std::unique_ptr<Gradient<DensityOnGrid<T>>>& densityGradientOnGrid,
                                  std::unique_ptr<Hessian<DensityOnGrid<T>>>& densityHessianOnGrid,
                                  std::function<void(DensityOnGrid<T>& densityOnGrid, Gradient<DensityOnGrid<T>>& densityGradientOnGrid,
                                                     Hessian<DensityOnGrid<T>>& densityHessianOnGrid)>
                                      externalUpdateFunction);
  /**
   * @brief Default destructor.
   */
  virtual ~ExternalDensityOnGridController();
  /**
   * @return The density on a grid.
   */
  virtual const DensityOnGrid<T>& getDensityOnGrid() override final;
  /**
   * @return The density gradient on a gird.
   */
  virtual const Gradient<DensityOnGrid<T>>& getDensityGradientOnGrid() override final;
  /**
   * @return The hessian of the density on a grid.
   */
  virtual const Hessian<DensityOnGrid<T>>& getDensityHessianOnGrid() override final;
  /**
   * @brief Delete old data and update it.
   */
  virtual void notify() override final;
  /**
   * @brief Setter for the highest derivative.
   * @param The order of the highest derivative.
   */
  virtual void setHighestDerivative(unsigned int) override final {
    assert(false);
  }

 private:
  /**
   * @brief The external update function.
   */
  std::function<void(DensityOnGrid<T>& densityOnGrid, Gradient<DensityOnGrid<T>>& densityGradientOnGrid,
                     Hessian<DensityOnGrid<T>>& densityHessianOnGrid)>
      _externalUpdateFunction;
};

} /* namespace Serenity */
#endif /* EXTERNALDENSITYONGRIDCONTROLLER_H */
