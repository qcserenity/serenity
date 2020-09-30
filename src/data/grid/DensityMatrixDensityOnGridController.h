/**
 * @file   DensityMatrixDensityOnGridController.h
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   30. Dezember 2014, 16:05
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
#ifndef DENSITYMATRIXDENSITYONGRIDCONTROLLER_H
#define DENSITYMATRIXDENSITYONGRIDCONTROLLER_H
/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGridController.h"
#include "data/matrices/DensityMatrix.h"
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <Eigen/SparseCore>

namespace Serenity {
/* Forward Declarations */
template<Options::SCF_MODES T>
class DensityMatrixController;
template<Options::SCF_MODES T>
class DensityOnGridCalculator;
/**
 * @class DensityMatrixDensityOnGridController DensityMatrixDensityOnGridController.h
 * @brief Calculates the DensityOnGrid based on a density matrix. This is the most common way.
 *
 * Wraps a DensityOnGridCalculator, makes sure densities are up-to-date, avoids recalculation.
 */
template<Options::SCF_MODES T>
class DensityMatrixDensityOnGridController : public ObjectSensitiveClass<DensityMatrix<T>>,
                                             public DensityOnGridController<T> {
 public:
  /**
   * @brief Constructor
   * @param densOnGridCalculator A DensityOnGridCalculator to calculate the needed data. The highest available
   * derivative is determined by this object.
   * @param densityMatrix The density to be represented on a grid is given here in matrix form.
   */
  DensityMatrixDensityOnGridController(std::shared_ptr<DensityOnGridCalculator<T>> densOnGridCalculator,
                                       const std::shared_ptr<DensityMatrixController<T>> densityMatrixController);
  /**
   * @brief Constructor
   * @param densOnGridCalculator see above
   * @param densityMatrix see above
   * @param highestDerivative limits the actually calculated derivative. Must not be larger than
   *                          allowed by the densOnGridCalculator.
   */
  DensityMatrixDensityOnGridController(std::shared_ptr<DensityOnGridCalculator<T>> densOnGridCalculator,
                                       const std::shared_ptr<DensityMatrixController<T>> densityMatrixController,
                                       const unsigned int highestDerivative);

  virtual ~DensityMatrixDensityOnGridController() = default;

  /**
   * @brief Getter for the density on a grid. Uses lazy evaluation.
   * @return The density on a grid.
   */
  virtual const DensityOnGrid<T>& getDensityOnGrid() override final;

  /**
   * @brief Getter for the density gradient on a grid. Uses lazy evaluation.
   * @return The density gradient on a grid.
   */
  virtual const Gradient<DensityOnGrid<T>>& getDensityGradientOnGrid() override final;

  /**
   * @brief Getter for the density Hessian on a grid. Uses lazy evaluation.
   * @return The density Hessian on a grid.
   */
  virtual const Hessian<DensityOnGrid<T>>& getDensityHessianOnGrid() override final;

  /**
   * @brief Sets a custom density on a grid.
   */
  void setDensityOnGrid(std::unique_ptr<DensityOnGrid<T>> densityOnGrid);

  /**
   * @brief Sets a custom density gradient on a grid.
   */
  void setDensityGradientOnGrid(std::unique_ptr<Gradient<DensityOnGrid<T>>> densityGradientOnGrid);

  /**
   * @brief Sets a custom density Hessian on a grid.
   */
  void setDensityHessianOnGrid(std::unique_ptr<Hessian<DensityOnGrid<T>>> densityHessianOnGrid);

  /**
   * @brief Is triggered by the notification system to delete outdated data and induce a recalculation.
   */
  virtual void notify() override final;

  /**
   * @brief Sets the highest derivative of the density to be calculated on the grid.
   */
  virtual void setHighestDerivative(unsigned int newHighestDerivative) override final;

  /**
   * @brief Getter for the blocks with entries different from zero.
   *        Note that the grid may become extremely sparsely occupied during
   *        embedding calculations.
   */
  virtual Eigen::SparseVector<int> getNonNegligibleBlocks() override {
    return _densOnGridCalculator->getNonNegligibleBlocks();
  }
  /**
   * @brief Getter for the number of blocks.
   */
  virtual unsigned int getNBlocks() override {
    return _densOnGridCalculator->getBasisFunctionOnGridController()->getNBlocks();
  }
  /**
   * @brief Getter for the maximum block size.
   */
  virtual unsigned int getMaxBlockSize() override {
    return _densOnGridCalculator->getBasisFunctionOnGridController()->getMaxBlockSize();
  }

 private:
  void calculateData();

  const std::shared_ptr<DensityOnGridCalculator<T>> _densOnGridCalculator;
  std::shared_ptr<DensityMatrixController<T>> _densityMatrixController;
};

} /* namespace Serenity */
#endif /* DENSITYMATRIXDENSITYONGRIDCONTROLLER_H */
