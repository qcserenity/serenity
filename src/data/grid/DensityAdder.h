/**
 * @file DensityAdder.h
 *
 * @date Sep 5, 2016
 * @author Michael Boeckers
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

#ifndef BASICS_GRID_DATAONGRID_DENSITYADDER_H_
#define BASICS_GRID_DATAONGRID_DENSITYADDER_H_

/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/matrices/DensityMatrix.h"
#include "settings/Options.h"

namespace Serenity {
/**
 * @brief Class to add densities which are not defined on the same grid.
 */
template<Options::SCF_MODES SCFMode>
class DensityAdder {
 public:
  DensityAdder() = default;
  virtual ~DensityAdder() = default;
  /**
   *
   * @param densityToBeAddedTo     Density to which shall be added
   * @param densityMatrix          The density matrix from which the density which shall be added
   *                               is calculated
   * @param basisController        The basis controller corresponding to the density matrix
   * @param blockAverageThreshold  Grid blocks of basis function values \f$ \chi_i(r) \f$ are neglected if
   *                               the average value for all points of this block is lower than this threshold.
   * @brief                        Calculates a density from densityMatrix and  it's corresponding
   *                               basis functions on the grid given by densityToBeAddedTo. Use this
   *                               function  only if the densities you want to add are not defined
   *                               on the same grid.
   */
  static void add(DensityOnGrid<SCFMode>& densityToBeAddedTo, const DensityMatrix<SCFMode>& densityMatrix,
                  std::shared_ptr<BasisController> basisController, const double blockAverageThreshold);
};

} /* namespace Serenity */

#endif /* BASICS_GRID_DATAONGRID_DENSITYADDER_H_ */
