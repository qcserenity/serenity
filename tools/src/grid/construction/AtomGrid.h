/**
 * @file   AtomGrid.h
 * @author Thomas Dresselhaus
 *
 * @date   14. März 2014, 23:21
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
#ifndef ATOMGRID_H
#define ATOMGRID_H
/* Include Serenity Internal Headers */
#include "geometry/Point.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <vector>

namespace Serenity {
/* Forward declarations */
class GridFactory;
/**
 * @class AtomGrid AtomGrid.h
 * @brief An atom-centered grid.
 *
 * From all AtomGrids in the Atoms of a System the actual integration grid will be produced.
 */
class AtomGrid {
 public:
  /**
   * @param gridPoints
   * @param weights
   */
  AtomGrid(Eigen::Matrix3Xd gridPoints, Eigen::VectorXd weights) : _gridPoints(gridPoints), _weights(weights) {
  }
  virtual ~AtomGrid() = default;
  /**
   * @returns the size of the grid, i.e. length of the gridPoint and weight vectors.
   */
  unsigned int nPoints() const {
    return _weights.size();
  }
  /**
   * @returns the grid points
   */
  const Eigen::Matrix3Xd& getGridPoints() const {
    return _gridPoints;
  }
  /**
   * @returns the grid weights
   */
  const Eigen::VectorXd& getWeights() const {
    return _weights;
  }

 private:
  Eigen::Matrix3Xd _gridPoints;
  Eigen::VectorXd _weights;
};

} // namespace Serenity

#endif /* ATOMGRID_H */
