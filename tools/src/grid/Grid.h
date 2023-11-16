/**
 * @file   Grid.h
 * @author Thomas Dresselhaus
 *
 * @date 14. März 2014, 21:05
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
#ifndef GRID_H
#define GRID_H
/* Include Serenity Internal Headers */
#include "geometry/Point.h"
#include "grid/HilbertRTreeSorting.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <cassert>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class GridController;
/**
 * @class Grid Grid.h
 * @brief A grid usable e.g. for numerical integration or putting 3-D data into easily plottable data
 */
class Grid {
  friend GridController;

 public:
  /**
   * @param gridPoints real points in space which define the grid
   * @param weights basically the volume that each grid point shall represent; in the end we want
   *        to represent the whole R^3 space
   */
  Grid(std::unique_ptr<Eigen::Matrix3Xd> gridPoints, std::unique_ptr<Eigen::VectorXd> weights)
    : // These should call the move constructors of the vectors.
      _gridPoints(*std::move(gridPoints)),
      _weights(*std::move(weights)),
      _nPoints(_weights.size()),
      _sorted(false) {
    assert(_weights.size() == (int)_gridPoints.cols());
  };

  virtual ~Grid() = default;

  /**
   * @brief Uses a Hilbert R Tree to sort the grid.
   */
  void sort() {
    Timings::takeTime("Tech. -           Grid Sorting");
    HilbertRTreeSorting sorter(_gridPoints, _weights);
    sorter.sort();
    _sorted = true;
    Timings::timeTaken("Tech. -           Grid Sorting");
  }

 private:
  Eigen::Matrix3Xd _gridPoints;
  /*
   * The weights are not of type GridData, because they are not 'defined on a grid', but rather are
   * a part of the grid definition.
   */
  Eigen::VectorXd _weights;
  const unsigned int _nPoints;

 protected:
  bool _sorted;
};

} // namespace Serenity

#endif /* GRID_H */
