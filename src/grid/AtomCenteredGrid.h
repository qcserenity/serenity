/**
 * @file   AtomCenteredGrid.h
 * @author Thomas Dresselhaus
 *
 * @date 09. März 2016, 15:39
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
#ifndef ATOMCENTEREDGRID_H
#define ATOMCENTEREDGRID_H
/* Include Serenity Internal Headers */
#include "grid/Grid.h"

namespace Serenity {
/* Forward declarations */
class AtomCenteredGridController;
/**
 * @class AtomCenteredGrid AtomCenteredGrid.h
 * @brief Atom-centered version of a numerical integration grid.
 *
 * Most integration grids used in calculations are constructed as a combination of grids for
 * the atoms of the system. However, also other grids (e.g. cubical grids for plotting) are
 * used sometimes. This AtomCenteredGrid stores the contributions of the different atoms to
 * this grid additionally to the normal grid information.\n
 *
 * A similar design exists for the Basis.
 */
class AtomCenteredGrid : public Grid {
  friend AtomCenteredGridController;

 public:
  /**
   * @brief See the constructor of Grid
   * @param gridPoints
   * @param weights
   * @param gridIndicesOfAtoms for each atom used in the construction of the grid the first
   *                           and the end (one after the last) grid point index which belongs
   *                           to that atom.
   */
  AtomCenteredGrid(std::unique_ptr<Eigen::Matrix3Xd> gridPoints, std::unique_ptr<Eigen::VectorXd> weights,
                   const std::vector<std::pair<unsigned int, unsigned int>>& gridIndicesOfAtoms)
    : Grid(std::move(gridPoints), std::move(weights)), _gridIndicesOfAtoms(gridIndicesOfAtoms) {
  }

  virtual ~AtomCenteredGrid() = default;

  const std::vector<std::pair<unsigned int, unsigned int>>& getGridIndicesOfAtoms() const {
    assert(!_sorted);
    return _gridIndicesOfAtoms;
  }

 private:
  /*
   * See the constructor
   */
  const std::vector<std::pair<unsigned int, unsigned int>> _gridIndicesOfAtoms;
};

} // namespace Serenity

#endif /* ATOMCENTEREDGRID_H */
