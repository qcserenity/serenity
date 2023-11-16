/**
 * @file   HilbertRTreeSorting.h
 *
 * @date    Sep 29, 2017
 * @author  Jan Unsleber
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

#ifndef GRID_HILBERTRTREESORTING_H_
#define GRID_HILBERTRTREESORTING_H_
/* Include Serenity Internal Headers */
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {
/**
 * @class HilbertRTreeSorting HilbertRTreeSorting.h
 * @brief An implementation of a sorting algorithm for points in space, using a 3D Hilbert R-Tree.
 *
 * The algorithm for n points (4n doubles)
 *  addditional 4n doubles and n integers of memory
 * The generatrion of hilbert indices scales O(d*n) where d is the tree depth.
 * The following sorting algorithm should scale O(n*2) and ist parallel at the moment.
 * (-JU, Sep 29, 2017)
 */
class HilbertRTreeSorting {
 public:
  HilbertRTreeSorting(Eigen::Matrix3Xd& points, Eigen::VectorXd& weights);

  /// @brief Default destructor
  virtual ~HilbertRTreeSorting() = default;

  /// @brief Sorts the points and weights in place.
  void sort();

 private:
  // the coordinates to be sorted
  Eigen::Matrix3Xd& _points;
  // the weight to be sorted
  Eigen::VectorXd& _weights;
  // the minimal value for each cartesian dimension
  Eigen::Vector3d _min;
  // the conversions from the actual spread to the spread in
  //   the space the hilbert function is defined in
  Eigen::Vector3d _spread;
  // the depth of the hilbert tree
  int _depth;
  // number of points in the coordinates
  int _nPoints;
  // the number of vertices in the hilbert curve
  int _nVert;
  // the length of each side of the cube the curve is defined in
  int _length;
};

} /* namespace Serenity */

#endif /* GRID_HILBERTRTREESORTING_H_ */
