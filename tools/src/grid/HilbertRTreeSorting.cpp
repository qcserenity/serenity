/**
 * @file   HilbertRTreeSorting.cpp
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

/* Include Class Header*/
#include "grid/HilbertRTreeSorting.h"
/* Include Std and External Headers */
#include <array>
#include <iostream>

namespace Serenity {

HilbertRTreeSorting::HilbertRTreeSorting(Eigen::Matrix3Xd& points, Eigen::VectorXd& weights)
  : _points(points), _weights(weights), _nPoints(points.cols()) {
  assert(_nPoints == weights.size());
  _depth = 1;
  _length = 2;
  _nVert = 8;
  while (_nPoints * 8 > _nVert) {
    _depth += 1;
    _length *= 2;
    _nVert *= 8;
  }
  for (unsigned int i = 0; i < 3; i++) {
    _min[i] = points.row(i).minCoeff();
    _spread[i] = double(_length) / (points.row(i).maxCoeff() - _min[i]);
  }
}

void HilbertRTreeSorting::sort() {
  Eigen::VectorXi index(_nPoints);

  double* pts = _points.data();

  /*
   *  the transformation matrix for all hilbert cubes to the
   *   next level of depth
   */
  constexpr std::array<std::array<int, 8>, 8> _trans{{
      {{0, 7, 6, 1, 2, 5, 4, 3}},
      {{0, 3, 4, 6, 7, 5, 2, 1}},
      {{0, 3, 4, 6, 7, 5, 2, 1}},
      {{2, 3, 0, 1, 6, 7, 4, 5}},
      {{2, 3, 0, 1, 6, 7, 4, 5}},
      {{6, 5, 2, 1, 0, 3, 4, 7}},
      {{6, 5, 2, 1, 0, 3, 4, 7}},
      {{4, 3, 2, 5, 6, 1, 0, 7}},
  }};

  /*
   * the number of the current point in the first hilbert cube based
   *   on its relation to the center
   *  number = _val[>0.5x][>0.5y][>0.5z]
   */
  constexpr std::array<std::array<std::array<int, 2>, 2>, 2> _val{{
      {{{{5, 6}}, {{4, 7}}}},
      {{{{2, 1}}, {{3, 0}}}},
  }};

#pragma omp parallel for
  for (int i = 0; i < _nPoints; i++) {
    Eigen::Vector3i ipt;
    int* pt = ipt.data();
    int shift(i * 3);
    pt[0] = int((pts[shift + 0] - _min[0]) * _spread[0]);
    pt[1] = int((pts[shift + 1] - _min[1]) * _spread[1]);
    pt[2] = int((pts[shift + 2] - _min[2]) * _spread[2]);

    int l = _length / 2;
    int gthx(pt[0] > l);
    int gthy(pt[1] > l);
    int gthz(pt[2] > l);
    int v = _val[gthx][gthy][gthz];
    int idx = v;
    while (l > 1) {
      idx *= 8;
      if (gthx)
        pt[0] -= l;
      if (gthy)
        pt[1] -= l;
      if (gthz)
        pt[2] -= l;
      l /= 2;
      gthx = pt[0] > l;
      gthy = pt[1] > l;
      gthz = pt[2] > l;
      int x = _val[gthx][gthy][gthz];
      v = _trans[v][x];
      idx += v;
    }
    index.data()[i] = idx;
  }

  /* =======================
   *   Parallel Merge Sort
   * =======================*/

  // initial setup
  // number of cores
  int nprocs = omp_get_max_threads();
  // number of nodes to start sorting is one for each core
  int nnodes = std::min((int)(_nPoints / 2), nprocs);
  // size per node
  Eigen::VectorXi size(nnodes);
  for (int i = 0; i < nnodes; i++)
    size[i] = int(_nPoints / nnodes);
  size[0] += _nPoints % nnodes;
  // starting offset per node
  Eigen::VectorXi start = Eigen::VectorXi::Zero(nnodes);
  for (int i = 1; i < nnodes; i++) {
    start[i] = start[i - 1] + size[i - 1];
  }
  // list to hold the sorted grid indices
  //   this list assures that the points and weights only
  //   have to be move once at the end in a final parallel
  //   O(n) scaling operation.
  Eigen::VectorXi sortedidx(_nPoints);

  // distribute and sort blocks/nodes of points
#pragma omp parallel for schedule(static, 1)
  for (int i = 0; i < nnodes; i++) {
    Eigen::VectorXi newidx(size[i]);
    for (int pt = 0; pt < size[i]; pt++) {
      int idx;
      newidx.data()[pt] = index.segment(start[i], size[i]).maxCoeff(&idx);
      sortedidx.data()[pt + start[i]] = idx + start[i];
      index.data()[idx + start[i]] = -1;
    }
    // regenerate initial indix list  with new ones
    index.segment(start[i], size[i]) = newidx;
  }
  // collapse/merge tree
  while (size.size() > 1) {
    int newnnodes = int(nnodes / 2) + nnodes % 2;
    Eigen::VectorXi newstart(newnnodes);
    Eigen::VectorXi newsize(newnnodes);
    if (nnodes % 2 == 1) {
      newstart[newnnodes - 1] = start[start.size() - 1];
      newsize[newnnodes - 1] = size[size.size() - 1];
    }

#pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < int(nnodes / 2); i++) {
      int counta = 0;
      int countb = 0;
      int sizea = size[2 * i];
      int sizeb = size[2 * i + 1];
      int starta = start[2 * i];
      int startb = start[2 * i + 1];
      newstart[i] = starta;
      newsize[i] = sizea + sizeb;
      Eigen::VectorXi tmpsort(newsize[i]);
      Eigen::VectorXi tmpidx(newsize[i]);
      /*
       * merge
       */
      for (int pt = 0; pt < newsize[i]; pt++) {
        if (index.data()[starta + counta] > index.data()[startb + countb]) {
          tmpsort.data()[pt] = sortedidx.data()[starta + counta];
          tmpidx.data()[pt] = index.data()[starta + counta];
          counta++;
        }
        else {
          tmpsort.data()[pt] = sortedidx.data()[startb + countb];
          tmpidx.data()[pt] = index.data()[startb + countb];
          countb++;
        }
        // if one list is empty, append the rest of the other and stop
        if (counta == sizea) {
          int remaining = sizeb - countb;
          tmpsort.segment(sizea + countb, remaining) = sortedidx.segment(startb + countb, remaining);
          tmpidx.segment(sizea + countb, remaining) = index.segment(startb + countb, remaining);
          break;
        }
        else if (countb == sizeb) {
          int remaining = sizea - counta;
          tmpsort.segment(sizeb + counta, remaining) = sortedidx.segment(starta + counta, remaining);
          tmpidx.segment(sizeb + counta, remaining) = index.segment(starta + counta, remaining);
          break;
        }
      }
      sortedidx.segment(newstart[i], newsize[i]) = tmpsort;
      index.segment(newstart[i], newsize[i]) = tmpidx;
    }
    size = newsize;
    start = newstart;
    nnodes = newnnodes;
  }

  // move points to final destination
  //   TODO this could be done without a copy (there should be chains to follow)
  Eigen::Matrix3Xd pcopy(_points);
  Eigen::VectorXd wcopy(_weights);
  for (int i = 0; i < _nPoints; i++) {
    _points.col(i) = pcopy.col(sortedidx[i]);
    _weights[i] = wcopy[sortedidx[i]];
  }
}

} /* namespace Serenity */
