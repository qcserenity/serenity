/**
 * @file CholeskyDecomposer.cpp
 *
 * @date Jan 25, 2017
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

/* Include Class Header*/
#include "math/linearAlgebra/CholeskyDecomposer.h"
/* Include Serenity Internal Headers */
#include "memory/MemoryManager.h"
#include "misc/SerenityError.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

CholeskyDecomposer::CholeskyDecomposer(
    Eigen::VectorXd& diagonal,
    std::function<std::unique_ptr<Eigen::MatrixXd>(std::vector<unsigned int>& qualifiedSet)> matrixColumnCalculator,
    const double decompositionThreshold, const double screeningDamping, const double spanFactor, const unsigned int nMaxQual,
    const double negativeThreshold, const double negativeFailThreshold, const double pruningThreshold)
  : _diagonal(diagonal),
    _matrixColumnCalculator(matrixColumnCalculator),
    _decompositionThreshold(decompositionThreshold),
    _screeningDamping(screeningDamping),
    _spanFactor(spanFactor),
    _nMaxQual(nMaxQual),
    _negativeThreshold(negativeThreshold),
    _negativeFailThreshold(negativeFailThreshold),
    _pruningThreshold(pruningThreshold),
    _hasBeenCalculated(false) {
  assert(_screeningDamping >= 1.0);
}

std::unique_ptr<Eigen::SparseMatrix<double>> CholeskyDecomposer::getCholeskyVectors() {
  _choleskyVectors =
      std::unique_ptr<Eigen::SparseMatrix<double>>(new Eigen::SparseMatrix<double>(_diagonal.rows(), _diagonal.rows()));
  decompose();
  return std::move(_choleskyVectors);
}

void CholeskyDecomposer::decompose() {
  // Triplet list for the construction of the sparse Cholesky matrix (see Eigen manual)
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  //  tripletList.reserve(_diagonal.rows() * _diagonal.rows());
  // 1. Find dMax
  double dMax = std::abs(_diagonal.maxCoeff());

  // 2. Compute reduced set of significant diagonal elements
  //   We do not work with a reduced set because it prevents the usage
  //   of eigen block operations which are necessary to make the update
  //   of the qualifiedCols efficient. Furthermore, it makes the routine
  //   which calculates the qualifiedCols more complex.

  // 3. Initialize vector counter n
  unsigned int n = 0;

  // 4. Initialize iteration counter
  unsigned int i = 0;

  // 5. While Dmax^(i+1) > tau
  do {
    // 5.a.Update iteration counter
    i += 1;

    // 5.b. Compute smallest diagonal that may be treated
    double dMin = std::max(_spanFactor * dMax, _decompositionThreshold);

    // 5.c. Compute set of qualified diagonals
    std::vector<unsigned int> qualifiedSet = calculateQualifiedSet(dMin);
    // 5.d. Get matrix columns corresponding to qualified diagonals
    auto qualifiedCols = _matrixColumnCalculator(qualifiedSet);

    // 5.e. Subtract contributions from previous vectors (if any)
    if (n > 0) {
      for (unsigned int q = 0; q < qualifiedSet.size(); ++q) {
        Eigen::VectorXd vec = (*_choleskyVectors) * (*_choleskyVectors).transpose().col(qualifiedSet[q]);
        (*qualifiedCols).col(q) -= vec;
      }
    }

    // 5.f. Compute largest diagonal among qualified
    double qMax;
    unsigned int qMaxIndex;
    findQmax(qMax, qMaxIndex, qualifiedSet);

    // 5.g. Initialize Counter for cholesky vectors
    unsigned int j = 0;

    // 5.h. While j , dim(D^(i)) and Qmax > Dmin^(i)
    do {
      // 5.h.i. Update counter j=j+1, J = n + j
      j += 1;
      unsigned int J = n + j - 1;
      // 5.h.iii Calculate Cholseky vector (and store in decomposed column of qualified set)
      double qMaxFactor = 1.0 / std::sqrt(qMax);
      for (unsigned int p = 0; p < _diagonal.rows(); ++p) {
        if (_diagonal(p) == 0.0) {
          (*qualifiedCols)(p, qMaxIndex) = 0.0;
        }
        else {
          (*qualifiedCols)(p, qMaxIndex) *= qMaxFactor;
        }
      }
      // 5.h.ii Assign Cholesky basis function h_J to [q]_J, the index corresponding to Q_max
      //[q]_J (qMaxIndex) has been calculated by function findQMax

      // 5.h.iv. Update
      _diagonal.array() -= (*qualifiedCols).col(qMaxIndex).array().square();
      // set treated diagonal element to zero
      _diagonal(qualifiedSet[qMaxIndex]) = 0.0;
      // set diagonal elements below threshold to zero
      for (unsigned int p = 0; p < _diagonal.rows(); ++p) {
        if (_diagonal(p) < _negativeThreshold) {
          if (_diagonal(p) < _negativeFailThreshold) {
            throw SerenityError("CholeskyDecomposer: Matrix not positive definite");
          }
          else {
            _diagonal(p) = 0.0;
          }
        }
      }

      for (unsigned int q = 0; q < qualifiedSet.size(); ++q) {
        if (_diagonal(qualifiedSet[q]) == 0.0)
          continue;
        (*qualifiedCols).col(q) -= (*qualifiedCols)(qualifiedSet[q], qMaxIndex) * (*qualifiedCols).col(qMaxIndex);
      }

      // save cholesky vector in sparse matrix
      for (unsigned int p = 0; p < _diagonal.rows(); ++p) {
        if ((*qualifiedCols)(p, qMaxIndex) == 0.0)
          continue;
        tripletList.push_back(T(p, J, (*qualifiedCols)(p, qMaxIndex)));
      }

      findQmax(qMax, qMaxIndex, qualifiedSet);

    } while (j < qualifiedSet.size() and qMax > dMin); /*5.h. While j , dim(D^(i)) and Qmax > Dmin^(i) */
    // 5.i. Update vector counter n = n + j
    n += j;

    // 5.j. Compute larget diagonal
    dMax = std::abs(_diagonal.maxCoeff());

    // 5.k Compute reduces set of significant diagonal elements
    //    Not used here, see comment 2.

    // Make Cholesky vectors from triplets
    // This copying step is necessary since we cannot work with the insert method (see eigen manual).
    // For the insert method, we would a priori need to know the number of non-zeros for each column
    // in order to reserve the memory. However, since we use a pivoting algorithm, it is not possible
    // to determine the number of non zeros before the calculation. One thus would have to reserve
    // memory for the maximum number of non-zeros which means that we would have a memory demand equal
    // to storing the complete matrix we want to decompose.
    // ToDo: Check if enough memory is available. If so, use insert mode.
    (*_choleskyVectors).setFromTriplets(tripletList.begin(), tripletList.end());
    // Remove too small entries
    (*_choleskyVectors).prune(_pruningThreshold, 1.0);
    // If not already compressed, compress Cholesky vectors
    (*_choleskyVectors).makeCompressed();

  } while (dMax > _decompositionThreshold); /* 5. While Dmax^(i+1) > tau */
  _hasBeenCalculated = true;
}

std::vector<unsigned int> CholeskyDecomposer::calculateQualifiedSet(double dMin) {
  // create pair-vector to sort diagonal elements and
  // retain original indices, i.e. the pivoting indices
  std::vector<std::pair<double, unsigned int>> pairVector;
  for (unsigned int i = 0; i < _diagonal.rows(); ++i) {
    pairVector.push_back(std::make_pair(_diagonal(i), i));
  }
  // sort in decending order
  std::stable_sort(pairVector.rbegin(), pairVector.rend());

  // calculate qualified set
  std::vector<unsigned int> qualifiedSet;
  for (unsigned int i = 0; i < _diagonal.rows(); ++i) {
    if (i > _nMaxQual)
      break;
    if (_diagonal(pairVector[i].second) > dMin) {
      qualifiedSet.push_back(pairVector[i].second);
    }
    else {
      break;
    }
  }
  return qualifiedSet;
}

void CholeskyDecomposer::findQmax(double& qMax, unsigned int& qMaxIndex, std::vector<unsigned int>& qualifiedSet) {
  qMax = _diagonal(qualifiedSet[0]);
  qMaxIndex = 0;
  for (unsigned int q = 1; q < qualifiedSet.size(); ++q) {
    if (_diagonal(qualifiedSet[q]) > qMax) {
      qMaxIndex = q;
      qMax = _diagonal(qualifiedSet[q]);
    }
  }
}

} /* namespace Serenity */
