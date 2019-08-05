/**
 * @file DirectCD.cpp
 *
 * @date Jun 1, 2016
 * @author Michael Boeckers
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "math/linearAlgebra/DirectCD.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <functional>


namespace Serenity {

DirectCD::DirectCD(
    VectorXd diagonal,
    std::function<Matrix<double>(VectorXi& qualifiedSet)> matrixColumnCalculator,
    double decompositionThreshold,
    double spanFactor) :
      _diagonal(diagonal),
      _matrixColumnCalculator(matrixColumnCalculator),
      _decompositionThreshold(decompositionThreshold),
      _spanFactor(spanFactor),
      _nDimension(diagonal.rows()) {
  _choleskyVectors.resize(_nDimension,_nDimension);
}

void DirectCD::decompose() {
  //Triplet list for the construction of the sparse Cholesky matrix
  typedef Triplet<double> T;
  std::vector<T> tripletList;
  //Largest diagonal element
  double dMax = std::abs(_diagonal.maxCoeff());
  //Iterate over sets of qualified diagonals until dMax < _decompositionThreshold
  //ToDo: It may be, that the first qualified consists of the whole matrix or a
  //      large block of the matrix. There should be a possibility to write the
  //      Cholesky vectors into a file. As suggested by Aquilante et. al., one can
  //      introduced a maximum dimension of 100 for the qualified set.
  do {
    //smallest diagonal element that may be treated
    double dMin = std::max(_spanFactor * dMax, _decompositionThreshold);
    //find qualified set
    VectorXi qualifiedSet;
    unsigned int nQualified;
    computeQualified(dMin, qualifiedSet, nQualified);
    //get matrix columns of qualified diagonals
    //Note: In principle, one could perform this algorithm column by column and never
    //      keep the complete qualified set in memory (i.e. each qualified set only
    //      consists of one column). However, I assume that this will be very slow.
    //ToDo: Build significant set for prescreening, see Aquilante et. al
    Matrix<double> qualifiedCols = _matrixColumnCalculator(qualifiedSet);
    //subtract contributions from previous vectors
    if (_nCholesky > 0) {
      for (unsigned int p = 0; p < _nDimension; ++p) {
        for (unsigned int q = 0; q < nQualified; ++q) {
          for (unsigned int j = 0; j < _nCholesky; ++j) {
           qualifiedCols(p,q) -= _choleskyVectors.coeffRef(p, j) *
                _choleskyVectors.coeffRef(qualifiedSet(q), j);
          }
        }
      }
    }
    unsigned int iAdd = 0;
    //compute largest diagonal among qualified
    double qMax;
    unsigned int qMaxIndex;
    findQmax(qMax, qMaxIndex, qualifiedSet);
    //decompose
    do {
      //calculate Cholesky vector
      double qMaxFactor = 1.0 / std::sqrt(qMax);
      unsigned int iCho = _nCholesky + iAdd;
      for (unsigned int p = 0; p < _nDimension; ++p) {
        tripletList.push_back(T(p,iCho,qMaxFactor*qualifiedCols(p,qMaxIndex)));
      }
      //make Cholesky vectors from triplets
      _choleskyVectors.setFromTriplets(tripletList.begin(), tripletList.end());
      Matrix<double> test = MatrixXd(_choleskyVectors);
      //update diagonal
      //ToDo: The function coeffRef involves a quite expensive binary search.
      //      There must be a better way to access the current Cholesky vector!
      for (unsigned int p = 0; p < _nDimension; ++p) {
        _diagonal(p) -= _choleskyVectors.coeffRef(p,iCho)
            * _choleskyVectors.coeffRef(p,iCho);
      }
      //update qualified columns
      for (unsigned int p = 0; p < _nDimension; ++p) {
        for (unsigned int q = 0; q < nQualified; ++q) {
          qualifiedCols(p,q) -= _choleskyVectors.coeffRef(p,iCho)*
              _choleskyVectors.coeffRef(qualifiedSet(q),iCho);
        }
      }
      //set treated diagonal element to zero
      _diagonal(qualifiedSet(qMaxIndex)) = 0.0;
      //set diagonal elements below threshold to zero
      for (unsigned int p = 0; p < _nDimension; ++p) {
        if (_diagonal(p) < _negativeThreshold) {
          if (_diagonal(p) < _negativeFailThreshold) {
            throw SerenityError("DirectCD: Matrix not positive definite");
          } else {
            _diagonal(p) = 0.0;
          }
        }
      }
      //update
      findQmax(qMax, qMaxIndex, qualifiedSet);
      dMax = _diagonal.maxCoeff();
      dMin = std::max(_spanFactor * dMax, _decompositionThreshold);
      iAdd += 1;
    } while (iAdd < nQualified and qMax > dMin);
    _nCholesky += iAdd;
  } while (dMax > _decompositionThreshold);
  //If not already compressed, compress Cholesky vectors
  _choleskyVectors.makeCompressed();
}

void DirectCD::findQmax(double& qMax, unsigned int& qMaxIndex, VectorXi& qualifiedSet) {
  qMax = _diagonal(qualifiedSet(0));
  qMaxIndex = 0;
  for (unsigned int q = 1; q < qualifiedSet.rows(); ++q) {
    if (_diagonal(qualifiedSet(q)) >= qMax) {
      qMaxIndex = q;
      qMax = _diagonal(qualifiedSet(q));
    }
  }
}

void DirectCD::computeQualified(double& dMin, VectorXi& qualifiedSet, unsigned int& nQualified) {
  //create pair-vector to sort diagonal elements and
  //retain original indices, i.e. the pivoting indices
  //and count number of qualified diagonals
  std::vector<std::pair<double, unsigned int> > pairVector;
  nQualified = 0;
  for (unsigned int i = 0; i < _nDimension; ++i) {
    pairVector.push_back(std::make_pair(_diagonal(i), i));
    if (_diagonal(i) > dMin) {
      nQualified += 1;
    }
  }
  //sort in decending order
  std::stable_sort(pairVector.rbegin(), pairVector.rend());
  //write original indices to pivoting vector
  //and calculate qualified set
  qualifiedSet.resize(nQualified);
  //reset counter and use as index for qualified set
  nQualified = 0;
  _pivotingIndices.resize(_nDimension);
  for (unsigned int i = 0; i < _nDimension; ++i) {
    _pivotingIndices(i) = pairVector[i].second;
    if (_diagonal(_pivotingIndices(i)) > dMin) {
      qualifiedSet(nQualified) = _pivotingIndices(i);
      nQualified += 1;
    } else {
      break;
    }
  }
}

SparseMatrix<double>& DirectCD::getCholeskyVectors() {
  return _choleskyVectors;
}

} /* namespace Serenity */
