/**
 * @file HelperFunctions.h
 *
 * @date Juli 11, 2015
 * @author Thomas Dresselhaus, Michael Boeckers
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
#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

/* Include Serenity Internal Headers */
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <Eigen/SparseCore>
#include <sstream>

namespace Serenity {
/**
 * @brief Returns a string of a number with a adjustable precision, from:
 * https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values/16606128#16606128
 * @param a_value the value
 * @param n decimal precision
 * @returns the value as string
 */
template<typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

/**
 * @brief a python-like zip based on the lambda function func. Two arguments (t, u).
 * @param t Argument 1
 * @param u Argument 2
 * @param func ziplike function
 */
template<class T, class U, class FuncType>
inline void zip(T& t, U& u, const FuncType& func) {
  auto tIter = t.begin();
  auto uIter = u.begin();
  _zipIts_(func, tIter, t.end(), uIter);
  assert(!(uIter != u.end()));
}
/**
 * @brief a python-like zip based on the lambda function func. Three arguments (t, u, v).
 * @param t Argument 1
 * @param u Argument 2
 * @param v Argument 3
 * @param func ziplike function
 */
template<class T, class U, class V, class FuncType>
inline void zip(T& t, U& u, V& v, const FuncType& func) {
  auto tIter = t.begin();
  auto uIter = u.begin();
  auto vIter = v.begin();
  _zipIts_(func, tIter, t.end(), uIter, vIter);
  assert(!(uIter != u.end()));
  assert(!(vIter != v.end()));
}
/**
 * @brief a python-like zip based on the lambda function func. Four arguments (t, u, v, w).
 * @param t Argument 1
 * @param u Argument 2
 * @param v Argument 3
 * @param w Argument 4
 * @param func ziplike function
 */
template<class T, class U, class V, class W, class FuncType>
inline void zip(T& t, U& u, V& v, W& w, const FuncType& func) {
  auto tIter = t.begin();
  auto uIter = u.begin();
  auto vIter = v.begin();
  auto wIter = w.begin();
  _zipIts_(func, tIter, t.end(), uIter, vIter, wIter);
  assert(!(uIter == u.end()));
  assert(!(vIter == v.end()));
  assert(!(wIter == w.end()));
}
/**
 * @brief A helper function that constructs a projection matrix from a prescreening
 *        vector <negligible>.\n
 *        The task of the projection matrix \f$ \mathbf{P} \f$ is to transform a large matrix
 *        \f$ \mathbf{M} \f$ with many redundant or negligible entries to a smaller matrix
 *        \f$ \mathbf{M}^\mathrm{sig}\f$ that only contains the important contributions.\n
 *        For example consider the basis function values on a grid for a given block
 *        of grid points. For the block only a small number of basis functions will
 *        have values different from zero. In the code only these values are calculated.
 *        However, the resulting matrix \f$ \mathbf{B} \f$ (block size x N-basis functions)
 *        will contain columns for all basis functions. The important columns can now
 *        be easily obtained by multiplying with the sparse projection matrix:\n\n
 *
 *               \f$ \mathbf{BP} = \mathbf{B}^\mathrm{sig} \f$\n\n
 *
 *        After evaluation in the reduced dimension the result can be easily projected back
 *        into the original dimension by multiplying with \f$ \mathbf{P}^\dagger \f$.\n\n
 *
 *        Note that this multiplication is extremely efficient since the sparse nature
 *        of the projection matrix is fully exploited.
 * @param negligible The prescreening vector.
 * @param allowXZero If true matrices with 0 columns may be generated, otherwise a nx1 zero vector is returned.
 * @return The projection matrix.
 */
inline Eigen::SparseMatrix<double> constructProjectionMatrix(Eigen::VectorXi negligible, bool allowXZero = false) {
  negligible -= Eigen::VectorXi::Constant(negligible.rows(), 1);
  negligible = negligible.array().abs();
  Eigen::SparseVector<int> blockToBasis(negligible.rows());
  blockToBasis = negligible.sparseView();
  unsigned int nBasisFunctions = negligible.size();

  std::vector<Eigen::Triplet<double>> projectionTriplets;
  unsigned int col = 0;
  for (Eigen::SparseVector<int>::InnerIterator itMu(blockToBasis); itMu; ++itMu) {
    projectionTriplets.push_back(Eigen::Triplet<double>(itMu.row(), col, 1.0));
    col++;
  } // for itMu
  Eigen::SparseMatrix<double> projectionMatrix(nBasisFunctions, blockToBasis.nonZeros());
  if (blockToBasis.nonZeros() == 0 && not allowXZero) {
    projectionMatrix.resize(nBasisFunctions, 1);
  }
  projectionMatrix.setFromTriplets(projectionTriplets.begin(), projectionTriplets.end());
  return projectionMatrix;
}

/**
 * @brief Constructs a projection matrix from a sparse-non-negligible vector.
 *        See constructProjectionMatrix.
 * @param nonNegligible The prescreenion vector.
 * @return The sparse projection.
 */
inline Eigen::SparseMatrix<double> constructProjectionMatrixFromSparse(const Eigen::SparseVector<int>& nonNegligible) {
  std::vector<Eigen::Triplet<double>> projectionTriplets;
  unsigned int col = 0;
  for (Eigen::SparseVector<int>::InnerIterator it(nonNegligible); it; ++it) {
    projectionTriplets.push_back(Eigen::Triplet<double>(it.row(), col, 1.0));
    col++;
  } // for it
  Eigen::SparseMatrix<double> projectionMatrix(nonNegligible.rows(), nonNegligible.nonZeros());
  projectionMatrix.setFromTriplets(projectionTriplets.begin(), projectionTriplets.end());
  return projectionMatrix;
}
/**
 * @brief Constructs a projection matrix from a sparse-non-negligible vector.
 *        See constructProjectionMatrix and constructProjectionMatrixFromSparse.
 *        The matrix is transposed with respect to constructProjectionMatrixFromSparse().
 * @param nonNegligible The prescreenion vector.
 * @return The sparse projection (transposed).
 */
inline Eigen::SparseMatrix<double> constructProjectionMatrixFromSparse_T(const Eigen::SparseVector<int>& nonNegligible) {
  std::vector<Eigen::Triplet<double>> projectionTriplets;
  unsigned int row = 0;
  for (Eigen::SparseVector<int>::InnerIterator it(nonNegligible); it; ++it) {
    projectionTriplets.push_back(Eigen::Triplet<double>(row, it.row(), 1.0));
    row++;
  } // for it
  Eigen::SparseMatrix<double> projectionMatrix(nonNegligible.nonZeros(), nonNegligible.rows());
  projectionMatrix.setFromTriplets(projectionTriplets.begin(), projectionTriplets.end());
  return projectionMatrix;
}

/**
 * @brief Performs Mulliken net population based prescreening and constructs the associated
 *        projection matrix.
 * @param coefficients  The coefficients to be prescreened.
 * @param mnpThreshold  The prescreening threshold.
 * @return The projection matrix.
 */
inline Eigen::SparseMatrix<double> constructSignificantMnPProjectionMatrix(Eigen::VectorXd coefficients, double mnpThreshold) {
  unsigned int nC2s = coefficients.size();
  Eigen::VectorXd c2s = coefficients.cwiseProduct(coefficients);
  std::vector<Eigen::Triplet<double>> projectionTriplets;
  unsigned int col = 0;
  for (unsigned int iC2 = 0; iC2 < nC2s; ++iC2) {
    if (c2s(iC2) > mnpThreshold) {
      projectionTriplets.push_back(Eigen::Triplet<double>(iC2, col, 1.0));
      col++;
    }
  } // for itMu
  Eigen::SparseMatrix<double> projectionMatrix(nC2s, col);
  projectionMatrix.setFromTriplets(projectionTriplets.begin(), projectionTriplets.end());
  return projectionMatrix;
}
/// @cond false
template<class It>
inline void _advance_(It& it) {
  ++it;
}
template<class It, class... Its>
inline void _advance_(It& it, Its&... its) {
  ++it;
  _advance_(its...);
}
template<class FuncType, class It, class... Its>
inline void _zipIts_(const FuncType& func, It& it, const It& end, Its&... its) {
  for (; it != end; _advance_(it, its...)) {
    func(*it, *its...);
  }
}
/// @endcond

} /* namespace Serenity */

#endif /* HELPERFUNCTIONS_H */
