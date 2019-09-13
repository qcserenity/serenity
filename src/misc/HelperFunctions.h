/**
 * @file HelperFunctions.h
 *
 * @date Juli 11, 2015
 * @author Thomas Dresselhaus, Michael Boeckers
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
#ifndef HELPERFUNCTIONS_H
#define	HELPERFUNCTIONS_H

/* Include Serenity Internal Headers */
#include "math/Matrix.h"
#include <Eigen/SparseCore>


namespace Serenity {
/**
 * @param target a container of the data from which the result is produced
 * @returns a container of reference_wrappers holding references to the objects in target
 */
template<class ContentT>
inline std::vector<std::reference_wrapper<ContentT> >
makeReferenceContainer(std::vector<ContentT>& target) {
  return std::vector<std::reference_wrapper<ContentT> >(target.begin(), target.end());
}
/**
 * @brief analog to makeReferenceContainer, which takes a vector of shared_pointers
 */
template<class ContentT>
inline std::vector<std::reference_wrapper<ContentT> >
makeReferenceContainerFromPtr(std::vector<std::shared_ptr<ContentT> >& target) {
  std::vector<std::reference_wrapper<ContentT> > result;
  for (const auto& element : target) result.push_back(*element);
  return result;
}
/**
 * @brief analog to makeReferenceContainer, which takes a vector of unique_pointers
 */
template<class ContentT>
inline std::vector<std::reference_wrapper<ContentT> >
makeReferenceContainerFromPtr(std::vector<std::unique_ptr<ContentT> >& target) {
  std::vector<std::reference_wrapper<ContentT> > result;
  for (const auto& element : target) result.push_back(*element);
  return result;
}


/// @cond false
template<class It>
inline void _advance_ (It& it) {
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
/**
 * @brief a python-like zip based on the lambda function func. Two arguments (t, u).
 * @param t
 * @param u
 * @param func
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
 * @param t
 * @param u
 * @param v
 * @param func
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
 * @param t
 * @param u
 * @param v
 * @param w
 * @param func
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
 * @brief  Writes a lower triagnular matrix to a vector
 * @param  matrixL A lower triagnular matrix
 * @return Returns the corresponding vector. Note that this routine will not check
 *         whether matrixL is a upper triangular matrix or if it is symmetric or not.
 *         It will just take the lower triangular matrix of any matrixL you give to it and
 *         store this in a vector of appropriate length.
 */
inline Eigen::VectorXd triangularMatrix2Vector(Matrix<double>& matrixL) {
  Eigen::VectorXd triagonalVector(matrixL.cols()*(matrixL.cols()+1)/2);
  for (int i = 0,ij=0; i < matrixL.cols(); ++i) {
    for (int j = 0; j <= i; ++j,++ij) {
      triagonalVector(ij) = matrixL(i,j);
    }
  }
  return triagonalVector;
}
 /**
  *
  * @param triangularVector A vector
  * @return                 Returns the corresponding lower triangular matrix
  */
inline Matrix<double> vector2triangularMatrix(Eigen::VectorXd triangularVector) {
  int nDimension = int(-0.5 + std::sqrt(0.25 + 2.0 * triangularVector.rows()));
  Matrix<double> triangularMatrix(nDimension,nDimension);
  triangularMatrix.setZero();
  for (int i = 0,ij=0; i < nDimension; ++i) {
    for (int j = 0; j <= i;++j,++ij) {
      triangularMatrix(i,j) = triangularVector(ij);
    }
  }
  return triangularMatrix;
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
 * @return The projection matrix.
 */
inline Eigen::SparseMatrix<double> constructProjectionMatrix(
    Eigen::VectorXi negligible) {
  negligible -= Eigen::VectorXi::Constant(negligible.rows(),1);
  negligible = negligible.array().abs();
  Eigen::SparseVector<int> blockToBasis (negligible.rows());
  blockToBasis = negligible.sparseView();
  unsigned int nBasisFunctions = negligible.size();

  std::vector<Eigen::Triplet<double> > projectionTriplets;
  unsigned int col = 0;
  for(Eigen::SparseVector<int>::InnerIterator itMu(blockToBasis); itMu; ++itMu) {
    projectionTriplets.push_back(Eigen::Triplet<double>(itMu.row(),col,1.0));
    col++;
  }// for itMu
  Eigen::SparseMatrix<double> projectionMatrix(
      nBasisFunctions,blockToBasis.nonZeros());
  projectionMatrix.setFromTriplets(projectionTriplets.begin(),projectionTriplets.end());
  return projectionMatrix;
}

/**
 * @brief Performs Mulliken net population based prescreening and constructs the associated
 *        projection matrix.
 * @param coefficients  The coefficients to be prescreened.
 * @param mnpThreshold  The prescreening threshold.
 * @return The projection matrix.
 */
inline Eigen::SparseMatrix<double> constructSignificantMnPProjectionMatrix(
    Eigen::VectorXd coefficients,
    double mnpThreshold) {
  unsigned int nC2s = coefficients.size();
  Eigen::VectorXd c2s = coefficients.cwiseProduct(coefficients);
  std::vector<Eigen::Triplet<double> > projectionTriplets;
  unsigned int col = 0;
  for(unsigned int iC2 = 0; iC2 < nC2s; ++iC2) {
    if(c2s(iC2) > mnpThreshold) {
      projectionTriplets.push_back(Eigen::Triplet<double>(iC2,col,1.0));
      col++;
    }
  }// for itMu
  Eigen::SparseMatrix<double> projectionMatrix(
      nC2s,col);
  projectionMatrix.setFromTriplets(projectionTriplets.begin(),projectionTriplets.end());
  return projectionMatrix;
}


} /* namespace Serenity */

#endif	/* HELPERFUNCTIONS_H */
