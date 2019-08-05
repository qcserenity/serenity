/**
 * @file   Orthogonalization.h
 *
 * @date   Dec 7, 2015
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

#ifndef ORTHOGONALIZATION_H_
#define ORTHOGONALIZATION_H_

/* Include Std and External Headers */
#include <Eigen/Dense>


namespace Serenity {
/**
 * @class  Orthogonalization Orthogonalization.h
 * @brief  Orthogonalization routines
 */
class Orthogonalization {
public:
  /**
   *
   * @param A Matrix with column vectors to be orthogonalized
   * @brief Modified Gram-Schmidt orthogonalization. Scales with 2*n*k, where n is their
   *        dimension and k their number. Note that the column vectors will not be normalized.
   */
  static void modifiedGramSchmidt(Eigen::MatrixXd& A) {
    const unsigned int nVector = A.cols();
    Eigen::VectorXd projectionOperator(A.rows());
    for (unsigned int i = 0; i < nVector; ++i) {
      for (unsigned int j = i + 1; j < nVector; ++j) {
        projectionOperator = A.col(i).dot(A.col(j)) / A.col(i).dot(A.col(i)) * A.col(i);
        A.col(j) -= projectionOperator;
      }
    }
  }
  /**
   *
   * @param A Matrix with column vectors to be orthogonalized
   * @brief Orthogonalization by QR decomposition using Eigen3. Note thate this orthogonalization
   *        routine is more expensive than Gram-Schmidt.
   */
  static void qr(Eigen::MatrixXd& A) {
    A = A.householderQr().householderQ();
  }
  /**
   *
   * @param A
   * @param B
   * @brief Impose biorthogonalization condition between column vectors of matrix A and B,
   *        i.e. the usual biorthogonality condition for left and right eigenvectors of
   *        non-Hermitian matrices.
   *        \f[
   *           \left( \vec{a}_i, \vec{b}_i \right)= \delta_{ij}
   *        \f]
   *        Note that the column vectors will not be normalized.
   */
  static void modifiedGramSchmidtBiorthogonalization(Eigen::MatrixXd& A, Eigen::MatrixXd& B) {
    const unsigned int nVector = 2 * A.cols();
    Eigen::VectorXd projectionOperator(A.rows());
    //combine matrices A and B
    Eigen::MatrixXd combinedMatrix(A.rows(),nVector);
    combinedMatrix.block(0,0,A.rows(),A.cols()) = A;
    combinedMatrix.block(0,A.cols(),B.rows(),B.cols()) = B;
    //biorthonormalize column vectors of A and B
    for (unsigned int i = 0; i < nVector; ++i) {
      for (unsigned int j = A.cols() + i + 1; j < nVector; ++j) {
        projectionOperator = combinedMatrix.col(i).dot(combinedMatrix.col(j)) / combinedMatrix.col(i).dot(combinedMatrix.col(i)) * combinedMatrix.col(i);
        combinedMatrix.col(j) -= projectionOperator;
      }
    }
    //get biorthogonalized matrices A and B from combined matrix
    A = combinedMatrix.block(0,0,A.rows(),A.cols());
    B = combinedMatrix.block(0,A.cols(),B.rows(),B.cols());
  }

  static void modifiedGramSchmidtHessianProjection(Eigen::MatrixXd& A, Eigen::MatrixXd& B) {
    //Gram Schmidt procedure
    unsigned int dim1 = A.cols();
    unsigned int dim2 = A.rows();
    Eigen::MatrixXd c;
    Eigen::MatrixXd alpha(dim1,dim2);
    Eigen::VectorXd tmpVec;
    Eigen::MatrixXd a2(dim1,dim2);

    alpha.setZero();
    a2 = A;
    B = A;


    for(unsigned int i=0; i<dim2;i++){
      c.setZero();
      c = -1 * alpha * (B.block(0,i,dim1,dim2-i)) + B.block(0,i,dim1,dim2-i);
      for(unsigned int j=0; j<c.cols();j++){
        for(unsigned int k=0; k<c.rows();k++){
          B(k,j+i)=c(k,j);
        }
      }

      Eigen::Map<Eigen::VectorXd> colVec (B.col(i).data(),B.rows());
      auto max_norm = colVec.dot(colVec);
      auto max_j = i;

      for(unsigned int ii=i+1; ii<dim2; ii++){
        Eigen::Map<Eigen::VectorXd> colVec (B.col(ii).data(),B.rows());
        auto norms = colVec.dot(colVec);
        if(norms > max_norm){
          max_norm = norms;
          max_j=ii;
        }
      }
      if(i != max_j){
        tmpVec=B.col(i);
        B.col(i) = B.col(max_j);
        B.col(max_j) = tmpVec;
        tmpVec = a2.col(i);
        a2.col(i) = a2.col(max_j);
        a2.col(max_j) = tmpVec;
      }
      if(max_norm > 0.00001){
        auto normfac = 1.0 / sqrt(max_norm);
        B.col(i) *= normfac;
      } else {
        B.col(i).setZero();
      }
      for (unsigned int k=0; k<dim1;k++){
        for (unsigned int l=0; l<=k; l++){
          alpha(k,l) = B(k,i)*B(l,i);
          alpha(l,k) = alpha(k,l);
        }
      }
    }
  }

};
}

#endif /* ORHTOGONALIZATION_H_ */
