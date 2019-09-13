/**
 * @file MatrixFunctions.h
 *
 * @date Apr 17, 2019
 * @author: Moritz Bensberg
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

#ifndef MATH_LINEARALGEBRA_MATRIXFUNCTIONS_H_
#define MATH_LINEARALGEBRA_MATRIXFUNCTIONS_H_

/* Include Std and External Headers */
#include<Eigen/Dense>
#include<cmath>

namespace Serenity {
/**
 * @brief Calculate an arbitrary function of a symmetric matrix.
 * @param matrix The matrix.
 * @param f The function.
 * @return f(matrix).
 */
inline Eigen::MatrixXd mFunc_Sym(const Eigen::MatrixXd& matrix, std::function<double(const double&)> f) {
  assert(matrix.cols()==matrix.rows());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(
      matrix,
      Eigen::DecompositionOptions::ComputeEigenvectors);
  Eigen::VectorXd eigenvals = eigensolver.eigenvalues().eval();
  const Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors().eval();
  for (unsigned int col=0; col< matrix.cols(); ++col) {
    eigenvals(col) = f(eigenvals(col));
  }
  return (eigenvectors * eigenvals.asDiagonal() * eigenvectors.transpose()).eval();
}

/**
 * @brief Calculates the square root of the given matrix.
 * @param matrix The matrix.
 * @return sqrt(matrix).
 */
inline Eigen::MatrixXd mSqrt_Sym(const Eigen::MatrixXd& matrix) {
  auto const func = [&] (const double& e) {
    if(e < 0.0)throw SerenityError("You are trying to calculate the square root of a negative number!");
    return sqrt(e);
  };// func
  return mFunc_Sym(matrix,func);
}

/**
 * @brief Calculates the square root of the pseudo invers of the given matrix.
 * @param matrix The matrix.
 * @param threshold Threshold for setting entries in the invers to zero.
 * @return matrix^(-1/2).
 */
inline Eigen::MatrixXd pseudoInversSqrt_Sym(const Eigen::MatrixXd& matrix, double threshold = 1e-6) {
  auto const func = [&] (const double& e) {
    if(e < -1.0)throw SerenityError("Tolerance of negative eigenvalues in the pseudo inverse exceeded! You are trying to calculate the square root of a negative number!");
    return (e >= threshold)? 1.0/sqrt(e) : 0.0;
  };// func
  return mFunc_Sym(matrix,func);
}
} /* namespace Serenity */
#endif /* MATH_LINEARALGEBRA_MATRIXFUNCTIONS_H_ */