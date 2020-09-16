/**
 * @file SphericalHarmonicsRotations.cpp
 *
 * @author Moritz Bensberg
 * @date Feb 11, 2020
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
#include "basis/SphericalHarmonicsRotations.h"
/* Include Std and External Headers */
#include <math.h> //sin, cos, pow and sqrt

namespace Serenity {

/*
 * Map for the transformation matrices J.
 */
std::map<unsigned int, std::shared_ptr<Eigen::MatrixXd>> SphericalHarmonicsRotations::_jMatrices = {
    {0, std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Identity(1, 1))},
    {1, nullptr},
    {2, nullptr},
    {3, nullptr},
    {4, nullptr},
    {5, nullptr},
    {6, nullptr},
    {7, nullptr},
    {8, nullptr},
    {9, nullptr},
    {10, nullptr},
    {11, nullptr},
    {12, nullptr},
    {13, nullptr},
    {14, nullptr},
    {15, nullptr},
    {16, nullptr},
    {17, nullptr},
    {18, nullptr},
    {19, nullptr},
    {20, nullptr}};

Eigen::MatrixXd SphericalHarmonicsRotations::calculateHatGZMatrix(unsigned int l) {
  unsigned int l2Plus1 = 2 * l + 1;
  Eigen::MatrixXd gL = Eigen::MatrixXd::Zero(l2Plus1, l2Plus1);
  double denominator = sqrt((2 * l + 1) * (2 * l + 3));
  for (unsigned int k = 1; k <= l2Plus1; ++k) {
    gL(k - 1, k - 1) = sqrt(k * (2 * l + 2 - k)) / denominator;
  }          // for k
  return gL; // 2(l+1)+1x2l+1;
}

Eigen::MatrixXd SphericalHarmonicsRotations::calculateGYMatrix(unsigned int l) {
  unsigned int l2Plus1 = 2 * l + 1;
  Eigen::MatrixXd gL = Eigen::MatrixXd::Zero(2 * (l + 1) + 1, l2Plus1);
  double denominator = 2 * sqrt((2 * l + 1) * (2 * l + 3));
  for (unsigned int k = 1; k <= l - 1; ++k) {
    double coeff = sqrt(k * (k + 1)) / denominator;
    gL(2 * l + 2 - k - 1, k - 1) = coeff;
    gL(2 + k - 1, 2 * l + 2 - k - 1) = -coeff;
  } // for k
  for (unsigned int k = 1; k <= l; ++k) {
    double coeff = -sqrt((2 * l + 2 - k) * (2 * l + 3 - k)) / denominator;
    gL(k - 1, 2 * l + 2 - k - 1) = coeff;
    gL(2 * l + 4 - k - 1, k - 1) = -coeff;
  } // for k
  gL(l + 2 - 1, l - 1) = sqrt(2 * l * (l + 1)) / denominator;
  gL(l, l) = -sqrt(2 * (l + 1) * (l + 2)) / denominator;
  return gL;
}

Eigen::MatrixXd SphericalHarmonicsRotations::calculateXMatrix(unsigned int l, double alpha) {
  unsigned int matrixSize = 2 * l + 1;
  Eigen::MatrixXd x = Eigen::MatrixXd::Zero(matrixSize, matrixSize);
  unsigned int counter = 0;
  for (int k = l; k >= (int)-l; --k) {
    x(counter, matrixSize - counter - 1) = sin(k * alpha);
    x(counter, counter) = cos(k * alpha);
    ++counter;
  } // for k
  return x;
}

const Eigen::MatrixXd& SphericalHarmonicsRotations::getJMatrix(unsigned int l) {
  if (!_jMatrices[l]) {
    _jMatrices[l] = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(2 * l + 1, 2 * l + 1));
    Eigen::MatrixXd& jL = *_jMatrices[l];
    if (l == 1) {
      /*
       * The form of this matrix can be deduced from Fig.1 of
       * Journal of Physics A: Mathematical and Theoretical 40, 1597 (2007).
       * j*j has to be an identity matrix, since x(a)*j*x(b)*j*x(c) = delta
       * which holds for any choice of a,b and c. If a=b=c=0, x reduces to the identity
       * matrix and delta has to be the identity as well, since no rotation occurred.
       *
       * With, j(2l,2l) = 2^-(l-1) and j = j^T it follows that j(0,1) and j(1,0) have
       * to have an absolute value of 1. The correct reference for j(l=3) is reproduced
       * for j(0,1)=-1.
       *
       * However, it is probably the easiest to picture this as the exchange between
       * px and py, which if of course the intention of the matrix.
       */
      jL << 0, -1, 0, -1, 0, 0, 0, 0, 1;
      return jL;
    }
    // Set the block according to FIg. 1.
    const Eigen::MatrixXd& j_lMinus1 = getJMatrix(l - 1);
    Eigen::MatrixXd hatGZ_lMinus1 = calculateHatGZMatrix(l - 1);
    Eigen::MatrixXd gY_lMinus1 = calculateGYMatrix(l - 1);
    jL.block(0, 1, 2 * l + 1, 2 * l - 1) = gY_lMinus1 * j_lMinus1 * hatGZ_lMinus1.inverse(); // 2l+1 x 2l+1
    jL.block(1, 0, 2 * l - 1, 1) = jL.block(0, 1, 1, 2 * l - 1).transpose();
    jL.block(1, 2 * l, 2 * l - 1, 1) = jL.block(2 * l, 1, 1, 2 * l - 1).transpose();
    jL(2 * l, 2 * l) = pow(2.0, (int)-(l - 1));
  }
  return *_jMatrices[l];
}

Eigen::MatrixXd SphericalHarmonicsRotations::getTransformationMatrix(unsigned int l, double alpha, double beta, double gamma) {
  const Eigen::MatrixXd& jL = getJMatrix(l);
  Eigen::MatrixXd xAlpha = calculateXMatrix(l, alpha);
  Eigen::MatrixXd xBeta = calculateXMatrix(l, beta);
  Eigen::MatrixXd xGamma = calculateXMatrix(l, gamma);
  return xAlpha * jL * xBeta * jL * xGamma;
}

} /* namespace Serenity */
