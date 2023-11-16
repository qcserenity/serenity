/**
 * @file   Ao2MoHalfTransformer.cpp
 * @author Moritz Bensberg
 * @date   17. December 2018
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
#include "integrals/transformer/Ao2MoHalfTransformer.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"                         //Definition of a basis controller.
#include "integrals/looper/ExchangeInteractionIntLooper.h" //Integral looping.

namespace Serenity {

void Ao2MoHalfTransformer::transformTwoElectronIntegrals(Eigen::MatrixXd& result, const Eigen::MatrixXd& pairDensityMatrix) {
  const unsigned int nBasisB = _basisControllerB->getNBasisFunctions();
  // resize result
  result.resize(nBasisB, nBasisB);
  result.setZero();

  // Build integral looper and calculate the integrals.
#ifdef _OPENMP
  // create a vector of matrices for each thread
  std::vector<Eigen::MatrixXd> eriContr;
  for (int j = 0; j < omp_get_max_threads(); ++j) {
    eriContr.push_back(Eigen::MatrixXd(nBasisB, nBasisB));
    eriContr[j].setZero();
  }
#else
  // or just one
  std::vector<Eigen::MatrixXd> eriContr(1, Eigen::MatrixXd(nBasisB, nBasisB));
  eriContr[0].setZero();
#endif
  ExchangeInteractionIntLooper excLooper(LIBINT_OPERATOR::coulomb, 0, _basisControllerB, _basisControllerA, 1E-10);

  auto const excLooperFunction = [&](const unsigned int& r, const unsigned int& i, const unsigned int& s,
                                     const unsigned int& j, double intValues, const unsigned int& threadID) {
    eriContr[threadID](r, s) += pairDensityMatrix(i, j) * intValues;
  };
  excLooper.loopNoDerivative(excLooperFunction);
#ifdef _OPENMP
  // sum over all threads
  for (int j = 0; j < omp_get_max_threads(); ++j) {
    result += eriContr[j];
  } // for j
#else
  result += eriContr[0];
#endif
}

} /* namespace Serenity */
