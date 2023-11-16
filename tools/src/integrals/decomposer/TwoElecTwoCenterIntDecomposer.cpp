/**
 * @file   TwoElecTwoCenterIntDecomposer.cpp
 *
 * @date   Jun 6, 2019
 * @author Lars Hellmann
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
#include "integrals/decomposer/TwoElecTwoCenterIntDecomposer.h"
/* Include Serenity Internal Headers */
#include "integrals/CDStorageController.h"
#include "integrals/decomposer/CholeskyDecomposer.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/wrappers/Libint.h"

namespace Serenity {

TwoElecTwoCenterIntDecomposer::TwoElecTwoCenterIntDecomposer(const Settings& settings,
                                                             std::shared_ptr<BasisController> basController,
                                                             std::shared_ptr<CDIntegralController> cdIntController,
                                                             std::string label, LIBINT_OPERATOR op, double mu)
  : _settings(settings),
    _basisController(basController),
    _cdIntController(cdIntController),
    _label(label),
    _decomposer(nullptr),
    _cdThresh(_settings.basis.cdThreshold),
    _op(op),
    _mu(mu) {
}

void TwoElecTwoCenterIntDecomposer::run() {
  unsigned int nbfs = _basisController->getNBasisFunctions();
  auto storageController = _cdIntController->getStorageController(_label);
  double decompositionThreshold = _cdThresh;

  bool normAux = !(_basisController->isAtomicCholesky());

  if (storageController->getNVectors())
    return;

  // Calc Two Center Diagonal
  Eigen::VectorXd diagonal(nbfs);
  diagonal.setZero();

  auto& libint = Libint::getInstance();
  libint.finalize(_op, 0, 2);
  libint.initialize(_op, 0, 2, std::vector<std::shared_ptr<Atom>>(0), _mu);
  Eigen::MatrixXd ints;
  auto shells = _basisController->getBasis();
  unsigned int nShells = shells.size();

  for (unsigned int I = 0; I < nShells; ++I) {
    auto shellI = *shells[I];
    if (shellI.getAngularMomentum() > AM_MAX)
      continue;
    bool significant = libint.compute(_op, 0, shellI, shellI, ints, normAux);
    if (!significant)
      continue;
    auto nI = shellI.getNContracted();
    auto firstI = _basisController->extendedIndex(I);
    for (unsigned int i = firstI, ii = 0; i < firstI + nI; i++, ii++) {
      diagonal[i] = ints.col(0)[ii * nI + ii];
    }
  }

  storageController->storeDiag(std::make_shared<Eigen::VectorXd>(diagonal));

  const auto& basis = _basisController->getBasis();

  /*
   * This function calculates columns of the integral super-matrix M_{ij,kl}
   * using a vector with column indices as input. The function works similarly
   * as the loop function in TwoElecFourCenterIntLooper.h. The column indices
   * are unsorted and the resulting integral matrix is strongly rectangular
   * which makes it difficult to fully exploit symmetry.
   *
   * For documentation on the general lambda function see CholeskyDecomposer.h
   */

  auto columnCalculator = [&](const std::vector<int>& qualIndices, const std::vector<int>& redIndices) {
    unsigned int redSize = 0;
    unsigned int qualSize = 0;
    for (int q : qualIndices)
      if (q > -1)
        qualSize++;
    for (int p : redIndices)
      if (p > -1)
        redSize++;

    Eigen::MatrixXd mat1 = Eigen::MatrixXd::Zero(redSize, qualSize);
    std::vector<std::vector<double>> mat2(qualSize, std::vector<double>(redSize, 0.0));

    // Get list of qualified shell-pairs, i.e. shells containing the index i
    std::vector<int> qualifiedShells;
    for (unsigned int k = 0; k < nShells; ++k) {
      const unsigned int firstK = _basisController->extendedIndex(k);
      const unsigned int nK = basis[k]->getNContracted();
      for (unsigned int kk = firstK; kk < nK + firstK; ++kk) {
        if (qualIndices[kk] > -1) {
          qualifiedShells.push_back(k);
          break;
        }
      }
    }

    // Get reduced list of shell-pairs, i.e. shells containing the index i
    std::vector<int> reducedShells;
    for (unsigned int k = 0; k < nShells; ++k) {
      const unsigned int firstK = _basisController->extendedIndex(k);
      const unsigned int nK = basis[k]->getNContracted();
      for (unsigned int kk = firstK; kk < nK + firstK; ++kk) {
        if (redIndices[kk] > -1) {
          reducedShells.push_back(k);
          break;
        }
      }
    }

#pragma omp parallel
    {
      Eigen::MatrixXd ints;

      double screeningThresh = _settings.basis.integralThreshold;
      if (decompositionThreshold * 0.1 < screeningThresh)
        screeningThresh = decompositionThreshold * 0.1;
#pragma omp for schedule(static, 1)
      for (unsigned int I = 0; I < reducedShells.size(); I++) {
        const unsigned int i = reducedShells[I];
        const auto& basI = *basis[i];
        const unsigned int nI = basis[i]->getNContracted();
        const unsigned int firstI = _basisController->extendedIndex(i);
        for (unsigned int J = 0; J < qualifiedShells.size(); J++) {
          const unsigned int j = qualifiedShells[J];
          const auto& basJ = *basis[j];
          const unsigned int nJ = basis[j]->getNContracted();
          const unsigned int firstJ = _basisController->extendedIndex(j);
          // calculate integrals
          bool significant = libint.compute(_op, 0, basI, basJ, ints, normAux);
          if (significant == false) {
            continue;
          }

          // try to use other integrals from this set of integrals
          for (unsigned int ii = 0; ii < nI; ii++) {
            const unsigned int iii = ii + firstI;
            int ri = redIndices[iii];
            if (ri < 0)
              continue;
            for (unsigned int jj = 0; jj < nJ; jj++) {
              const unsigned int jjj = jj + firstJ;
              int qj = qualIndices[jjj];
              if (qj < 0)
                continue;
              double integral = ints.col(0)(jj * nI + ii);
              // Eigen is col major -> second vec index is rows; first is cols
              mat2[qj][ri] = integral;
            }
          }
        } /* loop over qualified shell pairs kl */
      }   /* loop over shell pairs ij */
    }

    //    Map vec<vec> to Eigen::Matrix
    Eigen::MatrixXd mat3 = Eigen::MatrixXd::Zero(redSize, qualSize);

    for (unsigned int v = 0; v < mat2[0].size(); v++) {
      for (unsigned int w = 0; w < mat2.size(); w++) {
        mat3(v, w) = mat2[w][v];
      }
    }

    std::unique_ptr<Eigen::MatrixXd> mat = std::make_unique<Eigen::MatrixXd>(mat3);

    return mat;
  };

  // Input parameters according to Aquilante et. al.
  double screeningDamping = 1.0;
  if (decompositionThreshold > 1.0e-8)
    screeningDamping = 1.0e+9 * decompositionThreshold;

  // Build CholeskyDecomposer and run the decomposition
  _decomposer = std::make_shared<CholeskyDecomposer>(_label, _cdIntController, diagonal, columnCalculator,
                                                     decompositionThreshold, screeningDamping, 0.01, 500);

  _decomposer->run();

  libint.finalize(_op, 0, 2);

  return;
}

std::vector<unsigned int> TwoElecTwoCenterIntDecomposer::getCholeskyBasis() {
  if (!_decomposer)
    this->run();
  return _decomposer->getCholeskyBasis();
}

void TwoElecTwoCenterIntDecomposer::setThreshold(double cdThresh) {
  _cdThresh = cdThresh;
}

} /* namespace Serenity */
