/**
 * @file   TwoElecFourCenterIntDecomposer.cpp
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
#include "integrals/decomposer/TwoElecFourCenterIntDecomposer.h"
/* Include Serenity Internal Headers */
#include "integrals/CDStorageController.h"
#include "integrals/decomposer/CholeskyDecomposer.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/wrappers/Libint.h"

namespace Serenity {

TwoElecFourCenterIntDecomposer::TwoElecFourCenterIntDecomposer(const Settings& settings,
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

void TwoElecFourCenterIntDecomposer::run() {
  unsigned int nbfs = _basisController->getNBasisFunctions();
  auto storageController = _cdIntController->getStorageController(_label);
  double decompositionThreshold = _cdThresh;

  if (storageController->getNVectors())
    return;

  // Calculate and store the diagonal of the matrix that will be decomposed
  Eigen::VectorXd diagonal(nbfs * nbfs);
  diagonal.setZero();
  TwoElecFourCenterIntLooper looper(_op, 0, _basisController, 1.0e-10, _mu);

  auto const calcDiagonal = [&](const unsigned int& i, const unsigned int& j, const Eigen::VectorXd& integral,
                                const unsigned int threadId) {
    (void)threadId; // no warnings, please
    unsigned int ij = i * nbfs + j;
    unsigned int ji = j * nbfs + i;
    diagonal(ij) = integral(0);
    diagonal(ji) = integral(0);
  };

  // Calls a function that only loops over the diagonal elements of the ERI-matrix
  looper.loopDiagonal(calcDiagonal);

  const auto& shellPairs = _basisController->getShellPairData();
  const auto& basis = _basisController->getBasis();
  auto& libint = Libint::getInstance();
  libint.initialize(_op, 0, 4, std::vector<std::shared_ptr<Atom>>(0), _mu);
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
    for (int i : qualIndices)
      if (i > -1)
        qualSize++;
    for (int i : redIndices)
      if (i > -1)
        redSize++;

    std::unique_ptr<Eigen::MatrixXd> mat = std::make_unique<Eigen::MatrixXd>(redSize, qualSize);
    (*mat).setZero();

    // Get list of qualified shell-pairs, i.e. shells containing the combined indices ij
    std::vector<int> qualifiedShellPairs;
    for (unsigned int klPair = 0; klPair < (*shellPairs).size(); ++klPair) {
      auto& kl = (*shellPairs)[klPair];
      const unsigned int k = kl.bf1;
      const unsigned int l = kl.bf2;
      const unsigned int firstK = _basisController->extendedIndex(k);
      const unsigned int firstL = _basisController->extendedIndex(l);
      const unsigned int nK = basis[k]->getNContracted();
      const unsigned int nL = basis[l]->getNContracted();
      for (unsigned int kk = 0; kk < nK; ++kk) {
        const unsigned int kkk = kk + firstK;
        for (unsigned int ll = 0; ll < nL; ++ll) {
          const unsigned int lll = ll + firstL;
          const unsigned int klIndex = kkk * nbfs + lll;
          if (qualIndices[klIndex] > -1)
            goto isInQualifiedSet;
          const unsigned int lkIndex = lll * nbfs + kkk;
          if (qualIndices[lkIndex] > -1)
            goto isInQualifiedSet;
        }
      }
      if (false) {
      isInQualifiedSet:
        qualifiedShellPairs.push_back(klPair);
      }
    }

    // Get reduced list of shell-pairs, i.e. shells containing the combined indices ij
    std::vector<int> reducedShellPairs;
    for (unsigned int klPair = 0; klPair < (*shellPairs).size(); ++klPair) {
      auto& kl = (*shellPairs)[klPair];
      const unsigned int k = kl.bf1;
      const unsigned int l = kl.bf2;
      const unsigned int firstK = _basisController->extendedIndex(k);
      const unsigned int firstL = _basisController->extendedIndex(l);
      const unsigned int nK = basis[k]->getNContracted();
      const unsigned int nL = basis[l]->getNContracted();
      for (unsigned int kk = 0; kk < nK; ++kk) {
        const unsigned int kkk = kk + firstK;
        for (unsigned int ll = 0; ll < nL; ++ll) {
          const unsigned int lll = ll + firstL;
          const unsigned int klIndex = kkk * nbfs + lll;
          if (redIndices[klIndex] > -1)
            goto isInReducedSet;
          const unsigned int lkIndex = lll * nbfs + kkk;
          if (redIndices[lkIndex] > -1)
            goto isInReducedSet;
        }
      }
      if (false) {
      isInReducedSet:
        reducedShellPairs.push_back(klPair);
      }
    }

#pragma omp parallel
    {
      Eigen::MatrixXd ints;

      double screeningThresh = _settings.basis.integralThreshold;
      if (decompositionThreshold * 0.1 < screeningThresh)
        screeningThresh = decompositionThreshold * 0.1;
#pragma omp for schedule(static, 1)
      for (int ijPair = reducedShellPairs.size() - 1; ijPair >= 0; --ijPair) {
        auto& ij = (*shellPairs)[reducedShellPairs[ijPair]];
        const unsigned int i = ij.bf1;
        const unsigned int j = ij.bf2;
        const auto& basI = *basis[i];
        const auto& basJ = *basis[j];
        const unsigned int nI = basis[i]->getNContracted();
        const unsigned int nJ = basis[j]->getNContracted();
        const unsigned int firstI = _basisController->extendedIndex(i);
        const unsigned int firstJ = _basisController->extendedIndex(j);
        for (unsigned int klPair = 0; klPair < qualifiedShellPairs.size(); ++klPair) {
          auto& kl = (*shellPairs)[qualifiedShellPairs[klPair]];
          // Prescreening
          if (ij.factor * kl.factor < screeningThresh) {
            klPair = qualifiedShellPairs.size();
            continue;
          }
          const unsigned int k = kl.bf1;
          const unsigned int l = kl.bf2;
          const auto& basK = *basis[k];
          const auto& basL = *basis[l];
          const unsigned int nK = basis[k]->getNContracted();
          const unsigned int nL = basis[l]->getNContracted();
          const unsigned int firstK = _basisController->extendedIndex(k);
          const unsigned int firstL = _basisController->extendedIndex(l);
          // calculate integrals
          bool significant = libint.compute(_op, 0, basI, basJ, basK, basL, ints);
          if (!significant)
            continue;
          // try to use other integrals from this set of integrals
          for (unsigned int ii = 0; ii < nI; ++ii) {
            const unsigned int iii = ii + firstI;
            for (unsigned int jj = 0; jj < nJ; ++jj) {
              const unsigned int jjj = jj + firstJ;
              if (jjj > iii)
                continue;
              const unsigned int iijj = iii * nbfs + jjj;
              const unsigned int jjii = jjj * nbfs + iii;
              // get index in reduced set (if member of)
              int riijj = redIndices[iijj];
              int rjjii = redIndices[jjii];
              if (riijj < 0 and rjjii < 0)
                continue;
              for (unsigned int kk = 0; kk < nK; ++kk) {
                const unsigned int kkk = kk + firstK;
                for (unsigned int ll = 0; ll < nL; ++ll) {
                  const unsigned int lll = ll + firstL;
                  if (lll > kkk)
                    continue;
                  const unsigned int kkll = kkk * nbfs + lll;
                  const unsigned int llkk = lll * nbfs + kkk;
                  // get index in qualified set (if member of)
                  int qkkll = qualIndices[kkll];
                  int qllkk = qualIndices[llkk];
                  double integral = ints.col(0)(ii * nJ * nK * nL + jj * nK * nL + kk * nL + ll);
                  if (qkkll > -1) {
                    if (riijj > -1)
                      (*mat)(riijj, qkkll) = integral;
                    if (rjjii > -1)
                      (*mat)(rjjii, qkkll) = integral;
                  }
                  if (qllkk > -1) {
                    if (riijj > -1)
                      (*mat)(riijj, qllkk) = integral;
                    if (rjjii > -1)
                      (*mat)(rjjii, qllkk) = integral;
                  }
                } /* loop over ll */
              }   /* loop over kk */
            }     /* loop over jj */
          }       /* loop over ii */
        }         /* loop over qualified shell pairs kl */
      }           /* loop over shell pairs ij */
    }

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

  libint.finalize(_op, 0, 4);

  return;
}

std::vector<unsigned int> TwoElecFourCenterIntDecomposer::getCholeskyBasis() {
  if (!_decomposer)
    this->run();
  return _decomposer->getCholeskyBasis();
}

void TwoElecFourCenterIntDecomposer::setThreshold(double cdThresh) {
  _cdThresh = cdThresh;
}

} /* namespace Serenity */
