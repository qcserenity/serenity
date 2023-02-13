/**
 * @file   CholeskyDecomposer.cpp
 *
 * @date   Jun 28, 2018
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
#include "integrals/decomposer/CholeskyDecomposer.h"
/* Include Serenity Internal Headers */
#include "integrals/CDStorageController.h"
#include "integrals/decomposer/SimpleCholeskyDecomposer.h"
#include "io/FormattedOutputStream.h" //Filtered output streams.
#include "memory/MemoryManager.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Serenity {

CholeskyDecomposer::CholeskyDecomposer(
    std::string matrixName, std::shared_ptr<CDIntegralController> cdIntegralController, Eigen::VectorXd& diagonal,
    std::function<std::unique_ptr<Eigen::MatrixXd>(std::vector<int>& qualIndices, std::vector<int>& redIndices)> matrixColumnCalculator,
    const double decompositionThreshold, const double screeningDamping, const double spanFactor,
    const unsigned int nMaxQual, const double negativeThreshold, const double negativeFailThreshold)
  : _matrixName(matrixName),
    _decomposed(false),
    _cdIntegralController(cdIntegralController),
    _diagonal(diagonal),
    _matrixColumnCalculator(matrixColumnCalculator),
    _decompositionThreshold(decompositionThreshold),
    _screeningDamping(screeningDamping),
    _spanFactor(spanFactor),
    _nMaxQual(nMaxQual),
    _negativeThreshold(negativeThreshold),
    _negativeFailThreshold(negativeFailThreshold),
    _maxJ(0) {
  assert(screeningDamping >= 1.0);
  assert(_cdIntegralController);
  try {
    _cdIntegralController->getStorageController(_matrixName)->storeDiag(std::make_shared<Eigen::VectorXd>(_diagonal));
  }
  catch (...) {
  }
}

void CholeskyDecomposer::run() {
  OutputControl::vnOut << " Calculating (screened) Cholesky vectors for \"" + _matrixName + "\":" << std::endl;
  OutputControl::vnOut << " ---------------------------------------------------" << std::endl;
  OutputControl::vnOut << "  Decomposition Progress (maxDiag/cdThreshold): \n" << std::endl;

  Timings::takeTime("Math -  Cholesky decomposition");
  decompose(_matrixName);
  Timings::timeTaken("Math -  Cholesky decomposition");

  OutputControl::vnOut << "    Number of calculated Cholesky vectors: " << _maxJ << std::endl;

  OutputControl::vnOut << " ---------------------------------------------------" << std::endl << std::endl;
}

void CholeskyDecomposer::decompose(std::string label) {
  //   -----------------------------------------------------------------------------------------------
  //   Initialization
  //   -----------------------------------------------------------------------------------------------

  // Get storage Controller and obtain unmodified diagonal (of the matrix being decomposed)
  auto storageController = _cdIntegralController->getStorageController(label);
  _diagonal = (*storageController->loadDiag());
  _cBasis.clear();

  // Find dMax
  double dMax = _diagonal.maxCoeff();

  // 2. calc reduced set savings
  std::vector<int> redIndices(_diagonal.rows(), -1);
  std::vector<std::pair<double, unsigned int>> pairVector;
  std::vector<unsigned int> qualifiedSet;
  std::vector<int> qualIndices(_diagonal.rows(), -1);

  // Initialize vector counter n
  unsigned int n = 0;

  // While Dmax^(i+1) > tau
  do {
    takeTime("cycle");
    //   -----------------------------------------------------------------------------------------------
    //   Reduced Set
    //   -----------------------------------------------------------------------------------------------

    // 2. and 5.k Compute reduced set of significant diagonal elements
    ////compute reducedSet and reducedIndices for each thread
    std::vector<std::vector<unsigned int>> reducedSetSplit;
    std::vector<std::vector<int>> redIndicesSplit;
    std::vector<unsigned int> tmp;
    std::vector<int> tmp2(_diagonal.rows(), -1);
#ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    for (unsigned int i = 0; i < nThreads; i++) {
      reducedSetSplit.push_back(tmp);
      redIndicesSplit.push_back(tmp2);
    }
#else
    const unsigned int nThreads = 1;
    reducedSetSplit.push_back(tmp);
    redIndicesSplit.push_back(tmp2);
#endif
    unsigned int unassignedThreads = nThreads;

    std::vector<int> redIndices(_diagonal.rows(), -1);

    int reducedSetSize = 0;

    // Calculate the reduced set depending on the decomposition mode
    for (unsigned int i = 0; i < _diagonal.size(); i++) {
      if (_screeningDamping * sqrt(_diagonal[i] * dMax) > _decompositionThreshold) {
        reducedSetSize += 1;
      }
    }

    // Distribute the reduced set on available threads
    int diagDistrCounter = reducedSetSize;
    int k = 0;
    int l = 0;
    for (unsigned int threadId = 0; threadId < nThreads; threadId++) {
      int upper = ceil((double)diagDistrCounter / (double)unassignedThreads);
      int m = 0;
      do {
        if (_screeningDamping * sqrt(_diagonal[k] * dMax) > _decompositionThreshold) {
          reducedSetSplit[threadId].push_back(k);
          redIndicesSplit[threadId][k] = reducedSetSplit[threadId].size() - 1;
          redIndices[k] = l;
          m++;
          l++;
        }
        k++;
      } while (m < upper);

      unassignedThreads -= 1;
      diagDistrCounter -= upper;
      if (diagDistrCounter == 0)
        break;
    }

    //   -----------------------------------------------------------------------------------------------
    //   Qualified Set (and minimal diagonal element)
    //   -----------------------------------------------------------------------------------------------

    // Compute smallest diagonal that may be treated
    double dMin = std::max(_spanFactor * dMax, _decompositionThreshold);

    // Compute set of qualified diagonals
    // create pair-vector to sort diagonal elements and
    // retain original indices, i.e. the pivoting indices
    pairVector.clear();
    qualifiedSet.clear();
    for (unsigned int i = 0; i < qualIndices.size(); i++) {
      qualIndices[i] = -1;
    }

    for (unsigned int threadId = 0; threadId < nThreads; threadId++) {
      for (unsigned int i = 0; i < reducedSetSplit[threadId].size(); ++i) {
        pairVector.push_back(std::make_pair(_diagonal(reducedSetSplit[threadId][i]), reducedSetSplit[threadId][i]));
      }
    }

    // sort in descending order
    std::stable_sort(pairVector.rbegin(), pairVector.rend());
    // calculate qualified set
    for (unsigned int i = 0; i < pairVector.size(); ++i) {
      if (i >= _nMaxQual)
        break;
      if (pairVector[i].first > dMin) {
        qualifiedSet.push_back(pairVector[i].second);
        qualIndices[pairVector[i].second] = qualifiedSet.size() - 1;
      }
      else {
        qualIndices[pairVector[i].second] = -1;
      }
    }

    //   -----------------------------------------------------------------------------------------------
    //   Perform decomposition of the qual x qual matrix block to determine partial Cholesky basis
    //   -----------------------------------------------------------------------------------------------
    takeTime("CholsBasis");
    std::vector<unsigned int> partialCholBasis;
    if (qualifiedSet.size() == 0)
      throw SerenityError("Decomposition with insufficient accuracy chosen. No vectors remain to reconstruct the "
                          "matrix. Try choosing a smaller value for cdThreshold.");
    {
      // Calculate all integrals in the qual x qual block
      std::shared_ptr<Eigen::MatrixXd> qqCols(std::shared_ptr<Eigen::MatrixXd>(nullptr));
      qqCols = (_matrixColumnCalculator(qualIndices, qualIndices));
      // Subtract contributions of previous Cholesky vectors
      // Note that for this step all previous vectors have to be loaded resulting
      // in double the load times in disk mode. However, this is compensated by the
      // reduction of computation time.
      if (n > 0) {
        unsigned int countJ = 0;
        while (countJ < _maxJ) {
          unsigned int batchSize = storageController->loadBatch(countJ);
          std::shared_ptr<Eigen::RowVectorXd> cholVec;
          auto reducedVec = std::make_shared<Eigen::VectorXd>(qualifiedSet.size());
          auto qualVec = std::make_shared<Eigen::VectorXd>(qualifiedSet.size());
          unsigned int qualSize = qualifiedSet.size();

          unsigned int redSize = qualSize;

          for (unsigned int J = countJ; J < countJ + batchSize; J++) {
            cholVec = storageController->loadVector(J);

            for (unsigned int p = 0; p < redSize; p++) {
              (*reducedVec)[p] = (*cholVec)[qualifiedSet[p]];
            }
            for (unsigned int q = 0; q < qualSize; q++) {
              (*qualVec)[q] = (*cholVec)[qualifiedSet[q]];
            }
            auto qualData = qualVec->data();
            auto colData = qqCols->data();

            for (unsigned int q = 0; q < qualSize; q++, qualData++) {
              auto redData = reducedVec->data();
              for (unsigned int p = 0; p < redSize; p++, colData++, redData++) {
                (*colData) -= (*redData) * (*qualData);
              }
            }
          }
          reducedVec.reset();
          qualVec.reset();
          cholVec.reset();

          storageController->freeLastBatch();
          countJ += batchSize;
        }
      }

      // Reduce diagonal to the qual x qual block
      Eigen::VectorXd qualDiag(qualifiedSet.size());

      // calculate maximum diagonal element in the qualified diagonal
      qualDiag = qqCols->diagonal();
      double qMax;
      unsigned int qMaxIndex;
      qMax = qualDiag(0);
      qMaxIndex = 0;
      for (unsigned int q = 1; q < qualDiag.size(); ++q) {
        if (qualDiag(q) > qMax) {
          qMaxIndex = q;
          qMax = qualDiag(q);
        }
      }

      // Calculate Cholesky vectors to determine the Cholesky basis for the qual
      // x qual block. The vectors themself are discarded after their contribution
      // is subtracted from the remainder matrix.
      unsigned int j = 0;
      do {
        partialCholBasis.push_back(qualifiedSet[qMaxIndex]);
        j += 1;

        double qMaxFactor = 1.0 / std::sqrt(qMax);

        Eigen::VectorXd cholVec = qqCols->col(qMaxIndex) * qMaxFactor;

        (*qqCols) -= cholVec * cholVec.transpose();

        qualDiag -= cholVec.cwiseProduct(cholVec);

        // set treated diagonal element to zero
        qualDiag(qMaxIndex) = 0.0;
        (*qqCols)(qMaxIndex, qMaxIndex) = 0.0;

        // set diagonal elements below threshold to zero
        for (unsigned int p = 0; p < qualDiag.size(); ++p) {
          if (qualDiag(p) < _negativeThreshold) {
            if (qualDiag(p) < _negativeFailThreshold) {
              std::cout << "\nFailed during partial decomposition at element " << p << " with a value of "
                        << qualDiag(p) << " n: " << n << std::endl
                        << "The label is: " << _matrixName << std::endl;
              throw SerenityError("CholeskyDecomposer: Matrix not positive definite");
            }
            else {
              qualDiag(p) = 0.0;
            }
          }
        }

        qMax = qualDiag(0);
        qMaxIndex = 0;
        for (unsigned int q = 1; q < qualDiag.size(); ++q) {
          if (qualDiag(q) > qMax) {
            qMaxIndex = q;
            qMax = qualDiag(q);
          }
        }

      } while (j < qualifiedSet.size() and qMax > dMin);

    } // end of cholesky basis calculation.

    for (unsigned int i = 0; i < partialCholBasis.size(); i++) {
      _cBasis.push_back(partialCholBasis[i]);
    }

    qualifiedSet = partialCholBasis;
    for (unsigned int i = 0; i < qualIndices.size(); i++) {
      qualIndices[i] = -1;
    }
    for (unsigned int i = 0; i < qualifiedSet.size(); i++) {
      qualIndices[qualifiedSet[i]] = i;
    }

    timeTaken(3, "CholsBasis");

    //   -----------------------------------------------------------------------------------------------
    //   Calculate matrix elements corresponding to diagonal elements in the qualified set
    //   -----------------------------------------------------------------------------------------------

    takeTime("colCalc");

    // Get matrix columns corresponding to qualified diagonals
    std::vector<std::shared_ptr<Eigen::MatrixXd>> qualifiedColsSplit;
    for (unsigned int threadId = 0; threadId < nThreads; threadId++) {
      qualifiedColsSplit.push_back(std::shared_ptr<Eigen::MatrixXd>(nullptr));
      qualifiedColsSplit[threadId] = (_matrixColumnCalculator(qualIndices, redIndicesSplit[threadId]));
    }

    timeTaken(3, "colCalc");

    //   -----------------------------------------------------------------------------------------------
    //   Subtract contribution of previous vectors from the integral columns.
    //   -----------------------------------------------------------------------------------------------

    // Subtract contributions from previous vectors (if any)
    takeTime("subtraction");

    if (n > 0) {
      unsigned int countJ = 0;
      while (countJ < _maxJ) {
        unsigned int batchSize = storageController->loadBatch(countJ);
        Eigen::setNbThreads(1);
#pragma omp parallel
        {
          const unsigned int threadId = omp_get_thread_num();
          std::shared_ptr<Eigen::RowVectorXd> cholVec;
          auto reducedVec = std::make_shared<Eigen::VectorXd>(reducedSetSplit[threadId].size());
          auto qualVec = std::make_shared<Eigen::VectorXd>(qualifiedSet.size());
          unsigned int qualSize = qualifiedSet.size();

          unsigned int redSize = reducedSetSplit[threadId].size();

          for (unsigned int J = countJ; J < countJ + batchSize; J++) {
            cholVec = storageController->loadVector(J);

            for (unsigned int p = 0; p < redSize; p++) {
              (*reducedVec)[p] = (*cholVec)[reducedSetSplit[threadId][p]];
            }
            for (unsigned int q = 0; q < qualSize; q++) {
              (*qualVec)[q] = (*cholVec)[qualifiedSet[q]];
            }

            auto qualData = qualVec->data();
            auto colData = qualifiedColsSplit[threadId]->data();

            for (unsigned int q = 0; q < qualSize; q++, qualData++) {
              auto redData = reducedVec->data();
              for (unsigned int p = 0; p < redSize; p++, colData++, redData++) {
                (*colData) -= (*redData) * (*qualData);
              }
            }
          }
          reducedVec.reset();
          qualVec.reset();
          cholVec.reset();
        } // end of parallel section
        Eigen::setNbThreads(0);
        storageController->freeLastBatch();
        countJ += batchSize;
      }
    }

    timeTaken(3, "subtraction");

    //   -----------------------------------------------------------------------------------------------
    //   Calculation of the Cholesky vectors
    //   -----------------------------------------------------------------------------------------------

    takeTime("EvalChol");
    // Compute largest diagonal among qualified
    double qMax;
    unsigned int qMaxIndex;
    findQmax(qMax, qMaxIndex, qualifiedSet);

    // Initialize Counter for cholesky vectors
    unsigned int j = 0;

    // While j , dim(D^(i)) and Qmax > Dmin^(i)
    do {
      // Update counter j=j+1, J = n + j
      j += 1;
      unsigned int J = n + j - 1;

      double qMaxFactor = 1.0 / std::sqrt(qMax);

#pragma omp parallel
      {
        const unsigned int threadId = omp_get_thread_num();

        // Calculate Cholseky vector (and store in decomposed column of qualified set)
        for (unsigned int p = 0; p < reducedSetSplit[threadId].size(); ++p) {
          if (_diagonal(reducedSetSplit[threadId][p]) <= 0.0) {
            (*qualifiedColsSplit[threadId])(p, qMaxIndex) = 0.0;
          }
          else {
            (*qualifiedColsSplit[threadId])(p, qMaxIndex) *= qMaxFactor;
          }
        }

        // Assign Cholesky basis function h_J to [q]_J, the index corresponding to Q_max
        //[q]_J (qMaxIndex) has been calculated by function findQMax

        // Update

        for (unsigned int p = 0; p < reducedSetSplit[threadId].size(); p++) {
          _diagonal(reducedSetSplit[threadId][p]) -=
              ((*qualifiedColsSplit[threadId])(p, qMaxIndex) * (*qualifiedColsSplit[threadId])(p, qMaxIndex));
          if (_diagonal(reducedSetSplit[threadId][p]) < 0.0) {
          }
        }
      } // end of parallel section

      // set treated diagonal element to zero
      _diagonal(qualifiedSet[qMaxIndex]) = 0.0;

      // set diagonal elements below threshold to zero
      for (unsigned int p = 0; p < _diagonal.rows(); ++p) {
        if (_diagonal(p) < _negativeThreshold) {
          if (_diagonal(p) < _negativeFailThreshold) {
            std::cout << "\nFailed at element " << p << " with a value of " << _diagonal(p) << " n: " << n << std::endl;
            throw SerenityError("CholeskyDecomposer: Matrix not positive definite");
          }
          else {
            _diagonal(p) = 0.0;
          }
        }
      }

      // Obtain full reduced Cholesky vector by combining results from all threads.
      Eigen::VectorXd qMaxVector(reducedSetSize);
      qMaxVector.setZero();
      for (unsigned int i = 0, j = 0; i < nThreads; i++) {
        for (unsigned int k = 0; k < (*qualifiedColsSplit[i]).rows(); k++) {
          qMaxVector[j] += (*qualifiedColsSplit[i])(k, qMaxIndex);
          j++;
        }
      }

      auto tempCholVec = std::make_shared<Eigen::RowVectorXd>(_diagonal.size());
      tempCholVec->setZero();
#pragma omp parallel
      {
        const unsigned int threadId = omp_get_thread_num();

        for (unsigned int q = 0; q < qualifiedSet.size(); q++) {
          if (_diagonal(qualifiedSet[q]) <= 0.0)
            continue;
          if (redIndices[qualifiedSet[q]] > -1) {
            double qFactor = qMaxVector[redIndices[qualifiedSet[q]]];
            for (unsigned int p = 0; p < reducedSetSplit[threadId].size(); p++) {
              (*qualifiedColsSplit[threadId])(p, q) -= qFactor * (*qualifiedColsSplit[threadId])(p, qMaxIndex);
            }
          }
        }

        for (unsigned int p = 0; p < reducedSetSplit[threadId].size(); p++) {
          (*tempCholVec)(reducedSetSplit[threadId][p]) = (*qualifiedColsSplit[threadId])(p, qMaxIndex);
        }

      } // end of parallel section
      storageController->storeVector(n + j - 1, tempCholVec);

      findQmax(qMax, qMaxIndex, qualifiedSet);

      _maxJ = J + 1;

    } while (j < qualifiedSet.size() and qMax > dMin); /*5.h. While j < dim(Q^(i)) and Qmax > Dmin^(i) */

    timeTaken(3, "EvalChol");
    //   -----------------------------------------------------------------------------------------------
    //   Updates and clean-up
    //   -----------------------------------------------------------------------------------------------

    // free Memory from Cholesky vectors stored in integral Matrix
    for (unsigned int threadId = 0; threadId < qualifiedColsSplit.size(); threadId++) {
      qualifiedColsSplit[threadId].reset();
    }

    // Update vector counter n = n + j
    n += j;

    // 5.j. Compute largest diagonal
    dMax = 0.0;
    for (unsigned int threadId = 0; threadId < nThreads; threadId++) {
      for (unsigned int p = 0; p < reducedSetSplit[threadId].size(); p++) {
        if (_diagonal[reducedSetSplit[threadId][p]] > dMax) {
          dMax = _diagonal[reducedSetSplit[threadId][p]];
        }
      }
    }

    // Compute reduces set of significant diagonal elements
    // reduced set is calculated at the beginning of the loop

    // print decomposition progress
    std::ios_base::fmtflags f(std::cout.flags());
    OutputControl::vnOut << "    " << std::scientific << std::setprecision(6) << dMax << " / " << std::setprecision(0)
                         << _decompositionThreshold << std::endl;
    std::cout.flags(f);
    if (GLOBAL_PRINT_LEVEL == Options::GLOBAL_PRINT_LEVELS::VERBOSE or
        GLOBAL_PRINT_LEVEL == Options::GLOBAL_PRINT_LEVELS::DEBUGGING) {
      timeTaken(0, "cycle");
    }
    else {
      timeTaken(3, "cycle");
    }

  } while (dMax > _decompositionThreshold); /* 5. While Dmax^(i+1) > tau */

  _decomposed = true;
}

std::vector<unsigned int> CholeskyDecomposer::getCholeskyBasis() {
  if (!_decomposed)
    run();
  return _cBasis;
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
