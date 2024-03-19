/**
 * @file   CDHFPotential.cpp
 *
 * @date   Oct 26, 2018
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
#include "potentials/CDHFPotential.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/CDIntegralController.h"
#include "integrals/CDStorageController.h"
#include "integrals/decomposer/TwoElecFourCenterIntDecomposer.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CDHFPotential<SCFMode>::CDHFPotential(std::shared_ptr<SystemController> systemController,
                                      std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                                      const double exchangeRatio, const double prescreeningThreshold,
                                      double prescreeningIncrementStart, double prescreeningIncrementEnd)
  : HFPotential<SCFMode>(systemController, dMat, exchangeRatio, prescreeningThreshold, prescreeningIncrementStart,
                         prescreeningIncrementEnd, 0),
    _cdIntController(systemController->getCDIntegralController()) {
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);

  if (systemController->getSettings().basis.densFitJ != Options::DENS_FITS::CD ||
      systemController->getSettings().basis.densFitK != Options::DENS_FITS::CD) {
    throw SerenityError("A CDHFPotential was generated but density-fitting was not set to CD");
  }
  _storageLabel = "AO";
  // Generate vectors here
  TwoElecFourCenterIntDecomposer decomposer(systemController->getSettings(), systemController->getBasisController(),
                                            systemController->getCDIntegralController(), _storageLabel);
  decomposer.run();
};

template<>
void CDHFPotential<Options::SCF_MODES::RESTRICTED>::addToMatrix(FockMatrix<Options::SCF_MODES::RESTRICTED>& F,
                                                                const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densityMatrix) {
  auto systemController = _systemController.lock();
  auto cdStorageController = _cdIntController->getStorageController(_storageLabel);

  // Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fc;
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fx;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  for (unsigned int i = 0; i < nThreads; ++i) {
    fc.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
    fc[i]->setZero();
    fx.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
    fx[i]->setZero();
  }

#else
  fc.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
  fc[0]->setZero();
  fx.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
  fx[0]->setZero();
#endif

  Eigen::VectorXd esSigns;
  unsigned int nOcc;
  Eigen::MatrixXd occCoeff;

  // Check if orbitals are available. If not generate pseudo coefficients from the density matrix
  if (systemController->hasElectronicStructure<RESTRICTED>()) {
    const auto orbController = systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>();
    auto coefficientMatrix = orbController->getCoefficients();
    // Get number of occupied orbitals
    nOcc = systemController->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
    // Coefficient matrix of occupied orbitals
    occCoeff = coefficientMatrix.leftCols(nOcc);
    esSigns = Eigen::VectorXd::Ones(nOcc);
    esSigns *= 2;
  }
  else {
    Eigen::MatrixXd p = densityMatrix;
    auto pseudoCoeff = systemController->getCDIntegralController()->generatePseudoCoefficients(p);
    esSigns = pseudoCoeff.first;
    occCoeff = pseudoCoeff.second;
    nOcc = occCoeff.cols();
  }

  // Get number of basis functions
  const unsigned int n = _basis->getNBasisFunctions();
  // Get number of Cholesky vectors
  unsigned int M = cdStorageController->getNVectors();
  // number of currently loaded vectors
  unsigned int batchSize = 0;

  for (unsigned int J = 0; J < M; J++) {
#pragma omp parallel
    {
      // Determine ThreadId
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif

      // load the current batch of Cholesky vectors
#pragma omp master
      {
        batchSize = cdStorageController->loadBatch(J);
        if (batchSize <= 0) {
          throw SerenityError("Cholesky Vector batch size is less than 1.");
        }
      }
#pragma omp barrier

      Eigen::setNbThreads(1);
#pragma omp for schedule(dynamic)
      // loop over all Cholesky vectors
      for (unsigned int K = J; K < J + batchSize; K++) {
        // load the current vector
        std::shared_ptr<Eigen::RowVectorXd> cholVector = cdStorageController->loadVector(K);
        // Transform Two Index Vector into its Matrix representation
        Eigen::Map<Eigen::MatrixXd> cholVecMat(cholVector->data(), n, n);

        // Calculate Coulomb Contribution
        double densTimesChol = densityMatrix.cwiseProduct(cholVecMat).sum();
        *(fc[threadId]) += densTimesChol * cholVecMat;

        // Calculate Exchange contribution
        if (_xRatio > 0.0) {
          Eigen::MatrixXd half = cholVecMat * occCoeff;
          *(fx[threadId]) -= 0.5 * _xRatio * half * esSigns.asDiagonal() * half.transpose();
        }
      }
      Eigen::setNbThreads(0);

#pragma omp master
      {
        cdStorageController->freeLastBatch();
        J += batchSize - 1;
      }
#pragma omp barrier
    }
  }

#ifdef _OPENMP
  for (unsigned int i = 1; i < nThreads; ++i) {
    *fc[0] += *fc[i];
    *fx[0] += *fx[i];
    delete fx[i];
    delete fc[i];
  }
  *_fullXpotential = *fx[0];
#endif
  F = *fc[0] + *fx[0];
}

template<>
void CDHFPotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F,
                                                                  const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix) {
  auto systemController = _systemController.lock();
  auto cdStorageController = _cdIntController->getStorageController(_storageLabel);

  // Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>*> fc;
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>*> fx;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  // Let the first thread use the Fock matrix directly
  for (unsigned int i = 0; i < nThreads; ++i) {
    fc.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
    fc[i]->alpha.setZero();
    fc[i]->beta.setZero();
    fx.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
    fx[i]->alpha.setZero();
    fx[i]->beta.setZero();
  }

#else
  fc.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
  fc[0]->alpha.setZero();
  fc[0]->beta.setZero();
#endif

  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> occCoeff;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int> nOcc;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd> esSigns;

  if (systemController->hasElectronicStructure<UNRESTRICTED>()) {
    const auto orbController = systemController->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>();
    auto coefficientMatrix = orbController->getCoefficients();
    // Get number of occupied orbitals
    nOcc = systemController->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
    // Coefficient matrix of occupied orbitals
    for_spin(occCoeff, nOcc, esSigns, coefficientMatrix) {
      occCoeff_spin = coefficientMatrix_spin.leftCols(nOcc_spin);
      esSigns_spin = Eigen::VectorXd::Ones(nOcc_spin);
    };
  }
  else {
    for_spin(occCoeff, nOcc, densityMatrix, esSigns) {
      Eigen::MatrixXd p = densityMatrix_spin;
      auto pseudoCoeff = systemController->getCDIntegralController()->generatePseudoCoefficients(p);
      esSigns_spin = pseudoCoeff.first;
      occCoeff_spin = pseudoCoeff.second;
      nOcc_spin = occCoeff_spin.cols();
    };
  }

  // Get number of basis functions
  const unsigned int n = _basis->getNBasisFunctions();
  // Get number of Cholesky vectors
  unsigned int M = cdStorageController->getNVectors();
  // number of currently loaded vectors
  unsigned int batchSize = 0;

  for (unsigned int J = 0; J < M; J++) {
#pragma omp parallel
    {
      // Determine ThreadId
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif

      // load the current batch of Cholesky vectors
#pragma omp master
      {
        batchSize = cdStorageController->loadBatch(J);
        if (batchSize <= 0) {
          throw SerenityError("Cholesky Vector batch size is less than 1.");
        }
      }
#pragma omp barrier

      Eigen::setNbThreads(1);
#pragma omp for schedule(dynamic)
      // loop over all Cholesky vectors
      for (unsigned int K = J; K < J + batchSize; K++) {
        // load the current vector
        std::shared_ptr<Eigen::RowVectorXd> cholVector = cdStorageController->loadVector(K);
        // Transform Two Index Vector into its Matrix representation
        Eigen::Map<Eigen::MatrixXd> cholVecMat(cholVector->data(), n, n);

        // Calculate Coulomb Contribution
        double densTimesChol = densityMatrix.total().cwiseProduct(cholVecMat).sum();
        Eigen::MatrixXd tmp = densTimesChol * cholVecMat;
        fc[threadId]->alpha += tmp;
        fc[threadId]->beta += tmp;

        // Calculate Exchange contribution
        if (_xRatio > 0.0) {
          Eigen::MatrixXd half = cholVecMat * occCoeff.alpha;
          fx[threadId]->alpha -= _xRatio * half * esSigns.alpha.asDiagonal() * half.transpose();
          half = cholVecMat * occCoeff.beta;
          fx[threadId]->beta -= _xRatio * half * esSigns.beta.asDiagonal() * half.transpose();
        }
      }
      Eigen::setNbThreads(0);

#pragma omp master
      {
        cdStorageController->freeLastBatch();
        J += batchSize - 1;
      }
#pragma omp barrier
    }
  }

#ifdef _OPENMP
  for (unsigned int i = 1; i < nThreads; ++i) {
    fc[0]->alpha += fc[i]->alpha;
    fc[0]->beta += fc[i]->beta;
    fx[0]->alpha += fx[i]->alpha;
    fx[0]->beta += fx[i]->beta;
    delete fx[i];
    delete fc[i];
  }
  _fullXpotential->alpha = fx[0]->alpha;
  _fullXpotential->beta = fx[0]->beta;
#endif
  F.alpha = fc[0]->alpha + fx[0]->alpha;
  F.beta = fc[0]->beta + fx[0]->beta;
}

template class CDHFPotential<Options::SCF_MODES::RESTRICTED>;
template class CDHFPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
