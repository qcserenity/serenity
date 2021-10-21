/**
 * @file   CDExchangePotential.cpp
 *
 * @date   Apr 04, 2020
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
#include "potentials/CDExchangePotential.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/CDIntegralController.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "integrals/decomposer/SimpleCholeskyDecomposer.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "misc/SerenityError.h"
#include "misc/Timing.h"
#include "system/SystemController.h"

/* Include Std and External Headers */
#include <algorithm>

namespace Serenity {
template<Options::SCF_MODES SCFMode>
CDExchangePotential<SCFMode>::CDExchangePotential(std::shared_ptr<SystemController> systemController,
                                                  std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                                                  const double exchangeRatio, const double prescreeningThreshold,
                                                  LIBINT_OPERATOR op, const double mu)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _riints(nullptr),
    _systemController(systemController),
    _prescreeningThreshold(prescreeningThreshold),
    _exc(exchangeRatio),
    _dMatController(dMat),
    _excPotential(nullptr),
    _outOfDate(true),
    _op(op),
    _mu(mu) {
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& CDExchangePotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -   Exchange Pot.");
  if (_outOfDate) {
    DensityMatrix<SCFMode> densityMatrix = _dMatController->getDensityMatrix();
    _excPotential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot2 = *_excPotential;
    for_spin(pot2) {
      pot2_spin.setZero();
    };
    this->addToMatrix(*_excPotential, densityMatrix);
    _outOfDate = false;
  }
  Timings::timeTaken("Active System -   Exchange Pot.");
  return *_excPotential;
}

template<Options::SCF_MODES SCFMode>
double CDExchangePotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (_outOfDate)
    this->getMatrix();
  Timings::takeTime("Active System -   Exchange Pot.");
  auto& pot = *_excPotential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += 0.5 * pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("Active System -   Exchange Pot.");

  return energy;
};

template<>
void CDExchangePotential<Options::SCF_MODES::RESTRICTED>::addToMatrix(FockMatrix<Options::SCF_MODES::RESTRICTED>& F,
                                                                      const DensityMatrix<Options::SCF_MODES::RESTRICTED>& P) {
  if (_exc == 0.0)
    return;
  // Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fx;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  // Let the first thread use the Fock matrix directly
  //  fx.push_back(&F);
  for (unsigned int i = 0; i < nThreads; ++i) {
    fx.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
    fx[i]->setZero();
  }

#else
  fx.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
  fx[0]->setZero();
#endif

  takeTime("prepare");

  auto systemController = _systemController.lock();

  Eigen::VectorXd esSigns;
  unsigned int nOcc;
  Eigen::MatrixXd occCoeff;

  if (systemController->hasElectronicStructure<RESTRICTED>()) {
    //  if(false){
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
    // If there are no orbitals available generate pseudo coefficients from the density matrix
    Eigen::MatrixXd p = P;
    auto pseudoCoeff = systemController->getCDIntegralController()->generatePseudoCoefficients(p);
    esSigns = pseudoCoeff.first;
    occCoeff = pseudoCoeff.second;
    nOcc = occCoeff.cols();
  }

  auto& libint = Libint::getInstance();

  auto basisController = systemController->getBasisController();
  auto basis = basisController->getBasis();
  auto n = basisController->getNBasisFunctions();
  auto nShells = basis.size();

  std::shared_ptr<BasisController> auxBasisController;
  if (systemController->getSettings().basis.densityFitting == Options::DENS_FITS::ACD) {
    if (_op == LIBINT_OPERATOR::erf_coulomb) {
      auxBasisController = systemController->getBasisController(Options::BASIS_PURPOSES::ERF_ATOMIC_CHOLESKY);
    }
    else {
      auxBasisController = systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_CHOLESKY);
    }
  }
  else if (systemController->getSettings().basis.densityFitting == Options::DENS_FITS::ACCD) {
    if (_op == LIBINT_OPERATOR::erf_coulomb) {
      auxBasisController = systemController->getBasisController(Options::BASIS_PURPOSES::ERF_ATOMIC_COMPACT_CHOLESKY);
    }
    else {
      auxBasisController = systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY);
    }
  }
  else {
    auxBasisController = systemController->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  }

  auto auxBasis = auxBasisController->getBasis();
  auto m = auxBasisController->getNBasisFunctions();

  timeTaken(3, "prepare");

  libint.initialize(_op, 0, 3, std::vector<std::shared_ptr<Atom>>(0), _mu);

  takeTime("riints");
  // precalculate  two-center integrals and shellpairfactors
  if (!_riints)
    _riints = std::make_shared<RI_J_IntegralController>(basisController, auxBasisController, nullptr, _op, _mu);
  _riints->getInverseM();
  getShellPairFactor(0, 0);
  Eigen::MatrixXd minv = _riints->getInverseMSqrt();

  auxBasisController->getRIPrescreeningFactors();

  timeTaken(3, "riints");
  takeTime("memory");
  // Check available memory and set block size
  auto memManager = MemoryManager::getInstance();
  double availableMem = memManager->getAvailableSystemMemory();
  double additionalMemBuffer = 8 * n * m * omp_get_max_threads() * 28 * 20;
  additionalMemBuffer += 8 * (20 * 28) * (20 * 28) * (20 * 28);
  if (availableMem < 2e+9 + additionalMemBuffer) {
    _riints->clearCache();
    availableMem = memManager->getAvailableSystemMemory();
    if (availableMem < 2e+9 + additionalMemBuffer)
      throw SerenityError("Low Memory available");
  }
  availableMem = memManager->getAvailableSystemMemory();
  unsigned int maxBlockSize = std::floor((availableMem - 2e+9 - additionalMemBuffer) / (8 * nOcc * m));
  if (maxBlockSize > n)
    maxBlockSize = n;

  // Setup vector of matrices (nu,J) for each mu in the current block.
  std::vector<Eigen::MatrixXd> muMatrices;

  // while/for loop combination to loop over all mu-blocks
  unsigned int muShellCounter = 0;
  unsigned int muCounter = 0;
  timeTaken(3, "memory");
  // outer while loop over blocks (the number of iterations in this loop is the number of blocks)
  takeTime("loop over nShells");
  while (muShellCounter < nShells) {
    // get number of shells and basis functions for current Block
    unsigned int blockSize = 0;
    unsigned int nBlockShells = 0;
    for (unsigned int i = muShellCounter; i < nShells; i++) {
      auto nI = basis[i]->getNContracted();
      if (blockSize + nI > maxBlockSize)
        break;
      nBlockShells++;
      blockSize += nI;
    }
    for (unsigned int i = 0; i < blockSize; i++) {
      Eigen::MatrixXd tmp(m, nOcc);
      tmp.setZero();
      muMatrices.push_back(tmp);
    }

    takeTime("ints CDX");
    // loop over mu shells in the block
#pragma omp parallel for schedule(dynamic)
    for (unsigned int muShell = muShellCounter; muShell < muShellCounter + nBlockShells; muShell++) {
      const auto& basMu = *basis[muShell];
      const unsigned int nMu = basis[muShell]->getNContracted();
      const unsigned int firstMu = basisController->extendedIndex(muShell);

      // set up temporary matrix to store all integrals containing the current mu shell.
      Eigen::MatrixXd tmpIntegrals(n * m, nMu);
      tmpIntegrals.setZero();

      for (unsigned int nuShell = 0; nuShell < nShells; nuShell++) {
        Eigen::MatrixXd ints;

        // Very rudimentary prescreening
        double munuFactor = getShellPairFactor(muShell, nuShell);
        if (munuFactor < 1E-20)
          continue;

        const auto& basNu = *basis[nuShell];
        const unsigned int nNu = basis[nuShell]->getNContracted();
        const unsigned int firstNu = basisController->extendedIndex(nuShell);

        const auto& riPrescreeningFactors = auxBasisController->getRIPrescreeningFactors();

        // loop over shells from the auxiliary basis
        for (unsigned int jCount = 0; jCount < riPrescreeningFactors->size(); jCount++) {
          auto& jShell = (*riPrescreeningFactors)[jCount];
          auto& j = jShell.bf1;

          if (munuFactor * jShell.factor < _prescreeningThreshold)
            continue;

          const auto& basJ = *auxBasis[j];
          const unsigned int nJ = auxBasis[j]->getNContracted();
          const unsigned int firstJ = auxBasisController->extendedIndex(j);

          if (libint.compute(_op, 0, basJ, basMu, basNu, ints)) {
            // sort integrals into tmpIntegrals
            for (unsigned int mu = 0; mu < nMu; mu++) {
              for (unsigned int nu = 0; nu < nNu; nu++) {
                for (unsigned int J = 0; J < nJ; J++) {
                  tmpIntegrals((nu + firstNu) * m + (J + firstJ), mu) = ints.col(0)[J * nMu * nNu + mu * nNu + nu];
                }
              }
            }
          }

        } // end loop over jCount
      }   // end loop over nu
      // sort integrals into tmpIntegrals in n x m matrices for each mu and push to back of muMatrices
      for (unsigned int mu = 0; mu < nMu; mu++) {
        Eigen::MatrixXd tmp = Eigen::Map<Eigen::MatrixXd>(tmpIntegrals.col(mu).data(), m, n); // m x n
        tmp = tmp * occCoeff;                                                                 // m x i
        tmp = tmp.transpose() * minv;                                                         // m x i
        muMatrices[firstMu - muCounter + mu] = tmp;
      }
    } // end loop over mu in block
    timeTaken(3, "ints CDX");

    takeTime("eval CDX");
    // calculate Fock matrix contributions from muMatrices
    // (This means only contributions that can be calculated for mu and nu in the current muBlock)
#pragma omp parallel
    {
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      for (unsigned int mu = muCounter; mu < muCounter + blockSize; mu++) {
        auto localMu = esSigns.asDiagonal() * muMatrices[mu - muCounter];
#pragma omp for
        for (unsigned int nu = muCounter; nu <= mu; nu++) {
          (*fx[threadId])(mu, nu) = -0.5 * _exc * (localMu.cwiseProduct(muMatrices[nu - muCounter]).sum());
        }
      }
    }
    timeTaken(3, "eval CDX");

    /*********************************************
     *  Calculate integrals for all remaining mu
     *********************************************/
    //(If it is not possible to cache all integrals the remaining integrals will be calculated after one another
    // and will be evaluted with the integrals stored for the muBlock. When that is finished the next block
    // of muMatrices is calculated and the following step is repeated until all integral combinations are
    // accounted for. Warning for a number of blocks larger than 1 this results in many integral recalculations.)

    // loop over mu shells in the block
    for (unsigned int mu1Shell = muShellCounter + nBlockShells; mu1Shell < nShells; mu1Shell++) {
      Eigen::MatrixXd ints;

      const auto& basMu = *basis[mu1Shell];
      const unsigned int nMu = basis[mu1Shell]->getNContracted();
      const unsigned int firstMu = basisController->extendedIndex(mu1Shell);

      std::vector<Eigen::MatrixXd> mu1Matrices;

      // set up temporary matrix to store all integrals containing the current mu shell.
      Eigen::MatrixXd tmpIntegrals(n * m, nMu);
      tmpIntegrals.setZero();

      for (unsigned int nuShell = 0; nuShell < nShells; nuShell++) {
        double munuFactor = getShellPairFactor(mu1Shell, nuShell);
        if (munuFactor < 1E-14)
          continue;

        const auto& basNu = *basis[nuShell];
        const unsigned int nNu = basis[nuShell]->getNContracted();
        const unsigned int firstNu = basisController->extendedIndex(nuShell);

        const auto& riPrescreeningFactors = auxBasisController->getRIPrescreeningFactors();

        for (unsigned int jCount = 0; jCount < riPrescreeningFactors->size(); jCount++) {
          auto& jShell = (*riPrescreeningFactors)[jCount];
          auto& j = jShell.bf1;

          if (munuFactor * jShell.factor < _prescreeningThreshold)
            continue;

          const auto& basJ = *auxBasis[j];
          const unsigned int nJ = auxBasis[j]->getNContracted();
          const unsigned int firstJ = auxBasisController->extendedIndex(j);

          auto significant = libint.compute(_op, 0, basJ, basMu, basNu, ints);
          if (!significant)
            continue;

          // sort integrals into tmpIntegrals
          for (unsigned int mu = 0; mu < nMu; mu++) {
            for (unsigned int nu = 0; nu < nNu; nu++) {
              for (unsigned int J = 0; J < nJ; J++) {
                tmpIntegrals((nu + firstNu) * m + (J + firstJ), mu) = ints.col(0)[J * nMu * nNu + mu * nNu + nu];
              }
            }
          }
        }
      } // end of loop over nuShell

      // sort integrals int tmpIntegrals in n x m matrices for each mu and push to back of muMatrices
      for (unsigned int mu = 0; mu < nMu; mu++) {
        Eigen::MatrixXd tmp = Eigen::Map<Eigen::MatrixXd>(tmpIntegrals.col(mu).data(), m, n);
        mu1Matrices.push_back(tmp);
      }

      // half transform integrals in muMatrices and include two-center integrals
      for (unsigned int mu = 0; mu < nMu; mu++) {
        mu1Matrices[mu] = mu1Matrices[mu] * occCoeff;
        mu1Matrices[mu] = mu1Matrices[mu].transpose() * minv;
      }

      // calculate all mu nu contributions to the Fock matrix
      for (unsigned int mu = muCounter; mu < muCounter + blockSize; mu++) {
        for (unsigned int nu = firstMu; nu < firstMu + nMu; nu++) {
          (*fx[0])(mu, nu) =
              -0.5 * _exc *
              ((esSigns.asDiagonal() * muMatrices[mu - muCounter].cwiseProduct(mu1Matrices[nu - firstMu])).sum());
        }
      }
    } // end of loop over remaining mu

    // update general mu counter and adjust blocksize for the last block
    muShellCounter += nBlockShells;
    muCounter += blockSize;
    if (muCounter + maxBlockSize > n) {
      maxBlockSize = n - muCounter;
    }
    // reset old muMatrices
    muMatrices.clear();
  }
  timeTaken(3, "loop over nShells");
  takeTime("mirror CDX");
  // only one diagonal matrix was calculated, thus it has to be mirrored to obtain the full matrix
#pragma omp parallel
  {
#ifdef _OPENMP
    const unsigned int threadId = omp_get_thread_num();
#else
    const unsigned int threadId = 0;
#endif
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = i + 1; j < n; j++) {
        *(fx[threadId]->data() + (i + j * n)) = *(fx[threadId]->data() + (j + i * n));
      }
    }
  }
  timeTaken(3, "mirror CDX");

  libint.finalize(_op, 0, 3);

  takeTime("add fx CDX");

#ifdef _OPENMP
  for (unsigned int i = 1; i < nThreads; ++i) {
    *fx[0] += *fx[i];
    delete fx[i];
  }
#endif
  F = *fx[0];
  timeTaken(3, "add fx CDX");
}

template<>
void CDExchangePotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F,
                                                                        const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& P) {
  if (_exc == 0.0)
    return;
  // Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>*> fx;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  // Let the first thread use the Fock matrix directly
  for (unsigned int i = 0; i < nThreads; ++i) {
    fx.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
    fx[i]->alpha.setZero();
    fx[i]->beta.setZero();
  }

#else
  // Simply use the fock matrix directly
  fx.push_back(&F);
  fx[0]->alpha.setZero();
  fx[0]->beta.setZero();
#endif

  auto systemController = _systemController.lock();

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
    // If there are no orbitals available generate pseudo coefficients from the density matrix
    for_spin(occCoeff, nOcc, P, esSigns) {
      Eigen::MatrixXd p = P_spin;
      auto pseudoCoeff = systemController->getCDIntegralController()->generatePseudoCoefficients(p);
      esSigns_spin = pseudoCoeff.first;
      occCoeff_spin = pseudoCoeff.second;
      nOcc_spin = occCoeff_spin.cols();
    };
  }

  // Initialize controller and matrices
  //  const auto orbController = systemController->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>();
  //  auto coefficientMatrix = orbController->getCoefficients();
  //  // Get number of occupied orbitals
  //  auto nOcc = systemController->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  // Coefficient matrix of occupied orbitals
  Eigen::MatrixXd occCoeff_alpha = occCoeff.alpha;
  Eigen::MatrixXd occCoeff_beta = occCoeff.beta;

  auto& libint = Libint::getInstance();

  auto basisController = systemController->getBasisController();
  auto basis = basisController->getBasis();
  auto n = basisController->getNBasisFunctions();
  auto nShells = basis.size();

  std::shared_ptr<BasisController> auxBasisController;
  if (systemController->getSettings().basis.densityFitting == Options::DENS_FITS::ACD) {
    auxBasisController = systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_CHOLESKY);
  }
  else if (systemController->getSettings().basis.densityFitting == Options::DENS_FITS::ACCD) {
    auxBasisController = systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY);
  }
  else {
    auxBasisController = systemController->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  }

  auto auxBasis = auxBasisController->getBasis();
  auto m = auxBasisController->getNBasisFunctions();

  libint.initialize(_op, 0, 3, std::vector<std::shared_ptr<Atom>>(0), _mu);

  // precalculate  two-center integrals and shellpairfactors
  if (!_riints)
    _riints = std::make_shared<RI_J_IntegralController>(basisController, auxBasisController, nullptr, _op, _mu);
  _riints->getInverseM();
  getShellPairFactor(0, 0);
  Eigen::MatrixXd minv = _riints->getInverseMSqrt();
  auxBasisController->getRIPrescreeningFactors();

  // Check available memory and set block size
  auto memManager = MemoryManager::getInstance();
  auto availableMem = memManager->getAvailableSystemMemory();
  unsigned long long additionalMemBuffer = 8 * n * m * omp_get_max_threads() * 28 * 20;
  additionalMemBuffer += 8 * (20 * 28) * (20 * 28) * (20 * 28);
  if (availableMem < 2e+9 + additionalMemBuffer) {
    std::cout << "\ntry to free memory from cached 3C Integrals" << std::endl;
    _riints->clearCache();
    availableMem = memManager->getAvailableSystemMemory();
    if (availableMem < 2e+9 + additionalMemBuffer)
      throw SerenityError("Low Memory available");
  }

  unsigned int maxBlockSize = std::floor((availableMem - 2e+9 - additionalMemBuffer) / (8 * (nOcc.alpha + nOcc.beta) * m));
  if (maxBlockSize > n)
    maxBlockSize = n;

  // Setup vector of matrices (nu,J) for each mu in the current block.
  std::vector<Eigen::MatrixXd> muMatrices_alpha;
  std::vector<Eigen::MatrixXd> muMatrices_beta;

  // while/for loop combination to loop over all mu-blocks
  unsigned int muShellCounter = 0;
  unsigned int muCounter = 0;
  // outer while loop over blocks (the number of iterations in this loop is the number of blocks)
  while (muShellCounter < nShells) {
    // get number of shells and basis functions for current Block
    unsigned int blockSize = 0;
    unsigned int nBlockShells = 0;
    for (unsigned int i = muShellCounter; i < nShells; i++) {
      auto nI = basis[i]->getNContracted();
      if (blockSize + nI > maxBlockSize)
        break;
      nBlockShells++;
      blockSize += nI;
    }
    for (unsigned int i = 0; i < blockSize; i++) {
      Eigen::MatrixXd tmp(m, nOcc.alpha);
      tmp.setZero();
      muMatrices_alpha.push_back(tmp);
      Eigen::MatrixXd tmp1(m, nOcc.beta);
      tmp1.setZero();
      muMatrices_beta.push_back(tmp1);
    }

    // loop over mu shells in the block
#pragma omp parallel for schedule(dynamic)
    for (unsigned int muShell = muShellCounter; muShell < muShellCounter + nBlockShells; muShell++) {
      const auto& basMu = *basis[muShell];
      const unsigned int nMu = basis[muShell]->getNContracted();
      const unsigned int firstMu = basisController->extendedIndex(muShell);

      // set up temporary matrix to store all integrals containing the current mu shell.
      Eigen::MatrixXd tmpIntegrals(n * m, nMu);
      tmpIntegrals.setZero();

      for (unsigned int nuShell = 0; nuShell < nShells; nuShell++) {
        Eigen::MatrixXd ints;

        // Very rudimentary prescreening
        double munuFactor = getShellPairFactor(muShell, nuShell);
        if (munuFactor < 1E-14)
          continue;

        const auto& basNu = *basis[nuShell];
        const unsigned int nNu = basis[nuShell]->getNContracted();
        const unsigned int firstNu = basisController->extendedIndex(nuShell);

        const auto& riPrescreeningFactors = auxBasisController->getRIPrescreeningFactors();

        // loop over shells from the auxiliary basis
        for (unsigned int jCount = 0; jCount < riPrescreeningFactors->size(); jCount++) {
          auto& jShell = (*riPrescreeningFactors)[jCount];
          auto& j = jShell.bf1;

          if (munuFactor * jShell.factor < _prescreeningThreshold)
            continue;

          const auto& basJ = *auxBasis[j];
          const unsigned int nJ = auxBasis[j]->getNContracted();
          const unsigned int firstJ = auxBasisController->extendedIndex(j);

          if (libint.compute(_op, 0, basJ, basMu, basNu, ints)) {
            // sort integrals into tmpIntegrals
            for (unsigned int mu = 0; mu < nMu; mu++) {
              for (unsigned int nu = 0; nu < nNu; nu++) {
                for (unsigned int J = 0; J < nJ; J++) {
                  tmpIntegrals((nu + firstNu) * m + (J + firstJ), mu) = ints.col(0)[J * nMu * nNu + mu * nNu + nu];
                }
              }
            }
          }

        } // end loop over jCount
      }   // end loop over nu

      // sort integrals into tmpIntegrals in n x m matrices for each mu and push to back of muMatrices
      for (unsigned int mu = 0; mu < nMu; mu++) {
        Eigen::MatrixXd tmp = Eigen::Map<Eigen::MatrixXd>(tmpIntegrals.col(mu).data(), m, n);
        Eigen::MatrixXd tmp_alpha = tmp * occCoeff_alpha;
        tmp_alpha = tmp_alpha.transpose() * minv;
        muMatrices_alpha[firstMu - muCounter + mu] = tmp_alpha;
        Eigen::MatrixXd tmp_beta = tmp * occCoeff_beta;
        tmp_beta = tmp_beta.transpose() * minv;
        muMatrices_beta[firstMu - muCounter + mu] = tmp_beta;
      }
    } // end loop over mu in block

    // calculate Fock matrix contributions from muMatrices
    // (This means only contributions that can be calculated for mu and nu in the current muBlock)
#pragma omp parallel
    {
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      for (unsigned int mu = muCounter; mu < muCounter + blockSize; mu++) {
#pragma omp for
        for (unsigned int nu = muCounter; nu <= mu; nu++) {
          fx[threadId]->alpha(mu, nu) =
              -1 * _exc *
              ((esSigns.alpha.asDiagonal() * muMatrices_alpha[mu - muCounter].cwiseProduct(muMatrices_alpha[nu - muCounter]))
                   .sum());
          fx[threadId]->beta(mu, nu) =
              -1 * _exc *
              ((esSigns.beta.asDiagonal() * muMatrices_beta[mu - muCounter].cwiseProduct(muMatrices_beta[nu - muCounter]))
                   .sum());
        }
      }
    }

    /*********************************************
     *  Calculate integrals for all remaining mu
     *********************************************/
    //(If it is not possible to cache all integrals the remaining integrals will be calculated after one another
    // and will be evaluted with the integrals stored for the muBlock. When that is finished the next block
    // of muMatrices is calculated and the following step is repeated until all integral combinations are
    // accounted for. Warning for a number of blocks larger than 1 this results in many integral recalculations.)

    // loop over mu shells in the block
    for (unsigned int mu1Shell = muShellCounter + nBlockShells; mu1Shell < nShells; mu1Shell++) {
      Eigen::MatrixXd ints;

      const auto& basMu = *basis[mu1Shell];
      const unsigned int nMu = basis[mu1Shell]->getNContracted();
      const unsigned int firstMu = basisController->extendedIndex(mu1Shell);

      std::vector<Eigen::MatrixXd> innerMuMatrices_alpha;
      std::vector<Eigen::MatrixXd> innerMuMatrices_beta;

      // set up temporary matrix to store all integrals containing the current mu shell.
      Eigen::MatrixXd tmpIntegrals(n * m, nMu);
      tmpIntegrals.setZero();

      for (unsigned int nuShell = 0; nuShell < nShells; nuShell++) {
        double munuFactor = getShellPairFactor(mu1Shell, nuShell);
        if (munuFactor < 1E-14)
          continue;

        const auto& basNu = *basis[nuShell];
        const unsigned int nNu = basis[nuShell]->getNContracted();
        const unsigned int firstNu = basisController->extendedIndex(nuShell);

        const auto& riPrescreeningFactors = auxBasisController->getRIPrescreeningFactors();

        for (unsigned int jCount = 0; jCount < riPrescreeningFactors->size(); jCount++) {
          auto& jShell = (*riPrescreeningFactors)[jCount];
          auto& j = jShell.bf1;

          if (munuFactor * jShell.factor < _prescreeningThreshold)
            continue;

          const auto& basJ = *auxBasis[j];
          const unsigned int nJ = auxBasis[j]->getNContracted();
          const unsigned int firstJ = auxBasisController->extendedIndex(j);

          auto significant = libint.compute(_op, 0, basJ, basMu, basNu, ints);
          if (!significant)
            continue;

          // sort integrals into tmpIntegrals
          for (unsigned int mu = 0; mu < nMu; mu++) {
            for (unsigned int nu = 0; nu < nNu; nu++) {
              for (unsigned int J = 0; J < nJ; J++) {
                tmpIntegrals((nu + firstNu) * m + (J + firstJ), mu) = ints.col(0)[J * nMu * nNu + mu * nNu + nu];
              }
            }
          }
        }
      } // end of loop over nuShell

      // sort integrals int tmpIntegrals in n x m matrices for each mu and push to back of muMatrices
      // and half transform integrals in muMatrices and include two-center integrals
      for (unsigned int mu = 0; mu < nMu; mu++) {
        Eigen::MatrixXd tmp = Eigen::Map<Eigen::MatrixXd>(tmpIntegrals.col(mu).data(), m, n);

        Eigen::MatrixXd tmp_alpha = tmp * occCoeff_alpha;
        tmp_alpha = tmp_alpha.transpose() * minv;
        innerMuMatrices_alpha.push_back(tmp_alpha);
        Eigen::MatrixXd tmp_beta = tmp * occCoeff_beta;
        tmp_beta = tmp_beta.transpose() * minv;
        innerMuMatrices_beta.push_back(tmp_beta);
      }

      // calculate all mu nu contributions to the Fock matrix
      for (unsigned int mu = muCounter; mu < muCounter + blockSize; mu++) {
        for (unsigned int nu = firstMu; nu < firstMu + nMu; nu++) {
          fx[0]->alpha(nu, mu) =
              -1 * _exc *
              ((esSigns.alpha.asDiagonal() * muMatrices_alpha[mu - muCounter].cwiseProduct(innerMuMatrices_alpha[nu - firstMu]))
                   .sum());
          fx[0]->beta(nu, mu) =
              -1 * _exc *
              ((esSigns.beta.asDiagonal() * muMatrices_beta[mu - muCounter].cwiseProduct(innerMuMatrices_beta[nu - firstMu]))
                   .sum());
        }
      }
    } // end of loop over remaining mu

    // update general mu counter and adjust blocksize for the last block
    muShellCounter += nBlockShells;
    muCounter += blockSize;
    if (muCounter + maxBlockSize > n) {
      maxBlockSize = n - muCounter;
    }
    // reset old muMatrices
    muMatrices_alpha.clear();
    muMatrices_beta.clear();
  }

  // only one diagonal matrix was calculated, thus it has to be mirrored to obtain the full matrix
#pragma omp parallel
  {
#ifdef _OPENMP
    const unsigned int threadId = omp_get_thread_num();
#else
    const unsigned int threadId = 0;
#endif
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = i + 1; j < n; j++) {
        *(fx[threadId]->alpha.data() + (i + j * n)) = *(fx[threadId]->alpha.data() + (j + i * n));
        *(fx[threadId]->beta.data() + (i + j * n)) = *(fx[threadId]->beta.data() + (j + i * n));
      }
    }
  }

  libint.finalize(_op, 0, 3);

#ifdef _OPENMP
  for (unsigned int i = 1; i < nThreads; ++i) {
    fx[0]->alpha += fx[i]->alpha;
    fx[0]->beta += fx[i]->beta;
    delete fx[i];
  }
#endif
  F.alpha = fx[0]->alpha;
  F.beta = fx[0]->beta;
}

template<Options::SCF_MODES SCFMode>
double CDExchangePotential<SCFMode>::getShellPairFactor(unsigned int i, unsigned int j) {
#pragma omp critical
  {
    if (_shellPairFactors.size() == 0) {
      auto systemController = this->_systemController.lock();
      auto basisController = systemController->getBasisController();

      auto& basis = basisController->getBasis();
      auto& libint = Libint::getInstance();
      libint.initialize(_op, 0, 4, std::vector<std::shared_ptr<Atom>>(0), _mu);
      // loops over shells
      Eigen::MatrixXd integrals;
      const unsigned int nShells = basis.size();
      for (unsigned int i = 0; i < nShells; ++i) {
        std::vector<double> tmp;
        _shellPairFactors.push_back(tmp);
        const auto& shellI = *(basis)[i];
        for (unsigned int j = 0; j <= i; ++j) {
          const auto& shellJ = *(basis)[j];
          // calculate integrals
          if (libint.compute(_op, 0, shellI, shellJ, shellI, shellJ, integrals)) {
            double factor = std::sqrt(std::abs(integrals.col(0).maxCoeff()));
            _shellPairFactors[i].push_back(factor);
          }
          else {
            _shellPairFactors[i].push_back(0.0);
          }
        } /* j/shellJ */
      }   /* i/shellI */
      // finalize libint
      libint.finalize(_op, 0, 4);
    }
  }
  // change to shell indices
  return (i >= j ? _shellPairFactors[i][j] : _shellPairFactors[j][i]);
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd CDExchangePotential<SCFMode>::getGeomGradients() {
  throw SerenityError("Gradients not implemented yet");
}

template class CDExchangePotential<Options::SCF_MODES::RESTRICTED>;
template class CDExchangePotential<Options::SCF_MODES::UNRESTRICTED>;

} // namespace Serenity
