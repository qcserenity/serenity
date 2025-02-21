/**
 * @file   RIExchangePotential.cpp
 *
 * @date   Feb 23, 2022
 * @author Lars Hellmann
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
#include "potentials/RIExchangePotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/OrbitalController.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/CDIntegralController.h"
#include "integrals/RI_J_IntegralControllerFactory.h" //Construct RIJ integral controller.
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/looper/TwoElecThreeCenterCalculator.h"
#include "io/FormattedOutputStream.h" //Filtered output streams.
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
#include "parameters/Constants.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
RIExchangePotential<SCFMode>::RIExchangePotential(std::shared_ptr<SystemController> systemController,
                                                  std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                                                  std::shared_ptr<RI_J_IntegralController> ri_j_IntController,
                                                  const double exchangeRatio, const double prescreeningThreshold,
                                                  LIBINT_OPERATOR op, const double mu)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _systemController(systemController),
    _exc(exchangeRatio),
    _dMatController(dMat),
    _ri_j_IntController(ri_j_IntController),
    _fullpotential(nullptr),
    _outOfDate(true),
    _op(op),
    _mu(mu) {
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  (void)prescreeningThreshold;

  _fullpotential = std::make_shared<FockMatrix<SCFMode>>(FockMatrix<SCFMode>(this->_basis));
  auto& temp = *_fullpotential;
  for_spin(temp) {
    temp_spin.setZero();
  };

  if (_op == LIBINT_OPERATOR::coulomb) {
    _timingsLabel = "Active System -   Exchange Pot.";
  }
  else {
    _timingsLabel = "Active System -LR Exchange Pot.";
  }
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& RIExchangePotential<SCFMode>::getMatrix() {
  Timings::takeTime(_timingsLabel);
  if (_outOfDate) {
    DensityMatrix<SCFMode> densityMatrix = _dMatController->getDensityMatrix();
    _fullpotential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot2 = *_fullpotential;
    for_spin(pot2) {
      pot2_spin.setZero();
    };
    this->addToMatrix(*_fullpotential, densityMatrix);
    _outOfDate = false;
  }
  Timings::timeTaken(_timingsLabel);
  return *_fullpotential;
}

template<Options::SCF_MODES SCFMode>
double RIExchangePotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (_outOfDate)
    this->getMatrix();
  Timings::takeTime(_timingsLabel);
  auto& pot = *_fullpotential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += 0.5 * pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken(_timingsLabel);
  return energy;
};

template<>
void RIExchangePotential<Options::SCF_MODES::RESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::RESTRICTED>& F, const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densityMatrix) {
  OutputControl::vOut << " **** Performing RI-K ****" << std::endl;

  auto& libint = Libint::getInstance();
  libint.keepEngines(_op, 0, 3);
  auto systemController = _systemController.lock();

  const auto& basisController = _ri_j_IntController->getBasisController();
  auto auxBasisController = _ri_j_IntController->getAuxBasisController();

  // Call this once to initialize libint as needed
  {
    TwoElecThreeCenterCalculator intCalculator(_op, _mu, basisController, nullptr, auxBasisController,
                                               basisController->getPrescreeningThreshold());
  }

  // Obtain the number of basis functions in that system.
  const unsigned int nbfs = basisController->getNBasisFunctions();
  const unsigned int nAux = auxBasisController->getNBasisFunctions();

  Eigen::VectorXd esSigns;
  unsigned int nOcc;
  Eigen::MatrixXd occCoeff;

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
    // If there are no orbitals available generate pseudo coefficients from the density matrix
    Eigen::MatrixXd p = densityMatrix;
    auto pseudoCoeff = systemController->getCDIntegralController()->generatePseudoCoefficients(p);
    esSigns = pseudoCoeff.first;
    occCoeff = pseudoCoeff.second;
    nOcc = occCoeff.cols();
  }

  // Include factor of -1/2 and exchange ratio in esSigns
  //  esSigns *= -0.5 * _exc;

  // Initialization
  F.setZero();
  // AR: I guess this is necesssary so that InverseMSqrt is already calculated and can be accessed in a parallel region
  // which would otherwise cause problems
  _ri_j_IntController->getInverseMSqrt();
  unsigned nThreads = omp_get_max_threads();
  auto auxBas = auxBasisController->getBasis();

  // Get available system memory
  auto memManager = MemoryManager::getInstance();
  double availableMem = memManager->getAvailableSystemMemory();
  availableMem *= 0.9; // Include 10% Buffer

  // Initialize Blocking of occupied orbitals

  //  Calculate maximum number of auxiliary basis functions in a shell
  unsigned int maxShellBfs = 0;
  unsigned int maxAm = auxBasisController->getMaxAngularMomentum();
  if (auxBasisController->isPureSpherical())
    maxShellBfs = N_SHELL_SPH[maxAm];
  else
    maxShellBfs = N_SHELL_CART[maxAm];

  // Memory demand first step:
  double maxOcc1d = availableMem - sizeof(double) * nThreads * nbfs * nbfs * maxShellBfs;
  if (maxOcc1d < 0)
    throw SerenityError("Something went wrong in the memory allocation of the RIExchangePotential");
  maxOcc1d /= sizeof(double) * nbfs * nAux;
  unsigned int maxOcc1 = std::floor(maxOcc1d);

  // Memory demand second step:
  unsigned int maxOcc2 = std::floor(availableMem / (sizeof(double) * nbfs * nAux * 2));
  unsigned int occBlockSize = (maxOcc1 < maxOcc2) ? maxOcc1 : maxOcc2;
  if (occBlockSize > nOcc)
    occBlockSize = nOcc;

  // Details for verbose output
  OutputControl::vOut << " Exchange Ration:      " << _exc << std::endl;
  OutputControl::vOut << " Basis Functions:      " << nbfs << std::endl;
  OutputControl::vOut << " Aux. Basis Functions: " << nAux << std::endl;
  if (systemController->hasElectronicStructure<RESTRICTED>())
    OutputControl::vOut << " Occupied Orbitals:    " << nOcc << std::endl;
  else
    OutputControl::vOut << " Pseudo Orbitals:      " << nOcc << std::endl;
  OutputControl::vOut << " Blocking of occupied orbitals:" << std::endl;
  OutputControl::vOut << " \tAvailable Memory: " << 1e-9 * availableMem << " GB" << std::endl;
  OutputControl::vOut << "  \tBlocksize:        " << occBlockSize << std::endl;

  for (unsigned int iocc = 0; iocc < nOcc; iocc += occBlockSize) { // loop (in blocks) over all occupied orbitals
    if (iocc + occBlockSize >= nOcc)
      occBlockSize = nOcc - iocc; // Check and adjust block size
    /// Calculate integrals and perform half transformation
    Eigen::MatrixXd halfTransformedInts(Eigen::MatrixXd::Zero(nbfs * occBlockSize, nAux));
    {
      // Keep the calculator in this scope avoid unnecessary memory allocation
      TwoElecThreeCenterCalculator intCalculator(_op, _mu, basisController, nullptr, auxBasisController,
                                                 basisController->getPrescreeningThreshold());

#pragma omp parallel for schedule(dynamic)
      for (unsigned int pIndex = 0; pIndex < auxBas.size(); pIndex++) {
        auto integrals = intCalculator.calculateIntegrals(pIndex, omp_get_thread_num());
        auto firstP = auxBasisController->extendedIndex(pIndex);
        for (unsigned int p = 0; p < integrals.cols(); p++) {
          // Map column to matrix
          Eigen::Map<Eigen::MatrixXd> pIntegrals(integrals.col(p).data(), nbfs, nbfs);
          // Half transform
          Eigen::MatrixXd pHalfTransformed = pIntegrals * occCoeff.middleCols(iocc, occBlockSize);
          Eigen::Map<Eigen::VectorXd> pVector(pHalfTransformed.data(), nbfs * occBlockSize);
          halfTransformedInts.col(firstP + p) = pVector;
        }
      }
    }

    /// Normalization using the inverse sqrt of the 2-center matrix
    halfTransformedInts *= _ri_j_IntController->getInverseMSqrt();

    /// Evaluate contributions to the exchange matrix
    std::vector<MatrixInBasis<RESTRICTED>> fock_threads(nThreads, MatrixInBasis<RESTRICTED>(basisController));
#pragma omp parallel for schedule(dynamic)
    for (unsigned int p = 0; p < nAux; p++) {
      Eigen::Map<Eigen::MatrixXd> pIntegrals(halfTransformedInts.col(p).data(), nbfs, occBlockSize);
      fock_threads[omp_get_thread_num()] +=
          pIntegrals * esSigns.segment(iocc, occBlockSize).asDiagonal() * pIntegrals.transpose();
    }
    for (auto i : fock_threads) {
      F -= 0.5 * _exc * i;
    }
  }
  libint.freeEngines(_op, 0, 3);
}

template<>
void RIExchangePotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F, const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix) {
  OutputControl::vOut << " **** Performing RI-K ****" << std::endl;

  auto& libint = Libint::getInstance();
  libint.keepEngines(_op, 0, 3);
  auto systemController = _systemController.lock();

  const auto& basisController = _ri_j_IntController->getBasisController();
  auto auxBasisController = _ri_j_IntController->getAuxBasisController();

  // Call this once to initialize libint as needed
  {
    TwoElecThreeCenterCalculator intCalculator(_op, _mu, basisController, nullptr, auxBasisController,
                                               basisController->getPrescreeningThreshold());
  }

  // Obtain the number of basis functions in that system.
  const unsigned int nbfs = basisController->getNBasisFunctions();
  const unsigned int nAux = auxBasisController->getNBasisFunctions();

  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd> esSigns;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int> nOcc;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> occCoeff;

  if (systemController->hasElectronicStructure<UNRESTRICTED>()) {
    const auto orbController = systemController->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>();
    auto coefficientMatrix = orbController->getCoefficients();
    // Get number of occupied orbitals
    nOcc = systemController->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
    // Coefficient matrix of occupied orbitals
    for_spin(nOcc, occCoeff, esSigns, coefficientMatrix) {
      occCoeff_spin = coefficientMatrix_spin.leftCols(nOcc_spin);
      esSigns_spin = Eigen::VectorXd::Ones(nOcc_spin);
      // esSigns *= 2;
    };
  }
  else {
    // If there are no orbitals available generate pseudo coefficients from the density matrix
    DensityMatrix<Options::SCF_MODES::UNRESTRICTED> p = densityMatrix;
    for_spin(p, esSigns, occCoeff, nOcc) {
      auto pseudoCoeff = systemController->getCDIntegralController()->generatePseudoCoefficients(p_spin);
      esSigns_spin = pseudoCoeff.first;
      occCoeff_spin = pseudoCoeff.second;
      nOcc_spin = occCoeff_spin.cols();
    };
  }

  // Initialization
  for_spin(F) {
    F_spin.setZero();
  };
  _ri_j_IntController->getInverseMSqrt();
  unsigned nThreads = omp_get_max_threads();
  auto auxBas = auxBasisController->getBasis();

  // Get available system memory
  auto memManager = MemoryManager::getInstance();
  double availableMem = memManager->getAvailableSystemMemory();
  availableMem *= 0.9; // Include 10% Buffer

  // Initialize Blocking of occupied orbitals

  //  Calculate maximum number of auxiliary basis functions in a shell
  unsigned int maxShellBfs = 0;
  unsigned int maxAm = auxBasisController->getMaxAngularMomentum();
  if (auxBasisController->isPureSpherical())
    maxShellBfs = N_SHELL_SPH[maxAm];
  else
    maxShellBfs = N_SHELL_CART[maxAm];

  // Memory demand first step:
  double maxOcc1d = availableMem - sizeof(double) * nThreads * nbfs * nbfs * maxShellBfs;
  if (maxOcc1d < 0)
    throw SerenityError("Something went wrong in the memory allocation of the RIExchangePotential");
  maxOcc1d /= sizeof(double) * nbfs * nAux;
  unsigned int maxOcc1 = std::floor(maxOcc1d);

  // Memory demand second step:
  unsigned int maxOcc2 = std::floor(availableMem / (sizeof(double) * nbfs * nAux * 2));
  unsigned int occBlockSize = (maxOcc1 < maxOcc2) ? maxOcc1 : maxOcc2;
  if (occBlockSize > nOcc.alpha)
    occBlockSize = nOcc.alpha;

  // Details for verbose output
  OutputControl::vOut << " Exchange Ration:      " << _exc << std::endl;
  OutputControl::vOut << " Basis Functions:      " << nbfs << std::endl;
  OutputControl::vOut << " Aux. Basis Functions: " << nAux << std::endl;
  if (systemController->hasElectronicStructure<UNRESTRICTED>()) {
    OutputControl::vOut << " Occupied Orbitals (alpha):    " << nOcc.alpha << std::endl;
    OutputControl::vOut << " Occupied Orbitals (beta):     " << nOcc.beta << std::endl;
  }

  else {
    OutputControl::vOut << " Pseudo Orbitals (alpha):      " << nOcc.alpha << std::endl;
    OutputControl::vOut << " Pseudo Orbitals (beta):      " << nOcc.beta << std::endl;
  }

  OutputControl::vOut << " Blocking of occupied orbitals:" << std::endl;
  OutputControl::vOut << " \tAvailable Memory: " << 1e-9 * availableMem << " GB" << std::endl;
  OutputControl::vOut << "  \tBlocksize:        " << occBlockSize << std::endl;

  // Alpha

  for (unsigned int iocc = 0; iocc < nOcc.alpha; iocc += occBlockSize) { // loop (in blocks) over all occupied orbitals
    if (iocc + occBlockSize >= nOcc.alpha)
      occBlockSize = nOcc.alpha - iocc; // Check and adjust block size
    /// Calculate integrals and perform half transformation
    Eigen::MatrixXd halfTransformedInts(Eigen::MatrixXd::Zero(nbfs * occBlockSize, nAux));
    {
      // Keep the calculator in this scope avoid unnecessary memory allocation
      TwoElecThreeCenterCalculator intCalculator(_op, _mu, basisController, nullptr, auxBasisController,
                                                 basisController->getPrescreeningThreshold());

#pragma omp parallel for schedule(dynamic)
      for (unsigned int pIndex = 0; pIndex < auxBas.size(); pIndex++) {
        auto integrals = intCalculator.calculateIntegrals(pIndex, omp_get_thread_num());
        auto firstP = auxBasisController->extendedIndex(pIndex);
        for (unsigned int p = 0; p < integrals.cols(); p++) {
          // Map to vector
          Eigen::Map<Eigen::MatrixXd> pIntegrals(integrals.col(p).data(), nbfs, nbfs);
          // Half transform
          Eigen::MatrixXd pHalfTransformed = pIntegrals * occCoeff.alpha.middleCols(iocc, occBlockSize);
          Eigen::Map<Eigen::VectorXd> pVector(pHalfTransformed.data(), nbfs * occBlockSize);
          halfTransformedInts.col(firstP + p) = pVector;
        }
      }
    }

    /// Normalization using the inverse sqrt of the 2-center matrix
    halfTransformedInts *= _ri_j_IntController->getInverseMSqrt();

    /// Evaluate contributions to the exchange matrix
    std::vector<MatrixInBasis<RESTRICTED>> fock_threads(nThreads, MatrixInBasis<RESTRICTED>(basisController));
#pragma omp parallel for schedule(dynamic)
    for (unsigned int p = 0; p < nAux; p++) {
      Eigen::Map<Eigen::MatrixXd> pIntegrals(halfTransformedInts.col(p).data(), nbfs, occBlockSize);
      fock_threads[omp_get_thread_num()] +=
          pIntegrals * esSigns.alpha.segment(iocc, occBlockSize).asDiagonal() * pIntegrals.transpose();
    }
    for (auto i : fock_threads) {
      F.alpha -= 1.0 * _exc * i;
    }
  }

  // Beta

  occBlockSize = (maxOcc1 < maxOcc2) ? maxOcc1 : maxOcc2;
  if (occBlockSize > nOcc.beta)
    occBlockSize = nOcc.beta;

  for (unsigned int iocc = 0; iocc < nOcc.beta; iocc += occBlockSize) { // loop (in blocks) over all occupied orbitals
    if (iocc + occBlockSize >= nOcc.beta)
      occBlockSize = nOcc.beta - iocc; // Check and adjust block size
    /// Calculate integrals and perform half transformation
    Eigen::MatrixXd halfTransformedInts(Eigen::MatrixXd::Zero(nbfs * occBlockSize, nAux));
    {
      // Keep the calculator in this scope avoid unnecessary memory allocation
      TwoElecThreeCenterCalculator intCalculator(_op, _mu, basisController, nullptr, auxBasisController,
                                                 basisController->getPrescreeningThreshold());

#pragma omp parallel for schedule(dynamic)
      for (unsigned int pIndex = 0; pIndex < auxBas.size(); pIndex++) {
        auto integrals = intCalculator.calculateIntegrals(pIndex, omp_get_thread_num());
        auto firstP = auxBasisController->extendedIndex(pIndex);
        for (unsigned int p = 0; p < integrals.cols(); p++) {
          // Map to vector
          Eigen::Map<Eigen::MatrixXd> pIntegrals(integrals.col(p).data(), nbfs, nbfs);
          // Half transform
          Eigen::MatrixXd pHalfTransformed = pIntegrals * occCoeff.beta.middleCols(iocc, occBlockSize);
          Eigen::Map<Eigen::VectorXd> pVector(pHalfTransformed.data(), nbfs * occBlockSize);
          halfTransformedInts.col(firstP + p) = pVector;
        }
      }
    }

    /// Normalization using the inverse sqrt of the 2-center matrix
    halfTransformedInts *= _ri_j_IntController->getInverseMSqrt();

    /// Evaluate contributions to the exchange matrix
    std::vector<MatrixInBasis<RESTRICTED>> fock_threads(nThreads, MatrixInBasis<RESTRICTED>(basisController));
#pragma omp parallel for schedule(dynamic)
    for (unsigned int p = 0; p < nAux; p++) {
      Eigen::Map<Eigen::MatrixXd> pIntegrals(halfTransformedInts.col(p).data(), nbfs, occBlockSize);
      fock_threads[omp_get_thread_num()] +=
          pIntegrals * esSigns.beta.segment(iocc, occBlockSize).asDiagonal() * pIntegrals.transpose();
    }
    for (auto i : fock_threads) {
      F.beta -= 1.0 * _exc * i;
    }
  }

  libint.freeEngines(_op, 0, 3);
}

template<>
Eigen::MatrixXd RIExchangePotential<RESTRICTED>::getGeomGradients() {
  throw SerenityError("RIExchange Gradients not implemented yet");
}

template<>
Eigen::MatrixXd RIExchangePotential<UNRESTRICTED>::getGeomGradients() {
  throw SerenityError("RIExchange Gradients not implemented yet");
}

template class RIExchangePotential<Options::SCF_MODES::RESTRICTED>;
template class RIExchangePotential<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
