/**
 * @file ExchangePotential.cpp
 *
 * @author Moritz Bensberg
 * @date Sep 23, 2019
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
#include "potentials/ExchangePotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "io/FormattedOutputStream.h" //Filtered output streams.
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
#include "system/SystemController.h"

namespace Serenity {
template<Options::SCF_MODES SCFMode>
ExchangePotential<SCFMode>::ExchangePotential(std::shared_ptr<SystemController> systemController,
                                              std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                                              const double exchangeRatio, const double prescreeningThreshold,
                                              double prescreeningIncrementStart, double prescreeningIncrementEnd,
                                              unsigned int incrementSteps, bool clear4CenterCache, bool transitionDensity)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _systemController(systemController),
    _exc(exchangeRatio),
    _dMatController(dMat),
    _fullpotential(nullptr),
    _outOfDate(true),
    _incrementHelper(std::make_shared<IncrementalFockMatrix<SCFMode>>(
        dMat, (prescreeningThreshold != 0) ? prescreeningThreshold : this->_basis->getPrescreeningThreshold(),
        prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps, "Exact Exchange")),
    _clear4CenterCache(clear4CenterCache),
    _transitionDensity(transitionDensity) {
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);

  _fullpotential = std::make_shared<FockMatrix<SCFMode>>(FockMatrix<SCFMode>(this->_basis));
  auto& temp = *_fullpotential;
  for_spin(temp) {
    temp_spin.setZero();
  };
  _screening = prescreeningIncrementStart;
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& ExchangePotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -   Exchange Pot.");
  if (_outOfDate) {
    DensityMatrix<SCFMode> densityMatrix(this->_basis);
    std::vector<std::shared_ptr<FockMatrix<SCFMode>>> matrices = {_fullpotential};
    _incrementHelper->updateDensityAndThreshold(densityMatrix, _screening, matrices);
    if (_transitionDensity)
      this->addToMatrixForTransitionDensity(*_fullpotential, densityMatrix);
    else
      this->addToMatrix(*_fullpotential, densityMatrix);
    _outOfDate = false;
  }
  Timings::timeTaken("Active System -   Exchange Pot.");
  return *_fullpotential;
}

template<Options::SCF_MODES SCFMode>
double ExchangePotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (_outOfDate)
    this->getMatrix();
  Timings::takeTime("Active System -   Exchange Pot.");
  auto& pot = *_fullpotential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += 0.5 * pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("Active System -   Exchange Pot.");
  return energy;
};

template<>
void ExchangePotential<Options::SCF_MODES::RESTRICTED>::addToMatrix(FockMatrix<Options::SCF_MODES::RESTRICTED>& F,
                                                                    const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densityMatrix) {
  unsigned nb = _basis->getNBasisFunctions();
  unsigned ns = _basis->getReducedNBasisFunctions();
  unsigned nThreads = omp_get_max_threads();
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>> fx(nThreads, _basis);

  auto D = densityMatrix.data();

  auto distribute = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned iThread) {
    unsigned ik = i * nb + k;
    unsigned il = i * nb + l;
    unsigned jl = j * nb + l;
    unsigned jk = j * nb + k;

    auto Fx = fx[iThread].data();

    Fx[ik] -= D[jl] * integral;
    Fx[il] -= D[jk] * integral;
    Fx[jk] -= D[il] * integral;
    Fx[jl] -= D[ik] * integral;
  };

  auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  auto maxDensPtr = maxDensMat.data();
  auto maxDens = maxDensMat.maxCoeff();

  auto prescreen = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    double xschwarz = 0.5 * _exc * schwarz;
    if (maxDens * xschwarz < _screening) {
      return true;
    }
    double t1 = maxDensPtr[i * ns + k];
    double t2 = maxDensPtr[i * ns + l];
    double t3 = maxDensPtr[j * ns + k];
    double t4 = maxDensPtr[j * ns + l];
    double maxDBlock = std::max({t1, t2, t3, t4});
    return (maxDBlock * xschwarz < _screening);
  };

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, _incrementHelper->getPrescreeningThreshold());
  looper.loopNoDerivative(distribute, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);

  for (unsigned iThread = 1; iThread < nThreads; ++iThread) {
    fx[0] += fx[iThread];
  }

  // Symmetrize.
  fx[0] *= 0.5 * _exc;
  F += fx[0];
  F += fx[0].transpose();
}

template<>
void ExchangePotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F, const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix) {
  unsigned nb = _basis->getNBasisFunctions();
  unsigned ns = _basis->getReducedNBasisFunctions();
  unsigned nThreads = omp_get_max_threads();
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>> fx(nThreads, _basis);

  auto Da = densityMatrix.alpha.data();
  auto Db = densityMatrix.beta.data();

  auto distribute = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned ik = i * nb + k;
    unsigned il = i * nb + l;
    unsigned jl = j * nb + l;
    unsigned jk = j * nb + k;

    auto Fxa = fx[threadId].alpha.data();
    Fxa[ik] -= Da[jl] * integral;
    Fxa[il] -= Da[jk] * integral;
    Fxa[jk] -= Da[il] * integral;
    Fxa[jl] -= Da[ik] * integral;

    auto Fxb = fx[threadId].beta.data();
    Fxb[ik] -= Db[jl] * integral;
    Fxb[il] -= Db[jk] * integral;
    Fxb[jk] -= Db[il] * integral;
    Fxb[jl] -= Db[ik] * integral;
  };

  auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  auto maxDensPtr = maxDensMat.data();
  auto maxDens = maxDensMat.maxCoeff();

  auto prescreen = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    double xschwarz = _exc * schwarz;
    if (maxDens * xschwarz < _screening) {
      return true;
    }
    double t1 = maxDensPtr[i * ns + k];
    double t2 = maxDensPtr[i * ns + l];
    double t3 = maxDensPtr[j * ns + k];
    double t4 = maxDensPtr[j * ns + l];
    double maxDBlock = std::max({t1, t2, t3, t4});
    return (maxDBlock * xschwarz < _screening);
  };

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, _incrementHelper->getPrescreeningThreshold());
  looper.loopNoDerivative(distribute, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);

  for (unsigned iThread = 1; iThread < nThreads; ++iThread) {
    fx[0].alpha += fx[iThread].alpha;
    fx[0].beta += fx[iThread].beta;
  }

  // Symmetrize.
  fx[0].alpha *= _exc;
  F.alpha += fx[0].alpha;
  F.alpha += fx[0].alpha.transpose();
  fx[0].beta *= _exc;
  F.beta += fx[0].beta;
  F.beta += fx[0].beta.transpose();
}

template<>
void ExchangePotential<Options::SCF_MODES::RESTRICTED>::addToMatrixForTransitionDensity(
    FockMatrix<Options::SCF_MODES::RESTRICTED>& F, const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densityMatrix) {
  unsigned nb = _basis->getNBasisFunctions();
  unsigned ns = _basis->getReducedNBasisFunctions();
  unsigned nThreads = omp_get_max_threads();
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>> fx(nThreads, _basis);

  auto D = densityMatrix.data();

  auto distribute = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned iThread) {
    unsigned ik = i * nb + k;
    unsigned il = i * nb + l;
    unsigned jl = j * nb + l;
    unsigned jk = j * nb + k;
    unsigned ki = k * nb + i;
    unsigned li = l * nb + i;
    unsigned lj = l * nb + j;
    unsigned kj = k * nb + j;

    auto Fx = fx[iThread].data();

    Fx[ik] -= D[jl] * integral;
    Fx[il] -= D[jk] * integral;
    Fx[jk] -= D[il] * integral;
    Fx[jl] -= D[ik] * integral;
    Fx[ki] -= D[lj] * integral;
    Fx[li] -= D[kj] * integral;
    Fx[kj] -= D[li] * integral;
    Fx[lj] -= D[ki] * integral;
  };

  auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  auto maxDensPtr = maxDensMat.data();
  auto maxDens = maxDensMat.maxCoeff();

  auto prescreen = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    double xschwarz = 0.5 * _exc * schwarz;
    if (maxDens * xschwarz < _screening) {
      return true;
    }
    double t1 = maxDensPtr[i * ns + k];
    double t2 = maxDensPtr[i * ns + l];
    double t3 = maxDensPtr[j * ns + k];
    double t4 = maxDensPtr[j * ns + l];
    double t5 = maxDensPtr[k * ns + i];
    double t6 = maxDensPtr[l * ns + i];
    double t7 = maxDensPtr[k * ns + j];
    double t8 = maxDensPtr[l * ns + j];
    double maxDBlock = std::max({t1, t2, t3, t4, t5, t6, t7, t8});
    return (maxDBlock * xschwarz < _screening);
  };

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, _incrementHelper->getPrescreeningThreshold());
  looper.loopNoDerivative(distribute, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);

  for (unsigned iThread = 1; iThread < nThreads; ++iThread) {
    fx[0] += fx[iThread];
  }

  fx[0] *= 0.5 * _exc;
  F += fx[0];
}

template<>
void ExchangePotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrixForTransitionDensity(
    FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F, const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix) {
  unsigned nb = _basis->getNBasisFunctions();
  unsigned ns = _basis->getReducedNBasisFunctions();
  unsigned nThreads = omp_get_max_threads();
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>> fx(nThreads, _basis);

  auto Da = densityMatrix.alpha.data();
  auto Db = densityMatrix.beta.data();

  auto distribute = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned ik = i * nb + k;
    unsigned il = i * nb + l;
    unsigned jl = j * nb + l;
    unsigned jk = j * nb + k;
    unsigned ki = k * nb + i;
    unsigned li = l * nb + i;
    unsigned lj = l * nb + j;
    unsigned kj = k * nb + j;

    auto Fxa = fx[threadId].alpha.data();
    Fxa[ik] -= Da[jl] * integral;
    Fxa[il] -= Da[jk] * integral;
    Fxa[jk] -= Da[il] * integral;
    Fxa[jl] -= Da[ik] * integral;
    Fxa[ki] -= Da[lj] * integral;
    Fxa[li] -= Da[kj] * integral;
    Fxa[kj] -= Da[li] * integral;
    Fxa[lj] -= Da[ki] * integral;

    auto Fxb = fx[threadId].beta.data();
    Fxb[ik] -= Db[jl] * integral;
    Fxb[il] -= Db[jk] * integral;
    Fxb[jk] -= Db[il] * integral;
    Fxb[jl] -= Db[ik] * integral;
    Fxb[ki] -= Db[lj] * integral;
    Fxb[li] -= Db[kj] * integral;
    Fxb[kj] -= Db[li] * integral;
    Fxb[lj] -= Db[ki] * integral;
  };

  auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  auto maxDensPtr = maxDensMat.data();
  auto maxDens = maxDensMat.maxCoeff();

  auto prescreen = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    double xschwarz = _exc * schwarz;
    if (maxDens * xschwarz < _screening) {
      return true;
    }
    double t1 = maxDensPtr[i * ns + k];
    double t2 = maxDensPtr[i * ns + l];
    double t3 = maxDensPtr[j * ns + k];
    double t4 = maxDensPtr[j * ns + l];
    double t5 = maxDensPtr[k * ns + i];
    double t6 = maxDensPtr[l * ns + i];
    double t7 = maxDensPtr[k * ns + j];
    double t8 = maxDensPtr[l * ns + j];
    double maxDBlock = std::max({t1, t2, t3, t4, t5, t6, t7, t8});
    return (maxDBlock * xschwarz < _screening);
  };

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, _incrementHelper->getPrescreeningThreshold());
  looper.loopNoDerivative(distribute, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);

  for (unsigned iThread = 1; iThread < nThreads; ++iThread) {
    fx[0].alpha += fx[iThread].alpha;
    fx[0].beta += fx[iThread].beta;
  }

  fx[0].alpha *= _exc;
  F.alpha += fx[0].alpha;
  fx[0].beta *= _exc;
  F.beta += fx[0].beta;
}

template<>
Eigen::MatrixXd ExchangePotential<RESTRICTED>::getGeomGradients() {
  auto systemController = _systemController.lock();
  auto atoms = systemController->getAtoms();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _dMatController->getDensityMatrix();
  auto dMatPtr = densityMatrix.data();

  auto mapping = systemController->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

  const unsigned int nBFs = _basis->getNBasisFunctions();

#ifdef _OPENMP
  std::vector<Eigen::MatrixXd> priv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
#else
  std::vector<Eigen::MatrixXd> priv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
#endif

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _basis, _incrementHelper->getPrescreeningThreshold());
  const auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  const auto maxDens = maxDensMat.maxCoeff();
  auto const looperFunction = [&](const unsigned int& i, const unsigned int& j, const unsigned int& a,
                                  const unsigned int& b, const Eigen::VectorXd& intValues, const unsigned int& threadID) {
    double perm = (i == j) ? 1.0 : 2.0;
    perm *= (a == b) ? 1.0 : 2.0;
    perm *= (i == a) ? (j == b ? 1.0 : 2.0) : 2.0;
    perm *= 0.5;

    const unsigned int iAtom = mapping[i];
    const unsigned int jAtom = mapping[j];
    const unsigned int aAtom = mapping[a];
    const unsigned int bAtom = mapping[b];
    const double prefac = dMatPtr[i * nBFs + a] * dMatPtr[j * nBFs + b] + dMatPtr[i * nBFs + b] * dMatPtr[j * nBFs + a];
    priv[threadID].data()[0 * nAtoms + iAtom] += 0.25 * perm * prefac * intValues.data()[0];  // 0 * 3 + 0 = (atom, xyz)
    priv[threadID].data()[1 * nAtoms + iAtom] += 0.25 * perm * prefac * intValues.data()[1];  // 0 * 3 + 1 = (atom, xyz)
    priv[threadID].data()[2 * nAtoms + iAtom] += 0.25 * perm * prefac * intValues.data()[2];  // 0 * 3 + 2 = (atom, xyz)
    priv[threadID].data()[0 * nAtoms + jAtom] += 0.25 * perm * prefac * intValues.data()[3];  // 1 * 3 + 0 = (atom, xyz)
    priv[threadID].data()[1 * nAtoms + jAtom] += 0.25 * perm * prefac * intValues.data()[4];  // 1 * 3 + 1 = (atom, xyz)
    priv[threadID].data()[2 * nAtoms + jAtom] += 0.25 * perm * prefac * intValues.data()[5];  // 1 * 3 + 2 = (atom, xyz)
    priv[threadID].data()[0 * nAtoms + aAtom] += 0.25 * perm * prefac * intValues.data()[6];  // 2 * 3 + 0 = (atom, xyz)
    priv[threadID].data()[1 * nAtoms + aAtom] += 0.25 * perm * prefac * intValues.data()[7];  // 2 * 3 + 1 = (atom, xyz)
    priv[threadID].data()[2 * nAtoms + aAtom] += 0.25 * perm * prefac * intValues.data()[8];  // 2 * 3 + 2 = (atom, xyz)
    priv[threadID].data()[0 * nAtoms + bAtom] += 0.25 * perm * prefac * intValues.data()[9];  // 3 * 3 + 0 = (atom, xyz)
    priv[threadID].data()[1 * nAtoms + bAtom] += 0.25 * perm * prefac * intValues.data()[10]; // 3 * 3 + 1 = (atom, xyz)
    priv[threadID].data()[2 * nAtoms + bAtom] += 0.25 * perm * prefac * intValues.data()[11]; // 3 * 3 + 2 = (atom, xyz)
  };

  looper.loop(looperFunction, maxDens);

  Eigen::MatrixXd xGrad(nAtoms, 3);
  xGrad.setZero();
#ifdef _OPENMP
  // sum over all threads
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    xGrad -= priv[i] * _exc;
  }
#else
  xGrad -= priv[0] * _exc;
#endif
  return xGrad;
}

template<>
Eigen::MatrixXd ExchangePotential<UNRESTRICTED>::getGeomGradients() {
  auto systemController = _systemController.lock();
  auto atoms = systemController->getAtoms();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _dMatController->getDensityMatrix();
  auto dMatPtrA = densityMatrix.alpha.data();
  auto dMatPtrB = densityMatrix.beta.data();

  auto mapping = systemController->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

  const unsigned int nBFs = _basis->getNBasisFunctions();

#ifdef _OPENMP
  std::vector<Eigen::MatrixXd> priv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
#else
  std::vector<Eigen::MatrixXd> priv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
#endif

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _basis, _incrementHelper->getPrescreeningThreshold());
  const auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  const auto maxDens = maxDensMat.maxCoeff();
  auto const looperFunction = [&](const unsigned int& i, const unsigned int& j, const unsigned int& a,
                                  const unsigned int& b, const Eigen::VectorXd& intValues, const unsigned int& threadID) {
    double perm = (i == j) ? 1.0 : 2.0;
    perm *= (a == b) ? 1.0 : 2.0;
    perm *= (i == a) ? (j == b ? 1.0 : 2.0) : 2.0;
    perm *= 0.5;

    const unsigned int iAtom = mapping[i];
    const unsigned int jAtom = mapping[j];
    const unsigned int aAtom = mapping[a];
    const unsigned int bAtom = mapping[b];
    const double prefac =
        dMatPtrA[i * nBFs + a] * dMatPtrA[j * nBFs + b] + dMatPtrA[i * nBFs + b] * dMatPtrA[j * nBFs + a] +
        dMatPtrB[i * nBFs + a] * dMatPtrB[j * nBFs + b] + dMatPtrB[i * nBFs + b] * dMatPtrB[j * nBFs + a];
    priv[threadID].data()[0 * nAtoms + iAtom] += 0.5 * perm * prefac * intValues.data()[0];  // 0 * 3 + 0 = (atom, xyz)
    priv[threadID].data()[1 * nAtoms + iAtom] += 0.5 * perm * prefac * intValues.data()[1];  // 0 * 3 + 1 = (atom, xyz)
    priv[threadID].data()[2 * nAtoms + iAtom] += 0.5 * perm * prefac * intValues.data()[2];  // 0 * 3 + 2 = (atom, xyz)
    priv[threadID].data()[0 * nAtoms + jAtom] += 0.5 * perm * prefac * intValues.data()[3];  // 1 * 3 + 0 = (atom, xyz)
    priv[threadID].data()[1 * nAtoms + jAtom] += 0.5 * perm * prefac * intValues.data()[4];  // 1 * 3 + 1 = (atom, xyz)
    priv[threadID].data()[2 * nAtoms + jAtom] += 0.5 * perm * prefac * intValues.data()[5];  // 1 * 3 + 2 = (atom, xyz)
    priv[threadID].data()[0 * nAtoms + aAtom] += 0.5 * perm * prefac * intValues.data()[6];  // 2 * 3 + 0 = (atom, xyz)
    priv[threadID].data()[1 * nAtoms + aAtom] += 0.5 * perm * prefac * intValues.data()[7];  // 2 * 3 + 1 = (atom, xyz)
    priv[threadID].data()[2 * nAtoms + aAtom] += 0.5 * perm * prefac * intValues.data()[8];  // 2 * 3 + 2 = (atom, xyz)
    priv[threadID].data()[0 * nAtoms + bAtom] += 0.5 * perm * prefac * intValues.data()[9];  // 3 * 3 + 0 = (atom, xyz)
    priv[threadID].data()[1 * nAtoms + bAtom] += 0.5 * perm * prefac * intValues.data()[10]; // 3 * 3 + 1 = (atom, xyz)
    priv[threadID].data()[2 * nAtoms + bAtom] += 0.5 * perm * prefac * intValues.data()[11]; // 3 * 3 + 2 = (atom, xyz)
  };

  looper.loop(looperFunction, maxDens);

  Eigen::MatrixXd xGrad(nAtoms, 3);
  xGrad.setZero();
#ifdef _OPENMP
  // sum over all threads
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    xGrad -= priv[i] * _exc;
  }
#else
  xGrad -= priv[0] * _exc;
#endif
  return xGrad;
}

template<Options::SCF_MODES SCFMode>
ExchangePotential<SCFMode>::~ExchangePotential() {
  if (_systemController.lock() && _clear4CenterCache)
    _systemController.lock()->clear4CenterCache();
}

template class ExchangePotential<Options::SCF_MODES::RESTRICTED>;
template class ExchangePotential<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
