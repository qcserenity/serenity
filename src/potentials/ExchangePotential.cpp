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

namespace Serenity {
template<Options::SCF_MODES SCFMode>
ExchangePotential<SCFMode>::ExchangePotential(std::shared_ptr<SystemController> systemController,
                                              std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                                              const double exchangeRatio, const double prescreeningThreshold,
                                              double prescreeningIncrementStart, double prescreeningIncrementEnd,
                                              unsigned int incrementSteps, bool clear4CenterCache)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _systemController(systemController),
    _exc(exchangeRatio),
    _dMatController(dMat),
    _fullpotential(nullptr),
    _outOfDate(true),
    _incrementHelper(std::make_shared<IncrementalFockMatrix<SCFMode>>(
        dMat, (prescreeningThreshold != 0) ? prescreeningThreshold : this->_basis->getPrescreeningThreshold(),
        prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps, "Exact Exchange")),
    _clear4CenterCache(clear4CenterCache) {
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
  const unsigned int nBFs = _basis->getNBasisFunctions();

  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>> fx;
  const unsigned int nThreads = omp_get_max_threads();
  for (unsigned int i = 0; i < nThreads; ++i) {
    fx.emplace_back(_basis);
  }

  auto distribute = [&](unsigned i, unsigned j, unsigned k, const unsigned l, double integral, unsigned threadId) {
    unsigned long ik = i * nBFs + k;
    unsigned long il = i * nBFs + l;
    unsigned long jl = j * nBFs + l;
    unsigned long jk = j * nBFs + k;

    auto Fx = fx[threadId].data();
    auto D = densityMatrix.data();

    // Exchange.
    const double exc = 0.5 * _exc * integral;
    *(Fx + ik) -= *(D + jl) * exc;
    *(Fx + il) -= *(D + jk) * exc;
    *(Fx + jk) -= *(D + il) * exc;
    *(Fx + jl) -= *(D + ik) * exc;
  };

  const auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  const auto maxDens = maxDensMat.maxCoeff();
  auto prescreen = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwartz) {
    if (maxDens * schwartz < _screening) {
      return true;
    }
    double maxDBlock = 0.5 * maxDensMat(i, k);
    maxDBlock = std::max(maxDBlock, 0.5 * maxDensMat(i, l));
    maxDBlock = std::max(maxDBlock, 0.5 * maxDensMat(j, k));
    maxDBlock = std::max(maxDBlock, 0.5 * maxDensMat(j, l));
    if (maxDBlock * schwartz < _screening) {
      return true;
    }
    return false;
  };

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, _incrementHelper->getPrescreeningThreshold());
  looper.loopNoDerivative(distribute, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);

  for (unsigned int i = 1; i < nThreads; ++i) {
    fx[0] += fx[i];
  }

  // Symmetrize.
  Eigen::Ref<Eigen::MatrixXd> Fxa = fx[0];
  Eigen::MatrixXd tmp_x = Fxa + Fxa.transpose();
  F += tmp_x;
}

template<>
void ExchangePotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F, const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix) {
  const unsigned int nBFs = _basis->getNBasisFunctions();

  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>> fx;
  const unsigned int nThreads = omp_get_max_threads();
  for (unsigned int i = 0; i < nThreads; ++i) {
    fx.emplace_back(_basis);
  }

  auto distribute = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned long ik = i * nBFs + k;
    unsigned long il = i * nBFs + l;
    unsigned long jl = j * nBFs + l;
    unsigned long jk = j * nBFs + k;

    auto Fxa = fx[threadId].alpha.data();
    auto Fxb = fx[threadId].beta.data();

    auto Da = densityMatrix.alpha.data();
    auto Db = densityMatrix.beta.data();

    double exc = _exc * integral;
    *(Fxa + ik) -= *(Da + jl) * exc;
    *(Fxa + il) -= *(Da + jk) * exc;
    *(Fxa + jk) -= *(Da + il) * exc;
    *(Fxa + jl) -= *(Da + ik) * exc;

    *(Fxb + ik) -= *(Db + jl) * exc;
    *(Fxb + il) -= *(Db + jk) * exc;
    *(Fxb + jk) -= *(Db + il) * exc;
    *(Fxb + jl) -= *(Db + ik) * exc;
  };

  const auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  const auto maxDens = maxDensMat.maxCoeff();
  auto prescreen = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwartz) {
    if (maxDens * schwartz < _screening) {
      return true;
    }
    double maxDBlock = maxDensMat(i, k);
    maxDBlock = std::max(maxDBlock, maxDensMat(i, l));
    maxDBlock = std::max(maxDBlock, maxDensMat(j, k));
    maxDBlock = std::max(maxDBlock, maxDensMat(j, l));
    if (maxDBlock * schwartz < _screening) {
      return true;
    }
    return false;
  };

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, _incrementHelper->getPrescreeningThreshold());
  looper.loopNoDerivative(distribute, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);

  for (unsigned int i = 1; i < nThreads; ++i) {
    fx[0].alpha += fx[i].alpha;
    fx[0].beta += fx[i].beta;
  }

  // Symmetrize.
  Eigen::Ref<Eigen::MatrixXd> Fxa = fx[0].alpha;
  Eigen::Ref<Eigen::MatrixXd> Fxb = fx[0].beta;

  Eigen::MatrixXd tmp_x = Fxa + Fxa.transpose();
  F.alpha += tmp_x;

  tmp_x = Fxb + Fxb.transpose();
  F.beta += tmp_x;
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
