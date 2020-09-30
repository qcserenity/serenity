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
#include "misc/SerenityError.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"

namespace Serenity {
template<Options::SCF_MODES SCFMode>
ExchangePotential<SCFMode>::ExchangePotential(std::shared_ptr<SystemController> systemController,
                                              std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double exchangeRatio,
                                              const double prescreeningThreshold, double prescreeningIncrementStart,
                                              double prescreeningIncrementEnd, unsigned int incrementSteps)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    IncrementalFockMatrix<SCFMode>(dMat, prescreeningThreshold, prescreeningIncrementStart, prescreeningIncrementEnd,
                                   incrementSteps, "Exact Exchange"),
    _systemController(systemController),
    _exc(exchangeRatio),
    _dMatController(dMat),
    _fullpotential(nullptr),
    _outOfDate(true) {
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
    this->updateDensityAndThreshold(densityMatrix, _screening, *_fullpotential);
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

  /*
   * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
   */
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fx;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  // Let the first thread use the Fock matrix directly
  for (unsigned int i = 0; i < nThreads; ++i) {
    fx.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
    fx[i]->setZero();
  }

#else
  fx.push_back(&F);
  fx[0]->setZero();
#endif
  /*
   * Function which parses the integrals
   */
  auto distribute = [&](const unsigned& i, const unsigned& j, const unsigned& k, const unsigned& l,
                        const double& integral, const unsigned& threadId) {
    /*
     * Permutations
     */
    double perm = 2.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    /*
     * Exchange
     */
    const double exc = perm * integral * 0.25 * _exc;
    *(fx[threadId]->data() + j * nBFs + l) -= *(densityMatrix.data() + i * nBFs + k) * exc;
    *(fx[threadId]->data() + j * nBFs + k) -= *(densityMatrix.data() + i * nBFs + l) * exc;
    *(fx[threadId]->data() + i * nBFs + l) -= *(densityMatrix.data() + j * nBFs + k) * exc;
    *(fx[threadId]->data() + i * nBFs + k) -= *(densityMatrix.data() + j * nBFs + l) * exc;
  };
  /*
   * Detailed prescreening function
   */
  // Maximum absolute value in densityMatrix
  const double maxDens = densityMatrix.lpNorm<Eigen::Infinity>();
  auto prescreen = [&](const unsigned& i, const unsigned& j, const unsigned& k, const unsigned& l, const unsigned& nI,
                       const unsigned& nJ, const unsigned& nK, const unsigned& nL, const double& schwartz) {
    /*
     * Early return for insignificance based on the largest element in the whole density matrix
     */
    if (maxDens * schwartz < _screening)
      return true;
    double maxDBlock = densityMatrix.block(i, k, nI, nK).lpNorm<Eigen::Infinity>();
    maxDBlock = std::max(maxDBlock, densityMatrix.block(i, l, nI, nL).lpNorm<Eigen::Infinity>());
    maxDBlock = std::max(maxDBlock, densityMatrix.block(j, k, nJ, nK).lpNorm<Eigen::Infinity>());
    maxDBlock = std::max(maxDBlock, densityMatrix.block(j, l, nJ, nL).lpNorm<Eigen::Infinity>());
    if (maxDBlock * schwartz < _screening)
      return true;
    /*
     * Shell quadruple is significant.
     */
    return false;
  };
  /*
   * Construct the looper, which loops over all integrals
   */
  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 0, _basis, _screening);
  /*
   * Run
   */
  looper.loopNoDerivative(distribute, prescreen);
#ifdef _OPENMP
  for (unsigned int i = 1; i < nThreads; ++i) {
    *fx[0] += *fx[i];
    delete fx[i];
  }
#endif
  // Symmetrizing
  Eigen::MatrixXd& f_x = *fx[0];
  Eigen::MatrixXd temp_x = f_x + f_x.transpose();
  F += temp_x;
}

template<>
void ExchangePotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F, const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix) {
  /*
   * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
   */
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
  /*
   * Function which parses the integrals
   */
  auto distribute = [&](const unsigned& i, const unsigned& j, const unsigned& k, const unsigned& l,
                        const double& integral, const unsigned& threadId) {
    /*
     * Permutations
     */
    double perm = 2.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;

    /*
     * Exchange
     */

    const double exc = perm * integral * 0.5 * _exc;

    const double exc1a = densityMatrix.alpha(i, k) * exc;
    const double exc2a = densityMatrix.alpha(i, l) * exc;
    const double exc3a = densityMatrix.alpha(j, k) * exc;
    const double exc4a = densityMatrix.alpha(j, l) * exc;
    fx[threadId]->alpha(j, l) -= exc1a;
    fx[threadId]->alpha(j, k) -= exc2a;
    fx[threadId]->alpha(i, l) -= exc3a;
    fx[threadId]->alpha(i, k) -= exc4a;
    const double exc1b = densityMatrix.beta(i, k) * exc;
    const double exc2b = densityMatrix.beta(i, l) * exc;
    const double exc3b = densityMatrix.beta(j, k) * exc;
    const double exc4b = densityMatrix.beta(j, l) * exc;
    fx[threadId]->beta(j, l) -= exc1b;
    fx[threadId]->beta(j, k) -= exc2b;
    fx[threadId]->beta(i, l) -= exc3b;
    fx[threadId]->beta(i, k) -= exc4b;
  };
  /*
   * Detailed prescreening function
   */
  // Maximum absolute value in densityMatrix
  const double maxDens =
      std::max(densityMatrix.alpha.lpNorm<Eigen::Infinity>(), densityMatrix.beta.lpNorm<Eigen::Infinity>());
  auto prescreen = [&](const unsigned& i, const unsigned& j, const unsigned& k, const unsigned& l, const unsigned& nI,
                       const unsigned& nJ, const unsigned& nK, const unsigned& nL, const double& schwartz) {
    /*
     * Early return for insignificance based on the largest element in the whole density matrix
     */
    if (maxDens * schwartz < _screening)
      return true;
    double maxDBlock = 0.0;
    for_spin(densityMatrix) {
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(i, k, nI, nK).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(i, l, nI, nL).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(j, k, nJ, nK).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(j, l, nJ, nL).lpNorm<Eigen::Infinity>());
    };
    if (maxDBlock * schwartz < _screening)
      return true;
    /*
     * Shell quadruple is significant.
     */
    return false;
  };
  /*
   * Construct looper for loop over all integrals
   */
  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 0, _basis, _screening);
  /*
   * Run
   */
  looper.loopNoDerivative(distribute, prescreen);
#ifdef _OPENMP
  for (unsigned int i = 1; i < nThreads; ++i) {
    fx[0]->alpha += fx[i]->alpha;
    fx[0]->beta += fx[i]->beta;
  }
#endif
  // Symmetrizing
  Eigen::MatrixXd& f_x_a = fx[0]->alpha;
  Eigen::MatrixXd& f_x_b = fx[0]->beta;
  Eigen::MatrixXd temp_x = f_x_a + f_x_a.transpose();
  F.alpha += temp_x;
  temp_x = f_x_b + f_x_b.transpose();
  F.beta += temp_x;
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXi ExchangePotential<SCFMode>::createBasisToAtomMap(std::shared_ptr<Serenity::AtomCenteredBasisController> basis) {
  const unsigned int nBasisFunctions = basis->getNBasisFunctions();
  const unsigned int nBasisFunctionsRed = basis->getReducedNBasisFunctions();
  const auto basisIndicesRed = basis->getBasisIndicesRed();
  Eigen::VectorXi shellMap(nBasisFunctionsRed);
  Eigen::VectorXi bfMap(nBasisFunctions);
  // Vector to check whether ALL basis function shells are assigned to an atom index
  std::vector<bool> hasElementBeenSet(nBasisFunctionsRed, false);
  for (unsigned int iAtom = 0; iAtom < basisIndicesRed.size(); ++iAtom) {
    const unsigned int firstIndex = basisIndicesRed[iAtom].first;
    const unsigned int endIndex = basisIndicesRed[iAtom].second;
    for (unsigned int iShell = firstIndex; iShell < endIndex; ++iShell) {
      shellMap[iShell] = iAtom;
      hasElementBeenSet[iShell] = true;
    }
  }
  // Check
  for (bool x : hasElementBeenSet) {
    if (not x)
      throw SerenityError("ExchangePotential: Missed gradient element in gradient evaluation.");
  }
  for (unsigned int i = 0; i < nBasisFunctions; i++) {
    bfMap[i] = shellMap[basis->reducedIndex(i)];
  }
  return bfMap;
}

template<>
Eigen::MatrixXd ExchangePotential<RESTRICTED>::getGeomGradients() {
  auto systemController = _systemController.lock();
  auto atoms = systemController->getAtoms();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _dMatController->getDensityMatrix();
  auto dMatPtr = densityMatrix.data();

  auto mapping = createBasisToAtomMap(systemController->getAtomCenteredBasisController());

  const unsigned int nBFs = _basis->getNBasisFunctions();

#ifdef _OPENMP
  std::vector<Eigen::MatrixXd> priv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
#else
  std::vector<Eigen::MatrixXd> priv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
#endif

  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 1, _basis, _screening);

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

  looper.loop(looperFunction);

  Matrix<double> xGrad(nAtoms, 3);
  xGrad.setZero();
#ifdef _OPENMP
  // sum over all threads
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    xGrad -= priv[i] * _exc;
  }
#else
  xGrad -= priv[0] * _exc;
#endif
  return std::move(xGrad);
}

template<>
Eigen::MatrixXd ExchangePotential<UNRESTRICTED>::getGeomGradients() {
  auto systemController = _systemController.lock();
  auto atoms = systemController->getAtoms();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _dMatController->getDensityMatrix();
  auto dMatPtrA = densityMatrix.alpha.data();
  auto dMatPtrB = densityMatrix.beta.data();

  auto mapping = createBasisToAtomMap(systemController->getAtomCenteredBasisController());

  const unsigned int nBFs = _basis->getNBasisFunctions();

#ifdef _OPENMP
  std::vector<Eigen::MatrixXd> priv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
#else
  std::vector<Eigen::MatrixXd> priv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
#endif

  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 1, _basis, _screening);

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

  looper.loop(looperFunction);

  Matrix<double> xGrad(nAtoms, 3);
  xGrad.setZero();
#ifdef _OPENMP
  // sum over all threads
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    xGrad -= priv[i] * _exc;
  }
#else
  xGrad -= priv[0] * _exc;
#endif
  return std::move(xGrad);
}

template class ExchangePotential<Options::SCF_MODES::RESTRICTED>;
template class ExchangePotential<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
