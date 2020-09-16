/**
 * @file   HFPotential.cpp
 *
 * @date   Nov 24, 2016
 * @author Jan Unsleber
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
#include "potentials/HFPotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "io/FormattedOutputStream.h" //Filtered output streams.
#include "misc/Timing.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <algorithm>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HFPotential<SCFMode>::HFPotential(std::shared_ptr<SystemController> systemController,
                                  std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double xRatio,
                                  const double prescreeningThreshold, double prescreeningIncrementStart,
                                  double prescreeningIncrementEnd, unsigned int incrementSteps)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    IncrementalFockMatrix<SCFMode>(dMat, prescreeningThreshold, prescreeningIncrementStart, prescreeningIncrementEnd,
                                   incrementSteps),
    _systemController(systemController),
    _xRatio(xRatio),
    _dMatController(dMat),
    _fullpotential(nullptr),
    _fullXpotential(nullptr),
    _outOfDate(true) {
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);

  _fullpotential = std::make_shared<FockMatrix<SCFMode>>(FockMatrix<SCFMode>(this->_basis));
  auto& temp = *_fullpotential;
  for_spin(temp) {
    temp_spin.setZero();
  };
  _fullXpotential = std::make_shared<FockMatrix<SCFMode>>(FockMatrix<SCFMode>(this->_basis));
  auto& temp2 = *_fullXpotential;
  for_spin(temp2) {
    temp2_spin.setZero();
  };
  _screening = prescreeningIncrementStart;
};

template<Options::SCF_MODES SCFMode>
void HFPotential<SCFMode>::resetFockMatrix(FockMatrix<SCFMode>& f, double nextTreshold) {
  // No printing in the fist cycle.
  if (this->getCounter() != 0) {
    OutputControl::dOut << " ***** Reset Incremental Fock Matrix Build: Coulomb and Exact-Exchange *****" << std::endl;
    OutputControl::dOut << " ***** New Prescreening Threshold - " << nextTreshold << " ***** " << std::endl;
  }
  auto& xPot = *_fullXpotential;
  for_spin(f, xPot) {
    f_spin.setZero();
    xPot_spin.setZero();
  };
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& HFPotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -         HF Pot.");
  if (_outOfDate) {
    DensityMatrix<SCFMode> densityMatrix(this->_basis);
    this->updateDensityAndThreshold(densityMatrix, _screening, *_fullpotential);
    this->addToMatrix(*_fullpotential, densityMatrix);
    _outOfDate = false;
  }
  Timings::timeTaken("Active System -         HF Pot.");
  return *_fullpotential;
}

template<Options::SCF_MODES SCFMode>
double HFPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (_outOfDate)
    this->getMatrix();
  Timings::takeTime("Active System -         HF Pot.");
  auto& pot = *_fullpotential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += 0.5 * pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("Active System -         HF Pot.");
  return energy;
};

template<Options::SCF_MODES SCFMode>
double HFPotential<SCFMode>::getXEnergy(const DensityMatrix<SCFMode>& P) {
  if (_outOfDate)
    this->getMatrix();
  auto& pot = *_fullXpotential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += 0.5 * pot_spin.cwiseProduct(P_spin).sum();
  };
  return energy;
};

template<>
void HFPotential<Options::SCF_MODES::RESTRICTED>::addToMatrix(FockMatrix<Options::SCF_MODES::RESTRICTED>& F,
                                                              const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densityMatrix) {
  const unsigned int nBFs = _basis->getNBasisFunctions();

  /*
   * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
   */
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fc;
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fx;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  // Let the first thread use the Fock matrix directly
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
     * Coulomb
     */
    const double coul = perm * integral;
    *(fc[threadId]->data() + i * nBFs + j) += *(densityMatrix.data() + k * nBFs + l) * coul;
    *(fc[threadId]->data() + k * nBFs + l) += *(densityMatrix.data() + i * nBFs + j) * coul;
    /*
     * Exchange
     */
    const double exc = perm * integral * 0.25 * _xRatio;
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
    double maxDBlock = densityMatrix.block(i, j, nI, nJ).lpNorm<Eigen::Infinity>();
    maxDBlock = std::max(maxDBlock, densityMatrix.block(i, k, nI, nK).lpNorm<Eigen::Infinity>());
    maxDBlock = std::max(maxDBlock, densityMatrix.block(i, l, nI, nL).lpNorm<Eigen::Infinity>());
    maxDBlock = std::max(maxDBlock, densityMatrix.block(j, k, nJ, nK).lpNorm<Eigen::Infinity>());
    maxDBlock = std::max(maxDBlock, densityMatrix.block(j, l, nJ, nL).lpNorm<Eigen::Infinity>());
    maxDBlock = std::max(maxDBlock, densityMatrix.block(k, l, nK, nL).lpNorm<Eigen::Infinity>());
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
    *fc[0] += *fc[i];
    *fx[0] += *fx[i];
    delete fx[i];
    delete fc[i];
  }
#endif
  // Symmetrizing
  Eigen::MatrixXd& f_x = *fx[0];
  Eigen::MatrixXd& f_c = *fc[0];
  Eigen::MatrixXd temp_c = f_c + f_c.transpose();
  Eigen::MatrixXd temp_x = f_x + f_x.transpose();
  *_fullXpotential += temp_x;
  F += temp_c + temp_x;
}

template<>
void HFPotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F,
                                                                const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix) {
  /*
   * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
   */
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>*> fx;
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>*> fc;
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
  // Simply use the fock matrix directly
  fc.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
  fc[0]->alpha.setZero();
  fc[0]->beta.setZero();
  fx.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
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
     * Coulomb
     */
    const double coul = perm * integral;
    const double coul1 = (densityMatrix.alpha(k, l) + densityMatrix.beta(k, l)) * coul;
    const double coul2 = (densityMatrix.alpha(i, j) + densityMatrix.beta(i, j)) * coul;
    fc[threadId]->beta(i, j) += coul1;
    fc[threadId]->beta(k, l) += coul2;
    fc[threadId]->alpha(i, j) += coul1;
    fc[threadId]->alpha(k, l) += coul2;
    /*
     * Exchange
     */
    const double exc = perm * integral * 0.5 * _xRatio;
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
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(i, j, nI, nJ).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(i, k, nI, nK).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(i, l, nI, nL).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(j, k, nJ, nK).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(j, l, nJ, nL).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(k, l, nK, nL).lpNorm<Eigen::Infinity>());
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
    fc[0]->alpha += fc[i]->alpha;
    fc[0]->beta += fc[i]->beta;
    fx[0]->alpha += fx[i]->alpha;
    fx[0]->beta += fx[i]->beta;
  }
#endif
  // Symmetrizing
  Eigen::MatrixXd& f_x_a = fx[0]->alpha;
  Eigen::MatrixXd& f_x_b = fx[0]->beta;
  Eigen::MatrixXd& f_c_a = fc[0]->alpha;
  Eigen::MatrixXd& f_c_b = fc[0]->beta;

  Eigen::MatrixXd temp_c = f_c_a + f_c_a.transpose();
  Eigen::MatrixXd temp_x = f_x_a + f_x_a.transpose();
  F.alpha += temp_c + temp_x;
  _fullXpotential->alpha += temp_x;

  temp_c = f_c_b + f_c_b.transpose();
  temp_x = f_x_b + f_x_b.transpose();
  F.beta += temp_c + temp_x;
  _fullXpotential->beta += temp_x;
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXi HFPotential<SCFMode>::createBasisToAtomMap(std::shared_ptr<Serenity::AtomCenteredBasisController> basis) {
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
      throw SerenityError("HFPotential: Missed gradient element in gradient evaluation.");
  }
  for (unsigned int i = 0; i < nBasisFunctions; i++) {
    bfMap[i] = shellMap[basis->reducedIndex(i)];
  }
  return bfMap;
}

template<>
Eigen::MatrixXd HFPotential<RESTRICTED>::getGeomGradients() {
  auto systemController = _systemController.lock();
  auto atoms = systemController->getAtoms();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _dMatController->getDensityMatrix();
  auto dMatPtr = densityMatrix.data();

  auto mapping = createBasisToAtomMap(systemController->getAtomCenteredBasisController());

  const unsigned int nBFs = _basis->getNBasisFunctions();

#ifdef _OPENMP
  // create a vector of matrices for each thread
  std::vector<Eigen::MatrixXd> cPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
  std::vector<Eigen::MatrixXd> xPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
#else
  // or just one
  std::vector<Eigen::MatrixXd> cPriv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
  std::vector<Eigen::MatrixXd> xPriv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
#endif

  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 1, _basis, 1E-10);

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
    // Coulomb
    const double cprefac = dMatPtr[i * nBFs + j] * dMatPtr[a * nBFs + b];
    cPriv[threadID].data()[0 * nAtoms + iAtom] += perm * cprefac * intValues.data()[0];  // 0 * 3 + 0 = (at, xyz)
    cPriv[threadID].data()[1 * nAtoms + iAtom] += perm * cprefac * intValues.data()[1];  // 0 * 3 + 1 = (at, xyz)
    cPriv[threadID].data()[2 * nAtoms + iAtom] += perm * cprefac * intValues.data()[2];  // 0 * 3 + 2 = (at, xyz)
    cPriv[threadID].data()[0 * nAtoms + jAtom] += perm * cprefac * intValues.data()[3];  // 1 * 3 + 0 = (at, xyz)
    cPriv[threadID].data()[1 * nAtoms + jAtom] += perm * cprefac * intValues.data()[4];  // 1 * 3 + 1 = (at, xyz)
    cPriv[threadID].data()[2 * nAtoms + jAtom] += perm * cprefac * intValues.data()[5];  // 1 * 3 + 2 = (at, xyz)
    cPriv[threadID].data()[0 * nAtoms + aAtom] += perm * cprefac * intValues.data()[6];  // 2 * 3 + 0 = (at, xyz)
    cPriv[threadID].data()[1 * nAtoms + aAtom] += perm * cprefac * intValues.data()[7];  // 2 * 3 + 1 = (at, xyz)
    cPriv[threadID].data()[2 * nAtoms + aAtom] += perm * cprefac * intValues.data()[8];  // 2 * 3 + 2 = (at, xyz)
    cPriv[threadID].data()[0 * nAtoms + bAtom] += perm * cprefac * intValues.data()[9];  // 3 * 3 + 0 = (at, xyz)
    cPriv[threadID].data()[1 * nAtoms + bAtom] += perm * cprefac * intValues.data()[10]; // 3 * 3 + 1 = (at, xyz)
    cPriv[threadID].data()[2 * nAtoms + bAtom] += perm * cprefac * intValues.data()[11]; // 3 * 3 + 2 = (at, xyz)
    // Exchange
    const double xprefac = dMatPtr[i * nBFs + a] * dMatPtr[j * nBFs + b] + dMatPtr[i * nBFs + b] * dMatPtr[j * nBFs + a];
    xPriv[threadID].data()[0 * nAtoms + iAtom] += 0.25 * perm * xprefac * intValues.data()[0];  // 0 * 3 + 0 = (at, xyz)
    xPriv[threadID].data()[1 * nAtoms + iAtom] += 0.25 * perm * xprefac * intValues.data()[1];  // 0 * 3 + 1 = (at, xyz)
    xPriv[threadID].data()[2 * nAtoms + iAtom] += 0.25 * perm * xprefac * intValues.data()[2];  // 0 * 3 + 2 = (at, xyz)
    xPriv[threadID].data()[0 * nAtoms + jAtom] += 0.25 * perm * xprefac * intValues.data()[3];  // 1 * 3 + 0 = (at, xyz)
    xPriv[threadID].data()[1 * nAtoms + jAtom] += 0.25 * perm * xprefac * intValues.data()[4];  // 1 * 3 + 1 = (at, xyz)
    xPriv[threadID].data()[2 * nAtoms + jAtom] += 0.25 * perm * xprefac * intValues.data()[5];  // 1 * 3 + 2 = (at, xyz)
    xPriv[threadID].data()[0 * nAtoms + aAtom] += 0.25 * perm * xprefac * intValues.data()[6];  // 2 * 3 + 0 = (at, xyz)
    xPriv[threadID].data()[1 * nAtoms + aAtom] += 0.25 * perm * xprefac * intValues.data()[7];  // 2 * 3 + 1 = (at, xyz)
    xPriv[threadID].data()[2 * nAtoms + aAtom] += 0.25 * perm * xprefac * intValues.data()[8];  // 2 * 3 + 2 = (at, xyz)
    xPriv[threadID].data()[0 * nAtoms + bAtom] += 0.25 * perm * xprefac * intValues.data()[9];  // 3 * 3 + 0 = (at, xyz)
    xPriv[threadID].data()[1 * nAtoms + bAtom] += 0.25 * perm * xprefac * intValues.data()[10]; // 3 * 3 + 1 = (at, xyz)
    xPriv[threadID].data()[2 * nAtoms + bAtom] += 0.25 * perm * xprefac * intValues.data()[11]; // 3 * 3 + 2 = (at, xyz)
  };

  looper.loop(looperFunction);

  Matrix<double> hfGrad(nAtoms, 3);
  hfGrad.setZero();
#ifdef _OPENMP
  // sum over all threads
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    hfGrad += cPriv[i];
    hfGrad -= xPriv[i] * _xRatio;
  }
#else
  hfGrad += cPriv[0];
  hfGrad -= xPriv[0] * _xRatio;
#endif
  return hfGrad;
}

template<>
Eigen::MatrixXd HFPotential<UNRESTRICTED>::getGeomGradients() {
  auto systemController = _systemController.lock();
  auto atoms = systemController->getAtoms();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _dMatController->getDensityMatrix();
  auto dMatPtrA = densityMatrix.alpha.data();
  auto dMatPtrB = densityMatrix.beta.data();

  auto mapping = createBasisToAtomMap(systemController->getAtomCenteredBasisController());

  const unsigned int nBFs = _basis->getNBasisFunctions();

#ifdef _OPENMP
  // create a vector of matrices for each thread
  std::vector<Eigen::MatrixXd> cPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
  std::vector<Eigen::MatrixXd> xPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
#else
  // or just one
  std::vector<Eigen::MatrixXd> cPriv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
  std::vector<Eigen::MatrixXd> xPriv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
#endif

  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 1, _basis, 1E-10);

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
    // Coulomb
    const double cprefac =
        dMatPtrA[i * nBFs + j] * dMatPtrA[a * nBFs + b] + dMatPtrB[i * nBFs + j] * dMatPtrB[a * nBFs + b] +
        dMatPtrB[i * nBFs + j] * dMatPtrA[a * nBFs + b] + dMatPtrA[i * nBFs + j] * dMatPtrB[a * nBFs + b];
    cPriv[threadID].data()[0 * nAtoms + iAtom] += perm * cprefac * intValues.data()[0];  // 0 * 3 + 0 = (at, xyz)
    cPriv[threadID].data()[1 * nAtoms + iAtom] += perm * cprefac * intValues.data()[1];  // 0 * 3 + 1 = (at, xyz)
    cPriv[threadID].data()[2 * nAtoms + iAtom] += perm * cprefac * intValues.data()[2];  // 0 * 3 + 2 = (at, xyz)
    cPriv[threadID].data()[0 * nAtoms + jAtom] += perm * cprefac * intValues.data()[3];  // 1 * 3 + 0 = (at, xyz)
    cPriv[threadID].data()[1 * nAtoms + jAtom] += perm * cprefac * intValues.data()[4];  // 1 * 3 + 1 = (at, xyz)
    cPriv[threadID].data()[2 * nAtoms + jAtom] += perm * cprefac * intValues.data()[5];  // 1 * 3 + 2 = (at, xyz)
    cPriv[threadID].data()[0 * nAtoms + aAtom] += perm * cprefac * intValues.data()[6];  // 2 * 3 + 0 = (at, xyz)
    cPriv[threadID].data()[1 * nAtoms + aAtom] += perm * cprefac * intValues.data()[7];  // 2 * 3 + 1 = (at, xyz)
    cPriv[threadID].data()[2 * nAtoms + aAtom] += perm * cprefac * intValues.data()[8];  // 2 * 3 + 2 = (at, xyz)
    cPriv[threadID].data()[0 * nAtoms + bAtom] += perm * cprefac * intValues.data()[9];  // 3 * 3 + 0 = (at, xyz)
    cPriv[threadID].data()[1 * nAtoms + bAtom] += perm * cprefac * intValues.data()[10]; // 3 * 3 + 1 = (at, xyz)
    cPriv[threadID].data()[2 * nAtoms + bAtom] += perm * cprefac * intValues.data()[11]; // 3 * 3 + 2 = (at, xyz)
    // Exchange
    const double xprefac =
        dMatPtrA[i * nBFs + a] * dMatPtrA[j * nBFs + b] + dMatPtrA[i * nBFs + b] * dMatPtrA[j * nBFs + a] +
        dMatPtrB[i * nBFs + a] * dMatPtrB[j * nBFs + b] + dMatPtrB[i * nBFs + b] * dMatPtrB[j * nBFs + a];
    xPriv[threadID].data()[0 * nAtoms + iAtom] += 0.5 * perm * xprefac * intValues.data()[0];  // 0 * 3 + 0 = (at, xyz)
    xPriv[threadID].data()[1 * nAtoms + iAtom] += 0.5 * perm * xprefac * intValues.data()[1];  // 0 * 3 + 1 = (at, xyz)
    xPriv[threadID].data()[2 * nAtoms + iAtom] += 0.5 * perm * xprefac * intValues.data()[2];  // 0 * 3 + 2 = (at, xyz)
    xPriv[threadID].data()[0 * nAtoms + jAtom] += 0.5 * perm * xprefac * intValues.data()[3];  // 1 * 3 + 0 = (at, xyz)
    xPriv[threadID].data()[1 * nAtoms + jAtom] += 0.5 * perm * xprefac * intValues.data()[4];  // 1 * 3 + 1 = (at, xyz)
    xPriv[threadID].data()[2 * nAtoms + jAtom] += 0.5 * perm * xprefac * intValues.data()[5];  // 1 * 3 + 2 = (at, xyz)
    xPriv[threadID].data()[0 * nAtoms + aAtom] += 0.5 * perm * xprefac * intValues.data()[6];  // 2 * 3 + 0 = (at, xyz)
    xPriv[threadID].data()[1 * nAtoms + aAtom] += 0.5 * perm * xprefac * intValues.data()[7];  // 2 * 3 + 1 = (at, xyz)
    xPriv[threadID].data()[2 * nAtoms + aAtom] += 0.5 * perm * xprefac * intValues.data()[8];  // 2 * 3 + 2 = (at, xyz)
    xPriv[threadID].data()[0 * nAtoms + bAtom] += 0.5 * perm * xprefac * intValues.data()[9];  // 3 * 3 + 0 = (at, xyz)
    xPriv[threadID].data()[1 * nAtoms + bAtom] += 0.5 * perm * xprefac * intValues.data()[10]; // 3 * 3 + 1 = (at, xyz)
    xPriv[threadID].data()[2 * nAtoms + bAtom] += 0.5 * perm * xprefac * intValues.data()[11]; // 3 * 3 + 2 = (at, xyz)
  };

  looper.loop(looperFunction);

  Matrix<double> hfGrad(nAtoms, 3);
  hfGrad.setZero();
#ifdef _OPENMP
  // sum over all threads
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    hfGrad += cPriv[i];
    hfGrad -= xPriv[i] * _xRatio;
  }
#else
  hfGrad += cPriv[0];
  hfGrad -= xPriv[0] * _xRatio;
#endif
  return hfGrad;
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& HFPotential<SCFMode>::getXPotential() {
  if (_outOfDate)
    this->getMatrix();
  return *_fullXpotential;
}

template class HFPotential<Options::SCF_MODES::RESTRICTED>;
template class HFPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
