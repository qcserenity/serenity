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
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <algorithm>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HFPotential<SCFMode>::HFPotential(std::shared_ptr<SystemController> systemController,
                                  std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double xRatio,
                                  const double prescreeningThreshold, double prescreeningIncrementStart,
                                  double prescreeningIncrementEnd, unsigned int incrementSteps, bool clear4CenterCache)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _systemController(systemController),
    _xRatio(xRatio),
    _dMatController(dMat),
    _fullpotential(nullptr),
    _fullXpotential(nullptr),
    _outOfDate(true),
    _incrementHelper(std::make_shared<IncrementalFockMatrix<SCFMode>>(
        dMat, (prescreeningThreshold != 0) ? prescreeningThreshold : this->_basis->getPrescreeningThreshold(),
        prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps, "Coulomb and Exact-Exchange")),
    _clear4CenterCache(clear4CenterCache) {
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
FockMatrix<SCFMode>& HFPotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -         HF Pot.");
  if (_outOfDate) {
    DensityMatrix<SCFMode> densityMatrix(this->_basis);
    std::vector<std::shared_ptr<FockMatrix<SCFMode>>> matrices = {_fullpotential, _fullXpotential};
    _incrementHelper->updateDensityAndThreshold(densityMatrix, _screening, matrices);
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

  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>> fc;
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>> fx;
  const unsigned int nThreads = omp_get_max_threads();
  for (unsigned int i = 0; i < nThreads; ++i) {
    fc.emplace_back(_basis);
    fx.emplace_back(_basis);
  }

  auto distribute = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned long ij = i * nBFs + j;
    unsigned long kl = k * nBFs + l;
    unsigned long ik = i * nBFs + k;
    unsigned long il = i * nBFs + l;
    unsigned long jl = j * nBFs + l;
    unsigned long jk = j * nBFs + k;
    if (integral != integral)
      throw SerenityError("NAN detected in HF potential for 4-center integral");

    auto Fc = fc[threadId].data();
    auto Fx = fx[threadId].data();
    auto D = densityMatrix.data();

    /*
     * Restricted Hartree-Fock:
     * F(i,j) = [(ij|kl) - 0.5 * (ij|kl)] D(k,l)
     *
     * First, recall that we exploit the eigthfold integral symmetry in the looper
     * so we can expect the above contraction to be performed eight times.
     * However, both the density and Fock matrix are symmetric and so we can save
     * some time by only performing half of the actually needed operations and
     * recover the missing contribution by a symmetrization after the Fock build.
     * This is why we have four contractions for the exchange and Coulomb matrix
     * (actually only two because two of the four contractions are in turn
     * identical so we end up with a factor of 2.0).
     */

    // Coulomb.
    const double coul = 2.0 * integral;
    *(Fc + ij) += *(D + kl) * coul;
    *(Fc + kl) += *(D + ij) * coul;

    // Exchange.
    const double exc = 0.5 * _xRatio * integral;
    *(Fx + ik) -= *(D + jl) * exc;
    *(Fx + il) -= *(D + jk) * exc;
    *(Fx + jk) -= *(D + il) * exc;
    *(Fx + jl) -= *(D + ik) * exc;
  };

  auto distributeNoExc = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned long ij = i * nBFs + j;
    unsigned long kl = k * nBFs + l;
    if (integral != integral)
      throw SerenityError("NAN detected during Hartree-Fock Fock matrix build for 4-center integral");

    auto Fc = fc[threadId].data();
    auto D = densityMatrix.data();

    /*
     * Restricted Hartree-Fock Coulomb:
     * F(i,j) = [(ij|kl)] D(k,l)
     *
     * First, recall that we exploit the eigthfold integral symmetry in the looper
     * so we can expect the above contraction to be performed eight times.
     * However, both the density and Fock matrix are symmetric and so we can save
     * some time by only performing half of the actually needed operations and
     * recover the missing contribution by a symmetrization after the Fock build.
     * This is why we have four contractions for the exchange and Coulomb matrix
     * (actually only two because two of the four contractions are in turn
     * identical so we end up with a factor of 2.0).
     */

    // Coulomb.
    const double coul = 2.0 * integral;
    *(Fc + ij) += *(D + kl) * coul;
    *(Fc + kl) += *(D + ij) * coul;
  };

  const auto maxDensMat = densityMatrix.shellWiseAbsMax();
  const auto maxDens = maxDensMat.maxCoeff();
  auto prescreen = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwartz) {
    if (maxDens * schwartz < _screening) {
      return true;
    }
    double maxDBlock = 2.0 * maxDensMat(i, j);
    maxDBlock = std::max(maxDBlock, 2.0 * maxDensMat(k, l));
    maxDBlock = std::max(maxDBlock, 0.5 * maxDensMat(i, k));
    maxDBlock = std::max(maxDBlock, 0.5 * maxDensMat(i, l));
    maxDBlock = std::max(maxDBlock, 0.5 * maxDensMat(j, k));
    maxDBlock = std::max(maxDBlock, 0.5 * maxDensMat(j, l));
    if (maxDBlock * schwartz < _screening) {
      return true;
    }
    return false;
  };
  /*
   * We will use the prescreening threshold for the full matrix build to check the general Schwarz inequality.
   * The density matrix-based prescreening is not affected by this non-adapted threshold between SCF iterations.
   */
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, _incrementHelper->getPrescreeningThreshold());
  /*
   * Run
   */
  if (_xRatio != 0.0) {
    looper.loopNoDerivative(distribute, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);
  }
  else {
    looper.loopNoDerivative(distributeNoExc, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);
  }
  for (unsigned int i = 1; i < nThreads; ++i) {
    fc[0] += fc[i];
    fx[0] += fx[i];
  }

  // Symmetrize.
  Eigen::Ref<Eigen::MatrixXd> Fca = fc[0];
  Eigen::Ref<Eigen::MatrixXd> Fxa = fx[0];
  Eigen::MatrixXd tmp_c = Fca + Fca.transpose();
  Eigen::MatrixXd tmp_x = Fxa + Fxa.transpose();
  (*_fullXpotential) += tmp_x;
  F += tmp_c + tmp_x;
}

template<>
void HFPotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F,
                                                                const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix) {
  const unsigned int nBFs = _basis->getNBasisFunctions();

  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>> fx;
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>> fc;
  const unsigned int nThreads = omp_get_max_threads();
  for (unsigned int i = 0; i < nThreads; ++i) {
    fc.emplace_back(_basis);
    fx.emplace_back(_basis);
  }

  auto distribute = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned long ij = i * nBFs + j;
    unsigned long kl = k * nBFs + l;
    unsigned long ik = i * nBFs + k;
    unsigned long il = i * nBFs + l;
    unsigned long jl = j * nBFs + l;
    unsigned long jk = j * nBFs + k;

    auto Fca = fc[threadId].alpha.data();
    auto Fcb = fc[threadId].beta.data();
    auto Fxa = fx[threadId].alpha.data();
    auto Fxb = fx[threadId].beta.data();

    auto Da = densityMatrix.alpha.data();
    auto Db = densityMatrix.beta.data();

    // Coulomb.
    double coul1 = 2.0 * (*(Da + kl) + *(Db + kl)) * integral;
    double coul2 = 2.0 * (*(Da + ij) + *(Db + ij)) * integral;
    *(Fca + ij) += coul1;
    *(Fca + kl) += coul2;
    *(Fcb + ij) += coul1;
    *(Fcb + kl) += coul2;

    // Exchange.
    double exc = _xRatio * integral;
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
  const double maxDens = maxDensMat.maxCoeff();
  auto prescreen = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwartz) {
    if (maxDens * schwartz < _screening) {
      return true;
    }
    double maxDBlock = 2.0 * maxDensMat(i, j);
    maxDBlock = std::max(maxDBlock, 2.0 * maxDensMat(k, l));
    maxDBlock = std::max(maxDBlock, maxDensMat(i, k));
    maxDBlock = std::max(maxDBlock, maxDensMat(i, l));
    maxDBlock = std::max(maxDBlock, maxDensMat(j, k));
    maxDBlock = std::max(maxDBlock, maxDensMat(j, l));
    if (maxDBlock * schwartz < _screening) {
      return true;
    }
    return false;
  };

  /*
   * We will use the prescreening threshold for the full matrix build to check the general Schwarz inequality.
   * The density matrix-based prescreening is not affected by this non-adapted threshold between SCF iterations.
   */
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, _incrementHelper->getPrescreeningThreshold());
  looper.loopNoDerivative(distribute, prescreen, maxDens, _systemController.lock()->getIntegralCachingController(), true);

  for (unsigned int i = 1; i < nThreads; ++i) {
    fc[0].alpha += fc[i].alpha;
    fc[0].beta += fc[i].beta;
    fx[0].alpha += fx[i].alpha;
    fx[0].beta += fx[i].beta;
  }

  // Symmetrize.
  Eigen::Ref<Eigen::MatrixXd> Fca = fc[0].alpha;
  Eigen::Ref<Eigen::MatrixXd> Fxa = fx[0].alpha;
  Eigen::Ref<Eigen::MatrixXd> Fcb = fc[0].beta;
  Eigen::Ref<Eigen::MatrixXd> Fxb = fx[0].beta;

  Eigen::MatrixXd tmp_c = Fca + Fca.transpose();
  Eigen::MatrixXd tmp_x = Fxa + Fxa.transpose();
  _fullXpotential->alpha += tmp_x;
  F.alpha += tmp_c + tmp_x;

  tmp_c = Fcb + Fcb.transpose();
  tmp_x = Fxb + Fxb.transpose();
  _fullXpotential->beta += tmp_x;
  F.beta += tmp_c + tmp_x;
}

template<>
Eigen::MatrixXd HFPotential<RESTRICTED>::getGeomGradients() {
  auto systemController = _systemController.lock();
  auto atoms = systemController->getAtoms();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _dMatController->getDensityMatrix();
  auto dMatPtr = densityMatrix.data();

  auto mapping = systemController->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

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

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _basis, _incrementHelper->getPrescreeningThreshold());
  const auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  const double maxDens = maxDensMat.maxCoeff();
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

  looper.loop(looperFunction, maxDens);

  Eigen::MatrixXd hfGrad(nAtoms, 3);
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

  auto mapping = systemController->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

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

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _basis, _incrementHelper->getPrescreeningThreshold());
  const auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  const double maxDens = maxDensMat.maxCoeff();
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

  looper.loop(looperFunction, maxDens);

  Eigen::MatrixXd hfGrad(nAtoms, 3);
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

template<Options::SCF_MODES SCFMode>
HFPotential<SCFMode>::~HFPotential() {
  if (_systemController.lock() && _clear4CenterCache)
    _systemController.lock()->clear4CenterCache();
}

template class HFPotential<Options::SCF_MODES::RESTRICTED>;
template class HFPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
