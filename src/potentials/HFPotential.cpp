/**
 * @file   HFPotential.cpp
 *
 * @date   Nov 24, 2016
 * @author Jan Unsleber
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
#include "potentials/HFPotential.h"
/* Include Serenity Internal Headers */
#include "data/matrices/FockMatrix.h"
#include "misc/Timing.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
/* Include Std and External Headers */
#include <algorithm>


namespace Serenity {



template <Options::SCF_MODES SCFMode>
HFPotential<SCFMode>::HFPotential(std::shared_ptr<SystemController> systemController,
    std::shared_ptr<DensityMatrixController<SCFMode> > dMat,
    const double exchangeRatio,
    const double prescreeningThreshold):
    Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _systemController(systemController),
    _prescreeningThreshold(prescreeningThreshold),
    _exc(exchangeRatio),
    _dMatController(dMat),
    _potential(nullptr),
    _excPotential(nullptr){
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode> >::_self);
};

template <Options::SCF_MODES SCFMode> FockMatrix<SCFMode>&
HFPotential<SCFMode>::getMatrix(){
  Timings::takeTime("Active System -   Coul./XC Pot.");
  if (!_potential or !_excPotential){
    const auto& densityMatrix = _dMatController->getDensityMatrix();
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot){
      pot_spin.setZero();
    };
    _excPotential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot2 = *_excPotential;
    for_spin(pot2){
      pot2_spin.setZero();
    };
    this->addToMatrix(*_potential,densityMatrix);
  }
  Timings::timeTaken("Active System -   Coul./XC Pot.");
  return *_potential;
}

template <Options::SCF_MODES SCFMode> double
HFPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P){
  if (!_potential) this->getMatrix();
  Timings::takeTime("Active System -   Coul./XC Pot.");
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot,P){
    energy += 0.5*pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("Active System -   Coul./XC Pot.");
  return energy;
};

template <Options::SCF_MODES SCFMode> double
HFPotential<SCFMode>::getXEnergy(const DensityMatrix<SCFMode>& P){
  if (!_excPotential) this->getMatrix();
  auto& pot = *_excPotential;
  double energy = 0.0;
  for_spin(pot,P){
    energy += 0.5*pot_spin.cwiseProduct(P_spin).sum();
  };
  return energy;
};

template <> void
HFPotential<Options::SCF_MODES::RESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::RESTRICTED>& F,
    const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densityMatrix){

    const unsigned int nBFs = _basis->getNBasisFunctions();

    /*
     * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
     */
    std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fc;
    std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fx;
#ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    // Let the first thread use the Fock matrix directly
    fc.push_back(&F);
    fx.push_back(_excPotential.get());
    for (unsigned int i=1; i<nThreads; ++i) {
      fc.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
      fc[i]->setZero();
      fx.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
      fx[i]->setZero();
    }

#else
    fc.push_back(&F);
    fx.push_back(_excPotential.get());
#endif
    /*
     * Function which parses the integrals
     */
    auto distribute = [&](
        unsigned int i, unsigned int j, unsigned int k, unsigned int l,
        const Eigen::VectorXd integral,
        unsigned int threadId) {

      /*
       * Permutations
       */
      double perm =2.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;
      perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;

      /*
       * Coulomb
       */
      const double coul = perm * integral[0];
      const double coul1 = *(densityMatrix.data()+k*nBFs+l)  * coul;
      const double coul2 = *(densityMatrix.data()+i*nBFs+j)  * coul;
      *(fc[threadId]->data()+i*nBFs+j) += coul1;
      *(fc[threadId]->data()+j*nBFs+i) += coul1;
      *(fc[threadId]->data()+k*nBFs+l) += coul2;
      *(fc[threadId]->data()+l*nBFs+k) += coul2;
      /*
       * Exchange
       */
      const double exc = perm * integral[0] * 0.25 * _exc;
      const double exc1 = *(densityMatrix.data()+i*nBFs+k)  * exc;
      const double exc2 = *(densityMatrix.data()+i*nBFs+l)  * exc;
      const double exc3 = *(densityMatrix.data()+j*nBFs+k)  * exc;
      const double exc4 = *(densityMatrix.data()+j*nBFs+l)  * exc;
      *(fx[threadId]->data()+j*nBFs+l) -= exc1;
      *(fx[threadId]->data()+l*nBFs+j) -= exc1;
      *(fx[threadId]->data()+j*nBFs+k) -= exc2;
      *(fx[threadId]->data()+k*nBFs+j) -= exc2;
      *(fx[threadId]->data()+i*nBFs+l) -= exc3;
      *(fx[threadId]->data()+l*nBFs+i) -= exc3;
      *(fx[threadId]->data()+i*nBFs+k) -= exc4;
      *(fx[threadId]->data()+k*nBFs+i) -= exc4;

    };
    /*
     * Detailed prescreening function
     */
    // Maximum absolute value in densityMatrix
    const double maxDens = densityMatrix.lpNorm<Eigen::Infinity>();
    auto prescreeningFunc = [&](
        unsigned int i, unsigned int j, unsigned int k, unsigned int l,
        unsigned int nI, unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
      /*
       * Early return for insignificance based on the largest element in the whole density matrix
       */
      if (maxDens*schwartz < _prescreeningThreshold) return true;
      double maxDBlock = densityMatrix.block(i,j,nI,nJ).lpNorm<Eigen::Infinity>();
      maxDBlock = std::max(maxDBlock, densityMatrix.block(i,k,nI,nK).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix.block(i,l,nI,nL).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix.block(j,k,nJ,nK).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix.block(j,l,nJ,nL).lpNorm<Eigen::Infinity>());
      maxDBlock = std::max(maxDBlock, densityMatrix.block(k,l,nK,nL).lpNorm<Eigen::Infinity>());
      if (maxDBlock*schwartz < _prescreeningThreshold) return true;
      /*
       * Shell quadruple is significant.
       */
      return false;
    };
    /*
     * Construct the looper, which loops over all integrals
     */
    TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb,0,_basis, _prescreeningThreshold);
    /*
     * Run
     */
    looper.loop(distribute, prescreeningFunc);
#ifdef _OPENMP
    for (unsigned int i=1; i<nThreads; ++i) {
      *_potential += *fc[i];
      *_excPotential += *fx[i];
      delete fx[i] ;
      delete fc[i] ;
    }
    *_potential += *_excPotential;
#endif

}

template <> void
HFPotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F,
    const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix){
    /*
     * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
     */
    std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>* > fx;
    std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>* > fc;
#ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    // Let the first thread use the Fock matrix directly
    fc.push_back(&F);
    fx.push_back(_excPotential.get());
    for (unsigned int i=1; i<nThreads; ++i) {
      fc.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
      fc[i]->alpha.setZero();
      fc[i]->beta.setZero();
      fx.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
      fx[i]->alpha.setZero();
      fx[i]->beta.setZero();
    }

#else
    // Simply use the fock matrix directly
    fc.push_back(&F);
    fx.push_back(_excPotential.get());
#endif
    /*
     * Function which parses the integrals
     */
    auto distribute = [&](
        unsigned int i, unsigned int j, unsigned int k, unsigned int l,
        const Eigen::VectorXd integral,
        unsigned int threadId) {


      /*
       * Permutations
       */
      double perm =2.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;
      perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;

      /*
       * Coulomb
       */
      const double coul = perm * integral[0];
      const double coul1 = (densityMatrix.alpha(k,l)+densityMatrix.beta(k,l))  * coul;
      const double coul2 = (densityMatrix.alpha(i,j)+densityMatrix.beta(i,j))  * coul;
      fc[threadId]->beta(i,j) += coul1;
      fc[threadId]->beta(j,i) += coul1;
      fc[threadId]->beta(k,l) += coul2;
      fc[threadId]->beta(l,k) += coul2;
      fc[threadId]->alpha(i,j) += coul1;
      fc[threadId]->alpha(j,i) += coul1;
      fc[threadId]->alpha(k,l) += coul2;
      fc[threadId]->alpha(l,k) += coul2;
      /*
       * Exchange
       */

      const double exc = perm * integral[0] * 0.5 * _exc;

      const double exc1a = densityMatrix.alpha(i,k)  * exc;
      const double exc2a = densityMatrix.alpha(i,l)  * exc;
      const double exc3a = densityMatrix.alpha(j,k)  * exc;
      const double exc4a = densityMatrix.alpha(j,l)  * exc;
      fx[threadId]->alpha(j,l) -= exc1a;
      fx[threadId]->alpha(l,j) -= exc1a;
      fx[threadId]->alpha(j,k) -= exc2a;
      fx[threadId]->alpha(k,j) -= exc2a;
      fx[threadId]->alpha(i,l) -= exc3a;
      fx[threadId]->alpha(l,i) -= exc3a;
      fx[threadId]->alpha(i,k) -= exc4a;
      fx[threadId]->alpha(k,i) -= exc4a;
      const double exc1b = densityMatrix.beta(i,k)  * exc;
      const double exc2b = densityMatrix.beta(i,l)  * exc;
      const double exc3b = densityMatrix.beta(j,k)  * exc;
      const double exc4b = densityMatrix.beta(j,l)  * exc;
      fx[threadId]->beta(j,l) -= exc1b;
      fx[threadId]->beta(l,j) -= exc1b;
      fx[threadId]->beta(j,k) -= exc2b;
      fx[threadId]->beta(k,j) -= exc2b;
      fx[threadId]->beta(i,l) -= exc3b;
      fx[threadId]->beta(l,i) -= exc3b;
      fx[threadId]->beta(i,k) -= exc4b;
      fx[threadId]->beta(k,i) -= exc4b;

    };
    /*
     * Detailed prescreening function
     */
    // Maximum absolute value in densityMatrix
    const double maxDens = std::max(
        densityMatrix.alpha.lpNorm<Eigen::Infinity>(),
        densityMatrix.beta.lpNorm<Eigen::Infinity>());
    auto prescreeningFunc = [&](
        unsigned int i, unsigned int j, unsigned int k, unsigned int l,
        unsigned int nI, unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
      /*
       * Early return for insignificance based on the largest element in the whole density matrix
       */
      if (maxDens*schwartz < _prescreeningThreshold) return true;
      double maxDBlock=0.0;
      for_spin(densityMatrix) {
        maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(i,j,nI,nJ).lpNorm<Eigen::Infinity>());
        maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(i,k,nI,nK).lpNorm<Eigen::Infinity>());
        maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(i,l,nI,nL).lpNorm<Eigen::Infinity>());
        maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(j,k,nJ,nK).lpNorm<Eigen::Infinity>());
        maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(j,l,nJ,nL).lpNorm<Eigen::Infinity>());
        maxDBlock = std::max(maxDBlock, densityMatrix_spin.block(k,l,nK,nL).lpNorm<Eigen::Infinity>());
      };
      if (maxDBlock*schwartz < _prescreeningThreshold) return true;
      /*
       * Shell quadruple is significant.
       */
      return false;
    };
    /*
     * Construct looper for loop over all integrals
     */
    TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb,0,_basis, _prescreeningThreshold);
    /*
     * Run
     */
    looper.loop(distribute, prescreeningFunc);
#ifdef _OPENMP
    for (unsigned int i=1; i<nThreads; ++i) {
      _potential->alpha += fc[i]->alpha;
      _potential->beta += fc[i]->beta;
      _excPotential->alpha += fx[i]->alpha;
      _excPotential->beta += fx[i]->beta;
      delete fx[i] ;
      delete fc[i] ;
    }
    _potential->alpha += _excPotential->alpha;
    _potential->beta += _excPotential->beta;
#endif
}


template<Options::SCF_MODES T> std::vector<unsigned int>
HFPotential<T>::createBasisToAtomIndexMapping(
      const std::vector<std::pair<unsigned int, unsigned int> >& basisIndicesRed,
      unsigned int nBasisFunctionsRed) {
  std::vector<unsigned int> mapping(nBasisFunctionsRed);
  // Vector to check whether ALL basis function shells are assigned to an atom index
  std::vector<bool> hasElementBeenSet(nBasisFunctionsRed, false);
  for (unsigned int iAtom=0; iAtom < basisIndicesRed.size(); ++iAtom) {
    const unsigned int firstIndex = basisIndicesRed[iAtom].first;
    const unsigned int endIndex = basisIndicesRed[iAtom].second;
    for (unsigned int iShell=firstIndex; iShell < endIndex; ++iShell) {
      mapping[iShell] = iAtom;
      hasElementBeenSet[iShell] = true;
    }
  }
  // Check
  for (bool x : hasElementBeenSet){
    if (not x)
      throw SerenityError("HFPotential: Missed gradient element in gradient evaluation.");
  } 
  return mapping;
}


template <> Eigen::MatrixXd
HFPotential<RESTRICTED>::getGeomGradients(){
  auto atoms = _systemController->getAtoms();
  const auto& orbitalSet = _systemController->getActiveOrbitalController<RESTRICTED>();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _systemController->getElectronicStructure<RESTRICTED>()->getDensityMatrix();

  unsigned int nBasisFunctionsRed = orbitalSet->getBasisController()->getReducedNBasisFunctions();
  auto basisIndicesRed = _systemController->getAtomCenteredBasisController()->getBasisIndicesRed();
  auto mapping = createBasisToAtomIndexMapping(basisIndicesRed,nBasisFunctionsRed);

  auto basis = _systemController->getBasisController();
  Libint& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::coulomb,1,4);

#ifdef _OPENMP
  // create a vector of matrices for each thread
  std::vector<Eigen::MatrixXd> eriContrPriv(omp_get_max_threads(),Eigen::MatrixXd::Zero(nAtoms,3));
  std::vector<Eigen::MatrixXd> exchangeContrPriv(omp_get_max_threads(),Eigen::MatrixXd::Zero(nAtoms,3));
#else
  // or just one
  std::vector<Eigen::MatrixXd > eriContrPriv(1,Eigen::MatrixXd::Zero(nAtoms,3));
  std::vector<Eigen::MatrixXd > exchangeContrPriv(1,Eigen::MatrixXd::Zero(nAtoms,3));
#endif

  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 1, basis, 1E-10);

  auto const looperFunction = [&]
                               (const unsigned int&  i,
                                   const unsigned int&  j,
                                   const unsigned int&  a,
                                   const unsigned int&  b,
                                   const Eigen::VectorXd& intValues,
                                   const unsigned int& threadID) {

    double perm = (i == j) ? 1.0 : 2.0;
    perm *= (a == b) ? 1.0 : 2.0;
    perm *= (i == a) ? (j == b ? 1.0 : 2.0) : 2.0;
    perm *= 0.5;

    for (unsigned int iAtom = 0; iAtom < 4; ++iAtom) {



      unsigned int nAtom;
      switch(iAtom){
        case(0): nAtom = mapping[basis->reducedIndex(i)]; break;
        case(1): nAtom = mapping[basis->reducedIndex(j)]; break;
        case(2): nAtom = mapping[basis->reducedIndex(a)]; break;
        case(3): nAtom = mapping[basis->reducedIndex(b)]; break;
      }

      for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
          double contr = perm * densityMatrix(i,j) * densityMatrix(a,b) * intValues(iAtom*3 + iDirection);
          eriContrPriv[threadID](nAtom,iDirection) += contr;

          //TODO: THE FOLLOWING ONLY WORKS FOR RHF SO FAR!!!
          double prefac = 0.0;
          prefac += 0.5 * densityMatrix(i,a) * densityMatrix(j,b);
          prefac += 0.5 * densityMatrix(i,b) * densityMatrix(j,a);

          double exchContr = 0.5 * perm *prefac*intValues(iAtom*3 + iDirection);

          exchangeContrPriv[threadID](nAtom,iDirection) += exchContr;
      }
    }
  };

  looper.loop(looperFunction);

  libint.finalize(libint2::Operator::coulomb,1,4);

  Matrix<double> eriContr(nAtoms,3);
  eriContr.setZero();
#ifdef _OPENMP
  // sum over all threads
  for (unsigned int i=0; i<(unsigned int)omp_get_max_threads(); ++i) {
    eriContr += eriContrPriv[i];
    eriContr -= exchangeContrPriv[i]*_exc;
  }
#else
  eriContr += eriContrPriv[0];
  eriContr -= exchangeContrPriv[0]*_exc;
#endif
  return eriContr;
}

template <> Eigen::MatrixXd
HFPotential<UNRESTRICTED>::getGeomGradients(){
  auto atoms = _systemController->getAtoms();
  const auto& orbitalSet = _systemController->getActiveOrbitalController<UNRESTRICTED>();
  unsigned int nAtoms = atoms.size();

  auto densityMatrix = _systemController->getElectronicStructure<UNRESTRICTED>()->getDensityMatrix();

  unsigned int nBasisFunctionsRed = orbitalSet->getBasisController()->getReducedNBasisFunctions();
  auto basisIndicesRed = _systemController->getAtomCenteredBasisController()->getBasisIndicesRed();
  auto mapping = createBasisToAtomIndexMapping(basisIndicesRed,nBasisFunctionsRed);

  auto basis = _systemController->getBasisController();
  Libint& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::coulomb,1,4);

#ifdef _OPENMP
  // create a vector of matrices for each thread
  std::vector<Eigen::MatrixXd> eriContrPriv(omp_get_max_threads(),Eigen::MatrixXd::Zero(nAtoms,3));
  std::vector<Eigen::MatrixXd> exchangeContrPriv(omp_get_max_threads(),Eigen::MatrixXd::Zero(nAtoms,3));
#else
  // or just one
  std::vector<Eigen::MatrixXd > eriContrPriv(1,Eigen::MatrixXd::Zero(nAtoms,3));
  std::vector<Eigen::MatrixXd > exchangeContrPriv(1,Eigen::MatrixXd::Zero(nAtoms,3));
#endif

  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 1, basis, 1E-10);

  auto const looperFunction = [&]
                               (const unsigned int&  i,
                                   const unsigned int&  j,
                                   const unsigned int&  a,
                                   const unsigned int&  b,
                                   const Eigen::VectorXd& intValues,
                                   const unsigned int& threadID) {

    double perm = (i == j) ? 1.0 : 2.0;
    perm *= (a == b) ? 1.0 : 2.0;
    perm *= (i == a) ? (j == b ? 1.0 : 2.0) : 2.0;
    perm *= 0.5;

    for (unsigned int iAtom = 0; iAtom < 4; ++iAtom) {



      unsigned int nAtom;
      switch(iAtom){
        case(0): nAtom = mapping[basis->reducedIndex(i)]; break;
        case(1): nAtom = mapping[basis->reducedIndex(j)]; break;
        case(2): nAtom = mapping[basis->reducedIndex(a)]; break;
        case(3): nAtom = mapping[basis->reducedIndex(b)]; break;
      }

      for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
          double contr = densityMatrix.alpha(i,j) * densityMatrix.alpha(a,b);
          contr += densityMatrix.beta(i,j)  * densityMatrix.beta(a,b);
          contr += densityMatrix.alpha(i,j) * densityMatrix.beta(a,b);
          contr += densityMatrix.beta(i,j) * densityMatrix.alpha(a,b);
          contr *= perm * intValues(iAtom*3 + iDirection);
          eriContrPriv[threadID](nAtom,iDirection) += contr;

          //TODO: THE FOLLOWING ONLY WORKS FOR RHF SO FAR!!!
          double prefac = 0.0;
          prefac += densityMatrix.alpha(i,a) * densityMatrix.alpha(j,b);
          prefac += densityMatrix.alpha(i,b) * densityMatrix.alpha(j,a);
          prefac += densityMatrix.beta(i,a) * densityMatrix.beta(j,b);
          prefac += densityMatrix.beta(i,b) * densityMatrix.beta(j,a);

          double exchContr = 0.5 * perm *prefac*intValues(iAtom*3 + iDirection);

          exchangeContrPriv[threadID](nAtom,iDirection) += exchContr;
      }
    }
  };

  looper.loop(looperFunction);

  libint.finalize(libint2::Operator::coulomb,1,4);

  Matrix<double> eriContr(nAtoms,3);
  eriContr.setZero();
#ifdef _OPENMP
  // sum over all threads
  for (unsigned int i=0; i<(unsigned int)omp_get_max_threads(); ++i) {
    eriContr += eriContrPriv[i];
    eriContr -= exchangeContrPriv[i]*_exc;
  }
#else
  eriContr += eriContrPriv[0];
  eriContr -= exchangeContrPriv[0]*_exc;
#endif
  return eriContr;
}

template <Options::SCF_MODES SCFMode> FockMatrix<SCFMode>&
HFPotential<SCFMode>::getXPotential() {
  if (!_potential or !_excPotential) getMatrix();
  return *_excPotential;
}

template class HFPotential<Options::SCF_MODES::RESTRICTED>;
template class HFPotential<Options::SCF_MODES::UNRESTRICTED>;


} /* namespace Serenity */
