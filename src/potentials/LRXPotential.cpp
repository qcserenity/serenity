/**
 * @file LRXPotential.cpp
 *
 * @date Mar 31, 2017
 * @author M. Boeckers
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
#include "potentials/LRXPotential.h"
/* Include Serenity Internal Headers */
#include "data/matrices/FockMatrix.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
/* Include Std and External Headers */
#include <algorithm>


namespace Serenity {

template <Options::SCF_MODES SCFMode>
LRXPotential<SCFMode>::LRXPotential(
    std::shared_ptr<DensityMatrixController<SCFMode> > dMat,
    const double exchangeRatio,
    const double prescreeningThreshold,
    const double mu):
    Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _prescreeningThreshold(prescreeningThreshold),
    _exc(exchangeRatio),
    _dMatController(dMat),
    _potential(nullptr),
    _mu(mu){
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode> >::_self);
}

template <Options::SCF_MODES SCFMode> FockMatrix<SCFMode>&
LRXPotential<SCFMode>::getMatrix(){
  if (!_potential){
    const auto& densityMatrix = _dMatController->getDensityMatrix();
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot){
      pot_spin.setZero();
    };
    this->addToMatrix(*_potential,densityMatrix);
  }
  return *_potential;
}

template <Options::SCF_MODES SCFMode> double
LRXPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P){
  if (!_potential) this->getMatrix();
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot,P){
    energy += 0.5*pot_spin.cwiseProduct(P_spin).sum();
  };
  return energy;
};

template <> void
LRXPotential<Options::SCF_MODES::RESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::RESTRICTED>& F,
    const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densityMatrix){

    const unsigned int nBFs = _basis->getNBasisFunctions();

    /*
     * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
     */
    std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> f;
#ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    // Let the first thread use the Fock matrix directly
    f.push_back(&F);
    for (unsigned int i=1; i<nThreads; ++i) {
      f.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basis));
      f[i]->setZero();
    }

#else
    f.push_back(&F);
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
       * Exchange
       */
      const double exc = perm * integral[0] * 0.25 * _exc;
      const double exc1 = *(densityMatrix.data()+i*nBFs+k)  * exc;
      const double exc2 = *(densityMatrix.data()+i*nBFs+l)  * exc;
      const double exc3 = *(densityMatrix.data()+j*nBFs+k)  * exc;
      const double exc4 = *(densityMatrix.data()+j*nBFs+l)  * exc;
      *(f[threadId]->data()+j*nBFs+l) -= exc1;
      *(f[threadId]->data()+l*nBFs+j) -= exc1;
      *(f[threadId]->data()+j*nBFs+k) -= exc2;
      *(f[threadId]->data()+k*nBFs+j) -= exc2;
      *(f[threadId]->data()+i*nBFs+l) -= exc3;
      *(f[threadId]->data()+l*nBFs+i) -= exc3;
      *(f[threadId]->data()+i*nBFs+k) -= exc4;
      *(f[threadId]->data()+k*nBFs+i) -= exc4;

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
    TwoElecFourCenterIntLooper looper(libint2::Operator::erf_coulomb,0,_basis, _prescreeningThreshold,_mu);
    /*
     * Run
     */
    looper.loop(distribute, prescreeningFunc);
#ifdef _OPENMP
    for (unsigned int i=1; i<nThreads; ++i) {
      F += *f[i];
      delete f[i] ;
    }
#endif

}

template <> void
LRXPotential<Options::SCF_MODES::UNRESTRICTED>::addToMatrix(
    FockMatrix<Options::SCF_MODES::UNRESTRICTED>& F,
    const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densityMatrix){
    /*
     * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
     */
    std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>* > f;
#ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    // Let the first thread use the Fock matrix directly
    f.push_back(&F);
    for (unsigned int i=1; i<nThreads; ++i) {
      f.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basis));
      f[i]->alpha.setZero();
      f[i]->beta.setZero();
    }

#else
    // Simply use the fock matrix directly
    f.push_back(&F);
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
       * Exchange
       */

      const double exc = perm * integral[0] * 0.5 * _exc;

      const double exc1a = densityMatrix.alpha(i,k)  * exc;
      const double exc2a = densityMatrix.alpha(i,l)  * exc;
      const double exc3a = densityMatrix.alpha(j,k)  * exc;
      const double exc4a = densityMatrix.alpha(j,l)  * exc;
      f[threadId]->alpha(j,l) -= exc1a;
      f[threadId]->alpha(l,j) -= exc1a;
      f[threadId]->alpha(j,k) -= exc2a;
      f[threadId]->alpha(k,j) -= exc2a;
      f[threadId]->alpha(i,l) -= exc3a;
      f[threadId]->alpha(l,i) -= exc3a;
      f[threadId]->alpha(i,k) -= exc4a;
      f[threadId]->alpha(k,i) -= exc4a;
      const double exc1b = densityMatrix.beta(i,k)  * exc;
      const double exc2b = densityMatrix.beta(i,l)  * exc;
      const double exc3b = densityMatrix.beta(j,k)  * exc;
      const double exc4b = densityMatrix.beta(j,l)  * exc;
      f[threadId]->beta(j,l) -= exc1b;
      f[threadId]->beta(l,j) -= exc1b;
      f[threadId]->beta(j,k) -= exc2b;
      f[threadId]->beta(k,j) -= exc2b;
      f[threadId]->beta(i,l) -= exc3b;
      f[threadId]->beta(l,i) -= exc3b;
      f[threadId]->beta(i,k) -= exc4b;
      f[threadId]->beta(k,i) -= exc4b;

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
    TwoElecFourCenterIntLooper looper(libint2::Operator::erf_coulomb,0,_basis, _prescreeningThreshold,_mu);
    /*
     * Run
     */
    looper.loop(distribute, prescreeningFunc);
#ifdef _OPENMP
    for (unsigned int i=1; i<nThreads; ++i) {
      F.alpha += f[i]->alpha;
      F.beta += f[i]->beta;
      delete f[i] ;
    }
#endif
}


template <Options::SCF_MODES SCFMode> Eigen::MatrixXd
LRXPotential<SCFMode>::getGeomGradients(){
  Eigen::MatrixXd gradientContr(1,3);
  gradientContr.setZero();

  return gradientContr;
}

template class LRXPotential<Options::SCF_MODES::RESTRICTED>;
template class LRXPotential<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
