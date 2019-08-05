/**
 * @file GSigmaVector.cpp
 *
 * @date Oct 09, 2017
 * @author Michael Boeckers
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
#include "postHF/LRSCF/SigmaVector/GSigmaVector.h"
/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "input/FunctionalClassResolver.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"

namespace Serenity {

template<Options::SCF_MODES T> GSigmaVector<T>::GSigmaVector(
    std::shared_ptr<SystemController> system,
    Eigen::MatrixXd& guessVector):
  SigmaVector<T>(system,guessVector){
  //Set HF exchange ratio
  _hfExchangeRatio = -1.0;
  _lrhfExchangeRatio = 0.0;
  _mu = 0.0;
  //Get settings for TDDFT
  if (system->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    auto functional = FunctionalClassResolver::resolveFunctional(
        system->getSettings().dft.functional);
    _hfExchangeRatio = -1.0 * functional.getHfExchangeRatio();
    _lrhfExchangeRatio = -1.0 * functional.getLRExchangeRatio();
    _mu = functional.getRangeSeparationParameter();
  }
}

template<>
Eigen::MatrixXd GSigmaVector<Options::SCF_MODES::RESTRICTED>::calculateBlock(
    Eigen::MatrixXd& guess,
    std::vector<SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::MatrixXd> >& dens) {
  const unsigned int nBFs = this->_system->getBasisController()->getNBasisFunctions();

  //Fock like matrices for each test vector in block
  std::vector<SPMatrix<Options::SCF_MODES::RESTRICTED>  > fock(guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    auto& f = fock[iGuess];
    for_spin(f) {
     f_spin.resize(nBFs,nBFs);
     f_spin.setZero();
    };
  }

  //Thread safety
  std::vector<std::vector<SPMatrix<Options::SCF_MODES::RESTRICTED>  > *> f;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  f.push_back(&fock);
  for (unsigned int iThread=1; iThread < nThreads; ++iThread) {
    f.push_back(new std::vector<SPMatrix<Options::SCF_MODES::RESTRICTED> >);
    f[iThread]->resize(guess.cols());
    for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
      auto& tf = (*f[iThread])[iGuess];
      for_spin(tf) {
        tf_spin.resize(nBFs,nBFs);
        tf_spin.setZero();
      };
    }
  }
#else
  f.push_back(&fock);
#endif

  // Function to add only Coulomb integrals to pseudo-Fock matrix (needed for pure TDDFT)
  //
  // F_{ij} = c_{coul}\sum_{kl} (ij|kl) P_{kl}
  //
  // Here, each integral is multiplied with a density matrix element. For LRSCF problems,
  // the (pseudo) density matrix is non symmetric. Every integral (ij|kl) is needed for two
  // Fock-matrix elements, F_{ij} and F_{kl}. Since the resulting Coulomb Fock-like matrix
  // is symmetric, we can add P_{kl} to P_{lk} (and P_{ij} to P_{ji} to save some time.
  // For the diagonal elements of P and F, we than need to divide by a factor of two which is
  // done by the the factor perm.
  auto distributeCoulomb = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    const unsigned int kl = k*nBFs+l;
    const unsigned int lk = l*nBFs+k;
    const unsigned int ij = i*nBFs+j;
    const unsigned int ji = j*nBFs+i;
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    double coul = 4.0 * perm * integral[0];
    for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
      auto& pd = dens[iGuess];
      auto& pf = (*f[threadId])[iGuess];
      double coul1 = 0;
      double coul2 = 0;
      coul1 += *(pd.data()+kl);
      coul1 += *(pd.data()+lk);
      coul2 += *(pd.data()+ij);
      coul2 += *(pd.data()+ji);
      coul1 *= coul;
      coul2 *= coul;
      *(pf.data()+ij) += coul1;
      *(pf.data()+ji) += coul1;
      *(pf.data()+kl) += coul2;
      *(pf.data()+lk) += coul2;
    }
  };

  //Function to add coulomb and exchange integrals to pseudo-Fock matrix. For the exchange, we have
  //
  // F_{ij} =  \sum_{kl} ( c_{ilkj)(il|kj) + c_{ikjl}(ik|jl) )P_{kl} .
  //
  // Here, each integral is multiplied with a density matrix element. For LRSCF problems,
  // the (pseudo) density matrix is non symmetric. Every integral (il|kj) amd (ik|jl) is
  // needed for eight Fock-matrix elements, which can easily be seen from the symmetry of
  // the two-electron integrals (i.e. (il|kj)=(li|kj)=(li|jk)=(il|jk)=(kj|il)=(kj|li)
  // =(jk|li)=(jk|il)) which than need to be multiplied with the correct density matrix
  // element. Given an orbital(ij|kl) from integral[] and the symmetry relation for the
  // integrals you can directly identify the density matrix element which is needed for
  // the Fock- matrix elements.
  auto distributeCoulombAndExchange = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    const unsigned int kl = k*nBFs+l;
    const unsigned int lk = l*nBFs+k;
    const unsigned int ij = i*nBFs+j;
    const unsigned int ji = j*nBFs+i;
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    double coul = perm * integral[0];
    const double exc = coul * _hfExchangeRatio;
    coul *= 4;
    for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
      auto& pd = dens[iGuess];
      auto& pf = (*f[threadId])[iGuess];

      //Coulomb
      double coul1 = 0;
      double coul2 = 0;
      coul1 += *(pd.data()+kl);
      coul1 += *(pd.data()+lk);
      coul2 += *(pd.data()+ij);
      coul2 += *(pd.data()+ji);
      coul1 *= coul;
      coul2 *= coul;
      *(pf.data()+ij) += coul1;
      *(pf.data()+ji) += coul1;
      *(pf.data()+kl) += coul2;
      *(pf.data()+lk) += coul2;

      //Exchange
      //ikjl
      *(pf.data() + i * nBFs + k) += *(pd.data() + l * nBFs + j) *  exc;
      *(pf.data() + k * nBFs + i) += *(pd.data() + j * nBFs + l) *  exc;
      *(pf.data() + l * nBFs + i) += *(pd.data() + j * nBFs + k) *  exc;
      *(pf.data() + i * nBFs + l) += *(pd.data() + k * nBFs + j) *  exc;
      *(pf.data() + k * nBFs + j) += *(pd.data() + i * nBFs + l) *  exc;
      *(pf.data() + j * nBFs + k) += *(pd.data() + l * nBFs + i) *  exc;
      *(pf.data() + l * nBFs + j) += *(pd.data() + i * nBFs + k) *  exc;
      *(pf.data() + j * nBFs + l) += *(pd.data() + k * nBFs + i) *  exc;

      //ilkj
      *(pf.data() + i * nBFs + k) += *(pd.data() + j * nBFs + l) *  exc;
      *(pf.data() + k * nBFs + i) += *(pd.data() + l * nBFs + j) *  exc;
      *(pf.data() + l * nBFs + i) += *(pd.data() + k * nBFs + j) *  exc;
      *(pf.data() + i * nBFs + l) += *(pd.data() + j * nBFs + k) *  exc;
      *(pf.data() + k * nBFs + j) += *(pd.data() + l * nBFs + i) *  exc;
      *(pf.data() + j * nBFs + k) += *(pd.data() + i * nBFs + l) *  exc;
      *(pf.data() + l * nBFs + j) += *(pd.data() + k * nBFs + i) *  exc;
      *(pf.data() + j * nBFs + l) += *(pd.data() + i * nBFs + k) *  exc;
    }
  };

  //Function to distribute the integrals for long-range exchange (for comments see above)
  auto distributeLRExchange = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    const double exc = perm * integral[0] * _lrhfExchangeRatio;
    for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
        auto& pd = dens[iGuess];
        auto& pf = (*f[threadId])[iGuess];
        //ikjl
        *(pf.data() + i * nBFs + k) += *(pd.data() + l * nBFs + j) *  exc;
        *(pf.data() + k * nBFs + i) += *(pd.data() + j * nBFs + l) *  exc;
        *(pf.data() + l * nBFs + i) += *(pd.data() + j * nBFs + k) *  exc;
        *(pf.data() + i * nBFs + l) += *(pd.data() + k * nBFs + j) *  exc;
        *(pf.data() + k * nBFs + j) += *(pd.data() + i * nBFs + l) *  exc;
        *(pf.data() + j * nBFs + k) += *(pd.data() + l * nBFs + i) *  exc;
        *(pf.data() + l * nBFs + j) += *(pd.data() + i * nBFs + k) *  exc;
        *(pf.data() + j * nBFs + l) += *(pd.data() + k * nBFs + i) *  exc;

        //ilkj
        *(pf.data() + i * nBFs + k) += *(pd.data() + j * nBFs + l) *  exc;
        *(pf.data() + k * nBFs + i) += *(pd.data() + l * nBFs + j) *  exc;
        *(pf.data() + l * nBFs + i) += *(pd.data() + k * nBFs + j) *  exc;
        *(pf.data() + i * nBFs + l) += *(pd.data() + j * nBFs + k) *  exc;
        *(pf.data() + k * nBFs + j) += *(pd.data() + l * nBFs + i) *  exc;
        *(pf.data() + j * nBFs + k) += *(pd.data() + i * nBFs + l) *  exc;
        *(pf.data() + l * nBFs + j) += *(pd.data() + k * nBFs + i) *  exc;
        *(pf.data() + j * nBFs + l) += *(pd.data() + i * nBFs + k) *  exc;
    }
  };

  //Maximum absolute value of all density matrices
  //ToDo: make argument?
  double maxDens = 0.0;
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    auto& pd = dens[iGuess];
    for_spin(pd) {
      maxDens = std::max(maxDens,fabs(pd_spin.maxCoeff()));
    };
  }

  //Detailed prescreening function
  auto prescreeningFunc = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      unsigned int nI, unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
    //    // No warnings
    (void)i;(void)j;(void)k;(void)l;
    (void)nI;(void)nJ;(void)nK;(void)nL;
    //Early return for insignificance based on the largest element in the whole density matrix
    if (maxDens*schwartz <  this->_system->getSettings().basis.integralThreshold) return true;
    //ToDo: Implement detailed prescreening procedure
    //Shell quadruple is significant.
    return false;
  };


  //Calculate pseudo Fock matrices
  TwoElecFourCenterIntLooper looper(
      libint2::Operator::coulomb,0,this->_system->getBasisController(), this->_system->getSettings().basis.integralThreshold);

  if (_hfExchangeRatio != 0.0) {
    looper.loop(distributeCoulombAndExchange,prescreeningFunc);
  } else {
    looper.loop(distributeCoulomb,prescreeningFunc);
  }

  //Add LR-Exchange
  if (_lrhfExchangeRatio != 0.0) {
    TwoElecFourCenterIntLooper lrlooper(
        libint2::Operator::erf_coulomb,0,this->_system->getBasisController(), this->_system->getSettings().basis.integralThreshold,_mu);
    lrlooper.loop(distributeLRExchange,prescreeningFunc);
  }

//Sum over all threads
#ifdef _OPENMP
  for (unsigned int iThread=1; iThread<nThreads; ++iThread) {
    for (unsigned int iGuess=0; iGuess<guess.cols(); ++iGuess) {
      (*f[0])[iGuess] += (*f[iThread])[iGuess];
    }
    delete f[iThread] ;
  }
#endif

  //Transform to MO basis and obtain sigma vectors
  Eigen::MatrixXd sigma(guess.rows(),guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    sigma.col(iGuess) = ao2mo(fock[iGuess]);
  }
  return sigma;
}

template<>
Eigen::MatrixXd GSigmaVector<Options::SCF_MODES::UNRESTRICTED>::calculateBlock(
    Eigen::MatrixXd& guess,
    std::vector<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> >& dens) {
  const unsigned int nBFs = this->_system->getBasisController()->getNBasisFunctions();

  //Fock like matrices for each test vector in block
  std::vector<SPMatrix<Options::SCF_MODES::UNRESTRICTED> > fock(guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    auto& f = fock[iGuess];
    for_spin(f) {
     f_spin.resize(nBFs,nBFs);
     f_spin.setZero();
    };
  }

  //Thread safety
  std::vector<std::vector<SPMatrix<Options::SCF_MODES::UNRESTRICTED> > *> f;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  f.push_back(&fock);
  for (unsigned int iThread=1; iThread < nThreads; ++iThread) {
    f.push_back(new std::vector<SPMatrix<Options::SCF_MODES::UNRESTRICTED> >);
    f[iThread]->resize(guess.cols());
    for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
      auto& tf = (*f[iThread])[iGuess];
      for_spin(tf) {
        tf_spin.resize(nBFs,nBFs);
        tf_spin.setZero();
      };
    }
  }
#else
  f.push_back(&fock);
#endif

  // Function to add only Coulomb integrals to pseudo-Fock matrix (needed for pure TDDFT)
  //
  // F_{ij} = c_{coul}\sum_{kl} (ij|kl) P_{kl}
  //
  // Here, each integral is multiplied with a density matrix element. For LRSCF problems,
  // the (pseudo) density matrix is non symmetric. Every integral (ij|kl) is needed for two
  // Fock-matrix elements, F_{ij} and F_{kl}. Since the resulting Coulomb Fock-like matrix
  // is symmetric, we can add P_{kl} to P_{lk} (and P_{ij} to P_{ji} to save some time.
  // For the diagonal elements of P and F, we than need to divide by a factor of two which is
  // done by the the factor perm.
  auto distributeCoulomb = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    const unsigned int kl = k*nBFs+l;
    const unsigned int lk = l*nBFs+k;
    const unsigned int ij = i*nBFs+j;
    const unsigned int ji = j*nBFs+i;
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    double coul = 2.0 * perm * integral[0];
    for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
        auto& pd = dens[iGuess];
        auto& pf = (*f[threadId])[iGuess];
        double coul1 = 0;
        double coul2 = 0;
        coul1 += *(pd.alpha.data()+kl);
        coul1 += *(pd.alpha.data()+lk);
        coul1 += *(pd.beta.data()+kl);
        coul1 += *(pd.beta.data()+lk);
        coul2 += *(pd.alpha.data()+ij);
        coul2 += *(pd.alpha.data()+ji);
        coul2 += *(pd.beta.data()+ij);
        coul2 += *(pd.beta.data()+ji);

        coul1 *= coul;
        coul2 *= coul;
        *(pf.alpha.data()+ij) += coul1;
        *(pf.beta.data()+ij) += coul1;
        *(pf.alpha.data()+ji) += coul1;
        *(pf.beta.data()+ji) += coul1;
        *(pf.alpha.data()+kl) += coul2;
        *(pf.beta.data()+kl) += coul2;
        *(pf.alpha.data()+lk) += coul2;
        *(pf.beta.data()+lk) += coul2;
    }
  };

  //Function to add coulomb and exchange integrals to pseudo-Fock matrix.
  auto distributeCoulombAndExchange = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    const unsigned int kl = k*nBFs+l;
    const unsigned int lk = l*nBFs+k;
    const unsigned int ij = i*nBFs+j;
    const unsigned int ji = j*nBFs+i;
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    double coul = perm * integral[0];
    const double exc = coul * _hfExchangeRatio;
    coul *= 2.0;
    for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
      auto& pd = dens[iGuess];
      auto& pf = (*f[threadId])[iGuess];

      //Coulomb
      double coul1 = 0;
      double coul2 = 0;
      coul1 += *(pd.alpha.data()+kl);
      coul1 += *(pd.alpha.data()+lk);
      coul1 += *(pd.beta.data()+kl);
      coul1 += *(pd.beta.data()+lk);
      coul2 += *(pd.alpha.data()+ij);
      coul2 += *(pd.alpha.data()+ji);
      coul2 += *(pd.beta.data()+ij);
      coul2 += *(pd.beta.data()+ji);

      coul1 *= coul;
      coul2 *= coul;
      *(pf.alpha.data()+ij) += coul1;
      *(pf.beta.data()+ij) += coul1;
      *(pf.alpha.data()+ji) += coul1;
      *(pf.beta.data()+ji) += coul1;
      *(pf.alpha.data()+kl) += coul2;
      *(pf.beta.data()+kl) += coul2;
      *(pf.alpha.data()+lk) += coul2;
      *(pf.beta.data()+lk) += coul2;

      //Exchange
      //ikjl
      *(pf.alpha.data() +  i * nBFs +k) += *(pd.alpha.data() +  l * nBFs +j) *  exc;
      *(pf.alpha.data() +  k * nBFs +i) += *(pd.alpha.data() +  j * nBFs +l) *  exc;
      *(pf.alpha.data() +  l * nBFs +i) += *(pd.alpha.data() +  j * nBFs +k) *  exc;
      *(pf.alpha.data() +  i * nBFs +l) += *(pd.alpha.data() +  k * nBFs +j) *  exc;
      *(pf.alpha.data() +  k * nBFs +j) += *(pd.alpha.data() +  i * nBFs +l) *  exc;
      *(pf.alpha.data() +  j * nBFs +k) += *(pd.alpha.data() +  l * nBFs +i) *  exc;
      *(pf.alpha.data() +  l * nBFs +j) += *(pd.alpha.data() +  i * nBFs +k) *  exc;
      *(pf.alpha.data() +  j * nBFs +l) += *(pd.alpha.data() +  k * nBFs +i) *  exc;

      *(pf.beta.data() +  i * nBFs +k) += *(pd.beta.data() +  l * nBFs +j) *  exc;
      *(pf.beta.data() +  k * nBFs +i) += *(pd.beta.data() +  j * nBFs +l) *  exc;
      *(pf.beta.data() +  l * nBFs +i) += *(pd.beta.data() +  j * nBFs +k) *  exc;
      *(pf.beta.data() +  i * nBFs +l) += *(pd.beta.data() +  k * nBFs +j) *  exc;
      *(pf.beta.data() +  k * nBFs +j) += *(pd.beta.data() +  i * nBFs +l) *  exc;
      *(pf.beta.data() +  j * nBFs +k) += *(pd.beta.data() +  l * nBFs +i) *  exc;
      *(pf.beta.data() +  l * nBFs +j) += *(pd.beta.data() +  i * nBFs +k) *  exc;
      *(pf.beta.data() +  j * nBFs +l) += *(pd.beta.data() +  k * nBFs +i) *  exc;


      //ilkj
      *(pf.alpha.data() +  i * nBFs +k) += *(pd.alpha.data() +  j * nBFs +l) *  exc;
      *(pf.alpha.data() +  k * nBFs +i) += *(pd.alpha.data() +  l * nBFs +j) *  exc;
      *(pf.alpha.data() +  l * nBFs +i) += *(pd.alpha.data() +  k * nBFs +j) *  exc;
      *(pf.alpha.data() +  i * nBFs +l) += *(pd.alpha.data() +  j * nBFs +k) *  exc;
      *(pf.alpha.data() +  k * nBFs +j) += *(pd.alpha.data() +  l * nBFs +i) *  exc;
      *(pf.alpha.data() +  j * nBFs +k) += *(pd.alpha.data() +  i * nBFs +l) *  exc;
      *(pf.alpha.data() +  l * nBFs +j) += *(pd.alpha.data() +  k * nBFs +i) *  exc;
      *(pf.alpha.data() +  j * nBFs +l) += *(pd.alpha.data() +  i * nBFs +k) *  exc;

      *(pf.beta.data() +  i * nBFs +k) += *(pd.beta.data() +  j * nBFs +l) *  exc;
      *(pf.beta.data() +  k * nBFs +i) += *(pd.beta.data() +  l * nBFs +j) *  exc;
      *(pf.beta.data() +  l * nBFs +i) += *(pd.beta.data() +  k * nBFs +j) *  exc;
      *(pf.beta.data() +  i * nBFs +l) += *(pd.beta.data() +  j * nBFs +k) *  exc;
      *(pf.beta.data() +  k * nBFs +j) += *(pd.beta.data() +  l * nBFs +i) *  exc;
      *(pf.beta.data() +  j * nBFs +k) += *(pd.beta.data() +  i * nBFs +l) *  exc;
      *(pf.beta.data() +  l * nBFs +j) += *(pd.beta.data() +  k * nBFs +i) *  exc;
      *(pf.beta.data() +  j * nBFs +l) += *(pd.beta.data() +  i * nBFs +k) *  exc;
    }
  };

  //Function to distribute the integrals for long-range exchange (for comments see above)
  auto distributeLRExchange = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    const double exc = perm * integral[0] * _lrhfExchangeRatio;
    for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
        auto& pd = dens[iGuess];
        auto& pf = (*f[threadId])[iGuess];
        //ikjl
        *(pf.alpha.data() +  i * nBFs +k) += *(pd.alpha.data() +  l * nBFs +j) *  exc;
        *(pf.alpha.data() +  k * nBFs +i) += *(pd.alpha.data() +  j * nBFs +l) *  exc;
        *(pf.alpha.data() +  l * nBFs +i) += *(pd.alpha.data() +  j * nBFs +k) *  exc;
        *(pf.alpha.data() +  i * nBFs +l) += *(pd.alpha.data() +  k * nBFs +j) *  exc;
        *(pf.alpha.data() +  k * nBFs +j) += *(pd.alpha.data() +  i * nBFs +l) *  exc;
        *(pf.alpha.data() +  j * nBFs +k) += *(pd.alpha.data() +  l * nBFs +i) *  exc;
        *(pf.alpha.data() +  l * nBFs +j) += *(pd.alpha.data() +  i * nBFs +k) *  exc;
        *(pf.alpha.data() +  j * nBFs +l) += *(pd.alpha.data() +  k * nBFs +i) *  exc;

        *(pf.beta.data() +  i * nBFs +k) += *(pd.beta.data() +  l * nBFs +j) *  exc;
        *(pf.beta.data() +  k * nBFs +i) += *(pd.beta.data() +  j * nBFs +l) *  exc;
        *(pf.beta.data() +  l * nBFs +i) += *(pd.beta.data() +  j * nBFs +k) *  exc;
        *(pf.beta.data() +  i * nBFs +l) += *(pd.beta.data() +  k * nBFs +j) *  exc;
        *(pf.beta.data() +  k * nBFs +j) += *(pd.beta.data() +  i * nBFs +l) *  exc;
        *(pf.beta.data() +  j * nBFs +k) += *(pd.beta.data() +  l * nBFs +i) *  exc;
        *(pf.beta.data() +  l * nBFs +j) += *(pd.beta.data() +  i * nBFs +k) *  exc;
        *(pf.beta.data() +  j * nBFs +l) += *(pd.beta.data() +  k * nBFs +i) *  exc;


        //ilkj
        *(pf.alpha.data() +  i * nBFs +k) += *(pd.alpha.data() +  j * nBFs +l) *  exc;
        *(pf.alpha.data() +  k * nBFs +i) += *(pd.alpha.data() +  l * nBFs +j) *  exc;
        *(pf.alpha.data() +  l * nBFs +i) += *(pd.alpha.data() +  k * nBFs +j) *  exc;
        *(pf.alpha.data() +  i * nBFs +l) += *(pd.alpha.data() +  j * nBFs +k) *  exc;
        *(pf.alpha.data() +  k * nBFs +j) += *(pd.alpha.data() +  l * nBFs +i) *  exc;
        *(pf.alpha.data() +  j * nBFs +k) += *(pd.alpha.data() +  i * nBFs +l) *  exc;
        *(pf.alpha.data() +  l * nBFs +j) += *(pd.alpha.data() +  k * nBFs +i) *  exc;
        *(pf.alpha.data() +  j * nBFs +l) += *(pd.alpha.data() +  i * nBFs +k) *  exc;

        *(pf.beta.data() +  i * nBFs +k) += *(pd.beta.data() +  j * nBFs +l) *  exc;
        *(pf.beta.data() +  k * nBFs +i) += *(pd.beta.data() +  l * nBFs +j) *  exc;
        *(pf.beta.data() +  l * nBFs +i) += *(pd.beta.data() +  k * nBFs +j) *  exc;
        *(pf.beta.data() +  i * nBFs +l) += *(pd.beta.data() +  j * nBFs +k) *  exc;
        *(pf.beta.data() +  k * nBFs +j) += *(pd.beta.data() +  l * nBFs +i) *  exc;
        *(pf.beta.data() +  j * nBFs +k) += *(pd.beta.data() +  i * nBFs +l) *  exc;
        *(pf.beta.data() +  l * nBFs +j) += *(pd.beta.data() +  k * nBFs +i) *  exc;
        *(pf.beta.data() +  j * nBFs +l) += *(pd.beta.data() +  i * nBFs +k) *  exc;
    }
  };

  //Maximum absolute value of all density matrices
  //ToDo: make argument?
  double maxDens = 0.0;
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    auto& pd = dens[iGuess];
    for_spin(pd) {
      maxDens = std::max(maxDens,fabs(pd_spin.maxCoeff()));
    };
  }

  //Detailed prescreening function
  auto prescreeningFunc = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      unsigned int nI, unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
    //    // No warnings
    (void)i;(void)j;(void)k;(void)l;
    (void)nI;(void)nJ;(void)nK;(void)nL;
    //Early return for insignificance based on the largest element in the whole density matrix
    if (maxDens*schwartz <  this->_system->getSettings().basis.integralThreshold) return true;
    //ToDo: Implement detailed prescreening procedure
    //Shell quadruple is significant.
    return false;
  };


  //Calculate pseudo Fock matrices
  TwoElecFourCenterIntLooper looper(
      libint2::Operator::coulomb,0,this->_system->getBasisController(), this->_system->getSettings().basis.integralThreshold);

  if (_hfExchangeRatio != 0.0) {
    looper.loop(distributeCoulombAndExchange,prescreeningFunc);
  } else {
    looper.loop(distributeCoulomb,prescreeningFunc);
  }

  //Add LR-Exchange
  if (_lrhfExchangeRatio != 0.0) {
    TwoElecFourCenterIntLooper lrlooper(
        libint2::Operator::erf_coulomb,0,this->_system->getBasisController(), this->_system->getSettings().basis.integralThreshold,_mu);
    lrlooper.loop(distributeLRExchange,prescreeningFunc);
  }

//Sum over all threads
#ifdef _OPENMP
  for (unsigned int iThread=1; iThread<nThreads; ++iThread) {
    for (unsigned int iGuess=0; iGuess<guess.cols(); ++iGuess) {
      (*f[0])[iGuess] += (*f[iThread])[iGuess];
    }
    delete f[iThread] ;
  }
#endif

  //Transform to MO basis and obtain sigma vectors
  Eigen::MatrixXd sigma(guess.rows(),guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    sigma.col(iGuess) = ao2mo(fock[iGuess]);
  }
  return sigma;
}

template class GSigmaVector<Options::SCF_MODES::RESTRICTED> ;
template class GSigmaVector<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
