/**
 * @file KSigmaVector.cpp
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
#include "postHF/LRSCF/SigmaVector/KSigmaVector.h"
/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "input/FunctionalClassResolver.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"

namespace Serenity {

template<Options::SCF_MODES T> KSigmaVector<T>::KSigmaVector(
    std::shared_ptr<SystemController> system,
    Eigen::MatrixXd& guessVector,
    int pm):
  SigmaVector<T>(system,guessVector),
  _pm(pm){
  assert(fabs(_pm) == 1.0);
  //Set HF exchange ratio
  _hfExchangeRatio = 1.0;
  _lrhfExchangeRatio = 0.0;
  _mu = 0.0;
  //Get settings for TDDFT
  if (system->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    auto functional = FunctionalClassResolver::resolveFunctional(
        system->getSettings().dft.functional);
    _hfExchangeRatio = functional.getHfExchangeRatio();
    _lrhfExchangeRatio = functional.getLRExchangeRatio();
    _mu = functional.getRangeSeparationParameter();
  }
}

template<>
Eigen::MatrixXd KSigmaVector<Options::SCF_MODES::RESTRICTED>::calculateBlock(
    Eigen::MatrixXd& guess,
    std::vector<SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::MatrixXd> >& dens) {
  //For comments see GSigmaVector.cpp
  const unsigned int nBFs = this->_system->getBasisController()->getNBasisFunctions();

  std::vector<SPMatrix<Options::SCF_MODES::RESTRICTED> > fock(guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    auto& f = fock[iGuess];
    for_spin(f) {
     f_spin.resize(nBFs,nBFs);
     f_spin.setZero();
    };
  }

  //Thread safety
  std::vector<std::vector<SPMatrix<Options::SCF_MODES::RESTRICTED> > *> f;
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

  auto distributeExchange = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    const double exc = perm * integral[0] * _hfExchangeRatio;
    const double exc_ilkj = _pm * exc;
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
        *(pf.data() + i * nBFs + k) += *(pd.data() + j * nBFs + l) *  exc_ilkj;
        *(pf.data() + k * nBFs + i) += *(pd.data() + l * nBFs + j) *  exc_ilkj;
        *(pf.data() + l * nBFs + i) += *(pd.data() + k * nBFs + j) *  exc_ilkj;
        *(pf.data() + i * nBFs + l) += *(pd.data() + j * nBFs + k) *  exc_ilkj;
        *(pf.data() + k * nBFs + j) += *(pd.data() + l * nBFs + i) *  exc_ilkj;
        *(pf.data() + j * nBFs + k) += *(pd.data() + i * nBFs + l) *  exc_ilkj;
        *(pf.data() + l * nBFs + j) += *(pd.data() + k * nBFs + i) *  exc_ilkj;
        *(pf.data() + j * nBFs + l) += *(pd.data() + i * nBFs + k) *  exc_ilkj;
    }
  };

  auto distributeLRExchange = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    const double exc = perm * integral[0] * _lrhfExchangeRatio;
    const double exc_ilkj = _pm * exc;
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
        *(pf.data() + i * nBFs + k) += *(pd.data() + j * nBFs + l) *  exc_ilkj;
        *(pf.data() + k * nBFs + i) += *(pd.data() + l * nBFs + j) *  exc_ilkj;
        *(pf.data() + l * nBFs + i) += *(pd.data() + k * nBFs + j) *  exc_ilkj;
        *(pf.data() + i * nBFs + l) += *(pd.data() + j * nBFs + k) *  exc_ilkj;
        *(pf.data() + k * nBFs + j) += *(pd.data() + l * nBFs + i) *  exc_ilkj;
        *(pf.data() + j * nBFs + k) += *(pd.data() + i * nBFs + l) *  exc_ilkj;
        *(pf.data() + l * nBFs + j) += *(pd.data() + k * nBFs + i) *  exc_ilkj;
        *(pf.data() + j * nBFs + l) += *(pd.data() + i * nBFs + k) *  exc_ilkj;
    }
  };

  double maxDens = 0.0;
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    auto& pd = dens[iGuess];
    for_spin(pd) {
      maxDens = std::max(maxDens,fabs(pd_spin.maxCoeff()));
    };
  }

  auto prescreeningFunc = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      unsigned int nI, unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
    (void)i;(void)j;(void)k;(void)l;
    (void)nI;(void)nJ;(void)nK;(void)nL;
    if (maxDens*schwartz <  this->_system->getSettings().basis.integralThreshold) return true;
    return false;
  };

  if (_hfExchangeRatio != 0.0) {
    TwoElecFourCenterIntLooper looper(
        libint2::Operator::coulomb,0,this->_system->getBasisController(), this->_system->getSettings().basis.integralThreshold);
    looper.loop(distributeExchange,prescreeningFunc);
  }

  if (_lrhfExchangeRatio != 0.0) {
    TwoElecFourCenterIntLooper lrlooper(
        libint2::Operator::erf_coulomb,0,this->_system->getBasisController(), this->_system->getSettings().basis.integralThreshold,_mu);
    lrlooper.loop(distributeLRExchange,prescreeningFunc);
  }

#ifdef _OPENMP
  for (unsigned int iThread=1; iThread<nThreads; ++iThread) {
    for (unsigned int iGuess=0; iGuess<guess.cols(); ++iGuess) {
      (*f[0])[iGuess] += (*f[iThread])[iGuess];
    }
    delete f[iThread] ;
  }
#endif

  Eigen::MatrixXd sigma(guess.rows(),guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    sigma.col(iGuess) = ao2mo(fock[iGuess]);
  }
  return sigma;
}

template<>
Eigen::MatrixXd KSigmaVector<Options::SCF_MODES::UNRESTRICTED>::calculateBlock(
    Eigen::MatrixXd& guess,
    std::vector<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> >& dens) {
  //For comments see GSigmaVector.cpp
  const unsigned int nBFs = this->_system->getBasisController()->getNBasisFunctions();

  std::vector<SPMatrix<Options::SCF_MODES::UNRESTRICTED> >
  fock(guess.cols(),SPMatrix<Options::SCF_MODES::UNRESTRICTED>(nBFs,nBFs));

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

  auto distributeExchange = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    const double exc = perm * integral[0] * _hfExchangeRatio;
    const double exc_ilkj = _pm * exc;
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
        *(pf.alpha.data() +  i * nBFs +k) += *(pd.alpha.data() +  j * nBFs +l) *  exc_ilkj;
        *(pf.alpha.data() +  k * nBFs +i) += *(pd.alpha.data() +  l * nBFs +j) *  exc_ilkj;
        *(pf.alpha.data() +  l * nBFs +i) += *(pd.alpha.data() +  k * nBFs +j) *  exc_ilkj;
        *(pf.alpha.data() +  i * nBFs +l) += *(pd.alpha.data() +  j * nBFs +k) *  exc_ilkj;
        *(pf.alpha.data() +  k * nBFs +j) += *(pd.alpha.data() +  l * nBFs +i) *  exc_ilkj;
        *(pf.alpha.data() +  j * nBFs +k) += *(pd.alpha.data() +  i * nBFs +l) *  exc_ilkj;
        *(pf.alpha.data() +  l * nBFs +j) += *(pd.alpha.data() +  k * nBFs +i) *  exc_ilkj;
        *(pf.alpha.data() +  j * nBFs +l) += *(pd.alpha.data() +  i * nBFs +k) *  exc_ilkj;

        *(pf.beta.data() +  i * nBFs +k) += *(pd.beta.data() +  j * nBFs +l) *  exc_ilkj;
        *(pf.beta.data() +  k * nBFs +i) += *(pd.beta.data() +  l * nBFs +j) *  exc_ilkj;
        *(pf.beta.data() +  l * nBFs +i) += *(pd.beta.data() +  k * nBFs +j) *  exc_ilkj;
        *(pf.beta.data() +  i * nBFs +l) += *(pd.beta.data() +  j * nBFs +k) *  exc_ilkj;
        *(pf.beta.data() +  k * nBFs +j) += *(pd.beta.data() +  l * nBFs +i) *  exc_ilkj;
        *(pf.beta.data() +  j * nBFs +k) += *(pd.beta.data() +  i * nBFs +l) *  exc_ilkj;
        *(pf.beta.data() +  l * nBFs +j) += *(pd.beta.data() +  k * nBFs +i) *  exc_ilkj;
        *(pf.beta.data() +  j * nBFs +l) += *(pd.beta.data() +  i * nBFs +k) *  exc_ilkj;
    }
  };

  auto distributeLRExchange = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
    //Permutations
    double perm =1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;
    perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
    const double exc = perm * integral[0] * _lrhfExchangeRatio;
    const double exc_ilkj = _pm * exc;
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
        *(pf.alpha.data() +  i * nBFs +k) += *(pd.alpha.data() +  j * nBFs +l) *  exc_ilkj;
        *(pf.alpha.data() +  k * nBFs +i) += *(pd.alpha.data() +  l * nBFs +j) *  exc_ilkj;
        *(pf.alpha.data() +  l * nBFs +i) += *(pd.alpha.data() +  k * nBFs +j) *  exc_ilkj;
        *(pf.alpha.data() +  i * nBFs +l) += *(pd.alpha.data() +  j * nBFs +k) *  exc_ilkj;
        *(pf.alpha.data() +  k * nBFs +j) += *(pd.alpha.data() +  l * nBFs +i) *  exc_ilkj;
        *(pf.alpha.data() +  j * nBFs +k) += *(pd.alpha.data() +  i * nBFs +l) *  exc_ilkj;
        *(pf.alpha.data() +  l * nBFs +j) += *(pd.alpha.data() +  k * nBFs +i) *  exc_ilkj;
        *(pf.alpha.data() +  j * nBFs +l) += *(pd.alpha.data() +  i * nBFs +k) *  exc_ilkj;

        *(pf.beta.data() +  i * nBFs +k) += *(pd.beta.data() +  j * nBFs +l) *  exc_ilkj;
        *(pf.beta.data() +  k * nBFs +i) += *(pd.beta.data() +  l * nBFs +j) *  exc_ilkj;
        *(pf.beta.data() +  l * nBFs +i) += *(pd.beta.data() +  k * nBFs +j) *  exc_ilkj;
        *(pf.beta.data() +  i * nBFs +l) += *(pd.beta.data() +  j * nBFs +k) *  exc_ilkj;
        *(pf.beta.data() +  k * nBFs +j) += *(pd.beta.data() +  l * nBFs +i) *  exc_ilkj;
        *(pf.beta.data() +  j * nBFs +k) += *(pd.beta.data() +  i * nBFs +l) *  exc_ilkj;
        *(pf.beta.data() +  l * nBFs +j) += *(pd.beta.data() +  k * nBFs +i) *  exc_ilkj;
        *(pf.beta.data() +  j * nBFs +l) += *(pd.beta.data() +  i * nBFs +k) *  exc_ilkj;
    }
  };

  double maxDens = 0.0;
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    auto& pd = dens[iGuess];
    for_spin(pd) {
      maxDens = std::max(maxDens,fabs(pd_spin.maxCoeff()));
    };
  }

  auto prescreeningFunc = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      unsigned int nI, unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
    (void)i;(void)j;(void)k;(void)l;
    (void)nI;(void)nJ;(void)nK;(void)nL;
    if (maxDens*schwartz <  this->_system->getSettings().basis.integralThreshold) return true;
    return false;
  };


  if (_hfExchangeRatio != 0.0) {
    TwoElecFourCenterIntLooper looper(
        libint2::Operator::coulomb,0,this->_system->getBasisController(), this->_system->getSettings().basis.integralThreshold);
    looper.loop(distributeExchange,prescreeningFunc);
  }

  if (_lrhfExchangeRatio != 0.0) {
    TwoElecFourCenterIntLooper lrlooper(
        libint2::Operator::erf_coulomb,0,this->_system->getBasisController(), this->_system->getSettings().basis.integralThreshold,_mu);
    lrlooper.loop(distributeLRExchange,prescreeningFunc);
  }

#ifdef _OPENMP
  for (unsigned int iThread=1; iThread<nThreads; ++iThread) {
    for (unsigned int iGuess=0; iGuess<guess.cols(); ++iGuess) {
      (*f[0])[iGuess] += (*f[iThread])[iGuess];
    }
    delete f[iThread] ;
  }
#endif

  Eigen::MatrixXd sigma(guess.rows(),guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    sigma.col(iGuess) = ao2mo(fock[iGuess]);
  }
  return sigma;
}

template class KSigmaVector<Options::SCF_MODES::RESTRICTED> ;
template class KSigmaVector<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
