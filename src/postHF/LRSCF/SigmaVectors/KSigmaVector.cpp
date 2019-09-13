/**
 * @file KSigmaVector.cpp
 *
 * @date Dec 11, 2018
 * @author Michael Boeckers, Johannes Toelle
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
#include "postHF/LRSCF/SigmaVectors/KSigmaVector.h"

/* Include Serenity Internal Headers */
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/looper/ExchangeInteractionIntLooper.h"
#include "input/FunctionalClassResolver.h"
#include "dft/Functional.h"
#include "misc/Timing.h"


namespace Serenity {

template<Options::SCF_MODES SCFMode>
KSigmaVector<SCFMode>::KSigmaVector(
    std::vector<std::shared_ptr<LRSCFController<SCFMode> > > lrscf,
    std::vector<Eigen::MatrixXd> b,
    const double densityScreeningThreshold,
    const std::vector<int> pm):
    SigmaVector<SCFMode>(lrscf,b,densityScreeningThreshold),
    _pm(pm){}


template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > > > KSigmaVector<Options::SCF_MODES::RESTRICTED>::calcF(
    unsigned int I,
    unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > > > densityMatrices){


    Timings::takeTime("LRSCF -   Fock-like matrix: K");

    //Set exchange ratios
    double hfExchangeRatio = 1.0;
    double lrExchangeRatio = 0.0;
    double mu = 0.0;

    if (this->_lrscf[I]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT ||
        this->_lrscf[J]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      auto funcI = FunctionalClassResolver::resolveFunctional(
          this->_lrscf[I]->getSysSettings().dft.functional);
      auto funcJ = FunctionalClassResolver::resolveFunctional(
          this->_lrscf[J]->getSysSettings().dft.functional);
      //Can only use exact exchange, if the same amount of it is used in every subsystem. Otherwise, the
      //Response matrix becomes non-symmetric.
      if (funcI.getHfExchangeRatio() != funcJ.getHfExchangeRatio()) throw std::runtime_error("HF exchange ratio must be the same in all subsystems");
      if (funcI.getLRExchangeRatio() != funcJ.getLRExchangeRatio()) throw std::runtime_error("HF long range exchange ratio must be the same in all subsystems");
      if (funcI.getRangeSeparationParameter() != funcJ.getRangeSeparationParameter()) throw std::runtime_error("Range separation parameter must be the same in all subsystems");
      hfExchangeRatio = funcI.getHfExchangeRatio();
      lrExchangeRatio = funcI.getLRExchangeRatio();
      mu = funcI.getRangeSeparationParameter();
    }
    //IF a coupled Hartree-Fock calculation is performed only compute naddXC if it is specified in the input
    if(I!=J &&
       this->_lrscf[I]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF &&
       this->_lrscf[I]->getLRSCFSettings().embedding.naddXCFunc == Options::XCFUNCTIONALS::NONE){
         hfExchangeRatio = 0.0;
      }
    //Just defining the variable for lambda functions but will be set equal to hfExchangeRatio/ lrhfExchangeRatio later
    double exchangeRatio = 0.0;

  //Set dimensions for Fock like matrices
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > > > fock(
      new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>(this->_nSet));
  for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }
    //Thread safety
    std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > > *> f;
  #ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    f.push_back(&(*fock));
    for (unsigned int iThread=1; iThread < nThreads; ++iThread) {
      f.push_back(new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > >);
      (*f[iThread]).resize(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
        }
      }
    }
  #else
    f.push_back(&fock);
  #endif
    
    
    //Basisfunctions
    const unsigned int nBFs_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
    const unsigned int nBFs_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();

    // Function to calculate K pseudo-Fock matrix (for I=J)
    auto distributeExchangeII = [&] (
        unsigned int i, unsigned int j, unsigned int k, unsigned int l,
        const Eigen::VectorXd integral, unsigned int threadId){
      //Permutations
  
      const unsigned int ik = i * nBFs_I + k;
      const unsigned int ki = k * nBFs_I + i;
      const unsigned int li = l * nBFs_I + i;
      const unsigned int il = i * nBFs_I + l;
      const unsigned int kj = k * nBFs_I + j;
      const unsigned int jk = j * nBFs_I + k;
      const unsigned int lj = l * nBFs_I + j;
      const unsigned int jl = j * nBFs_I + l;
      
      double perm =1.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;
      perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
      const double exc = perm * integral[0] * exchangeRatio;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        const double exc_ilkj = _pm[iSet] * exc;
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& pp = (*densityMatrices)[iSet][iGuess];
          auto& pf = (*f[threadId])[iSet][iGuess];

          double exc_ik =  *(pp.data() + lj)  *  exc_ilkj;
          double exc_ki =  *(pp.data() + lj)  *  exc;

          exc_ik +=  *(pp.data() + jl)  *  exc;
          exc_ki +=  *(pp.data() + jl)  *  exc_ilkj;

          double exc_li =  *(pp.data() + jk)  *  exc_ilkj;
          double exc_il =  *(pp.data() + jk)  *  exc;

          exc_li +=  *(pp.data() + kj)  *  exc;
          exc_il +=  *(pp.data() + kj)  *  exc_ilkj;

          double exc_kj =  *(pp.data() + il)  *  exc_ilkj;        
          double exc_jk =  *(pp.data() + il)  *  exc;         
          
          exc_kj +=  *(pp.data() + li)  *  exc;        
          exc_jk +=  *(pp.data() + li)  *  exc_ilkj;         
          
          double exc_lj = *(pp.data() + ik)  *  exc_ilkj;         
          double exc_jl =  *(pp.data() + ik)  *  exc;

          exc_lj +=  *(pp.data() + ki)  *  exc; 
          exc_jl +=  *(pp.data() + ki)  *  exc_ilkj;         
                  
          *(pf.data() + ik) += exc_ik;
          *(pf.data() + ki) += exc_ki;
          *(pf.data() + li) += exc_li;
          *(pf.data() + il) += exc_il;
          *(pf.data() + kj) += exc_kj;
          *(pf.data() + jk) += exc_jk;
          *(pf.data() + lj) += exc_lj;
          *(pf.data() + jl) += exc_jl;
/*
          //ikjl
          *(pf.data() + i * nBFs_I + k) += *(pp.data() + l * nBFs_I + j) *  exc_ilkj;
          *(pf.data() + k * nBFs_I + i) += *(pp.data() + j * nBFs_I + l) *  exc_ilkj;
          *(pf.data() + l * nBFs_I + i) += *(pp.data() + j * nBFs_I + k) *  exc_ilkj;
          *(pf.data() + i * nBFs_I + l) += *(pp.data() + k * nBFs_I + j) *  exc_ilkj;
          *(pf.data() + k * nBFs_I + j) += *(pp.data() + i * nBFs_I + l) *  exc_ilkj;
          *(pf.data() + j * nBFs_I + k) += *(pp.data() + l * nBFs_I + i) *  exc_ilkj;
          *(pf.data() + l * nBFs_I + j) += *(pp.data() + i * nBFs_I + k) *  exc_ilkj;
          *(pf.data() + j * nBFs_I + l) += *(pp.data() + k * nBFs_I + i) *  exc_ilkj;

          //ilkj
          *(pf.data() + i * nBFs_I + k) += *(pp.data() + j * nBFs_I + l) *  exc;
          *(pf.data() + k * nBFs_I + i) += *(pp.data() + l * nBFs_I + j) *  exc;
          *(pf.data() + l * nBFs_I + i) += *(pp.data() + k * nBFs_I + j) *  exc;
          *(pf.data() + i * nBFs_I + l) += *(pp.data() + j * nBFs_I + k) *  exc;
          *(pf.data() + k * nBFs_I + j) += *(pp.data() + l * nBFs_I + i) *  exc;
          *(pf.data() + j * nBFs_I + k) += *(pp.data() + i * nBFs_I + l) *  exc;
          *(pf.data() + l * nBFs_I + j) += *(pp.data() + k * nBFs_I + i) *  exc;
          *(pf.data() + j * nBFs_I + l) += *(pp.data() + i * nBFs_I + k) *  exc; */
        }
      }
    };
    
  // Function to calculate K pseudo-Fock matrix (for I!=J)
    auto distributeExchangeIJ = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
      double exc = integral[0] * exchangeRatio;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& pp = (*densityMatrices)[iSet][iGuess];
          auto& pf = (*f[threadId])[iSet][iGuess];
          *(pf.data() + i * nBFs_I + k) +=  *(pp.data() + j * nBFs_J + l) *  exc ;
          *(pf.data() + i * nBFs_I + k) +=  *(pp.data() + l * nBFs_J + j) *  exc * _pm[iSet];
        }
      }
    };

    double maxDens = 0.0;
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto& P =  (*densityMatrices)[iSet][iGuess];
          maxDens = std::max(maxDens,P.array().abs().maxCoeff());
      }
    }

    //Use smalles prescreening threshold of subsystem I and J
    double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
        this->_lrscf[J]->getSysSettings().basis.integralThreshold);

    auto prescreeningFunc = [&](
        unsigned int i, unsigned int j, unsigned int k, unsigned int l,
        unsigned int nI, unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
      (void) i;(void) j;(void) k;(void) l;(void) nI;(void) nJ;(void) nK;(void) nL;
      if (maxDens*schwartz < prescreeningThreshold) return true;
      return false;
    };

    //Calculate pseudo Fock matrices
    if (hfExchangeRatio != 0.0) {
      exchangeRatio = hfExchangeRatio;
      if (I==J) {
        TwoElecFourCenterIntLooper looper(
            libint2::Operator::coulomb,0,this->_lrscf[I]->getBasisController(), prescreeningThreshold);
        looper.loop(distributeExchangeII,prescreeningFunc);
      } else if (I != J) {
        ExchangeInteractionIntLooper looper(
            libint2::Operator::coulomb,
            0,
            this->_lrscf[I]->getBasisController(),
            this->_lrscf[J]->getBasisController(),
            prescreeningThreshold);
        looper.loop(distributeExchangeIJ,prescreeningFunc);
      } else {
        assert(false);
      }
    }

    if (lrExchangeRatio != 0.0) {
      exchangeRatio = lrExchangeRatio;
      if (I==J) {
        TwoElecFourCenterIntLooper looper(
            libint2::Operator::erf_coulomb,0,this->_lrscf[I]->getBasisController(), prescreeningThreshold,mu);
        looper.loop(distributeExchangeII,prescreeningFunc);
      } else if (I != J) {
        ExchangeInteractionIntLooper looper(
            libint2::Operator::erf_coulomb,
            0,
            this->_lrscf[I]->getBasisController(),
            this->_lrscf[J]->getBasisController(),
            prescreeningThreshold,
            mu);
        looper.loop(distributeExchangeIJ,prescreeningFunc);
      } else {
        assert(false);
      }
    }


    //Sum over all threads
    #ifdef _OPENMP
      for (unsigned int iThread=1; iThread<nThreads; ++iThread) {
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess=0; iGuess< this->_nGuess; ++iGuess) {
            (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
          }
        }
        delete f[iThread];
      }
    #endif

      Timings::timeTaken("LRSCF -   Fock-like matrix: K");
    return fock;
  
  }


template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED> > > > KSigmaVector<Options::SCF_MODES::UNRESTRICTED>::calcF(
    unsigned int I,
    unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED> > > > densityMatrices){

    Timings::takeTime("LRSCF -   Fock-like matrix: K");

    //Set exchange ratios
    double hfExchangeRatio = 1.0;
    double lrExchangeRatio = 0.0;
    double mu = 0.0;

    if (this->_lrscf[I]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT ||
        this->_lrscf[I]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      auto funcI = FunctionalClassResolver::resolveFunctional(
          this->_lrscf[I]->getSysSettings().dft.functional);
      auto funcJ = FunctionalClassResolver::resolveFunctional(
          this->_lrscf[J]->getSysSettings().dft.functional);
      //Can only use exact exchange, if the same amount of it is used in every subsystem. Otherwise, the
      //Response matrix becomes non-symmetric.
      if (funcI.getHfExchangeRatio() != funcJ.getHfExchangeRatio()) throw std::runtime_error("HF exchange ratio must be the same in all subsystems");
      if (funcI.getLRExchangeRatio() != funcJ.getLRExchangeRatio()) throw std::runtime_error("HF long range exchange ratio must be the same in all subsystems");
      if (funcI.getRangeSeparationParameter() != funcJ.getRangeSeparationParameter()) throw std::runtime_error("Range separation parameter must be the same in all subsystems");
      hfExchangeRatio = funcI.getHfExchangeRatio();
      lrExchangeRatio = funcI.getLRExchangeRatio();
      mu = funcI.getRangeSeparationParameter();
    }
    double exchangeRatio = 0;

  //Set dimensions for Fock like matrices
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED> > > > fock(
      new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>(this->_nSet));
  for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }
  //Thread safety
    std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED> > > *> f;
  #ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    f.push_back(&(*fock));
    for (unsigned int iThread=1; iThread < nThreads; ++iThread) {
      f.push_back(new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED> > >);
      (*f[iThread]).resize(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
        }
      }
    }
  #else
    f.push_back(&fock);
  #endif
    
    
    //Basisfunctions
    const unsigned int nBFs_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
    const unsigned int nBFs_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();

  // Function to calculate K pseudo-Fock matrix (for I=J)
    auto distributeExchangeII = [&] (
        unsigned int i, unsigned int j, unsigned int k, unsigned int l,
        const Eigen::VectorXd integral, unsigned int threadId){
      //Permutations
      double perm =1.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;
      perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
      const double exc = perm * integral[0] * exchangeRatio;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        const double exc_ilkj = _pm[iSet] * exc;
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& pp = (*densityMatrices)[iSet][iGuess];
          auto& pf = (*f[threadId])[iSet][iGuess];
          //ikjl
          *(pf.alpha.data() +  i * nBFs_I +k) += *(pp.alpha.data() +  l * nBFs_I +j) *  exc_ilkj;
          *(pf.alpha.data() +  k * nBFs_I +i) += *(pp.alpha.data() +  j * nBFs_I +l) *  exc_ilkj;
          *(pf.alpha.data() +  l * nBFs_I +i) += *(pp.alpha.data() +  j * nBFs_I +k) *  exc_ilkj;
          *(pf.alpha.data() +  i * nBFs_I +l) += *(pp.alpha.data() +  k * nBFs_I +j) *  exc_ilkj;
          *(pf.alpha.data() +  k * nBFs_I +j) += *(pp.alpha.data() +  i * nBFs_I +l) *  exc_ilkj;
          *(pf.alpha.data() +  j * nBFs_I +k) += *(pp.alpha.data() +  l * nBFs_I +i) *  exc_ilkj;
          *(pf.alpha.data() +  l * nBFs_I +j) += *(pp.alpha.data() +  i * nBFs_I +k) *  exc_ilkj;
          *(pf.alpha.data() +  j * nBFs_I +l) += *(pp.alpha.data() +  k * nBFs_I +i) *  exc_ilkj;

          *(pf.beta.data() +  i * nBFs_I +k) += *(pp.beta.data() +  l * nBFs_I +j) *  exc_ilkj;
          *(pf.beta.data() +  k * nBFs_I +i) += *(pp.beta.data() +  j * nBFs_I +l) *  exc_ilkj;
          *(pf.beta.data() +  l * nBFs_I +i) += *(pp.beta.data() +  j * nBFs_I +k) *  exc_ilkj;
          *(pf.beta.data() +  i * nBFs_I +l) += *(pp.beta.data() +  k * nBFs_I +j) *  exc_ilkj;
          *(pf.beta.data() +  k * nBFs_I +j) += *(pp.beta.data() +  i * nBFs_I +l) *  exc_ilkj;
          *(pf.beta.data() +  j * nBFs_I +k) += *(pp.beta.data() +  l * nBFs_I +i) *  exc_ilkj;
          *(pf.beta.data() +  l * nBFs_I +j) += *(pp.beta.data() +  i * nBFs_I +k) *  exc_ilkj;
          *(pf.beta.data() +  j * nBFs_I +l) += *(pp.beta.data() +  k * nBFs_I +i) *  exc_ilkj;


          //ilkj
          *(pf.alpha.data() +  i * nBFs_I +k) += *(pp.alpha.data() +  j * nBFs_I +l) *  exc;
          *(pf.alpha.data() +  k * nBFs_I +i) += *(pp.alpha.data() +  l * nBFs_I +j) *  exc;
          *(pf.alpha.data() +  l * nBFs_I +i) += *(pp.alpha.data() +  k * nBFs_I +j) *  exc;
          *(pf.alpha.data() +  i * nBFs_I +l) += *(pp.alpha.data() +  j * nBFs_I +k) *  exc;
          *(pf.alpha.data() +  k * nBFs_I +j) += *(pp.alpha.data() +  l * nBFs_I +i) *  exc;
          *(pf.alpha.data() +  j * nBFs_I +k) += *(pp.alpha.data() +  i * nBFs_I +l) *  exc;
          *(pf.alpha.data() +  l * nBFs_I +j) += *(pp.alpha.data() +  k * nBFs_I +i) *  exc;
          *(pf.alpha.data() +  j * nBFs_I +l) += *(pp.alpha.data() +  i * nBFs_I +k) *  exc;

          *(pf.beta.data() +  i * nBFs_I +k) += *(pp.beta.data() +  j * nBFs_I +l) *  exc;
          *(pf.beta.data() +  k * nBFs_I +i) += *(pp.beta.data() +  l * nBFs_I +j) *  exc;
          *(pf.beta.data() +  l * nBFs_I +i) += *(pp.beta.data() +  k * nBFs_I +j) *  exc;
          *(pf.beta.data() +  i * nBFs_I +l) += *(pp.beta.data() +  j * nBFs_I +k) *  exc;
          *(pf.beta.data() +  k * nBFs_I +j) += *(pp.beta.data() +  l * nBFs_I +i) *  exc;
          *(pf.beta.data() +  j * nBFs_I +k) += *(pp.beta.data() +  i * nBFs_I +l) *  exc;
          *(pf.beta.data() +  l * nBFs_I +j) += *(pp.beta.data() +  k * nBFs_I +i) *  exc;
          *(pf.beta.data() +  j * nBFs_I +l) += *(pp.beta.data() +  i * nBFs_I +k) *  exc;
        }
      }
    };
    
  // Function to calculate K pseudo-Fock matrix (for I!=J)
    auto distributeExchangeIJ = [&] (
      unsigned int i, unsigned int j, unsigned int k, unsigned int l,
      const Eigen::VectorXd integral, unsigned int threadId){
      double exc = integral[0] * exchangeRatio;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& pp = (*densityMatrices)[iSet][iGuess];
          auto& pf = (*f[threadId])[iSet][iGuess];
          *(pf.alpha.data() + i * nBFs_I + k) +=  *(pp.alpha.data() + j * nBFs_J + l) *  exc ;
          *(pf.alpha.data() + i * nBFs_I + k) +=  *(pp.alpha.data() + l * nBFs_J + j) *  exc * _pm[iSet];
          *(pf.beta.data() + i * nBFs_I + k) +=  *(pp.beta.data() + j * nBFs_J + l) *  exc ;
          *(pf.beta.data() + i * nBFs_I + k) +=  *(pp.beta.data() + l * nBFs_J + j) *  exc * _pm[iSet];
        }
      }
    };

    double maxDens = 0.0;
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto& P =  (*densityMatrices)[iSet][iGuess];
        for_spin(P) {
          maxDens = std::max(maxDens,P_spin.array().abs().maxCoeff());
        };
      }
    }

    //Use smalles prescreening threshold of subsystem I and J
    double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
        this->_lrscf[J]->getSysSettings().basis.integralThreshold);

    auto prescreeningFunc = [&](
        unsigned int i, unsigned int j, unsigned int k, unsigned int l,
        unsigned int nI, unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
      (void) i;(void) j;(void) k;(void) l;(void) nI;(void) nJ;(void) nK;(void) nL;
      if (maxDens*schwartz < prescreeningThreshold) return true;
      return false;
    };

    //Calculate pseudo Fock matrices
    if (hfExchangeRatio != 0.0) {
      exchangeRatio = hfExchangeRatio;
      if (I==J) {
        TwoElecFourCenterIntLooper looper(
            libint2::Operator::coulomb,0,this->_lrscf[I]->getBasisController(), prescreeningThreshold);
        looper.loop(distributeExchangeII,prescreeningFunc);
      } else if (I != J) {
        ExchangeInteractionIntLooper looper(
            libint2::Operator::coulomb,
            0,
            this->_lrscf[I]->getBasisController(),
            this->_lrscf[J]->getBasisController(),
            prescreeningThreshold);
        looper.loop(distributeExchangeIJ,prescreeningFunc);
      } else {
        assert(false);
      }
    }

    if (lrExchangeRatio != 0.0) {
      exchangeRatio = lrExchangeRatio;
      if (I==J) {
        TwoElecFourCenterIntLooper looper(
            libint2::Operator::erf_coulomb,0,this->_lrscf[I]->getBasisController(), prescreeningThreshold,mu);
        looper.loop(distributeExchangeII,prescreeningFunc);
      } else if (I != J) {
        ExchangeInteractionIntLooper looper(
            libint2::Operator::erf_coulomb,
            0,
            this->_lrscf[I]->getBasisController(),
            this->_lrscf[J]->getBasisController(),
            prescreeningThreshold,
            mu);
        looper.loop(distributeExchangeIJ,prescreeningFunc);
      } else {
        assert(false);
      }
    }


    //Sum over all threads
    #ifdef _OPENMP
      for (unsigned int iThread=1; iThread<nThreads; ++iThread) {
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess=0; iGuess< this->_nGuess; ++iGuess) {
            (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
          }
        }
        delete f[iThread];
      }
    #endif

    Timings::timeTaken("LRSCF -   Fock-like matrix: K");
    return fock;
  }
template class KSigmaVector<Options::SCF_MODES::RESTRICTED> ;
template class KSigmaVector<Options::SCF_MODES::UNRESTRICTED>;
}


