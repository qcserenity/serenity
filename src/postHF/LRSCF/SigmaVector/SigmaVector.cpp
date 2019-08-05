/**
 * @file SigmaVector.cpp
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
#include "postHF/LRSCF/SigmaVector/SigmaVector.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"

namespace Serenity {

template<Options::SCF_MODES T>
SigmaVector<T>::SigmaVector(
    std::shared_ptr<SystemController> system,
    Eigen::MatrixXd& guessVector):
  _system(system),
  _guessVector(guessVector),
  _hasBeenCalculated(false){
  _sigmaVector.resize(_guessVector.rows(),_guessVector.cols());
  _sigmaVector.setZero();
}

template<>
Eigen::VectorXd SigmaVector<Options::SCF_MODES::RESTRICTED>::ao2mo(
    SPMatrix<Options::SCF_MODES::RESTRICTED>& pF) {
  auto coeff = _system->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  auto nOccupied = _system->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  auto nVirtual = _system->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  pF = coeff.block(0,0,coeff.rows(),nOccupied).transpose() * pF * coeff.block(0,nOccupied,coeff.rows(),nVirtual);
  pF.transposeInPlace();
  Eigen::VectorXd sigma = Eigen::Map<Eigen::VectorXd>(pF.data(),pF.cols()*pF.rows());
  pF.resize(0,0);
  return  sigma;
}

template<>
Eigen::VectorXd SigmaVector<Options::SCF_MODES::UNRESTRICTED>::ao2mo(
    SPMatrix<Options::SCF_MODES::UNRESTRICTED>& pF) {
  auto coeff = _system->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
  auto nOccupied = _system->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  auto nVirtual = _system->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  for_spin(pF,coeff,nOccupied,nVirtual) {
    pF_spin = coeff_spin.block(0,0,coeff_spin.rows(),nOccupied_spin).transpose() * pF_spin* coeff_spin.block(0,nOccupied_spin,coeff_spin.rows(),nVirtual_spin);
    pF_spin.transposeInPlace();
  };
  Eigen::VectorXd sigma(_guessVector.rows());
  sigma.block(0,0,nOccupied.alpha*nVirtual.alpha,1) = Eigen::Map<Eigen::VectorXd>(pF.alpha.data(),pF.alpha.cols()*pF.alpha.rows());
  sigma.block(nOccupied.alpha*nVirtual.alpha,0,nOccupied.beta*nVirtual.beta,1) = Eigen::Map<Eigen::VectorXd>(pF.beta.data(),pF.beta.cols()*pF.beta.rows());
  pF.alpha.resize(0,0);
  pF.beta.resize(0,0);
  return sigma;
}

template<>
SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::MatrixXd> SigmaVector<Options::SCF_MODES::RESTRICTED>::pDens(unsigned int iGuess) {
  auto coeff = _system->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  unsigned int nOccupied = _system->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  unsigned int nVirtual = _system->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  //ToDo: Use Eigen::Map to cast vector into matrix
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::MatrixXd> guess(nOccupied,nVirtual);
  for (unsigned int i = 0, ia = 0; i < nOccupied; ++i) {
    for (unsigned int a = 0; a < nVirtual; ++a, ++ia) {
      guess(i,a) = _guessVector(ia,iGuess);
    }
  }
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::MatrixXd> P;
  P = coeff.block(0,0,coeff.rows(),nOccupied) * guess * coeff.block(0,nOccupied,coeff.rows(),nVirtual).transpose();
  return P;
}

template<>
SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> SigmaVector<Options::SCF_MODES::UNRESTRICTED>::pDens(unsigned int iGuess) {
  auto coeff = _system->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
  auto nOccupied = _system->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  auto nVirtual = _system->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> guess;
  //ToDo: Use Eigen::Map to cast vector into matrix
  guess.alpha.resize(nOccupied.alpha,nVirtual.alpha);
  for (unsigned int i = 0, ia = 0; i < nOccupied.alpha; ++i) {
    for (unsigned int a = 0; a < nVirtual.alpha; ++a, ++ia) {
      guess.alpha(i,a) = _guessVector(ia,iGuess);
    }
  }
  guess.beta.resize(nOccupied.beta,nVirtual.beta);
  for (unsigned int i = 0, ia = nOccupied.alpha * nVirtual.alpha; i < nOccupied.beta; ++i) {
    for (unsigned int a = 0; a < nVirtual.beta; ++a, ++ia) {
      guess.beta(i,a) = _guessVector(ia,iGuess);
    }
  }
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> P;
  for_spin(P,coeff,nOccupied,nVirtual,guess) {
    P_spin = coeff_spin.block(0,0,coeff_spin.rows(),nOccupied_spin) * guess_spin * coeff_spin.block(0,nOccupied_spin,coeff_spin.rows(),nVirtual_spin).transpose();
  };
  return P;
}

template<Options::SCF_MODES T>
void SigmaVector<T>::calculate() {
  //ToDo: Check if data can be kept in memory, otherwise calculate blockwise!
  std::vector<SpinPolarizedData<T,Eigen::MatrixXd> > dens(_guessVector.cols());
  for (unsigned int iGuess = 0; iGuess < _guessVector.cols(); ++iGuess) {
    dens[iGuess] = pDens(iGuess);
  }
  Eigen::MatrixXd block = calculateBlock(_guessVector,dens);
  _sigmaVector.block(0,0,block.rows(),block.cols()) = block;
  _hasBeenCalculated = true;
}


template class SigmaVector<Options::SCF_MODES::RESTRICTED> ;
template class SigmaVector<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
