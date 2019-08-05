/**
 * @file DeltaESigmaVector.cpp
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
#include "postHF/LRSCF/SigmaVector/DeltaESigmaVector.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"

namespace Serenity {

template<Options::SCF_MODES T> DeltaESigmaVector<T>::DeltaESigmaVector(
    std::shared_ptr<SystemController> system,
    Eigen::MatrixXd& guessVector):
  SigmaVector<T>(system,guessVector){

}

template<>
Eigen::MatrixXd DeltaESigmaVector<Options::SCF_MODES::RESTRICTED>::calculateBlock(
    Eigen::MatrixXd& guess,
    std::vector<SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::MatrixXd> >& dens) {
  (void) dens;
  auto orbitalEnergies = _system->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getEigenvalues();
  auto nOccupied = _system->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  auto nVirtual = _system->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  //Calculate orbital energy differences, i.e. the leading diagonal contribution
  Eigen::VectorXd orbitalEnergyDifferences(_guessVector.rows());
  for (unsigned int ia = 0; ia < _guessVector.rows(); ++ia) {
    unsigned int i = floor(ia/nVirtual);
    unsigned int a = nOccupied + ia - i*nVirtual;
    orbitalEnergyDifferences(ia) = orbitalEnergies(a) - orbitalEnergies(i);
  }
  //Contract with guess vectors
  Eigen::MatrixXd sigma(_guessVector.rows(),guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    sigma.col(iGuess) = guess.col(iGuess).cwiseProduct(orbitalEnergyDifferences);
  }
  return sigma;
}

template<>
Eigen::MatrixXd DeltaESigmaVector<Options::SCF_MODES::UNRESTRICTED>::calculateBlock(
    Eigen::MatrixXd& guess,
    std::vector<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> >& dens) {
  (void) dens;
  auto orbitalEnergies = _system->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getEigenvalues();
  auto nOccupied = _system->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  auto nVirtual = _system->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  //Calculate orbital energy differences, i.e. the leading diagonal contribution
  unsigned int nAlpha = nOccupied.alpha*nVirtual.alpha;
  unsigned int nBeta = nOccupied.beta*nVirtual.beta;
  Eigen::VectorXd orbitalEnergyDifferences(_guessVector.rows());
  for (unsigned int ia = 0; ia < nAlpha; ++ia) {
    unsigned int i = floor(ia/nVirtual.alpha);
    unsigned int a = nOccupied.alpha + ia - i*nVirtual.alpha;
    orbitalEnergyDifferences(ia) = orbitalEnergies.alpha(a) - orbitalEnergies.alpha(i);
  }
  for (unsigned int ia = 0; ia < nBeta; ++ia) {
    unsigned int i = floor(ia/nVirtual.beta);
    unsigned int a = nOccupied.beta + ia - i*nVirtual.beta;
    orbitalEnergyDifferences(nAlpha + ia) = orbitalEnergies.beta(a) - orbitalEnergies.beta(i);
  }
  //Contract with guess vectors
  Eigen::MatrixXd sigma(_guessVector.rows(),guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    sigma.col(iGuess) = guess.col(iGuess).cwiseProduct(orbitalEnergyDifferences);
  }
  return sigma;
}

template class DeltaESigmaVector<Options::SCF_MODES::RESTRICTED> ;
template class DeltaESigmaVector<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
