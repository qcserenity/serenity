/**
 * @file LRSCFAnalysis.cpp
 *
 * @date Dec 1, 2016
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
#include "postHF/LRSCF/Analysis/LRSCFAnalysis.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/matrices/DensityMatrix.h"
#include "io/FormattedOutput.h"
#include "integrals/wrappers/Libint.h"
#include "data/matrices/MatrixInBasis.h"
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "data/OrbitalController.h"
#include "math/linearAlgebra/Orthogonalization.h"
/* Include Std and External Headers */
#include <iostream>


namespace Serenity {

template<Options::SCF_MODES T>
LRSCFAnalysis<T>::LRSCFAnalysis(
    std::shared_ptr<SystemController> activeSystem,
    std::vector<std::shared_ptr<SystemController> > environmentSystems,
    Eigen::VectorXd& eigenvalues,
    std::vector<Eigen::MatrixXd >& eigenvectors):
    _activeSystem(activeSystem),
    _environmentSystems(environmentSystems),
    _eigenvalues(eigenvalues),
    _eigenvectors(eigenvectors),
    _nOccupied(activeSystem->getNOccupiedOrbitals<T>()),
    _nVirtual(activeSystem->getNVirtualOrbitals<T>()),
    _nMolecularOrbitals(activeSystem->getActiveOrbitalController<T>()->getNOrbitals()){
  assert(activeSystem);
  assert (_eigenvectors.size() <= 2);
  _x = _eigenvectors[0];
  if(_eigenvectors.size() == 2) {
    _y = _eigenvectors[1];
  } else {
    //For TDA, y is zero and not calculated.
    _y = Eigen::MatrixXd::Zero(_eigenvectors[0].rows(),_eigenvectors[0].cols());
  }
}

template<>
void LRSCFAnalysis<Options::SCF_MODES::RESTRICTED>::printStateInfo(const unsigned int iState){
  Eigen::MatrixXd coefficients(_x.rows(),_eigenvalues.rows());
  for (unsigned int i = 0; i < _eigenvalues.rows(); ++i) {
    for (unsigned int j = 0; j < _x.rows(); ++j) {
      coefficients(j,i)= _x(j,i)*_x(j,i) - _y(j,i)*_y(j,i);
    }
  }
  coefficients *= 100;
  //print largest coefficients
  print((std::string)"Excitation energy           : " + _eigenvalues(iState));
  print((std::string)"Excitation energy / eV      : " + (_eigenvalues(iState) * HARTREE_TO_EV));
  print((std::string)"Excitation energy / nm      : " + (HARTREE_TO_NM / _eigenvalues(iState)));
  print((std::string)"Excitation energy / cm^(-1) : " + (_eigenvalues(iState) * HARTREE_TO_OOCM));
  print((std::string) "\n Dominant contributions: \n");
  printTableHead(" occ. orbital   virt. orbital    |c|^2*100 ");
  Eigen::VectorXd coefficientVector = coefficients.col(iState);
  for (unsigned int j = 0, jb = 0; j <  _nOccupied; ++j) {
    for (unsigned int b = _nOccupied; b < _nMolecularOrbitals; ++b, ++jb) {
      if (coefficientVector(jb) > 10.0) {
        std::cout.precision(1);
        std::cout << "     "
            << j + 1 <<  " alpha        "
            << b + 1 <<  " alpha           "
            << coefficientVector(jb) << std::endl;
      }
    }
  }
}

template<>
void LRSCFAnalysis<Options::SCF_MODES::UNRESTRICTED>::printStateInfo(const unsigned int iState){
  Eigen::MatrixXd coefficients(_x.rows(),_eigenvalues.rows());
  for (unsigned int i = 0; i < _eigenvalues.rows(); ++i) {
    for (unsigned int j = 0; j < _x.rows(); ++j) {
      coefficients(j,i)= _x(j,i)*_x(j,i) - _y(j,i)*_y(j,i);
    }
  }
  coefficients *= 100;
  //print largest coefficients
  print((std::string)"Excitation energy           : " + _eigenvalues(iState));
  print((std::string)"Excitation energy / eV      : " + (_eigenvalues(iState) * HARTREE_TO_EV));
  print((std::string)"Excitation energy / nm      : " + (HARTREE_TO_NM / _eigenvalues(iState)));
  print((std::string)"Excitation energy / cm^(-1) : " + (_eigenvalues(iState) * HARTREE_TO_OOCM));
  print((std::string) "\n Dominant contributions: \n");
  printTableHead(" occ. orbital   virt. orbital    |c|^2*100 ");
  Eigen::VectorXd coefficientVector = coefficients.col(iState);
  for (unsigned int j = 0, jb = 0; j < _nOccupied.alpha; ++j) {
    for (unsigned int b = _nOccupied.alpha; b < _nMolecularOrbitals; ++b, ++jb) {
      if (coefficientVector(jb) > 10) {
        std::cout.precision(1);
        std::cout << "     "
            << j + 1 <<  " alpha        "
            << b + 1 <<  " alpha           "
            << coefficientVector(jb) << std::endl;
      }
    }
  }
  for (unsigned int j = 0, jb = _nOccupied.alpha * _nVirtual.alpha; j < _nOccupied.beta; ++j) {
    for (unsigned int b = _nOccupied.beta; b < _nMolecularOrbitals; ++b, ++jb) {
      if (coefficientVector(jb) > 10) {
        std::cout.precision(1);
        std::cout << "     "
            << j + 1 <<  " beta         "
            << b + 1 <<  " beta            "
            << coefficientVector(jb) << std::endl;
      }
    }
  }
}


template<>
Eigen::VectorXd LRSCFAnalysis<Options::SCF_MODES::RESTRICTED>::Mo2Ao(Eigen::VectorXd moVec) {
  unsigned int nOcc = _activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  unsigned int nVirt = _activeSystem->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  //Write vector to matrix
  Eigen::MatrixXd tmp(nOcc,nVirt);
  tmp.setZero();
  for (unsigned int i = 0, ia = 0; i < nOcc; i++) {
    for (unsigned int a = 0; a < nVirt; a++, ++ia) {
      tmp(i, a) = moVec(ia);
    }
  }
  //Transform matrix to AO basis
  auto coeff = _activeSystem->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  tmp = coeff.block(0,0,coeff.rows(),nOcc) * tmp * coeff.block(0,nOcc,coeff.rows(),nVirt).transpose();
  //Build AO excitation vector
  Eigen::VectorXd aoVec(coeff.rows() * coeff.rows());
  aoVec.setZero();
  for (unsigned int mu = 0, munu = 0; mu < coeff.rows(); ++mu) {
    for (unsigned int nu = 0; nu < coeff.rows(); ++nu, ++munu) {
      aoVec(munu) = tmp(mu,nu);
    }
  }
  return aoVec;
}

template<>
Eigen::VectorXd LRSCFAnalysis<Options::SCF_MODES::UNRESTRICTED>::Mo2Ao(Eigen::VectorXd moVec) {
  //Get number of occupied and virtual orbitals
  auto nOcc = _activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  auto nVirt = _activeSystem->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  //Write vector to matrix
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> tmp;
  tmp.alpha.resize(nOcc.alpha,nVirt.alpha);
  tmp.alpha.setZero();
  tmp.beta.resize(nOcc.beta,nVirt.beta);
  tmp.beta.setZero();
  for (unsigned int i = 0, ia = 0; i < nOcc.alpha; i++) {
    for (unsigned int a = 0; a < nVirt.alpha; a++, ++ia) {
      tmp.alpha(i, a) = moVec(ia);
    }
  }
  for (unsigned int i = 0, ia = nOcc.alpha * nVirt.alpha; i < nOcc.beta; i++) {
    for (unsigned int a = 0; a < nVirt.beta; a++, ++ia) {
      tmp.beta(i, a) = moVec(ia);
    }
  }
  //Transform transition densities to AO basis
  auto coeff = _activeSystem->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
  tmp.alpha = coeff.alpha.block(0,0,coeff.alpha.rows(),nOcc.alpha) * tmp.alpha * coeff.alpha.block(0,nOcc.alpha,coeff.alpha.rows(),nVirt.alpha).transpose();
  tmp.beta = coeff.beta.block(0,0,coeff.beta.rows(),nOcc.beta) * tmp.beta * coeff.beta.block(0,nOcc.beta,coeff.beta.rows(),nVirt.beta).transpose();
  //Build AO excitation vectors
  Eigen::VectorXd aoVec(_activeSystem->getBasisController()->getNBasisFunctions() * _activeSystem->getBasisController()->getNBasisFunctions());
  aoVec.setZero();
  for (unsigned int mu = 0, munu = 0; mu < coeff.alpha.rows(); ++mu) {
    for (unsigned int nu = 0; nu < coeff.alpha.rows(); ++nu, ++munu) {
      aoVec(munu) = tmp.alpha(mu,nu) + tmp.beta(mu,nu);
    }
  }
  return aoVec;
}


template<>
void LRSCFAnalysis<Options::SCF_MODES::RESTRICTED>::mullikenPopulationAnalysis(const unsigned int iState){
  Eigen::MatrixXd xpy = _x+_y;
  //Get number of occupied and virtual orbitals
  unsigned int nOcc = _activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  unsigned int nVirt = _activeSystem->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  //Store transition density in matrix
  Eigen::MatrixXd tmp(nOcc,nVirt);
  for (unsigned int i = 0, ia = 0; i < nOcc; ++i) {
    for (unsigned int a = 0; a < nVirt; ++a, ++ia) {
      tmp(i,a) = xpy(ia,iState);
    }
  }
  //Build density matrix for excited electron
  DensityMatrix<Options::SCF_MODES::RESTRICTED> dElec(_activeSystem->getBasisController());
  dElec = tmp.transpose() * tmp;
  //Build density matrix for hole
  DensityMatrix<Options::SCF_MODES::RESTRICTED> dHole(_activeSystem->getBasisController());
  dHole = tmp * tmp.transpose();
  //Transform into AO basis
  auto coeff = _activeSystem->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  dElec = -1.0 * coeff.block(0,nOcc,coeff.rows(),nVirt) * dElec * coeff.block(0,nOcc,coeff.rows(),nVirt).transpose();
  dHole = coeff.block(0,0,coeff.rows(),nOcc) * dHole * coeff.block(0,0,coeff.rows(),nOcc).transpose();

  //Perform Mulliken analysis
  auto overlap = _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals();
  auto popElec = MullikenPopulationCalculator<Options::SCF_MODES::RESTRICTED>::calculateAtomPopulations(dElec,overlap,_activeSystem->getAtomCenteredBasisController()->getBasisIndices());
  auto popHole = MullikenPopulationCalculator<Options::SCF_MODES::RESTRICTED>::calculateAtomPopulations(dHole,overlap,_activeSystem->getAtomCenteredBasisController()->getBasisIndices());
  const auto& atoms = _activeSystem->getAtoms();
  unsigned int nAtoms = atoms.size();
  print((std::string) "\n Particle/Hole Population Analysis (Mulliken): \n");
  printf(" Transferred charge: %+8.5f \n",(dElec * overlap).trace());
  printf("\n");
  printTableHead("     No.      Atom   Electron    Hole       Change");
  for (unsigned int i=0; i<nAtoms; ++i) {
    printf( "%4s %5d %9s %11f %11f %11f\n",
        "",(i+1),atoms[i]->getAtomType()->getName().c_str(),popElec[i],popHole[i],popElec[i]+popHole[i]);
  }
  printf("\n");
}

template<>
void LRSCFAnalysis<Options::SCF_MODES::UNRESTRICTED>::mullikenPopulationAnalysis(const unsigned int iState){
  Eigen::MatrixXd xpy = _x+_y;
  //Get number of occupied and virtual orbitals
  auto nOcc = _activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  auto nVirt = _activeSystem->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  //Store transition density in matrix
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> tmp;
  tmp.alpha.resize(nOcc.alpha,nVirt.alpha);
  tmp.alpha.setZero();
  tmp.beta.resize(nOcc.beta,nVirt.beta);
  tmp.beta.setZero();
  for (unsigned int i = 0, ia = 0; i < nOcc.alpha; i++) {
    for (unsigned int a = 0; a < nVirt.alpha; a++, ++ia) {
      tmp.alpha(i, a) = xpy(ia,iState);
    }
  }
  for (unsigned int i = 0, ia = nOcc.alpha * nVirt.alpha; i < nOcc.beta; i++) {
    for (unsigned int a = 0; a < nVirt.beta; a++, ++ia) {
      tmp.beta(i, a) = xpy(ia,iState);
    }
  }
  //Build density matrix for excited electron
  DensityMatrix<Options::SCF_MODES::UNRESTRICTED> dElec(_activeSystem->getBasisController());
  dElec.alpha = tmp.alpha.transpose() * tmp.alpha;
  dElec.beta = tmp.beta.transpose() * tmp.beta;
  //Build density matrix for hole
  DensityMatrix<Options::SCF_MODES::UNRESTRICTED> dHole(_activeSystem->getBasisController());
  dHole.alpha = tmp.alpha * tmp.alpha.transpose();
  dHole.beta = tmp.beta * tmp.beta.transpose();
  //Transform into AO basis
  auto coeff = _activeSystem->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
  dElec.alpha = -1.0 * coeff.alpha.block(0,nOcc.alpha,coeff.alpha.rows(),nVirt.alpha) * dElec.alpha * coeff.alpha.block(0,nOcc.alpha,coeff.alpha.rows(),nVirt.alpha).transpose();
  dElec.beta = -1.0 * coeff.beta.block(0,nOcc.beta,coeff.beta.rows(),nVirt.beta) * dElec.beta * coeff.beta.block(0,nOcc.beta,coeff.beta.rows(),nVirt.beta).transpose();
  dHole.alpha = coeff.alpha.block(0,0,coeff.alpha.rows(),nOcc.alpha) * dHole.alpha * coeff.alpha.block(0,0,coeff.alpha.rows(),nOcc.alpha).transpose();
  dHole.beta = coeff.beta.block(0,0,coeff.beta.rows(),nOcc.beta) * dHole.beta * coeff.beta.block(0,0,coeff.beta.rows(),nOcc.beta).transpose();

  //Perform Mulliken analysis
  auto overlap = _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals();
  auto popElec = MullikenPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>::calculateAtomPopulations(dElec,overlap,_activeSystem->getAtomCenteredBasisController()->getBasisIndices());
  auto popHole = MullikenPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>::calculateAtomPopulations(dHole,overlap,_activeSystem->getAtomCenteredBasisController()->getBasisIndices());
  const auto& atoms = _activeSystem->getAtoms();
  unsigned int nAtoms = atoms.size();
  print((std::string) "\n Particle/Hole Population Analysis (Mulliken): \n");
  printf(" Transferred charge: %+8.5f \n",(dElec.alpha * overlap).trace() + (dElec.beta * overlap).trace());
  printf("\n");
  printTableHead("     No.      Atom   Electron    Hole       Total   ");
  for (unsigned int i=0; i<nAtoms; ++i) {
    printf( "%4s %5d %9s %11f %11f %11f\n",
        "",(i+1),atoms[i]->getAtomType()->getName().c_str(),popElec.alpha[i] + popElec.beta[i],
        popHole.alpha[i] + popHole.beta [i],popElec.alpha[i]+popHole.alpha[i] + popElec.beta[i]+popHole.beta[i]);
  }
  printf("\n");
}




template<Options::SCF_MODES T>
void LRSCFAnalysis<T>::writeExcitationVectors(Eigen::MatrixXd& aoVecs, std::string fname) {
  FILE * f;
  f = fopen (fname.c_str(),"w");
  auto& atoms = _activeSystem->getAtoms();
  int shellIndexI = -1;
  for (unsigned int iAt = 0; iAt < atoms.size(); ++iAt) {
    auto iBasis = atoms[iAt]->getBasisFunctions();
    for (unsigned int iShell = 0; iShell < iBasis.size(); ++iShell) {
      shellIndexI += 1;
      unsigned int firstIndexI = _activeSystem->getBasisController()->extendedIndex(shellIndexI);
      for (unsigned int iM = 0; iM < iBasis[iShell]->getNContracted(); ++iM) {
        unsigned int bfIndexI = firstIndexI + iM;

        int shellIndexJ = -1;
        for (unsigned int jAt = 0; jAt < atoms.size(); ++jAt) {
          auto jBasis = atoms[jAt]->getBasisFunctions();
          for (unsigned int jShell = 0; jShell < jBasis.size(); ++jShell) {
            shellIndexJ += 1;
            unsigned int firstIndexJ = _activeSystem->getBasisController()->extendedIndex(shellIndexJ);
            for (unsigned int jM = 0; jM < jBasis[jShell]->getNContracted(); ++jM) {
              unsigned int bfIndexJ = firstIndexJ + jM;
              unsigned int ij = bfIndexI * _activeSystem->getBasisController()->getNBasisFunctions() + bfIndexJ;
              fprintf(f,"%4i %4i %4i %4i %4i %4i |",iAt,iShell,iM,jAt,jShell,jM);
              for (unsigned int iState = 0; iState < aoVecs.cols(); ++iState) {
                fprintf(f," %+8.5f",aoVecs(ij,iState));
              }
              fprintf(f,"\n");
            }
          }
        }

      }
    }
  }


}

template<Options::SCF_MODES T>
void LRSCFAnalysis<T>::printAOExcitationVectors() {
  Eigen::MatrixXd aoVecsX(_activeSystem->getBasisController()->getNBasisFunctions() * _activeSystem->getBasisController()->getNBasisFunctions(),_eigenvectors[0].cols());
  //Transform to AO basis
  for (unsigned int iState = 0; iState < aoVecsX.cols(); ++iState) {
    aoVecsX.col(iState) = Mo2Ao(_x.col(iState));
  }

  Eigen::MatrixXd aoVecsY(_activeSystem->getBasisController()->getNBasisFunctions() * _activeSystem->getBasisController()->getNBasisFunctions(),_eigenvectors[0].cols());
  //Transform to AO basis
  for (unsigned int iState = 0; iState < aoVecsY.cols(); ++iState) {
    aoVecsY.col(iState) = Mo2Ao(_y.col(iState));
  }

  //Remove nan's for analysis
  auto nans = _eigenvalues.array().isNaN();
  Eigen::MatrixXd tmpX(aoVecsX.rows(),aoVecsX.cols() - nans.sum());
  Eigen::MatrixXd tmpY(aoVecsX.rows(),aoVecsX.cols() - nans.sum());
  int index = -1;
  for (unsigned int i = 0; i < aoVecsX.cols(); ++i) {
    if (nans(i) == 1) continue;
    index += 1;
    tmpX.col(index) = aoVecsX.col(i);
    tmpY.col(index) = aoVecsY.col(i);
  }

  //Orthogonalize
  Orthogonalization::modifiedGramSchmidt(tmpX);
  if (_eigenvectors.size() == 2) {
    //if TDDFT:
    Orthogonalization::modifiedGramSchmidt(tmpY);
    Orthogonalization::modifiedGramSchmidtBiorthogonalization(tmpX, tmpY);
  }

  //Normalize
  auto& x = tmpX;
  auto& y = tmpY;
  Eigen::MatrixXd tmp = x.transpose()*x - y.transpose()*y;
  Eigen::VectorXd norms = tmp.colwise().norm();
  for (unsigned int i = 0; i < norms.rows(); ++i) {
    norms(i) = 1.0 / std::sqrt(norms(i));
  }
  x *= norms.asDiagonal();
  y *= norms.asDiagonal();

  index = -1;
  for (unsigned int i = 0; i < aoVecsX.cols(); ++i) {
    if (nans(i) == 1) continue;
    index += 1;
    aoVecsX.col(i) = tmpX.col(index);
    aoVecsY.col(i) = tmpY.col(index);
  }

  std::string fname = _activeSystem->getSettings().path + "lrscf.X";
  writeExcitationVectors(aoVecsX, fname);
  fname = _activeSystem->getSettings().path + "lrscf.Y";
  writeExcitationVectors(aoVecsY, fname);
}


template class LRSCFAnalysis<Options::SCF_MODES::RESTRICTED> ;
template class LRSCFAnalysis<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
