/**
 * @file ExcitationSpectrum.cpp
 *
 * @date Nov 02, 2017
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
#include "postHF/LRSCF/Analysis/ExcitationSpectrum.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "io/FormattedOutput.h"
#include "integrals/wrappers/Libint.h"
#include "data/OrbitalController.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
/* Include Std and External Headers */
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES T> ExcitationSpectrum<T>::ExcitationSpectrum(
    std::shared_ptr<SystemController> systemController,
    std::vector<Eigen::MatrixXd >& eigenvectors,
    Eigen::VectorXd& eigenvalues):
      _systemController(systemController),
      _eigenvectors(eigenvectors),
      _eigenvalues(eigenvalues){

}

template<>
Eigen::MatrixXd ExcitationSpectrum<Options::SCF_MODES::RESTRICTED>::ao2mo(std::vector<Eigen::MatrixXd>& ao_xyz) {
  //Transform integrals
  auto coeff = _systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  std::vector<Eigen::MatrixXd> mo_xyz(3);
  auto nOccupied = _systemController->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  auto nVirtual = _systemController->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  const unsigned int nDimension = nOccupied * nVirtual;
  mo_xyz[0] = coeff.block(0,0,coeff.rows(),nOccupied).transpose() * ao_xyz[0] * coeff.block(0,nOccupied,coeff.rows(),nVirtual);
  mo_xyz[1] = coeff.block(0,0,coeff.rows(),nOccupied).transpose() * ao_xyz[1] * coeff.block(0,nOccupied,coeff.rows(),nVirtual);
  mo_xyz[2] = coeff.block(0,0,coeff.rows(),nOccupied).transpose() * ao_xyz[2] * coeff.block(0,nOccupied,coeff.rows(),nVirtual);

  //Store in vectors
  Eigen::MatrixXd result(nDimension,3);
  mo_xyz[0].transposeInPlace();
  mo_xyz[1].transposeInPlace();
  mo_xyz[2].transposeInPlace();
  result.col(0) = Eigen::Map<Eigen::VectorXd>(mo_xyz[0].data(),mo_xyz[0].cols()*mo_xyz[0].rows());
  result.col(1) = Eigen::Map<Eigen::VectorXd>(mo_xyz[1].data(),mo_xyz[1].cols()*mo_xyz[1].rows());
  result.col(2) = Eigen::Map<Eigen::VectorXd>(mo_xyz[2].data(),mo_xyz[2].cols()*mo_xyz[2].rows());

  return result;
}

template<>
Eigen::MatrixXd ExcitationSpectrum<Options::SCF_MODES::UNRESTRICTED>::ao2mo(std::vector<Eigen::MatrixXd>& ao_xyz) {
  //Transform integrals
  auto coeff = _systemController->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
  std::vector<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> > mo_xyz(3);
  auto nOccupied = _systemController->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  auto nVirtual = _systemController->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  const unsigned int nDimension = nOccupied.alpha * nVirtual.alpha + nOccupied.beta * nVirtual.beta;
  mo_xyz[0].alpha = coeff.alpha.block(0,0,coeff.alpha.rows(),nOccupied.alpha).transpose() * ao_xyz[0] * coeff.alpha.block(0,nOccupied.alpha,coeff.alpha.rows(),nVirtual.alpha);
  mo_xyz[1].alpha = coeff.alpha.block(0,0,coeff.alpha.rows(),nOccupied.alpha).transpose() * ao_xyz[1] * coeff.alpha.block(0,nOccupied.alpha,coeff.alpha.rows(),nVirtual.alpha);
  mo_xyz[2].alpha = coeff.alpha.block(0,0,coeff.alpha.rows(),nOccupied.alpha).transpose() * ao_xyz[2] * coeff.alpha.block(0,nOccupied.alpha,coeff.alpha.rows(),nVirtual.alpha);
  mo_xyz[0].beta = coeff.beta.block(0,0,coeff.beta.rows(),nOccupied.beta).transpose() * ao_xyz[0] * coeff.beta.block(0,nOccupied.beta,coeff.beta.rows(),nVirtual.beta);
  mo_xyz[1].beta = coeff.beta.block(0,0,coeff.beta.rows(),nOccupied.beta).transpose() * ao_xyz[1] * coeff.beta.block(0,nOccupied.beta,coeff.beta.rows(),nVirtual.beta);
  mo_xyz[2].beta = coeff.beta.block(0,0,coeff.beta.rows(),nOccupied.beta).transpose() * ao_xyz[2] * coeff.beta.block(0,nOccupied.beta,coeff.beta.rows(),nVirtual.beta);

  //Store in vectors
  Eigen::MatrixXd result(nDimension,3);
  mo_xyz[0].alpha.transposeInPlace();
  mo_xyz[1].alpha.transposeInPlace();
  mo_xyz[2].alpha.transposeInPlace();
  mo_xyz[0].beta.transposeInPlace();
  mo_xyz[1].beta.transposeInPlace();
  mo_xyz[2].beta.transposeInPlace();
  result.block(0,0,nOccupied.alpha*nVirtual.alpha,1)= Eigen::Map<Eigen::VectorXd>(mo_xyz[0].alpha.data(),mo_xyz[0].alpha.cols()*mo_xyz[0].alpha.rows());
  result.block(nOccupied.alpha*nVirtual.alpha,0,nOccupied.beta*nVirtual.beta,1)= Eigen::Map<Eigen::VectorXd>(mo_xyz[0].beta.data(),mo_xyz[0].beta.cols()*mo_xyz[0].beta.rows());
  result.block(0,1,nOccupied.alpha*nVirtual.alpha,1)= Eigen::Map<Eigen::VectorXd>(mo_xyz[1].alpha.data(),mo_xyz[1].alpha.cols()*mo_xyz[1].alpha.rows());
  result.block(nOccupied.alpha*nVirtual.alpha,1,nOccupied.beta*nVirtual.beta,1)= Eigen::Map<Eigen::VectorXd>(mo_xyz[1].beta.data(),mo_xyz[1].beta.cols()*mo_xyz[1].beta.rows());
  result.block(0,2,nOccupied.alpha*nVirtual.alpha,1)= Eigen::Map<Eigen::VectorXd>(mo_xyz[2].alpha.data(),mo_xyz[2].alpha.cols()*mo_xyz[2].alpha.rows());
  result.block(nOccupied.alpha*nVirtual.alpha,2,nOccupied.beta*nVirtual.beta,1)= Eigen::Map<Eigen::VectorXd>(mo_xyz[2].beta.data(),mo_xyz[2].beta.cols()*mo_xyz[2].beta.rows());

  return result;
}



template<Options::SCF_MODES T>
Eigen::MatrixXd ExcitationSpectrum<T>::dipoleIntegrals() {
  //Calculate AO integrals
  auto& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::emultipole1,0,2);

  auto basisController = _systemController->getBasisController();
  auto basis = basisController->getBasis();
  const unsigned int nBFs = basisController->getNBasisFunctions();
  std::vector<Eigen::MatrixXd > dipolesAO(3,Eigen::MatrixXd::Zero(nBFs,nBFs));

  for (unsigned int i = 0; i < basis.size(); i++) {
    for (unsigned int j = 0; j < basis.size(); j++) {
      Eigen::MatrixXd dipoleInts;
      if(libint.compute(libint2::Operator::emultipole1,0,*basis[i], *basis[j], dipoleInts)){
        for (unsigned int k = 0; k < basis[i]->getNContracted(); k++) {
          auto mu = basisController->extendedIndex(i) + k;
          for (unsigned int l = 0; l < basis[j]->getNContracted(); l++) {
            auto nu = basisController->extendedIndex(j) + l;
            const unsigned int nj = basis[j]->getNContracted();
            dipolesAO[0](mu,nu) = dipoleInts((nj * k + l),1);
            dipolesAO[1](mu,nu) = dipoleInts((nj * k + l),2);
            dipolesAO[2](mu,nu) = dipoleInts((nj * k + l),3);
          }
        }
      }
    }
  };
  libint.finalize(libint2::Operator::emultipole1,0,2);

  return ao2mo(dipolesAO);
}

template<Options::SCF_MODES T>
Eigen::MatrixXd ExcitationSpectrum<T>::momentumIntegrals() {
  //Calculate AO momentum integrals
  auto gridController = _systemController->getGridController();
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      this->_systemController->getSettings().grid.blocksize,
      this->_systemController->getSettings().grid.basFuncRadialThreshold,
      1,
      this->_systemController->getBasisController(),
      gridController);
  //Get number of grid blocks
  const unsigned int nBlocks = basisFunctionOnGridController->getNBlocks();
  //Get weights
  const auto& weights = gridController->getWeights();
  const unsigned int nBasisFunc = _systemController->getBasisController()->getNBasisFunctions();
  std::vector<Eigen::MatrixXd > momentumAO(3,Eigen::MatrixXd::Zero(nBasisFunc,nBasisFunc));
  //Perform blockwise numerical integration
  for (unsigned int iBlock=0; iBlock < nBlocks; ++iBlock) {
    //Data for this block
    auto& thisBlockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
    //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
    const auto& basisFunctionValues = thisBlockData->functionValues;
    //gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
    const auto& gradBasisFunctionValues = thisBlockData->derivativeValues;
    //basis function negligebility
    const auto& negligible = thisBlockData->negligible;
    //number of grid points in this block
    const unsigned int blockSize = thisBlockData->functionValues.rows();
    //Get first index of this Block
    const unsigned int iGridStart = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
    //Add to matrix
    for (unsigned int nu = 0; nu < nBasisFunc; ++nu) {
      if (negligible[nu]) continue;
      Eigen::VectorXd nuW(basisFunctionValues.col(nu).cwiseProduct(weights.segment(iGridStart,blockSize)));
      for (unsigned int mu = 0; mu < nBasisFunc; ++mu) {
        if (negligible[mu]) continue;
        momentumAO[0](nu,mu) += nuW.dot((*gradBasisFunctionValues).x.col(mu));
        momentumAO[1](nu,mu) += nuW.dot((*gradBasisFunctionValues).y.col(mu));
        momentumAO[2](nu,mu) += nuW.dot((*gradBasisFunctionValues).z.col(mu));
      }
    }
  }

  return ao2mo(momentumAO);
}


template<Options::SCF_MODES T>
void ExcitationSpectrum<T>::printSpectrum() {
  //Calculate dipole integrals
  auto dipoles = dipoleIntegrals();

  //Calculate momentum integrals
  auto momentum = momentumIntegrals();

  //Calculate X+Y and X-Y excitation vector
  Eigen::MatrixXd xpy = _eigenvectors[0];
  Eigen::MatrixXd xmy = _eigenvectors[0];
  if (_eigenvectors.size() ==2) {
    xpy += _eigenvectors[1];
    xmy -= _eigenvectors[1];
  }

  //Calculate transition dipole moments
  Eigen::MatrixXd tdm_l(_eigenvectors[0].cols(),3);
  Eigen::MatrixXd tdm_v(_eigenvectors[0].cols(),3);
  double factor = (T==Options::SCF_MODES::RESTRICTED) ? 1.0 : 1/std::sqrt(2);
  for (unsigned int iState = 0; iState < tdm_l.rows(); ++iState) {
    for (unsigned int ix = 0; ix < 3; ++ix) {
      tdm_l(iState,ix) = factor * -1.0 * std::sqrt(2) * xpy.col(iState).dot(dipoles.col(ix));
      tdm_v(iState,ix) = factor * -1.0 / _eigenvalues(iState) * std::sqrt(2) * xmy.col(iState).dot(momentum.col(ix));
    }
  }

  //Calculate oscillator strengths
  Eigen::VectorXd os_l(tdm_l.rows());
  Eigen::VectorXd os_v(tdm_l.rows());
  for (unsigned int iState = 0; iState < tdm_l.rows(); ++iState) {
    os_l(iState) = 2.0/3.0 * _eigenvalues(iState) * tdm_l.row(iState).array().square().sum();
    os_v(iState) = 2.0/3.0 * _eigenvalues(iState) * tdm_v.row(iState).array().square().sum();
  }
  printBigCaption("LRSCF Excitation Spectrum");
  //Print spectra
  printf("--------------------------------------------------------------------------------\n");
  printf("                      Absorption Spectrum (dipole-length)                       \n");
  printf("--------------------------------------------------------------------------------\n");
  printf("state      energy      wavelength  f_osc      mu**2     mu_x     mu_y     mu_z \n");
  printf("       (eV)    (cm-1)     (nm)               (au**2)    (au)     (au)     (au) \n");
  printf("--------------------------------------------------------------------------------\n");
  for (unsigned int iState = 0; iState < _eigenvalues.rows(); ++iState) {
    printf(" %3i %6.2f %10.1f %8.1f  %10.7f %8.4f %8.4f %8.4f %8.4f \n",
        iState + 1,
        _eigenvalues(iState) * HARTREE_TO_EV,
        _eigenvalues(iState) * HARTREE_TO_OOCM,
        HARTREE_TO_NM / _eigenvalues(iState),
        os_l(iState),
        tdm_l.row(iState).array().square().sum(),
        tdm_l(iState,0),tdm_l(iState,1),tdm_l(iState,2));
  }
  printf("\n--------------------------------------------------------------------------------\n");
  printf("                      Absorption Spectrum (dipole-velocity)                       \n");
  printf("--------------------------------------------------------------------------------\n");
  printf("state      energy      wavelength  f_osc      mu**2     mu_x     mu_y     mu_z \n");
  printf("       (eV)    (cm-1)     (nm)               (au**2)    (au)     (au)     (au) \n");
  printf("--------------------------------------------------------------------------------\n");
  for (unsigned int iState = 0; iState < _eigenvalues.rows(); ++iState) {
    printf(" %3i %6.2f %10.1f %8.1f  %10.7f %8.4f %8.4f %8.4f %8.4f \n",
        iState + 1,
        _eigenvalues(iState) * HARTREE_TO_EV,
        _eigenvalues(iState) * HARTREE_TO_OOCM,
        HARTREE_TO_NM / _eigenvalues(iState),
        os_v(iState),
        tdm_v.row(iState).array().square().sum(),
        tdm_v(iState,0),tdm_v(iState,1),tdm_v(iState,2));
  }
  printf("\n");
}



template class ExcitationSpectrum<Options::SCF_MODES::RESTRICTED> ;
template class ExcitationSpectrum<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
