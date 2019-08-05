/**
 * @file KernelSigmaVector.cpp
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
#include "postHF/LRSCF/SigmaVector/KernelSigmaVector.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/DoublySpinPolarizedData.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/SPMatrix.h"

namespace Serenity {

template<Options::SCF_MODES T> KernelSigmaVector<T>::KernelSigmaVector(
    std::shared_ptr<SystemController> system,
    Eigen::MatrixXd& guessVector,
    std::shared_ptr<Kernel<T> > kernel):
        SigmaVector<T>(system,guessVector),
        _kernel(kernel),
        _gga(false),
        _nAddGGA(false){
  if (kernel) {
    //Get derivatives of integral kernel on grid. If quantity is not calculated for the used functional,
    //the getters return a nullptr.
    _d2FdRho2 = _kernel->getD2F_dRho2();
    _dFdSigma = _kernel->getDF_dSigma();
    _d2FdSigma2 = _kernel->getD2F_dSigma2();
    _d2FdRhodSigma = _kernel->getD2F_dRhodSigma();
    _densityGradient = _kernel->getActiveDensityGradient();
    //derivatives of GGA part for non-additive functionals evaluated for the total density. The non-additive
    //contribution from the active density is contained in the elements above.
    _totalD2FdSigma2 = _kernel->getTotalNaddD2F_dSigma2();
    _totalD2FdRhodSigma = _kernel->getTotalNaddD2F_dRhodSigma();
    _totalDensityGradient = _kernel->getTotalDensityGradient();

    if (_dFdSigma) _gga = true;
    if (_totalD2FdSigma2) _nAddGGA = true;

    //LDA Kernel needed in all cases
    assert(_d2FdRho2);
    //Either all GGA objects have to be there or none
    assert((_dFdSigma && _d2FdSigma2 && _d2FdRhodSigma && _densityGradient) ||
        (!_dFdSigma && !_d2FdSigma2 && !_d2FdRhodSigma && !_densityGradient));
    //Either all FDE objects have to be there or none
    assert((_totalD2FdSigma2 && _totalD2FdRhodSigma && _totalDensityGradient) ||
        (!_totalD2FdSigma2 && !_totalD2FdRhodSigma && !_totalDensityGradient));

    unsigned int derivativeLevel = (_gga || _nAddGGA) ? 1 : 0;
    _basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
        this->_system->getSettings().grid.blocksize,
        this->_system->getSettings().grid.basFuncRadialThreshold,
        derivativeLevel,
        this->_system->getBasisController(),
        _kernel->getGridController());
  }
}

template<> void KernelSigmaVector<Options::SCF_MODES::RESTRICTED>::contractGridBlock(
    unsigned int iBlock,
    const std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockData,
    std::vector<SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::MatrixXd> >& dens,
    std::vector<SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::VectorXd> >& scalarPart,
    std::vector<Gradient<SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::VectorXd> > >* gradientPart) {
  //Factor
  double factor = 4.0;
  //Number of basis functions
  const unsigned int nBasisFunc = this->_system->getBasisController()->getNBasisFunctions();
  //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValues = blockData->functionValues;
  //gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
  const auto& gradBasisFunctionValues = blockData->derivativeValues;
  //basis function negligebility
  const auto& negligible = blockData->negligible;
  //number of grid points in this block
  const unsigned int blockSize = blockData->functionValues.rows();
  //Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
  //set dimensions and null
  for (unsigned int iGuess = 0; iGuess < scalarPart.size(); ++iGuess) {
    scalarPart[iGuess].resize(blockSize);
    scalarPart[iGuess].setZero();
    if (_gga) {
      auto grad = makeGradient<SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::VectorXd> >(blockSize);
      (*gradientPart)[iGuess] = grad;
      (*gradientPart)[iGuess].x.setZero();
      (*gradientPart)[iGuess].y.setZero();
      (*gradientPart)[iGuess].z.setZero();
    }
  }
  if (!_gga) {
    //
    //LDA
    //
    //Loop over all guess vectors
    for (unsigned int iGuess = 0; iGuess < scalarPart.size(); ++iGuess) {
      //Get density Matrix for current guess vector
      const auto& density = dens[iGuess];
      //Premultiply basis functions with pseudo density
      Eigen::MatrixXd pbf(blockSize,nBasisFunc);
      pbf.setZero();
      for (unsigned int kappa = 0; kappa < nBasisFunc; ++kappa) {
        // Simple significance screening
        if (negligible[kappa]) continue;
        for (unsigned int lambda = 0; lambda <= kappa; ++lambda) {
          double perm = (kappa == lambda) ? 0.5 : 1.0;
          double dens = perm * (density(kappa,lambda) + density(lambda,kappa));
          pbf.col(kappa) += dens * basisFunctionValues.col(lambda);
        }
      }
      //Premultiply basisfunctions
      Eigen::VectorXd bfbf(blockSize);
      bfbf.setZero();
      for (unsigned int kappa = 0; kappa < nBasisFunc; ++kappa) {
        bfbf += basisFunctionValues.col(kappa).cwiseProduct(pbf.col(kappa));
      }
      //contraction
      scalarPart[iGuess] = factor * (*_d2FdRho2).block(iGridStart,0,blockSize,1).cwiseProduct(bfbf);
    }
  } else if (_gga || _nAddGGA) {
    //
    //GGA
    //
    //Loop over all guess vectors
    for (unsigned int iGuess = 0; iGuess < scalarPart.size(); ++iGuess) {
      //Get density Matrix for current guess vector
      const auto& density = dens[iGuess];
      //Premultiply basis functions and derivatives with pseudo density
      Eigen::MatrixXd pbf(blockSize,nBasisFunc);
      auto pgbf = makeGradient<Eigen::MatrixXd>(blockSize,nBasisFunc);
      pbf.setZero();
      pgbf.x.setZero();
      pgbf.y.setZero();
      pgbf.z.setZero();
      for (unsigned int kappa = 0; kappa < nBasisFunc; ++kappa) {
        // Simple significance screening
        if (negligible[kappa]) continue;
        for (unsigned int lambda = 0; lambda <= kappa; ++lambda) {
          double perm = (kappa == lambda) ? 0.5 : 1.0;
          double dens = perm * (density(kappa,lambda) + density(lambda,kappa));
          pbf.col(kappa) += dens * basisFunctionValues.col(lambda);
          pgbf.x.col(kappa) += dens * (*gradBasisFunctionValues).x.col(lambda);
          pgbf.y.col(kappa) += dens * (*gradBasisFunctionValues).y.col(lambda);
          pgbf.z.col(kappa) += dens * (*gradBasisFunctionValues).z.col(lambda);
        }
      }
      //Norm based prescreening
      if (pbf.norm() < this->_system->getSettings().grid.blockAveThreshold) continue;
      //Calculate \sum_{\kappa \lambda} P_{\kappa \lambda} (\phi_\kappa \phi_\lambda)
      //and \sum_{\kappa \lambda} P_{\kappa \lambda} \nabla(\phi_\kappa \phi_\lambda)
      Eigen::VectorXd pbfbf(blockSize);
      auto pgbfbf = makeGradient<Eigen::VectorXd>(blockSize);
      pbfbf.setZero();
      pgbfbf.x.setZero();
      pgbfbf.y.setZero();
      pgbfbf.z.setZero();
      for (unsigned int kappa = 0; kappa < nBasisFunc; ++kappa) {
        pbfbf += basisFunctionValues.col(kappa).cwiseProduct(pbf.col(kappa));
        pgbfbf.x += (*gradBasisFunctionValues).x.col(kappa).cwiseProduct(pbf.col(kappa));
        pgbfbf.y += (*gradBasisFunctionValues).y.col(kappa).cwiseProduct(pbf.col(kappa));
        pgbfbf.z += (*gradBasisFunctionValues).z.col(kappa).cwiseProduct(pbf.col(kappa));
        pgbfbf.x += basisFunctionValues.col(kappa).cwiseProduct(pgbf.x.col(kappa));
        pgbfbf.y += basisFunctionValues.col(kappa).cwiseProduct(pgbf.y.col(kappa));
        pgbfbf.z += basisFunctionValues.col(kappa).cwiseProduct(pgbf.z.col(kappa));
      }
      //Norm based prescreening
      if (pbfbf.norm() < this->_system->getSettings().grid.blockAveThreshold) continue;
      //Contract LDA part with kernel
      scalarPart[iGuess] += factor * (*_d2FdRho2).block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf);
      if (_gga) {
        //Calculate ...\nabla \rho (\phi_\kappa \phi_\lambda)
        auto grhopbfbf = makeGradient<Eigen::VectorXd>(blockSize);
        grhopbfbf.x = (*_densityGradient).x.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf);
        grhopbfbf.y = (*_densityGradient).y.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf);
        grhopbfbf.z = (*_densityGradient).z.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf);
        //Calculate ...\nabla \rho \nabla(\phi_\kappa \phi_\lambda)
        Eigen::VectorXd  grhopgbfbf(blockSize);
        grhopgbfbf.setZero();
        grhopgbfbf += (*_densityGradient).x.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x);
        grhopgbfbf += (*_densityGradient).y.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y);
        grhopgbfbf += (*_densityGradient).z.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z);
        //Calculate ...\nabla rho \nabla \rho \nabla(\phi_\kappa \phi_\lambda)
        auto grhogrhopgbfbf = makeGradient<Eigen::VectorXd>(blockSize);
        grhogrhopgbfbf.x = (*_densityGradient).x.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf);
        grhogrhopgbfbf.y = (*_densityGradient).y.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf);
        grhogrhopgbfbf.z = (*_densityGradient).z.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf);
        //Contract with kernel
        scalarPart[iGuess] += 2.0 * factor * (*_d2FdRhodSigma).block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf);
        (*gradientPart)[iGuess].x += 2.0 * factor * (*_d2FdRhodSigma).block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x);
        (*gradientPart)[iGuess].y += 2.0 * factor * (*_d2FdRhodSigma).block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y);
        (*gradientPart)[iGuess].z += 2.0 * factor * (*_d2FdRhodSigma).block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z);
        (*gradientPart)[iGuess].x += 2.0 * factor * (*_dFdSigma).block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x);
        (*gradientPart)[iGuess].y += 2.0 * factor * (*_dFdSigma).block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y);
        (*gradientPart)[iGuess].z += 2.0 * factor * (*_dFdSigma).block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z);
        (*gradientPart)[iGuess].x += 4.0 * factor * (*_d2FdSigma2).block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf.x);
        (*gradientPart)[iGuess].y += 4.0 * factor * (*_d2FdSigma2).block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf.y);
        (*gradientPart)[iGuess].z += 4.0 * factor * (*_d2FdSigma2).block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf.z);
      }
      if (_nAddGGA) {
        //Calculate ...\nabla \rho (\phi_\kappa \phi_\lambda)
        auto grhopbfbf = makeGradient<Eigen::VectorXd>(blockSize);
        grhopbfbf.x = (*_totalDensityGradient).x.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf);
        grhopbfbf.y = (*_totalDensityGradient).y.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf);
        grhopbfbf.z = (*_totalDensityGradient).z.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf);
        //Calculate ...\nabla \rho \nabla(\phi_\kappa \phi_\lambda)
        Eigen::VectorXd  grhopgbfbf(blockSize);
        grhopgbfbf.setZero();
        grhopgbfbf += (*_totalDensityGradient).x.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x);
        grhopgbfbf += (*_totalDensityGradient).y.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y);
        grhopgbfbf += (*_totalDensityGradient).z.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z);
        //Calculate ...\nabla rho \nabla \rho \nabla(\phi_\kappa \phi_\lambda)
        auto grhogrhopgbfbf = makeGradient<Eigen::VectorXd>(blockSize);
        grhogrhopgbfbf.x = (*_totalDensityGradient).x.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf);
        grhogrhopgbfbf.y = (*_totalDensityGradient).y.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf);
        grhogrhopgbfbf.z = (*_totalDensityGradient).z.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf);
        //Contract with kernel
        scalarPart[iGuess] += 2.0 * factor * (*_totalD2FdRhodSigma).block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf);
        (*gradientPart)[iGuess].x += 2.0 * factor * (*_totalD2FdRhodSigma).block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x);
        (*gradientPart)[iGuess].y += 2.0 * factor * (*_totalD2FdRhodSigma).block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y);
        (*gradientPart)[iGuess].z += 2.0 * factor * (*_totalD2FdRhodSigma).block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z);
        (*gradientPart)[iGuess].x += 4.0 * factor * (*_totalD2FdSigma2).block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf.x);
        (*gradientPart)[iGuess].y += 4.0 * factor * (*_totalD2FdSigma2).block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf.y);
        (*gradientPart)[iGuess].z += 4.0 * factor * (*_totalD2FdSigma2).block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf.z);
      }
    }
  } else {
    assert (false);
  }
}

template<> void KernelSigmaVector<Options::SCF_MODES::UNRESTRICTED>::contractGridBlock(
    unsigned int iBlock,
    const std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockData,
    std::vector<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> >& dens,
    std::vector<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> >& scalarPart,
    std::vector<Gradient<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> > >* gradientPart) {
//  Timings::takeTime("LRSCF: Contraction");
  //Number of basis functions
  const unsigned int nBasisFunc = this->_system->getBasisController()->getNBasisFunctions();
  //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValues = blockData->functionValues;
  //gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
  const auto& gradBasisFunctionValues = blockData->derivativeValues;
  //basis function negligebility
  const auto& negligible = blockData->negligible;
  //number of grid points in this block
  const unsigned int blockSize = blockData->functionValues.rows();
  //set dimensions and null
  for (unsigned int iGuess = 0; iGuess < scalarPart.size(); ++iGuess) {
    scalarPart[iGuess].alpha.resize(blockSize);
    scalarPart[iGuess].alpha.setZero();
    scalarPart[iGuess].beta.resize(blockSize);
    scalarPart[iGuess].beta.setZero();
    if (_gga) {
      auto grad = makeGradient<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> >(blockSize);
      (*gradientPart)[iGuess] = grad;
      (*gradientPart)[iGuess].x.alpha.setZero();
      (*gradientPart)[iGuess].x.beta.setZero();
      (*gradientPart)[iGuess].y.alpha.setZero();
      (*gradientPart)[iGuess].y.beta.setZero();
      (*gradientPart)[iGuess].z.alpha.setZero();
      (*gradientPart)[iGuess].z.beta.setZero();
    }
  }
  //Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);

  if (!_gga) {
    //
    //LDA
    //
    //Loop over all guess vectors
    for (unsigned int iGuess = 0; iGuess < scalarPart.size(); ++iGuess) {
      //Get density Matrix for current guess vector
      const auto& density = dens[iGuess];
      //Premultiply basis functions with pseudo density
      SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> pbf(blockSize,nBasisFunc);
      pbf.alpha.setZero();
      pbf.beta.setZero();
      for (unsigned int kappa = 0; kappa < nBasisFunc; ++kappa) {
        // Simple significance screening
        if (negligible[kappa]) continue;
        for (unsigned int lambda = 0; lambda <= kappa; ++lambda) {
          SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,double> dens;
          double perm = (kappa == lambda) ? 0.5 : 1.0;
          dens.alpha = perm * (density.alpha(kappa,lambda) + density.alpha(lambda,kappa));
          dens.beta = perm * (density.beta(kappa,lambda) + density.beta(lambda,kappa));
          pbf.alpha.col(kappa) += dens.alpha * basisFunctionValues.col(lambda);
          pbf.beta.col(kappa) += dens.beta * basisFunctionValues.col(lambda);
        }
      }
      //Premultiply basisfunctions
      SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> bfbf(blockSize);
      bfbf.alpha.setZero();
      bfbf.beta.setZero();
      for (unsigned int kappa = 0; kappa < nBasisFunc; ++kappa) {
        bfbf.alpha += basisFunctionValues.col(kappa).cwiseProduct(pbf.alpha.col(kappa));
        bfbf.beta += basisFunctionValues.col(kappa).cwiseProduct(pbf.beta.col(kappa));
      }
      //contraction
      scalarPart[iGuess].alpha += 2.0 * (*_d2FdRho2).aa.block(iGridStart,0,blockSize,1).cwiseProduct(bfbf.alpha);
      scalarPart[iGuess].alpha += 2.0 * (*_d2FdRho2).ab.block(iGridStart,0,blockSize,1).cwiseProduct(bfbf.beta);
      scalarPart[iGuess].beta += 2.0 * (*_d2FdRho2).ab.block(iGridStart,0,blockSize,1).cwiseProduct(bfbf.alpha);
      scalarPart[iGuess].beta += 2.0 * (*_d2FdRho2).bb.block(iGridStart,0,blockSize,1).cwiseProduct(bfbf.beta);
    }
  } else if (_gga || _nAddGGA) {
    //
    //GGA
    //
    //Loop over all guess vectors
    for (unsigned int iGuess = 0; iGuess < scalarPart.size(); ++iGuess) {
      //Get density Matrix for current guess vector
      const auto& density = dens[iGuess];
      //Premultiply basis functions with pseudo density
      SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> pbf(blockSize,nBasisFunc);
      auto pgbf = makeGradient<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::MatrixXd> >(blockSize,nBasisFunc);
      pbf.alpha.setZero();
      pbf.beta.setZero();
      pgbf.x.alpha.setZero();
      pgbf.x.beta.setZero();
      pgbf.y.alpha.setZero();
      pgbf.y.beta.setZero();
      pgbf.z.alpha.setZero();
      pgbf.z.beta.setZero();
      for (unsigned int kappa = 0; kappa < nBasisFunc; ++kappa) {
        // Simple significance screening
        if (negligible[kappa]) continue;
        for (unsigned int lambda = 0; lambda <= kappa; ++lambda) {
          SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,double> dens;
          double perm = (kappa == lambda) ? 0.5 : 1.0;
          dens.alpha = perm * (density.alpha(kappa,lambda) + density.alpha(lambda,kappa));
          dens.beta = perm * (density.beta(kappa,lambda) + density.beta(lambda,kappa));
          pbf.alpha.col(kappa) += dens.alpha * basisFunctionValues.col(lambda);
          pbf.beta.col(kappa) += dens.beta * basisFunctionValues.col(lambda);
          pgbf.x.alpha.col(kappa) += dens.alpha * (*gradBasisFunctionValues).x.col(lambda);
          pgbf.x.beta.col(kappa) += dens.beta * (*gradBasisFunctionValues).x.col(lambda);
          pgbf.y.alpha.col(kappa) += dens.alpha * (*gradBasisFunctionValues).y.col(lambda);
          pgbf.y.beta.col(kappa) += dens.beta * (*gradBasisFunctionValues).y.col(lambda);
          pgbf.z.alpha.col(kappa) += dens.alpha * (*gradBasisFunctionValues).z.col(lambda);
          pgbf.z.beta.col(kappa) += dens.beta * (*gradBasisFunctionValues).z.col(lambda);
        }
      }
      //Norm based prescreening
      double maxNorm = std::max(pbf.alpha.norm(),pbf.beta.norm());
      if (maxNorm < this->_system->getSettings().grid.blockAveThreshold) continue;
      //Calculate \sum_{\kappa \lambda} P_{\kappa \lambda} (\phi_\kappa \phi_\lambda)
      //and \sum_{\kappa \lambda} P_{\kappa \lambda} \nabla(\phi_\kappa \phi_\lambda)
      SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> pbfbf(blockSize);
      auto pgbfbf = makeGradient<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> >(blockSize);
      pbfbf.alpha.setZero();
      pbfbf.beta.setZero();
      pgbfbf.x.alpha.setZero();
      pgbfbf.x.beta.setZero();
      pgbfbf.y.alpha.setZero();
      pgbfbf.y.beta.setZero();
      pgbfbf.z.alpha.setZero();
      pgbfbf.z.beta.setZero();
      for (unsigned int kappa = 0; kappa < nBasisFunc; ++kappa) {
        pbfbf.alpha += basisFunctionValues.col(kappa).cwiseProduct(pbf.alpha.col(kappa));
        pbfbf.beta += basisFunctionValues.col(kappa).cwiseProduct(pbf.beta.col(kappa));
        pgbfbf.x.alpha += (*gradBasisFunctionValues).x.col(kappa).cwiseProduct(pbf.alpha.col(kappa));
        pgbfbf.x.beta += (*gradBasisFunctionValues).x.col(kappa).cwiseProduct(pbf.beta.col(kappa));
        pgbfbf.y.alpha += (*gradBasisFunctionValues).y.col(kappa).cwiseProduct(pbf.alpha.col(kappa));
        pgbfbf.y.beta += (*gradBasisFunctionValues).y.col(kappa).cwiseProduct(pbf.beta.col(kappa));
        pgbfbf.z.alpha += (*gradBasisFunctionValues).z.col(kappa).cwiseProduct(pbf.alpha.col(kappa));
        pgbfbf.z.beta += (*gradBasisFunctionValues).z.col(kappa).cwiseProduct(pbf.beta.col(kappa));
        pgbfbf.x.alpha += basisFunctionValues.col(kappa).cwiseProduct(pgbf.x.alpha.col(kappa));
        pgbfbf.x.beta += basisFunctionValues.col(kappa).cwiseProduct(pgbf.x.beta.col(kappa));
        pgbfbf.y.alpha += basisFunctionValues.col(kappa).cwiseProduct(pgbf.y.alpha.col(kappa));
        pgbfbf.y.beta += basisFunctionValues.col(kappa).cwiseProduct(pgbf.y.beta.col(kappa));
        pgbfbf.z.alpha += basisFunctionValues.col(kappa).cwiseProduct(pgbf.z.alpha.col(kappa));
        pgbfbf.z.beta += basisFunctionValues.col(kappa).cwiseProduct(pgbf.z.beta.col(kappa));
      }
      //Norm based prescreening
      maxNorm = std::max(pbfbf.alpha.norm(),pbfbf.beta.norm());
      if (maxNorm < this->_system->getSettings().grid.blockAveThreshold) continue;
      //Contract LDA part
      scalarPart[iGuess].alpha += 2.0 * (*_d2FdRho2).aa.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
      scalarPart[iGuess].alpha += 2.0 * (*_d2FdRho2).ab.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
      scalarPart[iGuess].beta += 2.0 * (*_d2FdRho2).ab.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
      scalarPart[iGuess].beta += 2.0 * (*_d2FdRho2).bb.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
      if (_gga) {
        //Calculate ...\nabla \rho (\phi_\kappa \phi_\lambda)
        auto grhopbfbf = makeGradient<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> >(blockSize);
        grhopbfbf.x.aa = (*_densityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.x.ab = (*_densityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.x.ba = (*_densityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.x.bb = (*_densityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.y.aa = (*_densityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.y.ab = (*_densityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.y.ba = (*_densityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.y.bb = (*_densityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.z.aa = (*_densityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.z.ab = (*_densityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.z.ba = (*_densityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.z.bb = (*_densityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        //Calculate ...\nabla \rho \nabla(\phi_\kappa \phi_\lambda)
        DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd>  grhopgbfbf(blockSize);
        grhopgbfbf.aa.setZero();
        grhopgbfbf.ab.setZero();
        grhopgbfbf.ba.setZero();
        grhopgbfbf.bb.setZero();
        grhopgbfbf.aa += (*_densityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.alpha);
        grhopgbfbf.aa += (*_densityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.alpha);
        grhopgbfbf.aa += (*_densityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.alpha);
        grhopgbfbf.ab += (*_densityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.beta);
        grhopgbfbf.ab += (*_densityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.beta);
        grhopgbfbf.ab += (*_densityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.beta);
        grhopgbfbf.ba += (*_densityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.alpha);
        grhopgbfbf.ba += (*_densityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.alpha);
        grhopgbfbf.ba += (*_densityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.alpha);
        grhopgbfbf.bb += (*_densityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.beta);
        grhopgbfbf.bb += (*_densityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.beta);
        grhopgbfbf.bb += (*_densityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.beta);
        //Calculate ...\nabla rho \nabla \rho \nabla(\phi_\kappa \phi_\lambda)
        auto grhogrhopgbfbf_aaa = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_bbb = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_bab = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_aba = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_bba = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_aab = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_baa = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_abb = makeGradient<Eigen::VectorXd>(blockSize);
        grhogrhopgbfbf_aaa.x = (*_densityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_aaa.y = (*_densityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_aaa.z = (*_densityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_bbb.x = (*_densityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_bbb.y = (*_densityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_bbb.z = (*_densityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_bab.x = (*_densityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_bab.y = (*_densityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_bab.z = (*_densityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_aba.x = (*_densityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_aba.y = (*_densityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_aba.z = (*_densityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_bba.x = (*_densityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_bba.y = (*_densityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_bba.z = (*_densityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_aab.x = (*_densityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_aab.y = (*_densityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_aab.z = (*_densityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_baa.x = (*_densityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_baa.y = (*_densityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_baa.z = (*_densityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_abb.x = (*_densityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_abb.y = (*_densityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_abb.z = (*_densityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        //Contract with kernel
        scalarPart[iGuess].alpha += 4.0 * (*_d2FdRhodSigma).aaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        scalarPart[iGuess].alpha += 2.0 * (*_d2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        scalarPart[iGuess].alpha += 4.0 * (*_d2FdRhodSigma).abb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        scalarPart[iGuess].alpha += 2.0 * (*_d2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        scalarPart[iGuess].beta += 4.0 * (*_d2FdRhodSigma).bbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        scalarPart[iGuess].beta += 2.0 * (*_d2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        scalarPart[iGuess].beta += 4.0 * (*_d2FdRhodSigma).baa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        scalarPart[iGuess].beta += 2.0 * (*_d2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_d2FdRhodSigma).aaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.aa);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_d2FdRhodSigma).aaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.aa);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_d2FdRhodSigma).aaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.aa);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_d2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.ba);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_d2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.ba);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_d2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.ba);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_d2FdRhodSigma).baa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.ab);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_d2FdRhodSigma).baa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.ab);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_d2FdRhodSigma).baa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.ab);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_d2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.bb);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_d2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.bb);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_d2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.bb);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_d2FdRhodSigma).bbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.bb);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_d2FdRhodSigma).bbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.bb);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_d2FdRhodSigma).bbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.bb);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_d2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.ab);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_d2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.ab);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_d2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.ab);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_d2FdRhodSigma).abb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.ba);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_d2FdRhodSigma).abb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.ba);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_d2FdRhodSigma).abb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.ba);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_d2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.aa);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_d2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.aa);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_d2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.aa);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_dFdSigma).aa.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.alpha);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_dFdSigma).aa.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.alpha);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_dFdSigma).aa.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.alpha);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_dFdSigma).ab.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.beta);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_dFdSigma).ab.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.beta);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_dFdSigma).ab.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.beta);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_dFdSigma).bb.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.beta);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_dFdSigma).bb.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.beta);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_dFdSigma).bb.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.beta);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_dFdSigma).ab.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.alpha);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_dFdSigma).ab.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.alpha);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_dFdSigma).ab.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.alpha);
        (*gradientPart)[iGuess].x.alpha += 8.0 * (*_d2FdSigma2).aaaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.x);
        (*gradientPart)[iGuess].y.alpha += 8.0 * (*_d2FdSigma2).aaaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.y);
        (*gradientPart)[iGuess].z.alpha += 8.0 * (*_d2FdSigma2).aaaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.z);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.x);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.y);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.z);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.x);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.y);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.z);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.x);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.y);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.z);
        (*gradientPart)[iGuess].x.alpha += 8.0 * (*_d2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.x);
        (*gradientPart)[iGuess].y.alpha += 8.0 * (*_d2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.y);
        (*gradientPart)[iGuess].z.alpha += 8.0 * (*_d2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.z);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.x);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.y);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.z);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.x);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.y);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.z);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.x);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.y);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.z);
        (*gradientPart)[iGuess].x.beta += 8.0 * (*_d2FdSigma2).bbbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.x);
        (*gradientPart)[iGuess].y.beta += 8.0 * (*_d2FdSigma2).bbbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.y);
        (*gradientPart)[iGuess].z.beta += 8.0 * (*_d2FdSigma2).bbbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.z);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.x);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.y);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.z);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.x);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.y);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.z);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.x);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.y);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.z);
        (*gradientPart)[iGuess].x.beta += 8.0 * (*_d2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.x);
        (*gradientPart)[iGuess].y.beta += 8.0 * (*_d2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.y);
        (*gradientPart)[iGuess].z.beta += 8.0 * (*_d2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.z);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.x);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.y);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_d2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.z);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.x);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.y);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_d2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.z);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.x);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.y);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_d2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.z);
      }
      if (_nAddGGA) {
        //Calculate ...\nabla \rho (\phi_\kappa \phi_\lambda)
        auto grhopbfbf = makeGradient<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> >(blockSize);
        grhopbfbf.x.aa = (*_totalDensityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.x.ab = (*_totalDensityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.x.ba = (*_totalDensityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.x.bb = (*_totalDensityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.y.aa = (*_totalDensityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.y.ab = (*_totalDensityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.y.ba = (*_totalDensityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.y.bb = (*_totalDensityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.z.aa = (*_totalDensityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.z.ab = (*_totalDensityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        grhopbfbf.z.ba = (*_totalDensityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.alpha);
        grhopbfbf.z.bb = (*_totalDensityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pbfbf.beta);
        //Calculate ...\nabla \rho \nabla(\phi_\kappa \phi_\lambda)
        DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd>  grhopgbfbf(blockSize);
        grhopgbfbf.aa.setZero();
        grhopgbfbf.ab.setZero();
        grhopgbfbf.ba.setZero();
        grhopgbfbf.bb.setZero();
        grhopgbfbf.aa += (*_totalDensityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.alpha);
        grhopgbfbf.aa += (*_totalDensityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.alpha);
        grhopgbfbf.aa += (*_totalDensityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.alpha);
        grhopgbfbf.ab += (*_totalDensityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.beta);
        grhopgbfbf.ab += (*_totalDensityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.beta);
        grhopgbfbf.ab += (*_totalDensityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.beta);
        grhopgbfbf.ba += (*_totalDensityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.alpha);
        grhopgbfbf.ba += (*_totalDensityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.alpha);
        grhopgbfbf.ba += (*_totalDensityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.alpha);
        grhopgbfbf.bb += (*_totalDensityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.x.beta);
        grhopgbfbf.bb += (*_totalDensityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.y.beta);
        grhopgbfbf.bb += (*_totalDensityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(pgbfbf.z.beta);
        //Calculate ...\nabla rho \nabla \rho \nabla(\phi_\kappa \phi_\lambda)
        auto grhogrhopgbfbf_aaa = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_bbb = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_bab = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_aba = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_bba = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_aab = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_baa = makeGradient<Eigen::VectorXd>(blockSize);
        auto grhogrhopgbfbf_abb = makeGradient<Eigen::VectorXd>(blockSize);
        grhogrhopgbfbf_aaa.x = (*_totalDensityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_aaa.y = (*_totalDensityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_aaa.z = (*_totalDensityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_bbb.x = (*_totalDensityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_bbb.y = (*_totalDensityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_bbb.z = (*_totalDensityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_bab.x = (*_totalDensityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_bab.y = (*_totalDensityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_bab.z = (*_totalDensityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_aba.x = (*_totalDensityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_aba.y = (*_totalDensityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_aba.z = (*_totalDensityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_bba.x = (*_totalDensityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_bba.y = (*_totalDensityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_bba.z = (*_totalDensityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        grhogrhopgbfbf_aab.x = (*_totalDensityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_aab.y = (*_totalDensityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_aab.z = (*_totalDensityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        grhogrhopgbfbf_baa.x = (*_totalDensityGradient).x.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_baa.y = (*_totalDensityGradient).y.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_baa.z = (*_totalDensityGradient).z.beta.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        grhogrhopgbfbf_abb.x = (*_totalDensityGradient).x.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_abb.y = (*_totalDensityGradient).y.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        grhogrhopgbfbf_abb.z = (*_totalDensityGradient).z.alpha.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        //Contract with kernel
        scalarPart[iGuess].alpha += 4.0 * (*_totalD2FdRhodSigma).aaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        scalarPart[iGuess].alpha += 2.0 * (*_totalD2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        scalarPart[iGuess].alpha += 4.0 * (*_totalD2FdRhodSigma).abb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        scalarPart[iGuess].alpha += 2.0 * (*_totalD2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        scalarPart[iGuess].beta += 4.0 * (*_totalD2FdRhodSigma).bbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.bb);
        scalarPart[iGuess].beta += 2.0 * (*_totalD2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ab);
        scalarPart[iGuess].beta += 4.0 * (*_totalD2FdRhodSigma).baa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.aa);
        scalarPart[iGuess].beta += 2.0 * (*_totalD2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopgbfbf.ba);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_totalD2FdRhodSigma).aaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.aa);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_totalD2FdRhodSigma).aaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.aa);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_totalD2FdRhodSigma).aaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.aa);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_totalD2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.ba);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_totalD2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.ba);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_totalD2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.ba);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_totalD2FdRhodSigma).baa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.ab);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_totalD2FdRhodSigma).baa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.ab);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_totalD2FdRhodSigma).baa.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.ab);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_totalD2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.bb);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_totalD2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.bb);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_totalD2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.bb);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_totalD2FdRhodSigma).bbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.bb);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_totalD2FdRhodSigma).bbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.bb);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_totalD2FdRhodSigma).bbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.bb);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_totalD2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.ab);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_totalD2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.ab);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_totalD2FdRhodSigma).bab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.ab);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_totalD2FdRhodSigma).abb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.ba);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_totalD2FdRhodSigma).abb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.ba);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_totalD2FdRhodSigma).abb.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.ba);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_totalD2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.x.aa);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_totalD2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.y.aa);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_totalD2FdRhodSigma).aab.block(iGridStart,0,blockSize,1).cwiseProduct(grhopbfbf.z.aa);
        (*gradientPart)[iGuess].x.alpha += 8.0 * (*_totalD2FdSigma2).aaaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.x);
        (*gradientPart)[iGuess].y.alpha += 8.0 * (*_totalD2FdSigma2).aaaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.y);
        (*gradientPart)[iGuess].z.alpha += 8.0 * (*_totalD2FdSigma2).aaaa.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.z);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.x);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.y);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.z);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.x);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.y);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.z);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.x);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.y);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.z);
        (*gradientPart)[iGuess].x.alpha += 8.0 * (*_totalD2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.x);
        (*gradientPart)[iGuess].y.alpha += 8.0 * (*_totalD2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.y);
        (*gradientPart)[iGuess].z.alpha += 8.0 * (*_totalD2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.z);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.x);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.y);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.z);
        (*gradientPart)[iGuess].x.alpha += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.x);
        (*gradientPart)[iGuess].y.alpha += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.y);
        (*gradientPart)[iGuess].z.alpha += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.z);
        (*gradientPart)[iGuess].x.alpha += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.x);
        (*gradientPart)[iGuess].y.alpha += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.y);
        (*gradientPart)[iGuess].z.alpha += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.z);
        (*gradientPart)[iGuess].x.beta += 8.0 * (*_totalD2FdSigma2).bbbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.x);
        (*gradientPart)[iGuess].y.beta += 8.0 * (*_totalD2FdSigma2).bbbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.y);
        (*gradientPart)[iGuess].z.beta += 8.0 * (*_totalD2FdSigma2).bbbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bbb.z);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.x);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.y);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bab.z);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.x);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.y);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_abb.z);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.x);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.y);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aab.z);
        (*gradientPart)[iGuess].x.beta += 8.0 * (*_totalD2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.x);
        (*gradientPart)[iGuess].y.beta += 8.0 * (*_totalD2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.y);
        (*gradientPart)[iGuess].z.beta += 8.0 * (*_totalD2FdSigma2).aabb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_baa.z);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.x);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.y);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_totalD2FdSigma2).abbb.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_bba.z);
        (*gradientPart)[iGuess].x.beta += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.x);
        (*gradientPart)[iGuess].y.beta += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.y);
        (*gradientPart)[iGuess].z.beta += 4.0 * (*_totalD2FdSigma2).aaab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aaa.z);
        (*gradientPart)[iGuess].x.beta += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.x);
        (*gradientPart)[iGuess].y.beta += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.y);
        (*gradientPart)[iGuess].z.beta += 2.0 * (*_totalD2FdSigma2).abab.block(iGridStart,0,blockSize,1).cwiseProduct(grhogrhopgbfbf_aba.z);
      }
    }
  } else {
    assert(false);
  }
//  Timings::timeTaken("LRSCF: Contraction");
}




template<Options::SCF_MODES T>
Eigen::MatrixXd KernelSigmaVector<T>::calculateBlock(
    Eigen::MatrixXd& guess,
    std::vector<SpinPolarizedData<T,Eigen::MatrixXd> >& dens) {
  if (!_kernel) return Eigen::MatrixXd::Zero(guess.rows(),guess.cols());
  //Numerical integrator
  ScalarOperatorToMatrixAdder<T> numInt(
      _basisFunctionOnGridController,
      this->_system->getSettings().grid.blockAveThreshold);
  //Number of grid blocks
  const unsigned int nBlocks = _basisFunctionOnGridController->getNBlocks();
  //Number of basis functions
  const unsigned int nBasisFunc = this->_system->getBasisController()->getNBasisFunctions();

  //Build matrix for parallel computation
  unsigned int nThreads = 1.0;
  nThreads = (unsigned int)omp_get_max_threads();
  std::vector<std::vector<SPMatrix<T> > > threadMatrices(nThreads);
  for (unsigned int iThread=0; iThread< nThreads; ++iThread) {
    threadMatrices[iThread].resize(guess.cols(),SPMatrix<T>(nBasisFunc,nBasisFunc));
  }
#pragma omp parallel
  {
    unsigned int threadID = omp_get_thread_num();
    std::vector<SPMatrix<T> >& privateMatrix =  threadMatrices[threadID];
#pragma omp for schedule(static,1)
    //Loop over all blocks
    for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
      //calculate data for this block
      const auto& blockData = _basisFunctionOnGridController->getBlockOnGridData(iBlock);
      //Contract kernel with density
      std::vector<SpinPolarizedData<T,Eigen::VectorXd> > scalarPart(guess.cols());
      std::vector<Gradient<SpinPolarizedData<T,Eigen::VectorXd> > > gradientPart((_gga || _nAddGGA) ? guess.cols() : 0);
      if (!_gga && !_nAddGGA) {
        //LDA
        contractGridBlock(iBlock,blockData,dens,scalarPart);
      } else {
        //GGA
        contractGridBlock(iBlock,blockData,dens,scalarPart,&gradientPart);
      }
      //Integrate block and add to matrix
      for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
        if (!_gga && !_nAddGGA) {
          //LDA
          numInt.addBlock(iBlock,blockData,privateMatrix[iGuess],scalarPart[iGuess]);
        } else {
          //GGA
          numInt.addBlock(iBlock,blockData,privateMatrix[iGuess],scalarPart[iGuess],gradientPart[iGuess]);
        }
      }/* loop over all guess vectors */
    }/* loop over all blocks */
  }/* omp parallel */
  //Add thread matrices to Fock matrix
  auto& fock = threadMatrices[0];
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    SPMatrix<T>& f = fock[iGuess];
    for (unsigned int iThread=1; iThread< nThreads; ++iThread) {
      SPMatrix<T>& matrix = threadMatrices[iThread][iGuess];
      for_spin(f,matrix){
        f_spin += matrix_spin;
      };
    }
  }

  //Transform to MO basis and obtain sigma vectors
  Eigen::MatrixXd sigma(guess.rows(),guess.cols());
  for (unsigned int iGuess = 0; iGuess < guess.cols(); ++iGuess) {
    sigma.col(iGuess) = this->ao2mo(fock[iGuess]);
  }
  return sigma;
}

template class KernelSigmaVector<Options::SCF_MODES::RESTRICTED> ;
template class KernelSigmaVector<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
