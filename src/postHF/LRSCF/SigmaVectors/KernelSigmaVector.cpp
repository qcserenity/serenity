/**
 * @file KernelSigmaVector.cpp
 *
 * @date Dec 07, 2018
 * @author Johannes Tölle
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
#include "postHF/LRSCF/SigmaVectors/KernelSigmaVector.h"

/* Include Serenity Internal Headers */
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "misc/Timing.h"
#include "misc/HelperFunctions.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unistd.h>
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
KernelSigmaVector<SCFMode>::KernelSigmaVector(
    std::vector<std::shared_ptr<LRSCFController<SCFMode> > > lrscf,
    std::vector<Eigen::MatrixXd> b,
    const double densityScreeningThreshold,
    std::shared_ptr<Kernel<SCFMode> > kernel):
        SigmaVector<SCFMode>(lrscf,b,densityScreeningThreshold),
        _kernel(kernel)
    {
    }

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode> > > > KernelSigmaVector<SCFMode>::calcF(
    unsigned int I,
    unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode> > > > densityMatrices) {

  Timings::takeTime("LRSCF - Fock-like matrix: fxc");
  
  //This is needed for coupling pattern where the same system is with coupled 
  //and uncoupled vectors. The Kernel has only system information while I and J are
  //related to the lrscf Controller, but the systems behind I and J could be indentical
  for (unsigned int ilrscf = 0; ilrscf < this->_lrscf.size(); ilrscf++){
    if(this->_lrscf[ilrscf]->getSys() == this->_lrscf[I]->getSys()) {
      I = ilrscf;
      break;
    }
  }
  for (unsigned int ilrscf = 0; ilrscf < this->_lrscf.size(); ilrscf++){
    if(this->_lrscf[ilrscf]->getSys() == this->_lrscf[J]->getSys()) {
      J = ilrscf;
      break;
    }
  }
  
  //Set dimensions for Fock like matrices
  //Final dimensions are the dimensions of subsystem I
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode> > > > fock(
      new std::vector<std::vector<MatrixInBasis<SCFMode>>>(this->_nSet));
  for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
  }

  //Return zero if kernel is nullptr
  if (!_kernel) return fock;

  //GGA or not
  bool gga = _kernel->isGGA();
  //Dervative Level required depends on the adiabatic functional
  unsigned int derivativeLevel = (gga) ? 1 : 0;
  //BasisfunctionOnGridController:
  std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerI = nullptr;
  std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerJ = nullptr;

  basisFunctionOnGridControllerI = BasisFunctionOnGridControllerFactory::produce(
                      this->_kernel->getBlocksize(I),
                      this->_kernel->getbasFuncRadialThreshold(I),
                      derivativeLevel,
                      this->_lrscf[I]->getBasisController(),
                      _kernel->getGridController());
  //Two controller needed if the Sigma Vector of two different subsystems is calculated
  if(I != J) {
    basisFunctionOnGridControllerJ = BasisFunctionOnGridControllerFactory::produce(
        this->_kernel->getBlocksize(J),
        this->_kernel->getbasFuncRadialThreshold(J),
        derivativeLevel,
        this->_lrscf[J]->getBasisController(),
        _kernel->getGridController());
  }
  //Initialize scalar and gradient part for integration
  std::vector<std::vector<GridPotential<SCFMode> > > scalarContr(this -> _nSet);
  for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < this ->_nGuess; ++iGuess) {
      scalarContr[iSet].emplace_back(_kernel->getGridController());
    }
  }
  std::vector<std::vector<Gradient<GridPotential<SCFMode> > > > gradientContr(this -> _nSet);
  if (gga){
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        gradientContr[iSet].emplace_back(makeGradient<GridPotential<SCFMode> >(_kernel->getGridController()));
      }
    }
  }

  //Contract kernel with density and basis functions and adds nummerical integration weights
  Timings::takeTime("LRSCF -    Kernel Contraction");
  if (I == J){
    contractKernel((*densityMatrices), basisFunctionOnGridControllerI ,scalarContr,gradientContr, I, J);
  }else{
    contractKernel((*densityMatrices), basisFunctionOnGridControllerJ ,scalarContr,gradientContr, I, J);
  }
  Timings::timeTaken("LRSCF -    Kernel Contraction");

  
  std::vector<std::vector<std::vector<MatrixInBasis<SCFMode> > > *> f;
#ifdef _OPENMP
  const unsigned int nThreads = omp_get_max_threads();
  f.push_back(&(*fock));
  for (unsigned int iThread=1; iThread < nThreads; ++iThread) {
    f.push_back(new std::vector<std::vector<MatrixInBasis<SCFMode> > >);
    (*f[iThread]).resize(this->_nSet);
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
    }
  }
  Eigen::setNbThreads(1);
#else
  f.push_back(&fock);
#endif

  Timings::takeTime("LRSCF - Numerical Integration");
  //Nummerical Integration
  if(gga){
    numIntSigma(f,basisFunctionOnGridControllerI,scalarContr,gradientContr,I);
  }else{
    numIntSigma(f,basisFunctionOnGridControllerI,scalarContr,I);
  }
  Timings::timeTaken("LRSCF - Numerical Integration");

//Sum over all threads
#ifdef _OPENMP
  Eigen::setNbThreads(0);
  for (unsigned int iThread=1; iThread<nThreads; ++iThread) {
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess=0; iGuess< this->_nGuess; ++iGuess) {
        (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
      }
    }
    delete f[iThread];
  }
#endif

  Timings::timeTaken("LRSCF - Fock-like matrix: fxc");
  
  return fock;
}

template<Options::SCF_MODES SCFMode>
void KernelSigmaVector<SCFMode>::contractKernel(
    std::vector<std::vector<MatrixInBasis<SCFMode> > >& dens,
    std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
    std::vector<std::vector<GridPotential<SCFMode> > >& scalarPart,
    std::vector<std::vector<Gradient<GridPotential<SCFMode> > > >& gradientPart,
    unsigned int I,
    unsigned int J) {

  //Number of grid blocks
  const unsigned int nBlocks = basisFunctionOnGridController->getNBlocks();
  //Number of basis functions
  const unsigned int nBasisFunc = basisFunctionOnGridController->getNBasisFunctions();
  bool gga = _kernel->isGGA();
  //Set threads for Eigen in parrallel region
  Eigen::setNbThreads(1);

#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
      //calculate data for this block
      const auto& blockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
      //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
      Eigen::MatrixXd& basisFunctionValues = blockData->functionValues;
      //gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
      const auto& gradBasisFunctionValues = blockData->derivativeValues;
      //number of grid points in this block
      const unsigned int blockSize = blockData->functionValues.rows();
      //Get first index of this Block
      const unsigned int iGridStart = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
      //weights for numerical integration
      const Eigen::VectorXd& weights = basisFunctionOnGridController->getGridController()->getWeights();
      // Negligible vector 1: uniportant 0: important
      const Eigen::VectorXi& negligible = blockData->negligible;   
      //calculate projector for non neglible basis functions
      const auto projector = constructProjectionMatrix(negligible);
      //Project into smaller basis
      const Eigen::MatrixXd basisFuncProj = basisFunctionValues * projector;
      Eigen::MatrixXd gradBasisFuncProj_x;
      Eigen::MatrixXd gradBasisFuncProj_y;
      Eigen::MatrixXd gradBasisFuncProj_z;
      if (gga) {
        gradBasisFuncProj_x = gradBasisFunctionValues->x * projector;
        gradBasisFuncProj_y = gradBasisFunctionValues->y * projector;
        gradBasisFuncProj_z = gradBasisFunctionValues->z * projector;
      }
      
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {

          SpinPolarizedData<SCFMode,Eigen::VectorXd>  pbb(blockSize);
          Gradient<SpinPolarizedData<SCFMode,Eigen::VectorXd> >  pnbb =  makeGradient<SpinPolarizedData<SCFMode,Eigen::VectorXd> >((gga) ? blockSize : 0);
          //contract basis functions with density of guess vector iGuess
          auto& p = dens[iSet][iGuess];
          SpinPolarizedData<SCFMode,Eigen::MatrixXd> pb(blockSize,nBasisFunc);
        
          for_spin(p,pb,pbb) {
            pb_spin.setZero();
            //This needs to be done twice because the density matrix is not symmetric
            pb_spin = basisFuncProj*p_spin.transpose();
            pb_spin += basisFuncProj*p_spin;
            pbb_spin = 0.5 * basisFuncProj.cwiseProduct(pb_spin).rowwise().sum();
            //Now contract with weights
            pbb_spin = pbb_spin.cwiseProduct(weights.segment(iGridStart,blockSize));
          };

          if (gga) {
            auto& pnbbx = pnbb.x;
            auto& pnbby = pnbb.y;
            auto& pnbbz = pnbb.z;
            for_spin(pb,pnbbx,pnbby,pnbbz) {
              pnbbx_spin = gradBasisFuncProj_x.cwiseProduct(pb_spin).rowwise().sum();
              pnbbx_spin = pnbbx_spin.cwiseProduct(weights.segment(iGridStart,blockSize));
              pnbby_spin = gradBasisFuncProj_y.cwiseProduct(pb_spin).rowwise().sum();
              pnbby_spin = pnbby_spin.cwiseProduct(weights.segment(iGridStart,blockSize));
              pnbbz_spin = gradBasisFuncProj_z.cwiseProduct(pb_spin).rowwise().sum();
              pnbbz_spin = pnbbz_spin.cwiseProduct(weights.segment(iGridStart,blockSize));
            };
          }
          //Contract with kernel (Note that restricted and unrestricted cannot be handled in
          //a single function since there is not yet a for_spin looper function for DoublySpinPolarized data)
          contractBlock(iGridStart,blockSize,pbb,pnbb,scalarPart[iSet][iGuess],gradientPart[iSet][iGuess], I, J);
        }/* Loop over iGuess */
      }/* Loop over iSet */
    }/* Loop over iBlock */
  }
  //reset threads for Eigen
  Eigen::setNbThreads(0);
}

template<Options::SCF_MODES SCFMode>
void KernelSigmaVector<SCFMode>::numIntSigma(
      std::vector<std::vector<std::vector<MatrixInBasis<SCFMode> > > *>& focklikeMatrix,
      std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
      std::vector<std::vector<GridPotential<SCFMode> > >& scalarPart,
      unsigned int I){

  //Number of grid blocks
  const unsigned int nBlocks = basisFunctionOnGridController->getNBlocks();
  //Block average threshold
  double blockAveThreshold = _kernel->getblockAveThreshold(I);
  //Set threads for Eigen in parrallel region
  Eigen::setNbThreads(1);
#pragma omp parallel
  {
  //ThreadID
  const unsigned int threadID = omp_get_thread_num();
#pragma omp for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
      //calculate data for this block
      const auto& blockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
      //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
      Eigen::MatrixXd& basisFunctionValues = blockData->functionValues;
      //number of grid points in this block
      const unsigned int blockSize = blockData->functionValues.rows();
      //Get first index of this Block
      const unsigned int iGridStart = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
      // Negligible vector 1: uniportant 0: important
      const Eigen::VectorXi& negligible = blockData->negligible;   
      //calculate projector for non neglible basis functions
      const auto projector = constructProjectionMatrix(negligible);
      //Project into smaller basis
      const Eigen::MatrixXd basisFuncProj = basisFunctionValues * projector;
      std::cout<<"here"<<std::endl;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          
          // Initialize data
          const GridPotential<SCFMode>& scalarpart = scalarPart[iSet][iGuess];
          MatrixInBasis<SCFMode>& f = (*focklikeMatrix[threadID])[iSet][iGuess];

          for_spin(f,scalarpart) {
            const Eigen::VectorXd scalar = scalarpart_spin.segment(iGridStart, blockSize);
            //Test of significance
            double average = scalar.cwiseAbs().sum();
            average /= blockSize;
            if (average < blockAveThreshold) return;
            //Evaluate quantities
            Eigen::MatrixXd temp = (basisFuncProj.array().colwise() * scalar.array()).matrix().transpose() * basisFuncProj;
      
            f_spin += projector * temp * projector.transpose();
          };
        }/*End Guess*/
      }/*End Set*/
    }/*End Block*/
  }/*End Parallel*/
  //reset threads for Eigen
  Eigen::setNbThreads(0);
  return;
}



template<Options::SCF_MODES SCFMode>
void KernelSigmaVector<SCFMode>::numIntSigma(
      std::vector<std::vector<std::vector<MatrixInBasis<SCFMode> > > *>& focklikeMatrix,
      std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
      std::vector<std::vector<GridPotential<SCFMode> >  >& scalarPart,
      std::vector<std::vector<Gradient<GridPotential<SCFMode> > > >& gradientPart,
      unsigned int I){

  //Number of grid blocks
  const unsigned int nBlocks = basisFunctionOnGridController->getNBlocks();
  //Block average threshold
  double blockAveThreshold = _kernel->getblockAveThreshold(I);
  //Set threads for Eigen in parrallel region
  Eigen::setNbThreads(1);
#pragma omp parallel
  {
  //ThreadID
  const unsigned int threadID = omp_get_thread_num();
#pragma omp for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
      //calculate data for this block
      const auto& blockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
      //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
      const Eigen::MatrixXd& basisFunctionValues = blockData->functionValues;
      //gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
      const auto& gradBasisFunctionValues = blockData->derivativeValues;
      //number of grid points in this block
      const unsigned int blockSize = blockData->functionValues.rows();
      //Get first index of this Block
      const unsigned int iGridStart = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
      // Negligible vector 1: uniportant 0: important
      const Eigen::VectorXi& negligible = blockData->negligible;   
      //calculate projector for non neglible basis functions
      const auto projector = constructProjectionMatrix(negligible);
      //Project into smaller basis
      const Eigen::MatrixXd basisFuncProj = basisFunctionValues * projector;
      const Eigen::MatrixXd gradBasisFuncProj_x = gradBasisFunctionValues->x * projector;
      const Eigen::MatrixXd gradBasisFuncProj_y = gradBasisFunctionValues->y * projector;
      const Eigen::MatrixXd gradBasisFuncProj_z = gradBasisFunctionValues->z * projector;

      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          
          // Initialize data
          const GridPotential<SCFMode>& scalarpart = scalarPart[iSet][iGuess];
          const GridData<SCFMode>& gradientPartX = gradientPart[iSet][iGuess].x;
          const GridData<SCFMode>& gradientPartY = gradientPart[iSet][iGuess].y;
          const GridData<SCFMode>& gradientPartZ = gradientPart[iSet][iGuess].z;
          
          MatrixInBasis<SCFMode>& f = (*focklikeMatrix[threadID])[iSet][iGuess];

          for_spin(f,scalarpart,gradientPartX,gradientPartY,gradientPartZ) {
            //Test of significance
            double average = scalarpart_spin.segment(iGridStart,blockSize).cwiseAbs().sum();
            average += gradientPartX_spin.segment(iGridStart,blockSize).cwiseAbs().sum();
            average += gradientPartY_spin.segment(iGridStart,blockSize).cwiseAbs().sum();
            average += gradientPartZ_spin.segment(iGridStart,blockSize).cwiseAbs().sum();
            average /= blockSize;
            if (average < blockAveThreshold) return;
          
            //Evaluate quantities
            Eigen::MatrixXd temp = basisFuncProj.array().colwise() * scalarpart_spin.segment(iGridStart,blockSize).array();
            Eigen::MatrixXd temp2 = temp.transpose() * basisFuncProj;
            temp = gradBasisFuncProj_x.array().colwise() * gradientPartX_spin.segment(iGridStart,blockSize).array()
                 + gradBasisFuncProj_y.array().colwise() * gradientPartY_spin.segment(iGridStart,blockSize).array()
                 + gradBasisFuncProj_z.array().colwise() * gradientPartZ_spin.segment(iGridStart,blockSize).array();
            temp2 += basisFuncProj.transpose() * temp + temp.transpose() * basisFuncProj; 
            //Project back into bigger basis         
            f_spin += projector * temp2 * projector.transpose();
          };
        }/*End Guess*/
      }/*End Set*/
    }/*End Block*/
  }/*End Parallel*/
  //reset threads for Eigen
  Eigen::setNbThreads(0);
  return;
}


template<>
void KernelSigmaVector<Options::SCF_MODES::RESTRICTED>::contractBlock(
    const unsigned int iGridStart,
    const unsigned int blockSize,
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::VectorXd>& pbb,
    Gradient<SpinPolarizedData<Options::SCF_MODES::RESTRICTED,Eigen::VectorXd> >&  pnbb,
    GridPotential<Options::SCF_MODES::RESTRICTED>& scalarPart,
    Gradient<GridPotential<Options::SCF_MODES::RESTRICTED> >& gradientPart,
    unsigned int I,
    unsigned int J) {
  auto ppptr = _kernel->getPP(I,J,blockSize,iGridStart);
  auto& pp = (*ppptr);
  scalarPart.segment(iGridStart, blockSize) += pp.cwiseProduct(pbb);
  //scalarPart.segment(iGridStart, blockSize) += pp.segment(iGridStart,blockSize).cwiseProduct(pbb);
  if (_kernel->isGGA()) {
    auto pgptr = _kernel->getPG(I,J,blockSize,iGridStart);
    auto ggptr = _kernel->getGG(I,J,blockSize,iGridStart);
    auto& pg = (*pgptr);
    auto& gg = (*ggptr);
    scalarPart.segment(iGridStart, blockSize) += pg.x.cwiseProduct(pnbb.x);
    scalarPart.segment(iGridStart, blockSize) += pg.y.cwiseProduct(pnbb.y);
    scalarPart.segment(iGridStart, blockSize) += pg.z.cwiseProduct(pnbb.z);
    gradientPart.x.segment(iGridStart, blockSize) += pg.x.cwiseProduct(pbb);
    gradientPart.y.segment(iGridStart, blockSize) += pg.y.cwiseProduct(pbb);
    gradientPart.z.segment(iGridStart, blockSize) += pg.z.cwiseProduct(pbb);
    gradientPart.x.segment(iGridStart, blockSize) += gg.xx.cwiseProduct(pnbb.x);
    gradientPart.x.segment(iGridStart, blockSize) += gg.xy.cwiseProduct(pnbb.y);
    gradientPart.x.segment(iGridStart, blockSize) += gg.xz.cwiseProduct(pnbb.z);
    gradientPart.y.segment(iGridStart, blockSize) += gg.xy.cwiseProduct(pnbb.x);
    gradientPart.y.segment(iGridStart, blockSize) += gg.yy.cwiseProduct(pnbb.y);
    gradientPart.y.segment(iGridStart, blockSize) += gg.yz.cwiseProduct(pnbb.z);
    gradientPart.z.segment(iGridStart, blockSize) += gg.xz.cwiseProduct(pnbb.x);
    gradientPart.z.segment(iGridStart, blockSize) += gg.yz.cwiseProduct(pnbb.y);
    gradientPart.z.segment(iGridStart, blockSize) += gg.zz.cwiseProduct(pnbb.z);
  }
}

template<>
void KernelSigmaVector<Options::SCF_MODES::UNRESTRICTED>::contractBlock(
    const unsigned int iGridStart,
    const unsigned int blockSize,
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd>& pbb,
    Gradient<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,Eigen::VectorXd> >&  pnbb,
    GridPotential<Options::SCF_MODES::UNRESTRICTED>& scalarPart,
    Gradient<GridPotential<Options::SCF_MODES::UNRESTRICTED> >& gradientPart,
    unsigned int I,
    unsigned int J) {
  auto ppptr = _kernel->getPP(I,J,blockSize,iGridStart);
  auto& pp = (*ppptr);
  scalarPart.alpha.segment(iGridStart, blockSize) += pp.aa.cwiseProduct(pbb.alpha);
  scalarPart.alpha.segment(iGridStart, blockSize) += pp.ab.cwiseProduct(pbb.beta);
  scalarPart.beta.segment(iGridStart, blockSize) += pp.bb.cwiseProduct(pbb.beta);
  scalarPart.beta.segment(iGridStart, blockSize) += pp.ab.cwiseProduct(pbb.alpha);
  if (_kernel->isGGA()) {
    auto pgptr = _kernel->getPG(I,J,blockSize,iGridStart);
    auto ggptr = _kernel->getGG(I,J,blockSize,iGridStart);
    auto& pg = (*pgptr);
    auto& gg = (*ggptr);
    scalarPart.alpha.segment(iGridStart, blockSize) += pg.x.aa.cwiseProduct(pnbb.x.alpha);
    scalarPart.alpha.segment(iGridStart, blockSize) += pg.y.aa.cwiseProduct(pnbb.y.alpha);
    scalarPart.alpha.segment(iGridStart, blockSize) += pg.z.aa.cwiseProduct(pnbb.z.alpha);
    scalarPart.alpha.segment(iGridStart, blockSize) += pg.x.ab.cwiseProduct(pnbb.x.beta);
    scalarPart.alpha.segment(iGridStart, blockSize) += pg.y.ab.cwiseProduct(pnbb.y.beta);
    scalarPart.alpha.segment(iGridStart, blockSize) += pg.z.ab.cwiseProduct(pnbb.z.beta);
    scalarPart.beta.segment(iGridStart, blockSize) += pg.x.ba.cwiseProduct(pnbb.x.alpha);
    scalarPart.beta.segment(iGridStart, blockSize) += pg.y.ba.cwiseProduct(pnbb.y.alpha);
    scalarPart.beta.segment(iGridStart, blockSize) += pg.z.ba.cwiseProduct(pnbb.z.alpha);
    scalarPart.beta.segment(iGridStart, blockSize) += pg.x.bb.cwiseProduct(pnbb.x.beta);
    scalarPart.beta.segment(iGridStart, blockSize) += pg.y.bb.cwiseProduct(pnbb.y.beta);
    scalarPart.beta.segment(iGridStart, blockSize) += pg.z.bb.cwiseProduct(pnbb.z.beta);
    gradientPart.x.alpha.segment(iGridStart, blockSize) += pg.x.aa.cwiseProduct(pbb.alpha);
    gradientPart.y.alpha.segment(iGridStart, blockSize) += pg.y.aa.cwiseProduct(pbb.alpha);
    gradientPart.z.alpha.segment(iGridStart, blockSize) += pg.z.aa.cwiseProduct(pbb.alpha);
    gradientPart.x.alpha.segment(iGridStart, blockSize) += pg.x.ba.cwiseProduct(pbb.beta);
    gradientPart.y.alpha.segment(iGridStart, blockSize) += pg.y.ba.cwiseProduct(pbb.beta);
    gradientPart.z.alpha.segment(iGridStart, blockSize) += pg.z.ba.cwiseProduct(pbb.beta);
    gradientPart.x.beta.segment(iGridStart, blockSize) += pg.x.ab.cwiseProduct(pbb.alpha);
    gradientPart.y.beta.segment(iGridStart, blockSize) += pg.y.ab.cwiseProduct(pbb.alpha);
    gradientPart.z.beta.segment(iGridStart, blockSize) += pg.z.ab.cwiseProduct(pbb.alpha);
    gradientPart.x.beta.segment(iGridStart, blockSize) += pg.x.bb.cwiseProduct(pbb.beta);
    gradientPart.y.beta.segment(iGridStart, blockSize) += pg.y.bb.cwiseProduct(pbb.beta);
    gradientPart.z.beta.segment(iGridStart, blockSize) += pg.z.bb.cwiseProduct(pbb.beta);
    gradientPart.x.alpha.segment(iGridStart, blockSize) += gg.xx.aa.cwiseProduct(pnbb.x.alpha);
    gradientPart.x.alpha.segment(iGridStart, blockSize) += gg.xy.aa.cwiseProduct(pnbb.y.alpha);
    gradientPart.x.alpha.segment(iGridStart, blockSize) += gg.xz.aa.cwiseProduct(pnbb.z.alpha);
    gradientPart.x.alpha.segment(iGridStart, blockSize) += gg.xx.ab.cwiseProduct(pnbb.x.beta);
    gradientPart.x.alpha.segment(iGridStart, blockSize) += gg.xy.ab.cwiseProduct(pnbb.y.beta);
    gradientPart.x.alpha.segment(iGridStart, blockSize) += gg.xz.ab.cwiseProduct(pnbb.z.beta);
    gradientPart.x.beta.segment(iGridStart, blockSize) += gg.xx.bb.cwiseProduct(pnbb.x.beta);
    gradientPart.x.beta.segment(iGridStart, blockSize) += gg.xy.bb.cwiseProduct(pnbb.y.beta);
    gradientPart.x.beta.segment(iGridStart, blockSize) += gg.xz.bb.cwiseProduct(pnbb.z.beta);
    gradientPart.x.beta.segment(iGridStart, blockSize) += gg.xx.ba.cwiseProduct(pnbb.x.alpha);
    gradientPart.x.beta.segment(iGridStart, blockSize) += gg.xy.ba.cwiseProduct(pnbb.y.alpha);
    gradientPart.x.beta.segment(iGridStart, blockSize) += gg.xz.ba.cwiseProduct(pnbb.z.alpha);
    gradientPart.y.alpha.segment(iGridStart, blockSize) += gg.xy.aa.cwiseProduct(pnbb.x.alpha);
    gradientPart.y.alpha.segment(iGridStart, blockSize) += gg.yy.aa.cwiseProduct(pnbb.y.alpha);
    gradientPart.y.alpha.segment(iGridStart, blockSize) += gg.yz.aa.cwiseProduct(pnbb.z.alpha);
    gradientPart.y.alpha.segment(iGridStart, blockSize) += gg.xy.ab.cwiseProduct(pnbb.x.beta);
    gradientPart.y.alpha.segment(iGridStart, blockSize) += gg.yy.ab.cwiseProduct(pnbb.y.beta);
    gradientPart.y.alpha.segment(iGridStart, blockSize) += gg.yz.ab.cwiseProduct(pnbb.z.beta);
    gradientPart.y.beta.segment(iGridStart, blockSize) += gg.xy.bb.cwiseProduct(pnbb.x.beta);
    gradientPart.y.beta.segment(iGridStart, blockSize) += gg.yy.bb.cwiseProduct(pnbb.y.beta);
    gradientPart.y.beta.segment(iGridStart, blockSize) += gg.yz.bb.cwiseProduct(pnbb.z.beta);
    gradientPart.y.beta.segment(iGridStart, blockSize) += gg.xy.ba.cwiseProduct(pnbb.x.alpha);
    gradientPart.y.beta.segment(iGridStart, blockSize) += gg.yy.ba.cwiseProduct(pnbb.y.alpha);
    gradientPart.y.beta.segment(iGridStart, blockSize) += gg.yz.ba.cwiseProduct(pnbb.z.alpha);
    gradientPart.z.alpha.segment(iGridStart, blockSize) += gg.xz.aa.cwiseProduct(pnbb.x.alpha);
    gradientPart.z.alpha.segment(iGridStart, blockSize) += gg.yz.aa.cwiseProduct(pnbb.y.alpha);
    gradientPart.z.alpha.segment(iGridStart, blockSize) += gg.zz.aa.cwiseProduct(pnbb.z.alpha);
    gradientPart.z.alpha.segment(iGridStart, blockSize) += gg.xz.ab.cwiseProduct(pnbb.x.beta);
    gradientPart.z.alpha.segment(iGridStart, blockSize) += gg.yz.ab.cwiseProduct(pnbb.y.beta);
    gradientPart.z.alpha.segment(iGridStart, blockSize) += gg.zz.ab.cwiseProduct(pnbb.z.beta);
    gradientPart.z.beta.segment(iGridStart, blockSize) += gg.xz.bb.cwiseProduct(pnbb.x.beta);
    gradientPart.z.beta.segment(iGridStart, blockSize) += gg.yz.bb.cwiseProduct(pnbb.y.beta);
    gradientPart.z.beta.segment(iGridStart, blockSize) += gg.zz.bb.cwiseProduct(pnbb.z.beta);
    gradientPart.z.beta.segment(iGridStart, blockSize) += gg.xz.ba.cwiseProduct(pnbb.x.alpha);
    gradientPart.z.beta.segment(iGridStart, blockSize) += gg.yz.ba.cwiseProduct(pnbb.y.alpha);
    gradientPart.z.beta.segment(iGridStart, blockSize) += gg.zz.ba.cwiseProduct(pnbb.z.alpha);
  }
}

template class KernelSigmaVector<Options::SCF_MODES::RESTRICTED>;
template class KernelSigmaVector<Options::SCF_MODES::UNRESTRICTED>;
}