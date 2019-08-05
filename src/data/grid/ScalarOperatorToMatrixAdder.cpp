/**
 * @file   ScalarOperatorToMatrixAdder.cpp
 *
 * @date   Mar 22, 2014, Apr 20, 2017
 * @author Thomas Dresselhaus, M. Boeckers
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
#include "data/grid/ScalarOperatorToMatrixAdder.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "math/FloatMaths.h"
#include "data/matrices/FockMatrix.h"
#include "grid/GridController.h"
#include "misc/HelperFunctions.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <cassert>


namespace Serenity {
using namespace std;

template<Options::SCF_MODES SCFMode>ScalarOperatorToMatrixAdder<SCFMode>::ScalarOperatorToMatrixAdder(
      shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
      const double blockAveThreshold) :
          _basisFunctionOnGridControllerA(basisFunctionOnGridController),
          _basisFunctionOnGridControllerB(basisFunctionOnGridController),
          _blockAveThreshold(blockAveThreshold) {
  assert(_basisFunctionOnGridControllerA);
#ifdef _OPENMP
  omp_init_lock(&_lock);
#endif
}

template<Options::SCF_MODES SCFMode>ScalarOperatorToMatrixAdder<SCFMode>::ScalarOperatorToMatrixAdder(
    std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerA,
    std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerB,
    double blockAveThreshold) :
        _basisFunctionOnGridControllerA(basisFunctionOnGridControllerA),
        _basisFunctionOnGridControllerB(basisFunctionOnGridControllerB),
        _blockAveThreshold(blockAveThreshold) {
  assert(_basisFunctionOnGridControllerA);
  assert(_basisFunctionOnGridControllerB);
  assert(isDefinedOnSameGrid(*_basisFunctionOnGridControllerA,*_basisFunctionOnGridControllerB) &&
      "Both basis sets have to be evaluated on the same grid!");
#ifdef _OPENMP
  omp_init_lock(&_lock);
#endif
}

template<Options::SCF_MODES SCFMode>void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToMatrix(
    SPMatrix<SCFMode>& matrix,
    const GridPotential<SCFMode>& scalarOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  //checks
  assert(scalarOperator.isValid());
  assert(isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA));
//  assert(_basisFunctionOnGridControllerA && _basisFunctionOnGridControllerB);
  //Get number of grid blocks
  const unsigned int nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  //Build matrix for parallel computation
  unsigned int nThreads = (unsigned int)omp_get_max_threads();
  std::vector<std::unique_ptr<SPMatrix<SCFMode> > > threadMatrices(nThreads);

#pragma omp parallel
  {
    unsigned int threadID = omp_get_thread_num();
    threadMatrices[threadID].reset(new SPMatrix<SCFMode>(
        _basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions(),
        _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions()));

    //Loop over blocks of grid points
#pragma omp for schedule(static,1)
    for (unsigned int blockNumber=0; blockNumber < nBlocks; ++blockNumber) {
      //Data for this block
      auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(blockNumber);
      auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(blockNumber);
      //Integrate and add to Fock matrix for this thread
//      if (isDefinedInSameBasis(*_basisFunctionOnGridControllerA,*_basisFunctionOnGridControllerB)) {
//        addBlock(blockNumber,blockDataA,*threadMatrices[threadID],scalarOperator);
//      } else {
        addBlock(blockNumber,blockDataA,blockDataB,*threadMatrices[threadID],scalarOperator);
//      }
    }/* loop over all blocks of grid points */
  } // end pragma omp parallel
  for (auto& mat : threadMatrices){
    if (mat != nullptr){
      matrix += *mat;
    }
  }
  Timings::timeTaken("Tech. -    Grid to Matrix Int.");
  return;
}

template<Options::SCF_MODES SCFMode>void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToMatrix(
    // TODO maybe a different type?
    SPMatrix<SCFMode>& matrix,
    const GridPotential<SCFMode>& scalarOperator,
    const Gradient<GridPotential<SCFMode> >& gradientOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  //checks
  assert(scalarOperator.isValid());
  assert(isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA));
  assert (_basisFunctionOnGridControllerA->getHighestDerivative() >= 1);
  assert (_basisFunctionOnGridControllerB->getHighestDerivative() >= 1);
  for (const auto& component : gradientOperator) {
    if (!component.isValid())
      throw SerenityError("ScalarOperatorToMatrixAdder: Component is invalid.");
    if (!isDefinedOnSameGrid(component, *_basisFunctionOnGridControllerA))
      throw SerenityError("ScalarOperatorToMatrixAdder: Components are not defined on the same grid.");
  }
  //Get number of grid blocks
  const unsigned int nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  //Build matrix for parallel computation
  unsigned int nThreads = (unsigned int)omp_get_max_threads();
  std::vector<std::unique_ptr<SPMatrix<SCFMode> > > threadMatrices(nThreads);

#pragma omp parallel
  {
    unsigned int threadID = omp_get_thread_num();
    threadMatrices[threadID].reset(new SPMatrix<SCFMode>(
        _basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions(),
        _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions()));

    //Loop over blocks of grid points
#pragma omp for schedule(static,1)
    for (unsigned int blockNumber=0; blockNumber < nBlocks; ++blockNumber) {
      //Data for this block
      auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(blockNumber);
      auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(blockNumber);
      //Integrate and add to Fock matrix for this thread
      addBlock(blockNumber,blockDataA,blockDataB,*threadMatrices[threadID],scalarOperator,gradientOperator);
    }/* loop over all blocks of grid points */
  } // end pragma omp parallel
  for (auto& mat : threadMatrices){
    if (mat != nullptr){
      matrix += *mat;
    }
  }
  Timings::timeTaken("Tech. -    Grid to Matrix Int.");
  return;
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addBlock(
    unsigned int iBlock,
    std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataA,
    std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataB,
    SPMatrix<SCFMode>& m_AB,
    const GridPotential<SCFMode>& scalarPart) {
  bool useSym = isDefinedInSameBasis(*_basisFunctionOnGridControllerA,*_basisFunctionOnGridControllerB);
  //Get weights
  const auto& weights = _basisFunctionOnGridControllerA->getGridController()->getWeights();
  //Number of basis functions
  const unsigned int nBasisFuncA = _basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions();
  const unsigned int nBasisFuncB = _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions();
  // Sanity checks
  assert(m_AB.rows() == nBasisFuncA);
  assert(m_AB.cols() == nBasisFuncB);
  //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValuesA = blockDataA->functionValues;
  const auto& basisFunctionValuesB = blockDataB->functionValues;
  //basis function negligebility
  const auto& negligibleA = blockDataA->negligible;
  const auto& negligibleB = blockDataB->negligible;
  //number of grid points in this block
  const unsigned int blockSize = blockDataA->functionValues.rows();
  //Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridControllerA->getFirstIndexOfBlock(iBlock);
  for_spin(m_AB,scalarPart) {
    //Pre-multiply operator and weights
    Eigen::VectorXd scalarW(blockSize);
    scalarW.array() = scalarPart_spin.segment(iGridStart,blockSize).array() * weights.segment(iGridStart,blockSize).array();
    double average = scalarW.cwiseAbs().sum();
    average /= blockSize;
    //cycle for_spin macro if non-significant
    if (average < _blockAveThreshold) return;
    //Contract with AO basis functions
    for (unsigned int nu_b = 0; nu_b < nBasisFuncB; ++nu_b) {
      if (negligibleB[nu_b]) continue;
      Eigen::VectorXd scalarWAO( scalarW.cwiseProduct(basisFunctionValuesB.col(nu_b)) );
      for (unsigned int mu_a = 0; mu_a < nBasisFuncA; ++mu_a) {
        // only loop over one half of the matrix if it is symmetrical
        if (useSym and mu_a > nu_b) break;
        if (negligibleA[mu_a]) continue;
        double preCalc = basisFunctionValuesA.col(mu_a).transpose()*scalarWAO;
        m_AB_spin(mu_a,nu_b) += preCalc;
        if(useSym and mu_a != nu_b) m_AB_spin(nu_b,mu_a) += preCalc;
      }
    }
  }; /* for_spin */
  return;
}

/* Out-dated with next version of Michael's code. */
template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addBlock(
    unsigned int iBlock,
    std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockData,
    SPMatrix<SCFMode>& matrix,
    const SpinPolarizedData<SCFMode,Eigen::VectorXd>& scalarPart,
    const Gradient<SpinPolarizedData<SCFMode,Eigen::VectorXd> >& gradientPart) {
  //Get weights
  const auto& weights = _basisFunctionOnGridControllerA->getGridController()->getWeights();
  //Number of basis functions
  unsigned int nBasisFunc;
  for_spin(matrix) {
    nBasisFunc = matrix_spin.rows();
  };
  //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValues = blockData->functionValues;
  //gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
  const auto& gradBasisFunctionValues = blockData->derivativeValues;
  //basis function negligebility
  const auto& negligible = blockData->negligible;
  //number of grid points in this block
  const unsigned int blockSize = blockData->functionValues.rows();
  //Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridControllerA->getFirstIndexOfBlock(iBlock);
  auto& gradientPartX = gradientPart.x;
  auto& gradientPartY = gradientPart.y;
  auto& gradientPartZ = gradientPart.z;
  for_spin(matrix,scalarPart,gradientPartX,gradientPartY,gradientPartZ) {
    //Pre-multiply operator and weights
    Eigen::VectorXd scalarW(blockSize);
    auto gradW = makeGradient<Eigen::VectorXd>(blockSize);
    scalarW.array() = scalarPart_spin.array() * weights.segment(iGridStart,blockSize).array();
    gradW.x.array() = gradientPartX_spin.array() * weights.segment(iGridStart,blockSize).array();
    gradW.y.array() = gradientPartY_spin.array() * weights.segment(iGridStart,blockSize).array();
    gradW.z.array() = gradientPartZ_spin.array() * weights.segment(iGridStart,blockSize).array();
    double average = scalarW.cwiseAbs().sum();
    average += gradW.x.cwiseAbs().sum();
    average += gradW.y.cwiseAbs().sum();
    average += gradW.z.cwiseAbs().sum();
    average /= blockSize;

    //cycle for_spin macro if non-significant
    if (average < _blockAveThreshold) return;
    //Contract with AO basis functions
    for (unsigned int nu = 0; nu < nBasisFunc; ++nu) {
      if (negligible[nu]) continue;
      Eigen::VectorXd scalarWAO( scalarW.cwiseProduct(basisFunctionValues.col(nu)) );
      Eigen::VectorXd gradWAO(  gradW.x.cwiseProduct((*gradBasisFunctionValues).x.col(nu)) );
      gradWAO += gradW.y.cwiseProduct((*gradBasisFunctionValues).y.col(nu));
      gradWAO += gradW.z.cwiseProduct((*gradBasisFunctionValues).z.col(nu));

      for (unsigned int mu = 0; mu < nBasisFunc; ++mu) {
        if (negligible[mu]) continue;
        matrix_spin(mu,nu) += basisFunctionValues.col(mu).transpose()*scalarWAO;
        matrix_spin(nu,mu) += gradWAO.transpose() * basisFunctionValues.col(mu);
        matrix_spin(mu,nu) += basisFunctionValues.col(mu).transpose() * gradWAO;
      }
    }
  }; /* for_spin */
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addBlock(
      unsigned int iBlock,
      std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataA,
      std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataB,
      SPMatrix<SCFMode>& m_AB,
      const GridPotential<SCFMode>& scalarPart,
      const Gradient<GridPotential<SCFMode> >& gradientPart) {
  bool useSym = isDefinedInSameBasis(*_basisFunctionOnGridControllerA,*_basisFunctionOnGridControllerB);
  //Get weights
  const auto& weights = _basisFunctionOnGridControllerA->getGridController()->getWeights();
  //Number of basis functions
  const unsigned int nBasisFuncA = _basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions();
  const unsigned int nBasisFuncB = _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions();
  //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValuesA = blockDataA->functionValues;
  const auto& basisFunctionValuesB = blockDataB->functionValues;
  //gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
  const auto& gradBasisFunctionValuesA = blockDataA->derivativeValues;
  const auto& gradBasisFunctionValuesB = blockDataB->derivativeValues;
  //basis function negligebility
  const auto& negligibleA = blockDataA->negligible;
  const auto& negligibleB = blockDataB->negligible;
  //number of grid points in this block
  const unsigned int blockSize = blockDataA->functionValues.rows();
  //Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridControllerA->getFirstIndexOfBlock(iBlock);
  auto& gradientPartX = gradientPart.x;
  auto& gradientPartY = gradientPart.y;
  auto& gradientPartZ = gradientPart.z;
  for_spin(m_AB,scalarPart,gradientPartX,gradientPartY,gradientPartZ) {
    //Pre-multiply operator and weights
    Eigen::VectorXd scalarW(blockSize);
    auto gradW = makeGradient<Eigen::VectorXd>(blockSize);
    scalarW.array() = scalarPart_spin.segment(iGridStart,blockSize).array() * weights.segment(iGridStart,blockSize).array();
    gradW.x.array() = gradientPartX_spin.segment(iGridStart,blockSize).array() * weights.segment(iGridStart,blockSize).array();
    gradW.y.array() = gradientPartY_spin.segment(iGridStart,blockSize).array() * weights.segment(iGridStart,blockSize).array();
    gradW.z.array() = gradientPartZ_spin.segment(iGridStart,blockSize).array() * weights.segment(iGridStart,blockSize).array();
    // Calculate average for prescreening.
    double average = scalarW.cwiseAbs().sum();
    average += gradW.x.cwiseAbs().sum();
    average += gradW.y.cwiseAbs().sum();
    average += gradW.z.cwiseAbs().sum();
    average /= blockSize;
    //cycle for_spin macro if non-significant
    if (average < _blockAveThreshold) return;
    //Contract with AO basis functions
    for (unsigned int nu_b = 0; nu_b < nBasisFuncB; ++nu_b) {
      if (negligibleB[nu_b]) continue;
      Eigen::VectorXd scalarWAO( scalarW.cwiseProduct(basisFunctionValuesB.col(nu_b)) );
      // weight * nabla AO_b
      Eigen::VectorXd gradWAO_nu_b(  gradW.x.cwiseProduct((*gradBasisFunctionValuesB).x.col(nu_b)) );
      gradWAO_nu_b += gradW.y.cwiseProduct((*gradBasisFunctionValuesB).y.col(nu_b));
      gradWAO_nu_b += gradW.z.cwiseProduct((*gradBasisFunctionValuesB).z.col(nu_b));

      for (unsigned int mu_a = 0; mu_a < nBasisFuncA; ++mu_a) {
        if (negligibleA[mu_a]) continue;
        if (useSym and mu_a > nu_b) break;
        // weight * nabla AO_a
        Eigen::VectorXd gradWAO_mu_a( gradW.x.cwiseProduct((*gradBasisFunctionValuesA).x.col(mu_a)) );
        gradWAO_mu_a += gradW.y.cwiseProduct((*gradBasisFunctionValuesA).y.col(mu_a));
        gradWAO_mu_a += gradW.z.cwiseProduct((*gradBasisFunctionValuesA).z.col(mu_a));

        // weight * AO_a * scalar * AO_b
        double wAOaSAOb =  basisFunctionValuesA.col(mu_a).transpose()*scalarWAO;
        m_AB_spin(mu_a,nu_b) += wAOaSAOb;
        // weight * AO_a * grad nabla AO_b
        double wAOaGNablaAOb = basisFunctionValuesA.col(mu_a).transpose() * gradWAO_nu_b;
        m_AB_spin(mu_a,nu_b) += wAOaGNablaAOb;
        // weight * nabla AO_a * grad AO_b
        double wNablaAOaGAOb = gradWAO_mu_a.transpose()*basisFunctionValuesB.col(nu_b);
        m_AB_spin(mu_a,nu_b) += wNablaAOaGAOb;
        if (useSym and mu_a != nu_b) m_AB_spin(nu_b,mu_a)+=wAOaSAOb+wAOaGNablaAOb+wNablaAOaGAOb;
      }
    }
  }; /* for_spin */
  return;
}

/* Out-dated with next version of Michael's code. */
template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addBlock(
    unsigned int iBlock,
    std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockData,
    SPMatrix<SCFMode>& matrix,
    const SpinPolarizedData<SCFMode,Eigen::VectorXd>& scalarPart) {
  //Get weights
  const auto& weights = _basisFunctionOnGridControllerA->getGridController()->getWeights();
  //Number of basis functions
  unsigned int nBasisFunc;
  for_spin(matrix) {
    nBasisFunc = matrix_spin.rows();
  };
  //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValues = blockData->functionValues;
  //basis function negligebility
  const auto& negligible = blockData->negligible;
  //number of grid points in this block
  const unsigned int blockSize = blockData->functionValues.rows();
  //Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridControllerA->getFirstIndexOfBlock(iBlock);
  for_spin(matrix,scalarPart) {
    //Pre-multiply operator and weights
    Eigen::VectorXd scalarW(blockSize);
    scalarW.array() = scalarPart_spin.array() * weights.segment(iGridStart,blockSize).array();
    double average = scalarW.cwiseAbs().sum();
    average /= blockSize;
    //cycle for_spin macro if non-significant
    if (average < _blockAveThreshold) return;
    for (unsigned int nu = 0; nu < nBasisFunc; ++nu) {
      if (negligible[nu]) continue;
      Eigen::VectorXd scalarWAO( scalarW.cwiseProduct(basisFunctionValues.col(nu)) );
      for (unsigned int mu = 0; mu < nBasisFunc; ++mu) {
        if (negligible[mu]) continue;
        matrix_spin(mu,nu) += basisFunctionValues.col(mu).transpose()*scalarWAO;
      }
    }
  }; /* for_spin */
}

template class ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED>;
template class ScalarOperatorToMatrixAdder<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
