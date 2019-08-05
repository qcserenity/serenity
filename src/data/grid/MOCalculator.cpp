/**
 * @file   MOCalculator.cpp
 *
 * @date   Apr 22, 2014
 * @author Thomas Dresselhaus
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
#include "data/grid/MOCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "grid/GridController.h"
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <iostream>


namespace Serenity {
using namespace std;

MOCalculator::MOCalculator(
    std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController):
  _basisFuncOnGridController(basisFuncOnGridController){

}


Eigen::VectorXd MOCalculator::calcMOValuesOnGrid(
    Eigen::VectorXd& coefficients) {

  //get some useful things
  unsigned int nBlocks=_basisFuncOnGridController->getNBlocks();
  unsigned int nGridPts=_basisFuncOnGridController->getNGridPoints();
  unsigned int nBasisFunc=_basisFuncOnGridController->getNBasisFunctions();

  //prepare return value
  Eigen::VectorXd moOnGrid(nGridPts);
  moOnGrid.setZero();

  //loop over all blocks and calculate the MO value on every gridpoint
 #pragma omp parallel for schedule(dynamic)
    for (unsigned int block=0; block<nBlocks; block++){
      const auto& basisFuncOnGrid = _basisFuncOnGridController->getBlockOnGridData(block);
      unsigned int blockEnd=0;
      if(block==nBlocks-1){
        blockEnd=nGridPts;
      }
      else{
        blockEnd=_basisFuncOnGridController->getFirstIndexOfBlock(block+1);
      }
      unsigned int firstOfBlock=_basisFuncOnGridController->getFirstIndexOfBlock(block);
      for (unsigned int gridPt=firstOfBlock; gridPt<blockEnd; gridPt++){
        for (unsigned int ao=0; ao<nBasisFunc; ao++){
          moOnGrid[gridPt]+=basisFuncOnGrid->functionValues(gridPt-firstOfBlock, ao)
                    *coefficients[ao];
        }
      }
    }
  return moOnGrid;

}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Matrix<double>>MOCalculator::calcOccMOValuesOnGrid(
    CoefficientMatrix<SCFMode>& coefficientMatrix,
    SpinPolarizedData<SCFMode, unsigned int>& nOccOrbs){

  //get some useful values
  unsigned int nGridPts=_basisFuncOnGridController->getNGridPoints();

  //prepare return value
  SpinPolarizedData<SCFMode, Matrix<double>> mosOnGrid;

  for_spin(mosOnGrid,nOccOrbs){
    mosOnGrid_spin.resize(nGridPts,nOccOrbs_spin);
    mosOnGrid_spin.setZero();
  };

  //call the function for every occupied MO
  for_spin(mosOnGrid, nOccOrbs, coefficientMatrix){
    for (unsigned int mo=0; mo<nOccOrbs_spin; mo++){
      Eigen::VectorXd coeffs=coefficientMatrix_spin.col(mo);
      auto tmpMoOnGrid=calcMOValuesOnGrid(coeffs);
      for (unsigned int gridPt=0; gridPt<nGridPts; gridPt++){
        mosOnGrid_spin(gridPt,mo)=tmpMoOnGrid[gridPt];
      }
    }
  };

  return mosOnGrid;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Matrix<double>>MOCalculator::calcAllMOValuesOnGrid(
      CoefficientMatrix<SCFMode>& coefficientMatrix){

  //call the function above with nBasisFunc as nOccOrbs -> all MOs are returned
  SpinPolarizedData<SCFMode, unsigned int> orbs(_basisFuncOnGridController->getNBasisFunctions());
  SpinPolarizedData<SCFMode, Matrix<double> > mosOnGrid=calcOccMOValuesOnGrid<SCFMode>(coefficientMatrix,orbs);

  return mosOnGrid;
}


/*
 * Calculation of the kinetic energy density.
 */
template<Options::SCF_MODES SCFMode>
GridData<SCFMode> MOCalculator::calcKineticEnergyDensityOnGrid(
    CoefficientMatrix<SCFMode> coefficientMatrix,
    SpinPolarizedData<SCFMode, unsigned int> nOccOrbs) {

  //get useful values
  unsigned int nBlocks=_basisFuncOnGridController->getNBlocks();
  unsigned int nGridPts=_basisFuncOnGridController->getNGridPoints();
  unsigned int nBasisFunc=_basisFuncOnGridController->getNBasisFunctions();

  //prepare return value
  GridData<SCFMode> kineticEnergyDensity(_basisFuncOnGridController->
      getGridController());

  //actual calculation
  for_spin(coefficientMatrix,nOccOrbs,kineticEnergyDensity){
    //Loop over the occupied orbitals
    for (unsigned int i=0; i<nOccOrbs_spin; ++i){
      std::vector<double> gradPhi2(nGridPts, 0.0);
      std::vector<double> gradPhipartX(nGridPts, 0.0);
      std::vector<double> gradPhipartY(nGridPts, 0.0);
      std::vector<double> gradPhipartZ(nGridPts, 0.0);

      // Loop over Grid. Data is accessed blockwise.
      for (unsigned int blockIndex = 0 ; blockIndex < nBlocks; ++blockIndex) {
        const auto& blockOnGridData = _basisFuncOnGridController->getBlockOnGridData(blockIndex);
        const auto& bfDerivatives = *blockOnGridData->derivativeValues;
        const unsigned int blockSize = blockOnGridData->functionValues.rows();

        // Loop over data entries in the current block
        for (unsigned int blockData = 0; blockData < blockSize; ++blockData) {
          const unsigned int pointIndex = _basisFuncOnGridController->
              getFirstIndexOfBlock(blockIndex)+blockData;

          //Loop over basis functions
          for (unsigned int mu = 0; mu < nBasisFunc; ++mu) {
            gradPhipartX[pointIndex]+=coefficientMatrix_spin(mu,i)*bfDerivatives.x(blockData,mu);
            gradPhipartY[pointIndex]+=coefficientMatrix_spin(mu,i)*bfDerivatives.y(blockData,mu);
            gradPhipartZ[pointIndex]+=coefficientMatrix_spin(mu,i)*bfDerivatives.z(blockData,mu);
          }

          gradPhi2[pointIndex]=pow(gradPhipartX[pointIndex],2.0)+
              pow(gradPhipartY[pointIndex],2.0)+pow(gradPhipartZ[pointIndex],2.0);

          kineticEnergyDensity_spin[pointIndex]+=gradPhi2[pointIndex];

          if(SCFMode==Options::SCF_MODES::RESTRICTED){
            kineticEnergyDensity_spin[pointIndex]*=2.0;
          }
        } /* Loop over the entries in the current block */
      } /* Loop over the grid blocks.*/
    } /* Loop over the occupied orbitals */
  };
  return kineticEnergyDensity;
}

template
SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Matrix<double>>MOCalculator::calcOccMOValuesOnGrid<Options::SCF_MODES::RESTRICTED>(
    CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& coefficientMatrix,
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>& nOccOrbs);
template
SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Matrix<double>>MOCalculator::calcOccMOValuesOnGrid<Options::SCF_MODES::UNRESTRICTED>(
    CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>& coefficientMatrix,
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int>& nOccOrbs);
template
SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Matrix<double>>MOCalculator::calcAllMOValuesOnGrid<Options::SCF_MODES::RESTRICTED>(
      CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& coefficientMatrix);
template
SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Matrix<double>>MOCalculator::calcAllMOValuesOnGrid<Options::SCF_MODES::UNRESTRICTED>(
      CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>& coefficientMatrix);
template
GridData<Options::SCF_MODES::RESTRICTED> MOCalculator::calcKineticEnergyDensityOnGrid<Options::SCF_MODES::RESTRICTED>(
    CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficientMatrix,
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int> nOccOrbs);
template
GridData<Options::SCF_MODES::UNRESTRICTED> MOCalculator::calcKineticEnergyDensityOnGrid<Options::SCF_MODES::UNRESTRICTED>(
    CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coefficientMatrix,
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int> nOccOrbs);

} /* namespace Serenity */
