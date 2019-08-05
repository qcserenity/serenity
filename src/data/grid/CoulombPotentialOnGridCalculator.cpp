/**
 * @file CoulombPotentialOnGridCalculator.cpp
 *
 * @date Mar 31, 2016
 * @author David Schnieders
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
#include "data/grid/CoulombPotentialOnGridCalculator.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/matrices/DensityMatrix.h"
#include "grid/Grid.h"
#include "grid/GridController.h"
#include "data/grid/GridData.h"
#include "integrals/wrappers/Libint.h"
#include "misc/SerenityError.h"
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <Eigen/Core>
#include <vector>


namespace Serenity {

CoulombPotentialOnGridCalculator::CoulombPotentialOnGridCalculator(
    std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController):
        _basisFuncOnGridController(basisFuncOnGridController){

}


template<Options::SCF_MODES SCF_MODE>
void CoulombPotentialOnGridCalculator::calculateElectronElectron(GridPotential<RESTRICTED>& result,
    const DensityMatrix<SCF_MODE>& densMat){

  auto nGridPoints=_basisFuncOnGridController->getNGridPoints();
  auto nBlocks=_basisFuncOnGridController->getNBlocks();
  auto basisController=_basisFuncOnGridController->getBasisController();
  auto basis=basisController->getBasis();
  const auto& gridPoints = _basisFuncOnGridController->getGridController()->getGridPoints();


  auto& libint = Libint::getInstance();
  libint.keepEngines(libint2::Operator::nuclear,0,2);

  //go through grid blockwise
#pragma omp for schedule(dynamic)
  for (unsigned int block=0; block<nBlocks; block++){
    unsigned int blockEnd=0;
    if(block==nBlocks-1){
      blockEnd=nGridPoints;
    }
    else{
      blockEnd=_basisFuncOnGridController->getFirstIndexOfBlock(block+1);
    }
    unsigned int firstOfBlock=_basisFuncOnGridController->getFirstIndexOfBlock(block);
    for(unsigned int gridpoint=firstOfBlock; gridpoint<blockEnd; gridpoint++){
      /*
       * integrals should be calculated to a point charge of -1
       * at the position of the current grid point
       */
      std::vector<std::pair<double,std::array<double,3>>> point =
      {{-1.0,{{gridPoints(0,gridpoint),gridPoints(1,gridpoint),gridPoints(2,gridpoint)}} }};

      Eigen::MatrixXd ints = libint.compute1eInts(libint2::Operator::nuclear,basisController,point);

      double pot = 0.0;
      for_spin(densMat){
        pot+=densMat_spin.cwiseProduct(ints).sum();
      };
      //... and store in the GridPotential
      result[gridpoint]+=pot;
    }//gridPoint
  }//block


}

void CoulombPotentialOnGridCalculator::calculateElectronNuclei(
    GridPotential<RESTRICTED>& result,
    const std::vector<std::shared_ptr<Atom> >& atoms){

  auto nGridPoints=_basisFuncOnGridController->getNGridPoints();
  auto nBlocks=_basisFuncOnGridController->getNBlocks();
  auto gridController=_basisFuncOnGridController->getGridController();
  const auto& gridPoints = gridController->getGridPoints();



  //go through grid blockwise
#pragma omp for schedule(dynamic)
  for (unsigned int block=0; block<nBlocks; block++){
    unsigned int blockEnd=0;
    if(block==nBlocks-1){
      blockEnd=nGridPoints;
    }
    else{
      blockEnd=_basisFuncOnGridController->getFirstIndexOfBlock(block+1);
    }
    unsigned int firstOfBlock=_basisFuncOnGridController->getFirstIndexOfBlock(block);
    for(unsigned int gridpoint=firstOfBlock; gridpoint<blockEnd; gridpoint++){
      //Eigen vector class for easy vector calculation
      Eigen::Vector3d gridCoord = gridPoints.col(gridpoint);
      //calculate sum_A [-charge_A/(|r-R_A|)] on every gridpoint
      for (auto atom : atoms){
        //Eigen vector class for easy vector calculation
        Eigen::Vector3d atomCoord(
            atom->getX(),
            atom->getY(),
            atom->getZ());
        /*
         * precalculate the potential before storing it in the
         * spin polarized GridPotential (to prevent double calculation)
         */
        /*
         * TODO I think the effective core potential must be added here as well, but I am not sure
         * about this. To be on the safe side, I thus throw an exception here.
         */
        if (atom->usesECP()) throw SerenityError(
          "The Coulomb potential on a grid is to be calculated in the presence of an Atom using"
          " an effective core potential. This is probably not done correctly and thus forbidden"
          " until the issue has been fixed.");
        double tmpPot=-atom->getEffectiveCharge()/((gridCoord-atomCoord).norm());
        //store
        result[gridpoint]+=tmpPot;
      }
    }//gridPoint
  }//block

}

template
void CoulombPotentialOnGridCalculator::calculateElectronElectron<Options::SCF_MODES::RESTRICTED>(
    GridPotential<Options::SCF_MODES::RESTRICTED>& result,
    const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densMat);
template
void CoulombPotentialOnGridCalculator::calculateElectronElectron<Options::SCF_MODES::UNRESTRICTED>(
    GridPotential<Options::SCF_MODES::RESTRICTED>& result,
    const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densMat);

} /* namespace Serenity */
