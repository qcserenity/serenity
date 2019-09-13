/**
 * @file DipoleIntegrals.cpp
 * @author Niklas Niemeyer
 *
 * @date Dec. 17, 2018
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
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"

/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
DipoleIntegrals<SCFMode>::DipoleIntegrals(
    std::vector<std::shared_ptr<LRSCFController<SCFMode> > > lrscf,
    Point gaugeOrigin) :
  _lrscf(lrscf),
  _gaugeOrigin(gaugeOrigin){
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<const Eigen::MatrixXd> DipoleIntegrals<SCFMode>::getLengths(Options::INTEGRAL_TYPE type){
  if(type == Options::INTEGRAL_TYPE::ANALYTICAL){
    if(!_lengthsAnalytical){
      this->computeIntegralsAnalytically();
    }
    return _lengthsAnalytical;
  }else if(type == Options::INTEGRAL_TYPE::NUMERICAL){
    if(!_lengthsNumerical){
      this->computeIntegralsNumerically();
    }
    return _lengthsNumerical;
  }else{
      assert(false);
      return nullptr;
  }
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<const Eigen::MatrixXd> DipoleIntegrals<SCFMode>::getVelocities(Options::INTEGRAL_TYPE type){
  if(type == Options::INTEGRAL_TYPE::ANALYTICAL){
    if(!_velocitiesAnalytical){
      this->computeIntegralsAnalytically();
    }
    return _velocitiesAnalytical;
  }else if(type == Options::INTEGRAL_TYPE::NUMERICAL){
    if(!_velocitiesNumerical){
      this->computeIntegralsNumerically();
    }
    return _velocitiesNumerical;
  }else{
      assert(false);
      return nullptr;
  }
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<const Eigen::MatrixXd> DipoleIntegrals<SCFMode>::getMagnetics(Options::INTEGRAL_TYPE type){
  if(type == Options::INTEGRAL_TYPE::ANALYTICAL){
    if(!_magneticsAnalytical){
      this->computeIntegralsAnalytically();
    }
    return _magneticsAnalytical;
  }else if(type == Options::INTEGRAL_TYPE::NUMERICAL){
    if(!_magneticsNumerical){
      this->computeIntegralsNumerically();
    }
    return _magneticsNumerical;
  }else{
      assert(false);
      return nullptr;
  }
}

template<Options::SCF_MODES SCFMode>
void DipoleIntegrals<SCFMode>::computeIntegralsAnalytically() {

  std::vector<std::vector<Eigen::MatrixXd> > len;
  std::vector<std::vector<Eigen::MatrixXd> > mag;
  std::vector<std::vector<Eigen::MatrixXd> > vel;

  Timings::takeTime("LRSCF - Dip Ints (Analytical)"); 

  //Compute all integrals analytically (electric dipole (len and vel rep) and magnetic dipole (mag))
  auto &libint = Libint::getInstance();
  auto &libintDeriv = Libint::getInstance();
  libint.initialize(libint2::Operator::emultipole1, 0, 2, _gaugeOrigin);
  libintDeriv.initialize(libint2::Operator::emultipole1, 1, 2, _gaugeOrigin);
  for(unsigned int iSys = 0; iSys < _lrscf.size(); ++iSys) {
    auto basisController = _lrscf[iSys]->getBasisController();
    auto basis = basisController->getBasis();
    const unsigned int nBFs = basisController->getNBasisFunctions();
    //Append a set of x,y,z matrices for each system
    len.push_back(std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(nBFs, nBFs)));
    vel.push_back(std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(nBFs, nBFs)));
    mag.push_back(std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(nBFs, nBFs)));
    for (unsigned int i = 0; i < basis.size(); i++) {
      for (unsigned int j = 0; j < basis.size(); j++) {
        Eigen::MatrixXd ints;
        if (libint.compute(libint2::Operator::emultipole1, 0, *basis[i], *basis[j], ints)) {
          for (unsigned int k = 0; k < basis[i]->getNContracted(); k++) {
            auto mu = basisController->extendedIndex(i) + k;
            for (unsigned int l = 0; l < basis[j]->getNContracted(); l++) {
              auto nu = basisController->extendedIndex(j) + l;
              const unsigned int nj = basis[j]->getNContracted();
              //Matrix ints has 4 columns:
              // (i|1|j) -> (:,0),  (i|x|j) -> (:,1),   (i|y|j) -> (:,2),   (i|z|j) -> (:,3)
              //Electric dipole integrals in length rep
              len[iSys][0](mu, nu) = ints((nj * k + l), 1);
              len[iSys][1](mu, nu) = ints((nj * k + l), 2);
              len[iSys][2](mu, nu) = ints((nj * k + l), 3);
            }
          }
        }
        Eigen::MatrixXd intsDeriv;
        if (libintDeriv.compute(libint2::Operator::emultipole1, 1, *basis[i], *basis[j], intsDeriv)) {
          for (unsigned int k = 0; k < basis[i]->getNContracted(); k++) {
            auto mu = basisController->extendedIndex(i) + k;
            for (unsigned int l = 0; l < basis[j]->getNContracted(); l++) {
              auto nu = basisController->extendedIndex(j) + l;
              const unsigned int nj = basis[j]->getNContracted();
              //Matrix intsDeriv has 24 columns (Note: derivative-level 1!)
              // (i|1|dx j) -> (:, 0),   (i|x|dx j) -> (:, 1),   (i|y|dx j) -> (:, 2),   (i|z|dx j) -> (:, 3),   
              // (i|1|dy j) -> (:, 4),   (i|x|dy j) -> (:, 5),   (i|y|dy j) -> (:, 6),   (i|z|dy j) -> (:, 7),   
              // (i|1|dz j) -> (:, 8),   (i|x|dz j) -> (:, 9),   (i|y|dz j) -> (:,10),   (i|z|dz j) -> (:,11),   

              // (dx i|1|j) -> (:,12),   (dx i|x|j) -> (:,13),   (dx i|y|j) -> (:,14),   (dx i|z|j) -> (:,15),   
              // (dy i|1|j) -> (:,16),   (dy i|x|j) -> (:,17),   (dy i|y|j) -> (:,18),   (dy i|z|j) -> (:,19),   
              // (dz i|1|j) -> (:,20),   (dz i|x|j) -> (:,21),   (dz i|y|j) -> (:,22),   (dz i|z|j) -> (:,23).

              //Electric dipole integrals in velocity rep (i|nabla|j)
              vel[iSys][0](mu, nu) = intsDeriv((nj * k + l), 0);
              vel[iSys][1](mu, nu) = intsDeriv((nj * k + l), 4);
              vel[iSys][2](mu, nu) = intsDeriv((nj * k + l), 8);
              //Magnetic dipole integrals (i|r x nabla|j)
              mag[iSys][0](mu, nu) = intsDeriv((nj * k + l), 10) - intsDeriv((nj * k + l), 7);
              mag[iSys][1](mu, nu) = intsDeriv((nj * k + l),  3) - intsDeriv((nj * k + l), 9);
              mag[iSys][2](mu, nu) = intsDeriv((nj * k + l),  5) - intsDeriv((nj * k + l), 2);
            }
          }
        }
      }
    }
  }
  libint.finalize(libint2::Operator::emultipole1, 0, 2);
  libintDeriv.finalize(libint2::Operator::emultipole1, 1, 2);
  /* Compute all integrals analytically */

  Timings::timeTaken("LRSCF - Dip Ints (Analytical)");

  //Perform AO -> MO transformation and assign shared pointer
  _lengthsAnalytical    = std::make_shared<const Eigen::MatrixXd>(ao2mo(len));
  _velocitiesAnalytical = std::make_shared<const Eigen::MatrixXd>(ao2mo(vel));
  _magneticsAnalytical  = std::make_shared<const Eigen::MatrixXd>(ao2mo(mag));
}

template<Options::SCF_MODES SCFMode>
void DipoleIntegrals<SCFMode>::computeIntegralsNumerically() {

  std::vector<std::vector<Eigen::MatrixXd> > len;
  std::vector<std::vector<Eigen::MatrixXd> > mag;
  std::vector<std::vector<Eigen::MatrixXd> > vel;

  Timings::takeTime("LRSCF -  Dip Ints (Numerical)"); 
  
  //Compute all integrals numerically (electric dipole (length & velocity rep) and magnetic dipole (mag))
  //(Note) I'm gonna keep this functionality here for testing purposes (NN)
  for(unsigned int iSys = 0; iSys < _lrscf.size(); ++iSys) {

    auto gridController = _lrscf[iSys]->getGridController();
    auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
            _lrscf[iSys]->getSysSettings().grid.blocksize,
            _lrscf[iSys]->getSysSettings().grid.basFuncRadialThreshold,
            1,
            _lrscf[iSys]->getBasisController(),
            _lrscf[iSys]->getGridController());
    //Get number of grid blocks
    const unsigned int nBlocks = basisFunctionOnGridController->getNBlocks();
    //Get weights
    const auto &weights = gridController->getWeights();
    //Get all grid points
    const auto &gridPoints = gridController->getGridPoints();
    //Get number of basis functions
    const unsigned int nBFs = _lrscf[iSys]->getBasisController()->getNBasisFunctions();
    //Append a set of x,y,z matrices for each system
    len.push_back(std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(nBFs, nBFs)));
    vel.push_back(std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(nBFs, nBFs)));
    mag.push_back(std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(nBFs, nBFs)));
    //Perform blockwise numerical integration
    //ToDo PARALLELIZE!!
    for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
      //Data for this block
      auto &thisBlockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
      //function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
      const auto &basisFunctionValues = thisBlockData->functionValues;
      //gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
      const auto &gradBasisFunctionValues = thisBlockData->derivativeValues;
      //basis function negligibility
      const auto &negligible = thisBlockData->negligible;
      //number of grid points in this block
      const unsigned int blockSize = thisBlockData->functionValues.rows();
      //Get first index of this Block
      const unsigned int iGridStart = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
      //Cut the gridPoints matrix into blockPoints segments
      const auto blockPoints = gridPoints.block(0, iGridStart, 3, blockSize).transpose();
      //Add to matrix
      for (unsigned int nu = 0; nu < nBFs; ++nu) {
        if (negligible[nu]) continue;
        Eigen::VectorXd nuW(basisFunctionValues.col(nu).cwiseProduct(weights.segment(iGridStart, blockSize)));
        for (unsigned int mu = 0; mu < nBFs; ++mu) {
          if (negligible[mu]) continue;

          const auto& gradBFx = (*gradBasisFunctionValues).x.col(mu);
          const auto& gradBFy = (*gradBasisFunctionValues).y.col(mu);
          const auto& gradBFz = (*gradBasisFunctionValues).z.col(mu);

          //Electric integrals length
          len[iSys][0](nu, mu) += nuW.dot(basisFunctionValues.col(mu).cwiseProduct(blockPoints.col(0)));
          len[iSys][1](nu, mu) += nuW.dot(basisFunctionValues.col(mu).cwiseProduct(blockPoints.col(1)));
          len[iSys][2](nu, mu) += nuW.dot(basisFunctionValues.col(mu).cwiseProduct(blockPoints.col(2)));
          /* Electric integrals length */
          
          //Electric integrals velocity
          vel[iSys][0](nu, mu) += nuW.dot(gradBFx);
          vel[iSys][1](nu, mu) += nuW.dot(gradBFy);
          vel[iSys][2](nu, mu) += nuW.dot(gradBFz);
          /* Electric integrals velocity */

          //Magnetic Integrals
          mag[iSys][0](nu, mu) += nuW.dot(blockPoints.col(1).cwiseProduct(gradBFz) - blockPoints.col(2).cwiseProduct(gradBFy));
          mag[iSys][1](nu, mu) += nuW.dot(blockPoints.col(2).cwiseProduct(gradBFx) - blockPoints.col(0).cwiseProduct(gradBFz));
          mag[iSys][2](nu, mu) += nuW.dot(blockPoints.col(0).cwiseProduct(gradBFy) - blockPoints.col(1).cwiseProduct(gradBFx));
          /* Magnetic Integrals */
        }
      }
    }
  }
  /* Numerical integration for electric dipole integrals (length & velocity rep) and magnetic integrals */

  Timings::timeTaken("LRSCF -  Dip Ints (Numerical)");

  //Perform AO -> MO transformation and assign shared pointer
  _lengthsNumerical    = std::make_shared<const Eigen::MatrixXd>(ao2mo(len));
  _velocitiesNumerical = std::make_shared<const Eigen::MatrixXd>(ao2mo(vel));
  _magneticsNumerical  = std::make_shared<const Eigen::MatrixXd>(ao2mo(mag));
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd DipoleIntegrals<SCFMode>::ao2mo(std::vector<std::vector<Eigen::MatrixXd> >& ao_xyz) {

  //Determine length of dipole matrices
  unsigned int nDimension = 0;
  for (unsigned int iSys = 0; iSys < _lrscf.size(); ++iSys) {
    auto nOccupied = _lrscf[iSys]->getNOccupied();
    auto nVirtual = _lrscf[iSys]->getNVirtual();
    for_spin(nOccupied,nVirtual) {nDimension += nOccupied_spin*nVirtual_spin;};
  }

  Eigen::MatrixXd dipoles = Eigen::MatrixXd::Zero(nDimension,3);

  unsigned int iStartSys = 0;
  for (unsigned int iSys = 0; iSys < _lrscf.size(); ++iSys) {
    auto coeff = _lrscf[iSys]->getCoefficients();
    auto nOccupied = _lrscf[iSys]->getNOccupied();
    auto nVirtual = _lrscf[iSys]->getNVirtual();
    for (unsigned int iMu = 0; iMu < 3; ++iMu) {
      unsigned int iStartSpin = 0;
      for_spin (coeff,nOccupied,nVirtual) {
        Eigen::MatrixXd tmp = coeff_spin.block(0,0,coeff_spin.rows(),nOccupied_spin).transpose()
                            * ao_xyz[iSys][iMu]
                            * coeff_spin.block(0,nOccupied_spin,coeff_spin.rows(),nVirtual_spin);
        tmp.transposeInPlace();
        dipoles.block(iStartSys+iStartSpin,iMu,nOccupied_spin*nVirtual_spin,1)
            = Eigen::Map<Eigen::VectorXd>(tmp.data(),tmp.cols()*tmp.rows());
        iStartSpin += nOccupied_spin*nVirtual_spin;
      };
    }
    for_spin(nOccupied,nVirtual) {iStartSys += nOccupied_spin*nVirtual_spin;};
  }

  return dipoles;
}

template class DipoleIntegrals<Options::SCF_MODES::RESTRICTED>;
template class DipoleIntegrals<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
