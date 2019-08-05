/**
 * @file NAddFuncPotential.cpp
 *
 * @date Nov 24, 2016
 * @author Jan Unsleber
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
#include "potentials/NAddFuncPotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridController.h"
#include "data/grid/DensityOnGridFactory.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "energies/EnergyContributions.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/grid/SupersystemDensityOnGridController.h"
#include "misc/Timing.h"
#include "dft/functionals/wrappers/XCFun.h"

namespace Serenity {

template <Options::SCF_MODES SCFMode>
NAddFuncPotential<SCFMode>::NAddFuncPotential(std::shared_ptr<SystemController> system,
    std::shared_ptr<DensityMatrixController<SCFMode> > activeDMat,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > envDMats,
    std::shared_ptr<GridController> grid,
    Functional functional,
    std::pair<bool, std::vector<std::shared_ptr<EnergyComponentController> > > eCon,
    bool evaluateEnergy,
    bool evaluateExactX):
    Potential<SCFMode>(activeDMat->getDensityMatrix().getBasisController()),
    _system(system),
    _actDMatController(activeDMat),
    _envDMatController(envDMats),
    _grid(grid),
    _functional(functional),
    _potential(nullptr),
    _densOnGridControllers(),
    _supersysDensOnGridController(nullptr),
    _gridToMatrix(nullptr),
    _basisFunctionOnGridController(nullptr),
    _excPot(nullptr),
    _energy(0.0),
    _helper(),
    _energyController(eCon.second),
    _isXC(eCon.first),
    _evaluateEnergy(evaluateEnergy),
    _evaluateExactX(evaluateExactX){
  assert(_energyController.size()==0 or _energyController.size()==(envDMats.size()+1));
  this->_grid->addSensitiveObject(ObjectSensitiveClass<Grid>::_self);
  this->_actDMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode> >::_self);
  for(auto& envMat : _envDMatController){
    envMat->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode> >::_self);
  }
  _basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _system->getSettings(), this->_basis, this->_grid);
  _densOnGridControllers.push_back(DensityOnGridFactory<SCFMode>::produce(
      _actDMatController,_grid,1,_system->getSettings()));
  _gridToMatrix = std::make_shared<ScalarOperatorToMatrixAdder<SCFMode> >(
      _basisFunctionOnGridController, _system->getSettings().grid.blockAveThreshold);

  unsigned int count =0;
  for(auto& envMat : _envDMatController){
    auto dens = DensityOnGridFactory<SCFMode>::produce(
        envMat,_grid,1,_system->getSettings());
    _densOnGridControllers.push_back(dens);
    _helper.push_back(nullptr);
    _helper[count].reset(new NAddEnergyHelper<SCFMode>(_functional,dens));
    count++;
  }

  _supersysDensOnGridController = std::make_shared<SupersystemDensityOnGridController<SCFMode> >(_densOnGridControllers);
}


template <Options::SCF_MODES SCFMode> FockMatrix<SCFMode>&
NAddFuncPotential<SCFMode>::getMatrix(){

  Timings::takeTime("FDE - Non-Add. Func. Pot.");
  if (!_potential){

    XCFun<SCFMode> xcfun(128);
    auto superFuncDat = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENTS,
        _functional,
        _supersysDensOnGridController,1);

    auto activeFuncDat = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENTS,
        _functional,
        _densOnGridControllers[0],1);

    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot){
      pot_spin.setZero();
    };

    if (_functional.isHybrid() and _evaluateExactX){
      if (!_excPot){
        _excPot.reset(new ExchangeInteractionPotential<SCFMode>(this->_basis,_envDMatController,
            _functional.getHfExchangeRatio(),_system->getSettings().basis.integralThreshold));
      }
      auto& xcMat = _excPot->getMatrix();
      for_spin(pot,xcMat){
        pot_spin += xcMat_spin;
      };
    }
    if (_functional.getFunctionalClass() == FUNCTIONAL_CLASSES::LDA){
      *superFuncDat.dFdRho -= *activeFuncDat.dFdRho;
      _gridToMatrix->addScalarOperatorToMatrix(pot,*superFuncDat.dFdRho);
    } else if (_functional.getFunctionalClass() == FUNCTIONAL_CLASSES::GGA) {
      *superFuncDat.dFdRho -= *activeFuncDat.dFdRho;
      superFuncDat.dFdGradRho->x -= activeFuncDat.dFdGradRho->x;
      superFuncDat.dFdGradRho->y -= activeFuncDat.dFdGradRho->y;
      superFuncDat.dFdGradRho->z -= activeFuncDat.dFdGradRho->z;
      _gridToMatrix->addScalarOperatorToMatrix(pot,*superFuncDat.dFdRho,*superFuncDat.dFdGradRho);
    } else if (_functional.getFunctionalClass() == FUNCTIONAL_CLASSES::NONE){
    } else {
      assert(false && "Unsupported functional type as nadd. functional!");
    }

    /*
     * This handling of the energies is not the finest piece of code, but
     * it makes FaT run efficiently.
     * The idea is the following:
     * If there is a supersystem grid present (no cut off) the energies
     * that are subtracted for each environment system can be reused rather
     * than recalculated.
     * The fact that the grid is a supersystem grid is handled in the FDE/FaT
     * tasks, only if this is the case, the variables used here are actually present.
     * (I.e. list of energy controllers is longer than 0.)
     * Calculated energies are then stored in the energy controllers.
     * If this can not be done the energy of each subsystem on the supersystem
     * grid has to be evaluated for each FDE run in the FaT run, this is costly
     * for large systems or the case of many systems.
     */
    if(_evaluateEnergy){
      _energy = superFuncDat.energy - activeFuncDat.energy;
      if (_energyController.size()>0){
        /*
         *  clear the energies of this subsystem so they are evaluated
         *  anew when the final density is calculated
         */
        if(_isXC){
          _energyController[0]->clearNAddXCSub();
        } else {
          _energyController[0]->clearNAddKinSub();
        }
        // subtract all other energies if present, if not then calculate them
        for (unsigned int i=1; i<_energyController.size();i++){
          if(_isXC){
            if (_energyController[i]->checkNAddXCSub()){
              _energy -= _energyController[i]->getNAddXCSub();
            } else {
              _energyController[i]->setNAddXCSub(_helper[i-1]->getEnergy());
              _energy -= _energyController[i]->getNAddXCSub();
            }
          } else {
            if (_energyController[i]->checkNAddKinSub()){
              _energy -= _energyController[i]->getNAddKinSub();
            } else {
              _energyController[i]->setNAddKinSub(_helper[i-1]->getEnergy());
              _energy -= _energyController[i]->getNAddKinSub();
            }
          }
        }
      } else {
        for(auto& h : _helper){
          _energy -= h->getEnergy();
        }
      }
    }
  }
  Timings::timeTaken("FDE - Non-Add. Func. Pot.");
  return *_potential;
}

template <Options::SCF_MODES SCFMode> double
NAddFuncPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P){
  assert(P.getBasisController() == this->_basis);

  if(!_potential){
    (void)this->getMatrix();
  }

  if (_functional.isHybrid() and _evaluateExactX){
    if (!_excPot){
      _excPot.reset(new ExchangeInteractionPotential<SCFMode>(this->_basis,_envDMatController,
          _functional.getHfExchangeRatio(),_system->getSettings().basis.integralThreshold));
    }
    this->_system->template getElectronicStructure<SCFMode>()->getEnergyComponentController()->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_EXACT_EXCHANGE,_excPot->getEnergy(P));
  }else{
    this->_system->template getElectronicStructure<SCFMode>()->getEnergyComponentController()->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_EXACT_EXCHANGE,0.0);
  }

  return _energy;
};

template <Options::SCF_MODES SCFMode> Eigen::MatrixXd
NAddFuncPotential<SCFMode>::getGeomGradients(){

  XCFun<SCFMode> xcfun(128);
  auto funcData = xcfun.calcData(
      FUNCTIONAL_DATA_TYPE::GRADIENTS,
      _functional,
      _supersysDensOnGridController,
      1);
  auto activeFuncDat = xcfun.calcData(
      FUNCTIONAL_DATA_TYPE::GRADIENTS,
      _functional,
      _densOnGridControllers[0],
      1);

  unsigned int derivorder = 0;
  *(funcData.dFdRho) -= *(activeFuncDat.dFdRho);
  if (_functional.getFunctionalClass() == FUNCTIONAL_CLASSES::GGA) {
    funcData.dFdGradRho->x -= activeFuncDat.dFdGradRho->x;
    funcData.dFdGradRho->y -= activeFuncDat.dFdGradRho->y;
    funcData.dFdGradRho->z -= activeFuncDat.dFdGradRho->z;
    derivorder = _basisFunctionOnGridController->getHighestDerivative();
    if (derivorder<2) _basisFunctionOnGridController->setHighestDerivative(2);
  }

  unsigned int nAtoms = _system->getAtoms().size();
  auto basisController = _system->getAtomCenteredBasisController();
  const unsigned int nBasisFunctions = basisController->getNBasisFunctions();

  // Create mapping
  auto orbitalSet = _system->getActiveOrbitalController<SCFMode>();
  std::vector<unsigned int> mapping(nBasisFunctions);
  std::vector<bool> hasElementBeenSet(nBasisFunctions, false);
  const auto& basisIndices = basisController->getBasisIndices();
  // Vector to check whether ALL basis function shells are assigned to an atom index
  for (unsigned int iAtom=0; iAtom < basisIndices.size(); ++iAtom) {
    const unsigned int firstIndex = basisIndices[iAtom].first;
    const unsigned int endIndex = basisIndices[iAtom].second;
    for (unsigned int mu=firstIndex; mu < endIndex; ++mu) {
      mapping[mu] = iAtom;
      hasElementBeenSet[mu] = true;
    }
  }

  // Check
  for (bool x : hasElementBeenSet){
    if (not x)
      throw SerenityError("NAddFuncPotential: Missed gradient element in gradient evaluation.");
  } 

  // Get potential and (in case of GGAs) potential gradients
  const auto& potential = *(funcData.dFdRho);
  auto potentialGradients = funcData.dFdGradRho;

  // Get grid and BFs on it TODO check whether derivative level needs to be modified
  unsigned int nBlocks = _basisFunctionOnGridController->getNBlocks();
  const auto& weights = _basisFunctionOnGridController->getGridController()->getWeights();

  // Get density Matrix
  auto densityMatrix(_system->getElectronicStructure<SCFMode>()->getDensityMatrix());

  Eigen::MatrixXd gradientContr = Eigen::MatrixXd::Zero(nAtoms,3);
#pragma omp parallel
  {
    Eigen::MatrixXd gradientContrPriv(nAtoms,3);
    gradientContrPriv.setZero();
    auto gpp = gradientContrPriv.data();
#pragma omp for schedule(dynamic)
    /*
     * Loop over Grid. Data is accessed blockwise.
     */
    for (unsigned int blockIndex = 0 ; blockIndex < nBlocks; ++blockIndex) {
      const auto& blockOnGridData = _basisFunctionOnGridController->getBlockOnGridData(blockIndex);
      const auto& bfValues = blockOnGridData->functionValues;
      assert(blockOnGridData->derivativeValues);
      const auto& bfDerivatives = *blockOnGridData->derivativeValues;
      const unsigned int blockSize = blockOnGridData->functionValues.rows();
      /*
       * Check ExchangeCorrelationDerivativeCalculator.h for formulas and reference
       */
      for (unsigned int mu = 0; mu < nBasisFunctions; ++mu) {
        auto muAtom = mapping[mu];
        if (blockOnGridData->negligible[mu]) continue;
        for (unsigned int nu = 0; nu <=mu; ++nu) {
          if (blockOnGridData->negligible[nu]) continue;
          auto nuAtom = mapping[nu];
          const unsigned int pointStartIndex = _basisFunctionOnGridController->getFirstIndexOfBlock(blockIndex);
          // A common prefactor for the calculations below
          for_spin(potential,densityMatrix){
            Eigen::VectorXd prefac = (mu==nu?1.0:2.0) * densityMatrix_spin(mu,nu) * weights.segment(pointStartIndex,blockSize);
            Eigen::VectorXd pot = potential_spin.segment(pointStartIndex,blockSize).array() * prefac.array();
            gpp[muAtom+(0*nAtoms)] -= pot.dot(Eigen::VectorXd(bfValues.col(nu).array() * bfDerivatives.x.col(mu).array()));
            gpp[muAtom+(1*nAtoms)] -= pot.dot(Eigen::VectorXd(bfValues.col(nu).array() * bfDerivatives.y.col(mu).array()));
            gpp[muAtom+(2*nAtoms)] -= pot.dot(Eigen::VectorXd(bfValues.col(nu).array() * bfDerivatives.z.col(mu).array()));
            gpp[nuAtom+(0*nAtoms)] -= pot.dot(Eigen::VectorXd(bfValues.col(mu).array() * bfDerivatives.x.col(nu).array()));
            gpp[nuAtom+(1*nAtoms)] -= pot.dot(Eigen::VectorXd(bfValues.col(mu).array() * bfDerivatives.y.col(nu).array()));
            gpp[nuAtom+(2*nAtoms)] -= pot.dot(Eigen::VectorXd(bfValues.col(mu).array() * bfDerivatives.z.col(nu).array()));
          };
          if (_functional.getFunctionalClass() == FUNCTIONAL_CLASSES::GGA) {

            const auto& bfSecondDer = *blockOnGridData->secondDerivativeValues;

            // variables for the spin macro
            const auto& vdx = potentialGradients->x;
            const auto& vdy = potentialGradients->y;
            const auto& vdz = potentialGradients->z;
            for_spin(vdx,vdy,vdz,densityMatrix){


              Eigen::VectorXd prefac = (mu==nu?1.0:2.0) * densityMatrix_spin(mu,nu) * weights.segment(pointStartIndex,blockSize);

              Eigen::VectorXd potentialx = vdx_spin.segment(pointStartIndex,blockSize).array() * prefac.array();
              Eigen::VectorXd potentialy = vdy_spin.segment(pointStartIndex,blockSize).array() * prefac.array();
              Eigen::VectorXd potentialz = vdz_spin.segment(pointStartIndex,blockSize).array() * prefac.array();
              Eigen::VectorXd XmuXnu = bfDerivatives.x.col(mu).array() * bfDerivatives.x.col(nu).array();
              Eigen::VectorXd XmuYnu = bfDerivatives.x.col(mu).array() * bfDerivatives.y.col(nu).array();
              Eigen::VectorXd XmuZnu = bfDerivatives.x.col(mu).array() * bfDerivatives.z.col(nu).array();
              Eigen::VectorXd YmuXnu = bfDerivatives.y.col(mu).array() * bfDerivatives.x.col(nu).array();
              Eigen::VectorXd YmuYnu = bfDerivatives.y.col(mu).array() * bfDerivatives.y.col(nu).array();
              Eigen::VectorXd YmuZnu = bfDerivatives.y.col(mu).array() * bfDerivatives.z.col(nu).array();
              Eigen::VectorXd ZmuXnu = bfDerivatives.z.col(mu).array() * bfDerivatives.x.col(nu).array();
              Eigen::VectorXd ZmuYnu = bfDerivatives.z.col(mu).array() * bfDerivatives.y.col(nu).array();
              Eigen::VectorXd ZmuZnu = bfDerivatives.z.col(mu).array() * bfDerivatives.z.col(nu).array();

              gpp[muAtom+(0*nAtoms)] -=  (potentialx.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xx.col(mu).array()))+ XmuXnu))
                                              +(potentialy.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xy.col(mu).array()))+ XmuYnu))
                                              +(potentialz.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xz.col(mu).array()))+ XmuZnu));

              gpp[muAtom+(1*nAtoms)] -=  (potentialx.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xy.col(mu).array()))+ YmuXnu))
                                              +(potentialy.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.yy.col(mu).array()))+ YmuYnu))
                                              +(potentialz.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.yz.col(mu).array()))+ YmuZnu));

              gpp[muAtom+(2*nAtoms)] -=  (potentialx.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xz.col(mu).array()))+ ZmuXnu))
                                              +(potentialy.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.yz.col(mu).array()))+ ZmuYnu))
                                              +(potentialz.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.zz.col(mu).array()))+ ZmuZnu));

              gpp[nuAtom+(0*nAtoms)] -=  (potentialx.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xx.col(nu).array()))+ XmuXnu))
                                              +(potentialy.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xy.col(nu).array()))+ YmuXnu))
                                              +(potentialz.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xz.col(nu).array()))+ ZmuXnu));

              gpp[nuAtom+(1*nAtoms)] -=  (potentialx.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xy.col(nu).array()))+ XmuYnu))
                                              +(potentialy.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.yy.col(nu).array()))+ YmuYnu))
                                              +(potentialz.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.yz.col(nu).array()))+ ZmuYnu));

              gpp[nuAtom+(2*nAtoms)] -=  (potentialx.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xz.col(nu).array()))+ XmuZnu))
                                              +(potentialy.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.yz.col(nu).array()))+ YmuZnu))
                                              +(potentialz.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.zz.col(nu).array()))+ ZmuZnu));
            };
          }
        }

      }
    }
#pragma omp critical
    {
      gradientContr += gradientContrPriv;
    }
  } /* END OpenMP parallel */
  if (_functional.getFunctionalClass() == FUNCTIONAL_CLASSES::GGA) {
    _basisFunctionOnGridController->setHighestDerivative(derivorder);
  }
  return gradientContr;
}


template class NAddFuncPotential<Options::SCF_MODES::RESTRICTED>;
template class NAddFuncPotential<Options::SCF_MODES::UNRESTRICTED>;


template <Options::SCF_MODES SCFMode>
NAddEnergyHelper<SCFMode>::NAddEnergyHelper(Functional func,
    std::shared_ptr<DensityOnGridController<SCFMode> > ctrs):
    _functional(func),
    _densOnGridControllers(ctrs){
  _densOnGridControllers->addSensitiveObject(ObjectSensitiveClass<DensityOnGrid<SCFMode> >::_self);
}

template <Options::SCF_MODES SCFMode>
double NAddEnergyHelper<SCFMode>::getEnergy(){
  if (_old){
    XCFun<SCFMode> xcfun(128);
    _energy = xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        _functional,
        _densOnGridControllers,
        0).energy;
    _old=false;
  }
  return _energy;
}
template class NAddEnergyHelper<Options::SCF_MODES::RESTRICTED>;
template class NAddEnergyHelper<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
