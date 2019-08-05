/**
 * @file TDReconstructionPotential.cpp
 *
 * @date Apr 03, 2018
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

#include "potentials/TDReconstructionPotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/CoulombPotentialOnGridCalculator.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/ElectronicStructure.h"
#include "input/FunctionalClassResolver.h"
#include "grid/GridController.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/OneIntControllerFactory.h"
#include "potentials/OptEffPotential.h"
#include "potentials/bundles/PotentialBundle.h"
#include "potentials/HFPotential.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "dft/functionals/wrappers/XCFun.h"
#include "energies/EnergyContributions.h"
#include "potentials/ZeroPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/bundles/ESIPotentials.h"
#include "potentials/ECPInteractionPotential.h"
#include "potentials/bundles/PBEPotentials.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
TDReconstructionPotential<SCFMode>::TDReconstructionPotential(
    std::shared_ptr<SystemController> actSys,
    std::shared_ptr<SystemController> supSys,
    std::vector<std::shared_ptr<SystemController>> envSystems,
    double smoothFactor,
    std::string potBasisLabel,
    const double singValThreshold,
    double lbDamping,
    unsigned int lbCycles,
    unsigned int carterCycles,
    bool noSupRec):
    Potential<SCFMode>(actSys->getBasisController()),
    _potential(nullptr),
    _actSys(actSys),
    _supSys(supSys),
    _envSystems(envSystems),
    _smoothFactor(smoothFactor),
    _potBasisLabel(potBasisLabel),
    _singValThreshold(singValThreshold),
    _lbDamping(lbDamping),
    _lbCycles(lbCycles),
    _carterCycles(carterCycles),
    _noSupRec(noSupRec),
    _supXEnergy(0){

}

template<Options::SCF_MODES SCFMode>
void TDReconstructionPotential<SCFMode>::calculatePotential(){
  /*
   * Prepare output data
   */
  this->_potential.reset(new FockMatrix<SCFMode>(this->_basis));
  auto& pot=*this->_potential;

  /*
   * LibInt
   */
  auto& libint = Libint::getInstance();

  /*
   * Get grid for reconstruction
   */
  auto superSystemGrid=_supSys->getGridController();

  /*
   * Active system objects
   */
  auto actNOccOrbs=_actSys->template getNOccupiedOrbitals<SCFMode>();
  const auto& actNEl=_actSys->template getNElectrons<SCFMode>();
  auto actSysBasisController=_actSys->getBasisController();
  auto actSysGeom=_actSys->getGeometry();
  auto actBasFuncOnGridController=BasisFunctionOnGridControllerFactory::produce(
      128,0.0,2,actSysBasisController,superSystemGrid);
  auto actOneEIntController=OneIntControllerFactory::getInstance().produce(_actSys->getBasisController(),actSysGeom);
  auto actDensOnGridCalc=std::make_shared<DensityOnGridCalculator<SCFMode> >(
      actBasFuncOnGridController, 0.0);
  const auto& actSysDensMat=_actSys->template getElectronicStructure<SCFMode>()->getDensityMatrix();
  auto actDensOnGrid=actDensOnGridCalc->calcDensityOnGrid(actSysDensMat);
  auto orbsAct=_actSys->template getElectronicStructure<SCFMode>()->getMolecularOrbitals();

  /*
   * Supersystem objects
   */
  auto supNOccOrbs=_supSys->template getNOccupiedOrbitals<SCFMode>();
  const auto& supNEl=_supSys->template getNElectrons<SCFMode>();
  auto supSysBasisController=_supSys->getBasisController();
  auto supSysGeom=_supSys->getGeometry();
  auto supBasFuncOnGridController=BasisFunctionOnGridControllerFactory::produce(
      128,0.0,2,supSysBasisController,superSystemGrid);
  auto supOneEIntController=OneIntControllerFactory::getInstance().produce(_supSys->getBasisController(),supSysGeom);
  auto supDensOnGridCalc=std::make_shared<DensityOnGridCalculator<SCFMode> >(
      supBasFuncOnGridController, 0.0);
  const auto& supSysDensMat=_supSys->template getElectronicStructure<SCFMode>()->getDensityMatrix();
  auto supDensOnGrid=supDensOnGridCalc->calcDensityOnGrid(supSysDensMat);


  /*
   * Hybrid functional?
   */
  auto functional=FunctionalClassResolver::resolveFunctional(_supSys->getSettings().dft.functional);
  double exc=_carterCycles!=0? -1.0 : functional.isHybrid()? functional.getHfExchangeRatio() : -1.0;

  /*
   * Calculate energy from supersystem orbitals
   */
  auto supKin=libint.compute1eInts(libint2::Operator::kinetic, supSysBasisController);
  _supEnergy=0.0;
  for_spin(supSysDensMat){
    _supEnergy+=supSysDensMat_spin.cwiseProduct(supKin).sum();
  };


  /*
   * Get basis to express the potential in
   */

  /*
   * Generate geometry without dummyAtoms, since supersystem
   * potentialbases introduce instable potential reconstructions
   */

  std::vector<std::shared_ptr<Atom>> onlyActAtoms;
  for(auto atom : actSysGeom->getAtoms()){
    if(!atom->isDummy())onlyActAtoms.push_back(atom);
  }

  auto onlyActGeom=std::make_shared<Geometry>(onlyActAtoms);

  auto actPotBasisController=AtomCenteredBasisControllerFactory::produce(
      onlyActGeom,
      _actSys->getSettings().basis.basisLibPath,
      _actSys->getSettings().basis.makeSphericalBasis,
      false,
      _actSys->getSettings().basis.firstECP,
      _potBasisLabel);

  auto actPotBasFuncOnGridController=BasisFunctionOnGridControllerFactory::produce(
      128,0.0,2,actPotBasisController,superSystemGrid);

  auto supPotBasisController=AtomCenteredBasisControllerFactory::produce(
      supSysGeom,
      _supSys->getSettings().basis.basisLibPath,
      _supSys->getSettings().basis.makeSphericalBasis,
      false,
      _supSys->getSettings().basis.firstECP,
      _potBasisLabel);

  auto supPotBasFuncOnGridController=BasisFunctionOnGridControllerFactory::produce(
      128,0.0,2,supPotBasisController,superSystemGrid);

  OptEffPotential<SCFMode> oepCalcAct(actBasFuncOnGridController,actPotBasFuncOnGridController,actOneEIntController,actDensOnGridCalc,actNOccOrbs,_smoothFactor,_singValThreshold,exc);
  OptEffPotential<SCFMode> oepCalc(supBasFuncOnGridController,supPotBasFuncOnGridController,supOneEIntController,supDensOnGridCalc,supNOccOrbs,_smoothFactor,_singValThreshold,exc);

  if(_carterCycles==0){
    /*
     * Run Wu-Yang and van Leeuwen-Baerends
     */

    /*
     * Active system initial guess: Fermi-Amaldi potential + Nuclear potential
     */
    CoulombPotentialOnGridCalculator coulOnGridCalcAct(actBasFuncOnGridController);
    GridPotential<SCFMode> initialGuessAct(superSystemGrid);
    GridPotential<SCFMode> negCoulOnGrid(superSystemGrid);
    GridPotential<RESTRICTED> tmp(superSystemGrid);
    coulOnGridCalcAct.calculateElectronElectron<SCFMode>(tmp,actSysDensMat);
    for_spin(initialGuessAct,actNEl){
      initialGuessAct_spin = tmp;
      initialGuessAct_spin *= (1.0-(1.0/actNEl_spin));
    };
    tmp.setZero();
    coulOnGridCalcAct.calculateElectronNuclei(tmp,actSysGeom->getAtoms());
    for_spin(initialGuessAct){
      initialGuessAct_spin += tmp;
    };

    /*
     * Supersystem initial guess: Fermi-Amaldi potential + Nuclear potential
     */
    CoulombPotentialOnGridCalculator coulOnGridCalc(supBasFuncOnGridController);
    GridPotential<SCFMode> initialGuess(superSystemGrid);
    tmp.setZero();
    coulOnGridCalc.calculateElectronElectron<SCFMode>(tmp,supSysDensMat);
    for_spin(initialGuess,supNEl,negCoulOnGrid){
      initialGuess_spin += tmp;
      negCoulOnGrid_spin -= tmp;
      initialGuess_spin *= (1.0-(1.0/supNEl_spin));
    };
    tmp.setZero();
    coulOnGridCalc.calculateElectronNuclei(tmp,supSysGeom->getAtoms());
    for_spin(initialGuess,negCoulOnGrid){
      negCoulOnGrid_spin -= tmp;
      initialGuess_spin += tmp;
    };

    /*
     * Perform potential reconstruction to obtain missing parts of the embedding potential
     */
    std::cout << "Active System" << std::endl;

    /*
     * First Wu-Yang
     */
    std::cout << "Running Wu-Yang potential reconstruction" << std::endl;
    oepCalcAct.calculateOEP(actDensOnGrid,orbsAct,initialGuessAct);
    initialGuessAct+=oepCalcAct.getPotentialOnGrid();


    /*
     * Then van Leeuwen--Baerends
     */
    if(_lbCycles!=0){
      std::cout << "Running van Leeuwen-Baerends potential reconstruction" << std::endl;
      oepCalcAct.calculateOEPLB(actDensOnGrid,orbsAct,initialGuessAct,_lbDamping,_lbCycles);
      initialGuessAct+=oepCalcAct.getPotentialOnGrid();
    }


    if(!_noSupRec){
      /*
       * Get the supersystem potential via potential reconstruction
       */
      std::cout << "Supersystem" << std::endl;
      auto dummyOrbs=std::make_shared<OrbitalController<SCFMode>>(supSysBasisController);
      /*
       * First Wu-Yang
       */
      std::cout << "Running Wu-Yang potential reconstruction" << std::endl;
      oepCalc.calculateOEP(supDensOnGrid,dummyOrbs,initialGuess);
      initialGuess+=oepCalc.getPotentialOnGrid();
      /*
       * Then van Leeuwen--Baerends
       */
      if(_lbCycles!=0){
        std::cout << "Running van Leeuwen-Baerends potential reconstruction" << std::endl;
        oepCalc.calculateOEPLB(supDensOnGrid,dummyOrbs,initialGuess,_lbDamping,_lbCycles);
        initialGuess+=oepCalc.getPotentialOnGrid();
      }

    }


    /*
     * Construct non-additive kin pot as difference between the
     * supersystem potential and every other part of the total
     * potential
     */
    ScalarOperatorToMatrixAdder<SCFMode> scalarOpToMat(actBasFuncOnGridController,0.0);

    GridPotential<SCFMode> potOnGrid(superSystemGrid);

    if(!_noSupRec)potOnGrid-=initialGuess;
    potOnGrid+=initialGuessAct;
    /*
     * Get the matrix representation
     */
    scalarOpToMat.addScalarOperatorToMatrix(pot, potOnGrid);


    if(_noSupRec){
      scalarOpToMat.addScalarOperatorToMatrix(pot, negCoulOnGrid);

      auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<SCFMode> >(
          supDensOnGridCalc, _supSys->template getElectronicStructure<SCFMode>()->getDensityMatrixController());

      XCFun<SCFMode> xcFun(128);
      MatrixInBasis<SCFMode> tmp(this->_basis);
      auto funcData = xcFun.calcData(
          FUNCTIONAL_DATA_TYPE::GRADIENTS,
          functional,
          densOnGridController);

      if (functional.getFunctionalClass() == FUNCTIONAL_CLASSES::LDA){
        scalarOpToMat.addScalarOperatorToMatrix(tmp,*funcData.dFdRho);
      } else if (functional.getFunctionalClass() == FUNCTIONAL_CLASSES::GGA) {
        scalarOpToMat.addScalarOperatorToMatrix(tmp,*funcData.dFdRho,*funcData.dFdGradRho);
      } else {
        assert(false && "Unsupported functional type as functional!");
      }
      pot-=tmp;
    }

  }else{
    /*
     * Run Zhang-Carter reconstruction
     */

    std::cout << "Active system" << std::endl;
    /*
     * Copy of actSysDensMat to use as target density (actSysDensMat will get
     * altered during the reconstruction procedure)
     */
    DensityMatrix<SCFMode> targetDens(actSysDensMat);

    /*
     * Prepare initial guess: The use potential that is requested for the active system
     */
    FockMatrix<SCFMode> initialGuessMat(actSysBasisController);

    std::shared_ptr<PotentialBundle<SCFMode> > potentials;
    if (_supSys->getSettings().method==Options::ELECTRONIC_STRUCTURE_THEORIES::HF){
      potentials = _actSys->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
    } else if (_supSys->getSettings().method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
      potentials = _actSys->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(
          Options::GRID_PURPOSES::DEFAULT);
    } else {
      std::cout << "ERROR: None existing electronicStructureTheory requested." << std::endl;
      assert(false);
    }

    std::shared_ptr<PotentialBundle<SCFMode> > esiPot;
    if (_supSys->getSettings().dft.densityFitting!=Options::DENS_FITS::RI){
      esiPot = std::shared_ptr<PotentialBundle<SCFMode> >(
          new ESIPotentials<SCFMode>(
              _actSys,
              {_envSystems},
              _actSys->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
              _actSys->getGeometry(),
              {_envSystems[0]->getElectronicStructure<SCFMode>()->getDensityMatrixController()},
              {_envSystems[0]->getGeometry()}));
    } else {
      std::vector<std::shared_ptr<BasisController> > envAuxBasis;

      esiPot = std::shared_ptr<PotentialBundle<SCFMode> >(
          new ESIPotentials<SCFMode>(
              _actSys,
              {_envSystems},
              _actSys->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
              _actSys->getGeometry(),
              {_envSystems[0]->getElectronicStructure<SCFMode>()->getDensityMatrixController()},
              {_envSystems[0]->getGeometry()},
              _actSys->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB),
              {_envSystems[0]->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)}));
    }
    // ECP TODO: Check consistency!
    std::shared_ptr<Potential<SCFMode> > ecpInt
      (new ECPInteractionPotential<SCFMode>(
          _actSys,
          _actSys->getGeometry()->getAtoms(),
          _envSystems[0]->getGeometry()->getAtoms(),
          {_envSystems[0]->getElectronicStructure<SCFMode>()->getDensityMatrixController()},
          _actSys->getBasisController()));

    auto naddXCPot = std::shared_ptr<NAddFuncPotential<SCFMode> >(new NAddFuncPotential<SCFMode>(
        _actSys,
        _actSys->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        {_envSystems[0]->template getElectronicStructure<SCFMode>()->getDensityMatrixController()},
        superSystemGrid,
        functional,
        {true , std::vector<std::shared_ptr<EnergyComponentController> >(0)},
        true,
        false));

    auto zero = std::shared_ptr<Potential<SCFMode> >(
          new ZeroPotential<SCFMode>(_actSys->getBasisController()));

    /*
     * Bundle potentials
     */
     auto guessPot = std::shared_ptr<PotentialBundle<SCFMode> >(
         new PBEPotentials<SCFMode>(
             potentials,
             esiPot,
             naddXCPot,
             zero,
             ecpInt));

     auto eController=std::make_shared<EnergyComponentController>();

    initialGuessMat+=guessPot->getFockMatrix(actSysDensMat,eController);

    /*
     * Subtract kinetic part
     */
     auto actKin=libint.compute1eInts(libint2::Operator::kinetic, actSysBasisController);
     for_spin(initialGuessMat){
       initialGuessMat_spin-=actKin;
     };

     std::cout << "Running Carter potential reconstruction" << std::endl;
     oepCalcAct.calculateOEPCarter(targetDens,orbsAct,initialGuessMat,_carterCycles);
     pot+=oepCalcAct.getMatrix();

  }

  /*
   * Exact exchange part
   */
  if(functional.isHybrid()){

    auto eCont=_actSys->template getElectronicStructure<SCFMode>()->getEnergyComponentController();

    _supXEnergy=0.0;

    double thresh = _supSys->getSettings().basis.integralThreshold;
    HFPotential<SCFMode> KSup(_supSys,
        _supSys->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        functional.getHfExchangeRatio(),
        thresh);
    _supXEnergy+=KSup.getXEnergy(supSysDensMat);

    if(_noSupRec and _carterCycles==0){
      pot+=KSup.getXPotential();
    }

  }

}

template <Options::SCF_MODES SCFMode>
double TDReconstructionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P){
  if(!this->_potential){
    calculatePotential();
  };
  double energy=_supEnergy;
  auto& libint = Libint::getInstance();
  auto kin=libint.compute1eInts(libint2::Operator::kinetic, _actSys->getBasisController());
  for_spin(P){
    double kinE=P_spin.cwiseProduct(kin).sum();
    energy-=kinE;
  };

  for(auto sys : this->_envSystems){
    auto kin=libint.compute1eInts(libint2::Operator::kinetic, sys->getBasisController());
    auto mat=sys->template getElectronicStructure<SCFMode>()->getDensityMatrix();
    for_spin(mat){
      double kinE=mat_spin.cwiseProduct(kin).sum();
      energy-=kinE;
    };
  }

  if(_supXEnergy!=0.0){
    double xEnergy=_supXEnergy;
    xEnergy-=_actSys->getElectronicStructure<SCFMode>()->getEnergyComponentController()->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXACT_EXCHANGE);
    for(auto sys : this->_envSystems){
      xEnergy-=sys->template getElectronicStructure<SCFMode>()->getEnergyComponentController()->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXACT_EXCHANGE);
    }
    _actSys->getElectronicStructure<SCFMode>()->getEnergyComponentController()->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_EXACT_EXCHANGE,xEnergy);
  }


  return energy;
};

template <Options::SCF_MODES SCFMode> Eigen::MatrixXd
TDReconstructionPotential<SCFMode>::getGeomGradients(){

  throw SerenityError(
      (string)"Geometrical gradients not available for potential reconstruction methods!");

  Eigen::MatrixXd gradientContr(1,3);
  gradientContr.setZero();

  return gradientContr;
}

template
class TDReconstructionPotential<Options::SCF_MODES::RESTRICTED>;
template
class TDReconstructionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
