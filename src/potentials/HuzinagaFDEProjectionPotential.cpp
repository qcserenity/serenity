/**
 * @file HuzinagaFDEProjectionPotential.cpp
 *
 * @date Nov 23, 2017
 * @author Moritz Bensberg
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
#include "potentials/HuzinagaFDEProjectionPotential.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/OrbitalController.h"
#include "dft/Functional.h"
#include "dft/functionals/wrappers/XCFun.h"
#include "grid/GridControllerFactory.h"
#include "math/Derivatives.h"
#include "misc/Timing.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/bundles/ESIPotentials.h"
#include "potentials/ZeroPotential.h"
#include "system/SystemController.h"
#include "data/matrices/SPMatrix.h"
#include "potentials/ECPInteractionPotential.h"
#include "data/ElectronicStructure.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "misc/SystemSplittingTools.h"
#include "tasks/LocalizationTask.h"
#include "io/IOOptions.h"

/* ABFockMatrixConstruction */
#include "potentials/ABFockMatrixConstruction/ABCoreHamiltonian.h"
#include "potentials/ABFockMatrixConstruction/ABCoulombInteractionPotential.h"
#include "potentials/ABFockMatrixConstruction/ABFuncPotential.h"
#include "potentials/ABFockMatrixConstruction/ABNAddFuncPotential.h"
#include "potentials/ABFockMatrixConstruction/ABExchangePotential.h"
#include "potentials/ABFockMatrixConstruction/ABHFPotential.h"
#include "potentials/ABFockMatrixConstruction/ABZeroPotential.h"

namespace Serenity {

template <Options::SCF_MODES SCFMode>
HuzinagaFDEProjectionPotential<SCFMode>::HuzinagaFDEProjectionPotential(
    std::shared_ptr<SystemController> activeSystem,
    std::vector<std::shared_ptr<SystemController> > environmentSystems,
    Functional naddXCFunc,
    bool truncatedProjection,
    double projectionOverlapThreshold,
    bool buildHoffmannsOperator,
    bool topDown,
    bool distantKinFunc,
    double basisFunctionRatio,
    double borderAtomThreshold,
    std::shared_ptr<GridController> supersystemgrid,
    Options::KINFUNCTIONALS naddKinFunc,
    double gridCutOff,
    std::vector<std::shared_ptr<EnergyComponentController> > allEConts):
    Potential<SCFMode>(activeSystem->getBasisController()),
    _activeSystem(activeSystem),
    _environmentSystems(environmentSystems),
    _naddXCFunc(naddXCFunc),
    _truncatedProjection(truncatedProjection),
    _projectionOverlapThreshold(projectionOverlapThreshold),
    _buildHoffmannsOperator(buildHoffmannsOperator),
    _topDown(topDown),
    _distantKinFunc(distantKinFunc),
    _supersystemgrid(supersystemgrid){
  takeTime("Huzinaga -- pre-calculations");
  // Expects that an electronic structure is already in place for the active system!
  _activeSystem->getElectronicStructure<SCFMode>()
      ->getDensityMatrixController()
      ->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode> >::_self);
  // the basis controller of the active system.
  auto basisContA = _activeSystem->getBasisController();
  // A projection does only make sense if there is something to project
  assert(_environmentSystems.size() > 0 && "There is no environment system given!");
  /*
   * 1. Calculate S^AB.
   * 2. Determine which systems should be projected.
   * 3. Calculate H^AB, J^AB[env], rho[env] and build grids grid^AB for systems which should be projected
   *    If Hoffmann's operator is supposed to be constructed, save the potential bundle to calculate f_BB.
   */
  unsigned int nEnv = _environmentSystems.size();
  // Default: Project all systems.
  _projectedEnvSystems = Eigen::VectorXi::Constant(nEnv,1);
  /* ==================== */
  /*  1. Calculate S^AB   */
  /* ==================== */
  for (unsigned int i = 0; i < _environmentSystems.size(); ++i) {
    // Container for AB-Aux basis controller
    _abAuxBasisController.push_back(nullptr);
    auto& libint = Libint::getInstance();
    auto s_AB = std::make_shared<Eigen::MatrixXd> (
        libint.compute1eInts(libint2::Operator::overlap, _environmentSystems[i]->getBasisController(),basisContA));
    assert(s_AB->rows()==basisContA->getNBasisFunctions() && s_AB->cols()==_environmentSystems[i]->getBasisController()->getNBasisFunctions());
    /* ================================= */
    /*  2. Determine projected systems   */
    /* ================================= */
    if(_truncatedProjection) {
      /*
       * 1. Evaluate overlap between the basis sets of the subsystems
       * 2. Include all subsystems with a total overlap over a given threshold
       */
      double totalInterSystemBasisOverlap = 0.0;
      for(unsigned int k = 0; k < s_AB->rows(); ++k) {
        for (unsigned int l = 0; l < s_AB->cols(); ++l) {
          totalInterSystemBasisOverlap += std::abs((*s_AB)(k,l));
        }
      }
      assert(totalInterSystemBasisOverlap > -1e-6);
      std::cout << "Total basis overlap ("<<
          _activeSystem->getSystemName()<<" and "<< _environmentSystems[i]->getSystemName() <<"): " <<
          totalInterSystemBasisOverlap << std::endl;
      if (totalInterSystemBasisOverlap <= _projectionOverlapThreshold) {
        _projectedEnvSystems[i] = 0;
        // If the system does not get projected, forget the basis overlap.
        s_AB = nullptr;
      }
    }
    _s_ABs.push_back(s_AB);
  }
  /* ============================================== */
  /*  3. Pre-calculate constant env contributions   */
  /* ============================================== */
  // Build supersystem geometry
  _supersystemGeometry = std::make_shared<Geometry>();
  auto actGeom = _activeSystem->getGeometry();
  *_supersystemGeometry += *actGeom;
  for (unsigned int i = 0; i < _environmentSystems.size(); ++i){
    *_supersystemGeometry += *_environmentSystems[i]->getGeometry();
  }
  _supersystemGeometry->deleteIdenticalAtoms();
  // Build the density matrix controllers of the environment systems in the same
  // spin-polarization as the active system.
  // These are needed for the Coulomb and exchange contribution to the supersystem fock operator of systems,
  // which are not in the system pair.
  // In a top-down calculation, the previouse calculation will always have happened with the right SCFMode and has to be used!
  for (auto sys : _environmentSystems){
    if (sys->getSettings().scfMode == SCFMode or _topDown) {
      assert(sys->hasElectronicStructure<SCFMode>());
      _envDensityCont.push_back(sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController());
    } else {
      if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        //Build unrestricted DensityMatrixController
        DensityMatrix<SCFMode> uDensMat(sys->getBasisController());
        for_spin(uDensMat) {
          uDensMat_spin = 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
        };
        _envDensityCont.push_back(std::make_shared<DensityMatrixController<SCFMode>>(uDensMat));
      } else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        //Build restricted DensityMatrixController
        DensityMatrix<SCFMode> rDensMat(sys->getBasisController());
        for_spin(rDensMat) {
          rDensMat_spin = sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
        };
        _envDensityCont.push_back(std::make_shared<DensityMatrixController<SCFMode>>(rDensMat));
      } else {
        assert(false);
      }
    }
    // If it is a hybrid within each environment system. Construct the not projected
    // density matrix.
    if(_distantKinFunc) {
      // Localization (already done for top-down calculations)
      if(!_topDown) {
        LocalizationTask locTask(sys);
        locTask.run();
      }
      const auto& envCoeff = sys->getActiveOrbitalController<SCFMode>()->getCoefficients();
      auto orbitalPopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
          envCoeff,
          sys->getOneElectronIntegralController()->getOverlapIntegrals(),
          sys->getAtomCenteredBasisController()->getBasisIndices());
      auto distantOrbitals = SystemSplittingTools<SCFMode>::selectDistantOrbitals(
          orbitalPopulations,
          _activeSystem,
          sys,
          basisFunctionRatio,
          borderAtomThreshold);
      auto nonOrthogonalDensityMatrix = SystemSplittingTools<SCFMode>::buildNonOrthogonalDensityMatrix(
          sys,distantOrbitals);
      _notProjectedEnvDensities.push_back(std::make_shared<DensityMatrixController<SCFMode> >(
          *nonOrthogonalDensityMatrix));
    }
  }
  if (_distantKinFunc) {
    assert(supersystemgrid);
    assert(naddKinFunc != Options::KINFUNCTIONALS::NONE);
    // Build non-additive kinetic energy potential for not projected systems
    _naddKinPot = std::make_shared<NAddFuncPotential<SCFMode> >(
        _activeSystem,
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        _notProjectedEnvDensities,
        _supersystemgrid,
        FunctionalClassResolver::resolveFunctional(naddKinFunc),
        std::make_pair(false,(gridCutOff<0)?allEConts:std::vector<std::shared_ptr<EnergyComponentController> >(0)),
        (gridCutOff<0.0)?true:false);
  }
  // Calculate H^AB and J^AB[env]
  std::vector<std::shared_ptr<BasisController> > envAuxBasis;
  if (_activeSystem->getSettings().dft.densityFitting == Options::DENS_FITS::RI) {
    for (const auto& envSys : _environmentSystems) {
      envAuxBasis.push_back(envSys->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));
    }
  } else {
    for (unsigned int i = 0; i < _environmentSystems.size();++i) {
      envAuxBasis.push_back(nullptr);
    }
  }
  for (unsigned int i=0; i < _environmentSystems.size();++i) {
    if (_projectedEnvSystems[i]) {
      _coreAB_Potentials.push_back(
          std::make_shared<ABCoreHamiltonian<SCFMode> >(
              basisContA,
              _environmentSystems[i]->getBasisController(),
              _supersystemGeometry));

      /*
       * Active Coulomb and in the case of HF/Hybrid exchange contribution.
       */
      auto combinedGeometry = std::make_shared<Geometry>();
      auto actGeom = _activeSystem->getGeometry();
      auto envGeom = _environmentSystems[i]->getGeometry();
      *combinedGeometry+= *actGeom;
      *combinedGeometry+= *envGeom;
      combinedGeometry->deleteIdenticalAtoms();
      if(!_abAuxBasisController[i]) {
        _abAuxBasisController[i] = AtomCenteredBasisControllerFactory::produce(
            combinedGeometry,
            _activeSystem->getSettings().basis.basisLibPath,
            _activeSystem->getSettings().basis.makeSphericalBasis,
            false,
            _activeSystem->getSettings().basis.firstECP,
            _activeSystem->getSettings().basis.auxJLabel);
      }
      const auto activeFunctional = FunctionalClassResolver::resolveFunctional(_activeSystem->getSettings().dft.functional);
      double exchangeRatioActive = (_activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF)?
          1.0 : activeFunctional.getHfExchangeRatio();
      std::vector<std::shared_ptr<BasisController> > actAuxVec = {_abAuxBasisController[i]};
      std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > actDensityMatrixController = {
          _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()};
      _activeCoulombAB_Potentials.push_back(
          std::make_shared<ABHFPotential<SCFMode> >(
              _activeSystem,
              basisContA,
              _environmentSystems[i]->getBasisController(),
              actDensityMatrixController,
              exchangeRatioActive,
              _topDown,
              _activeSystem->getSettings().dft.densityFitting,
              _abAuxBasisController[i],
              actAuxVec));
      double exchangeRatioNAdd = _naddXCFunc.getHfExchangeRatio();
    _coulombAB_Potentials.push_back(
        std::make_shared<ABHFPotential<SCFMode> >(
            _environmentSystems[i],
            basisContA,
            _environmentSystems[i]->getBasisController(),
            _envDensityCont,
            exchangeRatioNAdd,
            _topDown,
            _activeSystem->getSettings().dft.densityFitting,
            _abAuxBasisController[i],
            envAuxBasis));
      // Exchange--correlation contributions
      // Grid construction
      //Build A+B grid controller
      std::shared_ptr<GridController> grid_AB;
      grid_AB = GridControllerFactory::produce(combinedGeometry,_activeSystem->getSettings());
      std::shared_ptr<ABPotential<SCFMode> > activeXCPotential;
      if (_activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
        auto actDensityMatrixController = {_activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()};
        activeXCPotential = std::make_shared<ABFuncPotential<SCFMode> > (
            _activeSystem,
            basisContA,
            _environmentSystems[i]->getBasisController(),
            grid_AB,
            actDensityMatrixController,
            activeFunctional);
      } else if (_activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
        activeXCPotential = std::make_shared<ABZeroPotential<SCFMode> > (
            basisContA,
            _environmentSystems[i]->getBasisController());
      } else {
        assert(false && "Unknwon electronic structure theory!");
      }
      _exchangeAB_Potentials.push_back(activeXCPotential);
      _naddExchangeAB_Potentials.push_back(std::make_shared<ABNAddFuncPotential<SCFMode> > (
          _activeSystem,
          basisContA,
          _environmentSystems[i]->getBasisController(),
          _envDensityCont,
          grid_AB,
          _naddXCFunc));
      // build non-additive kinetic energy contribution in case of an hybrid method
      if(_distantKinFunc) {
        _naddKinAB_Potentials.push_back(std::make_shared<ABNAddFuncPotential<SCFMode> >(
            _activeSystem,
            basisContA,
            _environmentSystems[i]->getBasisController(),
            _notProjectedEnvDensities,
            grid_AB,
            FunctionalClassResolver::resolveFunctional(naddKinFunc)));
      }
      if (_buildHoffmannsOperator) {
        // Calculate the core fock matrix for each projected environment system.
        std::shared_ptr<PotentialBundle<SCFMode> > potBundle;
        if (_environmentSystems[i]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
          potBundle=_environmentSystems[i]->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
        } else if (_environmentSystems[i]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
          potBundle=_environmentSystems[i]->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
        } else {
          assert(false && "Unknwon electronic structure theory!");
        }
        std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > densMatControllers;
        std::vector<std::shared_ptr<SystemController> > otherSystemControllers;
        std::vector<std::shared_ptr<const Geometry> > otherGeometries;
        densMatControllers.push_back(_activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController());
        otherSystemControllers.push_back(_activeSystem);
        otherGeometries.push_back(_activeSystem->getGeometry());
        for(unsigned int j = 0; j < _environmentSystems.size();++j) {
          if (i!=j){
            densMatControllers.push_back(_envDensityCont[j]);
            otherSystemControllers.push_back(_environmentSystems[j]);
            otherGeometries.push_back(_environmentSystems[j]->getGeometry());
          }
        }
        // build coulomb interaction potentials
        auto esiPot = std::make_shared<ESIPotentials<SCFMode> >(
            _environmentSystems[i],
            otherSystemControllers,
            _envDensityCont[i],
            _environmentSystems[i]->getGeometry(),
            densMatControllers,
            otherGeometries
        );
        // build non-additive exchange--correlation potential
        auto nonAddExchange = std::make_shared<NAddFuncPotential<SCFMode> >(
            _environmentSystems[i],
            _envDensityCont[i],
            densMatControllers,
            _environmentSystems[i]->getGridController(),
            _naddXCFunc);
        // ECP interaction
        std::vector<std::shared_ptr<Atom> > envAtoms;
        for (const auto& atom : _activeSystem->getGeometry()->getAtoms()) {
          envAtoms.push_back(atom);
        }
        for (unsigned int j = 0; j < _environmentSystems.size(); ++j) {
          if (i != j) {
            for (const auto& atom : _environmentSystems[j]->getGeometry()->getAtoms()) {
              envAtoms.push_back(atom);
            }
          }
        }
        std::shared_ptr<Potential<SCFMode> > ecpInt
        (new ECPInteractionPotential<SCFMode>(
            _environmentSystems[i],
            _environmentSystems[i]->getGeometry()->getAtoms(),
            envAtoms,
            densMatControllers,
            _environmentSystems[i]->getBasisController()));
        // Save it all in a bundle
        _env_i_fdePot.push_back(std::make_shared<FDEPotentials<SCFMode> >(
            potBundle,
            esiPot,
            nonAddExchange,
            std::make_shared<ZeroPotential<SCFMode> >(
                _environmentSystems[i]->getBasisController()),
                ecpInt));
      } // if _buildHoffmannsOperator
    } else {
      _coreAB_Potentials.push_back(nullptr);
      _coulombAB_Potentials.push_back(nullptr);
      _activeCoulombAB_Potentials.push_back(nullptr);
      _exchangeAB_Potentials.push_back(nullptr);
      _naddExchangeAB_Potentials.push_back(nullptr);
      if (_buildHoffmannsOperator){
        _env_i_fdePot.push_back(nullptr);
      }
    }
  }
  if(iOOptions.printDebugInfos) {
    writeInterSubsystemOccOverlap();
  }
  timeTaken(2,"Huzinaga -- pre-calculations");
}

template <Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& HuzinagaFDEProjectionPotential<SCFMode>::getMatrix(){
  if (!_potential) {
    if(iOOptions.printDebugInfos) {
      writeInterSubsystemOccOverlap();
    }
    // reset the fock matrix
    _potential.reset(new FockMatrix<SCFMode> (this->_basis));
    auto& pot = *_potential;
    for_spin(pot) {
      pot_spin.setZero();
    };
    // Loop over all environmental systems.
    for(unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
      if (_projectedEnvSystems[iEnv]) {
        auto environmentSystem = _environmentSystems[iEnv];
        /* ======================================================= *
         *   1. Build outer diagonal fock and get overlap matrix   *
         * ======================================================= */
        // build the outer diagonal fock matrix f_AB
        const SPMatrix<SCFMode> f_AB= buildOuterDiagonalFockMatrix(iEnv); // N(act) x N(env)
        // get the outer diagonal overlap matrix s_AB
        const Eigen::MatrixXd& s_AB = *_s_ABs[iEnv]; // N(act) x N(env)
        /* ========================== *
         *   2. Calculate projector   *
         * ========================== */
        DensityMatrix<SCFMode> projectedMatrix = _envDensityCont[iEnv]->getDensityMatrix();
        if(_naddKinPot) {
          const auto& nonOrthogonalDensityMatrix = _notProjectedEnvDensities[iEnv]->getDensityMatrix();
          projectedMatrix -= nonOrthogonalDensityMatrix;
        }
        double preFactor = (SCFMode == Options::SCF_MODES::RESTRICTED)? 0.5 : 1.0;
        for_spin(f_AB,pot,projectedMatrix) {
          auto matrix =  f_AB_spin * projectedMatrix_spin * s_AB.transpose();
          pot_spin += - preFactor * (matrix+matrix.transpose());
        };
        if (_buildHoffmannsOperator) {
          const FockMatrix<SCFMode> f_BB = _env_i_fdePot[iEnv]->getFockMatrix(
              _envDensityCont[iEnv]->getDensityMatrix(),
              std::make_shared<EnergyComponentController>());
          for_spin(f_BB,pot,projectedMatrix) {
            pot_spin += preFactor * preFactor
                * (s_AB * projectedMatrix_spin * f_BB_spin * projectedMatrix_spin * s_AB.transpose());
          };
        }
      } // if _projectedEnvSystems[iEnv]
    } // for iEnv
    if (_naddKinPot) {
      *_potential+= _naddKinPot->getMatrix();
    }
    if(iOOptions.printDebugInfos) {
      writeInterSubsystemOccOverlap();
    }
  } // if !_potential
  return *_potential;
}

template <Options::SCF_MODES SCFMode>
double HuzinagaFDEProjectionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  double energy = 0.0;
  if (!_potential) getMatrix();
  if (_naddKinPot) energy += _naddKinPot->getEnergy(P);
  return energy;
}

template <Options::SCF_MODES SCFMode>
Eigen::MatrixXd HuzinagaFDEProjectionPotential<SCFMode>::getGeomGradients(){
  Eigen::MatrixXd gradientContr(1,3);
  gradientContr.setZero();

  return gradientContr;
}

template <Options::SCF_MODES SCFMode>
SPMatrix<SCFMode> HuzinagaFDEProjectionPotential<SCFMode>::buildOuterDiagonalFockMatrix(
    unsigned int iEnv) {
  takeTime("Huzinaga -- Build core potential");
  const unsigned int nBasisAct = _activeSystem->getBasisController()->getNBasisFunctions();
  const unsigned int nBasisEnv = _environmentSystems[iEnv]->getBasisController()->getNBasisFunctions();
  // initialize the Fock matrix
  SPMatrix<SCFMode> f_AB(Eigen::MatrixXd::Zero(nBasisAct,nBasisEnv));

  // 1. Add the core hamiltonian
  f_AB += _coreAB_Potentials[iEnv]->getMatrix();
  // 2. Add the exchange--correlation potential
  f_AB += _exchangeAB_Potentials[iEnv]->getMatrix();
  f_AB += _naddExchangeAB_Potentials[iEnv]->getMatrix();

  // 3. Add the coulomb contribution
  f_AB += _activeCoulombAB_Potentials[iEnv]->getMatrix();
  f_AB += _coulombAB_Potentials[iEnv]->getMatrix();
  if(_naddKinPot) f_AB += _naddKinAB_Potentials[iEnv]->getMatrix();
  timeTaken(2,"Huzinaga -- Build core potential");
  return f_AB;
}

template <Options::SCF_MODES SCFMode>
void HuzinagaFDEProjectionPotential<SCFMode>::writeInterSubsystemOccOverlap() {

  auto nOccAct = _activeSystem->getNOccupiedOrbitals<SCFMode>();
  const auto& coeffAct = _activeSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  for (unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
    if(_projectedEnvSystems[iEnv]){
      const auto& coeffEnv = _environmentSystems[iEnv]->getActiveOrbitalController<SCFMode>()->getCoefficients();
      auto nOccEnv = _environmentSystems[iEnv]->getNOccupiedOrbitals<SCFMode>();
      assert(_s_ABs[iEnv]);
      const auto& s_AB = *_s_ABs[iEnv];
      // Calculate the average overlap between the occupied orbitals of the subsystems.
      for_spin(coeffAct,coeffEnv,nOccAct,nOccEnv) {
        double totalInterSystemOccOverlap = 0.0;
        Eigen::MatrixXd totalOverlap = (
            coeffAct_spin.leftCols(nOccAct_spin).transpose()*s_AB*coeffEnv_spin.leftCols(nOccEnv_spin)).array().abs().matrix();
        assert(totalOverlap.rows()==nOccAct_spin&&totalOverlap.cols()==nOccEnv_spin);
        for(unsigned int i = 0; i < nOccAct_spin; ++i) {
          for (unsigned int j = 0; j < nOccEnv_spin; ++j) {
            double overlap =totalOverlap(i,j);
            totalInterSystemOccOverlap += overlap;
          }
        }
        unsigned int loopEnd = (nOccAct_spin*nOccEnv_spin >= 10)? 10 : nOccAct_spin*nOccEnv_spin;
        std::cout << loopEnd<<" largest overlaps "<< _activeSystem->getSystemName()<< " "
            <<_environmentSystems[iEnv]->getSystemName() << std::endl;
        for(unsigned int i=0; i < loopEnd; ++i) {
          int indexAct=0;
          int indexEnv=0;
          totalOverlap.maxCoeff(&indexAct,&indexEnv);
          std::cout << "Act "<<indexAct<<"<->"<<indexEnv<<"  "<<totalOverlap(indexAct,indexEnv)<<std::endl;
          totalOverlap(indexAct,indexEnv) = 0.0;
        }
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "Total intersystem overlap of occ. orbitals("<<_activeSystem->getSystemName() <<" and "
            << _environmentSystems[iEnv]->getSystemName() << "): " << totalInterSystemOccOverlap <<std::endl;
      };
    }
  }
}

template class HuzinagaFDEProjectionPotential<Options::SCF_MODES::RESTRICTED>;
template class HuzinagaFDEProjectionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
