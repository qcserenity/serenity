/**
 * @file BUReconstructionPotential.cpp
 *
 * @date Oct 13, 2016
 * @author David Schnieders
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "potentials/BUReconstructionPotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/CoulombPotentialOnGridCalculator.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/grid/SupersystemDensityOnGridController.h"
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h"
#include "integrals/OneIntControllerFactory.h"
#include "integrals/wrappers/Libint.h"
#include "potentials/ERIPotential.h"
#include "potentials/OptEffPotential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
BUReconstructionPotential<SCFMode>::BUReconstructionPotential(
    std::shared_ptr<SystemController> actSys, std::shared_ptr<DensityMatrixController<SCFMode>> activeDMat,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDMats,
    std::shared_ptr<GridController> superSystemGrid, Functional functional,
    std::vector<std::shared_ptr<SystemController>> envSystems, double smoothFactor, std::string potBasisLabel,
    const double singValThreshold, double lbDamping, unsigned int lbCycles, unsigned int carterCycles)
  : NAddFuncPotential<SCFMode>(actSys, activeDMat, envDMats, superSystemGrid, functional),
    _basis(actSys->getBasisController()),
    _smoothFactor(smoothFactor),
    _potBasisLabel(potBasisLabel),
    _singValThreshold(singValThreshold),
    _lbDamping(lbDamping),
    _lbCycles(lbCycles),
    _carterCycles(carterCycles),
    _potentialOnGrid(nullptr) {
  for (auto e : envSystems)
    _envSystems.push_back(e);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& BUReconstructionPotential<SCFMode>::getMatrix() {
  if (!this->_potential) {
    calculatePotential();
  };
  return *this->_potential;
};

template<Options::SCF_MODES SCFMode>
double BUReconstructionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (!this->_potential) {
    calculatePotential();
  };
  double energy = _supEnergy;
  auto& libint = Libint::getInstance();
  auto system = this->_system.lock();
  auto kin = libint.compute1eInts(LIBINT_OPERATOR::kinetic, system->getBasisController());
  for_spin(P) {
    double kinE = P_spin.cwiseProduct(kin).sum();
    energy -= kinE;
  };

  for (auto weak : this->_envSystems) {
    auto sys = weak.lock();
    auto kin = libint.compute1eInts(LIBINT_OPERATOR::kinetic, sys->getBasisController());
    auto mat = sys->template getElectronicStructure<SCFMode>()->getDensityMatrix();
    for_spin(mat) {
      double kinE = mat_spin.cwiseProduct(kin).sum();
      energy -= kinE;
    };
  }

  return energy;
}

template<Options::SCF_MODES SCFMode>
void BUReconstructionPotential<SCFMode>::calculatePotential() {
  /*
   * Prepare output data
   */
  this->_potential.reset(new FockMatrix<SCFMode>(_basis));
  this->_potentialOnGrid.reset(new GridPotential<SCFMode>(this->_grid));
  auto& pot = *this->_potential;

  /*
   * LibInt
   */
  auto& libint = Libint::getInstance();

  /*
   * Supersystem geometry and nOccupiedOrbitals
   */
  auto system = this->_system.lock();
  auto supSysGeom = std::make_shared<Geometry>();
  *supSysGeom += *(system->getGeometry());
  auto supNOccOrbs = system->template getNOccupiedOrbitals<SCFMode>();
  auto supNEl = system->template getNElectrons<SCFMode>();
  for (auto weak : this->_envSystems) {
    auto sys = weak.lock();
    *supSysGeom += *(sys->getGeometry());
    supNOccOrbs += sys->template getNOccupiedOrbitals<SCFMode>();
    supNEl += sys->template getNElectrons<SCFMode>();
  }

  /*
   * Active system objects
   */
  auto actNEl = system->template getNElectrons<SCFMode>();
  auto actNOccOrbs = system->template getNOccupiedOrbitals<SCFMode>();
  auto actSysBasisController = system->getBasisController();
  auto actSysGeom = system->getGeometry();

  auto actBasFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(system->getSettings().grid.blocksize, 0.0,
                                                                                  2, actSysBasisController, this->_grid);
  auto oneEIntControllerAct = system->getOneElectronIntegralController();
  auto actDensOnGridCalc = std::make_shared<DensityOnGridCalculator<SCFMode>>(actBasFuncOnGridController, 0.0);
  auto dummyOrbsAct = std::make_shared<OrbitalController<SCFMode>>(actSysBasisController, 0);

  /*
   * Supersystem basis
   */
  std::shared_ptr<BasisController> supSysBasisController = AtomCenteredBasisControllerFactory::produce(
      supSysGeom, system->getSettings().basis.basisLibPath, system->getSettings().basis.makeSphericalBasis, false,
      system->getSettings().basis.firstECP, system->getSettings().basis.label);

  /*
   * If subsystem is expressed in supersystem basis: Set supersystem basis to subsystem basis
   * (to have the correct order of basis functions)
   */
  bool superSystemBasis = supSysBasisController->getNBasisFunctions() == actSysBasisController->getNBasisFunctions();
  if (superSystemBasis)
    supSysBasisController = actSysBasisController;

  /*
   * Supersystem density matrix (block diagonal if subsystem bases are used)
   */
  DensityMatrix<SCFMode> supSysDensMat(supSysBasisController);
  auto actDMat = this->_actDMatController->getDensityMatrix();
  unsigned int blockStart = 0;
  for_spin(supSysDensMat, actDMat) {
    supSysDensMat_spin.block(0, 0, actDMat_spin.rows(), actDMat_spin.cols()) += actDMat_spin;
  };
  if (!superSystemBasis)
    blockStart += actSysBasisController->getNBasisFunctions();
  for (auto dmat : this->_envDMatController) {
    auto envDMat = dmat->getDensityMatrix();
    for_spin(supSysDensMat, envDMat) {
      supSysDensMat_spin.block(blockStart, blockStart, envDMat_spin.rows(), envDMat_spin.cols()) += envDMat_spin;
    };
    if (!superSystemBasis)
      blockStart += envDMat.getBasisController()->getNBasisFunctions();
  }

  /*
   * Supersystem objects
   */
  auto oneEIntController = OneIntControllerFactory::getInstance().produce(supSysBasisController, supSysGeom,
                                                                          system->getExternalChargeController());
  auto supBasFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(system->getSettings().grid.blocksize, 0.0,
                                                                                  2, supSysBasisController, this->_grid);
  auto supDensOnGridCalc = std::make_shared<DensityOnGridCalculator<SCFMode>>(supBasFuncOnGridController, 0.0);
  auto dummyOrbs = std::make_shared<OrbitalController<SCFMode>>(supSysBasisController, system->getNCoreElectrons() / 2);

  /*
   * Hybrid functional?
   */
  auto envSystemZero = _envSystems[0].lock();
  auto functional = envSystemZero->getSettings().customFunc.basicFunctionals.size()
                        ? Functional(envSystemZero->getSettings().customFunc)
                        : resolveFunctional(envSystemZero->getSettings().dft.functional);
  double exc = _carterCycles != 0 ? -1.0 : functional.isHybrid() ? functional.getHfExchangeRatio() : -1.0;

  /*
   * Get basis to express the potential in
   */
  auto actPotBasisController = AtomCenteredBasisControllerFactory::produce(
      actSysGeom, system->getSettings().basis.basisLibPath, system->getSettings().basis.makeSphericalBasis, false,
      system->getSettings().basis.firstECP, _potBasisLabel);

  auto actPotBasFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
      system->getSettings().grid.blocksize, 0.0, 2, actPotBasisController, this->_grid);

  auto supPotBasisController = AtomCenteredBasisControllerFactory::produce(
      supSysGeom, envSystemZero->getSettings().basis.basisLibPath, envSystemZero->getSettings().basis.makeSphericalBasis,
      false, envSystemZero->getSettings().basis.firstECP, _potBasisLabel);

  auto supPotBasFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
      system->getSettings().grid.blocksize, 0.0, 2, supPotBasisController, this->_grid);

  OptEffPotential<SCFMode> oepCalcAct(actBasFuncOnGridController, actPotBasFuncOnGridController, oneEIntControllerAct,
                                      actDensOnGridCalc, actNOccOrbs, _smoothFactor, _singValThreshold, exc);
  OptEffPotential<SCFMode> oepCalc(supBasFuncOnGridController, supPotBasFuncOnGridController, oneEIntController,
                                   supDensOnGridCalc, supNOccOrbs, _smoothFactor, _singValThreshold, exc);

  /*
   * Run Wu-Yang and van Leeuwen-Baerends
   */

  /*
   * Active system initial guess: Fermi-Amaldi potential + Nuclear potential
   */
  GridPotential<SCFMode> initialGuessAct(this->_grid);
  GridPotential<SCFMode> negCoulOnGrid(this->_grid);
  GridPotential<RESTRICTED> tmp(this->_grid);
  CoulombPotentialOnGridCalculator::calculateElectronElectron<SCFMode>(tmp, this->_actDMatController->getDensityMatrix());
  for_spin(initialGuessAct, actNEl) {
    initialGuessAct_spin = tmp;
    initialGuessAct_spin *= (1.0 - (1.0 / actNEl_spin));
  };
  tmp.setZero();
  CoulombPotentialOnGridCalculator::calculateElectronNuclei(tmp, actSysGeom->getAtoms());
  for_spin(initialGuessAct) {
    initialGuessAct_spin += tmp;
  };

  /*
   * Supersystem initial guess: Fermi-Amaldi potential + Nuclear potential
   */
  GridPotential<SCFMode> initialGuess(this->_grid);
  tmp.setZero();
  CoulombPotentialOnGridCalculator::calculateElectronElectron<SCFMode>(tmp, supSysDensMat);
  for_spin(initialGuess, supNEl, negCoulOnGrid) {
    initialGuess_spin += tmp;
    negCoulOnGrid_spin -= tmp;
    initialGuess_spin *= (1.0 - (1.0 / supNEl_spin));
  };
  tmp.setZero();
  CoulombPotentialOnGridCalculator::calculateElectronNuclei(tmp, supSysGeom->getAtoms());
  for_spin(initialGuess, negCoulOnGrid) {
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
  oepCalcAct.calculateOEP(this->_densOnGridControllers[0]->getDensityOnGrid(), dummyOrbsAct, initialGuessAct);
  initialGuessAct += oepCalcAct.getPotentialOnGrid();

  /*
   * Then van Leeuwen--Baerends
   */
  if (_lbCycles != 0) {
    std::cout << "Running van Leeuwen-Baerends potential reconstruction" << std::endl;
    oepCalcAct.calculateOEPLB(this->_densOnGridControllers[0]->getDensityOnGrid(), dummyOrbsAct, initialGuessAct,
                              _lbDamping, _lbCycles);
    initialGuessAct += oepCalcAct.getPotentialOnGrid();
  }

  /*
   * Perform potential reconstruction to obtain missing parts of the embedding potential
   */
  std::cout << "Supersystem" << std::endl;

  /*
   * First Wu-Yang
   */
  std::cout << "Running Wu-Yang potential reconstruction" << std::endl;
  oepCalc.calculateOEP(this->_supersysDensOnGridController->getDensityOnGrid(), dummyOrbs, initialGuess);
  initialGuess += oepCalc.getPotentialOnGrid();

  /*
   * Then van Leeuwen--Baerends
   */
  if (_lbCycles != 0) {
    std::cout << "Running van Leeuwen-Baerends potential reconstruction" << std::endl;
    oepCalc.calculateOEPLB(this->_supersysDensOnGridController->getDensityOnGrid(), dummyOrbs, initialGuess, _lbDamping,
                           _lbCycles);
    initialGuess += oepCalc.getPotentialOnGrid();
  }

  /*
   * Construct non-additive kin pot as difference between the
   * supersystem potential and every other part of the total
   * potential
   */
  ScalarOperatorToMatrixAdder<SCFMode> scalarOpToMat(actBasFuncOnGridController, 0.0);

  GridPotential<SCFMode> potOnGrid(this->_grid);

  potOnGrid -= initialGuess;
  potOnGrid += initialGuessAct;
  /*
   * Get the matrix representation
   */
  scalarOpToMat.addScalarOperatorToMatrix(pot, potOnGrid);

  /*
   * Exact exchange part
   */
  if (exc > 0.0) {
    auto eCont = system->template getElectronicStructure<SCFMode>()->getEnergyComponentController();

    double exactXEnergy = 0.0;

    DensityMatrix<SCFMode> oepDMat(supSysBasisController);
    int factor = 2;
    CoefficientMatrix<SCFMode> supCoeffs = dummyOrbs->getCoefficients();
    if (SCFMode == Options::SCF_MODES::UNRESTRICTED)
      factor = 1;
    for (unsigned int i = 0; i < supSysBasisController->getNBasisFunctions(); i++) {
      for (unsigned int j = 0; j < supSysBasisController->getNBasisFunctions(); j++) {
        for_spin(supCoeffs, supNOccOrbs, oepDMat) {
          for (unsigned int k = 0; k < supNOccOrbs_spin; k++) {
            oepDMat_spin(i, j) += supCoeffs_spin(i, k) * supCoeffs_spin(j, k) * factor;
          }
        };
      }
    }

    auto oepDMatController = std::make_shared<DensityMatrixController<SCFMode>>(oepDMat);

    double thresh = envSystemZero->getSettings().basis.integralThreshold;
    ERIPotential<SCFMode> KSup(system, oepDMatController, functional.getHfExchangeRatio(), thresh,
                               system->getSettings().basis.integralIncrementThresholdStart,
                               system->getSettings().basis.integralIncrementThresholdEnd,
                               system->getSettings().basis.incrementalSteps);
    exactXEnergy += KSup.getEnergy(oepDMatController->getDensityMatrix());

    ERIPotential<SCFMode> K(system, system->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
                            functional.getHfExchangeRatio(), thresh, system->getSettings().basis.integralIncrementThresholdStart,
                            system->getSettings().basis.integralIncrementThresholdEnd,
                            system->getSettings().basis.incrementalSteps);
    exactXEnergy -= K.getEnergy(this->_actDMatController->getDensityMatrix());

    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_EXACT_EXCHANGE, exactXEnergy);
  }

  /*
   * Calculate energy from supersystem and subsystem orbitals
   */

  auto supKin = libint.compute1eInts(LIBINT_OPERATOR::kinetic, supSysBasisController);

  _supEnergy = 0.0;
  int factor = 2;
  CoefficientMatrix<SCFMode> supCoeffs = dummyOrbs->getCoefficients();
  if (SCFMode == Options::SCF_MODES::UNRESTRICTED)
    factor = 1;
  for (unsigned int i = 0; i < supSysBasisController->getNBasisFunctions(); i++) {
    for (unsigned int j = 0; j < supSysBasisController->getNBasisFunctions(); j++) {
      for_spin(supCoeffs, supNOccOrbs) {
        for (unsigned int k = 0; k < supNOccOrbs_spin; k++) {
          _supEnergy += supKin(i, j) * supCoeffs_spin(i, k) * supCoeffs_spin(j, k) * factor;
        }
      };
    }
  }
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd BUReconstructionPotential<SCFMode>::getGeomGradients() {
  Eigen::MatrixXd gradientContr(1, 3);
  gradientContr.setZero();
  throw SerenityError("No gradients for BUreconstructionPotentials available!");
  return gradientContr;
}

template class BUReconstructionPotential<Options::SCF_MODES::RESTRICTED>;
template class BUReconstructionPotential<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
