/**
 * @file   EDAPotentials.cpp
 *
 * @date   Aug 18, 2017
 * @author Moritz Bensberg, Jan Unsleber
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
#include "potentials/bundles/EDAPotentials.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "potentials/ERIPotential.h"
#include "potentials/HCorePotential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
EDAPotentials<SCFMode>::EDAPotentials(std::shared_ptr<SystemController> systemA, std::shared_ptr<SystemController> systemB,
                                      std::shared_ptr<SystemController> super, EDAEnergyContributions contribution)
  : _a(systemA),
    _b(systemB),
    _super(super),
    _pot(super->getBasisController()),
    _AO2MO(super->getBasisController()),
    _MO2AO(super->getBasisController()),
    _hcorepot(new HCorePotential<SCFMode>(_super)),
    _hfpot(new ERIPotential<SCFMode>(
        _super, _super->getElectronicStructure<SCFMode>()->getDensityMatrixController(), 1.0,
        _super->getSettings().basis.integralThreshold, _super->getSettings().basis.integralIncrementThresholdStart,
        _super->getSettings().basis.integralIncrementThresholdEnd, _super->getSettings().basis.incrementalSteps)),
    _energyContr(contribution) {
  assert(_a);
  assert(_b);
  assert(_super);
  _a->getElectronicStructure<SCFMode>();
  _b->getElectronicStructure<SCFMode>();

  // Generate transformation matrices
  // 1. Coefficient matrix for AO-MO and MO-AO transformation
  auto coeffA(_a->getActiveOrbitalController<SCFMode>()->getCoefficients());
  auto coeffB(_b->getActiveOrbitalController<SCFMode>()->getCoefficients());
  for_spin(_AO2MO, _MO2AO, coeffA, coeffB) {
    _AO2MO_spin.setZero();
    _AO2MO_spin.topLeftCorner(coeffA_spin.rows(), coeffA_spin.cols()) = coeffA_spin;
    _AO2MO_spin.bottomRightCorner(coeffB_spin.rows(), coeffB_spin.cols()) = coeffB_spin;
    _MO2AO_spin = _AO2MO_spin.inverse();
  };

  // 2. Orthogonalization Matrices using modified overlap
  // transform to MO basis
  MatrixInBasis<SCFMode> S(_super->getBasisController());
  for_spin(S, _AO2MO) {
    S_spin = _AO2MO_spin.adjoint() * _super->getOneElectronIntegralController()->getOverlapIntegrals() * _AO2MO_spin;
  };

  // blocking
  this->setBlocksZero(S, _energyContr);

  // transform back to AO basis
  for_spin(S, _MO2AO) {
    S_spin = _MO2AO_spin.adjoint() * S_spin * _MO2AO_spin;
  };
  _super->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->useCustomOverlap(S);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> EDAPotentials<SCFMode>::getFockMatrix(const DensityMatrix<SCFMode>& P,
                                                          std::shared_ptr<EnergyComponentController> energies) {
  // Set up initial potential
  for_spin(_pot) {
    _pot_spin.setZero();
  };
  _pot += _hcorepot->getMatrix();
  _pot += _hfpot->getMatrix();

  if (_energyContr == EDAEnergyContributions::ES or _energyContr == EDAEnergyContributions::ESPL) {
    // subtract supersystem exchange
    _pot -= _hfpot->getXPotential();
    // calculate subsystem exchange for a and b
    DensityMatrix<SCFMode> dMatA(_a->getBasisController());
    DensityMatrix<SCFMode> dMatB(_b->getBasisController());
    for_spin(dMatA, dMatB, P) {
      dMatA_spin = P_spin.topLeftCorner(dMatA_spin.rows(), dMatA_spin.cols());
      dMatB_spin = P_spin.bottomRightCorner(dMatB_spin.rows(), dMatB_spin.cols());
    };
    auto dMatCA = std::make_shared<DensityMatrixController<SCFMode>>(dMatA);
    auto dMatCB = std::make_shared<DensityMatrixController<SCFMode>>(dMatB);
    auto hfpota = std::make_shared<ERIPotential<SCFMode>>(_a, dMatCA, 1.0, _super->getSettings().basis.integralThreshold,
                                                          _super->getSettings().basis.integralIncrementThresholdStart,
                                                          _super->getSettings().basis.integralIncrementThresholdEnd,
                                                          _super->getSettings().basis.incrementalSteps);
    auto hfpotb = std::make_shared<ERIPotential<SCFMode>>(_b, dMatCB, 1.0, _super->getSettings().basis.integralThreshold,
                                                          _super->getSettings().basis.integralIncrementThresholdStart,
                                                          _super->getSettings().basis.integralIncrementThresholdEnd,
                                                          _super->getSettings().basis.incrementalSteps);
    auto xa(hfpota->getXPotential());
    auto xb(hfpotb->getXPotential());
    // add subsystem exchange potential of a and b
    for_spin(_pot, xa, xb) {
      _pot_spin.topLeftCorner(xa_spin.rows(), xa_spin.cols()) += xa_spin;
      _pot_spin.bottomRightCorner(xb_spin.rows(), xb_spin.cols()) += xb_spin;
    };
  }

  // transform to MO basis
  FockMatrix<SCFMode> h(_hcorepot->getMatrix());
  for_spin(h, _pot, _AO2MO) {
    h_spin = _AO2MO_spin.adjoint() * h_spin * _AO2MO_spin;
    _pot_spin = _AO2MO_spin.adjoint() * _pot_spin * _AO2MO_spin;
  };

  // blocking
  this->setBlocksZero(h, _energyContr);
  this->setBlocksZero(_pot, _energyContr);

  // transform back to AO basis
  for_spin(h, _pot, _MO2AO) {
    h_spin = _MO2AO_spin.adjoint() * h_spin * _MO2AO_spin;
    _pot_spin = _MO2AO_spin.adjoint() * _pot_spin * _MO2AO_spin;
  };

  // EDA energy
  double energy = _super->getGeometry()->getCoreCoreRepulsion();
  for_spin(_pot, h, P) {
    energy += 0.5 * (_pot_spin + h_spin).cwiseProduct(P_spin).sum();
  };
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::HF_ENERGY, energy);
  _EDAEnergy = energy;

  return _pot;
}

template<Options::SCF_MODES SCFMode>
void EDAPotentials<SCFMode>::setBlocksZero(MatrixInBasis<SCFMode>& matrix, EDAEnergyContributions x) {
  // Get the block selection matrix according to the choice in x.
  Eigen::Matrix<bool, 4, 4> blockSelectionMatrix;
  switch (x) {
    case EDAEnergyContributions::ESXCT:
      blockSelectionMatrix << 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1;
      break;
    case EDAEnergyContributions::ESXPLX:
    case EDAEnergyContributions::ESPL:
      blockSelectionMatrix << 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1;
      break;
    case EDAEnergyContributions::ESXEX:
      blockSelectionMatrix << 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1;
      break;
    case EDAEnergyContributions::ESX:
    case EDAEnergyContributions::ES:
      blockSelectionMatrix << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
      break;
  }

  auto nOCCA = _a->getNOccupiedOrbitals<SCFMode>();
  auto nOCCB = _b->getNOccupiedOrbitals<SCFMode>();

  // Set the blocks to 0.
  for_spin(matrix, nOCCA, nOCCB) {
    // get the number of the basis functions for the systems.
    const unsigned int nA(_a->getBasisController()->getNBasisFunctions());
    const unsigned int nB(_b->getBasisController()->getNBasisFunctions());
    const unsigned int nOcc_A = nOCCA_spin;
    const unsigned int nOcc_B = nOCCB_spin;
    const unsigned int nVirt_A(nA - nOcc_A);
    const unsigned int nVirt_B(nB - nOcc_B);

    // diagonal elements are not checked, they are never zero
    if (!blockSelectionMatrix(0, 1))
      matrix_spin.block(nOcc_A, 0, nVirt_A, nOcc_A).setZero();
    if (!blockSelectionMatrix(0, 2))
      matrix_spin.block(nA, 0, nOcc_B, nOcc_A).setZero();
    if (!blockSelectionMatrix(0, 3))
      matrix_spin.block(nA + nOcc_B, 0, nVirt_B, nOcc_A).setZero();
    if (!blockSelectionMatrix(1, 0))
      matrix_spin.block(0, nOcc_A, nOcc_A, nVirt_A).setZero();
    if (!blockSelectionMatrix(1, 2))
      matrix_spin.block(nA, nOcc_A, nOcc_B, nVirt_A).setZero();
    if (!blockSelectionMatrix(1, 3))
      matrix_spin.block(nA + nOcc_B, nOcc_A, nVirt_B, nVirt_A).setZero();
    if (!blockSelectionMatrix(2, 0))
      matrix_spin.block(0, nA, nOcc_A, nOcc_B).setZero();
    if (!blockSelectionMatrix(2, 1))
      matrix_spin.block(nOcc_A, nA, nVirt_A, nOcc_B).setZero();
    if (!blockSelectionMatrix(2, 3))
      matrix_spin.block(nA + nOcc_B, nA, nVirt_B, nOcc_B).setZero();
    if (!blockSelectionMatrix(3, 0))
      matrix_spin.block(0, nA + nOcc_B, nOcc_A, nVirt_B).setZero();
    if (!blockSelectionMatrix(3, 1))
      matrix_spin.block(nOcc_A, nA + nOcc_B, nVirt_A, nVirt_B).setZero();
    if (!blockSelectionMatrix(3, 2))
      matrix_spin.block(nA, nA + nOcc_B, nOcc_B, nVirt_B).setZero();
  };
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd EDAPotentials<SCFMode>::getGradients() {
  throw SerenityError("No geometry gradients available for the EDA potentials.");
  Eigen::MatrixXd gradientContr(1, 3);
  return gradientContr;
}

template class EDAPotentials<Options::SCF_MODES::RESTRICTED>;
template class EDAPotentials<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
