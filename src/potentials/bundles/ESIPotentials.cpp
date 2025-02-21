/**
 * @file ESIPotentials.cpp
 *
 * @date Nov 24, 2016
 * @author: Jan Unsleber
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
#include "potentials/bundles/ESIPotentials.h"
/* Include Serenity Internal Headers */
#include "energies/EnergyContributions.h"
#include "integrals/wrappers/Libint.h"
#include "potentials/ChargeElectronInteractionPotential.h"
#include "potentials/CoulombInteractionPotential.h"
#include "potentials/NEInteractionPotential.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ESIPotentials<SCFMode>::ESIPotentials(std::shared_ptr<SystemController> actSystem,
                                      std::vector<std::shared_ptr<SystemController>> envSystems,
                                      std::shared_ptr<DensityMatrixController<SCFMode>> activeDMat,
                                      std::shared_ptr<const Geometry> activeGeom,
                                      std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDMats,
                                      std::vector<std::shared_ptr<const Geometry>> envGeoms, unsigned int firstPassiveSystemIndex,
                                      bool useCharges, Options::POPULATION_ANALYSIS_ALGORITHMS chargeModel)
  : _actSystem(actSystem),
    _envSystems(envSystems),
    _activeDMat(activeDMat),
    _activeGeom(activeGeom),
    _envDMats(envDMats),
    _envGeoms(envGeoms),
    _nePot(nullptr),
    _cePot(nullptr),
    _coulPot(nullptr),
    _enAttr(nullptr),
    _nnRep(nullptr) {
  Timings::takeTime("FDE -            ESI Pot.");
  _nePot.reset(new NEInteractionPotential<SCFMode>(_actSystem, _envSystems,
                                                   _activeDMat->getDensityMatrix().getBasisController(), _envGeoms));
  _cePot.reset(new ChargeElectronInteractionPotential<SCFMode>(
      _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), chargeModel));
  if (firstPassiveSystemIndex < envSystems.size()) {
    std::vector<std::shared_ptr<SystemController>> passiveSystems;
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> passiveDensityMatrices;
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> environmentDensityMatrices;
    _envSystems = {};
    for (unsigned int envIndex = 0; envIndex < envSystems.size(); ++envIndex) {
      if (envIndex < firstPassiveSystemIndex) {
        _envSystems.push_back(envSystems[envIndex]);
        environmentDensityMatrices.push_back(_envDMats[envIndex]);
      }
      else {
        passiveSystems.push_back(envSystems[envIndex]);
        passiveDensityMatrices.push_back(_envDMats[envIndex]);
      }
    }
    if (useCharges) {
      _coulPot.reset(new ChargeElectronInteractionPotential<SCFMode>(
          _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), chargeModel));
      _passiveCoulPot.reset(new ChargeElectronInteractionPotential<SCFMode>(
          _actSystem, passiveSystems, _activeDMat->getDensityMatrix().getBasisController(), chargeModel));
    }
    else {
      _coulPot.reset(new CoulombInteractionPotential<SCFMode>(
          _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), environmentDensityMatrices));
      _passiveCoulPot.reset(new CoulombInteractionPotential<SCFMode>(
          _actSystem, passiveSystems, _activeDMat->getDensityMatrix().getBasisController(), passiveDensityMatrices, true));
    }
  }
  else {
    if (useCharges) {
      _coulPot.reset(new ChargeElectronInteractionPotential<SCFMode>(
          _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), chargeModel));
    }
    else {
      _coulPot.reset(new CoulombInteractionPotential<SCFMode>(
          _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), _envDMats));
    }
  }
  for (auto& atom : _activeGeom->getAtoms()) {
    atom->addSensitiveObject(ObjectSensitiveClass<Atom>::_self);
  }
  for (auto& geo : _envGeoms) {
    for (auto& atom : geo->getAtoms()) {
      atom->addSensitiveObject(ObjectSensitiveClass<Atom>::_self);
    }
  }
  for (auto& mat : _envDMats) {
    mat->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  }
  Timings::timeTaken("FDE -            ESI Pot.");
}

template<Options::SCF_MODES SCFMode>
ESIPotentials<SCFMode>::ESIPotentials(
    std::shared_ptr<SystemController> actSystem, std::vector<std::shared_ptr<SystemController>> envSystems,
    std::shared_ptr<DensityMatrixController<SCFMode>> activeDMat, std::shared_ptr<const Geometry> activeGeom,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDMats,
    std::vector<std::shared_ptr<const Geometry>> envGeoms, const std::shared_ptr<BasisController> actAuxBasis,
    std::vector<std::shared_ptr<BasisController>> envAuxBasis, unsigned int firstPassiveSystemIndex, bool useCharges,
    Options::POPULATION_ANALYSIS_ALGORITHMS chargeModel)
  : _actSystem(actSystem),
    _envSystems(envSystems),
    _activeDMat(activeDMat),
    _activeGeom(activeGeom),
    _envDMats(envDMats),
    _envGeoms(envGeoms),
    _nePot(nullptr),
    _cePot(nullptr),
    _coulPot(nullptr),
    _enAttr(nullptr),
    _nnRep(nullptr) {
  Timings::takeTime("FDE -            ESI Pot.");
  _nePot.reset(new NEInteractionPotential<SCFMode>(_actSystem, _envSystems,
                                                   _activeDMat->getDensityMatrix().getBasisController(), _envGeoms));
  _cePot.reset(new ChargeElectronInteractionPotential<SCFMode>(
      _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), chargeModel));

  if (firstPassiveSystemIndex < envSystems.size()) {
    std::vector<std::shared_ptr<SystemController>> passiveSystems;
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> passiveDensityMatrices;
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> environmentDensityMatrices;
    std::vector<std::shared_ptr<BasisController>> envAuxBasisController;
    std::vector<std::shared_ptr<BasisController>> passiveAuxBasisController;
    _envSystems = {};
    for (unsigned int envIndex = 0; envIndex < envSystems.size(); ++envIndex) {
      if (envIndex < firstPassiveSystemIndex) {
        _envSystems.push_back(envSystems[envIndex]);
        environmentDensityMatrices.push_back(_envDMats[envIndex]);
        envAuxBasisController.push_back(envAuxBasis[envIndex]);
      }
      else {
        passiveSystems.push_back(envSystems[envIndex]);
        passiveDensityMatrices.push_back(_envDMats[envIndex]);
        passiveAuxBasisController.push_back(envAuxBasis[envIndex]);
      }
    }

    if (useCharges) {
      _coulPot.reset(new ChargeElectronInteractionPotential<SCFMode>(
          _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), chargeModel));
      _passiveCoulPot.reset(new ChargeElectronInteractionPotential<SCFMode>(
          _actSystem, passiveSystems, _activeDMat->getDensityMatrix().getBasisController(), chargeModel));
    }
    else {
      _coulPot.reset(new CoulombInteractionPotential<SCFMode>(
          _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), environmentDensityMatrices,
          actAuxBasis, envAuxBasisController));
      _passiveCoulPot.reset(new CoulombInteractionPotential<SCFMode>(
          _actSystem, passiveSystems, _activeDMat->getDensityMatrix().getBasisController(), passiveDensityMatrices,
          actAuxBasis, passiveAuxBasisController, true));
    }
  }
  else {
    if (useCharges) {
      _coulPot.reset(new ChargeElectronInteractionPotential<SCFMode>(
          _actSystem, _envSystems, _activeDMat->getDensityMatrix().getBasisController(), chargeModel));
    }
    else {
      _coulPot.reset(new CoulombInteractionPotential<SCFMode>(_actSystem, _envSystems,
                                                              _activeDMat->getDensityMatrix().getBasisController(),
                                                              _envDMats, actAuxBasis, envAuxBasis));
    }
  }
  for (auto& atom : _activeGeom->getAtoms()) {
    atom->addSensitiveObject(ObjectSensitiveClass<Atom>::_self);
  }
  for (auto& geo : _envGeoms) {
    for (auto& atom : geo->getAtoms()) {
      atom->addSensitiveObject(ObjectSensitiveClass<Atom>::_self);
    }
  }
  for (auto& mat : _envDMats) {
    mat->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  }
  Timings::timeTaken("FDE -            ESI Pot.");
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> ESIPotentials<SCFMode>::getFockMatrix(const DensityMatrix<SCFMode>& P,
                                                          std::shared_ptr<EnergyComponentController> energies) {
  Timings::takeTime("FDE -            ESI Pot.");
  // electron(act) - electron(env) Coulomb
  //               and
  // electron(act) - nucleii(env) coulomb
  const auto& Ven = _nePot->getMatrix();
  FockMatrix<SCFMode> Vee = _coulPot->getMatrix();
  if (_passiveCoulPot)
    Vee = Vee + _passiveCoulPot->getMatrix();
  double en = _nePot->getEnergy(P);
  double ee = _coulPot->getEnergy(P);
  if (_passiveCoulPot)
    ee += _passiveCoulPot->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_ELECTRONS_ENV_NUCLEI_COULOMB, en);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_ELECTRONS_ENV_ELECTRONS_COULOMB, ee);

  // nucleii(act) - electron(env) Coulomb
  //               and
  // nucleii(act) - nucleii(env) coulomb
  if (!_enAttr || !_nnRep) {
    _enAttr.reset(new double(0.0e0));
    auto& libint = Libint::getInstance();
    for (unsigned int env = 0; env < _envDMats.size(); ++env) {
      auto envbasisController = _envDMats[env]->getDensityMatrix().getBasisController();
      auto envInts = libint.compute1eInts(LIBINT_OPERATOR::nuclear, envbasisController, _activeGeom->getAtoms());
      auto dmat = _envDMats[env]->getDensityMatrix();
      for_spin(dmat) {
        *_enAttr += envInts.cwiseProduct(dmat_spin).sum();
      };
    }

    _nnRep.reset(new double(0.0e0));
    for (const auto& envGeo : _envGeoms) {
      for (const auto& envAtom : envGeo->getAtoms()) {
        if (envAtom->isDummy())
          continue;
        for (const auto& actAtom : _activeGeom->getAtoms()) {
          if (actAtom->isDummy())
            continue;
          *_nnRep += (envAtom->getEffectiveCharge() * actAtom->getEffectiveCharge()) / distance(*envAtom, *actAtom);
        }
      }
    }
  }
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NUCLEI_ENV_ELECTRONS_COULOMB, *_enAttr);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NUCLEI_ENV_NUCLEI_COULOMB, *_nnRep);
  Timings::timeTaken("FDE -            ESI Pot.");
  return Ven + Vee;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd ESIPotentials<SCFMode>::getGradients() {
  auto gradients = _nePot->getGeomGradients();
  gradients += _coulPot->getGeomGradients();
  if (_passiveCoulPot)
    _passiveCoulPot->getGeomGradients();

  return gradients;
}

template class ESIPotentials<Options::SCF_MODES::RESTRICTED>;
template class ESIPotentials<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
