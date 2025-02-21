/**
 * @file ChargeElectronInteractionPotential.cpp
 *
 * @date Mai 20 2021
 * @author LarsHellmann
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
#include "potentials/ChargeElectronInteractionPotential.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/BeckePopulationCalculator.h"
#include "analysis/populationAnalysis/CHELPGPopulationCalculator.h"
#include "analysis/populationAnalysis/CM5PopulationCalculator.h"
#include "analysis/populationAnalysis/HirshfeldPopulationCalculator.h"
#include "analysis/populationAnalysis/IAOPopulationCalculator.h"
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/Basis.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "geometry/Atom.h"
#include "integrals/wrappers/Libint.h"
#include "misc/Timing.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ChargeElectronInteractionPotential<SCFMode>::ChargeElectronInteractionPotential(
    std::shared_ptr<SystemController> activeSystem, std::vector<std::shared_ptr<SystemController>> environmentSystems,
    std::shared_ptr<BasisController> basis, Options::POPULATION_ANALYSIS_ALGORITHMS chargeModel)
  : Potential<SCFMode>(basis), _actSystem(activeSystem), _chargeModel(chargeModel) {
  for (auto e : environmentSystems) {
    _envSystems.push_back(e);
    for (auto& atom : e->getGeometry()->getAtoms()) {
      atom->addSensitiveObject(ObjectSensitiveClass<Atom>::_self);
    }
  }
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& ChargeElectronInteractionPotential<SCFMode>::getMatrix() {
  Timings::takeTime("FDE -     Charge-Elec-Int Pot.");
  if (!_potential) {
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot) {
      pot_spin.setZero();
    };

    auto& libint = Libint::getInstance();

    // generate charges vector
    std::vector<std::pair<double, std::array<double, 3>>> pointCharges;
    for (auto e : _envSystems) {
      auto envCont = e.lock();

      auto geom = envCont->getGeometry();
      Eigen::VectorXd totPop(geom->getNAtoms());
      totPop.setZero();

      if (_chargeModel == Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN) {
        MullikenPopulationCalculator<SCFMode> calculator;
        auto pop = calculator.calculateMullikenPopulations(envCont);
        for_spin(pop) {
          totPop += pop_spin;
        };
      }
      if (_chargeModel == Options::POPULATION_ANALYSIS_ALGORITHMS::HIRSHFELD) {
        auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
            128, 0.0, 2, envCont->getBasisController(), envCont->getGridController());
        auto densOnGridCalc = std::make_shared<DensityOnGridCalculator<SCFMode>>(basFuncOnGridController, 0.0);
        auto densMatController = envCont->getElectronicStructure<SCFMode>()->getDensityMatrixController();
        auto densOnGridController =
            std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densOnGridCalc, densMatController);

        HirshfeldPopulationCalculator<SCFMode> hirshFeld(envCont, densOnGridController);
        auto pop = hirshFeld.getAtomPopulations();
        for_spin(pop) {
          totPop += pop_spin;
        };
      }
      if (_chargeModel == Options::POPULATION_ANALYSIS_ALGORITHMS::CM5) {
        auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
            128, 0.0, 2, envCont->getBasisController(), envCont->getGridController());
        auto densOnGridCalc = std::make_shared<DensityOnGridCalculator<SCFMode>>(basFuncOnGridController, 0.0);
        auto densMatController = envCont->getElectronicStructure<SCFMode>()->getDensityMatrixController();
        auto densOnGridController =
            std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densOnGridCalc, densMatController);
        auto hirshFeld = std::make_shared<HirshfeldPopulationCalculator<SCFMode>>(envCont, densOnGridController);
        CM5PopulationCalculator<SCFMode> cm5(envCont, hirshFeld);
        auto pop = cm5.getAtomPopulations();
        for_spin(pop) {
          totPop += pop_spin;
        };
      }
      if (_chargeModel == Options::POPULATION_ANALYSIS_ALGORITHMS::BECKE) {
        auto densPtr =
            std::make_shared<DensityMatrix<SCFMode>>(envCont->getElectronicStructure<SCFMode>()->getDensityMatrix());
        auto beckeCalculator = BeckePopulationCalculator<SCFMode>(envCont, densPtr);
        totPop = beckeCalculator.getAtomPopulations();
      }
      if (_chargeModel == Options::POPULATION_ANALYSIS_ALGORITHMS::CHELPG) {
        CHELPGPopulationCalculator<SCFMode> chelpg(envCont);
        auto atoms = envCont->getAtoms();
        auto pop = chelpg.getAtomPopulations();
        for (unsigned int iAtom = 0; iAtom < atoms.size(); ++iAtom) {
          totPop[iAtom] = atoms[iAtom]->getEffectiveCharge() - pop[iAtom];
        }
      }
      if (_chargeModel == Options::POPULATION_ANALYSIS_ALGORITHMS::IAO ||
          _chargeModel == Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell) {
        envCont->getElectronicStructure<SCFMode>()->getEnergy();
        auto pop = IAOPopulationCalculator<SCFMode>::calculateIAOPopulations(envCont);
        for_spin(pop) {
          totPop += pop_spin;
        };
      }

      for (unsigned int i = 0; i < geom->getNAtoms(); i++) {
        Eigen::Vector3d coords = (*geom)[i]->coords();
        std::array<double, 3> point = {coords(0), coords(1), coords(2)};
        pointCharges.push_back(std::make_pair(totPop(i), point));
      }
    }
    auto chargeInts = libint.compute1eInts(LIBINT_OPERATOR::nuclear, this->_basis, pointCharges);

    for_spin(pot) {
      pot_spin -= chargeInts;
    };
  }
  Timings::timeTaken("FDE -     Charge-Elec-Int Pot.");
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
double ChargeElectronInteractionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (!_potential)
    this->getMatrix();
  Timings::takeTime("FDE -     Charge-Elec-Int Pot.");
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("FDE -     Charge-Elec-Int Pot.");
  return energy;
};

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd ChargeElectronInteractionPotential<SCFMode>::getGeomGradients() {
  throw SerenityError(
      "Gradients for Charge Electron Interaction not implemented yet. Use full CoulombInteractionPotential");
}

template class ChargeElectronInteractionPotential<Options::SCF_MODES::RESTRICTED>;
template class ChargeElectronInteractionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
