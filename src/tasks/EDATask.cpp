/**
 * @file   EDATask.cpp
 *
 * @date   Mar 29, 2016
 * @author Moritz Bensberg
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
#include "tasks/EDATask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutput.h"
#include "parameters/Constants.h"
#include "potentials/bundles/EDAPotentials.h"
#include "scf/Scf.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <cmath>
#include <iostream>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
EDATask<SCFMode>::EDATask(std::vector<std::shared_ptr<SystemController>> systems) {
  if (systems.size() != 2)
    throw SerenityError("ERROR: Two active systems need to be given for an EDA calculation.");
  _systemA = systems[0];
  _systemB = systems[1];
};

template<Options::SCF_MODES SCFMode>
void EDATask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, {_systemA, _systemB});
  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);

  _systemA->getElectronicStructure<SCFMode>();
  _systemB->getElectronicStructure<SCFMode>();

  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;

  double eES, eESX, eESXCT, eESXPLX, eESPL, eESXEX, eInt, e0;
  {
    auto super = (*_systemA) + (*_systemB);
    auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ES);
    const auto& P(super->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
    pot->getFockMatrix(P, super->template getElectronicStructure<SCFMode>()->getEnergyComponentController());
    eES = pot->getEDAEnergyContribution();
  }
  {
    auto super = (*_systemA) + (*_systemB);
    auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESX);
    const auto& P(super->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
    pot->getFockMatrix(P, super->template getElectronicStructure<SCFMode>()->getEnergyComponentController());
    eESX = pot->getEDAEnergyContribution();
  }
  {
    auto super = (*_systemA) + (*_systemB);
    auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXCT);
    Scf<SCFMode>::perform(settings, super->template getElectronicStructure<SCFMode>(), pot);
    eESXCT = pot->getEDAEnergyContribution();
  }
  {
    auto super = (*_systemA) + (*_systemB);
    auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXPLX);
    Scf<SCFMode>::perform(settings, super->template getElectronicStructure<SCFMode>(), pot);
    eESXPLX = pot->getEDAEnergyContribution();
  }
  {
    auto super = (*_systemA) + (*_systemB);
    auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESPL);
    Scf<SCFMode>::perform(settings, super->template getElectronicStructure<SCFMode>(), pot);
    eESPL = pot->getEDAEnergyContribution();
  }
  {
    auto super = (*_systemA) + (*_systemB);
    auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXEX);
    const auto& P = super->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
    super->template getElectronicStructure<SCFMode>()->getMolecularOrbitals()->updateOrbitals(
        pot->getFockMatrix(P, super->template getElectronicStructure<SCFMode>()->getEnergyComponentController()),
        super->getOneElectronIntegralController());
    pot->getFockMatrix(P, super->template getElectronicStructure<SCFMode>()->getEnergyComponentController());
    eESXEX = pot->getEDAEnergyContribution();
  }
  {
    auto super = (*_systemA) + (*_systemB);
    ScfTask<SCFMode> scf(super);
    scf.run();
    e0 = _systemA->template getElectronicStructure<SCFMode>()->getEnergy() +
         _systemB->template getElectronicStructure<SCFMode>()->getEnergy();
    eInt = super->template getElectronicStructure<SCFMode>()->getEnergy() - e0;
  }
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);

  // Calculate and print final energies.
  double eMix(eInt + e0 - eESXEX - eESPL - eESXCT - eESXPLX + 2.0 * eESX + eESPL);

  printf("%11s\n", "-------------");
  printf("%12s\n", "Summary");
  printf("%11s\n", "-------------");
  printf("%12s %9f %8s %14f %8s %14f %6s\n", "E_{ES}= ", (eES - e0), "Hartree", (eES - e0) * HARTREE_TO_KCAL_PER_MOL,
         " kcal/mol", (eES - e0) * HARTREE_TO_KJ_PER_MOL, "kJ/mol");
  printf("%12s %9f %8s %14f %8s %14f %6s\n", "E_{EX}= ", (eESXEX - eES), "Hartree",
         (eESXEX - eES) * HARTREE_TO_KCAL_PER_MOL, " kcal/mol", (eESXEX - eES) * HARTREE_TO_KJ_PER_MOL, "kJ/mol");
  printf("%12s %9f %8s %14f %8s %14f %6s\n", "E_{PL}= ", (eESPL - eES), "Hartree",
         (eESPL - eES) * HARTREE_TO_KCAL_PER_MOL, " kcal/mol", (eESPL - eES) * HARTREE_TO_KJ_PER_MOL, "kJ/mol");
  printf("%12s %9f %8s %14f %8s %14f %6s\n", "E_{CT}= ", (eESXCT - eESX), "Hartree",
         (eESXCT - eESX) * HARTREE_TO_KCAL_PER_MOL, " kcal/mol", (eESXCT - eESX) * HARTREE_TO_KJ_PER_MOL, "kJ/mol");
  printf("%12s %9f %8s %14f %8s %14f %6s\n", "E_{EXPL}= ", (eES + eESXPLX - eESX - eESPL), "Hartree",
         (eES + eESXPLX - eESX - eESPL) * HARTREE_TO_KCAL_PER_MOL, " kcal/mol",
         (eES + eESXPLX - eESX - eESPL) * HARTREE_TO_KJ_PER_MOL, "kJ/mol");
  printf("%12s %9f %8s %14f %8s %14f %6s\n", "E_{MIX'}= ", eMix, "Hartree", eMix * HARTREE_TO_KCAL_PER_MOL, " kcal/mol",
         eMix * HARTREE_TO_KJ_PER_MOL, "kJ/mol");
  printf("%11s\n", "-------------");
  printf("%12s %9f %8s %14f %8s %14f %6s\n", "E_{inter}= ", eInt, "Hartree", eInt * HARTREE_TO_KCAL_PER_MOL,
         " kcal/mol", eInt * HARTREE_TO_KJ_PER_MOL, "kJ/mol");
  printf("%11s\n", "-------------");
}

template class EDATask<Options::SCF_MODES::RESTRICTED>;
template class EDATask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
