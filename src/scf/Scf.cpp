/**
 * @file   Scf.cpp
 *
 * @date   last rework Nov 29. 2016
 * @author Thomas Dresselhaus, Jan Unsleber
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
#include "scf/Scf.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "dft/Functional.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutput.h"
#include "misc/SerenityError.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
#include "potentials/bundles/PotentialBundle.h"
#include "scf/ConvergenceController.h"
#include "scf/ROHF.h"
#include "settings/Settings.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void Scf<SCFMode>::perform(const Settings& settings, std::shared_ptr<ElectronicStructure<SCFMode>> es,
                           std::shared_ptr<PotentialBundle<SCFMode>> potentials, bool allowNotConverged,
                           std::shared_ptr<SPMatrix<SCFMode>> momMatrix, unsigned int momCycles, bool useALMO) {
  allowNotConverged = (allowNotConverged || settings.scf.allowNotConverged);
  es->setDiskMode(false, "", "");
  if (useALMO)
    es->getDensityMatrixController()->setALMO(es->getOneElectronIntegralController());
  auto energyComponentController = es->getEnergyComponentController();
  auto orbitalController = es->getMolecularOrbitals();
  orbitalController->setCanOrthTh(settings.scf.canOrthThreshold);
  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  // Check if range seperate hybrid is used; then it is more efficient to keep those libint engines as well
  auto functional = settings.customFunc.basicFunctionals.size() ? Functional(settings.customFunc)
                                                                : resolveFunctional(settings.dft.functional);
  if (functional.isRSHybrid()) {
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 4);
  }
  /*
   * Perform one Roothan-Hall step without pre-diagonalization to
   * savely convert MOs and density matrix into the current basis
   * (important when the basis changes during the lifetime of
   * the ElectronicStructure)
   */
  es->attachPotentials(potentials);
  auto F(potentials->getFockMatrix(es->getDensityMatrix(), energyComponentController));
  ConvergenceController<SCFMode> convergenceController(settings, es->getDensityMatrixController(), orbitalController,
                                                       es->getOneElectronIntegralController(), energyComponentController);

  convergenceController.accelerateConvergence(F, es->getDensityMatrix());
  orbitalController->updateOrbitals(convergenceController.getLevelshift(), F, es->getOneElectronIntegralController());
  unsigned int counter = 0;
  bool converged = false;

  /*
   * SCF loop
   */
  while (true) {
    counter++;
    takeTime("the fock matrix formation");
    F = potentials->getFockMatrix(es->getDensityMatrix(), energyComponentController);

    if (settings.scf.rohf != Options::ROHF_TYPES::NONE) {
      ROHF<SCFMode>::addConstraint(F, es, settings.scf.rohf, settings.scf.suhfLambda);
    }

    timeTaken(2, "the fock matrix formation");

    takeTime("the generation of new MOs (possibly with conv. acc.)");
    /*
     * Hook in for convergence acceleration
     */
    convergenceController.accelerateConvergence(F, es->getDensityMatrix());
    orbitalController->updateOrbitals(convergenceController.getLevelshift(), F, es->getOneElectronIntegralController(),
                                      momMatrix);
    // Update MOM matrix (MOM and not IMOM procedure).
    if ((counter - 1 < momCycles) && momMatrix) {
      auto& mom = (*momMatrix);
      auto C = orbitalController->getCoefficients();
      auto nocc = es->getNOccupiedOrbitals();
      for_spin(mom, C, nocc) {
        mom_spin = C_spin.leftCols(nocc_spin);
      };
    }

    if (counter % settings.scf.writeRestart == 0) {
      orbitalController->toHDF5(settings.path + "tmp", settings.identifier);
    }
    timeTaken(2, "the generation of new MOs (possibly with conv. acc.)");
    takeTime("the convergence check");
    converged = convergenceController.checkConvergence();
    timeTaken(2, "the convergence check");
    if (converged and counter > 2) {
      if (iOOptions.printSCFCycleInfo) {
        printf("    Converged after %3i cycles. Loop exited.\n\n", counter);
      }
      es->state = ES_STATE::CONVERGED;
      /*
       * Replace temporary orbital files with the converged electronic structure
       */
      std::remove((settings.path + "tmp.orbs.res.h5").c_str());
      std::remove((settings.path + "tmp.orbs.unres.h5").c_str());
      es->setFockMatrix(F);
      es->toHDF5(settings.path + settings.name, settings.identifier);
      break;
    }
    if (counter >= settings.scf.maxCycles and allowNotConverged) {
      es->state = ES_STATE::CONVERGED;
      libint.clearAllEngines();
      WarningTracker::printWarning((std::string) "WARNING: SCF did NOT converge after " + settings.scf.maxCycles +
                                       " cycles!!! Continuing...",
                                   iOOptions.printSCFCycleInfo);
      break;
    }
    else if (counter >= settings.scf.maxCycles and !allowNotConverged) {
      es->state = ES_STATE::FAILED;
      libint.clearAllEngines();
      energyComponentController->printAllComponents();
      if (!allowNotConverged)
        throw SerenityError((std::string) "Cancelling SCF after " + settings.scf.maxCycles + " cycles. NOT CONVERGED!!!");
      break;
    }
  }

  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  if (functional.isRSHybrid()) {
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 2);
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 3);
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 4);
  }
}

template class Scf<Options::SCF_MODES::RESTRICTED>;
template class Scf<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
