/**
 * @file   ScfTask.h
 *
 * @date   Mar 7, 2014
 * @author Thomas Dresselhaus
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
#ifndef SCFTASK_H_
#define SCFTASK_H_
/* Include Serenity Internal Headers */
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //Local correlation settings for local MP2.
#include "settings/Options.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class ElectronicStructureCalculator;
class SystemController;
template<Options::SCF_MODES SCFMode>
class PotentialBundle;

using namespace Serenity::Reflection;

struct ScfTaskSettings {
  ScfTaskSettings()
    : restart(false),
      mp2Type(Options::MP2_TYPES::RI),
      maxResidual(1e-5),
      maxCycles(100),
      fractionalDegeneracy(false),
      skipSCF(false),
      allowNotConverged(false) {
    lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
    lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
  }
  REFLECTABLE((bool)restart, (Options::MP2_TYPES)mp2Type, (double)maxResidual, (int)maxCycles,
              (bool)fractionalDegeneracy, (bool)skipSCF, (bool)allowNotConverged)
 public:
  LocalCorrelationSettings lcSettings;
};

/**
 * @class ScfTask ScfTask.h
 * @brief Performs an SCF calculation, i.e. an electronic structure calculation.
 */
template<Options::SCF_MODES SCFMode>
class ScfTask : public Task {
 public:
  /**
   * @param systemController
   */
  ScfTask(std::shared_ptr<SystemController> systemController);

  virtual ~ScfTask() = default;

  void run();

  /**
   * @brief The settings/keywords for SCFTask:
   *        -restart :              Uses old orbitals to restart SCF (tries to restart from h5
   *                                file if available)
   *        -mp2Type:               Type of MP2 used for double-hybrid correlation part.
   *        -maxResidual:           Maximum residual threshold for local MP2.
   *        -maxCycles:             Maximum number of iterations before cancelling the amplitude optimization
   *                                in local MP2.
   *        -fractionalDegeneracy:  Allows fractional occupations for degenerate orbitals in order
   *                                to achieve spherical symmetry in atom-SCFs.
   *        -skipSCF:               Skip the SCF procedure and calculate only the energy.
   *        -allowNotConverged:     If the maximum number of SCF cycles is reached Serenity will continue
   *                                even with non-converged orbitals.
   *        -lcSettings:            Local correlation settings for local MP2.
   */
  ScfTaskSettings settings;

  /**
   * @brief Getter for the electronic structure.
   * @return Returns the electronic structure.
   */
  std::shared_ptr<ElectronicStructureCalculator<SCFMode>> getElectronicStructureCalculator();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(ScfTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
    }
    else if (!c.lcSettings.visitSettings(v, blockname)) {
      throw SerenityError((std::string) "Unknown settings block in ScfTaskSettings: " + blockname);
    }
  }

 private:
  const std::shared_ptr<SystemController> _systemController;
  /*
   * Helper Functions.
   */
  void performSCF(std::shared_ptr<PotentialBundle<SCFMode>> potentials);
  std::shared_ptr<PotentialBundle<SCFMode>> getPotentialBundle();
  void printResults();
  void printHeader();
  void loadRestartFiles();
  void finalDFTEnergyEvaluation();
  void calculateMP2Contribution();
  void calculateDispersionCorrection();
};

} /* namespace Serenity */

#endif /* SCFTASK_H_ */
