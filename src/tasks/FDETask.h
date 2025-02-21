/**
 * @file   FDETask.h
 * @author Jan Unsleber, Kevin Klahr, Thomas Dresselhaus
 *
 * @date   last reworked on Nov 11. 2016
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
#ifndef FDETASK_H
#define FDETASK_H

/* Include Serenity Internal Headers */
#include "data/grid/GridPotential.h"
#include "dft/Functional.h"
#include "misc/SerenityError.h"
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "settings/EmbeddingSettings.h"
#include "settings/Reflection.h"
#include "tasks/LocalizationTask.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory> //smart ptr.
#include <vector>

namespace Serenity {
class SystemController;
template<Options::SCF_MODES SCFMode>
class PotentialBundle;
class EnergyComponentController;

using namespace Serenity::Reflection;

struct FDETaskSettings {
  FDETaskSettings()
    : gridCutOff(-1.0),
      smallSupersystemGrid(false),
      finalGrid(true),
      calculateMP2Energy(true),
      calculateUnrelaxedMP2Density(false),
      maxResidual(1e-5),
      maxCycles(100),
      calculateEnvironmentEnergy(false),
      mp2Type(Options::MP2_TYPES::LOCAL),
      calculateSolvationEnergy(false),
      skipSCF(false) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
    lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
    lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
    loc.splitValenceAndCore = true;
  }
  REFLECTABLE((double)gridCutOff, (bool)smallSupersystemGrid, (bool)finalGrid, (bool)calculateMP2Energy,
              (bool)calculateUnrelaxedMP2Density, (double)maxResidual, (int)maxCycles, (bool)calculateEnvironmentEnergy,
              (Options::MP2_TYPES)mp2Type, (bool)calculateSolvationEnergy, (bool)skipSCF)
 public:
  LocalCorrelationSettings lcSettings;
  EmbeddingSettings embedding;
  LocalizationTaskSettings loc;
  ///@brief The index of the first passive environment subsystem that will never
  ///       be changed during a freeze-and-thaw procedure.
  unsigned int firstPassiveSystemIndex = 9999999;
  bool initializeSuperMolecularSurface = true;
};

/**
 * @class FDETask FDETask.h
 * @brief A task for an FDE run (1 active, x environment systems).
 *
 * References:
 * C.R. Jacob, J. Neugebauer, Subsystem Density-Functional Theory, WIREs Comput. Mol. Sci. 4 (2014), 325.
 */
template<Options::SCF_MODES SCFMode>
class FDETask : public Task {
 public:
  /**
   * @brief Constructor.
   *
   * @param activeSystem           The active system.
   * @param environmentSystems     A list of all the environment systems.
   */
  FDETask(std::shared_ptr<SystemController> activeSystem, std::vector<std::shared_ptr<SystemController>> environmentSystems);
  /**
   * @brief Default destructor.
   */
  virtual ~FDETask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(FDETaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    if (c.lcSettings.visitAsBlockSettings(v, blockname))
      return;
    if (c.loc.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the keyword is unknown.
    throw SerenityError((std::string) "Unknown block in FDETaskSettings: " + blockname);
  }

  /**
   * @brief The settings/keywords for the FDETask.
   */
  FDETaskSettings settings;

  /**
   * @brief Getter for the supersystem grid used in this run.
   *
   * This function should be called after the run() function is called.
   *
   * @return After a run, returns the supersystem grid (incl. cut off) if one is present.
   */
  std::shared_ptr<GridController> getSuperSystemGrid() {
    if (!_supersystemgrid)
      throw SerenityError("getSuperSystemGrid() was called before running the FDETask!");
    return _supersystemgrid;
  };
  /**
   * @brief Setter for a custom supersystem grid to use in this run.
   *
   * This function should be called before the run() function is called.
   *
   * @param grid The supersystem grid to use in this run.
   */
  void setSuperSystemGrid(std::shared_ptr<GridController> grid) {
    _supersystemgrid = grid;
  };

 private:
  std::shared_ptr<SystemController> _activeSystem;
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  std::shared_ptr<GridController> _supersystemgrid;
  std::shared_ptr<GridController> _finalGrid;
  void calculateNonAdditiveDispersionCorrection();
  void calculateMP2CorrelationContribution(std::shared_ptr<PotentialBundle<SCFMode>> fdePot);
  void calculateUnrelaxedMP2Density(std::shared_ptr<PotentialBundle<SCFMode>> fdePot,
                                    std::shared_ptr<EnergyComponentController> eCont);
  // prints the S*S values
  void printFaTAnalysis(std::vector<std::shared_ptr<SystemController>> systems,
                        std::shared_ptr<GridController> supersystemGrid = nullptr, bool actOnly = false);
};

} /* namespace Serenity */
#endif /* FDETASK_H */
