/**
 * @file DFTEmbeddedLocalCorrelationTask.cpp
 *
 * @date Apr 4, 2021
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
#include "DFTEmbeddedLocalCorrelationTask.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"        //Set up supersystem.
#include "io/FormattedOutputStream.h" //Filtered output streams.
#include "misc/SerenityError.h"       //Errors
#include "settings/Settings.h"        //Set up supersystem, name assignment.
#include "system/SystemController.h"
#include "tasks/FDETask.h"                    //Run FDE.
#include "tasks/FreezeAndThawTask.h"          //Run FaT.
#include "tasks/LocalCorrelationTask.h"       //Run local correlation calculation.
#include "tasks/ScfTask.h"                    //Set up supersystem.
#include "tasks/TopDownStaticEmbeddingTask.h" // SetUpSubsystems call.

namespace Serenity {

DFTEmbeddedLocalCorrelationTask::DFTEmbeddedLocalCorrelationTask(std::shared_ptr<SystemController> activeSystem,
                                                                 std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                                                 std::shared_ptr<SystemController> supersystem)
  : _activeSystem(activeSystem), _environmentSystems(environmentSystems), _supersystem(supersystem) {
}
DFTEmbeddedLocalCorrelationTask::~DFTEmbeddedLocalCorrelationTask() = default;

void DFTEmbeddedLocalCorrelationTask::run() {
  printSectionTitle("DFT-Embedded Local Correlation Calculation");
  std::string embeddingMode, naddXCFunc, naddKinFunc;
  Options::resolve<Options::KIN_EMBEDDING_MODES>(embeddingMode, this->settings.lcSettings.embeddingSettings.embeddingMode);
  Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(naddXCFunc, this->settings.lcSettings.embeddingSettings.naddXCFunc);
  Options::resolve<CompositeFunctionals::KINFUNCTIONALS>(naddKinFunc, this->settings.lcSettings.embeddingSettings.naddKinFunc);
  OutputControl::nOut << "  embedding mode:  " << embeddingMode << std::endl;
  OutputControl::nOut << "  naddXCFunc:      " << naddXCFunc << std::endl;
  if (this->settings.lcSettings.embeddingSettings.embeddingMode == Options::KIN_EMBEDDING_MODES::NADD_FUNC)
    OutputControl::nOut << "  naddKinFunc:     " << naddKinFunc << std::endl;
  // Set up subsystems if necessary
  if (this->settings.fromSupersystem)
    TopDownStaticEmbeddingTask<RESTRICTED>::setUpSubsystems(_supersystem, {_activeSystem}, _environmentSystems,
                                                            this->settings.loc, this->settings.split,
                                                            this->settings.add, this->settings.trunc);
  // Run FDE or FaT
  printSubSectionTitle("Embedded SCF Calculation");
  if (settings.runFaT) {
    std::vector<std::shared_ptr<SystemController>> allSystems = _environmentSystems;
    allSystems.insert(allSystems.begin(), _activeSystem);
    FreezeAndThawTask<RESTRICTED> fatTask(allSystems);
    fatTask.settings.embedding = this->settings.lcSettings.embeddingSettings;
    fatTask.settings.lcSettings = this->settings.lcSettings;
    fatTask.run();
  }
  FDETask<RESTRICTED> fdeTask(_activeSystem, _environmentSystems);
  fdeTask.settings.embedding = this->settings.lcSettings.embeddingSettings;
  fdeTask.settings.lcSettings = this->settings.lcSettings;
  fdeTask.settings.loc = this->settings.loc;
  fdeTask.settings.calculateEnvironmentEnergy = true;
  fdeTask.run();
  // Run localCorrelation task
  LocalCorrelationTask localCorrelationTask(_activeSystem, _environmentSystems);
  localCorrelationTask.settings.lcSettings = this->settings.lcSettings;
  if (this->settings.trunc.truncAlgorithm != Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE)
    localCorrelationTask.settings.lcSettings.useProjectedOccupiedOrbitals = true;
  localCorrelationTask.settings.normThreshold = this->settings.normThreshold;
  localCorrelationTask.settings.maxCycles = this->settings.maxCycles;
  localCorrelationTask.settings.writePairEnergies = this->settings.writePairEnergies;
  localCorrelationTask.settings.loc = this->settings.loc;
  localCorrelationTask.run();
}

void DFTEmbeddedLocalCorrelationTask::setUpSubsystems() {
  if (!_supersystem) {
    if (_environmentSystems.size() < 1)
      throw SerenityError("ERROR: The task settings suggest that an environment system should be available.\n"
                          "       However, none is given. This is probably an input error.");
    Settings supersystemSettings = _environmentSystems[0]->getSettings();
    supersystemSettings.name = "TMP_Supersystem";
    supersystemSettings.path = supersystemSettings.path + supersystemSettings.name + "/";
    supersystemSettings.charge = 0;
    supersystemSettings.spin = 0;
    _supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersystemSettings);
  }
  // Addition
  std::vector<std::shared_ptr<SystemController>> allSystems = _environmentSystems;
  allSystems.insert(allSystems.begin(), _activeSystem);
  SystemAdditionTask<RESTRICTED> additionTask(_supersystem, allSystems);
  additionTask.settings = this->settings.add;
  additionTask.run();
  // Scf
  printSubSectionTitle("Supersystem SCF Calculation");
  ScfTask<RESTRICTED> scfTask(_supersystem);
  scfTask.run();
  // Localization
  LocalizationTask localizationTask(_supersystem);
  localizationTask.settings = this->settings.loc;
  localizationTask.run();
  // Splitting
  SystemSplittingTask<RESTRICTED> splittingTask(_supersystem, allSystems);
  splittingTask.settings = this->settings.split;
  splittingTask.run();
  // Basis set truncation
  if (this->settings.trunc.truncAlgorithm != Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE) {
    BasisSetTruncationTask<RESTRICTED> truncationTask(_activeSystem);
    truncationTask.settings = this->settings.trunc;
    truncationTask.run();
  }
}

} /* namespace Serenity */
