/**
 * @file FCIDumpFileWriterTask.h
 *
 * @author Moritz Bensberg
 * @date Feb. 12, 2024
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

#ifndef SERENITY_FCIDUMPFILEWRITERTASK_H
#define SERENITY_FCIDUMPFILEWRITERTASK_H

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"     //Restricted/Unrestricted handling.
#include "misc/SerenityError.h"         //Errors.
#include "settings/EmbeddingSettings.h" //Embedding settings.
#include "settings/Reflection.h"        //Reflectable.
#include "tasks/Task.h"                 //Base class.
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {

/* Forward declarations */
class SystemController;

/**
 * @struct FCIDumpFileWriterTaskSettings
 * @brief The settings for the FCIDumpFileWriterTask.
 *
 * The following settings are available:
 * - mullikenThreshold:                 Prescreening threshold for the auxiliary function to orbital mapping.
 * - orbitalToShellThreshold:           Prescreening threshold for coefficient values.
 * - outputFilePath:                    The file path for the output file.
 * - orbitalRangeAlpha                  The orbital range for the alpha orbitals to be written to the FCI dump file.
 *                                      This is also the orbitals range, in the case restricted orbitals are used.
 * - orbitalRangeBeta                   The orbital range for the beta orbitals to be written to the FCI dump file.
 * - onlyValenceOrbitals                If true, and no orbital ranges are provided, the orbital ranges are assumed to
 * be the range of valence orbitals.
 * - valenceOrbitalsFromEnergyCutOff    If true, valence orbitals are determined by an energy cut off. If false,
 * tabulated values are used.
 * - energyCutOff                       Energy cut off to determine core orbitals.
 * - virtualEnergyCutOff                Energy cut off to determine virtual valence orbitals.
 * - integralSizeCutOff                 Integrals smaller than this value are not written to the output file.
 * - doiNetThreshold                    Prescreening threshold for the orbital coefficients during calculation of
 * differential overlap integrals.
 * - doiIntegralPrescreening            All orbital combinations ik/jl in the integral (ik|jl) are ignored during the
 * integral transformation that have a differential overlap smaller than this threshold.
 * - calculateCoreEnergy                If true, the energy of the core electrons is calculated and printed to the
 * output file.
 */
struct FCIDumpFileWriterTaskSettings {
  FCIDumpFileWriterTaskSettings()
    : mullikenThreshold(1e-4),
      orbitalToShellThreshold(1e-3),
      outputFilePath("./fcidump.txt"),
      orbitalRangeAlpha({}),
      orbitalRangeBeta({}),
      onlyValenceOrbitals(true),
      valenceOrbitalsFromEnergyCutOff(false),
      energyCutOff(-5.0),
      virtualEnergyCutOff(+1.0),
      integralSizeCutOff(1e-9),
      doiNetThreshold(1e-7),
      doiIntegralPrescreening(1e-7),
      calculateCoreEnergy(true) {
  }
  REFLECTABLE((double)mullikenThreshold, (double)orbitalToShellThreshold, (std::string)outputFilePath,
              (std::vector<unsigned int>)orbitalRangeAlpha, (std::vector<unsigned int>)orbitalRangeBeta,
              (bool)onlyValenceOrbitals, (bool)valenceOrbitalsFromEnergyCutOff, (double)energyCutOff,
              (double)virtualEnergyCutOff, (double)integralSizeCutOff, (double)doiNetThreshold,
              (double)doiIntegralPrescreening, (bool)calculateCoreEnergy)
 public:
  EmbeddingSettings embedding;
};

/**
 * @class FCIDumpFileWriterTask FCIDumpFileWriterTask.h
 * @brief Task to write FCI dump files. @see io/FCIDumpFileWriter.h for more information.
 */
template<Options::SCF_MODES SCFMode>
class FCIDumpFileWriterTask : public Task {
 public:
  /**
   * @brief Constructor.
   *
   * @param activeSystem           The active system.
   * @param environmentSystems     A list of all the environment systems.
   */
  FCIDumpFileWriterTask(std::shared_ptr<SystemController> activeSystem,
                        std::vector<std::shared_ptr<SystemController>> environmentSystems);
  /**
   * @brief Default destructor.
   */
  virtual ~FCIDumpFileWriterTask() = default;
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
  void visit(FCIDumpFileWriterTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the keyword is unknown.
    throw SerenityError((std::string) "Unknown block in FCIDumpFileWriterTaskSettings: " + blockname);
  }
  ///@brief The settings.
  FCIDumpFileWriterTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _activeSystem;
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> getOrbitalRanges();
};
} /* namespace Serenity */

#endif // SERENITY_FCIDUMPFILEWRITERTASK_H
