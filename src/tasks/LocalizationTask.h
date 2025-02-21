/**
 * @file   LocalizationTask.h
 *
 * @date   Apr 22, 2014
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
#ifndef LOCALIZATIONTASK_H_
#define LOCALIZATIONTASK_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "settings/LocalizationOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace Serenity::Reflection;
struct LocalizationTaskSettings {
  LocalizationTaskSettings()
    : locType(Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO),
      maxSweeps(100),
      alignExponent(4),
      useKineticAlign(false),
      splitValenceAndCore(false),
      useEnergyCutOff(true),
      energyCutOff(-5.0),
      nCoreOrbitals(std::numeric_limits<unsigned int>::infinity()),
      localizeVirtuals(false),
      splitVirtuals(true),
      virtualEnergyCutOff(+1.0),
      nRydbergOrbitals(std::numeric_limits<unsigned int>::infinity()),
      replaceVirtuals(false),
      separateSOMOs(false){};
  REFLECTABLE((Options::ORBITAL_LOCALIZATION_ALGORITHMS)locType, (unsigned int)maxSweeps, (unsigned int)alignExponent,
              (bool)useKineticAlign, (bool)splitValenceAndCore, (bool)useEnergyCutOff, (double)energyCutOff,
              (unsigned int)nCoreOrbitals, (bool)localizeVirtuals, (bool)splitVirtuals, (double)virtualEnergyCutOff,
              (unsigned int)nRydbergOrbitals, (bool)replaceVirtuals, (bool)separateSOMOs)
 public:
  /**
   * @brief Parse the settings from the input an instance of this class.
   * @param c The settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("LOC")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};
/**
 * @class  LocalizationTask LocalizationTask.h
 * @brief  Localize orbitals.
 *
 * In general that means: perform unitary transformations among the occupied orbitals
 * to optimize some (localization) criterion.
 */
class LocalizationTask : public Task {
 public:
  /**
   * @param system The system from which orbitals to localize are taken.
   * @param templateSystem A template to which the system is aligned.
   */
  LocalizationTask(std::shared_ptr<SystemController> systemController,
                   std::vector<std::shared_ptr<SystemController>> templateSystem = {});
  /**
   * @brief Default destructor.
   */
  virtual ~LocalizationTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief The settings/keywords for LocalizationTask:
   *        - locType: The localization algorithm. The following algorithm can be chosen:
   *          - IBO : Intrinsic Bond Orbitals (default)
   *          - PM : Pipek-Mezey
   *          - FB : Foster-Boys
   *          - ER : Edminston-Ruedenberg
   *          - ALIGN : Orbital alignment (needs environment/template system)
   *        - maxSweeps : Maximum number of micro-iterations (default: 1000)
   *        - alignExponent : Exponent used in the orbital alignment.
   *        - useKineticAlign : Use the kinetic energy as an additional criterion for the orbital alignment.
   *        - splitValenceAndCore : Localize valence and core separately.
   *        - useEnergyCutOff : Use an energy cut-off to determine core orbitals.
   *        - energyCutOff : Orbital eigenvalue threshold to select core orbitals.
   *        - nCoreOrbitals : Use a predefined number of core orbitals.
   *                          Needs useEnergyCutOff = false and splitValenceAndCore = true.
   *        - localizeVirtuals : If true, the virtual orbitals are localized as well. Default false.
   *                             This is only supported for the IBO and ALIGN schemes.
   *        - splitVirtuals : If true, the valence virtuals and diffuse virtuals are localized separately.
   *                          By default true.
   *        - virtualEnergyCutOff : Orbital eigenvalue threshold to select diffuse virtual orbitals. By default 1.0.
   *        - nRydbergOrbitals : Use a predefined number of diffuse virtual orbitals. Not used by default.
   *        - replaceVirtuals : Reconstruct the virtual orbitals before localization by projecting all occupied orbitals
   *                            and cleanly separating valence virtuals from diffuse virtuals. This should be used if
   *                            the IBO or ALIGN approaches are chosen and the same orbital set was not already
   *                            reconstructed in this manner before.
   *        - separateSOMOs : Run the orbital localization separately for singly occupied orbitals. This only makes
   *                          sense if the orbitals are (quasi) restricted.
   */
  LocalizationTaskSettings settings;

 private:
  template<Options::SCF_MODES SCFMode>
  void runByLastSCFMode();
  /**
   * @brief Split valence and core orbitals.
   * @return The indices of the valence and the core orbitals.
   */
  template<Options::SCF_MODES SCFMode>
  std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>>
  getValenceOrbitalIndices();
  /**
   * @brief Split virtual valence and rydberg orbitals.
   * @return The indices of the valence and the rydberg orbitals.
   */
  template<Options::SCF_MODES SCFMode>
  std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>>
  getVirtualValenceOrbitalIndices();
  /**
   * @brief Getter for spin polarized output praefix string.
   * @return The strings.
   */
  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode, std::string> getOutputPraefix();

  template<Options::SCF_MODES SCFMode>
  std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>>
  separateSOMOs(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& valenceOrbitalRange);

  template<Options::SCF_MODES SCFMode>
  std::vector<SpinPolarizedData<SCFMode, std::vector<unsigned int>>> separateOrbitalRanges();

  const std::shared_ptr<SystemController> _systemController;
  const std::vector<std::shared_ptr<SystemController>> _templateSystem;
};

} /* namespace Serenity */

#endif /* LOCALIZATIONTASK_H_ */
