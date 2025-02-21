/**
 * @file   VirtualOrbitalSpaceSelectionTask.h
 *
 * @date   Aug 5, 2019
 * @author Johannes Toelle
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
#ifndef VIRTUALORBITALSPACESELECTIONTASK_H_
#define VIRTUALORBITALSPACESELECTIONTASK_H_
/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "settings/EmbeddingSettings.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

using namespace Serenity::Reflection;

struct VirtualOrbitalSpaceSelectionTaskSettings {
  VirtualOrbitalSpaceSelectionTaskSettings()
    : excludeProjection(false),
      localCanonicalVirtuals(0.0),
      envCanonicalVirtuals(0.0),
      localizedVirtualorbitals(false),
      localizedEnvVirtualorbitals(false),
      recalculateFockMatrix(false),
      identifier(""),
      mixingOccAndVirtualOrbitals(false),
      relaxation(false),
      onlyOne(false) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::NONE;
    embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS ::NONE;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
  };
  REFLECTABLE((bool)excludeProjection, (double)localCanonicalVirtuals, (double)envCanonicalVirtuals,
              (bool)localizedVirtualorbitals, (bool)localizedEnvVirtualorbitals, (bool)recalculateFockMatrix,
              (std::string)identifier, (bool)mixingOccAndVirtualOrbitals, (bool)relaxation, (bool)onlyOne)
 public:
  EmbeddingSettings embedding;
};
/**
 * @class  VirtualOrbitalSpaceSelectionTask VirtualOrbitalSpaceSelectionTask.h
 * @brief  Performs different virtual orbital selection/localization schemes
 * Literature:
 * J. Chem. Phys., 153, 184113 (2020)
 */
template<Options::SCF_MODES SCFMode>
class VirtualOrbitalSpaceSelectionTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param activeSystems The active systems.
   * @param environmentSystems Environment systems.
   */
  VirtualOrbitalSpaceSelectionTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                                   const std::vector<std::shared_ptr<SystemController>>& passiveSystems);
  /**
   * @brief Default destructor.
   */
  virtual ~VirtualOrbitalSpaceSelectionTask() = default;
  /**
   * @see Task.h
   */
  void run();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(VirtualOrbitalSpaceSelectionTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the blockname is unknown.
    throw SerenityError((std::string) "Unknown block in VirtualOrbitalSpaceSelectionTaskSettings: " + blockname);
  }
  /**
   * @brief The settings/keywords for the LRSCFTask: \n
   *        -excludeProjection: Excludes artificially shifted envrionment in projection based embedding calculation
   *        -localCanonicalVirtuals: Select canonical virtual orbitals located on the active subsystem based on a
   * modified overlap criterion
   *        -envCanonicalVirtuals: Select canonical virtual orbitals located on the environment
   * subsystems based on a modified overlap criterion
   *        -localizedVirtualorbitals: Performs a localization of the virtual
   * orbitals and selects the virtual orbitals located on the subsystem
   *        -localizedEnvVirtualorbitals: Performs a
   * localization of the virtual orbitals and selects the virtual orbitals located on the subsystem
   *        -identifier:
   * Identifier for the file in which the new orbitals and orbital energies are stored (needed if the same system has
   * two sets for orbitals)
   *        -mixingOccAndVirtualOrbitals: Takes the occ orbitals of env subsystem 1 and the virtuals of
   * the environ subsystem 2 The structure of the supersystem needs to be given by the act system
   *        -relaxation: Performs an additional orbital space orthogonalization in case of mixingOccAndVirtualOrbitals
   */
  VirtualOrbitalSpaceSelectionTaskSettings settings;

  /**
   * @brief Writes the new CoefficientMatrix and orbital energies to disk
   */
  void writeOrbitalsToHDF5(CoefficientMatrix<SCFMode>& coefficients, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues,
                           std::shared_ptr<SystemController>& system);
  /**
   * @brief Updates the new coefficients and eigenvalues
   */
  void updateOrbitals(CoefficientMatrix<SCFMode>& coefs, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues,
                      SpinPolarizedData<SCFMode, std::vector<unsigned int>> indices);
  /**
   * @brief Prints the virtual orbitals with orbital energy
   */
  void printNewOrbitals(SpinPolarizedData<SCFMode, unsigned int>& nOcc, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues);
  ///@brief Check if input is correct
  void checkInput();

 private:
  // The systemcontroller of the active subsystems
  std::vector<std::shared_ptr<SystemController>> _act;
  // The systemcontroller of the environment subsystems
  std::vector<std::shared_ptr<SystemController>> _env;
};

} /* namespace Serenity */

#endif /* VIRTUALORBITALSPACESELECTIONTASK_H_ */