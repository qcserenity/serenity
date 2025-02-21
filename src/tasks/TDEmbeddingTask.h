/**
 * @file   TDEmbeddingTask.h
 *
 * @date   Apr 23, 2014
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
#ifndef TDEMBEDDINGTASK_H_
#define TDEMBEDDINGTASK_H_

/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "math/Matrix.h"
#include "misc/SerenityError.h"
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "settings/EmbeddingSettings.h"
#include "settings/LocalizationOptions.h"
#include "settings/MiscOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {

using namespace Serenity::Reflection;
class SystemController;
struct Settings;

struct TDEmbeddingTaskSettings {
  TDEmbeddingTaskSettings()
    : locType(Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO),
      orbitalThreshold(0.6),
      noSupRec(true),
      truncationFactor(0.0),
      truncAlgorithm(Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE),
      netThreshold(1e-4),
      load(""),
      name(""),
      maxResidual(1e-5),
      maxCycles(100),
      useFermiLevel(true),
      systemPartitioning(Options::SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD),
      mp2Type(Options::MP2_TYPES::LOCAL),
      splitValenceAndCore(false),
      addOrbitals(false) {
    lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
    lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
  }
  REFLECTABLE((Options::ORBITAL_LOCALIZATION_ALGORITHMS)locType, (double)orbitalThreshold, (bool)noSupRec,
              (double)truncationFactor, (Options::BASIS_SET_TRUNCATION_ALGORITHMS)truncAlgorithm, (double)netThreshold,
              (std::string)load, (std::string)name, (double)maxResidual, (unsigned int)maxCycles, (bool)useFermiLevel,
              (Options::SYSTEM_SPLITTING_ALGORITHM)systemPartitioning, (Options::MP2_TYPES)mp2Type,
              (bool)splitValenceAndCore, (bool)addOrbitals)
 public:
  LocalCorrelationSettings lcSettings;
  EmbeddingSettings embedding;
};

/**
 * @class  TDEmbeddingTask TDEmbeddingTask.h
 * @brief  An embedding scheme which starts from a solution for the supersystem.
 *
 * This task can use all embedding schemes listed in Options::KIN_EMBEDDING_MODES.\n\n
 *
 * Special case of projection-based embedding:\n
 * Using a projection operator orbitals of the active subsystem are forced to stay orthonormal
 * to orbitals of the environment. To get orbitals in the first place a calculation is performed
 * on the supersystem (e.g. a DFT calculation) followed by a localization method. Then another
 * calculation (e.g. WFT) on the active subsystem can be performed in a similar fashion to FDE,
 * i.e. the core Hamiltonian is modified by adding the Coulombic contributions of the environment
 * (exact) and the non-additive exchange-correlation contribution (calculated with a DFT
 * functional). Instead of the normally added non-additive kinetic energy potential a level shift
 * or a Huzinaga/Hoffmann's projection operator is added.
 */
template<Options::SCF_MODES SCFMode>
class TDEmbeddingTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param activeSystems      A vector containing all active systems (in the form of SystemControllers).
   * @param environmentSystems A vector containing all environment systems
   *                           (in the form of SystemControllers).
   */
  TDEmbeddingTask(std::shared_ptr<SystemController> activeSystem, std::shared_ptr<SystemController> environmentSystem);
  /**
   * @brief Default destructor.
   */
  virtual ~TDEmbeddingTask() = default;
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
  void visit(TDEmbeddingTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    if (c.lcSettings.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the blockname is unknown.
    throw SerenityError((std::string) "Unknown block in TDEmbeddingTaskSettings: " + blockname);
  }

  /**
   * @brief The settings/keywords for ProjectionBasedEmbTask:
   *        -locType : The localization method (see LocalizationTask, default: IBO)
   *        -orbitalThreshold : Threshold for until localized orbitals get assigned to active system (regarding Mulliken
   * population) -noSupRec : Only reconstructs subsystem potentials and evaluates supersystem potentials analytically
   *        -truncationFactor : The truncation factor used in the "primitive" truncation schemes
   *        -truncAlgorithm : The employed truncation algorithm.
   *        -netThreshold : The Mulliken net population threshold for the Mulliken net population truncation.
   *        -distantNonOrthogonal : True, if distant orbitals are not enforced to be orthogonal to the active orbitals.
   *        -nonOrthogonalCrit : The criterion to select non-orthogonal orbitals: NONE, DistantAtom
   *        -load : The path to the directory from which the supersystem is taken. If empty a supersystem calculation is
   * done. -name : The name of the system if it is loaded from disk. -maxResidual: Maximum residual/convergence
   * threshold for local MP2. -maxCycles: Maximum number of amplitude optimization cycles for local MP2. -useFermiLevel:
   * Use the Fermi level of the supersytem for the shift in the Fermi-shifted Huzinaga equation. -embeddingSettings: The
   * embedding settings. See settings/EmbeddingSettings.h for details. -lcSettings: Local correlation settings for local
   * MP2. -systemPartitioning: The system partitioning algorithm. -mp2Type: The type of MP2 used for the double hybrid
   * part. -splitValenceAndCore: Split valence and core orbitals during orbital localizations. -addOrbitals: Construct
   * the supersystem electronic structure from the subsystem electronic structures.
   */
  TDEmbeddingTaskSettings settings;

 private:
  /// @brief The active system
  std::shared_ptr<SystemController> _activeSystem;
  /// @brief The environment/embedding system
  std::shared_ptr<SystemController> _environmentSystem;
  /// @brief The indices of atoms which are not considered to be "distant" (optional)
  std::vector<unsigned int> _ghostIndices;
  /// @brief The supersystem indices of the active shells.
  std::vector<unsigned int> _activeShellIndices;
  /// @brief A bool for each atom in the system, which indicates whether it is "active" or not.
  std::vector<bool> _activeAtoms;

  /**
   * @brief Check the settings for logical errors. Nonsense is still allowed!
   */
  inline void checkInput();
  /**
   * @return Set up the supersystem from a supersystem calculation or fragment addition.
   */
  inline std::shared_ptr<SystemController> setUpSupersystem();
  /**
   * @param The supersystem controller.
   * @return The fermi-level of the supersystem or the largest occupied fock matrix element.
   */
  inline double getFermiLevel(std::shared_ptr<SystemController> supersystem);
  /**
   * @brief Run the FDE-like step for top-down potential reconstruction.
   * @param supersystem The supersystem controller.
   */
  inline void runTDPotentialReconstruction(std::shared_ptr<SystemController> supersystem);
};

} /* namespace Serenity */

#endif /* PROJECTIONBASEDEMBTASK_H_ */
