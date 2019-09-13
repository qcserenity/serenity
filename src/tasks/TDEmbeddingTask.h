/**
 * @file   TDEmbeddingTask.h
 *
 * @date   Apr 23, 2014
 * @author Thomas Dresselhaus
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
#ifndef TDEMBEDDINGTASK_H_
#define TDEMBEDDINGTASK_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "settings/Options.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
#include "math/Matrix.h"
#include "data/matrices/CoefficientMatrix.h"
#include "settings/EmbeddingSettings.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>


namespace Serenity {

using namespace Serenity::Reflection;
class SystemController;
class Settings;

struct TDEmbeddingTaskSettings {
  TDEmbeddingTaskSettings():
    locType(Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO),
    useEnvSys(false),
    orbitalThreshold(0.6),
    enforceCharges(false),
    noSupRec(true),
    truncationFactor(0.0),
    truncAlgorithm(Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE),
    netThreshold(1e-4),
    load(""),
    name(""),
    useFermiLevel(true)
  {}
  REFLECTABLE(
    (Options::ORBITAL_LOCALIZATION_ALGORITHMS) locType,
    (bool) useEnvSys,
    (double) orbitalThreshold,
    (bool) enforceCharges,
    (bool) noSupRec,
    (double) truncationFactor,
    (Options::BASIS_SET_TRUNCATION_ALGORITHMS) truncAlgorithm,
    (double) netThreshold,
    (std::string) load,
    (std::string) name,
    (bool) useFermiLevel
  )
public:
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
  TDEmbeddingTask(
      std::shared_ptr<SystemController>activeSystem,
      std::shared_ptr<SystemController> environmentSystem);
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
    if (!blockname.compare("")){
      visit_each(c, v);
    } else if(!c.embedding.visitSettings(v,blockname)){
      throw SerenityError((string)"Unknown block in TDEmbeddingTaskSettings: "+blockname);
    }
  }

  /**
   * @brief The settings/keywords for ProjectionBasedEmbTask:
   *        -locType : The localization method (see LocalizationTask, default: IBO)
   *        -useEnvSys : Relax the active system with respect to a previously calculated environment.
   *        -orbitalThreshold : Threshold for until localized orbitals get assigned to active system (regarding Mulliken population)
   *        -enforceCharges : If true, orbitals get assigned to match input charges
   *        -noSupRec : Only reconstructs subsystem potentials and evaluates supersystem potentials analytically
   *        -truncationFactor : The truncation factor used in the "primitive" truncation schemes
   *        -truncAlgorithm : The employed truncation algorithm.
   *        -netThreshold : The Mulliken net population threshold for the Mulliken net population truncation.
   *        -distantNonOrthogonal : True, if distant orbitals are not enforced to be orthogonal to the active orbitals.
   *        -nonOrthogonalCrit : The criterion to select non-orthogonal orbitals: NONE, DistantAtom
   *        -load : The path to the directory from which the supersystem is taken. If empty a supersystem calculation is done.
   *        -name : The name of the system if it is loaded from disk.
   *        -embeddingSettings: The embedding settings. See settings/EmbeddingSettings.h for details.
   */
  TDEmbeddingTaskSettings settings;

private:
  /// @brief The active system
  std::shared_ptr<SystemController>  _activeSystem;
  /// @brief The environment/embedding system
  std::shared_ptr<SystemController> _environmentSystem;
  /// @brief The indices of atoms which are not considered to be "distant" (optional)
  std::vector<unsigned int> _ghostIndices;
  /// @brief The supersystem indices of the active shells.
  std::vector<unsigned int> _activeShellIndices;
  /// @brief A bool for each atom in the system, which indicates whether it is "active" or not.
  std::vector<bool> _activeAtoms;

  /**
   * @brief Adjusts the basis set of the active system according to the choices in the settings.
   */
  void  setActiveSystemBasis();
  /**
   * @brief Creates a suffix for the naming of cube files.
   * @return "alpha,beta" if UNRESTRICTED, "" if RESTRICTED.
   */
  inline SpinPolarizedData<SCFMode, std::string> getDensityNameSuffix();
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
  inline double getFermiLevel(
      std::shared_ptr<SystemController> supersystem);
};

} /* namespace Serenity */

#endif /* PROJECTIONBASEDEMBTASK_H_ */
