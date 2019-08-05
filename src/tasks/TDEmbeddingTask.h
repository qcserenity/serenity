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
/* Include Std and External Headers */
#include <memory>
#include <vector>


namespace Serenity {

using namespace Serenity::Reflection;
class SystemController;
class Settings;

struct TDEmbeddingTaskSettings {
  TDEmbeddingTaskSettings():
    levelShiftParameter(1.0e6),
    locType(Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO),
    naddXCFunc(Options::XCFUNCTIONALS::BP86),
    naddKinFunc(Options::KINFUNCTIONALS::PW91K),
    longRangeNaddKinFunc(Options::KINFUNCTIONALS::NONE),
    embeddingMode(Options::KIN_EMBEDDING_MODES::LEVELSHIFT),
    dispersion(Options::DFT_DISPERSION_CORRECTIONS::NONE),
    useEnvSys(false),
    orbitalThreshold(0.6),
    enforceCharges(false),
    smoothFactor(0.0),
    potentialBasis(""),
    singValThreshold(0.0),
    lbDamping(0.995),
    lbCycles(0),
    carterCycles(0),
    noSupRec(true),
    truncationFactor(0.0),
    truncAlgorithm(Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE),
    netThreshold(1e-4),
    borderAtomThreshold(0.02),
    basisFunctionRatio(0.0),
    nonOrthogonalCrit(Options::NON_ORTHOGONAL_CRITERION::NONE),
    load(""),
    name("")
  {}
  REFLECTABLE(
    (double) levelShiftParameter,
    (Options::ORBITAL_LOCALIZATION_ALGORITHMS) locType,
    (Options::XCFUNCTIONALS) naddXCFunc,
    (Options::KINFUNCTIONALS) naddKinFunc,
    (Options::KINFUNCTIONALS) longRangeNaddKinFunc,
    (Options::KIN_EMBEDDING_MODES) embeddingMode,
    (Options::DFT_DISPERSION_CORRECTIONS) dispersion,
    (bool) useEnvSys,
    (double) orbitalThreshold,
    (bool) enforceCharges,
    (double) smoothFactor,
    (std::string) potentialBasis,
    (double) singValThreshold,
    (double) lbDamping,
    (unsigned int) lbCycles,
    (unsigned int) carterCycles,
    (bool) noSupRec,
    (double) truncationFactor,
    (Options::BASIS_SET_TRUNCATION_ALGORITHMS) truncAlgorithm,
    (double) netThreshold,
    (double) borderAtomThreshold,
    (double) basisFunctionRatio,
    (Options::NON_ORTHOGONAL_CRITERION) nonOrthogonalCrit,
    (std::string) load,
    (std::string) name
  )
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

  void run();

  /**
   * @brief The settings/keywords for ProjectionBasedEmbTask:
   *        -levelShiftParameter : Shift projected orbitals by levelShiftParameters (default: 1.0e-6 Eh)
   *        -locType : The localization method (see LocalizationTask, default: IBO)
   *        -naddXCFunc : The non-additive exchange correlation functional
   *        -naddKinFunc : The non-additive kinetic functional (used in for EMBEDDING_MODE::NADD_FUNC)
   *        -longRangeNaddKinFunc : The kinetic energy functional used to correct contributions from non orthogonal orbitals.
   *        -embeddingMode : The type of embedding to run (e.g., level shift, potential reconstruction... see Options class)
   *        -dispersion : The dispersion interaction between subsystems.
   *        -useEnvSys : Relax the active system with respect to a previously calculated environment.
   *        -orbitalThreshold : Threshold for until localized orbitals get assigned to active system (regarding Mulliken population)
   *        -enforceCharges : If true, orbitals get assigned to match input charges
   *        -smoothFactor : Smoothing to be used in potential reconstruction
   *        -potentialBasis : Basis to express the potential in during Wu-Yang reconstruction
   *        -singValThreshold : Threshold for singular value decomposition in Wu-Yang Newton-Raphson step
   *        -lbDamping : Damping to be used during the van Leeuwen-Baerends reconstruction
   *        -lbCycles : Maximum cycles for van Leeuwen-Baerends scheme
   *        -carterCycles : Maximum cycles for Zhang-Carter scheme
   *        -noSupRec : Only reconstructs subsystem potentials and evaluates supersystem potentials analytically
   *        -truncationFactor : The truncation factor used in the "primitive" truncation schemes
   *        -truncAlgorithm : The employed truncation algorithm.
   *        -netThreshold : The Mulliken net population threshold for the Mulliken net population truncation.
   *        -distantNonOrthogonal : True, if distant orbitals are not enforced to be orthogonal to the active orbitals.
   *        -borderAtomThreshold : The Mulliken population threshold used to determine if an orbital is considered "distant" or not.
   *                               The Mulliken population of the orbital on all not "distant" atoms has to exceed this threshold in order
   *                               to be included in the projector.
   *        -basisFunctionRatio : The minimum ratio of retained basis functions needed in order to consider a atom to be not "distant".
   *        -nonOrthogonalCrit : The criterion to select non-orthogonal orbitals: NONE, DistantAtom
   *        -load : The path to the directory from which the supersystem is taken. If empty a supersystem calculation is done.
   *        -name : The name of the system if it is loaded from disk.
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
};

} /* namespace Serenity */

#endif /* PROJECTIONBASEDEMBTASK_H_ */
