/**
 * @file   FDETask.h
 * @author Jan Unsleber, Kevin Klahr, Thomas Dresselhaus
 *
 * @date   last reworked on Nov 11. 2016
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
#ifndef FDETASK_H
#define	FDETASK_H

/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "data/grid/GridPotential.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>


namespace Serenity {
class SystemController;

using namespace Serenity::Reflection;

struct FDETaskSettings {
  FDETaskSettings():
    naddKinFunc(Options::KINFUNCTIONALS::PW91K),
    naddXCFunc(Options::XCFUNCTIONALS::PW91),
    dispersion(Options::DFT_DISPERSION_CORRECTIONS::NONE),
    embeddingMode(Options::KIN_EMBEDDING_MODES::NADD_FUNC),
    smoothFactor(0.0),
    potentialBasis(""),
    singValThreshold(0.0),
    lbDamping(0.995),
    lbCycles(0),
    carterCycles(0),
    locType(Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO),
    gridCutOff(-1.0),
    smallSupersystemGrid(false),
    finalGrid(true),
    truncateProjector(false),
    projecTruncThresh(1.0e+1),
    distantKinFunc(false),
    levelShiftParameter(1e+6),
    basisFunctionRatio(0.0),
    borderAtomThreshold(0.02){}
  REFLECTABLE(
      (Options::KINFUNCTIONALS) naddKinFunc,
      (Options::XCFUNCTIONALS) naddXCFunc,
      (Options::DFT_DISPERSION_CORRECTIONS) dispersion,
      (Options::KIN_EMBEDDING_MODES) embeddingMode,
      (double) smoothFactor,
      (std::string) potentialBasis,
      (double) singValThreshold,
      (double) lbDamping,
      (unsigned int) lbCycles,
      (unsigned int) carterCycles,
      (Options::ORBITAL_LOCALIZATION_ALGORITHMS) locType,
      (double) gridCutOff,
      (bool) smallSupersystemGrid,
      (bool) finalGrid,
      (bool) truncateProjector,
      (double) projecTruncThresh,
      (bool) distantKinFunc,
      (double) levelShiftParameter,
      (double) basisFunctionRatio,
      (double) borderAtomThreshold
  )

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
  FDETask(
      std::shared_ptr<SystemController> activeSystem,
      std::vector<std::shared_ptr<SystemController> > environmentSystems);
  /**
   * @brief Default destructor.
   */
  virtual ~FDETask() =default ;
  /**
   * @see Task
   */
  void run();

  /**
   * @brief The settings/keywords for the FDETask: \n
   *        -naddKinFunc: The non-additive kinetic energy functional (Default: TF)
   *        -naddXCFunc: The non-additive exchange-correlation functional (Default: BP86)
   *        -calculateGradients: If true, calculate FDE gradients (Default: false)
   *        -exactNaddKin: If true, use reconstructed non-additive kinetic FDE potential (ignores the naddKinFunc keyword, Default: false)
   *        -gridCutOff: Modifies the super-system grid used to evaluate the FDE potential.
   *                     Only atoms within gridCutOff distance (in Angstrom) from any atom of the active
   *                     system are used to construct the super-system grid (Default: 5.0).
   *        -truncateProjector: A flag whether the projector should be truncated (default: false). See HuzinagaFDEProjectionPotential.h.
   *        -projectionTruncThreshold: The projection truncation threshold (default: 1.0e+1). See HuzinagaFDEProjectionPotential.h.
   *        -distantKinFunc: A flag whether not projected subsystems should be treated with a non additive kin. energy func. See HuzinagaFDEProjectionPotential.h.
   *        -levelShiftParameter: The level-shift of environment orbitals.
   *        -basisFunctionRatio: The retained basis function ratio for hybrid-projection methods.
   *        -borderAtomThreshold: The mulliken population threshold for projecting an environment orbital.
   */
  FDETaskSettings settings;

  /**
   * @brief Getter for the supersystem grid used in this run.
   *
   * This function should be called after the run() function is called.
   *
   * @return After a run, returns the supersystem grid (incl. cut off) if one is present.
   */
  std::shared_ptr<GridController> getSuperSystemGrid(){
    assert(_supersystemgrid && "getSuperSystemGrid() was called before running the FDETask!");
    return _supersystemgrid;
  };
  /**
   * @brief Setter for a custom supersystem grid to use in this run.
   *
   * This function should be called before the run() function is called.
   *
   * @param grid The supersystem grid to use in this run.
   */
  void setSuperSystemGrid(std::shared_ptr<GridController> grid){
    _supersystemgrid=grid;
  };
  /**
   * @brief Getter for the grid on which the final energy is evaluated.
   * @return The final grid.
   */
  std::shared_ptr<GridController> getFinalGrid(){
    return _finalGrid;
  };
  /**
   * @brief Setter for the grid on which the final energy is evaluated.
   * @param grid The "final" grid controller.
   */
  void setFinalGrid(std::shared_ptr<GridController> grid){
    _finalGrid=grid;
  };
private:
  std::shared_ptr<SystemController>  _activeSystem;
  std::vector<std::shared_ptr<SystemController> > _environmentSystems;
  std::shared_ptr<GridController> _supersystemgrid;
  std::shared_ptr<GridController> _finalGrid;
};

} /* namespace Serenity */
#endif	/* FDETASK_H */
