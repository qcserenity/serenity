/**
 * @file FaTConvergenceAccelerator.h
 *
 * @date Mar 29, 2018
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

#ifndef MISC_FATCONVERGENCEACCELERATOR_H_
#define MISC_FATCONVERGENCEACCELERATOR_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/matrices/FockMatrix.h"
#include "misc/VectorOnDiskStorageController.h"
#include "scf/damper/Damper.h"
#include "settings/Options.h"
#include "tasks/FreezeAndThawTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

/* Forward Declarations */
class SystemController;
class DIIS;

/**
 * @class FaTConvergenceAccelerator FaTConvergenceAccelerator.h
 * @brief Manages the convergence acceleration methods for Freeze-and-Thaw (FaT) calculations.
 *
 * Currently implemented:\n
 *      - DIIS\n
 * \n
 *
 * DIIS:\n
 *   Perform DIIS procedure on top of the FaT calculation for the density matrices.
 *   The supersystem density matrix for a FaT calculation is given by a block diagonal matrix with the density matrices
 *   of the subsystems as the diagonal blocks. The error associated with a not fully relaxed density matrix is given by
 * the commutator \f$[F,P]\f$ which is a block diagonal matrix as well.\n Both matrices are mapped to a vector and
 * optimized using a standard DIIS (math/diis/DIIS.h) and then mapped back to their matrix formulation before the
 * density matrices of the subsystems are updated. For more information about the DIIS procedure see math/diis/DIIS.h.\n
 *   According to: Journal of Computational Chemistry 3 (4): 556–560; doi:10.1002/jcc.540030413
 */
template<Options::SCF_MODES SCFMode>
class FaTConvergenceAccelerator {
 public:
  /**
   * @brief Constructor.
   * @param maxStore The maximum number of density vectors, which are stored in the DIIS.
   * @param settings The FreezeAndThawTaskSettings, which are used for in the FaT calculations. Contains all settings
   * for the DIIS/Damping.
   * @param activeSystems The vector of the active systems.
   * @param environmentSystems The vector of the passive/pure environment systems.
   */
  FaTConvergenceAccelerator(unsigned int maxStore, const FreezeAndThawTaskSettings& settings,
                            std::vector<std::shared_ptr<SystemController>> activeSystems,
                            std::vector<std::shared_ptr<SystemController>> environmentSystems);

  /**
   * @brief Performs a DIIS step or damps.
   * @param energy The energy associated with the current density.
   */
  void accelerateConvergence();

  /**
   * @brief Default destructor.
   */
  virtual ~FaTConvergenceAccelerator() = default;

 private:
  ///@brief The vector of the active systems.
  std::vector<std::shared_ptr<SystemController>> _activeSystems;
  ///@brief The vector of the passive/pure environment systems.
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  ///@brief The entries of the block diagonal of [F,P].
  std::shared_ptr<VectorOnDiskStorageController> _fpsMinusSPF;
  ///@brief The entries of the block diagonal of P.
  std::shared_ptr<VectorOnDiskStorageController> _densityVector;
  ///@brief The old entries of the block diagonal of P.
  std::shared_ptr<VectorOnDiskStorageController> _oldDensityVector;
  ///@brief The FDETaskSettings, which are used for in the FaT calculations.
  const FreezeAndThawTaskSettings _settings;
  ///@brief The DIIS.
  std::shared_ptr<DIIS> _diis;
  ///@brief The cycle.
  unsigned int _cycle = 0;

  /* Helper Functions */
  /**
   * @brief Calculates [F,P] and saves it in _fpsMinusSPF.
   */
  void calcFPSminusSPF();
  /**
   * @brief Calculates F for the system with index "activeSystemIndex".
   * @param activeSystemIndex The index of the subsystem.
   * @return The embedded fock matrix for the subsystem.
   */
  FockMatrix<SCFMode> calcEmbeddedFockMatrix(unsigned int activeSystemIndex);

  /**
   * @brief Calculates the RMSD of the supersystem density matrix entries.
   * @return The max. RMSD.
   */
  double calcRMSDofDensity();
};

} /* namespace Serenity */

#endif /* MISC_FATCONVERGENCEACCELERATOR_H_ */
