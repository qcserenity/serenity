/**
 * @file ElectronicStructureCopyTask.h
 *
 * @author Moritz Bensberg
 * @date Feb 13, 2020
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

#ifndef TASKS_ELECTRONICSTRUCTURECOPYTASK_H_
#define TASKS_ELECTRONICSTRUCTURECOPYTASK_H_
/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h" //CoefficientMatrix definition.
#include "settings/Options.h"                //SCF_MODES
#include "settings/Reflection.h"             //Task settings.
#include "tasks/Task.h"                      //Inherits from.
/* Include Std and External Headers */
#include <Eigen/Dense> //Eigen::Vector3d.
#include <memory>      //smrt_ptr
#include <vector>      //std::vector.

namespace Serenity {
using namespace Serenity::Reflection;

/* Forward Declarations */
class SystemController;
class Geometry;
class BasisController;

struct ElectronicStructureCopyTaskSettings {
  ElectronicStructureCopyTaskSettings() : atomFrameIndices({}), orthogonalize(false), copyCharges(true) {
  }
  REFLECTABLE((std::vector<unsigned int>)atomFrameIndices, (bool)orthogonalize, (bool)copyCharges)
};

template<Options::SCF_MODES SCFMode>
class ElectronicStructureCopyTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param sourceSystem The systems which electronic structure should be copied.
   * @param targetSystems The systems to which the electronic structure is copied.
   */
  ElectronicStructureCopyTask(std::shared_ptr<SystemController> sourceSystem,
                              std::vector<std::shared_ptr<SystemController>> targetSystems);
  virtual ~ElectronicStructureCopyTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief The task settings.
   *   -- atomFrameIndices: Optional atom indices that may be prioritised in the frame construction.
   *   -- orthogonalize:    Orthogonalized the orbitals after rotation.
   *   -- copyCharges:      Adjust the charges/spins of the target systems to the charge/spin of the
   *                        source system.
   */
  ElectronicStructureCopyTaskSettings settings;

 private:
  /// @brief The active system
  std::shared_ptr<SystemController> _sourceSystem;
  /// @brief The environment/embedding system
  std::vector<std::shared_ptr<SystemController>> _targetSystems;
  ///@brief Construct the euler angles for rotating the source frame to the target frame.
  Eigen::Vector3d getEulerAngles(const std::vector<Eigen::Vector3d>& sourceFrame,
                                 const std::vector<Eigen::Vector3d>& targetFrame);
  ///@brief Construct an internal frame for the given geometry.
  std::vector<Eigen::Vector3d> getInternalFrame(std::shared_ptr<Geometry> geometry);
  ///@brief Rotate the coefficients from the source frame to the Cartesian frame by rotating the Cartesian axes.
  CoefficientMatrix<SCFMode> rotateMatrixIntoCartesianFrame(const std::vector<Eigen::Vector3d>& sourceFrame,
                                                            const CoefficientMatrix<SCFMode>& coefficients);
  ///@brief Rotate the coefficients from the Cartesian to the target frame by rotating the Cartesian axes.
  CoefficientMatrix<SCFMode> rotateMatrixIntoTargetFrame(const std::vector<Eigen::Vector3d>& targetFrame,
                                                         const CoefficientMatrix<SCFMode>& coefficients,
                                                         std::shared_ptr<BasisController> targetBasis);
  /**
   * @brief Check the source and target systems whether copying the electronic structure makes sense.
   *
   * Checked values:
   *   -- Spherical vs. Cartesian basis functions. Rotation is only allowed for purely spherical basis sets.
   *   -- Atom ordering. All atoms have to be ordered in the same way in order to ensure identical ordering of
   *      Basis functions.
   */
  void checkSystems();
  /**
   * @brief Orthogonalize the orbitals described via the coefficient matrix using a Löwdin orthogonalization.
   */
  void loewdingOrthogonalization(CoefficientMatrix<SCFMode>& coefficients, const MatrixInBasis<RESTRICTED>& S);
};

} /* namespace Serenity */

#endif /* TASKS_ELECTRONICSTRUCTURECOPYTASK_H_ */
