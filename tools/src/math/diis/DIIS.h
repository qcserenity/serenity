/**
 * @file   DIIS.h
 *
 * @date   Nov 18, 2013
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
#ifndef DIIS_H_
#define DIIS_H_
/* Include Serenity Internal Headers */
#include "math/optimizer/Optimizer.h"
#include "misc/VectorOnDiskStorageController.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward calculations */
/**
 * @class DIIS DIIS.h
 * @brief An implementation of the DIIS by Pulay.
 *
 * According to:
 *   [1] Chem. Phys. Lett. 73, 393 (1980) (original method)
 *   [2] J. Comput. Chem. 3, 556 (1982) (usage of FPS-SPF as error vector)
 *   [3] J. Chem. Phys. 84, 5728 (1986) (normalization and diagonal scaling)
 */
class DIIS {
 public:
  /**
   * @param maxStore                       How many error vectors and target matrices are stored at max.
   *                                       If reached (at least) the oldest ones are deleted in the
   *                                       following cycle.
   * @param diskMode                       Use of VectorOnDiskStorageController instead of conventional
   *                                       Eigen vectors.
   * @param scaling                        Scaling parameter for the diagonal of the B matrix; is used
   *                                       to suppress large coefficients when linear dependencies occur.
   */
  DIIS(unsigned maxStore, bool diskMode = false, double scaling = 1.0);

  ///@brief Default destructor.
  virtual ~DIIS() = default;

  /**
   * @brief Perform a DIIS step
   *
   * @param parameters  the target vector of the current cycle. It will be mutated to be
   *                    hopefully closer to the optimum as the incoming vector.
   * @param gradients   the error vector of the current optimization cycle. It is some measure
   *                    for the gradient.
   */
  void optimize(Eigen::Ref<Eigen::VectorXd> parameters, const Eigen::Ref<const Eigen::VectorXd>& gradients);

  /**
   * @brief Perform a DIIS step
   *
   * @param parameters  the target vector of the current cycle. It will be mutated to be
   *                    hopefully closer to the optimum as the incoming vector.
   * @param gradients   the error vector of the current optimization cycle. It is some measure
   *                    for the gradient.
   */
  void optimize(Eigen::MatrixXd& parameters, const Eigen::MatrixXd& gradients) {
    this->optimize(Eigen::Map<Eigen::VectorXd>(parameters.data(), parameters.cols() * parameters.rows()),
                   Eigen::Map<const Eigen::VectorXd>(gradients.data(), gradients.cols() * gradients.rows()));
  }

  /**
   * @brief Perform a DIIS step
   *
   * @param parameters  the target vector of the current cycle. It will be mutated to be
   *                    hopefully closer to the optimum as the incoming vector.
   * @param gradients   the error vector of the current optimization cycle. It is some measure
   *                    for the gradient.
   */
  void optimize(std::vector<double>& parameters, const std::vector<double>& gradients) {
    this->optimize(Eigen::Map<Eigen::VectorXd>(parameters.data(), parameters.size()),
                   Eigen::Map<const Eigen::VectorXd>(gradients.data(), gradients.size()));
  }
  /**
   * @brief Perform a DIIS step
   *
   * @param targetVector     the target vector of the current cycle. It will be mutated to be
   *                         hopefully closer to the optimum as the incoming vector.
   * @param newErrorVector   the error vector of the current optimization cycle. It is some measure
   *                         for the gradient.
   */
  void optimize(VectorOnDiskStorageController& targetVector, VectorOnDiskStorageController& newErrorVector);

  ///@brief (Re-)initializes this object. All stored data is erased.
  void reinit();
  /**
   * @brief Store a parameter vector without mutating it.
   * @param parameters  The parameters.
   * @param gradients   The errors.
   */
  void store(const Eigen::Ref<const Eigen::VectorXd>& parameters, const Eigen::Ref<const Eigen::VectorXd>& gradients);
  /**
   * @brief Store a parameter vector in a matrix representation without mutating it.
   * @param parameters  The parameters.
   * @param gradients   The errors.
   */
  void storeMatrix(const Eigen::MatrixXd& parameters, const Eigen::MatrixXd& gradients);
  /**
   * @brief Store a parameter vector without mutating it.
   * @param parameters  The parameters.
   * @param gradients   The errors.
   */
  void store(const std::vector<double>& parameters, const std::vector<double>& gradients);
  /**
   * @brief Store a parameter vector without mutating it.
   * @param parameters  The parameters.
   * @param gradients   The errors.
   */
  void store(VectorOnDiskStorageController& parameters, VectorOnDiskStorageController& gradients);
  /**
   * @brief Getter for the number of stored error/parameter vectors.
   * @return The number of stored vectors.
   */
  unsigned int getNVectorsStored();

 private:
  // Erase the oldest stuff and shift up the other entries in the vectors.
  void shiftVectors();

  // Initialize a new B matrix
  void initNewB();

  // Complete sweep.
  void cleanUp();

  // Checks whether the linear system is well conditioned.
  bool checkConditioning();

  // Excludes the oldest equation at that point from the B matrix construction.
  void excludeEquation();

  // How many error vectors and target matrices are stored at max.
  unsigned _maxStore;

  // The series of error vectors.
  std::vector<std::unique_ptr<Eigen::VectorXd>> _errorVectors;
  std::vector<std::unique_ptr<VectorOnDiskStorageController>> _errorDiskVectors;

  // The series of target vectors. Are linearly combined for the extrapolation.
  std::vector<std::unique_ptr<Eigen::VectorXd>> _targetVectors;
  std::vector<std::unique_ptr<VectorOnDiskStorageController>> _targetDiskVectors;

  // How many error vectors and target matrices are currently stored.
  unsigned _nStored;

  // Tracks how many vectors are to be excluded in in the current iteration.
  unsigned _nExcluded = 0;

  // Holds the scalar products of the error vectors.
  Eigen::MatrixXd _B;

  // Use of VectorOnDiskStorageController instead of regular vectors.
  bool _diskMode;

  // Parameter to dampen the coefficients. Is used to scale the diagonal of the B matrix.
  const double _scaling;

  // Keeps track of the iterations performed.
  unsigned _cycle;
};

} /* namespace Serenity */
#endif /* DIIS_H_ */
