/**
 * @file IncrementalFockMatrix.h
 *
 * @date   Jul 25, 2020
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

#ifndef POTENTIALS_INCREMENTALFOCKMATRIX_H_
#define POTENTIALS_INCREMENTALFOCKMATRIX_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
/**
 * @class IncrementalFockMatrix
 * @brief A class that handels the logic for incremental Fock matrix constructions.
 *        The density matrix increment and the associated prescreening threshold
 *        are selected automatically. Zero-changes in the density matrix due to
 *        unitary transformation of the occupied orbitals is checked for.
 *        The incremental Fock matrix build can be disabled by selecting 0 as
 *        incrementSteps.
 */
template<Options::SCF_MODES SCFMode>
class IncrementalFockMatrix : public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @brief Constructor.
   * @param densityMatrixController     The density matrix controller.
   * @param prescreeningThreshold       The presecreening threshold for full matrix constructions.
   * @param prescreeningIncrementStart  The final prescreening threshold for incremental constructions.
   * @param prescreeningIncrementEnd    The final prescreening threshold for incremental constructions.
   * @param incrementSteps              The interval in which the full matrix is constructed.
   * @param outputString                Additional output information on matrix reset.
   */
  IncrementalFockMatrix(std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrixController,
                        const double prescreeningThreshold, double prescreeningIncrementStart,
                        double prescreeningIncrementEnd, unsigned int incrementSteps, std::string outputString = "");
  /**
   * @brief Default destructor.
   */
  virtual ~IncrementalFockMatrix() = default;
  /**
   * @brief Assign density matrix (or its increment) and the associated threshold.
   * @param p            The density matrix. To be manipulated.
   * @param threshold    The threshold. To be manipulated.
   * @param f            The Fock matrix. May be set to zero.
   * @return True if the Fock matrix is reset and a full build is required.
   */
  bool updateDensityAndThreshold(DensityMatrix<SCFMode>& p, double& threshold,
                                 std::vector<std::shared_ptr<FockMatrix<SCFMode>>> fs);
  /**
   * @brief Reset the Fock matrix.
   * @param f              The Fock matrix to be resetted.
   * @param nextTreshold   The next threshold used for incremental construction.
   *
   * This function may be overwritten by derived classes.
   */
  void resetFockMatrix(std::vector<std::shared_ptr<FockMatrix<SCFMode>>> fs, double nextTreshold);
  /**
   * @brief Getter for the increment counter.
   * @return The increment counter.
   */
  int getCounter();
  /**
   * @brief  Getter for full construction interval.
   * @return The interval.
   */
  unsigned int getIncrementSteps();
  /**
   * @brief IncrementalFockMatrix is sensitive to any change in the basis set.
   */
  void notify() override final;
  /**
   * @brief Getter for the prescreening threshold for the full Fock-matrix construction.
   * @return The prescreening threshold.
   */
  double getPrescreeningThreshold() {
    return _prescreeningThreshold;
  };

 private:
  // The density matrix controller.
  std::shared_ptr<DensityMatrixController<SCFMode>> _densityMatrixController;
  // The presecreening threshold for full matrix constructions.
  double _prescreeningThreshold;
  // The initial prescreening threshold for incremental constructions.
  double _prescreeningIncrementStart;
  // The final prescreening threshold for incremental constructions.
  double _prescreeningIncrementEnd;
  // The interval in which the full matrix is constructed.
  unsigned int _incrementSteps;
  // Additional output information on matrix reset.
  std::string _outputString;
  // Flag to disable incremental Fock Matrix builds.
  bool _alwaysFullBuild = false;
  // Start with -1 in order to "hide" the first hidden SCF-step.
  int _counter = -1;
  // The basis changed since the last notfy call.
  bool _basisChanged;
  // The old density matrix.
  std::shared_ptr<DensityMatrix<SCFMode>> _oldDensityMatrix;
  // Stop any complete rebuild of the Fock matrix after reaching the final threshold.
  bool _reachedFinalThreshold = false;
  // Getter for the currently used prescreening threshold.
  double getCurrentThreshold();
  // Initialize the density matrix.
  void initializeOldDensityMatrix();
};

} /* namespace Serenity */

#endif /* POTENTIALS_INCREMENTALFOCKMATRIX_H_ */
