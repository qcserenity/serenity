/**
 * @file IncrementalFockMatrix.cpp
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

/* Include Class Header*/
#include "potentials/IncrementalFockMatrix.h"
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "io/FormattedOutputStream.h"
#include "misc/WarningTracker.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
IncrementalFockMatrix<SCFMode>::IncrementalFockMatrix(std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrixController,
                                                      const double prescreeningThreshold,
                                                      double prescreeningIncrementStart, double prescreeningIncrementEnd,
                                                      unsigned int incrementSteps, std::string outputString)
  : _densityMatrixController(densityMatrixController),
    _prescreeningThreshold(prescreeningThreshold),
    _prescreeningIncrementStart(prescreeningIncrementStart),
    _prescreeningIncrementEnd(prescreeningIncrementEnd != 0 ? prescreeningIncrementEnd : prescreeningThreshold),
    _incrementSteps(incrementSteps),
    _outputString(outputString),
    _alwaysFullBuild(incrementSteps == 0) {
  auto basis = densityMatrixController->getDensityMatrix().getBasisController();
  basis->addSensitiveObject(this->_self);
  initializeOldDensityMatrix();
  if (_incrementSteps == 0) {
    _incrementSteps = 1;
    WarningTracker::printWarning((std::string) " WARNING: Incremental steps must at least be 1, set incremental steps "
                                               "to 1 -> Rebuild Fock matrix "
                                               "every iteration!",
                                 true);
  }
}
template<Options::SCF_MODES SCFMode>
int IncrementalFockMatrix<SCFMode>::getCounter() {
  return _counter;
}
template<Options::SCF_MODES SCFMode>
unsigned int IncrementalFockMatrix<SCFMode>::getIncrementSteps() {
  return _incrementSteps;
}

template<Options::SCF_MODES SCFMode>
void IncrementalFockMatrix<SCFMode>::resetFockMatrix(std::vector<std::shared_ptr<FockMatrix<SCFMode>>> fs, double nextTreshold) {
  // No printing in the fist cycle.
  if (this->getCounter() != 0) {
    OutputControl::dOut << " ***** Reset Incremental Fock Matrix Build: " << _outputString << " *****" << std::endl;
    OutputControl::dOut << " ***** New Prescreening Threshold - " << nextTreshold << " ***** " << std::endl;
  }
  for (auto f_ptr : fs) {
    *f_ptr = FockMatrix<SCFMode>(_densityMatrixController->getDensityMatrix().getBasisController());
    auto& f = *f_ptr;
    for_spin(f) {
      f_spin.setZero();
    };
  }
}
template<Options::SCF_MODES SCFMode>
double IncrementalFockMatrix<SCFMode>::getCurrentThreshold() {
  double newThreshold = _prescreeningIncrementEnd;
  if (not _reachedFinalThreshold) {
    unsigned int nCycles = (_counter < 0) ? 0 : _counter / _incrementSteps;
    newThreshold = _prescreeningIncrementStart;
    for (unsigned int i = 0; i < nCycles; ++i)
      newThreshold /= 100;
    if (newThreshold <= _prescreeningIncrementEnd) {
      _reachedFinalThreshold = true;
      newThreshold = _prescreeningIncrementEnd;
    }
  }
  return newThreshold;
}

template<Options::SCF_MODES SCFMode>
bool IncrementalFockMatrix<SCFMode>::updateDensityAndThreshold(DensityMatrix<SCFMode>& p, double& threshold,
                                                               std::vector<std::shared_ptr<FockMatrix<SCFMode>>> fs) {
  if (!_oldDensityMatrix)
    initializeOldDensityMatrix();
  // Calculate current threshold and the density matrix change
  DensityMatrix<SCFMode> deltaP = _densityMatrixController->getDensityMatrix() - *_oldDensityMatrix;
  // Check for zero-changes.
  bool zeroChange = true;
  for_spin(deltaP) {
    double largestElement = deltaP_spin.array().abs().maxCoeff();
    if (largestElement > _prescreeningIncrementEnd * 0.1)
      zeroChange = false;
  };
  // do a full build
  //  - every _incrementSteps iterations
  //  - if the change is not zero
  //  - if the final threshold is not reached
  //  - if the alwaysFullBuild flag is set
  //  - if the basis set has changed
  bool fullBuild = (_counter % _incrementSteps == 0 && not zeroChange && not _reachedFinalThreshold) ||
                   _alwaysFullBuild || _basisChanged;
  double newThreshold = getCurrentThreshold();
  // Full build every few iterations for non-zero changes.
  // Prevent full construction due to zero/unitary change
  if (fullBuild) {
    // Full build. Reset Fock matrix, and assign default threshold and density matrix.
    p = _densityMatrixController->getDensityMatrix();
    this->resetFockMatrix(fs, newThreshold);
    threshold = _prescreeningThreshold;
  }
  else {
    // Incremental build. Set threshold and calculate density matrix increment.
    threshold = newThreshold;
    p = deltaP;
  }
  // reset old density matrix and increase the counter.
  *_oldDensityMatrix = _densityMatrixController->getDensityMatrix();
  if (not zeroChange)
    ++_counter;
  _basisChanged = false;
  return fullBuild;
}

template<Options::SCF_MODES SCFMode>
void IncrementalFockMatrix<SCFMode>::notify() {
  // the basis set has changed! Reset everything.
  _oldDensityMatrix = nullptr;
  _counter = 0;
  _basisChanged = true;
}

template<Options::SCF_MODES SCFMode>
void IncrementalFockMatrix<SCFMode>::initializeOldDensityMatrix() {
  _oldDensityMatrix = std::make_shared<DensityMatrix<SCFMode>>(_densityMatrixController->getDensityMatrix());
  DensityMatrix<SCFMode>& oldP = *_oldDensityMatrix;
  for_spin(oldP) {
    oldP_spin.setZero();
  };
}

template class IncrementalFockMatrix<Options::SCF_MODES::RESTRICTED>;
template class IncrementalFockMatrix<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
