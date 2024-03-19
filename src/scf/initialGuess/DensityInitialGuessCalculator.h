/**
 * @file   DensityInitialGuessCalculator.h
 * @author Thomas Dresselhaus
 *
 * @date   11. Juli 2014, 17:22
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
#ifndef DENSITYINITIALGUESSCALCULATOR_H
#define DENSITYINITIALGUESSCALCULATOR_H
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "settings/ElectronicStructureOptions.h"

namespace Serenity {
/* Forward declarations */
class SystemController;

/**
 * @class DensityInitialGuessCalculator DensityInitialGuessCalculator.h
 * @brief Initial guess orbitals are here created with a guessed density instead of a Fock matrix.
 *
 * Most initial guess strategies create the starting orbitals via diagonalization of a Fock matrix,
 * just as is done during the SCF procedure. The Fock matrix itself can be an approximate one (as
 * e.g. the HCoreGuess, which neglects the two-electron terms) or constructed using a guessed
 * density. Implementations of this interface are classes using the latter strategy. The guessed
 * density is used to construct an initial Fock matrix and corresponding orbitals.
 */
class DensityInitialGuessCalculator : public InitialGuessCalculator<Options::SCF_MODES::RESTRICTED> {
 public:
  /**
   * @brief Default constructor.
   */
  DensityInitialGuessCalculator() = default;
  /**
   * @brief Default destructor.
   */
  virtual ~DensityInitialGuessCalculator() = default;
  /**
   * @brief Creates a guessed (approximate) density matrix for the given system.
   *
   * Must be overridden in an actual implementation.
   *
   * @param   system for which a density is guessed.
   * @returns The guessed density matrix.
   */
  virtual std::unique_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>>
  calculateInitialDensity(std::shared_ptr<SystemController> systemController) = 0;

  std::unique_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>
  calculateInitialGuess(const std::shared_ptr<SystemController> systemController) override final;
};

} /* namespace Serenity */
#endif /* DENSITYINITIALGUESSCALCULATOR_H */
