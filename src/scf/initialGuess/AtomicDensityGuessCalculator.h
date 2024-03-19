/**
 * @file AtomicDensityGuessCalculator.h
 *
 * @date Jul 12, 2014
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

#ifndef ATOMICDENSITYGUESSCALCULATOR_H
#define ATOMICDENSITYGUESSCALCULATOR_H
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "scf/initialGuess/DensityInitialGuessCalculator.h"
/* Include Std and External Headers */
#include <memory> //smart ptr.

namespace Serenity {
/* Forward declarations */
class Atom;
struct Settings;

enum class GUESSMODES { SCF = 0, SCF_INPLACE = 12 };

/**
 * @class AtomicDensityGuessCalculator AtomicDensityGuessCalculator.h
 * @brief Calculates an initial guess based on atomic densities
 *
 * A spherical distribution of electrons around the atoms is assumed. The density in the basis
 * functions of the atoms is assumed to be diagonal and independent of the other atoms. With the
 * resulting (diagonal) density matrix an initial Fock matrix is formed and diagonalized to yield
 * the guess orbitals.
 */
class AtomicDensityGuessCalculator : public DensityInitialGuessCalculator {
 public:
  /**
   * @brief Default constructor.
   * @param scf Switch for HF-SCF densities or density guesses.
   */
  AtomicDensityGuessCalculator(GUESSMODES scf) : _scf(scf){};
  /**
   * @brief Default destructor.
   */
  virtual ~AtomicDensityGuessCalculator() = default;

  /**
   * @brief Creates a guessed (approximate) density matrix for the given system.
   * @returns The guessed density matrix.
   */
  std::unique_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>>
  calculateInitialDensity(std::shared_ptr<SystemController> systemController);

 private:
  GUESSMODES _scf;
  std::map<std::string, Eigen::MatrixXd> _atomDensities;
  /// @brief Calculates the atomic density matrix for a given atom.
  Eigen::MatrixXd performAtomInitialGuess(Settings settings, std::shared_ptr<Atom> atom);
};

} /* namespace Serenity */
#endif /* ATOMICDENSITYGUESSCALCULATOR_H */
