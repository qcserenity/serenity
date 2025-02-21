/**
 * @file   InitialGuessCalculator.h
 *
 * @date   Jul 5, 2013
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
#ifndef INITIALGUESSCALCULATOR_H_
#define INITIALGUESSCALCULATOR_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory> //smart ptr.

namespace Serenity {
/* Forward declarations */
class SystemController;
template<Options::SCF_MODES SCFMode>
class OrbitalController;
template<Options::SCF_MODES SCFMode>
class ElectronicStructure;

/**
 * @class InitialGuessCalculator InitialGuessCalculator.h
 * @brief Interface for calculating an initial guess
 *
 * for an electronic structure calculation.
 *
 * The SCF procedure iteratively updates the electron density, the Fock matrix and the
 * orbitals (through the orbital coefficients). Because these steps form a cycle, some
 * starting point is needed. Classes that implement this interface provide orbitals in
 * some way which are a hopefully good starting point for the SCF procedure.
 */
template<Options::SCF_MODES SCFMode>
class InitialGuessCalculator;

template<>
class InitialGuessCalculator<Options::SCF_MODES::RESTRICTED> {
 public:
  InitialGuessCalculator() = default;
  virtual ~InitialGuessCalculator() = default;
  /**
   * @param   system The guess will be produced for this system
   * @returns A set of starting orbitals suited for a subsequent electronic structure calculation.
   *          The basis set is already specified together with the system.
   */
  virtual std::unique_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>
  calculateInitialGuess(std::shared_ptr<SystemController> systemController) = 0;
};

template<>
class InitialGuessCalculator<Options::SCF_MODES::UNRESTRICTED> {
 public:
  InitialGuessCalculator() = default;
  virtual ~InitialGuessCalculator() = default;
  /**
   * @param   system The guess will be produced for this system
   * @returns a set of starting orbitals suited for a subsequent electronic structure calculation.
   *          The basis set is already specified together with the system.
   */
  virtual std::unique_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>>
  calculateInitialGuess(std::shared_ptr<SystemController> systemController) = 0;

 protected:
  /**
   * @brief Introduces an asymmetry between the alpha and beta orbitals.
   *
   * In unrestricted calculations one can quickly end up with the same result as in a restricted
   * calculation in case the number of alpha and beta electrons is the same. If the alpha and beta
   * starting orbitals are the same, this is (assuming no symmetry breaking due to numerical
   * inaccuracies) always the case, because the orbitals for both spins will always feel the same
   * potential. This scramble routine induces a breaking of the symmetry to try and avoid these
   * problems.
   */
  void scrambleOrbitals(OrbitalController<Options::SCF_MODES::UNRESTRICTED>& orbitals,
                        SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int> nElectrons);
};

/**
 * @class UnrestrictedFromRestrictedGuess InitialGuessCalculator.h
 * @brief Wraps a restricted guess to provide unrestricted guess orbitals using scrambling.
 */
class UnrestrictedFromRestrictedGuess : public InitialGuessCalculator<Options::SCF_MODES::UNRESTRICTED> {
 public:
  /**
   * @brief Constructor.
   * @param restrictedGuessCalculator The restricted InitialGuessCalculator.
   */
  UnrestrictedFromRestrictedGuess(std::shared_ptr<InitialGuessCalculator<Options::SCF_MODES::RESTRICTED>> restrictedGuessCalculator)
    : _restrictedGuessCalculator(restrictedGuessCalculator) {
  }
  virtual ~UnrestrictedFromRestrictedGuess() = default;

  /**
   * @brief Calculates the unrestricted initial guess from the restricted one using scrambling.
   * @param systemController The system controller.
   * @return The unrestricted initial electronic structure.
   */
  std::unique_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>>
  calculateInitialGuess(std::shared_ptr<SystemController> systemController) override final;

 private:
  std::shared_ptr<InitialGuessCalculator<Options::SCF_MODES::RESTRICTED>> _restrictedGuessCalculator;
};

} /* namespace Serenity */
#endif /* INITIALGUESSCALCULATOR_H_ */
