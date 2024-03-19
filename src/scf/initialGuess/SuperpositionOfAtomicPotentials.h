/**
 * @file SuperpositionOfAtomicPotentials.h
 *
 * @date   Aug 9, 2019
 * @author Jan Unsleber
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

#ifndef SUPERPOSITIONOFATOMICPOTENTIALS_H
#define SUPERPOSITIONOFATOMICPOTENTIALS_H
/* Include Serenity Internal Headers */
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory> //smart ptr.

namespace Serenity {

/**
 * @class SuperpositionOfAtomicPotentials
 * @brief Initial guess orbitals created using a grid based evaluation of tabulated atomic potentials.
 *
 * This guess has been described in:
 * "Assessment of Initial Guesses for Self-Consistent Field Calculations.
 *  Superposition of Atomic Potentials: Simple yet Efficient"
 *  - Susi Lehtola
 * J. Chem. Theory Comput.,2019, 15, 3
 * https://pubs.acs.org/doi/10.1021/acs.jctc.8b01089
 *
 * The data used is the non-relativistic CAPX SAP guess, the data used is taken from supplementary
 * material provided in the publication.
 */
class SuperpositionOfAtomicPotentials : public InitialGuessCalculator<RESTRICTED> {
 public:
  /// @brief Default constructor.
  SuperpositionOfAtomicPotentials() = default;
  /// @brief Default destructor.
  virtual ~SuperpositionOfAtomicPotentials() = default;
  std::unique_ptr<ElectronicStructure<RESTRICTED>>
  calculateInitialGuess(const std::shared_ptr<SystemController> systemController) override final;
};

} /* namespace Serenity */
#endif /* SUPERPOSITIONOFATOMICPOTENTIALS_H */
