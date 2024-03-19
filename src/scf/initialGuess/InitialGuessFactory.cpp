/**
 * @file InitialGuessFactory.cpp
 *
 * @date Jan 23, 2017
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

/* Include Class Header*/
#include "scf/initialGuess/InitialGuessFactory.h"
/* Include Serenity Internal Headers */
#include "scf/initialGuess/AtomicDensityGuessCalculator.h"
#include "scf/initialGuess/ExtendedHueckel.h"
#include "scf/initialGuess/HCoreGuessCalculator.h"
#include "scf/initialGuess/SuperpositionOfAtomicPotentials.h"
#include "settings/SCFOptions.h"

namespace Serenity {

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wswitch"

template<>
std::unique_ptr<InitialGuessCalculator<Options::SCF_MODES::RESTRICTED>>
InitialGuessFactory::produce<Options::SCF_MODES::RESTRICTED>(Options::INITIAL_GUESSES flavor) {
  switch (flavor) {
    case Options::INITIAL_GUESSES::H_CORE:
      return std::make_unique<HCoreGuessCalculator<RESTRICTED>>();
      break;
    case Options::INITIAL_GUESSES::EHT:
      return std::make_unique<ExtendedHueckel>();
      break;
    case Options::INITIAL_GUESSES::SAP:
      return std::make_unique<SuperpositionOfAtomicPotentials>();
      break;
    case Options::INITIAL_GUESSES::ATOM_SCF:
      return std::unique_ptr<AtomicDensityGuessCalculator>(new AtomicDensityGuessCalculator(GUESSMODES::SCF));
      break;
    case Options::INITIAL_GUESSES::ATOM_SCF_INPLACE:
      return std::unique_ptr<AtomicDensityGuessCalculator>(new AtomicDensityGuessCalculator(GUESSMODES::SCF_INPLACE));
      break;
      // NO default is given (intentional!). g++ will issue a warning/an error for missing cases.
  }
  throw SerenityError("Unknown flavor of initial guess.");
  return nullptr;
}

template<>
std::unique_ptr<InitialGuessCalculator<Options::SCF_MODES::UNRESTRICTED>>
InitialGuessFactory::produce<Options::SCF_MODES::UNRESTRICTED>(Options::INITIAL_GUESSES flavor) {
  switch (flavor) {
    case Options::INITIAL_GUESSES::H_CORE:
      return std::make_unique<HCoreGuessCalculator<UNRESTRICTED>>();
      break;
    case Options::INITIAL_GUESSES::EHT:
      return std::make_unique<UnrestrictedFromRestrictedGuess>(std::make_shared<ExtendedHueckel>());
      break;
    case Options::INITIAL_GUESSES::SAP:
      return std::make_unique<UnrestrictedFromRestrictedGuess>(std::make_shared<SuperpositionOfAtomicPotentials>());
      break;
    case Options::INITIAL_GUESSES::ATOM_SCF:
      return std::unique_ptr<UnrestrictedFromRestrictedGuess>(
          new UnrestrictedFromRestrictedGuess(std::make_shared<AtomicDensityGuessCalculator>(GUESSMODES::SCF)));
      break;
    case Options::INITIAL_GUESSES::ATOM_SCF_INPLACE:
      return std::unique_ptr<UnrestrictedFromRestrictedGuess>(
          new UnrestrictedFromRestrictedGuess(std::make_shared<AtomicDensityGuessCalculator>(GUESSMODES::SCF_INPLACE)));
      break;
    default:
      throw SerenityError("Unknown initial guess calculator specified for the unrestricted case.");
  }
}

#pragma GCC diagnostic pop

} /* namespace Serenity */
