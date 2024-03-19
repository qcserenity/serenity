/**
 * @file InitialGuessFactory.h
 *
 * @date   Jan 23, 2017
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

#ifndef INITIALGUESSFACTORY_H
#define INITIALGUESSFACTORY_H
/* Include Serenity Internal Headers */
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory> //smart ptr.

namespace Serenity {
/* Forward Declarations */
namespace Options {
enum class INITIAL_GUESSES;
}

/**
 * @class InitialGuessFactory InitialGuessFactory.h
 * @brief Simple factory to obtain a ready-to-use initial guess
 */
class InitialGuessFactory {
 private:
  /**
   * Private default constructor. Purely static class
   */
  InitialGuessFactory() = default;

 public:
  virtual ~InitialGuessFactory() = default;
  /**
   * @param flavor Determines the kind of requested initial guess
   * @returns an initial guess of the kind determined by flavor
   */
  template<Options::SCF_MODES SCFMode>
  static std::unique_ptr<InitialGuessCalculator<SCFMode>> produce(Options::INITIAL_GUESSES flavor);
};

} /* namespace Serenity */
#endif /* INITIALGUESSFACTORY_H */
