/**
 * @file   HCoreGuessCalculator.h
 *
 * @date   Nov 7, 2013
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
#ifndef HCOREGUESSCALCULATOR_H_
#define HCOREGUESSCALCULATOR_H_
/* Include Serenity Internal Headers */
#include "scf/initialGuess/InitialGuessCalculator.h"

namespace Serenity {
/**
 * @class HCoreGuessCalculator HCoreGuessCalculator.h
 * @brief The most simple way to calculate an initial guess for an electronic structure calculation, i.e. an initial
 * guess is created which assumes that there is no electron-electron interaction, thus the Fock matrix used in the
 * eigenvalue equation reduces to the one-electron integral matrix h.
 */
template<Options::SCF_MODES SCFMode>
class HCoreGuessCalculator final : public InitialGuessCalculator<SCFMode> {
 public:
  HCoreGuessCalculator() = default;
  virtual ~HCoreGuessCalculator() = default;

  std::unique_ptr<ElectronicStructure<SCFMode>>
  calculateInitialGuess(std::shared_ptr<SystemController> systemController) override final;
};

} /* namespace Serenity */
#endif /* HCOREGUESSCALCULATOR_H_ */
