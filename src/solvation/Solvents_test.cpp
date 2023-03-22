/**
 * @file Solvents_test.cpp
 *
 * @author Moritz Bensberg
 * @date Feb. 28, 2023
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

/* Include Serenity Internal Headers */
#include "solvation/Solvents.h"  //To be tested.
#include "settings/PCMOptions.h" //Solvent enum class.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class SolventsTest : public ::testing::Test {
 protected:
  SolventsTest() {
  }
  virtual ~SolventsTest() = default;
};
TEST_F(SolventsTest, all_parameters_present) {
  std::vector<Options::PCM_SOLVENTS> allSolvents = {
      Options::PCM_SOLVENTS::WATER,         Options::PCM_SOLVENTS::PROPYLENE_CARBONATE,
      Options::PCM_SOLVENTS::DMSO,          Options::PCM_SOLVENTS::NITROMETHANE,
      Options::PCM_SOLVENTS::ACETONITRILE,  Options::PCM_SOLVENTS::METHANOL,
      Options::PCM_SOLVENTS::ETHANOL,       Options::PCM_SOLVENTS::ACETONE,
      Options::PCM_SOLVENTS::DICHLORETHANE, Options::PCM_SOLVENTS::METHYLENECHLORIDE,
      Options::PCM_SOLVENTS::THF,           Options::PCM_SOLVENTS::ANILINE,
      Options::PCM_SOLVENTS::CHLOROBENZENE, Options::PCM_SOLVENTS::CHLOROFORM,
      Options::PCM_SOLVENTS::TOLUENE,       Options::PCM_SOLVENTS::DIOXANE,
      Options::PCM_SOLVENTS::BENZENE,       Options::PCM_SOLVENTS::CARBON_TETRACHLORIDE,
      Options::PCM_SOLVENTS::CYCLOHEXANE,   Options::PCM_SOLVENTS::N_HEPTANE};
  for (auto i : allSolvents) {
    ASSERT_NO_THROW(Solvents::getProbeRadius(i));
    ASSERT_NO_THROW(Solvents::getStaticPermittivity(i));
    ASSERT_NO_THROW(Solvents::getMolecularWeight(i));
  }
}
} /* namespace Serenity */
