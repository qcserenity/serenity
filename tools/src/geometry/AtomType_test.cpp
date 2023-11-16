/**
 * @file AtomType_test.cpp
 *
 * @date Oct 12, 2017
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

/* Include Serenity Internal Headers */
#include "geometry/AtomType.h"
#include "geometry/AtomTypeFactory.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>

namespace Serenity {
/**
 * @class AtomTypeTest AtomType_test.cpp
 * @brief Sets everything up for the tests of AtomType.h/.cpp .
 */
class AtomTypeTest : public ::testing::Test {};

/**
 * @test
 * @brief Tests AtomType.h/.cpp: Tests spin calculation.
 */
TEST_F(AtomTypeTest, Spin) {
  auto& fac = AtomTypeFactory::getInstance();
  auto f = fac.getAtomType("F");
  EXPECT_EQ(fabs(getAtomSpin(*f)), 1);
}
/**
 * @test
 * @brief Tests AtomType.h/.cpp: Tests occupation factors.
 */
TEST_F(AtomTypeTest, Occupations_R) {
  auto& fac = AtomTypeFactory::getInstance();
  auto f = fac.getAtomType("F");
  auto occ = getOccupationFactors<RESTRICTED>(*f);
  EXPECT_EQ(occ[0], 2);
  EXPECT_EQ(occ[1], 2);
  EXPECT_EQ(occ[2], 5.0 / 3.0);
  EXPECT_EQ(occ[3], 5.0 / 3.0);
  EXPECT_EQ(occ[4], 5.0 / 3.0);
}
/**
 * @test
 * @brief Tests AtomType.h/.cpp: Tests occupation factors.
 */
TEST_F(AtomTypeTest, Occupations_U) {
  auto& fac = AtomTypeFactory::getInstance();
  auto f = fac.getAtomType("F");
  auto occ = getOccupationFactors<UNRESTRICTED>(*f);
  EXPECT_EQ(occ.alpha[0], 1);
  EXPECT_EQ(occ.alpha[1], 1);
  EXPECT_EQ(occ.alpha[2], 1);
  EXPECT_EQ(occ.alpha[3], 1);
  EXPECT_EQ(occ.alpha[4], 1);
  EXPECT_EQ(occ.beta[0], 1);
  EXPECT_EQ(occ.beta[1], 1);
  EXPECT_EQ(occ.beta[2], 2.0 / 3.0);
  EXPECT_EQ(occ.beta[3], 2.0 / 3.0);
  EXPECT_EQ(occ.beta[4], 2.0 / 3.0);
}

/**
 * @test
 * @brief Tests AtomType.h/.cpp: Tests minimal basis size.
 */
TEST_F(AtomTypeTest, MinimalBasisSize) {
  auto& fac = AtomTypeFactory::getInstance();
  EXPECT_EQ(fac.getAtomType("F")->getMinimalBasisSize(), 5);
  EXPECT_EQ(fac.getAtomType("Ar")->getMinimalBasisSize(), 9);
  EXPECT_EQ(fac.getAtomType("Cs")->getMinimalBasisSize(), 31);
  EXPECT_EQ(fac.getAtomType("Na")->getMinimalBasisSize(), 9);
}

} /* namespace Serenity */
