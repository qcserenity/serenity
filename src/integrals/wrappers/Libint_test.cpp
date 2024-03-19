/**
 * @file   Libint_test.cpp
 *
 * @date   Apr 2, 2014
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
/* Include Serenity Internal Headers */
#include "integrals/wrappers/Libint.h"
#include "basis/Basis.h"
#include "basis/Shell.h"
#include "geometry/Atom.h"
#include "geometry/Point.h"
#include "integrals/Normalization.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * The test environment to test Libint
 *
 * The integral tests were modified to start from normalized basis functions (before only the
 * primitives were normalized, but not the contraction coefficients). The results have not been
 * recalculated fresh with Mathematica, but only corrected for the renormalization.
 */
class LibintTest : public testing::Test {
 protected:
  LibintTest() = default;
  virtual ~LibintTest() = default;
  /**
   * The object under test
   */
  Libint& libint = Libint::getInstance();
};
/**
 * @test
 * @brief Overlap integrals
 *
 * Limited accuracy (due to limited accuracy of the reference numbers created with numerical
 * integration using Mathematica).
 */
TEST_F(LibintTest, OverlapIntegrals) {
  const double expectedPrecision = 1e-7;
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  auto basis = basisController->getBasis();
  Eigen::MatrixXd results;

  /* Overlaps */
  libint.initialize(LIBINT_OPERATOR::overlap, 0, 2);
  libint.compute(LIBINT_OPERATOR::overlap, 0, *basis[0], *basis[0], results);
  ASSERT_EQ(results.size(), (unsigned int)1);
  EXPECT_NEAR(results(0, 0), 1.1930904 / 1.1930904, expectedPrecision);
  libint.compute(LIBINT_OPERATOR::overlap, 0, *basis[0], *basis[1], results);
  ASSERT_EQ(results.size(), (unsigned int)3);
  EXPECT_NEAR(results(0, 0), -0.71871395 / sqrt(1.1930904 * 1.1826193), expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.0, expectedPrecision);
  libint.compute(LIBINT_OPERATOR::overlap, 0, *basis[0], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)6);
  EXPECT_NEAR(results(0, 0), 0.049524003 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(1, 0), -0.12211283 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(2, 0), -0.048845134 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.19757795 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(4, 0), 0.12211283 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(5, 0), 0.049524003 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  libint.compute(LIBINT_OPERATOR::overlap, 0, *basis[1], *basis[1], results);
  ASSERT_EQ(results.size(), (unsigned int)9);
  EXPECT_NEAR(results(0, 0), 1.1826193 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(4, 0), 1.1826193 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(5, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(6, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(7, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(8, 0), 1.1826193 / 1.1826193, expectedPrecision);
  libint.compute(LIBINT_OPERATOR::overlap, 0, *basis[1], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)18);
  EXPECT_NEAR(results(0, 0), -0.025154920 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.074263027 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.029705211 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(3, 0), -0.094007735 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(4, 0), -0.058067663 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(5, 0), -0.023604436 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(6, 0), 0.079793615 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(7, 0), -0.11167132 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(8, 0), -0.058067663 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(9, 0), 0.069159704 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(10, 0), 0.055835659 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(11, 0), 0.029505544 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(12, 0), 0.031917446 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(13, 0), -0.058067663 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(14, 0), 0.010270773 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(15, 0), 0.047003867 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(16, 0), -0.012838467 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(17, 0), -0.0075377683 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  libint.compute(LIBINT_OPERATOR::overlap, 0, *basis[2], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)36);
  EXPECT_NEAR(results(0, 0), 1.2623580 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.42078602 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(4, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(5, 0), 0.42078602 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(6, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(7, 0), 1.2623580 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(8, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(9, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(10, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(11, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(12, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(13, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(14, 0), 1.2623580 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(15, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(16, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(17, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(18, 0), 0.42078602 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(19, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(20, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(21, 0), 1.2623580 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(22, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(23, 0), 0.42078602 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(24, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(25, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(26, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(27, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(28, 0), 1.2623580 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(29, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(30, 0), 0.42078602 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(31, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(32, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(33, 0), 0.42078602 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(34, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(35, 0), 1.2623580 / 1.2623580, expectedPrecision);
  libint.finalize(LIBINT_OPERATOR::overlap, 0, 2);
}
/**
 * @test
 * @brief Kinetic energy integrals
 *
 * Limited accuracy (due to limited accuracy of the reference numbers created with numerical
 * integration using Mathematica).
 */
TEST_F(LibintTest, KineticEnergyIntegrals) {
  const double expectedPrecision = 2e-7;
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  auto basis = basisController->getBasis();
  Eigen::MatrixXd results;
  libint.initialize(LIBINT_OPERATOR::kinetic, 0, 2);
  libint.compute(LIBINT_OPERATOR::kinetic, 0, *basis[0], *basis[0], results);
  ASSERT_EQ(results.size(), (unsigned int)1);
  EXPECT_NEAR(results(0, 0), 1.8961808 / 1.1930904, expectedPrecision);
  libint.compute(LIBINT_OPERATOR::kinetic, 0, *basis[0], *basis[1], results);
  ASSERT_EQ(results.size(), (unsigned int)3);
  EXPECT_NEAR(results(0, 0), -1.4901775 / sqrt(1.1930904 * 1.1826193), expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.0, expectedPrecision);
  libint.compute(LIBINT_OPERATOR::kinetic, 0, *basis[0], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)6);
  EXPECT_NEAR(results(0, 0), -0.044097900 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(1, 0), -0.039677301 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(2, 0), -0.015870920 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.0040082713 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(4, 0), 0.039677301 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(5, 0), -0.044097900 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  libint.compute(LIBINT_OPERATOR::kinetic, 0, *basis[1], *basis[1], results);
  ASSERT_EQ(results.size(), (unsigned int)9);
  EXPECT_NEAR(results(0, 0), 3.1253977 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(4, 0), 3.1253977 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(5, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(6, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(7, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(8, 0), 3.1253977 / 1.1826193, expectedPrecision);
  libint.compute(LIBINT_OPERATOR::kinetic, 0, *basis[1], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)18);
  EXPECT_NEAR(results(0, 0), -0.0024641700 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.047322082 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.018928833 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.0016494222 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(4, 0), -0.016729361 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(5, 0), 0.021932695 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(6, 0), -0.012927817 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(7, 0), -0.052914091 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(8, 0), -0.016729361 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(9, 0), 0.013946251 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(10, 0), 0.026457045 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(11, 0), -0.027415868 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(12, 0), -0.0051711267 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(13, 0), -0.016729361 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(14, 0), -0.017782432 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(15, 0), -0.00082471112 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(16, 0), 0.022228040 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(17, 0), -0.0045631357 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  libint.compute(LIBINT_OPERATOR::kinetic, 0, *basis[2], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)36);
  EXPECT_NEAR(results(0, 0), 1.7276239 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(3, 0), -0.22857367 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(4, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(5, 0), -0.22857367 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(6, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(7, 0), 2.9342964 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(8, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(9, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(10, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(11, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(12, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(13, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(14, 0), 2.9342964 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(15, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(16, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(17, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(18, 0), -0.22857367 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(19, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(20, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(21, 0), 1.7276239 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(22, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(23, 0), -0.22857367 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(24, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(25, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(26, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(27, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(28, 0), 2.9342964 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(29, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(30, 0), -0.22857367 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(31, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(32, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(33, 0), -0.22857367 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(34, 0), 0.0, expectedPrecision);
  EXPECT_NEAR(results(35, 0), 1.7276239 / 1.2623580, expectedPrecision);
  libint.finalize(LIBINT_OPERATOR::kinetic, 0, 2);
}
/**
 * @test
 * @brief Nuclear attraction integrals
 *
 * Limited accuracy (due to limited accuracy of the reference numbers created with numerical
 * integration using Mathematica).
 */
TEST_F(LibintTest, NuclearAttractionIntegrals) {
  const double expectedPrecision = 2e-7;
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  auto basis = basisController->getBasis();
  Eigen::MatrixXd results;

  std::vector<std::pair<double, std::array<double, 3>>> ptchrgs = {{11.0, {{1.7, -0.4, 3.3}}}, {3.0, {{-1.2, 0.1, -2.0}}}};

  libint.initialize(LIBINT_OPERATOR::nuclear, 0, 2, ptchrgs, 0.0);
  libint.compute(LIBINT_OPERATOR::nuclear, 0, *basis[0], *basis[0], results);
  ASSERT_EQ(results.size(), (unsigned int)1);
  EXPECT_NEAR(results(0, 0), -5.048262928 / 1.1930904, expectedPrecision);
  libint.compute(LIBINT_OPERATOR::nuclear, 0, *basis[0], *basis[1], results);
  ASSERT_EQ(results.size(), (unsigned int)3);
  EXPECT_NEAR(results(0, 0), 3.051721730 / sqrt(1.1930904 * 1.1826193), expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.028497005 / sqrt(1.1930904 * 1.1826193), expectedPrecision);
  EXPECT_NEAR(results(2, 0), -0.168618146 / sqrt(1.1930904 * 1.1826193), expectedPrecision);
  libint.compute(LIBINT_OPERATOR::nuclear, 0, *basis[0], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)6);
  EXPECT_NEAR(results(0, 0), -0.205663564 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.517657942 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.194045774 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(3, 0), -0.820713984 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(4, 0), -0.489827936 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(5, 0), -0.196849405 / sqrt(1.1930904 * 1.2623580), expectedPrecision);
  libint.compute(LIBINT_OPERATOR::nuclear, 0, *basis[1], *basis[1], results);
  ASSERT_EQ(results.size(), (unsigned int)9);
  EXPECT_NEAR(results(0, 0), -4.973913164 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.008098898 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(2, 0), -0.0947756602 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.008098898 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(4, 0), -4.914930511 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(5, 0), 0.029125858 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(6, 0), -0.0947756602 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(7, 0), 0.029125858 / 1.1826193, expectedPrecision);
  EXPECT_NEAR(results(8, 0), -5.177260422 / 1.1826193, expectedPrecision);
  libint.compute(LIBINT_OPERATOR::nuclear, 0, *basis[1], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)18);
  EXPECT_NEAR(results(0, 0), 0.104093519 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(1, 0), -0.3133169381 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(2, 0), -0.1182615205 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.396791179 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(4, 0), 0.2343042872 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(5, 0), 0.094260381 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(6, 0), -0.324561973 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(7, 0), 0.461937957 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(8, 0), 0.221420690 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(9, 0), -0.2863814663 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(10, 0), -0.2172143901 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(11, 0), -0.1129544822 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(12, 0), -0.154488863 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(13, 0), 0.2836515437 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(14, 0), -0.041894137 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(15, 0), -0.2260676973 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(16, 0), 0.0543626809 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  EXPECT_NEAR(results(17, 0), 0.026792177 / sqrt(1.1826193 * 1.2623580), expectedPrecision);
  libint.compute(LIBINT_OPERATOR::nuclear, 0, *basis[2], *basis[2], results);
  ASSERT_EQ(results.size(), (unsigned int)36);
  EXPECT_NEAR(results(0, 0), -3.967376346 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.102691822 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(2, 0), -0.089062072 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(3, 0), -1.337617893 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(4, 0), 0.016568569 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(5, 0), -1.332951704 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(6, 0), 0.102691822 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(7, 0), -4.012853679 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(8, 0), 0.028697603 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(9, 0), 0.103493002 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(10, 0), -0.062618008 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(11, 0), 0.036826278 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(12, 0), -0.089062072 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(13, 0), 0.028697603 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(14, 0), -3.998855112 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(15, 0), -0.036152524 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(16, 0), 0.063784985 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(17, 0), -0.0882305794 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(18, 0), -1.337617893 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(19, 0), 0.103493002 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(20, 0), -0.036152524 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(21, 0), -4.0402149262 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(22, 0), 0.027498029 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(23, 0), -1.3507669655 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(24, 0), 0.016568569 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(25, 0), -0.062618008 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(26, 0), 0.063784985 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(27, 0), 0.027498029 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(28, 0), -4.052300897 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(29, 0), 0.0195141542 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(30, 0), -1.332951704 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(31, 0), 0.036826278 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(32, 0), -0.0882305794 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(33, 0), -1.3507669655 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(34, 0), 0.0195141542 / 1.2623580, expectedPrecision);
  EXPECT_NEAR(results(35, 0), -4.026150420 / 1.2623580, expectedPrecision);
  libint.finalize(LIBINT_OPERATOR::nuclear, 0, 2);
}
/**
 * @test
 * @brief Integrals with f-Functions
 *
 * Limited accuracy (due to limited accuracy of the reference numbers created with numerical
 * integration using Mathematica).
 * Although after d-functions nothing fundamentally new happens there was a bug with
 * f-type functions. This test shall make sure it won't happen again.
 */
TEST_F(LibintTest, F_typeFunctionIntegrals) {
  const double expectedPrecision = 1e-8;
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  auto basis = basisController->getBasis();
  Shell fFunction({0.25}, {1.0}, 3, false, {{3.1, 1.2, -2.0}});
  Eigen::MatrixXd results;
  std::vector<std::pair<double, std::array<double, 3>>> ptchrgs = {{11.0, {{1.7, -0.4, 3.3}}}, {3.0, {{-1.2, 0.1, -2.0}}}};

  /* Overlaps */
  libint.initialize(LIBINT_OPERATOR::overlap, 0, 2);
  libint.compute(LIBINT_OPERATOR::overlap, 0, *basis[0], fFunction, results);
  ASSERT_EQ(results.size(), (unsigned int)10);
  EXPECT_NEAR(results(0, 0), -0.176338013 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(1, 0), -0.136513533 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.227522554 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(3, 0), -0.070545824 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(4, 0), 0.143540568 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(5, 0), -0.158943888 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(6, 0), -0.019421705 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(7, 0), 0.045513435 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(8, 0), -0.061526666 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(9, 0), 0.0578745502 / sqrt(1.1930904), expectedPrecision);
  libint.finalize(LIBINT_OPERATOR::overlap, 0, 2);
  /* Kinetic */
  libint.initialize(LIBINT_OPERATOR::kinetic, 0, 2);
  libint.compute(LIBINT_OPERATOR::kinetic, 0, *basis[0], fFunction, results);
  ASSERT_EQ(results.size(), (unsigned int)10);
  EXPECT_NEAR(results(0), -0.035495889 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(1, 0), -0.060828579 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(2, 0), 0.101380965 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.009511464 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(4, 0), 0.084793046 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(5, 0), -0.042707509 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(6, 0), 0.015109600 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(7, 0), -0.006136429 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(8, 0), -0.016531939 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(9, 0), -0.010116192 / sqrt(1.1930904), expectedPrecision);
  libint.finalize(LIBINT_OPERATOR::kinetic, 0, 2);
  /* Nuclear */
  libint.initialize(LIBINT_OPERATOR::nuclear, 0, 2, ptchrgs, 0.0);
  libint.compute(LIBINT_OPERATOR::nuclear, 0, *basis[0], fFunction, results);
  ASSERT_EQ(results.size(), (unsigned int)10);
  EXPECT_NEAR(results(0, 0), 0.7287011069 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(1, 0), 0.569204468 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(2, 0), -0.946276711 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(3, 0), 0.290752600 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(4, 0), -0.605218277 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(5, 0), 0.670843881 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(6, 0), 0.079272777 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(7, 0), -0.190720343 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(8, 0), 0.264209166 / sqrt(1.1930904), expectedPrecision);
  EXPECT_NEAR(results(9, 0), -0.2491949391 / sqrt(1.1930904), expectedPrecision);
  libint.finalize(LIBINT_OPERATOR::nuclear, 0, 2);
}
/**
 * @test
 * @brief Integration over two-center integrals (P|1/r|Q)
 *
 * Same and different centers, exponents and angular momenta are tested. Limited accuracy (due to
 * limited accuracy of the reference numbers created with numerical integration using Mathematica).
 * Caveat: Assumes a working normalization (not unit tested seperately).
 */
TEST_F(LibintTest, Test2centerIntegrals) {
  Shell dummyBasFunc({1.0}, {1.0}, 0, false, {{0.0, 0.0, 0.0}});
  Shell dummyBasFuncShifted({1.0}, {1.0}, 0, false, {{1.0, 1.0, 1.0}});
  Shell dummyBasFuncP({0.5}, {1.0}, 1, false, {{0.0, 0.0, 0.0}});
  Shell dummyBasFuncD({0.5}, {1.0}, 2, false, {{0.0, 0.0, 0.0}});
  Shell dummyBasFuncDAgain({0.5}, {1.0}, 2, false, {{3.0, -1.0, 2.0}});
  Shell dummyBasFuncShifted2({1.0}, {1.0}, 0, false, {{1.0, 2.0, 3.0}});
  Shell dummyBasFuncShifted2Again({1.0}, {1.0}, 0, false, {{4.0, 1.0, 5.0}});
  Shell dummyBasFuncShiftedP({1.0}, {1.0}, 1, false, {{1.0, 2.0, 3.0}});
  double normFactor = Normalization::normalizeTotalAngMom(0, 1.0);
  normFactor *= Normalization::finalNormalization(0, 0, 0);
  double normFactorP = Normalization::normalizeTotalAngMom(1, 0.5);
  normFactorP *= Normalization::finalNormalization(0, 1, 0);
  double normFactorD = Normalization::normalizeTotalAngMom(2, 0.5);
  normFactorD *= Normalization::finalNormalization(0, 1, 1);
  double normFactorP2 = Normalization::normalizeTotalAngMom(1, 1.0);
  normFactorP2 *= Normalization::finalNormalization(0, 1, 0);

  Eigen::MatrixXd results;
  libint.initialize(LIBINT_OPERATOR::coulomb, 0, 2);
  // s-function with itself, exp 1.0
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFunc, dummyBasFunc, results);
  EXPECT_NEAR(results(0, 0), 24.739 * normFactor * normFactor, 1e-3);
  // s-function with s-function, exp 1.0, different centers
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFunc, dummyBasFuncShifted, results);
  EXPECT_NEAR(results(0, 0), 16.412 * normFactor * normFactor, 1e-3);
  // p-function, exp 0.5 with s-function, exp 1.0, different centers
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncP, dummyBasFuncShifted, results);
  EXPECT_NEAR(results(1, 0), 7.2170 * normFactorP * normFactor, 1e-4);
  // d-function, exp 0.5 with s-function, exp 1.0, different centers
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncD, dummyBasFuncShifted2, results);
  EXPECT_NEAR(results(4, 0), 1.94482 * normFactorD * normFactor, 1e-4);
  // same as above, given in base functions the other way around
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncShifted2, dummyBasFuncD, results);
  EXPECT_NEAR(results(4, 0), 1.94482 * normFactorD * normFactor, 1e-4);
  // Same as above, now no function at origin
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncDAgain, dummyBasFuncShifted2Again, results);
  EXPECT_NEAR(results(4, 0), 1.94482 * normFactorD * normFactor, 1e-4);
  // d-function, exp 0.5 with p-function, exp 1.0, different centers
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncD, dummyBasFuncShiftedP, results);
  EXPECT_NEAR(results(13, 0), -0.10605 * normFactorD * normFactorP2, 1e-4);
  // same as above, permuted
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncShiftedP, dummyBasFuncD, results);
  EXPECT_NEAR(results(10, 0), -0.10605 * normFactorD * normFactorP2, 1e-4);
  libint.finalize(LIBINT_OPERATOR::coulomb, 0, 2);
}

/**
 * @test
 * @brief Integration over three-center integrals (P|1/r|kl)
 *
 * Same and different centers, exponents and angular momenta are tested. Limited accuracy (due to
 * limited accuracy of the reference numbers created with numerical integration using Mathematica).
 * Caveat: Assumes a working normalization (not unit tested seperately).
 */
TEST_F(LibintTest, Test3centerIntegrals) {
  Shell dummyBasFunc({1.0}, {1.0}, 0, false, {{0.0, 0.0, 0.0}});
  Shell dummyBasFunc2({0.5}, {1.0}, 0, false, {{0.0, 0.0, 0.0}});
  Shell dummyBasFuncShifted({1.0}, {1.0}, 0, false, {{1.0, 2.0, 3.0}});
  Shell dummyBasFuncP({1.0}, {1.0}, 1, false, {{1.0, 2.0, 3.0}});
  Shell dummyBasFuncPShifted({1.0}, {1.0}, 1, false, {{3.0, -1.0, 2.0}});
  Shell dummyBasFuncP2({0.5}, {1.0}, 1, false, {{0.5, 4.5, 0.5}});
  Shell dummyBasFuncP2Shifted({0.5}, {1.0}, 1, false, {{2.5, 1.5, -0.5}});
  Shell dummyBasFuncD2({0.5}, {1.0}, 2, false, {{0.0, 0.0, 0.0}});
  Shell dummyBasFuncD2Shifted({0.5}, {1.0}, 2, false, {{2.0, -3.0, -1.0}});
  double normFactor = Normalization::normalizeTotalAngMom(0, 1.0);
  normFactor *= Normalization::finalNormalization(0, 0, 0);
  double normFactor2 = Normalization::normalizeTotalAngMom(0, 0.5);
  normFactor2 *= Normalization::finalNormalization(0, 0, 0);
  double normFactorP = Normalization::normalizeTotalAngMom(1, 1.0);
  normFactorP *= Normalization::finalNormalization(0, 1, 0);
  double normFactorP2 = Normalization::normalizeTotalAngMom(1, 0.5);
  normFactorP2 *= Normalization::finalNormalization(0, 1, 0);
  double normFactorD2 = Normalization::normalizeTotalAngMom(2, 0.5);
  normFactorD2 *= Normalization::finalNormalization(0, 1, 1);

  Eigen::MatrixXd results;
  // All s-functions, same center
  libint.initialize(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFunc, dummyBasFunc2, dummyBasFunc2, results);
  EXPECT_NEAR(results(0, 0), 24.739 * normFactor * normFactor2 * normFactor2, 1e-4);
  // all s-functions, different centers
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFunc, dummyBasFuncShifted, dummyBasFuncShifted, results);
  EXPECT_NEAR(results(0, 0), 2.92827 * normFactor * normFactor * normFactor, 1e-3);
  // d and p-functions, all different centers (checking D_yz, P_y, P_x)
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncD2, dummyBasFuncP, dummyBasFuncP2, results);
  EXPECT_NEAR(results(39, 0), 0.00287007 * normFactorD2 * normFactorP * normFactorP2, 1e-5);
  // same as above, shifted
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncD2Shifted, dummyBasFuncPShifted, dummyBasFuncP2Shifted, results);
  EXPECT_NEAR(results(39, 0), 0.00287007 * normFactorD2 * normFactorP * normFactorP2, 1e-5);
  // Same as above permuted
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncP, dummyBasFuncD2, dummyBasFuncP2, results);
  EXPECT_NEAR(results(30, 0), -0.000223906 * normFactorP * normFactorD2 * normFactorP2, 1e-5);
  // Same as above permuted again
  libint.compute(LIBINT_OPERATOR::coulomb, 0, dummyBasFuncP, dummyBasFuncP2, dummyBasFuncD2, results);
  EXPECT_NEAR(results(22, 0), -0.000223906 * normFactorP * normFactorP2 * normFactorD2, 1e-5);
  libint.finalize(LIBINT_OPERATOR::coulomb, 0, 3);
}

/**
 * @test
 * @brief Test the 1eInt function for cross basis terms works.
 */
TEST_F(LibintTest, 1eInts_CrossBasis_Overlap) {
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  /* Overlaps */
  auto res1 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController, basisController);
  auto res2 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController);
  for (unsigned int i = 0; i < res1.cols(); i++) {
    for (unsigned int j = 0; j < res1.cols(); j++) {
      EXPECT_DOUBLE_EQ(res1(i, j), res2(i, j));
    }
  }
}

/**
 * @test
 * @brief Test the 1eInt function for cross basis terms works.
 */
TEST_F(LibintTest, 1eInts_CrossBasis_Nuclear1) {
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  /* Overlaps */
  auto res1 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController, basisController);
  auto res2 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController);
  for (unsigned int i = 0; i < res1.cols(); i++) {
    for (unsigned int j = 0; j < res1.cols(); j++) {
      EXPECT_DOUBLE_EQ(res1(i, j), res2(i, j));
    }
  }
}

/**
 * @test
 * @brief Test the 1eInt function for cross basis terms works.
 */
TEST_F(LibintTest, 1eInts_CrossBasis_Nuclear2) {
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  std::vector<std::pair<double, std::array<double, 3>>> test = {{1.0, std::array<double, 3>{0.0, 0.0, 0.0}}};
  /* Overlaps */
  auto res1 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController, basisController, {test});
  auto res2 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController, test);
  for (unsigned int i = 0; i < res1.cols(); i++) {
    for (unsigned int j = 0; j < res1.cols(); j++) {
      EXPECT_DOUBLE_EQ(res1(i, j), res2(i, j));
    }
  }
}

/**
 * @test
 * @brief Test the nuclear attraction for both ways of lading nucleii.
 */
TEST_F(LibintTest, 1eInts_NuclearVersions1) {
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  auto test1 = std::make_shared<Atom>("H", 0.0, 0.0, 0.0);
  std::vector<std::pair<double, std::array<double, 3>>> test2 = {{1.0, std::array<double, 3>{0.0, 0.0, 0.0}}};
  /* Overlaps */
  auto res1 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController, {test1});
  auto res2 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController, test2);
  for (unsigned int i = 0; i < res1.cols(); i++) {
    for (unsigned int j = 0; j < res1.cols(); j++) {
      EXPECT_DOUBLE_EQ(res1(i, j), res2(i, j));
    }
  }
}

/**
 * @test
 * @brief Test the nuclear attraction for both ways of lading nucleii.
 */
TEST_F(LibintTest, 1eInts_NuclearVersions2) {
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  auto test1 = std::make_shared<Atom>("H", 0.0, 0.0, 0.0);
  std::pair<double, std::array<double, 3>> test2 = {1.0, {{0.0, 0.0, 0.0}}};
  /* Overlaps */
  auto res1 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController, basisController, {test1});
  auto res2 = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisController, basisController, {test2});
  for (unsigned int i = 0; i < res1.cols(); i++) {
    for (unsigned int j = 0; j < res1.cols(); j++) {
      EXPECT_DOUBLE_EQ(res1(i, j), res2(i, j));
    }
  }
}

/**
 * @test
 * @brief Test if the clearing of engines crashes.
 */
TEST_F(LibintTest, ClearAllEngines) {
  libint.clearAllEngines();
}

} // namespace Serenity
