/**
 * @file AtomTypeFactory_test.cpp
 *
 * @date Mar 17, 2017
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
#include "geometry/AtomTypeFactory.h"
#include "geometry/AtomType.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class AtomTypeFactoryTest AtomTypeFactory_test.cpp
 * @brief Sets everything up for the tests of AtomTypeFactory.h/.cpp .
 */
class AtomTypeFactoryTest : public ::testing::Test {};

/**
 * @test
 * @brief Tests AtomTypeFactory.h/.cpp: Test cases.
 */
TEST_F(AtomTypeFactoryTest, Cases) {
  auto& fac = AtomTypeFactory::getInstance();
  auto he = fac.getAtomType("he");
  auto He = fac.getAtomType("He");
  auto HE = fac.getAtomType("HE");
  auto hE = fac.getAtomType("hE");
  EXPECT_TRUE(he == He);
  EXPECT_TRUE(he == HE);
  EXPECT_TRUE(he == hE);
}
/**
 * @test
 * @brief Tests AtomTypeFactory.h/.cpp: Test the PSE.
 */
TEST_F(AtomTypeFactoryTest, PSE_charges) {
  auto& fac = AtomTypeFactory::getInstance();
  std::vector<std::string> pse = {
      "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
      "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
      "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
      "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
      "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
      "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
      "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

  for (unsigned int i = 0; i < pse.size(); ++i) {
    auto atomt = fac.getAtomType(pse[i]);
    EXPECT_EQ((int)i + 1, atomt->getNuclearCharge());
  }
}
/**
 * @test
 * @brief Tests AtomTypeFactory.h/.cpp: Test the isotopes of hydrogen.
 */
TEST_F(AtomTypeFactoryTest, DT_charges) {
  auto& fac = AtomTypeFactory::getInstance();
  auto d = fac.getAtomType("D");
  EXPECT_EQ(1, d->getNuclearCharge());
  EXPECT_EQ(2.01410177812, d->getMass());
  auto t = fac.getAtomType("T");
  EXPECT_EQ(1, t->getNuclearCharge());
  EXPECT_EQ(3.0160492779, t->getMass());
}
/**
 * @test
 * @brief Tests AtomTypeFactory.h/.cpp: Test ghost atoms.
 */
TEST_F(AtomTypeFactoryTest, Ghost) {
  auto& fac = AtomTypeFactory::getInstance();
  auto he = fac.getAtomType("he:");
  EXPECT_EQ(0, he->getMass());
  EXPECT_EQ(0, he->getNuclearCharge());
  EXPECT_EQ(2, he->getPSEPosition());
}
} /*namespace Serenity*/
