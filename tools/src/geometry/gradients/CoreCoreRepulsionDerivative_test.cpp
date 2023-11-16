/**
 * @file   CoreCoreRepulsionDerivative_test.cpp
 *
 * @date   Nov 21, 2014
 * @author k_klah01
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
#include "geometry/gradients/CoreCoreRepulsionDerivative.h"
#include "geometry/Atom.h"
#include "geometry/AtomType.h"
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @test
 * @brief Tests CoreCoreRepulsionDerivative.h: calculation of the derivative
 *        of the core core repulsion w.r.t nuclear coordinates
 */
TEST(CoreCoreRepulsionDerivativeTest, CalculateDerivative) {
  std::shared_ptr<AtomType> atomtypeOne(new AtomType(std::string("H"), 1, 1.0, 1.0, 1.0, 1.0, 0, {}, 1.0));
  std::shared_ptr<AtomType> atomtypeTwo(new AtomType(std::string("He"), 2, 1.0, 1.0, 1.0, 1.0, 0, {}, 1.0));
  auto atomOne = std::make_shared<Atom>(atomtypeOne, 0.0, 0.0, 0.0);
  auto atomTwo = std::make_shared<Atom>(atomtypeTwo, 0.0, 0.0, 2.0);
  std::vector<std::shared_ptr<Atom>> fakeAtoms(2, nullptr);
  fakeAtoms[0] = atomOne;
  fakeAtoms[1] = atomTwo;
  auto deriv = CoreCoreRepulsionDerivative::calculateDerivative(fakeAtoms);
  EXPECT_NE(deriv(0, 2), 2.0 / 3.0);
  EXPECT_NE(deriv(0, 1), 42);
  EXPECT_EQ(deriv(1, 1), 0);
  EXPECT_EQ(deriv(0, 2), 0.5);
  EXPECT_EQ(deriv(1, 2), -0.5);
}

} // namespace Serenity
