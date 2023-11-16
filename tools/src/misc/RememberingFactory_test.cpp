/**
 * @file   RememberingFactory_test.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   25. November 2015, 11:07
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
#include "misc/RememberingFactory.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

/* Further Includes*/

namespace Serenity {
class RememberingFactoryTest : public ::testing::Test {
 protected:
  typedef int Dependency;
  class ProducableClass {
   public:
    ProducableClass(Dependency dep) : _dep(dep) {
    }
    Dependency _dep;
  };

  class Factory : public RememberingFactory<ProducableClass, Dependency> {
   public:
    std::shared_ptr<ProducableClass> produce(Dependency dep) {
      return getOrProduce(dep);
    }

   protected:
    std::unique_ptr<ProducableClass> produceNew(Dependency dep) {
      return std::unique_ptr<ProducableClass>(new ProducableClass(dep));
    }
  };
  Factory fac;
};

/**
 * @test
 * @brief Make sure repeated construction gives the same object if the identifiers are equal
 */
TEST_F(RememberingFactoryTest, RepeatedConstruction) {
  Dependency dep1(1);
  Dependency dep2(2);
  auto product1 = fac.produce(dep1);
  auto product2 = fac.produce(dep2);
  auto product3 = fac.produce(dep1);
  EXPECT_EQ(product1, product3);
  EXPECT_NE(product1, product2);
  EXPECT_NE(product2, product3);
}

/**
 * @test
 * @brief Make sure the automatic cleanUp works correctly.
 *
 * I.e. things are destroyed, but only if not held somewhere else.
 */
TEST_F(RememberingFactoryTest, CleanUp) {
  Dependency dep(1);
  auto product1 = fac.produce(dep);
  // This should not do anything, because product1 is still in use.
  auto product2 = fac.produce(dep);
  EXPECT_EQ(product1, product2);
  std::weak_ptr<ProducableClass> p1weak = product1;
  product1.reset();
  EXPECT_FALSE(p1weak.expired());
  product2.reset();
  // Now product1 and product2 (they are the same anyway) should be deleted.
  EXPECT_TRUE(p1weak.expired());
  /*
   * To make sure we actually get a new product now that the old one is forgot, we make an equality
   * check of the expired weak_ptr and the new shared_ptr. This is a bit tricky (see below), so I
   * use the solution suggested here:
   * http://stackoverflow.com/questions/12301916/equality-compare-stdweak-ptr
   *
   * Comment for the old version:
   * Hope that the old memory address will not be reused for product1. We try to achieve this by
   * first constructing a different product which probably uses up that memory. This procedure is,
   * however, NOT GUARANTEED TO WORK!!
   */
  product1 = fac.produce(dep);
  EXPECT_TRUE(p1weak.owner_before(product1) || product1.owner_before(p1weak));
}

} /* namespace Serenity */
