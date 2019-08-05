/**
 * @file XCFun_test.cpp
 * @author: Jan Unsleber
 *
 * @date Sep 20, 2017
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "dft/functionals/wrappers/XCFun.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <xcfun.h>

namespace Serenity {

class XCFunTest : public ::testing::Test {
protected:
  XCFunTest(){
  }
  virtual ~XCFunTest() = default;
};

/**
 * @test XCFunTest
 * @brief Tests the functional aliases of the implemented functionals.
 */
TEST_F(XCFunTest, FunctionalAliases) {
  XCFun<RESTRICTED> xcfun(1);
  for (unsigned int i=0;i<static_cast<unsigned int>(BASIC_FUNCTIONALS::SAOP);i++){
    xc_functional func = xc_new_functional();
    int iErr = xc_set(func,xcfun.getAlias(static_cast<BASIC_FUNCTIONALS>(i)),1.0);
    EXPECT_TRUE(iErr==0);
  }
}

} /* namespace Serenity */
