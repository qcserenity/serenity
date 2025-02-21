/**
 * @file   BasicFunctionals_test.cpp
 *
 * @date   Sep 3, 2020
 * @author Jan P. Unsleber
 *
 * IMPORTANT:\n
 * This file was automatically generated, please do not alter it.
 * Any required changes should be made to the generating Python script
 * which should be located close by.
 *
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
#include "dft/functionals/BasicFunctionals.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
namespace BasicFunctionals {

#ifdef SERENITY_USE_LIBXC
TEST(BasicFunctionalsTest, BasicFunctionalsAliases_LibXC) {
  const unsigned int nBFuncs = 556;
  for (unsigned int i = 0; i < nBFuncs; i++) {
    auto func = static_cast<BASIC_FUNCTIONALS>(i);
    if (func == BASIC_FUNCTIONALS::NONE || func == BASIC_FUNCTIONALS::XC_SAOP)
      continue;
    if (getImplementation[i] == IMPLEMENTATIONS::LIBXC || getImplementation[i] == IMPLEMENTATIONS::BOTH) {
      ASSERT_NO_THROW(getLibXCAlias(func));
    }
    else {
      ASSERT_THROW(getLibXCAlias(func), SerenityError);
    }
  }
}
#endif /* SERENITY_USE_LIBXC */

#ifdef SERENITY_USE_XCFUN
TEST(BasicFunctionalsTest, BasicFunctionalsAliases_XCFun) {
  const unsigned int nBFuncs = 556;
  for (unsigned int i = 0; i < nBFuncs; i++) {
    auto func = static_cast<BASIC_FUNCTIONALS>(i);
    if (func == BASIC_FUNCTIONALS::NONE || func == BASIC_FUNCTIONALS::XC_SAOP)
      continue;
    if (getImplementation[i] == IMPLEMENTATIONS::XCFUN || getImplementation[i] == IMPLEMENTATIONS::BOTH) {
      ASSERT_NO_THROW(getXCFunAlias(func));
    }
    else {
      ASSERT_THROW(getXCFunAlias(func), SerenityError);
    }
  }
}
#endif /* SERENITY_USE_XCFUN */

TEST(BasicFunctionalsTest, BasicFunctionalsPurposes) {
  const unsigned int nBFuncs = 556;
  for (unsigned int i = 0; i < nBFuncs; i++) {
    ASSERT_NO_THROW(getPurpose[i]);
    int p = static_cast<int>(getPurpose[i]);
    ASSERT_TRUE(p < 5);
    ASSERT_TRUE(p >= 0);
  }
}

TEST(BasicFunctionalsTest, BasicFunctionalsClass) {
  const unsigned int nBFuncs = 556;
  for (unsigned int i = 0; i < nBFuncs; i++) {
    ASSERT_NO_THROW(getClass[i]);
    int c = static_cast<int>(getClass[i]);
    ASSERT_TRUE(c < 5);
    ASSERT_TRUE(c >= 0);
  }
}

} /* namespace BasicFunctionals */
} /* namespace Serenity */
