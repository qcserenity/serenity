/**
 * @file   CompositeFunctionals_test.cpp
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
#include "dft/functionals/CompositeFunctionals.h"
/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "dft/functionals/BasicFunctionals.h"
#include "misc/SerenityError.h"
#include "settings/DFTOptions.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
namespace CompositeFunctionals {

TEST(CompositeFunctionalsTest, ResolveEnum_Functional) {
  const unsigned int nFuncs = 57;
  for (unsigned int i = 0; i < nFuncs; i++) {
    auto func = static_cast<FUNCTIONALS>(i);
    if (getImplementation[i] == IMPLEMENTATIONS::EITHER_OR) {
      ASSERT_NO_THROW(resolveFunctional(func));
    }
    else if (getImplementation[i] == IMPLEMENTATIONS::LIBXC) {
#ifdef SERENITY_USE_LIBXC
      ASSERT_NO_THROW(resolveFunctional(func));
#else  /* SERENITY_USE_LIBXC */
      ASSERT_THROW(resolveFunctional(func), SerenityError);
#endif /* SERENITY_USE_LIBXC */
    }
    else if (getImplementation[i] == IMPLEMENTATIONS::XCFUN) {
#ifdef SERENITY_USE_XCFUN
      ASSERT_NO_THROW(resolveFunctional(func));
#else  /* SERENITY_USE_XCFUN */
      ASSERT_THROW(resolveFunctional(func), SerenityError);
#endif /* SERENITY_USE_XCFUN */
    }
    else {
      ASSERT_THROW(resolveFunctional(func), SerenityError);
    }
  }
}

TEST(CompositeFunctionalsTest, ResolveEnum_XCFunctional) {
  const unsigned int nFuncs = 57;
  for (unsigned int i = 0; i < nFuncs; i++) {
    auto func = static_cast<FUNCTIONALS>(i);
    if (getPurpose[i] != PURPOSES::EXCHANGE_CORRELATION) {
      continue;
    }
    auto xcfunc = static_cast<XCFUNCTIONALS>(i);
    if (getImplementation[i] == IMPLEMENTATIONS::EITHER_OR) {
      ASSERT_EQ(resolveFunctional(func), resolveFunctional(xcfunc));
    }
    else if (getImplementation[i] == IMPLEMENTATIONS::LIBXC) {
#ifdef SERENITY_USE_LIBXC
      ASSERT_EQ(resolveFunctional(func), resolveFunctional(xcfunc));
#else  /* SERENITY_USE_LIBXC */
      ASSERT_THROW(resolveFunctional(xcfunc), SerenityError);
#endif /* SERENITY_USE_LIBXC */
    }
    else if (getImplementation[i] == IMPLEMENTATIONS::XCFUN) {
#ifdef SERENITY_USE_XCFUN
      ASSERT_EQ(resolveFunctional(func), resolveFunctional(xcfunc));
#else  /* SERENITY_USE_XCFUN */
      ASSERT_THROW(resolveFunctional(xcfunc), SerenityError);
#endif /* SERENITY_USE_XCFUN */
    }
    else {
      ASSERT_THROW(resolveFunctional(func), SerenityError);
    }
  }
}

TEST(CompositeFunctionalsTest, ResolveEnum_KinFunctional) {
  const unsigned int nFuncs = 57;
  for (unsigned int i = 0; i < nFuncs; i++) {
    auto func = static_cast<FUNCTIONALS>(i);
    if (getPurpose[i] != PURPOSES::KINETIC) {
      continue;
    }
    auto kinfunc = static_cast<KINFUNCTIONALS>(i);
    if (getImplementation[i] == IMPLEMENTATIONS::EITHER_OR) {
      ASSERT_EQ(resolveFunctional(func), resolveFunctional(kinfunc));
    }
    else if (getImplementation[i] == IMPLEMENTATIONS::LIBXC) {
#ifdef SERENITY_USE_LIBXC
      ASSERT_EQ(resolveFunctional(func), resolveFunctional(kinfunc));
#else  /* SERENITY_USE_LIBXC */
      ASSERT_THROW(resolveFunctional(kinfunc), SerenityError);
#endif /* SERENITY_USE_LIBXC */
    }
    else if (getImplementation[i] == IMPLEMENTATIONS::XCFUN) {
#ifdef SERENITY_USE_XCFUN
      ASSERT_EQ(resolveFunctional(func), resolveFunctional(kinfunc));
#else  /* SERENITY_USE_XCFUN */
      ASSERT_THROW(resolveFunctional(kinfunc), SerenityError);
#endif /* SERENITY_USE_XCFUN */
    }
    else {
      ASSERT_THROW(resolveFunctional(func), SerenityError);
    }
  }
}

#ifdef SERENITY_USE_LIBXC
TEST(CompositeFunctionalsTest, ResolveEnum_Functional_LibXC) {
  const unsigned int nFuncs = 57;
  for (unsigned int i = 0; i < nFuncs; i++) {
    auto func = static_cast<FUNCTIONALS>(i);
    if (getImplementation[i] == IMPLEMENTATIONS::LIBXC || getImplementation[i] == IMPLEMENTATIONS::EITHER_OR) {
      ASSERT_NO_THROW(resolveLibXC(func));
    }
    else {
      ASSERT_THROW(resolveLibXC(func), SerenityError);
    }
  }
}
#endif /* SERENITY_USE_LIBXC */

#ifdef SERENITY_USE_XCFUN
TEST(CompositeFunctionalsTest, ResolveEnum_Functional_XCFun) {
  const unsigned int nFuncs = 57;
  for (unsigned int i = 0; i < nFuncs; i++) {
    auto func = static_cast<FUNCTIONALS>(i);
    if (getImplementation[i] == IMPLEMENTATIONS::XCFUN || getImplementation[i] == IMPLEMENTATIONS::EITHER_OR) {
      ASSERT_NO_THROW(resolveXCFun(func));
    }
    else {
      ASSERT_THROW(resolveXCFun(func), SerenityError);
    }
  }
}
#endif /* SERENITY_USE_XCFUN */

TEST(CompositeFunctionalsTest, ResolveString_Functional) {
  const unsigned int nFuncs = 57;
  for (unsigned int i = 0; i < nFuncs; i++) {
    auto func = static_cast<FUNCTIONALS>(i);
    std::string dummy;
    FUNCTIONALS test;
    Options::resolve<FUNCTIONALS>(dummy, func);
    Options::resolve<FUNCTIONALS>(dummy, test);
    ASSERT_EQ(func, test);
  }
}

TEST(CompositeFunctionalsTest, ResolveString_XCFunctional) {
  const unsigned int nFuncs = 57;
  for (unsigned int i = 0; i < nFuncs; i++) {
    if (getPurpose[i] != PURPOSES::EXCHANGE_CORRELATION) {
      continue;
    }
    auto func = static_cast<XCFUNCTIONALS>(i);
    std::string dummy;
    XCFUNCTIONALS test;
    Options::resolve<XCFUNCTIONALS>(dummy, func);
    Options::resolve<XCFUNCTIONALS>(dummy, test);
    ASSERT_EQ(func, test);
  }
}

TEST(CompositeFunctionalsTest, ResolveString_KinFunctional) {
  const unsigned int nFuncs = 57;
  for (unsigned int i = 0; i < nFuncs; i++) {
    if (getPurpose[i] != PURPOSES::KINETIC) {
      continue;
    }
    auto func = static_cast<KINFUNCTIONALS>(i);
    std::string dummy;
    KINFUNCTIONALS test;
    Options::resolve<KINFUNCTIONALS>(dummy, func);
    Options::resolve<KINFUNCTIONALS>(dummy, test);
    ASSERT_EQ(func, test);
  }
}

} /* namespace CompositeFunctionals */
} /* namespace Serenity */
