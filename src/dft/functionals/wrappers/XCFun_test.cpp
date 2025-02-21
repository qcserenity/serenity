/**
 * @file XCFun_test.cpp
 * @author: Jan Unsleber
 *
 * @date Sep 20, 2017
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

#ifdef SERENITY_USE_XCFUN
/* Include Serenity Internal Headers */
#include "dft/functionals/wrappers/XCFun.h"
/* Include Std and External Headers */
#include <XCFun/xcfun.h>
#include <gtest/gtest.h>

namespace Serenity {

class XCFunTest : public ::testing::Test {
 protected:
  XCFunTest() {
  }
  virtual ~XCFunTest() = default;
};

} /* namespace Serenity */
#endif /* SERENITY_USE_XCFUN */