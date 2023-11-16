/**
 * @file   SerenityError_test.cpp
 * @author Jan Unsleber
 *
 * @date   Oct 25, 2017
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
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @test
 * @brief Tests if the correct message is given as what().
 *
 * This test should also print the Serenity footer.
 */
TEST(SerenityErrorTest, Message) {
  SerenityError error("asdfghjkl");
  std::string ret(error.what());
  EXPECT_FALSE(ret.compare("asdfghjkl"));
  EXPECT_TRUE(ret.compare("asdfghjklm"));
}

} /* namespace Serenity */
