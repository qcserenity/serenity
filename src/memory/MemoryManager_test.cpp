/**
 * @file   MemoryManager_test.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   24. November 2015, 16:09
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
#include "memory/MemoryManager.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <omp.h>
#include <algorithm>
#include <vector>

namespace Serenity {
class MemoryManagerTest : public testing::Test {
 protected:
  std::shared_ptr<MemoryManager> mm = MemoryManager::getInstance();
  const size_t testMemSize = 128;
};

#if __APPLE__ || __MACH__
/**
 * @test
 * @brief Ensures that the available memory is smaller than the total RAM
 */
TEST_F(MemoryManagerTest, AvailableMemory) {
  EXPECT_NE(mm->getAvailableSystemMemory(), 0);
}
#elif __unix__ || __linux__ || __unix
/**
 * @test
 * @brief Ensures the MemoryManager works in parallel
 */
TEST_F(MemoryManagerTest, ParallelAccessWorking) {
#ifdef _OPENMP
  /*
   * To increase chances for problems with concurrency the number of threads is raised to a number
   * which is probably larger than the number of available cores.
   */
  const int testNThreads = 48;
  const int oldNThreads = omp_get_max_threads();
  omp_set_num_threads(testNThreads);
#endif
  unsigned int nRequests = 10;
#pragma omp parallel for
  for (unsigned int i = 0; i < nRequests; ++i) {
    /*
     * The program might crash if the tokens are wrong.
     */
    mm->getSerenityMemoryUsage();
  }
#ifdef _OPENMP
  // Clean up
  omp_set_num_threads(oldNThreads);
#endif
}
/**
 * @test
 * @brief Ensures that the available memory is smaller than the total RAM
 */
TEST_F(MemoryManagerTest, AvailableMemory) {
  EXPECT_TRUE(mm->getAvailableSystemMemory() < mm->getSystemPhysicalMemorySize());
}
/**
 * @test
 * @brief Ensures that the used memory is smaller than the total RAM
 */
TEST_F(MemoryManagerTest, SystemMemory) {
  EXPECT_TRUE(mm->getSystemMemoryUsage() < mm->getSystemPhysicalMemorySize());
}
/**
 * @test
 * @brief Ensures that the serenity used memory is smaller than the total RAM
 */
TEST_F(MemoryManagerTest, SerenityMemory) {
  EXPECT_TRUE(mm->getSerenityMemoryUsage() < mm->getSystemPhysicalMemorySize());
}
#endif
} /* namespace Serenity */
