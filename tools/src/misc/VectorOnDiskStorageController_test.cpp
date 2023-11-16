/**
 * @file VectorOnDiskStorageController_test.cpp
 *
 * @date Apr 26, 2018
 * @author Moritz Bensberg
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
/* Include Serenity Internal Headers */
#include "misc/VectorOnDiskStorageController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class VectorOnDiskStorageControllerTest
 * @brief Sets everything up for the tests of VectorOnDiskStorageController.h/.cpp .
 */
class VectorOnDiskStorageControllerTest : public ::testing::Test {
 protected:
  VectorOnDiskStorageControllerTest() {
  }

  virtual ~VectorOnDiskStorageControllerTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests the storage of vector segments in memory.
 */
TEST_F(VectorOnDiskStorageControllerTest, Store_In_Memory) {
  VectorOnDiskStorageController vectorController(2e+9, "Test.h5");
  auto segment = std::make_shared<Eigen::VectorXd>(4);
  *segment << 1.5, 99.3, 48.96, 0.0;
  vectorController.storeVectorSegment(segment, "TestLabel");
  auto storedSegment = vectorController.getVectorSegment("TestLabel");
  for (unsigned int i = 0; i < segment->size(); ++i) {
    EXPECT_NEAR((*segment)[i], (*storedSegment)[i], 1e-7);
  }
}

/**
 * @test
 * @brief Tests the storage of vector segments on hard disk.
 */
TEST_F(VectorOnDiskStorageControllerTest, Store_On_Disk) {
  VectorOnDiskStorageController vectorController(0, "Test.h5");
  auto segment = std::make_shared<Eigen::VectorXd>(4);
  *segment << 1.5, 99.3, 48.96, 0.0;
  vectorController.storeVectorSegment(segment, "TestLabel");
  auto storedSegment = vectorController.getVectorSegment("TestLabel");
  for (unsigned int i = 0; i < segment->size(); ++i) {
    EXPECT_NEAR((*segment)[i], (*storedSegment)[i], 1e-7);
  }
}

/**
 * @test
 * @brief Tests the storage of vector segments in memory and on hard disk.
 */
TEST_F(VectorOnDiskStorageControllerTest, Store_On_Disk_and_in_Memory) {
  VectorOnDiskStorageController vectorController(32, "Test.h5");
  auto segmentA = std::make_shared<Eigen::VectorXd>(4);
  auto segmentB = std::make_shared<Eigen::VectorXd>(4);
  *segmentA << 1.5, 99.3, 48.96, 0.0;
  *segmentB << 1.8, 6, 3, 1.0;
  vectorController.storeVectorSegment(segmentA, "TestA");
  vectorController.storeVectorSegment(segmentB, "TestB");
  auto storedSegmentA = vectorController.getVectorSegment("TestA");
  auto storedSegmentB = vectorController.getVectorSegment("TestB");
  for (unsigned int i = 0; i < segmentA->size(); ++i) {
    EXPECT_NEAR((*segmentA)[i], (*storedSegmentA)[i], 1e-7);
    EXPECT_NEAR((*segmentB)[i], (*storedSegmentB)[i], 1e-7);
  }
}

/**
 * @test
 * @brief Tests the storage of multiple vector segments in memory.
 */
TEST_F(VectorOnDiskStorageControllerTest, Store_Multiple_in_Memory) {
  VectorOnDiskStorageController vectorController(1e+9, "Test.h5");
  auto segmentA = std::make_shared<Eigen::VectorXd>(4);
  auto segmentB = std::make_shared<Eigen::VectorXd>(4);
  auto segmentC = std::make_shared<Eigen::VectorXd>(4);
  *segmentA << 1.5, 99.3, 48.96, 0.0;
  *segmentB << 1.8, 6, 3, 1.0;
  *segmentC << 8, 3, 3, 48.5;
  vectorController.storeVectorSegment(segmentA, "TestA");
  vectorController.storeVectorSegment(segmentB, "TestB");
  vectorController.storeVectorSegment(segmentC, "TestC");
  auto storedSegmentA = vectorController.getVectorSegment("TestA");
  auto storedSegmentB = vectorController.getVectorSegment("TestB");
  auto storedSegmentC = vectorController.getVectorSegment("TestC");
  for (unsigned int i = 0; i < segmentA->size(); ++i) {
    EXPECT_NEAR((*segmentA)[i], (*storedSegmentA)[i], 1e-7);
    EXPECT_NEAR((*segmentB)[i], (*storedSegmentB)[i], 1e-7);
    EXPECT_NEAR((*segmentC)[i], (*storedSegmentC)[i], 1e-7);
  }
}

/**
 * @test
 * @brief Tests the storage of multiple vector segments on disk.
 */
TEST_F(VectorOnDiskStorageControllerTest, Store_Multiple_On_Disk) {
  VectorOnDiskStorageController vectorController(0, "Test.h5");
  auto segmentA = std::make_shared<Eigen::VectorXd>(4);
  auto segmentB = std::make_shared<Eigen::VectorXd>(4);
  auto segmentC = std::make_shared<Eigen::VectorXd>(4);
  *segmentA << 1.5, 99.3, 48.96, 0.0;
  *segmentB << 1.8, 6, 3, 1.0;
  *segmentC << 8, 3, 3, 48.5;
  vectorController.storeVectorSegment(segmentA, "TestA");
  vectorController.storeVectorSegment(segmentB, "TestB");
  vectorController.storeVectorSegment(segmentC, "TestC");
  auto storedSegmentA = vectorController.getVectorSegment("TestA");
  auto storedSegmentB = vectorController.getVectorSegment("TestB");
  auto storedSegmentC = vectorController.getVectorSegment("TestC");
  for (unsigned int i = 0; i < segmentA->size(); ++i) {
    EXPECT_NEAR((*segmentA)[i], (*storedSegmentA)[i], 1e-7);
    EXPECT_NEAR((*segmentB)[i], (*storedSegmentB)[i], 1e-7);
    EXPECT_NEAR((*segmentC)[i], (*storedSegmentC)[i], 1e-7);
  }
}

/**
 * @test
 * @brief Tests the overwriting of a segment.
 */
TEST_F(VectorOnDiskStorageControllerTest, Overwrite_Segment) {
  VectorOnDiskStorageController vectorController(0, "Test.h5");
  auto segmentA = std::make_shared<Eigen::VectorXd>(4);
  auto segmentB = std::make_shared<Eigen::VectorXd>(4);
  *segmentA << 1.5, 99.3, 48.96, 0.0;
  *segmentB << 1.8, 6, 3, 1.0;
  vectorController.storeVectorSegment(segmentA, "TestA");
  auto storedSegmentA = vectorController.getVectorSegment("TestA");
  for (unsigned int i = 0; i < segmentA->size(); ++i) {
    EXPECT_NEAR((*segmentA)[i], (*storedSegmentA)[i], 1e-7);
  }
  vectorController.storeVectorSegment(segmentB, "TestA");
  auto storedSegmentB = vectorController.getVectorSegment("TestA");
  for (unsigned int i = 0; i < segmentB->size(); ++i) {
    EXPECT_NEAR((*segmentB)[i], (*storedSegmentB)[i], 1e-7);
  }
}

/**
 * @test
 * @brief Tests the operator*.
 */
TEST_F(VectorOnDiskStorageControllerTest, ScalarProduct_Operator) {
  VectorOnDiskStorageController vectorControllerA(0, "TestA.h5");
  VectorOnDiskStorageController vectorControllerB(0, "TestB.h5");
  auto segmentA = std::make_shared<Eigen::VectorXd>(4);
  auto segmentB = std::make_shared<Eigen::VectorXd>(4);
  auto segmentC = std::make_shared<Eigen::VectorXd>(4);
  auto segmentD = std::make_shared<Eigen::VectorXd>(4);
  *segmentA << 1.5, 99.3, 48.96, 0.0;
  *segmentB << 1.8, 6, 3, 1.0;
  *segmentC << 8, 3, 3, 48.5;
  *segmentD << 8, 4, 62, 84.5;
  vectorControllerA.storeVectorSegment(segmentA, "TestA");
  vectorControllerA.storeVectorSegment(segmentB, "TestB");
  vectorControllerB.storeVectorSegment(segmentC, "TestA");
  vectorControllerB.storeVectorSegment(segmentD, "TestB");
  double scalarProduct = vectorControllerA * vectorControllerB;
  EXPECT_NEAR(scalarProduct, 765.68, 1e-7);
}

/**
 * @test
 * @brief Tests the getLabelList.
 */
TEST_F(VectorOnDiskStorageControllerTest, getLabelList) {
  VectorOnDiskStorageController vectorController(1e+9, "Test.h5");
  auto segmentA = std::make_shared<Eigen::VectorXd>(4);
  auto segmentB = std::make_shared<Eigen::VectorXd>(4);
  auto segmentC = std::make_shared<Eigen::VectorXd>(4);
  auto segmentD = std::make_shared<Eigen::VectorXd>(4);
  *segmentA << 1.5, 99.3, 48.96, 0.0;
  *segmentB << 1.8, 6, 3, 1.0;
  *segmentC << 8, 3, 3, 48.5;
  *segmentD << 8, 4, 62, 84.5;
  vectorController.storeVectorSegment(segmentA, "TestA");
  vectorController.storeVectorSegment(segmentB, "TestB");
  vectorController.storeVectorSegment(segmentC, "TestC");
  vectorController.storeVectorSegment(segmentD, "TestD");
  auto labels = vectorController.getLabelList();
  EXPECT_TRUE(labels[0] == "TestA");
  EXPECT_TRUE(labels[1] == "TestB");
  EXPECT_TRUE(labels[2] == "TestC");
  EXPECT_TRUE(labels[3] == "TestD");
}

/**
 * @test
 * @brief Tests the copy constructor.
 */
TEST_F(VectorOnDiskStorageControllerTest, Copy_Constructor) {
  VectorOnDiskStorageController vectorController(1e+9, "Test.h5");
  auto segmentA = std::make_shared<Eigen::VectorXd>(4);
  auto segmentB = std::make_shared<Eigen::VectorXd>(4);
  auto segmentC = std::make_shared<Eigen::VectorXd>(4);
  auto segmentD = std::make_shared<Eigen::VectorXd>(4);
  *segmentA << 1.5, 99.3, 48.96, 0.0;
  *segmentB << 1.8, 6, 3, 1.0;
  *segmentC << 8, 3, 3, 48.5;
  *segmentD << 8, 4, 62, 84.5;
  vectorController.storeVectorSegment(segmentA, "TestA");
  vectorController.storeVectorSegment(segmentB, "TestB");
  vectorController.storeVectorSegment(segmentC, "TestC");
  vectorController.storeVectorSegment(segmentD, "TestD");
  VectorOnDiskStorageController copy(vectorController, "CopyOfTest.h5");
  auto labels = copy.getLabelList();
  EXPECT_TRUE(labels[0] == "TestA");
  EXPECT_TRUE(labels[1] == "TestB");
  EXPECT_TRUE(labels[2] == "TestC");
  EXPECT_TRUE(labels[3] == "TestD");
}

} /*namespace Serenity*/
