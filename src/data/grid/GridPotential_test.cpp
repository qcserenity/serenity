/**
 * @file   GridPotential_test.cpp
 * @author Jan Unsleber
 *
 * @date   May 9, 2017
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
#include "testsupply/GridController__TEST_SUPPLY.h"
#include "data/grid/GridPotential.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>


namespace Serenity {


class GridPotentialTest : public ::testing::Test {
protected:
  GridPotentialTest() :
      tGridController(
        GridController__TEST_SUPPLY::getGridController(TEST_GRID_CONTROLLERS::TINY)),
        tr(tGridController),tu(tGridController) {
    tr(0) = 0.0;
    tr(1) = 0.1;
    tr(2) = 1.0;
    tr(3) = 1.1;
    tu.alpha(0) = 0.0;
    tu.alpha(1) = 0.1;
    tu.alpha(2) = 1.0;
    tu.alpha(3) = 1.1;
    tu.beta(0) = 0.0;
    tu.beta(1) = 0.1;
    tu.beta(2) = 1.0;
    tu.beta(3) = 1.1;
  }
  virtual ~GridPotentialTest() = default;
  // The basis in which it is defined
  std::shared_ptr<GridController> tGridController;
  // The object under test
  GridPotential<RESTRICTED> tr;
  GridPotential<UNRESTRICTED> tu;
};
/**
 * @test
 * @brief Tests GridPotential.h: Normal construction.
 */
TEST_F(GridPotentialTest, GridConstruct_R) {
  GridPotential<RESTRICTED> vector(tGridController);
  // Check for zero initialization
  EXPECT_EQ(vector(0), 0.0);
  EXPECT_EQ(vector(1), 0.0);
  EXPECT_EQ(vector(2), 0.0);
  EXPECT_EQ(vector(3), 0.0);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: Copy construction.
 */
TEST_F(GridPotentialTest, CopyConstruct_R) {
  GridPotential<RESTRICTED> vector(tr);
  // Check for zero initialization
  EXPECT_EQ(vector(0), 0.0);
  EXPECT_EQ(vector(1), 0.1);
  EXPECT_EQ(vector(2), 1.0);
  EXPECT_EQ(vector(3), 1.1);
  EXPECT_EQ(tr(0), 0.0);
  EXPECT_EQ(tr(1), 0.1);
  EXPECT_EQ(tr(2), 1.0);
  EXPECT_EQ(tr(3), 1.1);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator= MatrixXd.
 */
TEST_F(GridPotentialTest, Eq_MatrixXd_R) {
  GridPotential<RESTRICTED> vector(tr);
  Eigen::VectorXd test (4);
  test << 1.0,2.0,2.0,3.0;
  vector=test;
  EXPECT_EQ(vector(0), 1.0);
  EXPECT_EQ(vector(1), 2.0);
  EXPECT_EQ(vector(2), 2.0);
  EXPECT_EQ(vector(3), 3.0);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator= GridPotential.
 */
TEST_F(GridPotentialTest, Eq_GD_1_R) {
  GridPotential<RESTRICTED> vector(tr);
  GridPotential<RESTRICTED> test(tGridController);
  test(0) = 1.0;
  test(2) = 2.0;
  test(1) = 2.0;
  test(3) = 3.0;
  vector=test;
  EXPECT_EQ(vector(0), 1.0);
  EXPECT_EQ(vector(1), 2.0);
  EXPECT_EQ(vector(2), 2.0);
  EXPECT_EQ(vector(3), 3.0);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator= GridPotential.
 */
TEST_F(GridPotentialTest, Eq_Ref_R) {
  GridPotential<RESTRICTED> vector(tr);
  vector(0)=500.0;
  const auto& test = &vector;
  EXPECT_EQ((*test)(0), 500.0);
  EXPECT_EQ((*test)(1), 0.1);
  EXPECT_EQ((*test)(2), 1.0);
  EXPECT_EQ((*test)(3), 1.1);
  EXPECT_EQ(test->getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator= Eigen3 computed expression.
 */
TEST_F(GridPotentialTest, Eq_MatProd_R) {
  GridPotential<RESTRICTED> vector(tr);
  Eigen::VectorXd test1 (4);
  test1 << 1.0,2.0,2.0,3.0;
  Eigen::VectorXd test2 (4);
  test2 << 1.0,0.0,0.0,1.0;
  vector=test1.cwiseProduct(test2);
  EXPECT_EQ(vector(0), 1.0);
  EXPECT_EQ(vector(1), 0.0);
  EXPECT_EQ(vector(2), 0.0);
  EXPECT_EQ(vector(3), 3.0);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator+= MatrixXd.
 */
TEST_F(GridPotentialTest, PE_VectorXd_R) {
  GridPotential<RESTRICTED> vector(tr);
  Eigen::VectorXd test (4);
  test << 1.0,2.0,2.0,3.0;
  vector+=test;
  EXPECT_EQ(vector(0), 1.0);
  EXPECT_EQ(vector(1), 2.1);
  EXPECT_EQ(vector(2), 3.0);
  EXPECT_EQ(vector(3), 4.1);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator-= MatrixXd.
 */
TEST_F(GridPotentialTest, ME_VectorXd_R) {
  GridPotential<RESTRICTED> vector(tr);
  Eigen::VectorXd test (4);
  test << 1.0,2.0,2.0,3.0;
  vector-=test;
  EXPECT_EQ(vector(0), -1.0);
  EXPECT_EQ(vector(1), -1.9);
  EXPECT_EQ(vector(2), -1.0);
  EXPECT_EQ(vector(3), -1.9);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator+= GridPotential.
 */
TEST_F(GridPotentialTest, PE_GD_R) {
  GridPotential<RESTRICTED> vector(tr);
  GridPotential<RESTRICTED> test (tGridController);
  test(0) = 1.0;
  test(2) = 2.0;
  test(1) = 2.0;
  test(3) = 3.0;
  vector+=test;
  EXPECT_EQ(vector(0), 1.0);
  EXPECT_EQ(vector(1), 2.1);
  EXPECT_EQ(vector(2), 3.0);
  EXPECT_EQ(vector(3), 4.1);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator-= GridPotential.
 */
TEST_F(GridPotentialTest, ME_GD_R) {
  GridPotential<RESTRICTED> vector(tr);
  GridPotential<RESTRICTED> test (tGridController);
  test(0) = 1.0;
  test(2) = 2.0;
  test(1) = 2.0;
  test(3) = 3.0;
  vector-=test;
  EXPECT_EQ(vector(0), -1.0);
  EXPECT_EQ(vector(1), -1.9);
  EXPECT_EQ(vector(2), -1.0);
  EXPECT_EQ(vector(3), -1.9);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: operator* mixed.
 */
TEST_F(GridPotentialTest, Multiply_Mixed_R) {
  GridPotential<RESTRICTED> vector(tr);
  GridPotential<RESTRICTED> test1 (tGridController);
  test1(0) = 1.0;
  test1(2) = 2.0;
  test1(1) = 2.0;
  test1(3) = 3.0;
  Eigen::VectorXd test2 (4);
  test2 << 1.0,0.0,0.0,1.0;
  Eigen::VectorXd result = test1.cwiseProduct(test2);
  EXPECT_EQ(result(0), 1.0);
  EXPECT_EQ(result(1), 0.0);
  EXPECT_EQ(result(2), 0.0);
  EXPECT_EQ(result(3), 3.0);
}
/**
 * @test
 * @brief Tests GridPotential.h: getGridController(), getNGridPoints()
 */
TEST_F(GridPotentialTest, GetGrid_R) {
  GridPotential<RESTRICTED> vector(tGridController);
  EXPECT_EQ(vector.getGridController(), tGridController);
  EXPECT_EQ(vector.getNGridPoints(), (unsigned int)4);
}
/**
 * @test
 * @brief Tests GridPotential.h: total()
 */
TEST_F(GridPotentialTest, Total_R) {
  GridPotential<RESTRICTED> vector(tr);
  auto total = vector.total();
  EXPECT_EQ(total(0), 0.0);
  EXPECT_EQ(total(1), 0.1);
  EXPECT_EQ(total(2), 1.0);
  EXPECT_EQ(total(3), 1.1);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: difference().
 */
TEST_F(GridPotentialTest, Difference_R) {
  GridPotential<RESTRICTED> vector(tr);
  auto diff = vector.difference();
  EXPECT_EQ(diff(0), 0.0);
  EXPECT_EQ(diff(1), 0.0);
  EXPECT_EQ(diff(2), 0.0);
  EXPECT_EQ(diff(3), 0.0);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: Move Construction
 */
TEST_F(GridPotentialTest, MoveConstruct_R) {
  EXPECT_NE(tr.data(), nullptr);
  GridPotential<RESTRICTED> result(std::move(tr));
  EXPECT_EQ(result(0), 0.0);
  EXPECT_EQ(result(1), 0.1);
  EXPECT_EQ(result(3), 1.1);
  EXPECT_NE(tr.data(), result.data());
  EXPECT_EQ(result.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: Move assignement operator
 */
TEST_F(GridPotentialTest, MoveAssignement_R) {
  GridPotential<RESTRICTED> result(tGridController);
  result = std::move(tr);
  EXPECT_EQ(result(0), 0.0);
  EXPECT_EQ(result(1), 0.1);
  EXPECT_EQ(result(3), 1.1);
  EXPECT_NE(tr.data(), result.data());
  EXPECT_EQ(0.0,tr(0));
  EXPECT_EQ(result.getGridController(), tGridController);
}

/**
 * @test
 * @brief Tests GridPotential.h: to cout
 */
TEST_F(GridPotentialTest, Stream_R) {
  std::cout << tr <<std::endl;
}

/* ==================
 *    UNRESTRICTED
 * ==================
 */


/**
 * @test
 * @brief Tests GridPotential.h: Normal construction
 */
TEST_F(GridPotentialTest, GridConstruct_U) {
  GridPotential<UNRESTRICTED> vector(tGridController);
  // Check for zero initialization
  EXPECT_EQ(vector.alpha(0), 0.0);
  EXPECT_EQ(vector.alpha(1), 0.0);
  EXPECT_EQ(vector.alpha(2), 0.0);
  EXPECT_EQ(vector.alpha(3), 0.0);
  EXPECT_EQ(vector.beta(0), 0.0);
  EXPECT_EQ(vector.beta(2), 0.0);
  EXPECT_EQ(vector.beta(1), 0.0);
  EXPECT_EQ(vector.beta(3), 0.0);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: Copy construction.
 */
TEST_F(GridPotentialTest, CopyConstruct_U) {
  GridPotential<UNRESTRICTED> vector(tu);
  // Check for zero initialization
  EXPECT_EQ(vector.alpha(0), 0.0);
  EXPECT_EQ(vector.alpha(1), 0.1);
  EXPECT_EQ(vector.alpha(2), 1.0);
  EXPECT_EQ(vector.alpha(3), 1.1);
  EXPECT_EQ(vector.beta(0), 0.0);
  EXPECT_EQ(vector.beta(2), 1.0);
  EXPECT_EQ(vector.beta(1), 0.1);
  EXPECT_EQ(vector.beta(3), 1.1);
  EXPECT_EQ(tu.alpha(0), 0.0);
  EXPECT_EQ(tu.alpha(1), 0.1);
  EXPECT_EQ(tu.alpha(2), 1.0);
  EXPECT_EQ(tu.alpha(3), 1.1);
  EXPECT_EQ(tu.beta(0), 0.0);
  EXPECT_EQ(tu.beta(2), 1.0);
  EXPECT_EQ(tu.beta(1), 0.1);
  EXPECT_EQ(tu.beta(3), 1.1);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: getGridController(), getNGridPoints()
 */
TEST_F(GridPotentialTest, GetGrid_U) {
  GridPotential<UNRESTRICTED> vector(tGridController);
  EXPECT_EQ(vector.getGridController(), tGridController);
  EXPECT_EQ(vector.getNGridPoints(), (unsigned int)4);
}
/**
 * @test
 * @brief Tests GridPotential.h: total()
 */
TEST_F(GridPotentialTest, Total_U) {
  GridPotential<UNRESTRICTED> vector(tu);
  auto total = vector.total();
  EXPECT_EQ(total(0), 0.0);
  EXPECT_EQ(total(1), 0.2);
  EXPECT_EQ(total(2), 2.0);
  EXPECT_EQ(total(3), 2.2);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: difference().
 */
TEST_F(GridPotentialTest, Difference_U) {
  GridPotential<UNRESTRICTED> vector(tu);
  vector.beta.setZero();
  auto diff = vector.difference();
  EXPECT_EQ(diff(0), 0.0);
  EXPECT_EQ(diff(1), 0.1);
  EXPECT_EQ(diff(2), 1.0);
  EXPECT_EQ(diff(3), 1.1);
  EXPECT_EQ(vector.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: Move assignement operator
 */
TEST_F(GridPotentialTest, MoveAssignement_U) {
  GridPotential<UNRESTRICTED> result(tGridController);
  result = std::move(tu);
  EXPECT_EQ(result.alpha(0), 0.0);
  EXPECT_EQ(result.alpha(1), 0.1);
  EXPECT_EQ(result.alpha(3), 1.1);
  EXPECT_NE(tu.alpha.data(), result.alpha.data());
  EXPECT_EQ(0.0,tu.alpha(0));
  EXPECT_EQ(result.getGridController(), tGridController);
}
/**
 * @test
 * @brief Tests GridPotential.h: Move Construction
 */
TEST_F(GridPotentialTest, MoveConstruct_U) {
  EXPECT_NE(tu.alpha.data(), nullptr);
  GridPotential<UNRESTRICTED> result(std::move(tu));
  EXPECT_EQ(result.alpha(0), 0.0);
  EXPECT_EQ(result.alpha(1), 0.1);
  EXPECT_EQ(result.alpha(3), 1.1);
  EXPECT_NE(tu.alpha.data(), result.alpha.data());
  EXPECT_EQ(result.getGridController(), tGridController);
}

/**
 * @test
 * @brief Tests GridPotential.h: to cout
 */
TEST_F(GridPotentialTest, Stream_U) {
  std::cout << tu <<std::endl;
}

} /*namespace Serenity*/
