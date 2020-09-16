/**
 * @file   SpinPolarizedData_test.cpp
 * @author Thomas Dresselhaus
 *
 * @date Aug 26, 2015
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
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <memory>
#include <vector>

namespace Serenity {

class SpinPolarizedDataTest : public ::testing::Test {
 protected:
  SpinPolarizedDataTest()
    : tRestrictedDouble(-2.3),
      tUnrestrictedDouble(1.4),
      tRestrictedIntVec({4, -3, 1}),
      tUnrestrictedIntVec(makeUnrestrictedFromPieces(std::vector<int>({-1, 3}), std::vector<int>({17, -4}))) {
    tUnrestrictedDouble.beta = -0.2;
  }
  virtual ~SpinPolarizedDataTest() = default;

  // The objects under test
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, double> tRestrictedDouble;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> tUnrestrictedDouble;
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::vector<int>> tRestrictedIntVec;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<int>> tUnrestrictedIntVec;
};
/**
 * @test
 * @brief Tests SpinPolarizedData.h: Normal Construction
 */
TEST_F(SpinPolarizedDataTest, Construct) {
  // double
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, double> rDouble(1.2);
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> uDouble(1.7);
  EXPECT_EQ(rDouble, 1.2);
  EXPECT_EQ(uDouble.alpha, 1.7);
  EXPECT_EQ(uDouble.beta, 1.7);
  // std::vector<int>
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::vector<int>> rIntVec(3, -5);
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<int>> uIntVec(4, 3);
  EXPECT_EQ(rIntVec.size(), (unsigned int)3);
  EXPECT_EQ(rIntVec[0], -5);
  EXPECT_EQ(uIntVec.alpha.size(), (unsigned int)4);
  EXPECT_EQ(uIntVec.alpha[3], 3);
  EXPECT_EQ(uIntVec.beta.size(), (unsigned int)4);
  EXPECT_EQ(uIntVec.beta[1], 3);
  // Make sure alpha and beta are independent
  uIntVec.alpha[1] = 7;
  EXPECT_EQ(uIntVec.alpha[1], 7);
  EXPECT_EQ(uIntVec.beta[1], 3);
  // Default construct
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<int>> uIntVec2;
  EXPECT_EQ(uIntVec2.alpha.size(), (unsigned int)0);
  uIntVec2.alpha.push_back(-3);
  EXPECT_EQ(uIntVec2.alpha[0], -3);
}
/**
 * @test
 * @brief Tests SpinPolarizedData.h: Unrestricted construction from (T&, T&)
 */
TEST_F(SpinPolarizedDataTest, MakeUnrestrictedFromPieces) {
  std::vector<int> a{-3, 5, 4};
  std::vector<int> b{1, 0};
  auto uIntVec = makeUnrestrictedFromPieces<std::vector<int>>(a, b);
  EXPECT_EQ(uIntVec.alpha.size(), (unsigned int)3);
  EXPECT_EQ(uIntVec.alpha[1], 5);
  EXPECT_EQ(uIntVec.beta.size(), (unsigned int)2);
  EXPECT_EQ(uIntVec.beta[0], 1);
  auto uIntVec2 = makeUnrestrictedFromPieces<std::vector<int>>({-3, 5, 4}, {1, 0});
  EXPECT_EQ(uIntVec2.alpha, uIntVec.alpha);
  EXPECT_EQ(uIntVec2.beta, uIntVec.beta);
}
/**
 * @test
 * @brief Tests SpinPolarizedData.h: Copy construction
 */
TEST_F(SpinPolarizedDataTest, CopyConstruction) {
  /*
   * Primitive
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, double> rCopy(tRestrictedDouble);
  EXPECT_EQ(rCopy, -2.3);
  // Make sure it's a copy
  rCopy = 0.7;
  EXPECT_EQ(tRestrictedDouble, -2.3);
  // Unresctricted
  tUnrestrictedDouble.beta = 9.6;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> uCopy(tUnrestrictedDouble);
  EXPECT_EQ(uCopy.alpha, 1.4);
  EXPECT_EQ(uCopy.beta, 9.6);
  // Make sure it's a copy
  uCopy.beta = 2.63;
  EXPECT_EQ(tUnrestrictedDouble.beta, 9.6);
  /*
   * vector<int>
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::vector<int>> rCopyVec(tRestrictedIntVec);
  EXPECT_EQ(rCopyVec[0], 4);
  // Make sure it's a copy
  rCopyVec[0] = -99;
  EXPECT_EQ(tRestrictedIntVec[0], 4);
  // Unresctricted
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<int>> uCopyVec(tUnrestrictedIntVec);
  EXPECT_EQ(uCopyVec.alpha[1], 3);
  EXPECT_EQ(uCopyVec.beta[1], -4);
  // Make sure it's a copy
  uCopyVec.alpha[1] = -99;
  EXPECT_EQ(tUnrestrictedIntVec.alpha[1], 3);
}
/**
 * @test
 * @brief Tests SpinPolarizedData.h: Move construction
 */
TEST_F(SpinPolarizedDataTest, MoveConstruction) {
  /*
   * Primitive
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, double> rCopy(std::move(tRestrictedDouble));
  EXPECT_EQ(rCopy, -2.3);
  // Unresctricted
  tUnrestrictedDouble.beta = 9.6;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> uCopy(std::move(tUnrestrictedDouble));
  EXPECT_EQ(uCopy.alpha, 1.4);
  EXPECT_EQ(uCopy.beta, 9.6);
  /*
   * vector<int>
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::vector<int>> rCopyVec(std::move(tRestrictedIntVec));
  EXPECT_EQ(rCopyVec[0], 4);
  // Unresctricted
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<int>> uCopyVec(std::move(tUnrestrictedIntVec));
  EXPECT_EQ(uCopyVec.alpha[1], 3);
  EXPECT_EQ(uCopyVec.beta[1], -4);
  /*
   * std::unique_ptr to ensure move instead of copy
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::unique_ptr<int>> rPtr(new int(1));
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::unique_ptr<int>> rPtrMove(std::move(rPtr));
  EXPECT_EQ(*rPtrMove, 1);
  EXPECT_EQ(rPtr, nullptr);
  auto uPtr =
      makeUnrestrictedFromPieces<std::unique_ptr<int>>(std::unique_ptr<int>(new int(2)), std::unique_ptr<int>(new int(3)));
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::unique_ptr<int>> uPtrMove(std::move(uPtr));
  EXPECT_EQ(*uPtrMove.alpha, 2);
  EXPECT_EQ(uPtr.alpha, nullptr);
  EXPECT_EQ(*uPtrMove.beta, 3);
  EXPECT_EQ(uPtr.beta, nullptr);
}
/**
 * @test
 * @brief Tests SpinPolarizedData.h: Assignment operator
 */
TEST_F(SpinPolarizedDataTest, Assignment) {
  /*
   * Primitive
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, double> rCopy(6.4);
  rCopy = tRestrictedDouble;
  EXPECT_EQ(rCopy, -2.3);
  // Make sure it's a copy
  rCopy = 0.7;
  EXPECT_EQ(tRestrictedDouble, -2.3);
  // Unresctricted
  tUnrestrictedDouble.beta = 9.6;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> uCopy(-1.7);
  uCopy = tUnrestrictedDouble;
  EXPECT_EQ(uCopy.alpha, 1.4);
  EXPECT_EQ(uCopy.beta, 9.6);
  // Make sure it's a copy
  uCopy.beta = 2.63;
  EXPECT_EQ(tUnrestrictedDouble.beta, 9.6);
  /*
   * vector<int>
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::vector<int>> rCopyVec(0);
  rCopyVec = tRestrictedIntVec;
  EXPECT_EQ(rCopyVec[0], 4);
  // Make sure it's a copy
  rCopyVec[0] = -99;
  EXPECT_EQ(tRestrictedIntVec[0], 4);
  // Unresctricted
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<int>> uCopyVec(0);
  uCopyVec = tUnrestrictedIntVec;
  EXPECT_EQ(uCopyVec.alpha[1], 3);
  EXPECT_EQ(uCopyVec.beta[1], -4);
  // Make sure it's a copy
  uCopyVec.alpha[1] = -99;
  EXPECT_EQ(tUnrestrictedIntVec.alpha[1], 3);
}
/**
 * @test
 * @brief Tests SpinPolarizedData.h: Move assignment
 */
TEST_F(SpinPolarizedDataTest, MoveAssignment) {
  /*
   * Primitive
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, double> rCopy(-7.44);
  rCopy = std::move(tRestrictedDouble);
  EXPECT_EQ(rCopy, -2.3);
  // Unresctricted
  tUnrestrictedDouble.beta = 9.6;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> uCopy(2.3);
  uCopy = std::move(tUnrestrictedDouble);
  EXPECT_EQ(uCopy.alpha, 1.4);
  EXPECT_EQ(uCopy.beta, 9.6);
  /*
   * vector<int>
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::vector<int>> rCopyVec(0);
  rCopyVec = std::move(tRestrictedIntVec);
  EXPECT_EQ(rCopyVec[0], 4);
  // Unresctricted
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::vector<int>> uCopyVec(0);
  uCopyVec = std::move(tUnrestrictedIntVec);
  EXPECT_EQ(uCopyVec.alpha[1], 3);
  EXPECT_EQ(uCopyVec.beta[1], -4);
  /*
   * std::unique_ptr to ensure move instead of copy
   */
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::unique_ptr<int>> rPtr(new int(1));
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::unique_ptr<int>> rPtrMove(new int(-7));
  rPtrMove = std::move(rPtr);
  EXPECT_EQ(*rPtrMove, 1);
  EXPECT_EQ(rPtr, nullptr);
  auto uPtr =
      makeUnrestrictedFromPieces<std::unique_ptr<int>>(std::unique_ptr<int>(new int(2)), std::unique_ptr<int>(new int(3)));
  auto uPtrMove = makeUnrestrictedFromPieces<std::unique_ptr<int>>(std::unique_ptr<int>(new int(66)),
                                                                   std::unique_ptr<int>(new int(666)));
  uPtrMove = std::move(uPtr);
  EXPECT_EQ(*uPtrMove.alpha, 2);
  EXPECT_EQ(uPtr.alpha, nullptr);
  EXPECT_EQ(*uPtrMove.beta, 3);
  EXPECT_EQ(uPtr.beta, nullptr);
}
/**
 * @test
 * @brief Tests SpinPolarizedData.h: total
 */
TEST_F(SpinPolarizedDataTest, Total) {
  EXPECT_EQ(tRestrictedDouble.total(), -2.3);
  EXPECT_DOUBLE_EQ(tUnrestrictedDouble.total(), 1.2);
}
/**
 * @test
 * @brief Tests SpinPolarizedData.h: arithmetic operators
 */
TEST_F(SpinPolarizedDataTest, ArithmeticOperators) {
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, double> rOther(2.0);
  EXPECT_DOUBLE_EQ(tRestrictedDouble + rOther, -0.3);
  EXPECT_DOUBLE_EQ(tRestrictedDouble - rOther, -4.3);
  EXPECT_DOUBLE_EQ(tRestrictedDouble * rOther, -4.6);
  EXPECT_DOUBLE_EQ(tRestrictedDouble / rOther, -1.15);
  // Unrestricted
  auto uOther = makeUnrestrictedFromPieces(-1.0, 2.0); // 1.4, -0.2
  EXPECT_DOUBLE_EQ((tUnrestrictedDouble + uOther).alpha, 0.4);
  EXPECT_DOUBLE_EQ((tUnrestrictedDouble + uOther).beta, 1.8);
  EXPECT_DOUBLE_EQ((tUnrestrictedDouble - uOther).alpha, 2.4);
  EXPECT_DOUBLE_EQ((tUnrestrictedDouble - uOther).beta, -2.2);
  EXPECT_DOUBLE_EQ((tUnrestrictedDouble * uOther).alpha, -1.4);
  EXPECT_DOUBLE_EQ((tUnrestrictedDouble * uOther).beta, -0.4);
  EXPECT_DOUBLE_EQ((tUnrestrictedDouble / uOther).alpha, -1.4);
  EXPECT_DOUBLE_EQ((tUnrestrictedDouble / uOther).beta, -0.1);
  // +=, *=, -=, /=
  rOther += tRestrictedDouble;
  EXPECT_DOUBLE_EQ(rOther, -0.3);
  rOther -= tRestrictedDouble;
  EXPECT_DOUBLE_EQ(rOther, 2.0);
  rOther *= tRestrictedDouble;
  EXPECT_DOUBLE_EQ(rOther, -4.6);
  rOther /= tRestrictedDouble;
  EXPECT_DOUBLE_EQ(rOther, 2.0);
  // Unrestricted
  uOther += tUnrestrictedDouble;
  EXPECT_DOUBLE_EQ(uOther.alpha, 0.4);
  EXPECT_DOUBLE_EQ(uOther.beta, 1.8);
  uOther -= tUnrestrictedDouble;
  EXPECT_DOUBLE_EQ(uOther.alpha, -1.0);
  EXPECT_DOUBLE_EQ(uOther.beta, 2.0);
  uOther *= tUnrestrictedDouble;
  EXPECT_DOUBLE_EQ(uOther.alpha, -1.4);
  EXPECT_DOUBLE_EQ(uOther.beta, -0.4);
  uOther /= tUnrestrictedDouble;
  EXPECT_DOUBLE_EQ(uOther.alpha, -1.0);
  EXPECT_DOUBLE_EQ(uOther.beta, 2.0);
}
/**
 * @test
 * @brief Tests SpinPolarizedData.h: for_spin macro
 */
TEST_F(SpinPolarizedDataTest, ForSpin) {
  for_spin(tRestrictedIntVec) {
    tRestrictedIntVec_spin[2] = 9;
  };
  EXPECT_EQ(tRestrictedIntVec[2], 9);
  for_spin(tUnrestrictedIntVec) {
    tUnrestrictedIntVec_spin[0] *= 2;
  };
  EXPECT_EQ(tUnrestrictedIntVec.alpha[0], -2);
  EXPECT_EQ(tUnrestrictedIntVec.beta[0], 34);
  // 2 Arguments
  for_spin(tRestrictedDouble, tRestrictedIntVec) {
    tRestrictedDouble_spin += tRestrictedIntVec_spin[0];
  };
  EXPECT_DOUBLE_EQ(tRestrictedDouble, 1.7);
  for_spin(tUnrestrictedDouble, tUnrestrictedIntVec) {
    tUnrestrictedDouble_spin += tUnrestrictedIntVec_spin[1];
  };
  EXPECT_DOUBLE_EQ(tUnrestrictedDouble.alpha, 4.4);
  EXPECT_DOUBLE_EQ(tUnrestrictedDouble.beta, -4.2);
  // 3-5 Arguments
  int foo = 0;
  // Restricted
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, int> r1(1);
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, int> r2(2);
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, int> r3(3);
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, int> r4(4);
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, int> r5(5);
  for_spin(r1, r2, r3) {
    foo += r1_spin + r2_spin + r3_spin;
  };
  EXPECT_EQ(foo, 6);
  foo = 0;
  for_spin(r1, r2, r3, r4) {
    foo += r1_spin + r2_spin + r3_spin + r4_spin;
  };
  EXPECT_EQ(foo, 10);
  foo = 0;
  for_spin(r1, r2, r3, r4, r5) {
    foo += r1_spin + r2_spin + r3_spin + r4_spin + r5_spin;
  };
  EXPECT_EQ(foo, 15);
  // Unrestricted
  foo = 0;
  auto u1 = makeUnrestrictedFromPieces<int>(1, 100);
  auto u2 = makeUnrestrictedFromPieces<int>(2, 200);
  auto u3 = makeUnrestrictedFromPieces<int>(3, 300);
  auto u4 = makeUnrestrictedFromPieces<int>(4, 400);
  auto u5 = makeUnrestrictedFromPieces<int>(5, 500);
  for_spin(u1, u2, u3) {
    foo += u1_spin + u2_spin + u3_spin;
  };
  EXPECT_EQ(foo, 606);
  foo = 0;
  for_spin(u1, u2, u3, u4) {
    foo += u1_spin + u2_spin + u3_spin + u4_spin;
  };
  EXPECT_EQ(foo, 1010);
  foo = 0;
  for_spin(u1, u2, u3, u4, u5) {
    foo += u1_spin + u2_spin + u3_spin + u4_spin + u5_spin;
  };
  EXPECT_EQ(foo, 1515);
  foo = 0;
}

} /* namespace Serenity */
