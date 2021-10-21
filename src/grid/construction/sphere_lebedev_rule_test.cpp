/**
 * @file   sphere_lebedev_rule_test.cpp
 *
 * @date   Mar 12, 2017
 * @author Jan Unsleber
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
#include "grid/construction/sphere_lebedev_rule.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>

namespace Serenity {

/**
 * @class AtomGridFactoryTest
 * @brief Sets everything up for the tests of the lebedev sphere .
 */
class LebedevSphereTest : public ::testing::Test {
 protected:
  LebedevSphereTest() {
  }

  virtual ~LebedevSphereTest() = default;

  std::vector<double> sphX;
  std::vector<double> sphY;
  std::vector<double> sphZ;
  std::vector<double> sphW;
};

/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0006) {
  unsigned int sphN = 6;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0006(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0014) {
  unsigned int sphN = 14;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0014(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0026) {
  unsigned int sphN = 26;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0026(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0038) {
  unsigned int sphN = 38;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0038(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0050) {
  unsigned int sphN = 50;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0050(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0086) {
  unsigned int sphN = 86;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0086(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0110) {
  unsigned int sphN = 110;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0110(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0146) {
  unsigned int sphN = 146;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0146(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0170) {
  unsigned int sphN = 170;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0170(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0194) {
  unsigned int sphN = 194;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0194(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0230) {
  unsigned int sphN = 230;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0230(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0266) {
  unsigned int sphN = 266;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0266(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0302) {
  unsigned int sphN = 302;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0302(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0350) {
  unsigned int sphN = 350;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0350(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0434) {
  unsigned int sphN = 434;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0434(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0590) {
  unsigned int sphN = 590;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0590(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0770) {
  unsigned int sphN = 770;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0770(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld0974) {
  unsigned int sphN = 974;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld0974(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld1202) {
  unsigned int sphN = 1202;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld1202(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld1454) {
  unsigned int sphN = 1454;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld1454(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld1730) {
  unsigned int sphN = 1730;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld1730(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld2030) {
  unsigned int sphN = 2030;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld2030(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld2354) {
  unsigned int sphN = 2354;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld2354(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld2702) {
  unsigned int sphN = 2702;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld2702(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld3074) {
  unsigned int sphN = 3074;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld3074(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld3470) {
  unsigned int sphN = 3470;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld3470(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld3890) {
  unsigned int sphN = 3890;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld3890(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld4334) {
  unsigned int sphN = 4334;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld4334(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld4802) {
  unsigned int sphN = 4802;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld4802(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld5294) {
  unsigned int sphN = 5294;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld5294(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
/**
 * @test
 * @brief Tests the sum of weigths of one lebedev type sphere.
 */
TEST_F(LebedevSphereTest, ld5810) {
  unsigned int sphN = 5810;
  sphX.resize(sphN);
  sphY.resize(sphN);
  sphZ.resize(sphN);
  sphW.resize(sphN);
  ld5810(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
  double sum = 0.0;
  for (unsigned int i = 0; i < sphN; ++i) {
    sum += sphW[i];
    EXPECT_NEAR(1.0, sqrt(sphX[i] * sphX[i] + sphY[i] * sphY[i] + sphZ[i] * sphZ[i]), 10e-8);
  }
  EXPECT_NEAR(1.0, sum, 10e-8);
}
} // namespace Serenity
