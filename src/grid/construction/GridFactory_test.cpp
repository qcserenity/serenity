/**
 * @file   GridFactory_test.cpp
 * @date Jun 19, 2017
 * @author Jan Unsleber
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
#include "grid/construction/GridFactory.h"
#include "grid/GridController.h"
#include "parameters/Constants.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class GridFactoryTest : public ::testing::Test {
 protected:
  inline double sphereFunction(const Eigen::VectorXd& point) const {
    return (point.norm() <= 7.0) ? 1.0 : 0.0;
  }

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
// This one just takes way too long
//                             - JU
// TEST_F(GridFactoryTest, Becke_Ahlrichs_3_7) {
//  GridFactory tGridFactory(
//          Options::GRID_TYPES::BECKE,
//          3,
//          Options::RADIAL_GRID_TYPES::AHLRICHS,
//          Options::SPHERICAL_GRID_TYPES::LEBEDEV,
//          7,
//          0.0);
//  auto sys  = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS,true);
//  GridController grid (tGridFactory.produce(sys->getGeometry()));
//  auto& pts = grid.getGridPoints();
//  auto& w = grid.getWeights();
//  double sum =0.0;
//  for (unsigned int i=0;i<w.size();i++){
//    sum += this->sphereFunction(pts[i])*w[i];
//  }
//  sum /= (343.0*4.0/3.0*PI);
//  EXPECT_NEAR(sum,1.0,4e-4);
//  EXPECT_GT(sum,1.0+3e-4);
//}

TEST_F(GridFactoryTest, Becke_Ahlrichs_3_4) {
  GridFactory tGridFactory(Options::GRID_TYPES::BECKE, 3, Options::RADIAL_GRID_TYPES::AHLRICHS,
                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 4, 0.0);
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, true);
  GridController grid(tGridFactory.produce(sys->getGeometry()));
  auto& pts = grid.getGridPoints();
  auto& w = grid.getWeights();
  double sum = 0.0;
  for (unsigned int i = 0; i < w.size(); i++) {
    sum += this->sphereFunction(pts.col(i)) * w[i];
  }
  sum /= (343.0 * 4.0 / 3.0 * PI);
  EXPECT_NEAR(sum, 1.0, 9e-4);
  EXPECT_GT(sum, 1.0 + 8e-4);
}

TEST_F(GridFactoryTest, Becke_Ahlrichs_3_2) {
  GridFactory tGridFactory(Options::GRID_TYPES::BECKE, 3, Options::RADIAL_GRID_TYPES::AHLRICHS,
                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 2, 0.0);
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, true);
  GridController grid(tGridFactory.produce(sys->getGeometry()));
  auto& pts = grid.getGridPoints();
  auto& w = grid.getWeights();
  double sum = 0.0;
  for (unsigned int i = 0; i < w.size(); i++) {
    sum += this->sphereFunction(pts.col(i)) * w[i];
  }
  sum /= (343.0 * 4.0 / 3.0 * PI);
  EXPECT_NEAR(sum, 1.0, 7e-3);
  EXPECT_GT(sum, 1.0 + 6e-3);
}

TEST_F(GridFactoryTest, Becke_Ahlrichs_3_1) {
  GridFactory tGridFactory(Options::GRID_TYPES::BECKE, 3, Options::RADIAL_GRID_TYPES::AHLRICHS,
                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 1, 0.0);
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, true);
  GridController grid(tGridFactory.produce(sys->getGeometry()));
  auto& pts = grid.getGridPoints();
  auto& w = grid.getWeights();
  double sum = 0.0;
  for (unsigned int i = 0; i < w.size(); i++) {
    sum += this->sphereFunction(pts.col(i)) * w[i];
  }
  sum /= (343.0 * 4.0 / 3.0 * PI);
  EXPECT_NEAR(sum, 1.0, 8e-3);
  EXPECT_GT(sum, 1.0 + 7e-3);
}

TEST_F(GridFactoryTest, Becke_Ahlrichs_2_4) {
  GridFactory tGridFactory(Options::GRID_TYPES::BECKE, 2, Options::RADIAL_GRID_TYPES::AHLRICHS,
                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 4, 0.0);
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, true);
  GridController grid(tGridFactory.produce(sys->getGeometry()));
  auto& pts = grid.getGridPoints();
  auto& w = grid.getWeights();
  double sum = 0.0;
  for (unsigned int i = 0; i < w.size(); i++) {
    sum += this->sphereFunction(pts.col(i)) * w[i];
  }
  sum /= (343.0 * 4.0 / 3.0 * PI);
  EXPECT_NEAR(sum, 1.0, 4e-4);
  EXPECT_GT(sum, 1.0 + 3e-4);
}

TEST_F(GridFactoryTest, Becke_Ahlrichs_1_4) {
  GridFactory tGridFactory(Options::GRID_TYPES::BECKE, 1, Options::RADIAL_GRID_TYPES::AHLRICHS,
                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 4, 0.0);
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, true);
  GridController grid(tGridFactory.produce(sys->getGeometry()));
  auto& pts = grid.getGridPoints();
  auto& w = grid.getWeights();
  double sum = 0.0;
  for (unsigned int i = 0; i < w.size(); i++) {
    sum += this->sphereFunction(pts.col(i)) * w[i];
  }
  sum /= (343.0 * 4.0 / 3.0 * PI);
  EXPECT_NEAR(sum, 1.0, 8e-5);
  EXPECT_GT(sum, 1.0 + 7e-5);
}

TEST_F(GridFactoryTest, Becke_Ahlrichs_0_4) {
  GridFactory tGridFactory(Options::GRID_TYPES::BECKE, 0, Options::RADIAL_GRID_TYPES::AHLRICHS,
                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 4, 0.0);
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, true);
  GridController grid(tGridFactory.produce(sys->getGeometry()));
  auto& pts = grid.getGridPoints();
  auto& w = grid.getWeights();
  double sum = 0.0;
  for (unsigned int i = 0; i < w.size(); i++) {
    sum += this->sphereFunction(pts.col(i)) * w[i];
  }
  sum /= (343.0 * 4.0 / 3.0 * PI);
  EXPECT_NEAR(sum, 1.0, 8e-5);
  EXPECT_GT(sum, 1.0 + 7e-5);
}

TEST_F(GridFactoryTest, Voronoi_Ahlrichs_1_4) {
  GridFactory tGridFactory(Options::GRID_TYPES::VORONOI, 0, Options::RADIAL_GRID_TYPES::AHLRICHS,
                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 4, 0.0);
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, true);
  GridController grid(tGridFactory.produce(sys->getGeometry()));
  auto& pts = grid.getGridPoints();
  auto& w = grid.getWeights();
  double sum = 0.0;
  for (unsigned int i = 0; i < w.size(); i++) {
    sum += this->sphereFunction(pts.col(i)) * w[i];
  }
  sum /= (343.0 * 4.0 / 3.0 * PI);
  EXPECT_NEAR(sum, 1.0, 2e-3);
  EXPECT_LT(sum, 1.0 - 1e-3);
}

TEST_F(GridFactoryTest, SSF_Ahlrichs_1_4) {
  GridFactory tGridFactory(Options::GRID_TYPES::SSF, 0, Options::RADIAL_GRID_TYPES::AHLRICHS,
                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 4, 0.0);
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, true);
  GridController grid(tGridFactory.produce(sys->getGeometry()));
  auto& pts = grid.getGridPoints();
  auto& w = grid.getWeights();
  double sum = 0.0;
  for (unsigned int i = 0; i < w.size(); i++) {
    sum += this->sphereFunction(pts.col(i)) * w[i];
  }
  sum /= (343.0 * 4.0 / 3.0 * PI);
  EXPECT_NEAR(sum, 1.0, 2e-3);
}

} // namespace Serenity
