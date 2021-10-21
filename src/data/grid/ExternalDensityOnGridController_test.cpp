/**
 * @file ExternalDensityOnGridController_test.cpp
 *
 * @date    Nov 9, 2017
 * @author: Jan Unsleber
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
#include "data/grid/ExternalDensityOnGridController.h"
#include "data/ElectronicStructure.h"
#include "data/grid/DensityOnGridFactory.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class ExternalDensityOnGridControllerTest
 * @brief Sets everything up for the tests of ExternalDensityOnGridController.h/.cpp .
 */
class ExternalDensityOnGridControllerTest : public ::testing::Test {
 protected:
  ExternalDensityOnGridControllerTest() {
  }

  virtual ~ExternalDensityOnGridControllerTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests ExternalDensityOnGridController.h/.cpp:
 */
TEST_F(ExternalDensityOnGridControllerTest, All) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatc = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  auto grid = systemController->getGridController();
  auto densContDMatc = DensityOnGridFactory<RESTRICTED>::produce(dMatc, grid, 2, systemController->getSettings());

  // get/generate originald ata
  auto& ref_dens = densContDMatc->getDensityOnGrid();
  auto& ref_dgrad = densContDMatc->getDensityGradientOnGrid();
  auto& ref_dhess = densContDMatc->getDensityHessianOnGrid();

  // test constructors
  std::unique_ptr<DensityOnGrid<RESTRICTED>> d0(new DensityOnGrid<RESTRICTED>(ref_dens));
  ExternalDensityOnGridController<RESTRICTED> deriv0(d0);

  std::unique_ptr<DensityOnGrid<RESTRICTED>> d1(new DensityOnGrid<RESTRICTED>(ref_dens));
  auto g1 = makeGradientPtr<DensityOnGrid<RESTRICTED>>(grid);
  g1->x = ref_dgrad.x;
  g1->y = ref_dgrad.y;
  g1->z = ref_dgrad.z;
  ExternalDensityOnGridController<RESTRICTED> deriv1(d1, g1);

  std::unique_ptr<DensityOnGrid<RESTRICTED>> d2(new DensityOnGrid<RESTRICTED>(ref_dens));
  auto g2 = makeGradientPtr<DensityOnGrid<RESTRICTED>>(grid);
  g2->x = ref_dgrad.x;
  g2->y = ref_dgrad.y;
  g2->z = ref_dgrad.z;
  auto h2 = makeHessianPtr<DensityOnGrid<RESTRICTED>>(grid);
  h2->xx = ref_dhess.xx;
  h2->xy = ref_dhess.xy;
  h2->xz = ref_dhess.xz;
  h2->yy = ref_dhess.yy;
  h2->yz = ref_dhess.yz;
  h2->zz = ref_dhess.zz;
  ExternalDensityOnGridController<RESTRICTED> deriv2(d2, g2, h2);

  // test getters and compare data
  auto& dens0 = deriv0.getDensityOnGrid();
  auto& dens1 = deriv1.getDensityOnGrid();
  auto& dens2 = deriv2.getDensityOnGrid();
  for (unsigned int iPt = 0; iPt < ref_dens.size(); iPt++) {
    EXPECT_EQ(dens0[iPt], ref_dens[iPt]);
    EXPECT_EQ(dens1[iPt], ref_dens[iPt]);
    EXPECT_EQ(dens2[iPt], ref_dens[iPt]);
  }
  auto& grad1 = deriv1.getDensityGradientOnGrid();
  auto& grad2 = deriv2.getDensityGradientOnGrid();
  for (unsigned int iPt = 0; iPt < ref_dens.size(); iPt++) {
    EXPECT_EQ(grad1.x[iPt], ref_dgrad.x[iPt]);
    EXPECT_EQ(grad1.y[iPt], ref_dgrad.y[iPt]);
    EXPECT_EQ(grad1.z[iPt], ref_dgrad.z[iPt]);
    EXPECT_EQ(grad2.x[iPt], ref_dgrad.x[iPt]);
    EXPECT_EQ(grad2.y[iPt], ref_dgrad.y[iPt]);
    EXPECT_EQ(grad2.z[iPt], ref_dgrad.z[iPt]);
  }
  auto& hess2 = deriv2.getDensityHessianOnGrid();
  for (unsigned int iPt = 0; iPt < ref_dens.size(); iPt++) {
    EXPECT_EQ(hess2.xx[iPt], ref_dhess.xx[iPt]);
    EXPECT_EQ(hess2.xy[iPt], ref_dhess.xy[iPt]);
    EXPECT_EQ(hess2.xz[iPt], ref_dhess.xz[iPt]);
    EXPECT_EQ(hess2.yy[iPt], ref_dhess.yy[iPt]);
    EXPECT_EQ(hess2.yz[iPt], ref_dhess.yz[iPt]);
    EXPECT_EQ(hess2.zz[iPt], ref_dhess.zz[iPt]);
  }

  // notify and check that data is still the same
  deriv0.notify();
  deriv1.notify();
  deriv2.notify();

  for (unsigned int iPt = 0; iPt < ref_dens.size(); iPt++) {
    EXPECT_EQ(dens0[iPt], ref_dens[iPt]);
    EXPECT_EQ(dens1[iPt], ref_dens[iPt]);
    EXPECT_EQ(dens2[iPt], ref_dens[iPt]);
  }
  for (unsigned int iPt = 0; iPt < ref_dens.size(); iPt++) {
    EXPECT_EQ(grad1.x[iPt], ref_dgrad.x[iPt]);
    EXPECT_EQ(grad1.y[iPt], ref_dgrad.y[iPt]);
    EXPECT_EQ(grad1.z[iPt], ref_dgrad.z[iPt]);
    EXPECT_EQ(grad2.x[iPt], ref_dgrad.x[iPt]);
    EXPECT_EQ(grad2.y[iPt], ref_dgrad.y[iPt]);
    EXPECT_EQ(grad2.z[iPt], ref_dgrad.z[iPt]);
  }
  for (unsigned int iPt = 0; iPt < ref_dens.size(); iPt++) {
    EXPECT_EQ(hess2.xx[iPt], ref_dhess.xx[iPt]);
    EXPECT_EQ(hess2.xy[iPt], ref_dhess.xy[iPt]);
    EXPECT_EQ(hess2.xz[iPt], ref_dhess.xz[iPt]);
    EXPECT_EQ(hess2.yy[iPt], ref_dhess.yy[iPt]);
    EXPECT_EQ(hess2.yz[iPt], ref_dhess.yz[iPt]);
    EXPECT_EQ(hess2.zz[iPt], ref_dhess.zz[iPt]);
  }
}

} /* namespace Serenity */
