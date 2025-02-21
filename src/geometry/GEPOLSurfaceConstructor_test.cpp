/**
 * @file GEPOLSurfaceConstructor_test.cpp
 *
 * @date Feb 28, 2017
 * @author David Schnieders
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
#include "geometry/Geometry.h"
#include "geometry/MolecularSurfaceController.h"
#include "geometry/MolecularSurfaceFactory.h"
#include "settings/PCMSettings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class GEPOLSurfaceConstructorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(GEPOLSurfaceConstructorTest, co) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);

  PCMSettings pcmSettings;
  pcmSettings.use = true;
  pcmSettings.cavity = Options::PCM_CAVITY_TYPES::GEPOL_SES;
  auto surfaceController = std::make_shared<MolecularSurfaceController>(systemController->getGeometry(), pcmSettings);

  auto surface = surfaceController->getGridPoints();

  unsigned int shoulbBe = 70;

  unsigned int size = surface.cols();

  EXPECT_EQ(shoulbBe, size);
  EXPECT_NEAR(199.07272856925294, surfaceController->getWeights().sum(), 1e-6);
}

TEST_F(GEPOLSurfaceConstructorTest, jacobsen) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::JACOBSEN_MINBAS);

  PCMSettings pcmSettings;
  pcmSettings.use = true;
  pcmSettings.cavity = Options::PCM_CAVITY_TYPES::GEPOL_SES;
  pcmSettings.minDistance = 0.1;
  auto surfaceController = std::make_shared<MolecularSurfaceController>(systemController->getGeometry(), pcmSettings);

  auto surface = surfaceController->getGridPoints();

  unsigned int shoulbBe = 1404;

  unsigned int size = surface.cols();

  EXPECT_EQ(shoulbBe, size);
  EXPECT_NEAR(2424.98273166, surfaceController->getWeights().sum(), 1e-6);
}

} // namespace Serenity
