/**
 * @file DelleySurfaceConstructor_test.cpp
 *
 * @author Moritz Bensberg
 * @date Jun 10, 2020
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
#include "geometry/MolecularSurfaceController.h"      //Surface is constructed via the controller.
#include "settings/PCMSettings.h"                     //Change PCMSettings.
#include "settings/Settings.h"                        //Change the settings of a test system.
#include "system/SystemController.h"                  //getSettings.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test resources
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class DelleySurfaceConstructorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests DelleySurfaceConstructor.h: surface construction for a single atom.
 */
TEST_F(DelleySurfaceConstructorTest, SurfaceConstruction_F) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs);
  PCMSettings pcmSettings;
  pcmSettings.use = true;
  pcmSettings.cavity = Options::PCM_CAVITY_TYPES::DELLEY;
  pcmSettings.lLarge = 3;
  pcmSettings.alpha = 50;
  auto surfaceController = std::make_shared<MolecularSurfaceController>(system->getGeometry(), pcmSettings);
  double sphereRadius = 3.33348;
  double sphericalArea = 4.0 * M_PI * sphereRadius * sphereRadius;
  // The total area covered by the points has to represent the total sphere surface.
  EXPECT_NEAR(surfaceController->getWeights().sum(), sphericalArea, 5e-3);
  // The first points are on the axes of the coordinate system.
  // Parallelization is done atom-wise. Thus, the ordering of points is
  // independent from the number of threads.
  const Eigen::Matrix3Xd& normals = surfaceController->getNormalVectors();
  EXPECT_NEAR(normals(0, 0), 1, 1e-12);
  EXPECT_NEAR(normals(0, 1), -1, 1e-12);
  EXPECT_NEAR(normals(1, 2), 1, 1e-12);
  EXPECT_NEAR(normals(1, 3), -1, 1e-12);
  EXPECT_NEAR(normals(2, 4), 1, 1e-12);
  EXPECT_NEAR(normals(2, 5), -1, 1e-12);
  unsigned int nPoints = normals.cols();
  const Eigen::Matrix3Xd& points = surfaceController->getGridPoints();
  // All points have to be on the sphere.
  for (unsigned int i = 0; i < nPoints; ++i)
    EXPECT_NEAR(points.col(i).norm(), sphereRadius, 5e-5);
}

/**
 * @test
 * @brief Tests DelleySurfaceConstructor.h: surface construction for a di-atomic molecule.
 */
TEST_F(DelleySurfaceConstructorTest, SurfaceConstruction_Ar2) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Ar2_6_31Gs);
  PCMSettings pcmSettings;
  pcmSettings.use = true;
  pcmSettings.cavity = Options::PCM_CAVITY_TYPES::DELLEY;
  pcmSettings.lLarge = 10;
  pcmSettings.cacheSize = 10000;
  pcmSettings.alpha = 50;
  Settings settings = system->getSettings();
  settings.pcm = pcmSettings;
  system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Ar2_6_31Gs, settings);
  auto surfaceController = std::make_shared<MolecularSurfaceController>(system->getGeometry(), pcmSettings);
  double sphereRadius = 4.26322;
  // The total area covered by the points has to represent the total sphere surface.
  EXPECT_NEAR(surfaceController->getWeights().sum(), 324.43963853260595, 5e-3);
  // The first points are on the axes of the coordinate system.
  unsigned int nThreads = 1;
#ifdef _OPENMP
  nThreads = omp_get_max_threads();
#endif
  if (nThreads == 1) {
    const Eigen::Matrix3Xd& normals = surfaceController->getNormalVectors();
    EXPECT_NEAR(normals(0, 0), 1, 1e-12);
    EXPECT_NEAR(normals(0, 1), -1, 1e-12);
    EXPECT_NEAR(normals(1, 2), 1, 1e-12);
    EXPECT_NEAR(normals(1, 3), -1, 1e-12);
    EXPECT_NEAR(normals(2, 4), -1, 1e-12);
    const Eigen::Matrix3Xd& points = surfaceController->getGridPoints();
    // First few points have to be on the sphere.
    for (unsigned int i = 0; i < 5; ++i)
      EXPECT_NEAR(points.col(i).norm(), sphereRadius, 5e-5);
  }
}

/**
 * @test
 * @brief Tests DelleySurfaceConstructor.h: surface construction for a di-atomic molecule.
 */
TEST_F(DelleySurfaceConstructorTest, SurfaceConstruction_water) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);
  Settings settings = system->getSettings();
  PCMSettings pcmSettings;
  pcmSettings.use = true;
  pcmSettings.cavity = Options::PCM_CAVITY_TYPES::DELLEY;
  pcmSettings.lLarge = 7;
  pcmSettings.alpha = 50;
  pcmSettings.cacheSize = 10000;
  pcmSettings.patchLevel = 2;
  settings.pcm = pcmSettings;
  system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  auto surfaceController = std::make_shared<MolecularSurfaceController>(system->getGeometry(), pcmSettings);
  const double weightDELLEY = surfaceController->getWeights().sum();
  // The total area covered by the points has to represent the total sphere surface.
  EXPECT_NEAR(weightDELLEY, 167.87040537607655, 5e-3);
  // Ensure that the point coordinates are the same!
  EXPECT_NEAR(surfaceController->getGridPoints().sum(), 22.895626579964713, 1e-4);

  // Compare to GEPOL surface.
  pcmSettings.cavity = Options::PCM_CAVITY_TYPES::GEPOL_SES;
  surfaceController = std::make_shared<MolecularSurfaceController>(system->getGeometry(), pcmSettings);
  const double weightGEPOL = surfaceController->getWeights().sum();
  EXPECT_NEAR(weightDELLEY, weightGEPOL, 15); // allow for roughly 10% deviation.
}

/**
 * @test
 * @brief Tests DelleySurfaceConstructor.h: surface construction for a di-atomic molecule.
 */
TEST_F(DelleySurfaceConstructorTest, SurfaceConstruction_C60) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS);
  Settings settings = system->getSettings();
  PCMSettings pcmSettings;
  pcmSettings.use = true;
  pcmSettings.cavity = Options::PCM_CAVITY_TYPES::DELLEY;
  pcmSettings.lLarge = 7;
  pcmSettings.alpha = 50;
  pcmSettings.minDistance = 0.1;
  pcmSettings.oneCavity = true;
  settings.pcm = pcmSettings;
  system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS, settings);
  auto surfaceController = std::make_shared<MolecularSurfaceController>(system->getGeometry(), pcmSettings);
  const double weightDELLEY = surfaceController->getWeights().sum();
  const unsigned int nPoints = surfaceController->getWeights().size();

  EXPECT_NEAR(weightDELLEY, 1190.23, 5e-2);
  EXPECT_EQ(nPoints, 1176);
}

} /* namespace Serenity */
