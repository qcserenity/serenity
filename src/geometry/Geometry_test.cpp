/**
 * @file Geometry_test.cpp
 *
 * @date Jun 8, 2016
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
#include "geometry/Geometry.h"
#include "io/HDF5.h"
#include "math/Matrix.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class GeometryTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests Geometry.h: constructors.
 */
TEST_F(GeometryTest, Constructors) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto& orig = *systemController->getGeometry();

  Geometry copy(orig.getAtomSymbols(), orig.getCoordinates());
  auto copycoords = copy.getCoordinates();
  EXPECT_NEAR(0.0, copycoords(0, 0), 1e-6);
  EXPECT_NEAR(0.0, copycoords(1, 0), 1e-6);
  EXPECT_NEAR(0.0, copycoords(0, 1), 1e-6);
  EXPECT_NEAR(0.0, copycoords(1, 1), 1e-6);
  EXPECT_NEAR(-0.7, copycoords(0, 2), 1e-6);
  EXPECT_NEAR(0.7, copycoords(1, 2), 1e-6);
  EXPECT_EQ("H", copy.getAtomSymbols()[0]);
  EXPECT_EQ("H", copy.getAtomSymbols()[1]);

  Geometry copy2(orig.getAtoms());
  auto copycoords2 = copy2.getCoordinates();
  EXPECT_EQ(orig.getAtoms()[0], copy2.getAtoms()[0]);
  EXPECT_EQ(orig.getAtoms()[1], copy2.getAtoms()[1]);
  EXPECT_NEAR(0.0, copycoords2(0, 0), 1e-6);
  EXPECT_NEAR(0.0, copycoords2(1, 0), 1e-6);
  EXPECT_NEAR(0.0, copycoords2(0, 1), 1e-6);
  EXPECT_NEAR(0.0, copycoords2(1, 1), 1e-6);
  EXPECT_NEAR(-0.7, copycoords2(0, 2), 1e-6);
  EXPECT_NEAR(0.7, copycoords2(1, 2), 1e-6);
  EXPECT_EQ("H", copy2.getAtomSymbols()[0]);
  EXPECT_EQ("H", copy2.getAtomSymbols()[1]);
}

/**
 * @test
 * @brief Tests Geometry.h: atom symbols.
 */
TEST_F(GeometryTest, AtomSymbols) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto& orig = *systemController->getGeometry();
  auto symb = orig.getAtomSymbols();
  EXPECT_EQ("H", symb[0]);
  EXPECT_EQ("H", symb[1]);
}

/**
 * @test
 * @brief Tests Geometry.h: core core repulsion.
 */
TEST_F(GeometryTest, CoreCoreRepulsion) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::JACOBSEN_MINBAS);
  EXPECT_NEAR(8575.78723528, systemController->getGeometry()->getCoreCoreRepulsion(), 1e-6);
}

/**
 * @test
 * @brief Tests Geometry.h: boundaries.
 */
TEST_F(GeometryTest, Boundaries) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::JACOBSEN_MINBAS);
  auto& orig = *systemController->getGeometry();
  EXPECT_NEAR(15.37423254018665, orig.getMaxX(), 1e-6);
  EXPECT_NEAR(-16.749656568374437, orig.getMinX(), 1e-6);
  EXPECT_NEAR(17.250441151481514, orig.getMaxY(), 1e-6);
  EXPECT_NEAR(-5.9900813974518581, orig.getMinY(), 1e-6);
  EXPECT_NEAR(8.0335852419249072, orig.getMaxZ(), 1e-6);
  EXPECT_NEAR(-8.2060858546196123, orig.getMinZ(), 1e-6);
}

/**
 * @test
 * @brief Tests Geometry.h: atoms.
 */
TEST_F(GeometryTest, Atoms) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::JACOBSEN_MINBAS);
  EXPECT_EQ((unsigned int)111, systemController->getGeometry()->getNAtoms());
  auto atoms = systemController->getGeometry()->getAtoms();
  for (unsigned int i = 0; i < 111; i++) {
    EXPECT_NE(nullptr, atoms[i]);
  }
}

/**
 * @test
 * @brief Tests Geometry.h: gradients.
 */
TEST_F(GeometryTest, Gradients) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto& orig = *systemController->getGeometry();
  Matrix<double> grads(2, 3);
  grads(0, 0) = 0.0;
  grads(1, 0) = 1.0;
  grads(0, 1) = 0.1;
  grads(1, 1) = 1.1;
  grads(0, 2) = 0.2;
  grads(1, 2) = 1.2;
  orig.setGradients(grads);
  auto copygrads = orig.getGradients();
  EXPECT_NEAR(grads(0, 0), copygrads(0, 0), 1e-6);
  EXPECT_NEAR(grads(1, 0), copygrads(1, 0), 1e-6);
  EXPECT_NEAR(grads(0, 1), copygrads(0, 1), 1e-6);
  EXPECT_NEAR(grads(1, 1), copygrads(1, 1), 1e-6);
  EXPECT_NEAR(grads(0, 2), copygrads(0, 2), 1e-6);
  EXPECT_NEAR(grads(1, 2), copygrads(1, 2), 1e-6);
}

/**
 * @test
 * @brief Tests Geometry.h: coordinates.
 */
TEST_F(GeometryTest, Coordinates) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto& orig = *systemController->getGeometry();
  auto coords = orig.getCoordinates();
  EXPECT_NEAR(0.0, coords(0, 0), 1e-6);
  EXPECT_NEAR(0.0, coords(1, 0), 1e-6);
  EXPECT_NEAR(0.0, coords(0, 1), 1e-6);
  EXPECT_NEAR(0.0, coords(1, 1), 1e-6);
  EXPECT_NEAR(-0.7, coords(0, 2), 1e-6);
  EXPECT_NEAR(0.7, coords(1, 2), 1e-6);
  coords(0, 0) = 0.0;
  coords(1, 0) = 1.0;
  coords(0, 1) = 0.1;
  coords(1, 1) = 1.1;
  coords(0, 2) = 0.2;
  coords(1, 2) = 1.2;
  orig.setCoordinates(coords);
  auto copycoords = orig.getCoordinates();
  EXPECT_NEAR(coords(0, 0), copycoords(0, 0), 1e-6);
  EXPECT_NEAR(coords(1, 0), copycoords(1, 0), 1e-6);
  EXPECT_NEAR(coords(0, 1), copycoords(0, 1), 1e-6);
  EXPECT_NEAR(coords(1, 1), copycoords(1, 1), 1e-6);
  EXPECT_NEAR(coords(0, 2), copycoords(0, 2), 1e-6);
  EXPECT_NEAR(coords(1, 2), copycoords(1, 2), 1e-6);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
}

/**
 * @test
 * @brief Tests Geometry.h: addition.
 */
TEST_F(GeometryTest, Addition) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto& orig = *systemController->getGeometry();

  Geometry test;
  test += orig;
  test += orig;

  auto coords = orig.getCoordinates();
  auto copycoords = test.getCoordinates();
  EXPECT_NEAR(coords(0, 0), copycoords(0, 0), 1e-6);
  EXPECT_NEAR(coords(1, 0), copycoords(1, 0), 1e-6);
  EXPECT_NEAR(coords(0, 1), copycoords(0, 1), 1e-6);
  EXPECT_NEAR(coords(1, 1), copycoords(1, 1), 1e-6);
  EXPECT_NEAR(coords(0, 2), copycoords(0, 2), 1e-6);
  EXPECT_NEAR(coords(1, 2), copycoords(1, 2), 1e-6);
  EXPECT_NEAR(coords(0, 0), copycoords(2, 0), 1e-6);
  EXPECT_NEAR(coords(1, 0), copycoords(3, 0), 1e-6);
  EXPECT_NEAR(coords(0, 1), copycoords(2, 1), 1e-6);
  EXPECT_NEAR(coords(1, 1), copycoords(3, 1), 1e-6);
  EXPECT_NEAR(coords(0, 2), copycoords(2, 2), 1e-6);
  EXPECT_NEAR(coords(1, 2), copycoords(3, 2), 1e-6);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
}

/**
 * @test
 * @brief Tests Geometry.h: notification.
 */
TEST_F(GeometryTest, Notification) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  auto& orig = *systemController->getGeometry();
  auto coords = orig.getCoordinates();

  orig.getAtoms()[0]->setZ(0.5);
  orig.getAtoms()[1]->setZ(-0.5);

  auto copycoords = orig.getCoordinates();
  EXPECT_NEAR(coords(0, 0), copycoords(0, 0), 1e-6);
  EXPECT_NEAR(coords(1, 0), copycoords(1, 0), 1e-6);
  EXPECT_NEAR(coords(0, 1), copycoords(0, 1), 1e-6);
  EXPECT_NEAR(coords(1, 1), copycoords(1, 1), 1e-6);
  EXPECT_NEAR(0.5, copycoords(0, 2), 1e-6);
  EXPECT_NEAR(-0.5, copycoords(1, 2), 1e-6);
  EXPECT_NEAR(1.0, orig.getCoreCoreRepulsion(), 1e-6);
  EXPECT_NEAR(0.5, orig.getMaxZ(), 1e-6);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
}

/**
 * @test
 * @brief Tests Geometry.h: aligned coordinates.
 */
TEST_F(GeometryTest, AlignedCoords) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);

  auto& orig = *systemController->getGeometry();
  auto coords = orig.getAlignedCoordinates();

  EXPECT_NEAR(coords(0, 0), 0.0, 1e-6);
  EXPECT_NEAR(coords(1, 0), 0.0, 1e-6);
  EXPECT_NEAR(coords(2, 0), 0.0, 1e-6);
  EXPECT_NEAR(coords(0, 1), 0.0, 1e-6);
  EXPECT_NEAR(coords(1, 1), -1.7731092365297005, 1e-6);
  EXPECT_NEAR(coords(2, 1), 0.0, 1e-6);
  EXPECT_NEAR(coords(0, 2), 1.8361467206227224, 1e-6);
  EXPECT_NEAR(coords(1, 2), -0.47699144050518677, 1e-6);
  EXPECT_NEAR(coords(2, 2), 0.0, 1e-6);
}

/**
 * @test
 * @brief Tests Geometry.h: add and remove ghost/dummy atoms.
 */
TEST_F(GeometryTest, AddAndDeleteGhosts) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);
  auto ghosts = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs);

  auto geometry = systemController->getGeometry();
  auto ghostGeometry = ghosts->getGeometry();
  unsigned int nInitialAtoms = geometry->getNAtoms();
  geometry->addAsDummy(*ghostGeometry);
  unsigned int nAtomsWithGhosts = geometry->getNAtoms();
  geometry->deleteGhostAtoms();
  unsigned int nFinalAtoms = geometry->getNAtoms();
  EXPECT_EQ(nInitialAtoms, nFinalAtoms);
  EXPECT_EQ(6, nAtomsWithGhosts);
}

/**
 * @test
 * @brief Tests Geometry.h: get correct number of core electrons.
 */
TEST_F(GeometryTest, CoreElectrons) {
  Settings settings;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS);
  EXPECT_EQ(60 * 2, systemController->getNCoreElectrons());
}

} // namespace Serenity
