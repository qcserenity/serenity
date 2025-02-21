/**
 * @file GeometryOptimizationTask_test.cpp
 *
 * @date Jan 31, 2017
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
#include "tasks/GeometryOptimizationTask.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h"
#include "parameters/Constants.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class GeometryOptimizationTest
 * @brief Sets everything up for the tests of GeometryOptimizationTask.h/.cpp .
 *
 */
class GeometryOptimizationTest : public ::testing::Test {
 protected:
  GeometryOptimizationTest() {
  }

  virtual ~GeometryOptimizationTest() = default;

  static void SetUpTestCase() {
    auto& libint = Libint::getInstance();
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 2);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 3);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 4);
  }
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
    auto& libint = Libint::getInstance();
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 2);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 3);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 4);
  }
};

/**
 * @test GeometryOptimizationTest
 * @brief Tests the geometry optimization routine with HF on an H2
 */
TEST_F(GeometryOptimizationTest, h2Hf) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.label = "DEF2-TZVP";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  std::vector<std::shared_ptr<SystemController>> activesystems;
  activesystems.push_back(systemController);
  auto geomTask = GeometryOptimizationTask<Options::SCF_MODES::RESTRICTED>(activesystems);
  geomTask.settings.rmsgradThresh = 1e-3;
  geomTask.run();
  auto path = systemController->getSystemPath();
  auto atoms = systemController->getGeometry()->getAtoms();
  double dist = (*atoms[0] - *atoms[1]).distanceToOrigin() * BOHR_TO_ANGSTROM;
  EXPECT_NEAR(dist, 0.735141, 1e-3);
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.xyz").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.trj").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.energies.res").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.settings").c_str()));
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test GeometryOptimizationTest
 * @brief Tests the geometry optimization routine with FaT on an H2 dimer
 */
TEST_F(GeometryOptimizationTest, H2DimerDftNori) {
  SystemController__TEST_SUPPLY::cleanUp();
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PW91;
  settings.basis.label = "STO-6G";
  settings.basis.makeSphericalBasis = true;
  auto activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  settings.name = "GeometryOptimizationTest_H2DimerDftNori_env";
  auto envSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, settings);

  GeometryOptimizationTask<RESTRICTED> geomTask({activeSystem}, {envSystem});
  geomTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  geomTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  geomTask.run();
  auto actCoM = activeSystem->getGeometry()->getCenterOfMass();
  auto envCoM = envSystem->getGeometry()->getCenterOfMass();
  double dist = (actCoM - envCoM).distanceToOrigin() * BOHR_TO_ANGSTROM;
  EXPECT_NEAR(dist, 2.7868424993460508, 1e-3);
  EXPECT_EQ(0, std::remove((activeSystem->getSystemPath() + "../opt.xyz").c_str()));
  EXPECT_EQ(0, std::remove((activeSystem->getSystemPath() + "../opt.trj").c_str()));
  EXPECT_EQ(0, std::remove((activeSystem->getSystemPath() + "SomeTestSystem.energies.res").c_str()));
  EXPECT_EQ(0, std::remove((activeSystem->getSystemPath() + "SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((activeSystem->getSystemPath() + "SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((activeSystem->getSystemPath() + "SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0, std::remove((activeSystem->getSystemPath() + "SomeTestSystem.xyz").c_str()));
  EXPECT_EQ(0, std::remove((envSystem->getSystemPath() + "GeometryOptimizationTest_H2DimerDftNori_env.energies.res").c_str()));
  EXPECT_EQ(0, std::remove((envSystem->getSystemPath() + "GeometryOptimizationTest_H2DimerDftNori_env.dmat.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((envSystem->getSystemPath() + "GeometryOptimizationTest_H2DimerDftNori_env.orbs.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((envSystem->getSystemPath() + "GeometryOptimizationTest_H2DimerDftNori_env.settings").c_str()));
  EXPECT_EQ(0, std::remove((envSystem->getSystemPath() + "GeometryOptimizationTest_H2DimerDftNori_env.xyz").c_str()));
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(activeSystem);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(envSystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test GeometryOptimizationTest
 * @brief Tests the geometry optimization routine with DFT/NORI on H2
 */

TEST_F(GeometryOptimizationTest, CO_DFT_no_RI) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "STO-6G";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS, settings);
  std::vector<std::shared_ptr<SystemController>> activesystems;
  activesystems.push_back(systemController);
  auto geomTask = GeometryOptimizationTask<Options::SCF_MODES::RESTRICTED>(activesystems);
  geomTask.settings.rmsgradThresh = 1e-3;
  geomTask.run();
  auto path = systemController->getSystemPath();
  auto atoms = systemController->getGeometry()->getAtoms();
  double dist = (*atoms[0] - *atoms[1]).distanceToOrigin() * BOHR_TO_ANGSTROM;
  EXPECT_NEAR(dist, 1.20372, 1e-3);
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.xyz").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.trj").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.energies.res").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.settings").c_str()));
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test GeometryOptimizationTest
 * @brief Tests the geometry optimization routine with DFT/RI on H2
 */

TEST_F(GeometryOptimizationTest, h2DftRi) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "STO-6G";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);
  std::vector<std::shared_ptr<SystemController>> activesystems;
  activesystems.push_back(systemController);
  auto geomTask = GeometryOptimizationTask<Options::SCF_MODES::RESTRICTED>(activesystems);
  geomTask.settings.rmsgradThresh = 1e-3;
  geomTask.run();
  auto path = systemController->getSystemPath();
  auto atoms = systemController->getGeometry()->getAtoms();
  double dist = (*atoms[0] - *atoms[1]).distanceToOrigin() * BOHR_TO_ANGSTROM;
  EXPECT_NEAR(dist, 0.7326501, 1e-3);
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.xyz").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.trj").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.energies.res").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0, std::remove((path + "SomeTestSystem.settings").c_str()));
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController);
  SystemController__TEST_SUPPLY::cleanUp();
}

}; /* NameSpace Serenity */
