/**
 * @file PCMPotential_test.cpp
 *
 * @author Moritz Bensberg
 * @date Apr 28, 2020
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
#include "potentials/PCMPotential.h"
#include "data/ElectronicStructure.h"
#include "data/grid/ElectrostaticPotentialOnGridController.h"
#include "geometry/MolecularSurfaceController.h" //Enums
#include "settings/Settings.h"                   //Member pcm.
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class PCMPotentialTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(PCMPotentialTest, WaterInWater) {
  auto water = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = water->getSettings();
  settings.pcm.use = true;
  settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  settings.pcm.patchLevel = 0;
  settings.pcm.cacheSize = 30;
  settings.pcm.solverType = Options::PCM_SOLVER_TYPES::CPCM;
  water = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);
  auto surface = water->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE);
  auto vdwSurface = water->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE_VDW);
  auto elecPot = water->getElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(MOLECULAR_SURFACE_TYPES::ACTIVE);
  PCMPotential<RESTRICTED> pcmCPCM(settings.pcm, water->getBasisController(), water->getGeometry(), surface, vdwSurface, elecPot);
  auto f_CPCM = pcmCPCM.getMatrix();

  settings.pcm.solverType = Options::PCM_SOLVER_TYPES::IEFPCM;
  PCMPotential<RESTRICTED> pcmIEFPCM(settings.pcm, water->getBasisController(), water->getGeometry(), surface,
                                     vdwSurface, elecPot);
  auto f_IEFPCM = pcmIEFPCM.getMatrix();

  auto P = water->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  EXPECT_NEAR(0.0, std::abs(pcmIEFPCM.getEnergy(P) - pcmCPCM.getEnergy(P)), 1e-4);
  EXPECT_NEAR(0.0, std::abs(pcmIEFPCM.getEnergy(P) - pcmIEFPCM.getTotalEnergy()), 1e-4);

  EXPECT_NEAR(f_CPCM(0, 0), 0.013779739269776362, 1e-6);
  EXPECT_NEAR(f_CPCM(0, 1), 0.0090576060584450033, 1e-6);
  EXPECT_NEAR(f_CPCM(0, 2), 0.00045871896963751513, 1e-6);
  EXPECT_NEAR(f_CPCM(0, 3), 0.0024422592475885131, 1e-6);
  EXPECT_NEAR(f_CPCM(0, 4), -0.00025688334810035818, 1e-6);
  EXPECT_NEAR(f_CPCM(0, 5), -5.6985703649095639e-05, 1e-6);
  EXPECT_NEAR(f_CPCM(0, 6), 0.0028535102484916414, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
TEST_F(PCMPotentialTest, WaterInWater_gradients) {
  auto water = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = water->getSettings();
  settings.pcm.use = true;
  settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  settings.pcm.patchLevel = 0;
  settings.pcm.cacheSize = 30;
  settings.pcm.solverType = Options::PCM_SOLVER_TYPES::CPCM;
  water = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);
  auto surface = water->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE);
  auto vdwSurface = water->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE_VDW);
  auto elecPot = water->getElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(MOLECULAR_SURFACE_TYPES::ACTIVE);
  PCMPotential<RESTRICTED> pcmCPCM(settings.pcm, water->getBasisController(), water->getGeometry(), surface, vdwSurface, elecPot);
  auto g_CPCM = pcmCPCM.getGeomGradients();

  EXPECT_NEAR(g_CPCM(0, 0), -0.0037746830207324952, 1e-6);
  EXPECT_NEAR(g_CPCM(0, 1), 3.299829468575995e-05, 1e-6);
  EXPECT_NEAR(g_CPCM(0, 2), -0.0076126377597855082, 1e-6);
  EXPECT_NEAR(g_CPCM(1, 0), -0.0083927455783924226, 1e-6);
  EXPECT_NEAR(g_CPCM(1, 1), 3.2251349925317272e-05, 1e-6);
  EXPECT_NEAR(g_CPCM(1, 2), -0.0016772755364749896, 1e-6);
  EXPECT_NEAR(g_CPCM(2, 0), 0.012167428599125193, 1e-6);
  EXPECT_NEAR(g_CPCM(2, 1), -6.5249644611078394e-05, 1e-6);
  EXPECT_NEAR(g_CPCM(2, 2), 0.0092899132962604367, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
} /* namespace Serenity */
