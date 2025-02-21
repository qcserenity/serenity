/**
 * @file ElectricFields_test.cpp
 *
 * @date   Oct 15, 2020
 * @author Niklas Niemeyer, Patrick Eschenbach
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
#include "data/ElectronicStructure.h"
#include "geometry/efields/EFieldPlates.h" //To be tested.
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

/**
 * @test
 * @brief Tests EField.h
 */
TEST(ElectricFieldTest, EFieldGroundstateEnergy) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = sys->getSettings();
  settings.basis.label = "def2-TZVP";

  auto sysa = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  double energya = sysa->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();

  settings.efield.use = true;
  auto sysb = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  double energyb = sysb->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();

  settings.efield.use = true;
  settings.efield.analytical = true;
  auto sysc = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  double energyc = sysc->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();

  // w/o electric field Turbomole Oct 2020 (uniform electric field, 1e-3 au along z-axis)
  EXPECT_NEAR(energya, -76.0577924498, 1e-6);
  // w/o electric field Serenity Oct 2020 (capacitor plates, default settings)
  EXPECT_NEAR(energya, -76.0577924528, 1e-6);
  // w/ electric field Turbomole Oct 2020 (uniform electric field, 1e-3 au along z-axis)
  EXPECT_NEAR(energyb, -76.0583133275, 1e-6);
  // w/ electric field Serenity Oct 2020 (capacitor plates, default settings)
  EXPECT_NEAR(energyb, -76.0583134053, 1e-6);
  // w/ electric field Turbomole Apr 2021 (uniform electric field, 1e-3 au along z-axis)
  EXPECT_NEAR(energyc, -76.0583133275, 1e-6);
  // w/ electric field Serenity Apr 2021 (uniform electric field, default settings)
  EXPECT_NEAR(energyc, -76.0583133281, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
