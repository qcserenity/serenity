/**
 * @file OrbitalAligner_test.cpp
 *
 *  @date Mar 14, 2019
 *  @author Moritz Bensberg
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
#include "analysis/orbitalLocalization/OrbitalAligner.h"         //To be tested.
#include "analysis/populationAnalysis/IAOPopulationCalculator.h" //Check alignment by population difference.
#include "data/OrbitalController.h"                              //Access to orbital coefficients.
#include "integrals/wrappers/Libint.h"                           //Check alignment by kinetic energy difference.
#include "system/SystemController.h"                             //Test system.
#include "tasks/LocalizationTask.h"                              //Test through localization task.
#include "testsupply/SystemController__TEST_SUPPLY.h"            //Test system.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class OrbitalAlignerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(OrbitalAlignerTest, testLocalizationRestricted) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto templateSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  system->getElectronicStructure<scfMode>();
  templateSystem->getElectronicStructure<scfMode>();
  LocalizationTask locTask(templateSystem);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  LocalizationTask alignTask(system, {templateSystem});
  alignTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN;
  alignTask.settings.splitValenceAndCore = true;
  alignTask.settings.localizeVirtuals = true;
  alignTask.run();
  auto popsTem = IAOPopulationCalculator<scfMode>::calculateShellwiseOrbitalPopulations(templateSystem);
  auto popsSys = IAOPopulationCalculator<scfMode>::calculateShellwiseOrbitalPopulations(system);

  auto diff = popsTem - popsSys;
  EXPECT_LE(diff.array().abs().sum(), 2e-1);
}

TEST_F(OrbitalAlignerTest, testLocalizationUnrestricted) {
  const auto scfMode = Options::SCF_MODES::UNRESTRICTED;
  auto templateSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  system->getElectronicStructure<scfMode>();
  templateSystem->getElectronicStructure<scfMode>();
  LocalizationTask locTask(templateSystem);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  LocalizationTask alignTask(system, {templateSystem});
  alignTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN;
  alignTask.settings.splitValenceAndCore = true;
  alignTask.run();
  auto popsTem = IAOPopulationCalculator<scfMode>::calculateShellwiseOrbitalPopulations(templateSystem);
  auto popsSys = IAOPopulationCalculator<scfMode>::calculateShellwiseOrbitalPopulations(system);

  auto diff = popsTem.alpha - popsSys.alpha;
  EXPECT_LE(diff.array().abs().sum(), 2e-1);
}

TEST_F(OrbitalAlignerTest, test_kineticAlign) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto templateSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  system->getElectronicStructure<scfMode>();
  templateSystem->getElectronicStructure<scfMode>();
  LocalizationTask locTask(templateSystem);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.run();
  auto libint = Libint::getSharedPtr();
  auto templateInts = libint->compute1eInts(LIBINT_OPERATOR::kinetic, templateSystem->getBasisController());
  auto ints = libint->compute1eInts(LIBINT_OPERATOR::kinetic, system->getBasisController());

  auto coeff_template = templateSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  auto coeff = system->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  Eigen::VectorXd initialTemplateEnergies = (coeff_template.transpose() * templateInts * coeff_template).diagonal();
  Eigen::VectorXd initialSystemEnergies = (coeff.transpose() * ints * coeff).diagonal();

  SpinPolarizedData<scfMode, Eigen::VectorXi> dummyVector;
  dummyVector = Eigen::VectorXi::Constant(5, 1);
  OrbitalAligner<scfMode> aligner(system, templateSystem, 4, false, false);
  aligner.alignOrbitals(*system->getActiveOrbitalController<scfMode>(), dummyVector, dummyVector, 100000, 1e-7);
  auto popsTem = IAOPopulationCalculator<scfMode>::calculateShellwiseOrbitalPopulations(templateSystem);
  auto popsSys = IAOPopulationCalculator<scfMode>::calculateShellwiseOrbitalPopulations(system);

  coeff_template = templateSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  coeff = system->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  Eigen::VectorXd templateEnergies = (coeff_template.transpose() * templateInts * coeff_template).diagonal();
  Eigen::VectorXd systemEnergies = (coeff.transpose() * ints * coeff).diagonal();

  auto diff = popsTem - popsSys;
  Eigen::VectorXd kinDiff = templateEnergies - systemEnergies;
  EXPECT_LE(diff.array().abs().sum(), 2e-1);

  OrbitalAligner<scfMode> secondAligner(system, templateSystem, 4, true, false);
  secondAligner.alignOrbitals(*system->getActiveOrbitalController<scfMode>(), dummyVector, dummyVector, 100000, 1e-7);
  coeff_template = templateSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  coeff = system->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  templateEnergies = (coeff_template.transpose() * templateInts * coeff_template).diagonal();
  systemEnergies = (coeff.transpose() * ints * coeff).diagonal();
  Eigen::VectorXd secondKinDiff = templateEnergies - systemEnergies;
  EXPECT_LE(secondKinDiff.array().abs().sum(), 2e-1);
}

} /* namespace Serenity */
