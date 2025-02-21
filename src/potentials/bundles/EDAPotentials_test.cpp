/**
 * @file   EDAPotentials_test.cpp
 *
 * @date   Sep 18, 2017
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
#include "potentials/bundles/EDAPotentials.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutputStream.h"
#include "scf/Scf.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class EDATaskTest
 * @brief Sets the EDATask test up.
 */
class EDAPotentialsTest : public ::testing::Test {
 protected:
  EDAPotentialsTest()
    : _systemA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs)),
      _systemB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs)) {
    _path = _systemA->getSystemPath() + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs/";
  }

  virtual ~EDAPotentialsTest() {
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.energies.res").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.orbs.res.h5").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.dmat.res.h5").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.energies.unres").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.orbs.unres.h5").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.dmat.unres.h5").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.settings").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.xyz").c_str());
    std::remove((_path).c_str());
  }

  static void SetUpTestCase() {
    auto& libint = Libint::getInstance();
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  }
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
    auto& libint = Libint::getInstance();
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  }

  /* The systems. */
  std::shared_ptr<SystemController> _systemA;
  std::shared_ptr<SystemController> _systemB;
  std::string _path;
};

/**
 * @test
 * @brief Tests the EDA energy component: ES
 *
 */
TEST_F(EDAPotentialsTest, EDA_rES) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ES);
  const auto& P(super->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
  pot->getFockMatrix(P, super->getElectronicStructure<SCFMode>()->getEnergyComponentController());
  EXPECT_NEAR(-152.03634092819178, pot->getEDAEnergyContribution(), 1e-6);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ES
 *
 */
TEST_F(EDAPotentialsTest, EDA_uES) {
  const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ES);
  const auto& P(super->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
  pot->getFockMatrix(P, super->getElectronicStructure<SCFMode>()->getEnergyComponentController());
  EXPECT_NEAR(-151.84251049065807, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESX
 *
 */
TEST_F(EDAPotentialsTest, EDA_rESX) {
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESX);
  const auto& P(super->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
  pot->getFockMatrix(P, super->getElectronicStructure<SCFMode>()->getEnergyComponentController());
  EXPECT_NEAR(-152.03634092719562, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESX
 *
 */
TEST_F(EDAPotentialsTest, EDA_uESX) {
  const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESX);
  const auto& P(super->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
  pot->getFockMatrix(P, super->getElectronicStructure<SCFMode>()->getEnergyComponentController());
  EXPECT_NEAR(-151.84251053401601, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESXCT
 *
 */
TEST_F(EDAPotentialsTest, EDA_rESXCT) {
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXCT);
  Scf<SCFMode>::perform(settings, super->getElectronicStructure<SCFMode>(), pot);
  EXPECT_NEAR(-152.03900546141278, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESXCT
 *
 */
TEST_F(EDAPotentialsTest, EDA_uESXCT) {
  const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXCT);
  Scf<SCFMode>::perform(settings, super->getElectronicStructure<SCFMode>(), pot);
  EXPECT_NEAR(-152.03900546035351, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESXPLX
 *
 */
TEST_F(EDAPotentialsTest, EDA_rESXPLX) {
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXPLX);
  Scf<SCFMode>::perform(settings, super->getElectronicStructure<SCFMode>(), pot);
  EXPECT_NEAR(-152.03808634721304, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESXPLX
 *
 */
TEST_F(EDAPotentialsTest, EDA_uESXPLX) {
  const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXPLX);
  Scf<SCFMode>::perform(settings, super->getElectronicStructure<SCFMode>(), pot);
  EXPECT_NEAR(-152.03808634719843, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESPL
 *
 */
TEST_F(EDAPotentialsTest, EDA_rESPL) {
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESPL);
  Scf<SCFMode>::perform(settings, super->getElectronicStructure<SCFMode>(), pot);
  EXPECT_NEAR(-152.03808634735159, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESPL
 *
 */
TEST_F(EDAPotentialsTest, EDA_uESPL) {
  const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESPL);
  Scf<SCFMode>::perform(settings, super->getElectronicStructure<SCFMode>(), pot);
  EXPECT_NEAR(-152.03808634739556, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESXEX
 *
 */
TEST_F(EDAPotentialsTest, EDA_rESXEX) {
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXEX);
  auto P(super->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
  super->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->updateOrbitals(
      pot->getFockMatrix(P, super->getElectronicStructure<SCFMode>()->getEnergyComponentController()),
      super->getOneElectronIntegralController());
  P = super->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
  pot->getFockMatrix(P, super->getElectronicStructure<SCFMode>()->getEnergyComponentController());
  EXPECT_NEAR(-152.02177018010286, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}
/**
 * @test
 * @brief Tests the EDA energy component: ESXEX
 *
 */
TEST_F(EDAPotentialsTest, EDA_uESXEX) {
  const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
  Settings settings(_systemA->getSettings());
  settings.scf.diisStartError = 0.0;
  auto super = (*_systemA) + (*_systemB);
  auto pot = std::make_shared<EDAPotentials<SCFMode>>(_systemA, _systemB, super, EDAEnergyContributions::ESXEX);
  auto P(super->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
  super->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->updateOrbitals(
      pot->getFockMatrix(P, super->getElectronicStructure<SCFMode>()->getEnergyComponentController()),
      super->getOneElectronIntegralController());
  P = super->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
  pot->getFockMatrix(P, super->getElectronicStructure<SCFMode>()->getEnergyComponentController());
  EXPECT_NEAR(-152.02177017825835, pot->getEDAEnergyContribution(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
}

} /* namespace Serenity */
