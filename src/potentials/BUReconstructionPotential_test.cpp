/**
 * @file BUReconstructionPotential_test.cpp
 *
 * @date Mar 16, 2017
 * @author David Schnieders
 * @copyright \n
 * This file is part of the program Serenity.\n\n
 * Serenity is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.\n\n
 * Serenity is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.\n\n
 * You should have received a copy of the GNU Lesser General
 * Public License along with Serenity.
 * If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "potentials/BUReconstructionPotential.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/ElectronicStructure.h"
#include "data/grid/GridPotential.h"
#include "data/matrices/FockMatrix.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "geometry/Geometry.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "grid/GridController.h"
#include "integrals/wrappers/Libint.h"
#include "settings/DFTOptions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <memory>
namespace Serenity {

class BUReconstructionPotentialTest : public ::testing::Test {
 protected:
  BUReconstructionPotentialTest() {
  }

  virtual ~BUReconstructionPotentialTest() = default;

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
};

TEST_F(BUReconstructionPotentialTest, H2FDE) {
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  auto superSystemGeometry = std::make_shared<Geometry>();
  *superSystemGeometry += *(_activeSystem->getGeometry());
  *superSystemGeometry += *(_environmentSystem->getGeometry());

  auto supersystemGrid = AtomCenteredGridControllerFactory::produce(superSystemGeometry, _activeSystem->getSettings().grid,
                                                                    Options::GRID_PURPOSES::DEFAULT);

  BUReconstructionPotential<Options::SCF_MODES::RESTRICTED> exactNaddKin(
      _activeSystem, _activeSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(),
      {_environmentSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController()},
      supersystemGrid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE0),
      {_environmentSystem}, 1e-4, "AUG-CC-PVQZ", 1e-4, 0.995, 0, 0);

  auto& F = exactNaddKin.getMatrix();

  EXPECT_NEAR(F(0, 0), 0.040870939685189782, 1e-3);
  EXPECT_NEAR(F(1, 0), -0.0011796262790443111, 1e-3);
  EXPECT_NEAR(F(2, 0), 0.033337935169037207, 1e-3);
  EXPECT_NEAR(F(0, 1), -0.0011796262790443111, 1e-3);
  EXPECT_NEAR(F(0, 2), 0.033337935169037207, 1e-3);
  EXPECT_NEAR(F(0, 3), 0.022709478413572917, 1e-3);

  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

TEST_F(BUReconstructionPotentialTest, H2FDEUnres) {
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  auto superSystemGeometry = std::make_shared<Geometry>();
  *superSystemGeometry += *(_activeSystem->getGeometry());
  *superSystemGeometry += *(_environmentSystem->getGeometry());

  auto supersystemGrid = AtomCenteredGridControllerFactory::produce(superSystemGeometry, _activeSystem->getSettings().grid,
                                                                    Options::GRID_PURPOSES::DEFAULT);

  BUReconstructionPotential<Options::SCF_MODES::UNRESTRICTED> exactNaddKin(
      _activeSystem, _activeSystem->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController(),
      {_environmentSystem->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController()},
      supersystemGrid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86),
      {_environmentSystem}, 1e-4, "AUG-CC-PVQZ", 1e-4, 0.995, 1000, 0);

  auto& F = exactNaddKin.getMatrix();

  EXPECT_NEAR(F.alpha(0, 0), -0.19815149498216425, 1e-3);
  EXPECT_NEAR(F.alpha(1, 0), -0.1706710237780926, 1e-3);
  EXPECT_NEAR(F.alpha(2, 0), -0.079019916938337736, 1e-3);
  EXPECT_NEAR(F.alpha(0, 1), -0.1706710237780926, 1e-3);
  EXPECT_NEAR(F.alpha(0, 2), -0.079019916938337736, 1e-3);
  EXPECT_NEAR(F.alpha(0, 3), -0.11476987778480806, 1e-3);

  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz")
                  .c_str());
  std::remove((_activeSystem->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

} // namespace Serenity
