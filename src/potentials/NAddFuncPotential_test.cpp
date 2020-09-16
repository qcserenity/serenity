/**
 * @file NAddFuncPotential_test.cpp
 *
 * @date Dec 5, 2016
 * @author: Kevin Klahr
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
#include "potentials/NAddFuncPotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/wrappers/Libint.h"
#include "potentials/Potential.h"
#include "settings/DFTOptions.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class NAddFuncPotentialTest : public ::testing::Test {
 protected:
  NAddFuncPotentialTest()
    : systemControllerAct(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerEnv(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)) {
  }

  virtual ~NAddFuncPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerAct;
  std::shared_ptr<SystemController> systemControllerEnv;
};

/**
 * @test NAddFuncPotentialTest
 * @brief Tests the LDA Fock Matrix of an H2 dimer
 */
TEST_F(NAddFuncPotentialTest, H2Dimer_FockMatrix_LDA) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMatAct =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMatEnv =
      systemControllerEnv->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> dMatsEnv;
  dMatsEnv.push_back(dMatEnv);

  std::shared_ptr<GridController> grid = systemControllerAct->getGridController();

  NAddFuncPotential<Options::SCF_MODES::RESTRICTED> NAddPot(
      systemControllerAct, dMatAct, dMatsEnv, grid,
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = NAddPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), -0.0032064208773879719, 1e-5);
  EXPECT_NEAR(F(1, 0), -0.0043578769550568185, 1e-5);
  EXPECT_NEAR(F(2, 0), -0.010168191220511096, 1e-5);
  EXPECT_NEAR(F(0, 1), -0.0043578769550568185, 1e-5);
  EXPECT_NEAR(F(0, 2), -0.010168191220511096, 1e-5);
  EXPECT_NEAR(F(0, 3), -0.0061493217326859771, 1e-5);
}
/**
 * @test NAddFuncPotentialTest
 * @brief Tests the LDA Fock Matrix of an unrestricted H2 dimer
 */
TEST_F(NAddFuncPotentialTest, H2Dimer_FockMatrix_LDA_UNRES) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatAct =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatEnv =
      systemControllerEnv->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> dMatsEnv;
  dMatsEnv.push_back(dMatEnv);

  std::shared_ptr<GridController> grid = systemControllerAct->getGridController();

  NAddFuncPotential<Options::SCF_MODES::UNRESTRICTED> NAddPot(
      systemControllerAct, dMatAct, dMatsEnv, grid,
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(NAddPot.getMatrix()));

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), -0.0033139248081612577, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), -0.0042599743144874028, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), -0.0092449849732655735, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), -0.0042599743144874028, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), -0.0092449849732655735, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), -0.0058283797137623104, 1e-5);
  EXPECT_NEAR(F.beta(0, 0), -0.0033139248081612577, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), -0.0042599743144874028, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), -0.0092449849732655735, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), -0.0042599743144874028, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), -0.0092449849732655735, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), -0.0058283797137623104, 1e-5);
}
/**
 * @test NAddFuncPotentialTest
 * @brief Tests the GGA Fock Matrix of an H2 dimer
 */
TEST_F(NAddFuncPotentialTest, H2Dimer_FockMatrix_GGA) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMatAct =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMatEnv =
      systemControllerEnv->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> dMatsEnv;
  dMatsEnv.push_back(dMatEnv);

  std::shared_ptr<GridController> grid = systemControllerAct->getGridController();

  NAddFuncPotential<Options::SCF_MODES::RESTRICTED> NAddPotGGA(
      systemControllerAct, dMatAct, dMatsEnv, grid,
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  FockMatrix<Options::SCF_MODES::RESTRICTED> FGGA = NAddPotGGA.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(FGGA(0, 0), -0.0021847667852869456, 1e-5);
  EXPECT_NEAR(FGGA(1, 0), -0.0033354523264519293, 1e-5);
  EXPECT_NEAR(FGGA(2, 0), -0.0080339760134858357, 1e-5);
  EXPECT_NEAR(FGGA(0, 1), -0.0033354523264519293, 1e-5);
  EXPECT_NEAR(FGGA(0, 2), -0.0080339760134858374, 1e-5);
  EXPECT_NEAR(FGGA(0, 3), -0.0049539349223450451, 1e-5);
}

/**
 * @test NAddFuncPotentialTest
 * @brief Tests the GGA Fock Matrix of an unrestricted H2 dimer
 */
TEST_F(NAddFuncPotentialTest, H2Dimer_FockMatrix_GGA_UNRES) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatAct =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatEnv =
      systemControllerEnv->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> dMatsEnv;
  dMatsEnv.push_back(dMatEnv);

  std::shared_ptr<GridController> grid = systemControllerAct->getGridController();

  NAddFuncPotential<Options::SCF_MODES::UNRESTRICTED> NAddPotGGA(
      systemControllerAct, dMatAct, dMatsEnv, grid,
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FGGA(std::move(NAddPotGGA.getMatrix()));

  // TODO create more data for the test
  EXPECT_NEAR(FGGA.alpha(0, 0), -0.002340764836631874, 1e-5);
  EXPECT_NEAR(FGGA.alpha(1, 0), -0.0033216945866467577, 1e-5);
  EXPECT_NEAR(FGGA.alpha(2, 0), -0.0073657451605040813, 1e-5);
  EXPECT_NEAR(FGGA.alpha(0, 1), -0.0033216945866467577, 1e-5);
  EXPECT_NEAR(FGGA.alpha(0, 2), -0.0073657451605040796, 1e-5);
  EXPECT_NEAR(FGGA.alpha(0, 3), -0.0047478778944491503, 1e-5);
  EXPECT_NEAR(FGGA.beta(0, 0), -0.002340764836631874, 1e-5);
  EXPECT_NEAR(FGGA.beta(1, 0), -0.0033216945866467577, 1e-5);
  EXPECT_NEAR(FGGA.beta(2, 0), -0.0073657451605040813, 1e-5);
  EXPECT_NEAR(FGGA.beta(0, 1), -0.0033216945866467577, 1e-5);
  EXPECT_NEAR(FGGA.beta(0, 2), -0.0073657451605040796, 1e-5);
  EXPECT_NEAR(FGGA.beta(0, 3), -0.0047478778944491503, 1e-5);
}
/**
 * @test CoulombInteractionPotentialTest
 * @brief Tests the LDA gradients of an H2 dimer
 */
TEST_F(NAddFuncPotentialTest, H2Dimer_Gradients_LDA) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMatAct =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMatEnv =
      systemControllerEnv->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> dMatsEnv;
  dMatsEnv.push_back(dMatEnv);

  std::shared_ptr<GridController> grid = systemControllerAct->getGridController();

  NAddFuncPotential<Options::SCF_MODES::RESTRICTED> NAddPot(
      systemControllerAct, dMatAct, dMatsEnv, grid,
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  auto result = NAddPot.getGeomGradients();

  // Reference for this test are XC Gradients as of 22.03.16
  EXPECT_NEAR(result(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(result(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(result(0, 2), -0.040343607583926161, 1e-5);
  EXPECT_NEAR(result(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(result(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(result(1, 2), -0.067425100964850535, 1e-5);
}

/**
 * @test CoulombInteractionPotentialTest
 * @brief Tests the LDA gradients of an unrestricted H2 dimer
 */
TEST_F(NAddFuncPotentialTest, H2Dimer_Gradients_LDA_UNRES) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatAct =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatEnv =
      systemControllerEnv->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> dMatsEnv;
  dMatsEnv.push_back(dMatEnv);

  std::shared_ptr<GridController> grid = systemControllerAct->getGridController();

  NAddFuncPotential<Options::SCF_MODES::UNRESTRICTED> NAddPot(
      systemControllerAct, dMatAct, dMatsEnv, grid,
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  auto result = NAddPot.getGeomGradients();

  EXPECT_NEAR(result(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(result(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(result(0, 2), -0.033562766887013505, 1e-5);
  EXPECT_NEAR(result(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(result(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(result(1, 2), -0.079896919910792702, 1e-5);
}
/**
 * @test NAddFuncPotentialTest
 * @brief Tests the GGA gradients of an unrestricted H2 dimer
 */
TEST_F(NAddFuncPotentialTest, H2Dimer_Gradients_GGA) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMatAct =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMatEnv =
      systemControllerEnv->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> dMatsEnv;
  dMatsEnv.push_back(dMatEnv);

  std::shared_ptr<GridController> grid = systemControllerAct->getGridController();

  NAddFuncPotential<Options::SCF_MODES::RESTRICTED> NAddPotGGA(
      systemControllerAct, dMatAct, dMatsEnv, grid,
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  auto result = NAddPotGGA.getGeomGradients();

  // Reference for this test are XC Gradients as of 22.03.16
  EXPECT_NEAR(result(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(result(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(result(0, 2), -0.035813083336678132, 1e-5);
  EXPECT_NEAR(result(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(result(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(result(1, 2), -0.06755537003605129, 1e-5);
}
/**
 * @test NAddFuncPotentialTest
 * @brief Tests the GGA gradients of an unrestricted H2 dimer
 */
TEST_F(NAddFuncPotentialTest, H2Dimer_Gradients_GGA_UNRES) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatAct =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatEnv =
      systemControllerEnv->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> dMatsEnv;
  dMatsEnv.push_back(dMatEnv);

  std::shared_ptr<GridController> grid = systemControllerAct->getGridController();

  NAddFuncPotential<Options::SCF_MODES::UNRESTRICTED> NAddPotGGA(
      systemControllerAct, dMatAct, dMatsEnv, grid,
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  auto result = NAddPotGGA.getGeomGradients();

  // Reference for this test are XC Gradients as of 22.03.16
  EXPECT_NEAR(result(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(result(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(result(0, 2), -0.029823046677914685, 1e-5);
  EXPECT_NEAR(result(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(result(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(result(1, 2), -0.07968569370351801, 1e-5);
}

} // namespace Serenity
