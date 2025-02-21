/**
 * @file DOIBasedSelector_test.cpp
 *
 * @date Apr 3, 2019
 * @author Moritz Bensberg
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

#ifndef ANALYSIS_PAOSELECTION_DOIBASEDSELECTOR_TEST_CPP_
#define ANALYSIS_PAOSELECTION_DOIBASEDSELECTOR_TEST_CPP_

/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/DOIBasedSelector.h"         //To be tested.
#include "basis/AtomCenteredBasisController.h"              //AtomCenteredBasisController definition.
#include "data/ElectronicStructure.h"                       //Density matrices for PAO controller.
#include "data/OrbitalController.h"                         //Occupied orbital coefficients.
#include "data/PAOController.h"                             //PAO controller.
#include "data/grid/BasisFunctionOnGridControllerFactory.h" //Basis function values on a grid.
#include "integrals/OneElectronIntegralController.h"        //Overlap matrix.
#include "system/SystemController.h"                        //Test systems.
#include "tasks/LocalizationTask.h"                         //Orbital localization.
#include "tasks/ScfTask.h"                                  //Fresh electronic structure from SCF.
#include "testsupply/SystemController__TEST_SUPPLY.h"       //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h>    //Testing framework
#include <Eigen/Dense>      //Dense matrices.
#include <Eigen/SparseCore> //Sparse matrices.

namespace Serenity {

class DOIBasedSelectorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Test the PAO selection.
 */
TEST_F(DOIBasedSelectorTest, selectPAOs) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto systemA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto systemB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto waterDimer = *systemA + *systemB;
  ScfTask<scfMode> scfTask(waterDimer);
  scfTask.run();

  LocalizationTask locTask(waterDimer);
  locTask.settings.splitValenceAndCore = true;
  locTask.run();

  unsigned int nOcc = waterDimer->getNOccupiedOrbitals<scfMode>();
  auto cOcc = std::make_shared<Eigen::MatrixXd>(
      waterDimer->getActiveOrbitalController<scfMode>()->getCoefficients().leftCols(nOcc).eval());
  auto densPtr = std::make_shared<DensityMatrix<scfMode>>(waterDimer->getElectronicStructure<scfMode>()->getDensityMatrix());
  auto sPtr =
      std::make_shared<MatrixInBasis<scfMode>>(waterDimer->getOneElectronIntegralController()->getOverlapIntegrals());
  auto paoController = std::make_shared<PAOController>(densPtr, sPtr, 1e-6);
  auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
      waterDimer->getSettings(), waterDimer->getBasisController(), waterDimer->getGridController());
  double doiThreshold = 1e-2;
  DOIBasedSelector selector(cOcc, paoController, basFuncOnGridController, waterDimer->getAtomCenteredBasisController(),
                            Eigen::VectorXd::Constant(10, doiThreshold), 1e-5);

  auto selection = selector.selectPAOs();
  Eigen::MatrixXi test(selection->rows(), selection->cols());
  test << 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1,
      1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  int diff = (test - Eigen::MatrixXi(*selection)).array().abs().sum();
  EXPECT_EQ(diff, 0);
  // clean up
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(waterDimer);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */

#endif /* ANALYSIS_PAOSELECTION_DOIBASEDSELECTOR_TEST_CPP_ */
