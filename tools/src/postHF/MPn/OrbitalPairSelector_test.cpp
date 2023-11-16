/**
 * @file OrbitalPairSelector_test.cpp
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

/* Include Serenity Internal Headers */
#include "postHF/MPn/OrbitalPairSelector.h"           //To be tested.
#include "data/ElectronicStructure.h"                 //Get density matrix.
#include "data/OrbitalController.h"                   //Coefficient matrix.
#include "data/OrbitalPair.h"                         //Check pair type.
#include "data/PAOController.h"                       //PAOs.
#include "integrals/OneElectronIntegralController.h"  //Overlap integrals.
#include "potentials/bundles/PotentialBundle.h"       //Fock matrix.
#include "system/SystemController.h"                  //Test systems.
#include "tasks/LocalizationTask.h"                   //Orbital localization.
#include "tasks/ScfTask.h"                            //Run SCF.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class OrbitalPairSelectorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(OrbitalPairSelectorTest, selectOrbitalPairs) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto systemA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto systemB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto waterDimer = *systemA + *systemB;
  ScfTask<scfMode> scfTask(waterDimer);
  scfTask.run();
  LocalizationTask locTask(waterDimer);
  locTask.run();

  unsigned int nOcc = waterDimer->getNOccupiedOrbitals<scfMode>();
  auto densPtr = std::make_shared<DensityMatrix<scfMode>>(waterDimer->getElectronicStructure<scfMode>()->getDensityMatrix());
  auto sPtr =
      std::make_shared<MatrixInBasis<scfMode>>(waterDimer->getOneElectronIntegralController()->getOverlapIntegrals());
  auto paoController = std::make_shared<PAOController>(densPtr, sPtr, 1e-6);

  double collTh = 1e-2;
  OrbitalPairSelector selector(waterDimer, paoController);
  std::vector<std::shared_ptr<OrbitalPair>> initialPairs;
  for (unsigned int i = 0; i < nOcc; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      auto newPair = std::make_shared<OrbitalPair>(i, j, 1e-7, 1e-3, collTh);
      initialPairs.push_back(newPair);
    } // for j
  }   // for i
  unsigned int nPAOs = paoController->getNPAOs();
  auto paoToOcc = std::make_shared<Eigen::SparseMatrix<int>>(Eigen::MatrixXi::Constant(nPAOs, nOcc, 1).sparseView());
  double doiTh = 1e-2;
  auto sysPot = waterDimer->getPotentials<scfMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  auto f =
      std::make_shared<FockMatrix<scfMode>>(sysPot->getFockMatrix(*densPtr, std::make_shared<EnergyComponentController>()));
  auto screenedPairs = selector.selectOrbitalPairs(initialPairs, 1e-5, 1e-6, f, paoToOcc, doiTh);
  for (const auto& pair : screenedPairs.first) {
    EXPECT_EQ(pair->type, OrbitalPairTypes::CLOSE);
  }
  for (const auto& pair : screenedPairs.second) {
    EXPECT_GE(pair->dipolePairEnergy, pair->dipoleCollinearPairEnergy);
  }
  EXPECT_EQ(screenedPairs.second.size(), 6);
  auto dPair1 = screenedPairs.second[0];
  auto dPair2 = screenedPairs.second[1];
  auto dPair3 = screenedPairs.second[2];
  EXPECT_EQ(dPair1->i, 1);
  EXPECT_EQ(dPair1->j, 0);
  EXPECT_EQ(dPair2->i, 2);
  EXPECT_EQ(dPair2->j, 1);
  EXPECT_EQ(dPair3->i, 4);
  EXPECT_EQ(dPair3->j, 1);
  // clean up
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(waterDimer);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
