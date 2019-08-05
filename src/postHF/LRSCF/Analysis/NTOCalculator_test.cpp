/**
 * @file   NTOCalculator_test.cpp
 *
 * @date   Oct 10, 2017
 * @author M. Boeckers
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "tasks/LRSCFTask.h"
#include "postHF/LRSCF/Analysis/NTOCalculator.h"
#include "data/OrbitalController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

TEST(NTOCalculator, RNTO){
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs,true);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  //Perform TDHF calculation for first excited state
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(systemController,{});
  lrscf.settings.nEigen = 1;
  lrscf.settings.rpa = true;
  lrscf.run();

  //Calculate NTOs for first state: pure HOME-LUMO transition, NTOs very similar to HOMO and LUMO
  NTOCalculator<Options::SCF_MODES::RESTRICTED> ntoCalc(systemController,0.1);
  Eigen::VectorXd occNTO = ntoCalc.getOccNTOs(0).col(4);
  Eigen::VectorXd virtNTO = ntoCalc.getVirtNTOs(0).col(0);
  auto coeff = systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  //Check
  for (unsigned int i = 0; i < occNTO.rows(); ++i) {
    EXPECT_NEAR(fabs(occNTO(i)), fabs(coeff(i,4)), 1e-5);
    EXPECT_NEAR(fabs(virtNTO(i)), fabs(coeff(i,5)), 0.3);
  }
  std::remove((systemController->getSettings().path+"/NTOS/NTO1/README").c_str());
  std::remove((systemController->getSettings().path+"/NTOS/NTO1").c_str());
  std::remove((systemController->getSettings().path+"/NTOS").c_str());
  std::remove((systemController->getSettings().path+"TestSystem_WaterMonOne_6_31Gs.lrscf.res.h5").c_str());
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);
}

TEST(NTOCalculator, UNTO){
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs,true);
  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(systemController,{});
  lrscf.settings.nEigen = 1;
  lrscf.settings.rpa = true;
  lrscf.run();
  NTOCalculator<Options::SCF_MODES::UNRESTRICTED> ntoCalc(systemController,0.1);
  Eigen::VectorXd occNTO = ntoCalc.getOccNTOs(0).alpha.col(4);
  Eigen::VectorXd virtNTO = ntoCalc.getVirtNTOs(0).alpha.col(0);
  auto coeff = systemController->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
  for (unsigned int i = 0; i < occNTO.rows(); ++i) {
    EXPECT_NEAR(fabs(occNTO(i)), fabs(coeff.alpha(i,4)), 1e-5);
    EXPECT_NEAR(fabs(virtNTO(i)), fabs(coeff.alpha(i,5)), 0.3);
  }
  std::remove((systemController->getSettings().path+"/NTOS/NTO1/README").c_str());
  std::remove((systemController->getSettings().path+"/NTOS/NTO1").c_str());
  std::remove((systemController->getSettings().path+"/NTOS").c_str());
  std::remove((systemController->getSettings().path+"TestSystem_WaterMonOne_6_31Gs.lrscf.unres.h5").c_str());
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);
}

}/* namespace Serenity */
