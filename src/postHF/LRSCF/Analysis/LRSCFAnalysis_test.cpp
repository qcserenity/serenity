/**
 * @file   LRSCFAnalysis_test.cpp
 * @author M. Boeckers
 *
 * @date   13. November 2017
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

#include <gtest/gtest.h>
#include "postHF/LRSCF/Analysis/LRSCFAnalysis.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"


namespace Serenity {

/**
 * @test
 * @brief Tests LRSCFAnalysis.h/.cpp: tests if all functions run clean with dummy arguments.
 */
TEST(LRSCFAnalysisTest,FailTest) {
  auto activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE,true);
  unsigned int rDimension = activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>() * activeSystem->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  unsigned int uDimension = activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>().alpha * activeSystem->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>().alpha;
  uDimension += activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>().beta * activeSystem->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>().beta;
  Eigen::VectorXd eigenvalues(1);
  eigenvalues(0)=1909;
  std::vector<Eigen::MatrixXd> rEigenvectors(2,Eigen::MatrixXd::Random(rDimension,1));
  std::vector<Eigen::MatrixXd> uEigenvectors(2,Eigen::MatrixXd::Random(uDimension,1));
  LRSCFAnalysis<Options::SCF_MODES::RESTRICTED> rAnalysis(activeSystem,{},eigenvalues,rEigenvectors);
  LRSCFAnalysis<Options::SCF_MODES::UNRESTRICTED> uAnalysis(activeSystem,{},eigenvalues,uEigenvectors);
  EXPECT_NO_THROW({
    rAnalysis.printStateInfo(0);
    uAnalysis.printStateInfo(0);
    rAnalysis.printAOExcitationVectors();
    uAnalysis.printAOExcitationVectors();
    rAnalysis.mullikenPopulationAnalysis(0);
    uAnalysis.mullikenPopulationAnalysis(0);
  });
  auto path=activeSystem->getSettings().path;
  std::remove((path+"lrscf.X").c_str());
  std::remove((path+"lrscf.Y").c_str());
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE);
}
}
