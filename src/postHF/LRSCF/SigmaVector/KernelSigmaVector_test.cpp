/**
 * @file KernelSigmaVector_test.cpp
 *
 * @date Nov 09, 2017
 * @author Michael Boeckers
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
#include "postHF/LRSCF/Kernel.h"
#include "postHF/LRSCF/SigmaVector/KernelSigmaVector.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>


namespace Serenity {

/**
 * @class KernelSigmaVectorTest
 * @brief Sets everything up for the tests of KernelSigmaVector.h/.cpp .
 */
class KernelSigmaVectorTest : public ::testing::Test {
 protected:
  KernelSigmaVectorTest() {
  }

  virtual ~KernelSigmaVectorTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(KernelSigmaVectorTest,RKernelSigmaVector) {
  //Setup test systems
  auto activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE);
  std::vector<std::shared_ptr<SystemController > > environmentSystem(1);
  environmentSystem[0] =  SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT);
  //Setup kernel
  auto kernel = std::make_shared<Kernel<Options::SCF_MODES::RESTRICTED> >(activeSystem,environmentSystem,true,false,
      Options::FUNCTIONALS::PW91,Options::FUNCTIONALS::PW91K,Options::FUNCTIONALS::PW91);
  //Setup guess vector
  Eigen::MatrixXd guess(1,1);
  guess.setOnes();
  //Calculate KernelSigmaVector
  KernelSigmaVector<Options::SCF_MODES::RESTRICTED> kernelSigma(activeSystem,guess,kernel);
  Eigen::MatrixXd sigma = kernelSigma.getSigma();
  EXPECT_NEAR(sigma(0,0), -0.33489367794665426, 5.0e-4);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(KernelSigmaVectorTest,UKernelSigmaVector) {
  //Setup test systems
  auto activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OH_MINBAS_PBE);
  std::vector<std::shared_ptr<SystemController > > environmentSystem(1);
  environmentSystem[0] =  SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT);
  //Setup kernel
  auto kernel = std::make_shared<Kernel<Options::SCF_MODES::UNRESTRICTED> >(activeSystem,environmentSystem,true,false,
      Options::FUNCTIONALS::PW91,Options::FUNCTIONALS::PW91K,Options::FUNCTIONALS::PW91);
  //Setup guess vector
  Eigen::MatrixXd guess(13,1);
  guess.setOnes();
  //Calculate KernelSigmaVector
  KernelSigmaVector<Options::SCF_MODES::UNRESTRICTED> kernelSigma(activeSystem,guess,kernel);
  Eigen::MatrixXd sigma = kernelSigma.getSigma();

  Eigen::VectorXd cSigma(13);
  cSigma << -3.502185046944e-02,
            -0.27488958342660813,
            -0.18958843574701995,
            -7.337999764487e-02,
            -7.002504978201e-02,
            -1.436720667941e-02,
            -3.762467519735e-02,
            -0.19784509797760275,
            -0.12515993262927669,
            -1.748585008936e-01,
            -0.0070690358790650998,
            -1.238100334696e-01,
            -8.188044278932e-02;

  for (unsigned int i = 0; i < 13; ++i) {
    EXPECT_NEAR(sigma(i,0), cSigma(i), 2.0e-4);
  }
  SystemController__TEST_SUPPLY::cleanUp();
}
}
