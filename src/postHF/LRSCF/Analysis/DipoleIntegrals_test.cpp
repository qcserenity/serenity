/**
 * @file DipoleIntegrals_test.cpp
 *
 * @date Apr 01, 2018
 * @author Niklas Niemeyer
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
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class DipoleIntegralsTest
 * @brief Sets everything up for the tests of DipoleIntegrals.h/.cpp .
 */
  class DipoleIntegralsTest : public ::testing::Test {
  protected:
    DipoleIntegralsTest() {
    }

    virtual ~DipoleIntegralsTest() = default;

    static void TearDownTestCase() {
      SystemController__TEST_SUPPLY::cleanUp();
    }
  };

/**
 * @test
 * @brief Tests DipoleIntegrals.h/.cpp.
 */
  TEST_F(DipoleIntegralsTest, res) {
    const auto SCFMode = Options::SCF_MODES::RESTRICTED;
    auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS, true);
    //Trigger SCF
    system->getElectronicStructure<SCFMode>();

    //Define settings
    LRSCFTaskSettings settings;
    std::vector<std::shared_ptr<LRSCFController<SCFMode> > > lrscf;
    lrscf.push_back(std::make_shared<LRSCFController<SCFMode> >(system,settings));
    Point origin(0.0,0.0,0.0);
    auto dipoles = std::make_shared<DipoleIntegrals<SCFMode> >(lrscf,origin);

    auto lenA = dipoles->getLengths(Options::INTEGRAL_TYPE::ANALYTICAL);
    auto lenN = dipoles->getLengths(Options::INTEGRAL_TYPE::NUMERICAL);
    Eigen::MatrixXd diffLen = (*lenA) - (*lenN);

    auto velA = dipoles->getVelocities(Options::INTEGRAL_TYPE::ANALYTICAL);
    auto velN = dipoles->getVelocities(Options::INTEGRAL_TYPE::NUMERICAL);
    Eigen::MatrixXd diffVel = (*velA) - (*velN);

    auto magA = dipoles->getMagnetics(Options::INTEGRAL_TYPE::ANALYTICAL);
    auto magN = dipoles->getMagnetics(Options::INTEGRAL_TYPE::NUMERICAL);
    Eigen::MatrixXd diffMag = (*magA) - (*magN);

    EXPECT_LT(diffLen.array().cwiseAbs().maxCoeff(),5e-6);
    EXPECT_LT(diffVel.array().cwiseAbs().maxCoeff(),5e-6);
    EXPECT_LT(diffMag.array().cwiseAbs().maxCoeff(),5e-6);
  }

    TEST_F(DipoleIntegralsTest, unres) {
    const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
    auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED, true);
    //Trigger SCF
    system->getElectronicStructure<SCFMode>();

    //Define settings
    LRSCFTaskSettings settings;
    std::vector<std::shared_ptr<LRSCFController<SCFMode> > > lrscf;
    lrscf.push_back(std::make_shared<LRSCFController<SCFMode> >(system,settings));
    Point origin(0.0,0.0,0.0);
    auto dipoles = std::make_shared<DipoleIntegrals<SCFMode> >(lrscf,origin);

    auto lenA = dipoles->getLengths(Options::INTEGRAL_TYPE::ANALYTICAL);
    auto lenN = dipoles->getLengths(Options::INTEGRAL_TYPE::NUMERICAL);
    Eigen::MatrixXd diffLen = (*lenA)- (*lenN);

    auto velA = dipoles->getVelocities(Options::INTEGRAL_TYPE::ANALYTICAL);
    auto velN = dipoles->getVelocities(Options::INTEGRAL_TYPE::NUMERICAL);
    Eigen::MatrixXd diffVel = (*velA) - (*velN);

    auto magA = dipoles->getMagnetics(Options::INTEGRAL_TYPE::ANALYTICAL);
    auto magN = dipoles->getMagnetics(Options::INTEGRAL_TYPE::NUMERICAL);
    Eigen::MatrixXd diffMag = (*magA) - (*magN);

    EXPECT_LT(diffLen.array().cwiseAbs().maxCoeff(),7.5e-6);
    EXPECT_LT(diffVel.array().cwiseAbs().maxCoeff(),7.5e-6);
    EXPECT_LT(diffMag.array().cwiseAbs().maxCoeff(),7.5e-6);
  }

} /* namespace Serenity */