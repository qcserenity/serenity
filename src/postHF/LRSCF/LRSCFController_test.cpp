/**
 * @file LRSCFController_test.cpp
 *
 * @date Dec 21, 2018
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
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class LRSCFControllerTest
 * @brief Sets everything up for the tests of LRSCFController.h/.cpp .
 */
  class LRSCFControllerTest : public ::testing::Test {
  protected:
    LRSCFControllerTest() {
    }

    virtual ~LRSCFControllerTest() = default;

    static void TearDownTestCase() {
      SystemController__TEST_SUPPLY::cleanUp();
    }
  };

/**
 * @test
 * @brief Tests LRSCFController.h/.cpp.
 */
  TEST_F(LRSCFControllerTest, unresROS) {
    const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
    auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

    //Define settings
    LRSCFTaskSettings settings;
    settings.setAlpha = {0,6,9,11};
    settings.setBeta = {0,3,9,10};

    //Build the LRSCFController
    auto lrscf = std::make_shared<LRSCFController<SCFMode> >(systemController,settings);

    auto orbitalEnergies = lrscf->getEigenvalues();
    auto coefficients = lrscf->getCoefficients();

    EXPECT_NEAR(coefficients.alpha(0,0),-1.847443758428e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(0,1), 3.861750062594e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(0,2), 6.232036305359e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(0,3), 1.085354432262e+00,1e-8);

    EXPECT_NEAR(coefficients.beta(1,0),-2.889703912574e-01,1e-8);
    EXPECT_NEAR(coefficients.beta(1,1), 1.966784325854e+00,1e-8);
    EXPECT_NEAR(coefficients.beta(1,2),-3.185164937658e+00,1e-8);
    EXPECT_NEAR(coefficients.beta(1,3), 9.483392887350e-01,1e-8);

    EXPECT_NEAR(orbitalEnergies.alpha(0),-5.944690830449e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(1), 1.888201704808e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(2), 2.271678208674e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(3), 3.795257888192e+00,1e-8);

    EXPECT_NEAR(orbitalEnergies.beta(0),-5.944690830449e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(1), 6.477099331389e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(2), 2.271678208674e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(3), 2.754630032187e+00,1e-8);

    //Redefine settings
    settings.setAlpha.clear();
    settings.setBeta.clear();
    settings.excludeAlpha = {1,2,3,4,5,7,8,10};
    settings.excludeBeta = {1,2,4,5,6,7,8,11};

    //Rebuild LRSCFController
    lrscf = std::make_shared<LRSCFController<SCFMode> >(systemController,settings);

    orbitalEnergies = lrscf->getEigenvalues();
    coefficients = lrscf->getCoefficients();

    EXPECT_NEAR(coefficients.alpha(0,0),-1.847443758428e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(0,1), 3.861750062594e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(0,2), 6.232036305359e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(0,3), 1.085354432262e+00,1e-8);

    EXPECT_NEAR(coefficients.beta(1,0),-2.889703912574e-01,1e-8);
    EXPECT_NEAR(coefficients.beta(1,1), 1.966784325854e+00,1e-8);
    EXPECT_NEAR(coefficients.beta(1,2),-3.185164937658e+00,1e-8);
    EXPECT_NEAR(coefficients.beta(1,3), 9.483392887350e-01,1e-8);

    EXPECT_NEAR(orbitalEnergies.alpha(0),-5.944690830449e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(1), 1.888201704808e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(2), 2.271678208674e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(3), 3.795257888192e+00,1e-8);

    EXPECT_NEAR(orbitalEnergies.beta(0),-5.944690830449e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(1), 6.477099331389e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(2), 2.271678208674e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(3), 2.754630032187e+00,1e-8);

    //Redefine settings
    settings.excludeAlpha.clear();
    settings.excludeBeta.clear();
    settings.energyInclusion = {-100.0, 0.0, 30.0, 60.0};

    //Reinitialize Orbitals
    lrscf = std::make_shared<LRSCFController<SCFMode> >(systemController,settings);

    orbitalEnergies = lrscf->getEigenvalues();
    coefficients = lrscf->getCoefficients();

    EXPECT_NEAR(coefficients.alpha(0,0),-1.847443758428e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(3,1), 4.830663075886e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(3,2), 3.315383886199e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(2,3), 3.663544556663e-01,1e-8);

    EXPECT_NEAR(coefficients.beta(0,0),-1.847443758428e-01,1e-8);
    EXPECT_NEAR(coefficients.beta(3,1), 4.830663075886e-01,1e-8);
    EXPECT_NEAR(coefficients.beta(3,2), 3.315383886199e-01,1e-8);
    EXPECT_NEAR(coefficients.beta(2,3), 3.663544556663e-01,1e-8);

    EXPECT_NEAR(orbitalEnergies.alpha(0),-5.944690830449e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(1), 1.420270210382e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(2), 1.420270210382e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(3), 1.888201704808e+00,1e-8);

    EXPECT_NEAR(orbitalEnergies.beta(0),-5.944690830449e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(1), 1.420270210382e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(2), 1.420270210382e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(3), 1.888201704808e+00,1e-8);

    settings.energyInclusion.clear();
    settings.energyExclusion = {0.0, 30.0, 60.0, 1000.0};

    //Reinitialize Orbitals
    lrscf = std::make_shared<LRSCFController<SCFMode> >(systemController,settings);

    orbitalEnergies = lrscf->getEigenvalues();
    coefficients = lrscf->getCoefficients();

    EXPECT_NEAR(coefficients.alpha(0,0),-1.847443758428e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(3,1), 4.830663075886e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(3,2), 3.315383886199e-01,1e-8);
    EXPECT_NEAR(coefficients.alpha(2,3), 3.663544556663e-01,1e-8);

    EXPECT_NEAR(coefficients.beta(0,0),-1.847443758428e-01,1e-8);
    EXPECT_NEAR(coefficients.beta(3,1), 4.830663075886e-01,1e-8);
    EXPECT_NEAR(coefficients.beta(3,2), 3.315383886199e-01,1e-8);
    EXPECT_NEAR(coefficients.beta(2,3), 3.663544556663e-01,1e-8);

    EXPECT_NEAR(orbitalEnergies.alpha(0),-5.944690830449e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(1), 1.420270210382e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(2), 1.420270210382e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.alpha(3), 1.888201704808e+00,1e-8);

    EXPECT_NEAR(orbitalEnergies.beta(0),-5.944690830449e-01,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(1), 1.420270210382e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(2), 1.420270210382e+00,1e-8);
    EXPECT_NEAR(orbitalEnergies.beta(3), 1.888201704808e+00,1e-8);

    //Reinitialize Orbitals
    lrscf = std::make_shared<LRSCFController<SCFMode> >(systemController,settings);
    
    orbitalEnergies = lrscf->getEigenvalues();
    auto orbitalIndices = lrscf->getReferenceOrbitals();
    for_spin(orbitalEnergies, orbitalIndices){
      unsigned int orbSize = orbitalEnergies_spin.size();
      unsigned int orbitalIndicesSize = orbitalIndices_spin.size();
      EXPECT_TRUE(orbSize == orbitalIndicesSize);
    };
  }

} /* namespace Serenity */


