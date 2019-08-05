/**
 * @file AtomicDensityGuessCalculator_test.cpp
 *
 * @date Mar 17, 2017
 * @author David Schnieders
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
#include "scf/initialGuess/AtomicDensityGuessCalculator.h"
#include "data/ElectronicStructure.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>



namespace Serenity {
class AtomicDensityGuessCalculatorTest : public ::testing::Test {
 protected:
  AtomicDensityGuessCalculatorTest() {
  }

  virtual ~AtomicDensityGuessCalculatorTest() = default;

  static void TearDownTestCase() {
   SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(AtomicDensityGuessCalculatorTest, O2MinBasRes) {

  auto sys =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_SING,true);

  AtomicDensityGuessCalculator guessCalc(GUESSMODES::OCCUPATIONS);

  auto elecStruct=guessCalc.calculateInitialGuess(sys);

  auto densMat=elecStruct->getDensityMatrix();


  EXPECT_NEAR(densMat(0,0),2.1049096637428,1e-5);
  EXPECT_NEAR(densMat(1,0),-0.45137903089973058,1e-5);
  EXPECT_NEAR(densMat(3,3),0.89482551947104338,1e-5);
  EXPECT_NEAR(densMat(5,0),0.0021267873100875665,1e-5);
  EXPECT_NEAR(densMat(0,2),0.0,1e-5);
}

TEST_F(AtomicDensityGuessCalculatorTest, FRes) {

  auto sys =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs,true);

  AtomicDensityGuessCalculator guessCalc(GUESSMODES::OCCUPATIONS);

  auto elecStruct=guessCalc.calculateInitialGuess(sys);

  auto densMat=elecStruct->getDensityMatrix();

  EXPECT_NEAR(densMat(0,0),2.092110866061641,1e-5);
  EXPECT_NEAR(densMat(1,0),-0.17616618711117527,1e-5);
  EXPECT_NEAR(densMat(3,3),0.81573683353275817,1e-5);
  EXPECT_NEAR(densMat(5,0),0.0,1e-5);
  EXPECT_NEAR(densMat(0,2),-0.30688103014939527,1e-5);
}



}
