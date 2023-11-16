/**
 * @file   DispersionCorrectionCalculator_test.cpp
 *
 * @date   Nov 26, 2015
 * @author Jan Unsleber
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
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "settings/Options.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class DispersionCorrectionCalculatorTest
 * @brief Sets everything up for the tests of DispersionCorrectionCalculator.h/.cpp .
 *
 * All reference data was calculated with the fortran tool provided on the webpage of
 * Prof. Stefan Grimme: http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=DFT-D3
 *
 * Due to differences in natural constants the tests for bigger systems have larger error bars allowed.
 *
 */
class DispersionCorrectionCalculatorTest : public ::testing::Test {
 protected:
  DispersionCorrectionCalculatorTest()
    : systemControllerOne(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS)),
      systemControllerTwo(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::JACOBSEN_MINBAS)),
      systemControllerThree(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H60_Ghost_MINBAS)) {
  }

  virtual ~DispersionCorrectionCalculatorTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// system
  std::shared_ptr<SystemController> systemControllerOne;
  std::shared_ptr<SystemController> systemControllerTwo;
  std::shared_ptr<SystemController> systemControllerThree;
};

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) energy correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3EnergyBP86CO) {
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3, systemControllerOne->getGeometry(), CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.00000074, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) energy correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3EnergyBP86COGhost) {
  auto geometry = *systemControllerOne->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3,
                                                                                 std::make_shared<Geometry>(geometry),
                                                                                 CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.00000074, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) energy correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3EnergyBP86Jacobsen) {
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3, systemControllerTwo->getGeometry(), CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.19527281, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) energy correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3EnergyBP86JacobsenGhost) {
  auto geometry = *systemControllerTwo->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3,
                                                                                 std::make_shared<Geometry>(geometry),
                                                                                 CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.19527281, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC energy correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3ABCEnergyBP86CO) {
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3ABC,
                                                                                 systemControllerOne->getGeometry(),
                                                                                 CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.00000074, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC energy correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3ABCEnergyBP86Jacobsen) {
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3ABC,
                                                                                 systemControllerTwo->getGeometry(),
                                                                                 CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.19303031, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC energy correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3ABCEnergyBP86JacobsenGhost) {
  auto geometry = *systemControllerTwo->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3ABC,
                                                                                 std::make_shared<Geometry>(geometry),
                                                                                 CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.19303031, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ energy correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJEnergyBP86CO) {
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3BJ, systemControllerOne->getGeometry(), CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.00062268, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC energy correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJEnergyBP86Jacobsen) {
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3BJ, systemControllerTwo->getGeometry(), CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.28100110, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC energy correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJEnergyBP86JacobsenGhost) {
  auto geometry = *systemControllerTwo->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJ,
                                                                                 std::make_shared<Geometry>(geometry),
                                                                                 CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.28100110, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ-ABC energy correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJABCEnergyBP86CO) {
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3BJABC, systemControllerOne->getGeometry(),
      CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.00062268, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC energy correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJABCEnergyBP86Jacobsen) {
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3BJABC, systemControllerTwo->getGeometry(),
      CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.27875860, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC energy correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJABCEnergyBP86JacobsenGhost) {
  auto geometry = *systemControllerTwo->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  double energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3BJABC, std::make_shared<Geometry>(geometry),
      CompositeFunctionals::XCFUNCTIONALS::BP86);

  EXPECT_NEAR(energy, -0.27875860, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) gradient correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3GradBP86CO) {
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3,
                                                                               systemControllerOne->getGeometry(),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 2), 0.27878862999524E-5, 1E-8);
  EXPECT_NEAR(grad(1, 2), -0.27878862999524E-5, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) gradient correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3GradBP86COGhost) {
  auto geometry = *systemControllerOne->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3,
                                                                               std::make_shared<Geometry>(geometry),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 2), 0.27878862999524E-5, 1E-8);
  EXPECT_NEAR(grad(1, 2), -0.27878862999524E-5, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) gradient correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3GradBP86Jacobsen) {
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3,
                                                                               systemControllerTwo->getGeometry(),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 0), -0.23939648114107E-03, 1E-6);
  EXPECT_NEAR(grad(0, 1), -0.29326368463196E-03, 1E-6);
  EXPECT_NEAR(grad(0, 2), 0.16643267153564E-03, 1E-6);
  EXPECT_NEAR(grad(1, 0), 0.28480927202525E-03, 1E-6);
  EXPECT_NEAR(grad(1, 1), -0.66020038784406E-03, 1E-6);
  EXPECT_NEAR(grad(1, 2), 0.38634083205742E-04, 1E-6);
  EXPECT_NEAR(grad(2, 0), -0.55805381243386E-04, 1E-6);
  EXPECT_NEAR(grad(2, 1), -0.14769958781335E-03, 1E-6);
  EXPECT_NEAR(grad(2, 2), 0.43079324928880E-03, 1E-6);
  EXPECT_NEAR(grad(3, 0), -0.10159865831420E-03, 1E-6);
  EXPECT_NEAR(grad(3, 1), -0.24115292587994E-03, 1E-6);
  EXPECT_NEAR(grad(3, 2), 0.39270485334179E-03, 1E-6);
  EXPECT_NEAR(grad(4, 0), 0.64648503102299E-03, 1E-6);
  EXPECT_NEAR(grad(4, 1), 0.11014200526491E-03, 1E-6);
  EXPECT_NEAR(grad(4, 2), 0.77315756624325E-03, 1E-6);
  EXPECT_NEAR(grad(5, 0), 0.15600742018761E-03, 1E-6);
  EXPECT_NEAR(grad(5, 1), -0.83701252002179E-04, 1E-6);
  EXPECT_NEAR(grad(5, 2), -0.10320055086246E-03, 1E-6);
  EXPECT_NEAR(grad(6, 0), -0.26369629280981E-03, 1E-6);
  EXPECT_NEAR(grad(6, 1), 0.43056225680949E-03, 1E-6);
  EXPECT_NEAR(grad(6, 2), 0.25897355496868E-03, 1E-6);
  EXPECT_NEAR(grad(7, 0), 0.23727668340770E-03, 1E-6);
  EXPECT_NEAR(grad(7, 1), -0.71058464700238E-03, 1E-6);
  EXPECT_NEAR(grad(7, 2), -0.22975462470042E-03, 1E-6);
  EXPECT_NEAR(grad(8, 0), 0.51233163024399E-04, 1E-6);
  EXPECT_NEAR(grad(8, 1), -0.20086917537838E-03, 1E-6);
  EXPECT_NEAR(grad(8, 2), -0.10567080208689E-03, 1E-6);
  EXPECT_NEAR(grad(9, 0), -0.25335844672692E-03, 1E-6);
  EXPECT_NEAR(grad(9, 1), -0.51486210410658E-03, 1E-6);
  EXPECT_NEAR(grad(9, 2), 0.11292206538659E-03, 1E-6);
  EXPECT_NEAR(grad(10, 0), -0.16750617929842E-03, 1E-6);
  EXPECT_NEAR(grad(10, 1), -0.13617186678462E-03, 1E-6);
  EXPECT_NEAR(grad(10, 2), 0.65781593998751E-04, 1E-6);
  EXPECT_NEAR(grad(11, 0), -0.64870113385215E-04, 1E-6);
  EXPECT_NEAR(grad(11, 1), -0.18880389933317E-03, 1E-6);
  EXPECT_NEAR(grad(11, 2), 0.12669866738940E-03, 1E-6);
  EXPECT_NEAR(grad(12, 0), -0.41101498125640E-04, 1E-6);
  EXPECT_NEAR(grad(12, 1), -0.91089636774285E-04, 1E-6);
  EXPECT_NEAR(grad(12, 2), -0.41780393924603E-05, 1E-6);
  EXPECT_NEAR(grad(13, 0), -0.21937055610206E-03, 1E-6);
  EXPECT_NEAR(grad(13, 1), -0.21022567376362E-03, 1E-6);
  EXPECT_NEAR(grad(13, 2), -0.53355798116362E-04, 1E-6);
  EXPECT_NEAR(grad(14, 0), 0.12096030488783E-03, 1E-6);
  EXPECT_NEAR(grad(14, 1), -0.33481498535850E-03, 1E-6);
  EXPECT_NEAR(grad(14, 2), -0.20951168017895E-04, 1E-6);
  EXPECT_NEAR(grad(15, 0), 0.29819435863964E-05, 1E-6);
  EXPECT_NEAR(grad(15, 1), -0.88856288253517E-04, 1E-6);
  EXPECT_NEAR(grad(15, 2), -0.87371658347790E-05, 1E-6);
  EXPECT_NEAR(grad(16, 0), 0.16722652455443E-03, 1E-6);
  EXPECT_NEAR(grad(16, 1), -0.41481076468878E-03, 1E-6);
  EXPECT_NEAR(grad(16, 2), -0.69862591390206E-04, 1E-6);
  EXPECT_NEAR(grad(17, 0), 0.10772916885965E-03, 1E-6);
  EXPECT_NEAR(grad(17, 1), -0.16632306472675E-03, 1E-6);
  EXPECT_NEAR(grad(17, 2), -0.27498132784997E-04, 1E-6);
  EXPECT_NEAR(grad(18, 0), 0.67637761977715E-03, 1E-6);
  EXPECT_NEAR(grad(18, 1), -0.21667916154471E-03, 1E-6);
  EXPECT_NEAR(grad(18, 2), 0.58562657565709E-04, 1E-6);
  EXPECT_NEAR(grad(19, 0), -0.74736649925546E-03, 1E-6);
  EXPECT_NEAR(grad(19, 1), -0.12040683859880E-04, 1E-6);
  EXPECT_NEAR(grad(19, 2), 0.72367216058471E-03, 1E-6);
  EXPECT_NEAR(grad(20, 0), 0.33510778156219E-03, 1E-6);
  EXPECT_NEAR(grad(20, 1), -0.12848194374529E-02, 1E-6);
  EXPECT_NEAR(grad(20, 2), 0.65719560841665E-03, 1E-6);
  EXPECT_NEAR(grad(21, 0), -0.51700318671666E-03, 1E-6);
  EXPECT_NEAR(grad(21, 1), -0.12995799212020E-02, 1E-6);
  EXPECT_NEAR(grad(21, 2), 0.44596927969442E-03, 1E-6);
  EXPECT_NEAR(grad(22, 0), 0.23229288681336E-03, 1E-6);
  EXPECT_NEAR(grad(22, 1), -0.46685255032151E-04, 1E-6);
  EXPECT_NEAR(grad(22, 2), 0.15731997236800E-03, 1E-6);
  EXPECT_NEAR(grad(23, 0), -0.12877068374006E-03, 1E-6);
  EXPECT_NEAR(grad(23, 1), -0.71406503902311E-04, 1E-6);
  EXPECT_NEAR(grad(23, 2), 0.45524638085392E-03, 1E-6);
  EXPECT_NEAR(grad(24, 0), -0.79414902533742E-03, 1E-6);
  EXPECT_NEAR(grad(24, 1), 0.13031633455296E-04, 1E-6);
  EXPECT_NEAR(grad(24, 2), 0.50053818798825E-03, 1E-6);
  EXPECT_NEAR(grad(25, 0), -0.70347565420089E-04, 1E-6);
  EXPECT_NEAR(grad(25, 1), 0.13634425309191E-03, 1E-6);
  EXPECT_NEAR(grad(25, 2), 0.99452415420100E-03, 1E-6);
  EXPECT_NEAR(grad(26, 0), -0.11104642475317E-03, 1E-6);
  EXPECT_NEAR(grad(26, 1), -0.46423862927161E-04, 1E-6);
  EXPECT_NEAR(grad(26, 2), 0.58987772591003E-03, 1E-6);
  EXPECT_NEAR(grad(27, 0), -0.23657197709089E-03, 1E-6);
  EXPECT_NEAR(grad(27, 1), -0.20104935787111E-03, 1E-6);
  EXPECT_NEAR(grad(27, 2), 0.47002007086454E-03, 1E-6);
  EXPECT_NEAR(grad(28, 0), -0.24782242124356E-04, 1E-6);
  EXPECT_NEAR(grad(28, 1), -0.19150823349044E-03, 1E-6);
  EXPECT_NEAR(grad(28, 2), 0.97555104070209E-03, 1E-6);
  EXPECT_NEAR(grad(29, 0), 0.12853758383436E-03, 1E-6);
  EXPECT_NEAR(grad(29, 1), 0.19637817812301E-03, 1E-6);
  EXPECT_NEAR(grad(29, 2), -0.14691683814867E-03, 1E-6);
  EXPECT_NEAR(grad(30, 0), 0.13602991186037E-03, 1E-6);
  EXPECT_NEAR(grad(30, 1), -0.20967995868233E-03, 1E-6);
  EXPECT_NEAR(grad(30, 2), -0.56941257610763E-04, 1E-6);
  EXPECT_NEAR(grad(31, 0), 0.57719700391863E-03, 1E-6);
  EXPECT_NEAR(grad(31, 1), -0.36817858906916E-03, 1E-6);
  EXPECT_NEAR(grad(31, 2), -0.54971762179040E-05, 1E-6);
  EXPECT_NEAR(grad(32, 0), 0.46534525505210E-04, 1E-6);
  EXPECT_NEAR(grad(32, 1), -0.92278239751090E-05, 1E-6);
  EXPECT_NEAR(grad(32, 2), -0.94781232424354E-04, 1E-6);
  EXPECT_NEAR(grad(33, 0), -0.17825669853133E-03, 1E-6);
  EXPECT_NEAR(grad(33, 1), -0.74080518499352E-04, 1E-6);
  EXPECT_NEAR(grad(33, 2), -0.15910188089558E-03, 1E-6);
  EXPECT_NEAR(grad(34, 0), -0.49119782519775E-03, 1E-6);
  EXPECT_NEAR(grad(34, 1), -0.54382238997521E-04, 1E-6);
  EXPECT_NEAR(grad(34, 2), 0.37033534453372E-03, 1E-6);
  EXPECT_NEAR(grad(35, 0), -0.28727571461093E-03, 1E-6);
  EXPECT_NEAR(grad(35, 1), 0.82022566482949E-04, 1E-6);
  EXPECT_NEAR(grad(35, 2), 0.63619548621831E-03, 1E-6);
  EXPECT_NEAR(grad(36, 0), -0.39231469999683E-03, 1E-6);
  EXPECT_NEAR(grad(36, 1), -0.14956241246950E-03, 1E-6);
  EXPECT_NEAR(grad(36, 2), 0.57041515808173E-03, 1E-6);
  EXPECT_NEAR(grad(37, 0), -0.70825954266066E-04, 1E-6);
  EXPECT_NEAR(grad(37, 1), 0.50026067941832E-03, 1E-6);
  EXPECT_NEAR(grad(37, 2), 0.34510037138103E-03, 1E-6);
  EXPECT_NEAR(grad(38, 0), 0.26496879921132E-04, 1E-6);
  EXPECT_NEAR(grad(38, 1), -0.37417878566112E-04, 1E-6);
  EXPECT_NEAR(grad(38, 2), -0.71129234809557E-05, 1E-6);
  EXPECT_NEAR(grad(39, 0), 0.83448184914835E-04, 1E-6);
  EXPECT_NEAR(grad(39, 1), -0.22008189141890E-03, 1E-6);
  EXPECT_NEAR(grad(39, 2), -0.24529271311915E-04, 1E-6);
  EXPECT_NEAR(grad(40, 0), 0.26708938471020E-03, 1E-6);
  EXPECT_NEAR(grad(40, 1), -0.27072292230860E-03, 1E-6);
  EXPECT_NEAR(grad(40, 2), 0.16428038551436E-04, 1E-6);
  EXPECT_NEAR(grad(41, 0), 0.30046832671835E-03, 1E-6);
  EXPECT_NEAR(grad(41, 1), 0.39242810715911E-03, 1E-6);
  EXPECT_NEAR(grad(41, 2), 0.84546460198716E-04, 1E-6);
  EXPECT_NEAR(grad(42, 0), -0.58858292946469E-03, 1E-6);
  EXPECT_NEAR(grad(42, 1), 0.27871026830993E-03, 1E-6);
  EXPECT_NEAR(grad(42, 2), 0.22211280085369E-03, 1E-6);
  EXPECT_NEAR(grad(43, 0), 0.56706507988895E-04, 1E-6);
  EXPECT_NEAR(grad(43, 1), 0.43582759454250E-03, 1E-6);
  EXPECT_NEAR(grad(43, 2), 0.27801828665525E-03, 1E-6);
  EXPECT_NEAR(grad(44, 0), 0.43743953900482E-03, 1E-6);
  EXPECT_NEAR(grad(44, 1), -0.99307370325658E-04, 1E-6);
  EXPECT_NEAR(grad(44, 2), 0.13036978975765E-03, 1E-6);
  EXPECT_NEAR(grad(45, 0), 0.14668060877444E-03, 1E-6);
  EXPECT_NEAR(grad(45, 1), 0.43042983525089E-03, 1E-6);
  EXPECT_NEAR(grad(45, 2), 0.10951569870183E-03, 1E-6);
  EXPECT_NEAR(grad(46, 0), -0.48384638486832E-03, 1E-6);
  EXPECT_NEAR(grad(46, 1), 0.54059541107157E-03, 1E-6);
  EXPECT_NEAR(grad(46, 2), 0.23151468012265E-03, 1E-6);
  EXPECT_NEAR(grad(47, 0), -0.31897390893856E-05, 1E-6);
  EXPECT_NEAR(grad(47, 1), -0.96230492759290E-04, 1E-6);
  EXPECT_NEAR(grad(47, 2), -0.10438134801079E-04, 1E-6);
  EXPECT_NEAR(grad(48, 0), -0.58833135719325E-03, 1E-6);
  EXPECT_NEAR(grad(48, 1), 0.76487599382896E-03, 1E-6);
  EXPECT_NEAR(grad(48, 2), 0.29974231333890E-03, 1E-6);
  EXPECT_NEAR(grad(49, 0), -0.69589688100785E-04, 1E-6);
  EXPECT_NEAR(grad(49, 1), 0.61987946322691E-04, 1E-6);
  EXPECT_NEAR(grad(49, 2), 0.37754781642052E-03, 1E-6);
  EXPECT_NEAR(grad(50, 0), -0.50042057497647E-03, 1E-6);
  EXPECT_NEAR(grad(50, 1), 0.23842094888645E-03, 1E-6);
  EXPECT_NEAR(grad(50, 2), 0.32693081986390E-03, 1E-6);
  EXPECT_NEAR(grad(51, 0), -0.13512127477070E-03, 1E-6);
  EXPECT_NEAR(grad(51, 1), 0.37754409652096E-04, 1E-6);
  EXPECT_NEAR(grad(51, 2), -0.10840049573994E-03, 1E-6);
  EXPECT_NEAR(grad(52, 0), 0.60883664576034E-03, 1E-6);
  EXPECT_NEAR(grad(52, 1), 0.18256601815862E-03, 1E-6);
  EXPECT_NEAR(grad(52, 2), 0.56297019445765E-04, 1E-6);
  EXPECT_NEAR(grad(53, 0), -0.54843471478183E-04, 1E-6);
  EXPECT_NEAR(grad(53, 1), -0.95069207015725E-04, 1E-6);
  EXPECT_NEAR(grad(53, 2), -0.47762032804462E-04, 1E-6);
  EXPECT_NEAR(grad(54, 0), 0.10067059644245E-02, 1E-6);
  EXPECT_NEAR(grad(54, 1), 0.18882544427757E-03, 1E-6);
  EXPECT_NEAR(grad(54, 2), 0.27680643620143E-03, 1E-6);
  EXPECT_NEAR(grad(55, 0), -0.16120892073712E-03, 1E-6);
  EXPECT_NEAR(grad(55, 1), 0.95888630818158E-04, 1E-6);
  EXPECT_NEAR(grad(55, 2), -0.42407326790284E-04, 1E-6);
  EXPECT_NEAR(grad(56, 0), 0.16679665659103E-03, 1E-6);
  EXPECT_NEAR(grad(56, 1), -0.85986842380693E-04, 1E-6);
  EXPECT_NEAR(grad(56, 2), -0.52242841541106E-04, 1E-6);
  EXPECT_NEAR(grad(57, 0), 0.10502431145237E-03, 1E-6);
  EXPECT_NEAR(grad(57, 1), -0.10873702095048E-03, 1E-6);
  EXPECT_NEAR(grad(57, 2), 0.11022059232849E-03, 1E-6);
  EXPECT_NEAR(grad(58, 0), 0.33834914873263E-03, 1E-6);
  EXPECT_NEAR(grad(58, 1), 0.25337847397295E-04, 1E-6);
  EXPECT_NEAR(grad(58, 2), -0.13951960351985E-03, 1E-6);
  EXPECT_NEAR(grad(59, 0), 0.10146858805040E-03, 1E-6);
  EXPECT_NEAR(grad(59, 1), -0.14292427521129E-04, 1E-6);
  EXPECT_NEAR(grad(59, 2), -0.76567065642898E-04, 1E-6);
  EXPECT_NEAR(grad(60, 0), 0.29118452525138E-04, 1E-6);
  EXPECT_NEAR(grad(60, 1), -0.22764740554297E-04, 1E-6);
  EXPECT_NEAR(grad(60, 2), -0.67860490987149E-04, 1E-6);
  EXPECT_NEAR(grad(61, 0), 0.18480864065011E-03, 1E-6);
  EXPECT_NEAR(grad(61, 1), -0.49275524140741E-04, 1E-6);
  EXPECT_NEAR(grad(61, 2), 0.31287187806142E-03, 1E-6);
  EXPECT_NEAR(grad(62, 0), -0.15234608856936E-04, 1E-6);
  EXPECT_NEAR(grad(62, 1), -0.44444645457712E-04, 1E-6);
  EXPECT_NEAR(grad(62, 2), 0.56313725109488E-04, 1E-6);
  EXPECT_NEAR(grad(63, 0), 0.36071929609932E-04, 1E-6);
  EXPECT_NEAR(grad(63, 1), -0.49416717891306E-04, 1E-6);
  EXPECT_NEAR(grad(63, 2), 0.11071147927040E-03, 1E-6);
  EXPECT_NEAR(grad(64, 0), -0.44719829563178E-05, 1E-6);
  EXPECT_NEAR(grad(64, 1), -0.13475590901461E-03, 1E-6);
  EXPECT_NEAR(grad(64, 2), -0.13023147443576E-03, 1E-6);
  EXPECT_NEAR(grad(65, 0), 0.74810712495518E-04, 1E-6);
  EXPECT_NEAR(grad(65, 1), -0.39211158329040E-04, 1E-6);
  EXPECT_NEAR(grad(65, 2), 0.22470527764129E-04, 1E-6);
  EXPECT_NEAR(grad(66, 0), -0.65821174667100E-04, 1E-6);
  EXPECT_NEAR(grad(66, 1), -0.16453999238782E-03, 1E-6);
  EXPECT_NEAR(grad(66, 2), 0.51530018769344E-04, 1E-6);
  EXPECT_NEAR(grad(67, 0), 0.51136003864819E-04, 1E-6);
  EXPECT_NEAR(grad(67, 1), 0.47256791711954E-04, 1E-6);
  EXPECT_NEAR(grad(67, 2), 0.12159305117493E-03, 1E-6);
  EXPECT_NEAR(grad(68, 0), 0.10822290357402E-03, 1E-6);
  EXPECT_NEAR(grad(68, 1), 0.76178857280038E-04, 1E-6);
  EXPECT_NEAR(grad(68, 2), -0.51330112084068E-04, 1E-6);
  EXPECT_NEAR(grad(69, 0), -0.67384873979368E-04, 1E-6);
  EXPECT_NEAR(grad(69, 1), 0.85634365497184E-04, 1E-6);
  EXPECT_NEAR(grad(69, 2), -0.27601519632602E-04, 1E-6);
  EXPECT_NEAR(grad(70, 0), 0.67297416968531E-03, 1E-6);
  EXPECT_NEAR(grad(70, 1), 0.22075467453647E-04, 1E-6);
  EXPECT_NEAR(grad(70, 2), 0.45379317113982E-03, 1E-6);
  EXPECT_NEAR(grad(71, 0), 0.60070738505901E-03, 1E-6);
  EXPECT_NEAR(grad(71, 1), 0.93614425119700E-04, 1E-6);
  EXPECT_NEAR(grad(71, 2), 0.29173525815529E-03, 1E-6);
  EXPECT_NEAR(grad(72, 0), 0.57388094216686E-04, 1E-6);
  EXPECT_NEAR(grad(72, 1), 0.24125830247974E-03, 1E-6);
  EXPECT_NEAR(grad(72, 2), 0.31489109430659E-03, 1E-6);
  EXPECT_NEAR(grad(73, 0), 0.66130579167566E-04, 1E-6);
  EXPECT_NEAR(grad(73, 1), 0.32145123257677E-03, 1E-6);
  EXPECT_NEAR(grad(73, 2), -0.21299802482742E-03, 1E-6);
  EXPECT_NEAR(grad(74, 0), 0.44755352732545E-03, 1E-6);
  EXPECT_NEAR(grad(74, 1), 0.12274240646114E-03, 1E-6);
  EXPECT_NEAR(grad(74, 2), -0.72619890743671E-04, 1E-6);
  EXPECT_NEAR(grad(75, 0), 0.44134258508854E-03, 1E-6);
  EXPECT_NEAR(grad(75, 1), 0.69810089287954E-04, 1E-6);
  EXPECT_NEAR(grad(75, 2), -0.15444003602011E-03, 1E-6);
  EXPECT_NEAR(grad(76, 0), 0.50794953953515E-04, 1E-6);
  EXPECT_NEAR(grad(76, 1), 0.14954410394429E-03, 1E-6);
  EXPECT_NEAR(grad(76, 2), 0.34718765594306E-03, 1E-6);
  EXPECT_NEAR(grad(77, 0), -0.35052846111756E-03, 1E-6);
  EXPECT_NEAR(grad(77, 1), 0.23690045086682E-03, 1E-6);
  EXPECT_NEAR(grad(77, 2), 0.24177273605188E-03, 1E-6);
  EXPECT_NEAR(grad(78, 0), -0.36711835134182E-03, 1E-6);
  EXPECT_NEAR(grad(78, 1), 0.34368266795802E-03, 1E-6);
  EXPECT_NEAR(grad(78, 2), 0.25139286706542E-03, 1E-6);
  EXPECT_NEAR(grad(79, 0), 0.51571849619032E-05, 1E-6);
  EXPECT_NEAR(grad(79, 1), 0.10873828545275E-03, 1E-6);
  EXPECT_NEAR(grad(79, 2), -0.63294402071644E-04, 1E-6);
  EXPECT_NEAR(grad(80, 0), -0.72217689699263E-04, 1E-6);
  EXPECT_NEAR(grad(80, 1), 0.10434438056613E-03, 1E-6);
  EXPECT_NEAR(grad(80, 2), 0.13031658651699E-03, 1E-6);
  EXPECT_NEAR(grad(81, 0), 0.97642354659945E-04, 1E-6);
  EXPECT_NEAR(grad(81, 1), 0.40438756251961E-04, 1E-6);
  EXPECT_NEAR(grad(81, 2), 0.29140205409474E-04, 1E-6);
  EXPECT_NEAR(grad(82, 0), -0.39200281549390E-03, 1E-6);
  EXPECT_NEAR(grad(82, 1), 0.52028337519622E-03, 1E-6);
  EXPECT_NEAR(grad(82, 2), -0.14003851372667E-03, 1E-6);
  EXPECT_NEAR(grad(83, 0), -0.40389620473747E-03, 1E-6);
  EXPECT_NEAR(grad(83, 1), 0.48176442571294E-03, 1E-6);
  EXPECT_NEAR(grad(83, 2), -0.43555156834671E-04, 1E-6);
  EXPECT_NEAR(grad(84, 0), 0.22015800174970E-03, 1E-6);
  EXPECT_NEAR(grad(84, 1), 0.49579514240542E-03, 1E-6);
  EXPECT_NEAR(grad(84, 2), 0.38723138307901E-04, 1E-6);
  EXPECT_NEAR(grad(85, 0), -0.41260096010419E-03, 1E-6);
  EXPECT_NEAR(grad(85, 1), -0.44591984297989E-04, 1E-6);
  EXPECT_NEAR(grad(85, 2), 0.46894670785220E-03, 1E-6);
  EXPECT_NEAR(grad(86, 0), -0.13470950532923E-03, 1E-6);
  EXPECT_NEAR(grad(86, 1), 0.34634413198267E-04, 1E-6);
  EXPECT_NEAR(grad(86, 2), 0.25743677376417E-04, 1E-6);
  EXPECT_NEAR(grad(87, 0), -0.84230600023913E-04, 1E-6);
  EXPECT_NEAR(grad(87, 1), -0.10641359062809E-03, 1E-6);
  EXPECT_NEAR(grad(87, 2), -0.66771174423974E-04, 1E-6);
  EXPECT_NEAR(grad(88, 0), -0.30776480888594E-03, 1E-6);
  EXPECT_NEAR(grad(88, 1), 0.85943994591127E-04, 1E-6);
  EXPECT_NEAR(grad(88, 2), 0.21326987237353E-03, 1E-6);
  EXPECT_NEAR(grad(89, 0), -0.10313170019803E-03, 1E-6);
  EXPECT_NEAR(grad(89, 1), -0.51259678857995E-05, 1E-6);
  EXPECT_NEAR(grad(89, 2), 0.74580805643801E-04, 1E-6);
  EXPECT_NEAR(grad(90, 0), -0.54172057048804E-04, 1E-6);
  EXPECT_NEAR(grad(90, 1), -0.35775842903383E-04, 1E-6);
  EXPECT_NEAR(grad(90, 2), 0.55319033649167E-04, 1E-6);
  EXPECT_NEAR(grad(91, 0), -0.45070026896180E-03, 1E-6);
  EXPECT_NEAR(grad(91, 1), 0.53885543656225E-03, 1E-6);
  EXPECT_NEAR(grad(91, 2), -0.37584689850402E-04, 1E-6);
  EXPECT_NEAR(grad(92, 0), -0.82631849638403E-04, 1E-6);
  EXPECT_NEAR(grad(92, 1), 0.82732510070554E-04, 1E-6);
  EXPECT_NEAR(grad(92, 2), -0.21839513125740E-05, 1E-6);
  EXPECT_NEAR(grad(93, 0), -0.89786851591123E-04, 1E-6);
  EXPECT_NEAR(grad(93, 1), 0.80813651863136E-04, 1E-6);
  EXPECT_NEAR(grad(93, 2), -0.68562172635354E-04, 1E-6);
  EXPECT_NEAR(grad(94, 0), 0.68781552639488E-03, 1E-6);
  EXPECT_NEAR(grad(94, 1), 0.86888375263331E-05, 1E-6);
  EXPECT_NEAR(grad(94, 2), -0.26447104074035E-03, 1E-6);
  EXPECT_NEAR(grad(95, 0), 0.45038436516607E-03, 1E-6);
  EXPECT_NEAR(grad(95, 1), 0.23851119441798E-03, 1E-6);
  EXPECT_NEAR(grad(95, 2), -0.12945730120987E-02, 1E-6);
  EXPECT_NEAR(grad(96, 0), 0.17754842185878E-04, 1E-6);
  EXPECT_NEAR(grad(96, 1), 0.29416639711523E-03, 1E-6);
  EXPECT_NEAR(grad(96, 2), -0.16110265629177E-02, 1E-6);
  EXPECT_NEAR(grad(97, 0), 0.22545324749821E-04, 1E-6);
  EXPECT_NEAR(grad(97, 1), -0.10707840976185E-03, 1E-6);
  EXPECT_NEAR(grad(97, 2), -0.12588619104665E-02, 1E-6);
  EXPECT_NEAR(grad(98, 0), 0.13357884094037E-03, 1E-6);
  EXPECT_NEAR(grad(98, 1), -0.86384364722966E-04, 1E-6);
  EXPECT_NEAR(grad(98, 2), -0.56435621164196E-03, 1E-6);
  EXPECT_NEAR(grad(99, 0), 0.34254467274320E-04, 1E-6);
  EXPECT_NEAR(grad(99, 1), 0.31499585837147E-03, 1E-6);
  EXPECT_NEAR(grad(99, 2), -0.11565146816106E-02, 1E-6);
  EXPECT_NEAR(grad(100, 0), 0.47138789630553E-04, 1E-6);
  EXPECT_NEAR(grad(100, 1), 0.21068582739932E-03, 1E-6);
  EXPECT_NEAR(grad(100, 2), -0.70182914437991E-03, 1E-6);
  EXPECT_NEAR(grad(101, 0), -0.48013535564164E-03, 1E-6);
  EXPECT_NEAR(grad(101, 1), 0.21506198743377E-03, 1E-6);
  EXPECT_NEAR(grad(101, 2), -0.14399730307689E-02, 1E-6);
  EXPECT_NEAR(grad(102, 0), 0.24615961411572E-03, 1E-6);
  EXPECT_NEAR(grad(102, 1), -0.36675005319376E-03, 1E-6);
  EXPECT_NEAR(grad(102, 2), -0.10964148957526E-02, 1E-6);
  EXPECT_NEAR(grad(103, 0), -0.51839241201035E-03, 1E-6);
  EXPECT_NEAR(grad(103, 1), 0.11651526428440E-03, 1E-6);
  EXPECT_NEAR(grad(103, 2), -0.12225089197687E-02, 1E-6);
  EXPECT_NEAR(grad(104, 0), 0.18221540271894E-03, 1E-6);
  EXPECT_NEAR(grad(104, 1), 0.34605389817808E-03, 1E-6);
  EXPECT_NEAR(grad(104, 2), -0.69304585970063E-03, 1E-6);
  EXPECT_NEAR(grad(105, 0), -0.11607514844081E-03, 1E-6);
  EXPECT_NEAR(grad(105, 1), 0.15597637073673E-03, 1E-6);
  EXPECT_NEAR(grad(105, 2), -0.54508449544137E-03, 1E-6);
  EXPECT_NEAR(grad(106, 0), -0.15183279265987E-03, 1E-6);
  EXPECT_NEAR(grad(106, 1), -0.24258823429891E-04, 1E-6);
  EXPECT_NEAR(grad(106, 2), -0.49306848074821E-03, 1E-6);
  EXPECT_NEAR(grad(107, 0), 0.17332992366834E-03, 1E-6);
  EXPECT_NEAR(grad(107, 1), -0.27852369391819E-03, 1E-6);
  EXPECT_NEAR(grad(107, 2), -0.59342359501126E-03, 1E-6);
  EXPECT_NEAR(grad(108, 0), 0.21946567983676E-03, 1E-6);
  EXPECT_NEAR(grad(108, 1), 0.15628395056612E-04, 1E-6);
  EXPECT_NEAR(grad(108, 2), -0.49867419140733E-03, 1E-6);
  EXPECT_NEAR(grad(109, 0), -0.77168461227281E-04, 1E-6);
  EXPECT_NEAR(grad(109, 1), -0.53417524868352E-04, 1E-6);
  EXPECT_NEAR(grad(109, 2), -0.11996635991775E-03, 1E-6);
  EXPECT_NEAR(grad(110, 0), -0.67202447662548E-04, 1E-6);
  EXPECT_NEAR(grad(110, 1), -0.83858691399475E-04, 1E-6);
  EXPECT_NEAR(grad(110, 2), -0.16068519102951E-03, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) gradient correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3GradBP86JacobsenGhost) {
  auto geometry = *systemControllerTwo->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3,
                                                                               std::make_shared<Geometry>(geometry),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 0), -0.23939648114107E-03, 1E-6);
  EXPECT_NEAR(grad(0, 1), -0.29326368463196E-03, 1E-6);
  EXPECT_NEAR(grad(0, 2), 0.16643267153564E-03, 1E-6);
  EXPECT_NEAR(grad(1, 0), 0.28480927202525E-03, 1E-6);
  EXPECT_NEAR(grad(1, 1), -0.66020038784406E-03, 1E-6);
  EXPECT_NEAR(grad(1, 2), 0.38634083205742E-04, 1E-6);
  EXPECT_NEAR(grad(2, 0), -0.55805381243386E-04, 1E-6);
  EXPECT_NEAR(grad(2, 1), -0.14769958781335E-03, 1E-6);
  EXPECT_NEAR(grad(2, 2), 0.43079324928880E-03, 1E-6);
  EXPECT_NEAR(grad(3, 0), -0.10159865831420E-03, 1E-6);
  EXPECT_NEAR(grad(3, 1), -0.24115292587994E-03, 1E-6);
  EXPECT_NEAR(grad(3, 2), 0.39270485334179E-03, 1E-6);
  EXPECT_NEAR(grad(4, 0), 0.64648503102299E-03, 1E-6);
  EXPECT_NEAR(grad(4, 1), 0.11014200526491E-03, 1E-6);
  EXPECT_NEAR(grad(4, 2), 0.77315756624325E-03, 1E-6);
  EXPECT_NEAR(grad(5, 0), 0.15600742018761E-03, 1E-6);
  EXPECT_NEAR(grad(5, 1), -0.83701252002179E-04, 1E-6);
  EXPECT_NEAR(grad(5, 2), -0.10320055086246E-03, 1E-6);
  EXPECT_NEAR(grad(6, 0), -0.26369629280981E-03, 1E-6);
  EXPECT_NEAR(grad(6, 1), 0.43056225680949E-03, 1E-6);
  EXPECT_NEAR(grad(6, 2), 0.25897355496868E-03, 1E-6);
  EXPECT_NEAR(grad(7, 0), 0.23727668340770E-03, 1E-6);
  EXPECT_NEAR(grad(7, 1), -0.71058464700238E-03, 1E-6);
  EXPECT_NEAR(grad(7, 2), -0.22975462470042E-03, 1E-6);
  EXPECT_NEAR(grad(8, 0), 0.51233163024399E-04, 1E-6);
  EXPECT_NEAR(grad(8, 1), -0.20086917537838E-03, 1E-6);
  EXPECT_NEAR(grad(8, 2), -0.10567080208689E-03, 1E-6);
  EXPECT_NEAR(grad(9, 0), -0.25335844672692E-03, 1E-6);
  EXPECT_NEAR(grad(9, 1), -0.51486210410658E-03, 1E-6);
  EXPECT_NEAR(grad(9, 2), 0.11292206538659E-03, 1E-6);
  EXPECT_NEAR(grad(10, 0), -0.16750617929842E-03, 1E-6);
  EXPECT_NEAR(grad(10, 1), -0.13617186678462E-03, 1E-6);
  EXPECT_NEAR(grad(10, 2), 0.65781593998751E-04, 1E-6);
  EXPECT_NEAR(grad(11, 0), -0.64870113385215E-04, 1E-6);
  EXPECT_NEAR(grad(11, 1), -0.18880389933317E-03, 1E-6);
  EXPECT_NEAR(grad(11, 2), 0.12669866738940E-03, 1E-6);
  EXPECT_NEAR(grad(12, 0), -0.41101498125640E-04, 1E-6);
  EXPECT_NEAR(grad(12, 1), -0.91089636774285E-04, 1E-6);
  EXPECT_NEAR(grad(12, 2), -0.41780393924603E-05, 1E-6);
  EXPECT_NEAR(grad(13, 0), -0.21937055610206E-03, 1E-6);
  EXPECT_NEAR(grad(13, 1), -0.21022567376362E-03, 1E-6);
  EXPECT_NEAR(grad(13, 2), -0.53355798116362E-04, 1E-6);
  EXPECT_NEAR(grad(14, 0), 0.12096030488783E-03, 1E-6);
  EXPECT_NEAR(grad(14, 1), -0.33481498535850E-03, 1E-6);
  EXPECT_NEAR(grad(14, 2), -0.20951168017895E-04, 1E-6);
  EXPECT_NEAR(grad(15, 0), 0.29819435863964E-05, 1E-6);
  EXPECT_NEAR(grad(15, 1), -0.88856288253517E-04, 1E-6);
  EXPECT_NEAR(grad(15, 2), -0.87371658347790E-05, 1E-6);
  EXPECT_NEAR(grad(16, 0), 0.16722652455443E-03, 1E-6);
  EXPECT_NEAR(grad(16, 1), -0.41481076468878E-03, 1E-6);
  EXPECT_NEAR(grad(16, 2), -0.69862591390206E-04, 1E-6);
  EXPECT_NEAR(grad(17, 0), 0.10772916885965E-03, 1E-6);
  EXPECT_NEAR(grad(17, 1), -0.16632306472675E-03, 1E-6);
  EXPECT_NEAR(grad(17, 2), -0.27498132784997E-04, 1E-6);
  EXPECT_NEAR(grad(18, 0), 0.67637761977715E-03, 1E-6);
  EXPECT_NEAR(grad(18, 1), -0.21667916154471E-03, 1E-6);
  EXPECT_NEAR(grad(18, 2), 0.58562657565709E-04, 1E-6);
  EXPECT_NEAR(grad(19, 0), -0.74736649925546E-03, 1E-6);
  EXPECT_NEAR(grad(19, 1), -0.12040683859880E-04, 1E-6);
  EXPECT_NEAR(grad(19, 2), 0.72367216058471E-03, 1E-6);
  EXPECT_NEAR(grad(20, 0), 0.33510778156219E-03, 1E-6);
  EXPECT_NEAR(grad(20, 1), -0.12848194374529E-02, 1E-6);
  EXPECT_NEAR(grad(20, 2), 0.65719560841665E-03, 1E-6);
  EXPECT_NEAR(grad(21, 0), -0.51700318671666E-03, 1E-6);
  EXPECT_NEAR(grad(21, 1), -0.12995799212020E-02, 1E-6);
  EXPECT_NEAR(grad(21, 2), 0.44596927969442E-03, 1E-6);
  EXPECT_NEAR(grad(22, 0), 0.23229288681336E-03, 1E-6);
  EXPECT_NEAR(grad(22, 1), -0.46685255032151E-04, 1E-6);
  EXPECT_NEAR(grad(22, 2), 0.15731997236800E-03, 1E-6);
  EXPECT_NEAR(grad(23, 0), -0.12877068374006E-03, 1E-6);
  EXPECT_NEAR(grad(23, 1), -0.71406503902311E-04, 1E-6);
  EXPECT_NEAR(grad(23, 2), 0.45524638085392E-03, 1E-6);
  EXPECT_NEAR(grad(24, 0), -0.79414902533742E-03, 1E-6);
  EXPECT_NEAR(grad(24, 1), 0.13031633455296E-04, 1E-6);
  EXPECT_NEAR(grad(24, 2), 0.50053818798825E-03, 1E-6);
  EXPECT_NEAR(grad(25, 0), -0.70347565420089E-04, 1E-6);
  EXPECT_NEAR(grad(25, 1), 0.13634425309191E-03, 1E-6);
  EXPECT_NEAR(grad(25, 2), 0.99452415420100E-03, 1E-6);
  EXPECT_NEAR(grad(26, 0), -0.11104642475317E-03, 1E-6);
  EXPECT_NEAR(grad(26, 1), -0.46423862927161E-04, 1E-6);
  EXPECT_NEAR(grad(26, 2), 0.58987772591003E-03, 1E-6);
  EXPECT_NEAR(grad(27, 0), -0.23657197709089E-03, 1E-6);
  EXPECT_NEAR(grad(27, 1), -0.20104935787111E-03, 1E-6);
  EXPECT_NEAR(grad(27, 2), 0.47002007086454E-03, 1E-6);
  EXPECT_NEAR(grad(28, 0), -0.24782242124356E-04, 1E-6);
  EXPECT_NEAR(grad(28, 1), -0.19150823349044E-03, 1E-6);
  EXPECT_NEAR(grad(28, 2), 0.97555104070209E-03, 1E-6);
  EXPECT_NEAR(grad(29, 0), 0.12853758383436E-03, 1E-6);
  EXPECT_NEAR(grad(29, 1), 0.19637817812301E-03, 1E-6);
  EXPECT_NEAR(grad(29, 2), -0.14691683814867E-03, 1E-6);
  EXPECT_NEAR(grad(30, 0), 0.13602991186037E-03, 1E-6);
  EXPECT_NEAR(grad(30, 1), -0.20967995868233E-03, 1E-6);
  EXPECT_NEAR(grad(30, 2), -0.56941257610763E-04, 1E-6);
  EXPECT_NEAR(grad(31, 0), 0.57719700391863E-03, 1E-6);
  EXPECT_NEAR(grad(31, 1), -0.36817858906916E-03, 1E-6);
  EXPECT_NEAR(grad(31, 2), -0.54971762179040E-05, 1E-6);
  EXPECT_NEAR(grad(32, 0), 0.46534525505210E-04, 1E-6);
  EXPECT_NEAR(grad(32, 1), -0.92278239751090E-05, 1E-6);
  EXPECT_NEAR(grad(32, 2), -0.94781232424354E-04, 1E-6);
  EXPECT_NEAR(grad(33, 0), -0.17825669853133E-03, 1E-6);
  EXPECT_NEAR(grad(33, 1), -0.74080518499352E-04, 1E-6);
  EXPECT_NEAR(grad(33, 2), -0.15910188089558E-03, 1E-6);
  EXPECT_NEAR(grad(34, 0), -0.49119782519775E-03, 1E-6);
  EXPECT_NEAR(grad(34, 1), -0.54382238997521E-04, 1E-6);
  EXPECT_NEAR(grad(34, 2), 0.37033534453372E-03, 1E-6);
  EXPECT_NEAR(grad(35, 0), -0.28727571461093E-03, 1E-6);
  EXPECT_NEAR(grad(35, 1), 0.82022566482949E-04, 1E-6);
  EXPECT_NEAR(grad(35, 2), 0.63619548621831E-03, 1E-6);
  EXPECT_NEAR(grad(36, 0), -0.39231469999683E-03, 1E-6);
  EXPECT_NEAR(grad(36, 1), -0.14956241246950E-03, 1E-6);
  EXPECT_NEAR(grad(36, 2), 0.57041515808173E-03, 1E-6);
  EXPECT_NEAR(grad(37, 0), -0.70825954266066E-04, 1E-6);
  EXPECT_NEAR(grad(37, 1), 0.50026067941832E-03, 1E-6);
  EXPECT_NEAR(grad(37, 2), 0.34510037138103E-03, 1E-6);
  EXPECT_NEAR(grad(38, 0), 0.26496879921132E-04, 1E-6);
  EXPECT_NEAR(grad(38, 1), -0.37417878566112E-04, 1E-6);
  EXPECT_NEAR(grad(38, 2), -0.71129234809557E-05, 1E-6);
  EXPECT_NEAR(grad(39, 0), 0.83448184914835E-04, 1E-6);
  EXPECT_NEAR(grad(39, 1), -0.22008189141890E-03, 1E-6);
  EXPECT_NEAR(grad(39, 2), -0.24529271311915E-04, 1E-6);
  EXPECT_NEAR(grad(40, 0), 0.26708938471020E-03, 1E-6);
  EXPECT_NEAR(grad(40, 1), -0.27072292230860E-03, 1E-6);
  EXPECT_NEAR(grad(40, 2), 0.16428038551436E-04, 1E-6);
  EXPECT_NEAR(grad(41, 0), 0.30046832671835E-03, 1E-6);
  EXPECT_NEAR(grad(41, 1), 0.39242810715911E-03, 1E-6);
  EXPECT_NEAR(grad(41, 2), 0.84546460198716E-04, 1E-6);
  EXPECT_NEAR(grad(42, 0), -0.58858292946469E-03, 1E-6);
  EXPECT_NEAR(grad(42, 1), 0.27871026830993E-03, 1E-6);
  EXPECT_NEAR(grad(42, 2), 0.22211280085369E-03, 1E-6);
  EXPECT_NEAR(grad(43, 0), 0.56706507988895E-04, 1E-6);
  EXPECT_NEAR(grad(43, 1), 0.43582759454250E-03, 1E-6);
  EXPECT_NEAR(grad(43, 2), 0.27801828665525E-03, 1E-6);
  EXPECT_NEAR(grad(44, 0), 0.43743953900482E-03, 1E-6);
  EXPECT_NEAR(grad(44, 1), -0.99307370325658E-04, 1E-6);
  EXPECT_NEAR(grad(44, 2), 0.13036978975765E-03, 1E-6);
  EXPECT_NEAR(grad(45, 0), 0.14668060877444E-03, 1E-6);
  EXPECT_NEAR(grad(45, 1), 0.43042983525089E-03, 1E-6);
  EXPECT_NEAR(grad(45, 2), 0.10951569870183E-03, 1E-6);
  EXPECT_NEAR(grad(46, 0), -0.48384638486832E-03, 1E-6);
  EXPECT_NEAR(grad(46, 1), 0.54059541107157E-03, 1E-6);
  EXPECT_NEAR(grad(46, 2), 0.23151468012265E-03, 1E-6);
  EXPECT_NEAR(grad(47, 0), -0.31897390893856E-05, 1E-6);
  EXPECT_NEAR(grad(47, 1), -0.96230492759290E-04, 1E-6);
  EXPECT_NEAR(grad(47, 2), -0.10438134801079E-04, 1E-6);
  EXPECT_NEAR(grad(48, 0), -0.58833135719325E-03, 1E-6);
  EXPECT_NEAR(grad(48, 1), 0.76487599382896E-03, 1E-6);
  EXPECT_NEAR(grad(48, 2), 0.29974231333890E-03, 1E-6);
  EXPECT_NEAR(grad(49, 0), -0.69589688100785E-04, 1E-6);
  EXPECT_NEAR(grad(49, 1), 0.61987946322691E-04, 1E-6);
  EXPECT_NEAR(grad(49, 2), 0.37754781642052E-03, 1E-6);
  EXPECT_NEAR(grad(50, 0), -0.50042057497647E-03, 1E-6);
  EXPECT_NEAR(grad(50, 1), 0.23842094888645E-03, 1E-6);
  EXPECT_NEAR(grad(50, 2), 0.32693081986390E-03, 1E-6);
  EXPECT_NEAR(grad(51, 0), -0.13512127477070E-03, 1E-6);
  EXPECT_NEAR(grad(51, 1), 0.37754409652096E-04, 1E-6);
  EXPECT_NEAR(grad(51, 2), -0.10840049573994E-03, 1E-6);
  EXPECT_NEAR(grad(52, 0), 0.60883664576034E-03, 1E-6);
  EXPECT_NEAR(grad(52, 1), 0.18256601815862E-03, 1E-6);
  EXPECT_NEAR(grad(52, 2), 0.56297019445765E-04, 1E-6);
  EXPECT_NEAR(grad(53, 0), -0.54843471478183E-04, 1E-6);
  EXPECT_NEAR(grad(53, 1), -0.95069207015725E-04, 1E-6);
  EXPECT_NEAR(grad(53, 2), -0.47762032804462E-04, 1E-6);
  EXPECT_NEAR(grad(54, 0), 0.10067059644245E-02, 1E-6);
  EXPECT_NEAR(grad(54, 1), 0.18882544427757E-03, 1E-6);
  EXPECT_NEAR(grad(54, 2), 0.27680643620143E-03, 1E-6);
  EXPECT_NEAR(grad(55, 0), -0.16120892073712E-03, 1E-6);
  EXPECT_NEAR(grad(55, 1), 0.95888630818158E-04, 1E-6);
  EXPECT_NEAR(grad(55, 2), -0.42407326790284E-04, 1E-6);
  EXPECT_NEAR(grad(56, 0), 0.16679665659103E-03, 1E-6);
  EXPECT_NEAR(grad(56, 1), -0.85986842380693E-04, 1E-6);
  EXPECT_NEAR(grad(56, 2), -0.52242841541106E-04, 1E-6);
  EXPECT_NEAR(grad(57, 0), 0.10502431145237E-03, 1E-6);
  EXPECT_NEAR(grad(57, 1), -0.10873702095048E-03, 1E-6);
  EXPECT_NEAR(grad(57, 2), 0.11022059232849E-03, 1E-6);
  EXPECT_NEAR(grad(58, 0), 0.33834914873263E-03, 1E-6);
  EXPECT_NEAR(grad(58, 1), 0.25337847397295E-04, 1E-6);
  EXPECT_NEAR(grad(58, 2), -0.13951960351985E-03, 1E-6);
  EXPECT_NEAR(grad(59, 0), 0.10146858805040E-03, 1E-6);
  EXPECT_NEAR(grad(59, 1), -0.14292427521129E-04, 1E-6);
  EXPECT_NEAR(grad(59, 2), -0.76567065642898E-04, 1E-6);
  EXPECT_NEAR(grad(60, 0), 0.29118452525138E-04, 1E-6);
  EXPECT_NEAR(grad(60, 1), -0.22764740554297E-04, 1E-6);
  EXPECT_NEAR(grad(60, 2), -0.67860490987149E-04, 1E-6);
  EXPECT_NEAR(grad(61, 0), 0.18480864065011E-03, 1E-6);
  EXPECT_NEAR(grad(61, 1), -0.49275524140741E-04, 1E-6);
  EXPECT_NEAR(grad(61, 2), 0.31287187806142E-03, 1E-6);
  EXPECT_NEAR(grad(62, 0), -0.15234608856936E-04, 1E-6);
  EXPECT_NEAR(grad(62, 1), -0.44444645457712E-04, 1E-6);
  EXPECT_NEAR(grad(62, 2), 0.56313725109488E-04, 1E-6);
  EXPECT_NEAR(grad(63, 0), 0.36071929609932E-04, 1E-6);
  EXPECT_NEAR(grad(63, 1), -0.49416717891306E-04, 1E-6);
  EXPECT_NEAR(grad(63, 2), 0.11071147927040E-03, 1E-6);
  EXPECT_NEAR(grad(64, 0), -0.44719829563178E-05, 1E-6);
  EXPECT_NEAR(grad(64, 1), -0.13475590901461E-03, 1E-6);
  EXPECT_NEAR(grad(64, 2), -0.13023147443576E-03, 1E-6);
  EXPECT_NEAR(grad(65, 0), 0.74810712495518E-04, 1E-6);
  EXPECT_NEAR(grad(65, 1), -0.39211158329040E-04, 1E-6);
  EXPECT_NEAR(grad(65, 2), 0.22470527764129E-04, 1E-6);
  EXPECT_NEAR(grad(66, 0), -0.65821174667100E-04, 1E-6);
  EXPECT_NEAR(grad(66, 1), -0.16453999238782E-03, 1E-6);
  EXPECT_NEAR(grad(66, 2), 0.51530018769344E-04, 1E-6);
  EXPECT_NEAR(grad(67, 0), 0.51136003864819E-04, 1E-6);
  EXPECT_NEAR(grad(67, 1), 0.47256791711954E-04, 1E-6);
  EXPECT_NEAR(grad(67, 2), 0.12159305117493E-03, 1E-6);
  EXPECT_NEAR(grad(68, 0), 0.10822290357402E-03, 1E-6);
  EXPECT_NEAR(grad(68, 1), 0.76178857280038E-04, 1E-6);
  EXPECT_NEAR(grad(68, 2), -0.51330112084068E-04, 1E-6);
  EXPECT_NEAR(grad(69, 0), -0.67384873979368E-04, 1E-6);
  EXPECT_NEAR(grad(69, 1), 0.85634365497184E-04, 1E-6);
  EXPECT_NEAR(grad(69, 2), -0.27601519632602E-04, 1E-6);
  EXPECT_NEAR(grad(70, 0), 0.67297416968531E-03, 1E-6);
  EXPECT_NEAR(grad(70, 1), 0.22075467453647E-04, 1E-6);
  EXPECT_NEAR(grad(70, 2), 0.45379317113982E-03, 1E-6);
  EXPECT_NEAR(grad(71, 0), 0.60070738505901E-03, 1E-6);
  EXPECT_NEAR(grad(71, 1), 0.93614425119700E-04, 1E-6);
  EXPECT_NEAR(grad(71, 2), 0.29173525815529E-03, 1E-6);
  EXPECT_NEAR(grad(72, 0), 0.57388094216686E-04, 1E-6);
  EXPECT_NEAR(grad(72, 1), 0.24125830247974E-03, 1E-6);
  EXPECT_NEAR(grad(72, 2), 0.31489109430659E-03, 1E-6);
  EXPECT_NEAR(grad(73, 0), 0.66130579167566E-04, 1E-6);
  EXPECT_NEAR(grad(73, 1), 0.32145123257677E-03, 1E-6);
  EXPECT_NEAR(grad(73, 2), -0.21299802482742E-03, 1E-6);
  EXPECT_NEAR(grad(74, 0), 0.44755352732545E-03, 1E-6);
  EXPECT_NEAR(grad(74, 1), 0.12274240646114E-03, 1E-6);
  EXPECT_NEAR(grad(74, 2), -0.72619890743671E-04, 1E-6);
  EXPECT_NEAR(grad(75, 0), 0.44134258508854E-03, 1E-6);
  EXPECT_NEAR(grad(75, 1), 0.69810089287954E-04, 1E-6);
  EXPECT_NEAR(grad(75, 2), -0.15444003602011E-03, 1E-6);
  EXPECT_NEAR(grad(76, 0), 0.50794953953515E-04, 1E-6);
  EXPECT_NEAR(grad(76, 1), 0.14954410394429E-03, 1E-6);
  EXPECT_NEAR(grad(76, 2), 0.34718765594306E-03, 1E-6);
  EXPECT_NEAR(grad(77, 0), -0.35052846111756E-03, 1E-6);
  EXPECT_NEAR(grad(77, 1), 0.23690045086682E-03, 1E-6);
  EXPECT_NEAR(grad(77, 2), 0.24177273605188E-03, 1E-6);
  EXPECT_NEAR(grad(78, 0), -0.36711835134182E-03, 1E-6);
  EXPECT_NEAR(grad(78, 1), 0.34368266795802E-03, 1E-6);
  EXPECT_NEAR(grad(78, 2), 0.25139286706542E-03, 1E-6);
  EXPECT_NEAR(grad(79, 0), 0.51571849619032E-05, 1E-6);
  EXPECT_NEAR(grad(79, 1), 0.10873828545275E-03, 1E-6);
  EXPECT_NEAR(grad(79, 2), -0.63294402071644E-04, 1E-6);
  EXPECT_NEAR(grad(80, 0), -0.72217689699263E-04, 1E-6);
  EXPECT_NEAR(grad(80, 1), 0.10434438056613E-03, 1E-6);
  EXPECT_NEAR(grad(80, 2), 0.13031658651699E-03, 1E-6);
  EXPECT_NEAR(grad(81, 0), 0.97642354659945E-04, 1E-6);
  EXPECT_NEAR(grad(81, 1), 0.40438756251961E-04, 1E-6);
  EXPECT_NEAR(grad(81, 2), 0.29140205409474E-04, 1E-6);
  EXPECT_NEAR(grad(82, 0), -0.39200281549390E-03, 1E-6);
  EXPECT_NEAR(grad(82, 1), 0.52028337519622E-03, 1E-6);
  EXPECT_NEAR(grad(82, 2), -0.14003851372667E-03, 1E-6);
  EXPECT_NEAR(grad(83, 0), -0.40389620473747E-03, 1E-6);
  EXPECT_NEAR(grad(83, 1), 0.48176442571294E-03, 1E-6);
  EXPECT_NEAR(grad(83, 2), -0.43555156834671E-04, 1E-6);
  EXPECT_NEAR(grad(84, 0), 0.22015800174970E-03, 1E-6);
  EXPECT_NEAR(grad(84, 1), 0.49579514240542E-03, 1E-6);
  EXPECT_NEAR(grad(84, 2), 0.38723138307901E-04, 1E-6);
  EXPECT_NEAR(grad(85, 0), -0.41260096010419E-03, 1E-6);
  EXPECT_NEAR(grad(85, 1), -0.44591984297989E-04, 1E-6);
  EXPECT_NEAR(grad(85, 2), 0.46894670785220E-03, 1E-6);
  EXPECT_NEAR(grad(86, 0), -0.13470950532923E-03, 1E-6);
  EXPECT_NEAR(grad(86, 1), 0.34634413198267E-04, 1E-6);
  EXPECT_NEAR(grad(86, 2), 0.25743677376417E-04, 1E-6);
  EXPECT_NEAR(grad(87, 0), -0.84230600023913E-04, 1E-6);
  EXPECT_NEAR(grad(87, 1), -0.10641359062809E-03, 1E-6);
  EXPECT_NEAR(grad(87, 2), -0.66771174423974E-04, 1E-6);
  EXPECT_NEAR(grad(88, 0), -0.30776480888594E-03, 1E-6);
  EXPECT_NEAR(grad(88, 1), 0.85943994591127E-04, 1E-6);
  EXPECT_NEAR(grad(88, 2), 0.21326987237353E-03, 1E-6);
  EXPECT_NEAR(grad(89, 0), -0.10313170019803E-03, 1E-6);
  EXPECT_NEAR(grad(89, 1), -0.51259678857995E-05, 1E-6);
  EXPECT_NEAR(grad(89, 2), 0.74580805643801E-04, 1E-6);
  EXPECT_NEAR(grad(90, 0), -0.54172057048804E-04, 1E-6);
  EXPECT_NEAR(grad(90, 1), -0.35775842903383E-04, 1E-6);
  EXPECT_NEAR(grad(90, 2), 0.55319033649167E-04, 1E-6);
  EXPECT_NEAR(grad(91, 0), -0.45070026896180E-03, 1E-6);
  EXPECT_NEAR(grad(91, 1), 0.53885543656225E-03, 1E-6);
  EXPECT_NEAR(grad(91, 2), -0.37584689850402E-04, 1E-6);
  EXPECT_NEAR(grad(92, 0), -0.82631849638403E-04, 1E-6);
  EXPECT_NEAR(grad(92, 1), 0.82732510070554E-04, 1E-6);
  EXPECT_NEAR(grad(92, 2), -0.21839513125740E-05, 1E-6);
  EXPECT_NEAR(grad(93, 0), -0.89786851591123E-04, 1E-6);
  EXPECT_NEAR(grad(93, 1), 0.80813651863136E-04, 1E-6);
  EXPECT_NEAR(grad(93, 2), -0.68562172635354E-04, 1E-6);
  EXPECT_NEAR(grad(94, 0), 0.68781552639488E-03, 1E-6);
  EXPECT_NEAR(grad(94, 1), 0.86888375263331E-05, 1E-6);
  EXPECT_NEAR(grad(94, 2), -0.26447104074035E-03, 1E-6);
  EXPECT_NEAR(grad(95, 0), 0.45038436516607E-03, 1E-6);
  EXPECT_NEAR(grad(95, 1), 0.23851119441798E-03, 1E-6);
  EXPECT_NEAR(grad(95, 2), -0.12945730120987E-02, 1E-6);
  EXPECT_NEAR(grad(96, 0), 0.17754842185878E-04, 1E-6);
  EXPECT_NEAR(grad(96, 1), 0.29416639711523E-03, 1E-6);
  EXPECT_NEAR(grad(96, 2), -0.16110265629177E-02, 1E-6);
  EXPECT_NEAR(grad(97, 0), 0.22545324749821E-04, 1E-6);
  EXPECT_NEAR(grad(97, 1), -0.10707840976185E-03, 1E-6);
  EXPECT_NEAR(grad(97, 2), -0.12588619104665E-02, 1E-6);
  EXPECT_NEAR(grad(98, 0), 0.13357884094037E-03, 1E-6);
  EXPECT_NEAR(grad(98, 1), -0.86384364722966E-04, 1E-6);
  EXPECT_NEAR(grad(98, 2), -0.56435621164196E-03, 1E-6);
  EXPECT_NEAR(grad(99, 0), 0.34254467274320E-04, 1E-6);
  EXPECT_NEAR(grad(99, 1), 0.31499585837147E-03, 1E-6);
  EXPECT_NEAR(grad(99, 2), -0.11565146816106E-02, 1E-6);
  EXPECT_NEAR(grad(100, 0), 0.47138789630553E-04, 1E-6);
  EXPECT_NEAR(grad(100, 1), 0.21068582739932E-03, 1E-6);
  EXPECT_NEAR(grad(100, 2), -0.70182914437991E-03, 1E-6);
  EXPECT_NEAR(grad(101, 0), -0.48013535564164E-03, 1E-6);
  EXPECT_NEAR(grad(101, 1), 0.21506198743377E-03, 1E-6);
  EXPECT_NEAR(grad(101, 2), -0.14399730307689E-02, 1E-6);
  EXPECT_NEAR(grad(102, 0), 0.24615961411572E-03, 1E-6);
  EXPECT_NEAR(grad(102, 1), -0.36675005319376E-03, 1E-6);
  EXPECT_NEAR(grad(102, 2), -0.10964148957526E-02, 1E-6);
  EXPECT_NEAR(grad(103, 0), -0.51839241201035E-03, 1E-6);
  EXPECT_NEAR(grad(103, 1), 0.11651526428440E-03, 1E-6);
  EXPECT_NEAR(grad(103, 2), -0.12225089197687E-02, 1E-6);
  EXPECT_NEAR(grad(104, 0), 0.18221540271894E-03, 1E-6);
  EXPECT_NEAR(grad(104, 1), 0.34605389817808E-03, 1E-6);
  EXPECT_NEAR(grad(104, 2), -0.69304585970063E-03, 1E-6);
  EXPECT_NEAR(grad(105, 0), -0.11607514844081E-03, 1E-6);
  EXPECT_NEAR(grad(105, 1), 0.15597637073673E-03, 1E-6);
  EXPECT_NEAR(grad(105, 2), -0.54508449544137E-03, 1E-6);
  EXPECT_NEAR(grad(106, 0), -0.15183279265987E-03, 1E-6);
  EXPECT_NEAR(grad(106, 1), -0.24258823429891E-04, 1E-6);
  EXPECT_NEAR(grad(106, 2), -0.49306848074821E-03, 1E-6);
  EXPECT_NEAR(grad(107, 0), 0.17332992366834E-03, 1E-6);
  EXPECT_NEAR(grad(107, 1), -0.27852369391819E-03, 1E-6);
  EXPECT_NEAR(grad(107, 2), -0.59342359501126E-03, 1E-6);
  EXPECT_NEAR(grad(108, 0), 0.21946567983676E-03, 1E-6);
  EXPECT_NEAR(grad(108, 1), 0.15628395056612E-04, 1E-6);
  EXPECT_NEAR(grad(108, 2), -0.49867419140733E-03, 1E-6);
  EXPECT_NEAR(grad(109, 0), -0.77168461227281E-04, 1E-6);
  EXPECT_NEAR(grad(109, 1), -0.53417524868352E-04, 1E-6);
  EXPECT_NEAR(grad(109, 2), -0.11996635991775E-03, 1E-6);
  EXPECT_NEAR(grad(110, 0), -0.67202447662548E-04, 1E-6);
  EXPECT_NEAR(grad(110, 1), -0.83858691399475E-04, 1E-6);
  EXPECT_NEAR(grad(110, 2), -0.16068519102951E-03, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC gradient correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3ABCGradBP86CO) {
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3,
                                                                               systemControllerOne->getGeometry(),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 2), 0.27878862999524E-5, 1E-8);
  EXPECT_NEAR(grad(1, 2), -0.27878862999524E-5, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC gradient correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3ABCGradBP86Jacobsen) {
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3ABC,
                                                                               systemControllerTwo->getGeometry(),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 0), -0.23717801640170E-03, 1E-6);
  EXPECT_NEAR(grad(0, 1), -0.28987042072207E-03, 1E-6);
  EXPECT_NEAR(grad(0, 2), 0.18600295030725E-03, 1E-6);
  EXPECT_NEAR(grad(1, 0), 0.27714377137876E-03, 1E-6);
  EXPECT_NEAR(grad(1, 1), -0.64601348693930E-03, 1E-6);
  EXPECT_NEAR(grad(1, 2), 0.22098175770654E-04, 1E-6);
  EXPECT_NEAR(grad(2, 0), -0.63165373673185E-04, 1E-6);
  EXPECT_NEAR(grad(2, 1), -0.13687906670926E-03, 1E-6);
  EXPECT_NEAR(grad(2, 2), 0.43815172780821E-03, 1E-6);
  EXPECT_NEAR(grad(3, 0), -0.12332675677029E-03, 1E-6);
  EXPECT_NEAR(grad(3, 1), -0.22738995146344E-03, 1E-6);
  EXPECT_NEAR(grad(3, 2), 0.38611094473551E-03, 1E-6);
  EXPECT_NEAR(grad(4, 0), 0.64575814860422E-03, 1E-6);
  EXPECT_NEAR(grad(4, 1), 0.12162856186441E-03, 1E-6);
  EXPECT_NEAR(grad(4, 2), 0.78228828377076E-03, 1E-6);
  EXPECT_NEAR(grad(5, 0), 0.16818541189921E-03, 1E-6);
  EXPECT_NEAR(grad(5, 1), -0.92153755939268E-04, 1E-6);
  EXPECT_NEAR(grad(5, 2), -0.97564433713353E-04, 1E-6);
  EXPECT_NEAR(grad(6, 0), -0.28044777249553E-03, 1E-6);
  EXPECT_NEAR(grad(6, 1), 0.44618601379171E-03, 1E-6);
  EXPECT_NEAR(grad(6, 2), 0.25092593198033E-03, 1E-6);
  EXPECT_NEAR(grad(7, 0), 0.24501087991183E-03, 1E-6);
  EXPECT_NEAR(grad(7, 1), -0.69880968688684E-03, 1E-6);
  EXPECT_NEAR(grad(7, 2), -0.24194202305893E-03, 1E-6);
  EXPECT_NEAR(grad(8, 0), 0.59200909543117E-04, 1E-6);
  EXPECT_NEAR(grad(8, 1), -0.20702836543792E-03, 1E-6);
  EXPECT_NEAR(grad(8, 2), -0.97174388308769E-04, 1E-6);
  EXPECT_NEAR(grad(9, 0), -0.25530222548126E-03, 1E-6);
  EXPECT_NEAR(grad(9, 1), -0.50582250253523E-03, 1E-6);
  EXPECT_NEAR(grad(9, 2), 0.11157145257822E-03, 1E-6);
  EXPECT_NEAR(grad(10, 0), -0.16496930249006E-03, 1E-6);
  EXPECT_NEAR(grad(10, 1), -0.12726313793206E-03, 1E-6);
  EXPECT_NEAR(grad(10, 2), 0.68257972726503E-04, 1E-6);
  EXPECT_NEAR(grad(11, 0), -0.76811744317813E-04, 1E-6);
  EXPECT_NEAR(grad(11, 1), -0.17754812228552E-03, 1E-6);
  EXPECT_NEAR(grad(11, 2), 0.11493032735417E-03, 1E-6);
  EXPECT_NEAR(grad(12, 0), -0.42612090021756E-04, 1E-6);
  EXPECT_NEAR(grad(12, 1), -0.81441515719382E-04, 1E-6);
  EXPECT_NEAR(grad(12, 2), -0.89513699952739E-05, 1E-6);
  EXPECT_NEAR(grad(13, 0), -0.21827099739717E-03, 1E-6);
  EXPECT_NEAR(grad(13, 1), -0.19835049585459E-03, 1E-6);
  EXPECT_NEAR(grad(13, 2), -0.52359199059803E-04, 1E-6);
  EXPECT_NEAR(grad(14, 0), 0.12315402180042E-03, 1E-6);
  EXPECT_NEAR(grad(14, 1), -0.32629077803331E-03, 1E-6);
  EXPECT_NEAR(grad(14, 2), -0.20763908944206E-04, 1E-6);
  EXPECT_NEAR(grad(15, 0), 0.46211874845382E-05, 1E-6);
  EXPECT_NEAR(grad(15, 1), -0.85718368010777E-04, 1E-6);
  EXPECT_NEAR(grad(15, 2), -0.51443606585082E-05, 1E-6);
  EXPECT_NEAR(grad(16, 0), 0.16315334500310E-03, 1E-6);
  EXPECT_NEAR(grad(16, 1), -0.41100812582001E-03, 1E-6);
  EXPECT_NEAR(grad(16, 2), -0.63022662401746E-04, 1E-6);
  EXPECT_NEAR(grad(17, 0), 0.10289249631693E-03, 1E-6);
  EXPECT_NEAR(grad(17, 1), -0.15747144737284E-03, 1E-6);
  EXPECT_NEAR(grad(17, 2), -0.34958185107614E-04, 1E-6);
  EXPECT_NEAR(grad(18, 0), 0.67992378875347E-03, 1E-6);
  EXPECT_NEAR(grad(18, 1), -0.20268825107959E-03, 1E-6);
  EXPECT_NEAR(grad(18, 2), 0.56995805361152E-04, 1E-6);
  EXPECT_NEAR(grad(19, 0), -0.75714494690225E-03, 1E-6);
  EXPECT_NEAR(grad(19, 1), -0.10222744775054E-05, 1E-6);
  EXPECT_NEAR(grad(19, 2), 0.73443376080339E-03, 1E-6);
  EXPECT_NEAR(grad(20, 0), 0.32739575379033E-03, 1E-6);
  EXPECT_NEAR(grad(20, 1), -0.12758497091448E-02, 1E-6);
  EXPECT_NEAR(grad(20, 2), 0.64175507129259E-03, 1E-6);
  EXPECT_NEAR(grad(21, 0), -0.52441799028590E-03, 1E-6);
  EXPECT_NEAR(grad(21, 1), -0.12853515735787E-02, 1E-6);
  EXPECT_NEAR(grad(21, 2), 0.44642213284109E-03, 1E-6);
  EXPECT_NEAR(grad(22, 0), 0.21670526023561E-03, 1E-6);
  EXPECT_NEAR(grad(22, 1), -0.28999225726610E-04, 1E-6);
  EXPECT_NEAR(grad(22, 2), 0.14676694325673E-03, 1E-6);
  EXPECT_NEAR(grad(23, 0), -0.12333261329891E-03, 1E-6);
  EXPECT_NEAR(grad(23, 1), -0.71931707092751E-04, 1E-6);
  EXPECT_NEAR(grad(23, 2), 0.44887653605974E-03, 1E-6);
  EXPECT_NEAR(grad(24, 0), -0.79657380099408E-03, 1E-6);
  EXPECT_NEAR(grad(24, 1), 0.67776107497750E-05, 1E-6);
  EXPECT_NEAR(grad(24, 2), 0.50121240439673E-03, 1E-6);
  EXPECT_NEAR(grad(25, 0), -0.65040106103026E-04, 1E-6);
  EXPECT_NEAR(grad(25, 1), 0.12650321767912E-03, 1E-6);
  EXPECT_NEAR(grad(25, 2), 0.99985078988799E-03, 1E-6);
  EXPECT_NEAR(grad(26, 0), -0.10432373823870E-03, 1E-6);
  EXPECT_NEAR(grad(26, 1), -0.55318799754943E-04, 1E-6);
  EXPECT_NEAR(grad(26, 2), 0.58744450419345E-03, 1E-6);
  EXPECT_NEAR(grad(27, 0), -0.24562200600174E-03, 1E-6);
  EXPECT_NEAR(grad(27, 1), -0.20140492475905E-03, 1E-6);
  EXPECT_NEAR(grad(27, 2), 0.47103035784947E-03, 1E-6);
  EXPECT_NEAR(grad(28, 0), -0.28488943965613E-04, 1E-6);
  EXPECT_NEAR(grad(28, 1), -0.20256025714324E-03, 1E-6);
  EXPECT_NEAR(grad(28, 2), 0.98037957845948E-03, 1E-6);
  EXPECT_NEAR(grad(29, 0), 0.11147795588742E-03, 1E-6);
  EXPECT_NEAR(grad(29, 1), 0.20795196148474E-03, 1E-6);
  EXPECT_NEAR(grad(29, 2), -0.15301660100037E-03, 1E-6);
  EXPECT_NEAR(grad(30, 0), 0.13451771804650E-03, 1E-6);
  EXPECT_NEAR(grad(30, 1), -0.21996564946691E-03, 1E-6);
  EXPECT_NEAR(grad(30, 2), -0.61757580869626E-04, 1E-6);
  EXPECT_NEAR(grad(31, 0), 0.57414127146254E-03, 1E-6);
  EXPECT_NEAR(grad(31, 1), -0.37171570269451E-03, 1E-6);
  EXPECT_NEAR(grad(31, 2), -0.11805456544775E-04, 1E-6);
  EXPECT_NEAR(grad(32, 0), 0.32530921494129E-04, 1E-6);
  EXPECT_NEAR(grad(32, 1), 0.61876148944902E-05, 1E-6);
  EXPECT_NEAR(grad(32, 2), -0.99019785973023E-04, 1E-6);
  EXPECT_NEAR(grad(33, 0), -0.17139895446928E-03, 1E-6);
  EXPECT_NEAR(grad(33, 1), -0.78032692433836E-04, 1E-6);
  EXPECT_NEAR(grad(33, 2), -0.15732827256839E-03, 1E-6);
  EXPECT_NEAR(grad(34, 0), -0.49341198602911E-03, 1E-6);
  EXPECT_NEAR(grad(34, 1), -0.61743444855760E-04, 1E-6);
  EXPECT_NEAR(grad(34, 2), 0.36648513219175E-03, 1E-6);
  EXPECT_NEAR(grad(35, 0), -0.26547658033976E-03, 1E-6);
  EXPECT_NEAR(grad(35, 1), 0.71697420789409E-04, 1E-6);
  EXPECT_NEAR(grad(35, 2), 0.64264297133084E-03, 1E-6);
  EXPECT_NEAR(grad(36, 0), -0.38465784590321E-03, 1E-6);
  EXPECT_NEAR(grad(36, 1), -0.15947120297169E-03, 1E-6);
  EXPECT_NEAR(grad(36, 2), 0.57096258492138E-03, 1E-6);
  EXPECT_NEAR(grad(37, 0), -0.54932625298560E-04, 1E-6);
  EXPECT_NEAR(grad(37, 1), 0.49249136281485E-03, 1E-6);
  EXPECT_NEAR(grad(37, 2), 0.34401752509547E-03, 1E-6);
  EXPECT_NEAR(grad(38, 0), 0.96791122176194E-05, 1E-6);
  EXPECT_NEAR(grad(38, 1), -0.25769272209720E-04, 1E-6);
  EXPECT_NEAR(grad(38, 2), -0.88815305236382E-05, 1E-6);
  EXPECT_NEAR(grad(39, 0), 0.77696687296514E-04, 1E-6);
  EXPECT_NEAR(grad(39, 1), -0.22277274905422E-03, 1E-6);
  EXPECT_NEAR(grad(39, 2), -0.24664568728692E-04, 1E-6);
  EXPECT_NEAR(grad(40, 0), 0.26974825652855E-03, 1E-6);
  EXPECT_NEAR(grad(40, 1), -0.28063870364769E-03, 1E-6);
  EXPECT_NEAR(grad(40, 2), 0.16343519077373E-04, 1E-6);
  EXPECT_NEAR(grad(41, 0), 0.28086755687923E-03, 1E-6);
  EXPECT_NEAR(grad(41, 1), 0.39351507801869E-03, 1E-6);
  EXPECT_NEAR(grad(41, 2), 0.76915390321853E-04, 1E-6);
  EXPECT_NEAR(grad(42, 0), -0.59962854513644E-03, 1E-6);
  EXPECT_NEAR(grad(42, 1), 0.28353335400055E-03, 1E-6);
  EXPECT_NEAR(grad(42, 2), 0.21538495303442E-03, 1E-6);
  EXPECT_NEAR(grad(43, 0), 0.62601743123119E-04, 1E-6);
  EXPECT_NEAR(grad(43, 1), 0.43957115723508E-03, 1E-6);
  EXPECT_NEAR(grad(43, 2), 0.26961378468819E-03, 1E-6);
  EXPECT_NEAR(grad(44, 0), 0.45470139826651E-03, 1E-6);
  EXPECT_NEAR(grad(44, 1), -0.10374796869838E-03, 1E-6);
  EXPECT_NEAR(grad(44, 2), 0.13857950034783E-03, 1E-6);
  EXPECT_NEAR(grad(45, 0), 0.14621405858143E-03, 1E-6);
  EXPECT_NEAR(grad(45, 1), 0.44049117441358E-03, 1E-6);
  EXPECT_NEAR(grad(45, 2), 0.11350077007003E-03, 1E-6);
  EXPECT_NEAR(grad(46, 0), -0.49403232081828E-03, 1E-6);
  EXPECT_NEAR(grad(46, 1), 0.53267080440346E-03, 1E-6);
  EXPECT_NEAR(grad(46, 2), 0.23135400719157E-03, 1E-6);
  EXPECT_NEAR(grad(47, 0), -0.12507698657128E-04, 1E-6);
  EXPECT_NEAR(grad(47, 1), -0.91913868078356E-04, 1E-6);
  EXPECT_NEAR(grad(47, 2), -0.15923484875254E-04, 1E-6);
  EXPECT_NEAR(grad(48, 0), -0.59920015451261E-03, 1E-6);
  EXPECT_NEAR(grad(48, 1), 0.76347893970001E-03, 1E-6);
  EXPECT_NEAR(grad(48, 2), 0.28459865567657E-03, 1E-6);
  EXPECT_NEAR(grad(49, 0), -0.54350476610919E-04, 1E-6);
  EXPECT_NEAR(grad(49, 1), 0.44517323767844E-04, 1E-6);
  EXPECT_NEAR(grad(49, 2), 0.37718403727515E-03, 1E-6);
  EXPECT_NEAR(grad(50, 0), -0.47841592399576E-03, 1E-6);
  EXPECT_NEAR(grad(50, 1), 0.23655916460177E-03, 1E-6);
  EXPECT_NEAR(grad(50, 2), 0.32110226559319E-03, 1E-6);
  EXPECT_NEAR(grad(51, 0), -0.11920946235248E-03, 1E-6);
  EXPECT_NEAR(grad(51, 1), 0.29878907176467E-04, 1E-6);
  EXPECT_NEAR(grad(51, 2), -0.96850931160190E-04, 1E-6);
  EXPECT_NEAR(grad(52, 0), 0.61330818744228E-03, 1E-6);
  EXPECT_NEAR(grad(52, 1), 0.17975767686066E-03, 1E-6);
  EXPECT_NEAR(grad(52, 2), 0.52938062027837E-04, 1E-6);
  EXPECT_NEAR(grad(53, 0), -0.44299839413577E-04, 1E-6);
  EXPECT_NEAR(grad(53, 1), -0.88609368332008E-04, 1E-6);
  EXPECT_NEAR(grad(53, 2), -0.41585309230315E-04, 1E-6);
  EXPECT_NEAR(grad(54, 0), 0.10150763871525E-02, 1E-6);
  EXPECT_NEAR(grad(54, 1), 0.18756779178569E-03, 1E-6);
  EXPECT_NEAR(grad(54, 2), 0.28681211423454E-03, 1E-6);
  EXPECT_NEAR(grad(55, 0), -0.17655992229250E-03, 1E-6);
  EXPECT_NEAR(grad(55, 1), 0.91438478076860E-04, 1E-6);
  EXPECT_NEAR(grad(55, 2), -0.48379822979946E-04, 1E-6);
  EXPECT_NEAR(grad(56, 0), 0.14763856630641E-03, 1E-6);
  EXPECT_NEAR(grad(56, 1), -0.80163757160533E-04, 1E-6);
  EXPECT_NEAR(grad(56, 2), -0.46900095211109E-04, 1E-6);
  EXPECT_NEAR(grad(57, 0), 0.95910512539158E-04, 1E-6);
  EXPECT_NEAR(grad(57, 1), -0.99726589290551E-04, 1E-6);
  EXPECT_NEAR(grad(57, 2), 0.96102877321659E-04, 1E-6);
  EXPECT_NEAR(grad(58, 0), 0.33520552531891E-03, 1E-6);
  EXPECT_NEAR(grad(58, 1), 0.29913547620154E-04, 1E-6);
  EXPECT_NEAR(grad(58, 2), -0.14124748328383E-03, 1E-6);
  EXPECT_NEAR(grad(59, 0), 0.10205051868591E-03, 1E-6);
  EXPECT_NEAR(grad(59, 1), -0.15363718400346E-04, 1E-6);
  EXPECT_NEAR(grad(59, 2), -0.73215738010412E-04, 1E-6);
  EXPECT_NEAR(grad(60, 0), 0.31206295463529E-04, 1E-6);
  EXPECT_NEAR(grad(60, 1), -0.22748496722329E-04, 1E-6);
  EXPECT_NEAR(grad(60, 2), -0.72337845666214E-04, 1E-6);
  EXPECT_NEAR(grad(61, 0), 0.18042733573728E-03, 1E-6);
  EXPECT_NEAR(grad(61, 1), -0.45779690392349E-04, 1E-6);
  EXPECT_NEAR(grad(61, 2), 0.31427285253189E-03, 1E-6);
  EXPECT_NEAR(grad(62, 0), -0.17226594573020E-04, 1E-6);
  EXPECT_NEAR(grad(62, 1), -0.45591584626299E-04, 1E-6);
  EXPECT_NEAR(grad(62, 2), 0.60752993934197E-04, 1E-6);
  EXPECT_NEAR(grad(63, 0), 0.38205654772853E-04, 1E-6);
  EXPECT_NEAR(grad(63, 1), -0.50010959823938E-04, 1E-6);
  EXPECT_NEAR(grad(63, 2), 0.10916974670623E-03, 1E-6);
  EXPECT_NEAR(grad(64, 0), -0.33748835115731E-05, 1E-6);
  EXPECT_NEAR(grad(64, 1), -0.13484587103423E-03, 1E-6);
  EXPECT_NEAR(grad(64, 2), -0.12487578225330E-03, 1E-6);
  EXPECT_NEAR(grad(65, 0), 0.79639567517741E-04, 1E-6);
  EXPECT_NEAR(grad(65, 1), -0.38404958353020E-04, 1E-6);
  EXPECT_NEAR(grad(65, 2), 0.25194631743046E-04, 1E-6);
  EXPECT_NEAR(grad(66, 0), -0.62313269442493E-04, 1E-6);
  EXPECT_NEAR(grad(66, 1), -0.16301505345603E-03, 1E-6);
  EXPECT_NEAR(grad(66, 2), 0.49945375463408E-04, 1E-6);
  EXPECT_NEAR(grad(67, 0), 0.46385309389828E-04, 1E-6);
  EXPECT_NEAR(grad(67, 1), 0.44811757577484E-04, 1E-6);
  EXPECT_NEAR(grad(67, 2), 0.11337491216874E-03, 1E-6);
  EXPECT_NEAR(grad(68, 0), 0.10020531487048E-03, 1E-6);
  EXPECT_NEAR(grad(68, 1), 0.72029779165086E-04, 1E-6);
  EXPECT_NEAR(grad(68, 2), -0.47447693661869E-04, 1E-6);
  EXPECT_NEAR(grad(69, 0), -0.71885733670794E-04, 1E-6);
  EXPECT_NEAR(grad(69, 1), 0.79197431349367E-04, 1E-6);
  EXPECT_NEAR(grad(69, 2), -0.29985666081673E-04, 1E-6);
  EXPECT_NEAR(grad(70, 0), 0.66657750351279E-03, 1E-6);
  EXPECT_NEAR(grad(70, 1), 0.80348112451879E-05, 1E-6);
  EXPECT_NEAR(grad(70, 2), 0.46509998380525E-03, 1E-6);
  EXPECT_NEAR(grad(71, 0), 0.59686857378042E-03, 1E-6);
  EXPECT_NEAR(grad(71, 1), 0.87156782141222E-04, 1E-6);
  EXPECT_NEAR(grad(71, 2), 0.28834722803016E-03, 1E-6);
  EXPECT_NEAR(grad(72, 0), 0.43457789773417E-04, 1E-6);
  EXPECT_NEAR(grad(72, 1), 0.22939963012866E-03, 1E-6);
  EXPECT_NEAR(grad(72, 2), 0.31347335281305E-03, 1E-6);
  EXPECT_NEAR(grad(73, 0), 0.53666331383611E-04, 1E-6);
  EXPECT_NEAR(grad(73, 1), 0.31076405924674E-03, 1E-6);
  EXPECT_NEAR(grad(73, 2), -0.21786753826681E-03, 1E-6);
  EXPECT_NEAR(grad(74, 0), 0.43999351148239E-03, 1E-6);
  EXPECT_NEAR(grad(74, 1), 0.11509127644925E-03, 1E-6);
  EXPECT_NEAR(grad(74, 2), -0.71534603996133E-04, 1E-6);
  EXPECT_NEAR(grad(75, 0), 0.43179822590886E-03, 1E-6);
  EXPECT_NEAR(grad(75, 1), 0.57500761210677E-04, 1E-6);
  EXPECT_NEAR(grad(75, 2), -0.16973737775986E-03, 1E-6);
  EXPECT_NEAR(grad(76, 0), 0.54492436234436E-04, 1E-6);
  EXPECT_NEAR(grad(76, 1), 0.13388535195938E-03, 1E-6);
  EXPECT_NEAR(grad(76, 2), 0.33878565294576E-03, 1E-6);
  EXPECT_NEAR(grad(77, 0), -0.34718879492883E-03, 1E-6);
  EXPECT_NEAR(grad(77, 1), 0.23069474191604E-03, 1E-6);
  EXPECT_NEAR(grad(77, 2), 0.23294024710741E-03, 1E-6);
  EXPECT_NEAR(grad(78, 0), -0.37065763710889E-03, 1E-6);
  EXPECT_NEAR(grad(78, 1), 0.32799122541849E-03, 1E-6);
  EXPECT_NEAR(grad(78, 2), 0.25589878604821E-03, 1E-6);
  EXPECT_NEAR(grad(79, 0), 0.83248427105235E-05, 1E-6);
  EXPECT_NEAR(grad(79, 1), 0.97885464966770E-04, 1E-6);
  EXPECT_NEAR(grad(79, 2), -0.64106444935496E-04, 1E-6);
  EXPECT_NEAR(grad(80, 0), -0.68012435609132E-04, 1E-6);
  EXPECT_NEAR(grad(80, 1), 0.97248355558922E-04, 1E-6);
  EXPECT_NEAR(grad(80, 2), 0.12085109824906E-03, 1E-6);
  EXPECT_NEAR(grad(81, 0), 0.10067005866307E-03, 1E-6);
  EXPECT_NEAR(grad(81, 1), 0.31985509955836E-04, 1E-6);
  EXPECT_NEAR(grad(81, 2), 0.24142441382561E-04, 1E-6);
  EXPECT_NEAR(grad(82, 0), -0.39592147417943E-03, 1E-6);
  EXPECT_NEAR(grad(82, 1), 0.51487090788829E-03, 1E-6);
  EXPECT_NEAR(grad(82, 2), -0.15697821907618E-03, 1E-6);
  EXPECT_NEAR(grad(83, 0), -0.40574768483625E-03, 1E-6);
  EXPECT_NEAR(grad(83, 1), 0.47384114765033E-03, 1E-6);
  EXPECT_NEAR(grad(83, 2), -0.46267285246846E-04, 1E-6);
  EXPECT_NEAR(grad(84, 0), 0.22385332802663E-03, 1E-6);
  EXPECT_NEAR(grad(84, 1), 0.47977995413444E-03, 1E-6);
  EXPECT_NEAR(grad(84, 2), 0.33560043398743E-04, 1E-6);
  EXPECT_NEAR(grad(85, 0), -0.41113959268474E-03, 1E-6);
  EXPECT_NEAR(grad(85, 1), -0.45939920943882E-04, 1E-6);
  EXPECT_NEAR(grad(85, 2), 0.46548649050262E-03, 1E-6);
  EXPECT_NEAR(grad(86, 0), -0.12969186143636E-03, 1E-6);
  EXPECT_NEAR(grad(86, 1), 0.37120288362057E-04, 1E-6);
  EXPECT_NEAR(grad(86, 2), 0.20356776238145E-04, 1E-6);
  EXPECT_NEAR(grad(87, 0), -0.76099512050028E-04, 1E-6);
  EXPECT_NEAR(grad(87, 1), -0.10620104574951E-03, 1E-6);
  EXPECT_NEAR(grad(87, 2), -0.70375515396661E-04, 1E-6);
  EXPECT_NEAR(grad(88, 0), -0.29991445279378E-03, 1E-6);
  EXPECT_NEAR(grad(88, 1), 0.84537364662448E-04, 1E-6);
  EXPECT_NEAR(grad(88, 2), 0.21364953246601E-03, 1E-6);
  EXPECT_NEAR(grad(89, 0), -0.10196937659892E-03, 1E-6);
  EXPECT_NEAR(grad(89, 1), -0.56884565505157E-05, 1E-6);
  EXPECT_NEAR(grad(89, 2), 0.68544773544331E-04, 1E-6);
  EXPECT_NEAR(grad(90, 0), -0.51941504567816E-04, 1E-6);
  EXPECT_NEAR(grad(90, 1), -0.38852948426937E-04, 1E-6);
  EXPECT_NEAR(grad(90, 2), 0.55508831642056E-04, 1E-6);
  EXPECT_NEAR(grad(91, 0), -0.44243048904823E-03, 1E-6);
  EXPECT_NEAR(grad(91, 1), 0.54027733043081E-03, 1E-6);
  EXPECT_NEAR(grad(91, 2), -0.38226138152842E-04, 1E-6);
  EXPECT_NEAR(grad(92, 0), -0.73391765815008E-04, 1E-6);
  EXPECT_NEAR(grad(92, 1), 0.76554117053676E-04, 1E-6);
  EXPECT_NEAR(grad(92, 2), -0.78239905949963E-05, 1E-6);
  EXPECT_NEAR(grad(93, 0), -0.85672561524956E-04, 1E-6);
  EXPECT_NEAR(grad(93, 1), 0.72740622131106E-04, 1E-6);
  EXPECT_NEAR(grad(93, 2), -0.69150038840129E-04, 1E-6);
  EXPECT_NEAR(grad(94, 0), 0.69080009754999E-03, 1E-6);
  EXPECT_NEAR(grad(94, 1), 0.22841697508313E-04, 1E-6);
  EXPECT_NEAR(grad(94, 2), -0.26794673664244E-03, 1E-6);
  EXPECT_NEAR(grad(95, 0), 0.45946479186405E-03, 1E-6);
  EXPECT_NEAR(grad(95, 1), 0.23171388498112E-03, 1E-6);
  EXPECT_NEAR(grad(95, 2), -0.12913238069659E-02, 1E-6);
  EXPECT_NEAR(grad(96, 0), 0.22278079178956E-04, 1E-6);
  EXPECT_NEAR(grad(96, 1), 0.29683680427934E-03, 1E-6);
  EXPECT_NEAR(grad(96, 2), -0.16013416646360E-02, 1E-6);
  EXPECT_NEAR(grad(97, 0), 0.41336335437682E-04, 1E-6);
  EXPECT_NEAR(grad(97, 1), -0.94003775644169E-04, 1E-6);
  EXPECT_NEAR(grad(97, 2), -0.12548312828361E-02, 1E-6);
  EXPECT_NEAR(grad(98, 0), 0.13074465597892E-03, 1E-6);
  EXPECT_NEAR(grad(98, 1), -0.85509764039333E-04, 1E-6);
  EXPECT_NEAR(grad(98, 2), -0.54758526353588E-03, 1E-6);
  EXPECT_NEAR(grad(99, 0), 0.38930609352040E-04, 1E-6);
  EXPECT_NEAR(grad(99, 1), 0.30477564547043E-03, 1E-6);
  EXPECT_NEAR(grad(99, 2), -0.11494674177346E-02, 1E-6);
  EXPECT_NEAR(grad(100, 0), 0.43022405511433E-04, 1E-6);
  EXPECT_NEAR(grad(100, 1), 0.20043850424416E-03, 1E-6);
  EXPECT_NEAR(grad(100, 2), -0.68636705558784E-03, 1E-6);
  EXPECT_NEAR(grad(101, 0), -0.46771954066317E-03, 1E-6);
  EXPECT_NEAR(grad(101, 1), 0.22179019006460E-03, 1E-6);
  EXPECT_NEAR(grad(101, 2), -0.14385576397208E-02, 1E-6);
  EXPECT_NEAR(grad(102, 0), 0.24048306696801E-03, 1E-6);
  EXPECT_NEAR(grad(102, 1), -0.36386043989416E-03, 1E-6);
  EXPECT_NEAR(grad(102, 2), -0.10931017273365E-02, 1E-6);
  EXPECT_NEAR(grad(103, 0), -0.50052743219923E-03, 1E-6);
  EXPECT_NEAR(grad(103, 1), 0.11400129721611E-03, 1E-6);
  EXPECT_NEAR(grad(103, 2), -0.12169828715277E-02, 1E-6);
  EXPECT_NEAR(grad(104, 0), 0.19425076423601E-03, 1E-6);
  EXPECT_NEAR(grad(104, 1), 0.34933249970421E-03, 1E-6);
  EXPECT_NEAR(grad(104, 2), -0.69382387119186E-03, 1E-6);
  EXPECT_NEAR(grad(105, 0), -0.11656252760373E-03, 1E-6);
  EXPECT_NEAR(grad(105, 1), 0.16226375168502E-03, 1E-6);
  EXPECT_NEAR(grad(105, 2), -0.53351140507341E-03, 1E-6);
  EXPECT_NEAR(grad(106, 0), -0.14616372557534E-03, 1E-6);
  EXPECT_NEAR(grad(106, 1), -0.17237919499979E-04, 1E-6);
  EXPECT_NEAR(grad(106, 2), -0.48309210802868E-03, 1E-6);
  EXPECT_NEAR(grad(107, 0), 0.18546506531439E-03, 1E-6);
  EXPECT_NEAR(grad(107, 1), -0.27371222604169E-03, 1E-6);
  EXPECT_NEAR(grad(107, 2), -0.57828050803819E-03, 1E-6);
  EXPECT_NEAR(grad(108, 0), 0.22531793276349E-03, 1E-6);
  EXPECT_NEAR(grad(108, 1), 0.20698949813177E-04, 1E-6);
  EXPECT_NEAR(grad(108, 2), -0.48491863713252E-03, 1E-6);
  EXPECT_NEAR(grad(109, 0), -0.74660458444787E-04, 1E-6);
  EXPECT_NEAR(grad(109, 1), -0.53619943953220E-04, 1E-6);
  EXPECT_NEAR(grad(109, 2), -0.98769664880385E-04, 1E-6);
  EXPECT_NEAR(grad(110, 0), -0.69053165797977E-04, 1E-6);
  EXPECT_NEAR(grad(110, 1), -0.81565364432975E-04, 1E-6);
  EXPECT_NEAR(grad(110, 2), -0.14032655954852E-03, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC gradient correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3ABCGradBP86JacobsenGhost) {
  auto geometry = *systemControllerTwo->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3ABC,
                                                                               std::make_shared<Geometry>(geometry),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 0), -0.23717801640170E-03, 1E-6);
  EXPECT_NEAR(grad(0, 1), -0.28987042072207E-03, 1E-6);
  EXPECT_NEAR(grad(0, 2), 0.18600295030725E-03, 1E-6);
  EXPECT_NEAR(grad(1, 0), 0.27714377137876E-03, 1E-6);
  EXPECT_NEAR(grad(1, 1), -0.64601348693930E-03, 1E-6);
  EXPECT_NEAR(grad(1, 2), 0.22098175770654E-04, 1E-6);
  EXPECT_NEAR(grad(2, 0), -0.63165373673185E-04, 1E-6);
  EXPECT_NEAR(grad(2, 1), -0.13687906670926E-03, 1E-6);
  EXPECT_NEAR(grad(2, 2), 0.43815172780821E-03, 1E-6);
  EXPECT_NEAR(grad(3, 0), -0.12332675677029E-03, 1E-6);
  EXPECT_NEAR(grad(3, 1), -0.22738995146344E-03, 1E-6);
  EXPECT_NEAR(grad(3, 2), 0.38611094473551E-03, 1E-6);
  EXPECT_NEAR(grad(4, 0), 0.64575814860422E-03, 1E-6);
  EXPECT_NEAR(grad(4, 1), 0.12162856186441E-03, 1E-6);
  EXPECT_NEAR(grad(4, 2), 0.78228828377076E-03, 1E-6);
  EXPECT_NEAR(grad(5, 0), 0.16818541189921E-03, 1E-6);
  EXPECT_NEAR(grad(5, 1), -0.92153755939268E-04, 1E-6);
  EXPECT_NEAR(grad(5, 2), -0.97564433713353E-04, 1E-6);
  EXPECT_NEAR(grad(6, 0), -0.28044777249553E-03, 1E-6);
  EXPECT_NEAR(grad(6, 1), 0.44618601379171E-03, 1E-6);
  EXPECT_NEAR(grad(6, 2), 0.25092593198033E-03, 1E-6);
  EXPECT_NEAR(grad(7, 0), 0.24501087991183E-03, 1E-6);
  EXPECT_NEAR(grad(7, 1), -0.69880968688684E-03, 1E-6);
  EXPECT_NEAR(grad(7, 2), -0.24194202305893E-03, 1E-6);
  EXPECT_NEAR(grad(8, 0), 0.59200909543117E-04, 1E-6);
  EXPECT_NEAR(grad(8, 1), -0.20702836543792E-03, 1E-6);
  EXPECT_NEAR(grad(8, 2), -0.97174388308769E-04, 1E-6);
  EXPECT_NEAR(grad(9, 0), -0.25530222548126E-03, 1E-6);
  EXPECT_NEAR(grad(9, 1), -0.50582250253523E-03, 1E-6);
  EXPECT_NEAR(grad(9, 2), 0.11157145257822E-03, 1E-6);
  EXPECT_NEAR(grad(10, 0), -0.16496930249006E-03, 1E-6);
  EXPECT_NEAR(grad(10, 1), -0.12726313793206E-03, 1E-6);
  EXPECT_NEAR(grad(10, 2), 0.68257972726503E-04, 1E-6);
  EXPECT_NEAR(grad(11, 0), -0.76811744317813E-04, 1E-6);
  EXPECT_NEAR(grad(11, 1), -0.17754812228552E-03, 1E-6);
  EXPECT_NEAR(grad(11, 2), 0.11493032735417E-03, 1E-6);
  EXPECT_NEAR(grad(12, 0), -0.42612090021756E-04, 1E-6);
  EXPECT_NEAR(grad(12, 1), -0.81441515719382E-04, 1E-6);
  EXPECT_NEAR(grad(12, 2), -0.89513699952739E-05, 1E-6);
  EXPECT_NEAR(grad(13, 0), -0.21827099739717E-03, 1E-6);
  EXPECT_NEAR(grad(13, 1), -0.19835049585459E-03, 1E-6);
  EXPECT_NEAR(grad(13, 2), -0.52359199059803E-04, 1E-6);
  EXPECT_NEAR(grad(14, 0), 0.12315402180042E-03, 1E-6);
  EXPECT_NEAR(grad(14, 1), -0.32629077803331E-03, 1E-6);
  EXPECT_NEAR(grad(14, 2), -0.20763908944206E-04, 1E-6);
  EXPECT_NEAR(grad(15, 0), 0.46211874845382E-05, 1E-6);
  EXPECT_NEAR(grad(15, 1), -0.85718368010777E-04, 1E-6);
  EXPECT_NEAR(grad(15, 2), -0.51443606585082E-05, 1E-6);
  EXPECT_NEAR(grad(16, 0), 0.16315334500310E-03, 1E-6);
  EXPECT_NEAR(grad(16, 1), -0.41100812582001E-03, 1E-6);
  EXPECT_NEAR(grad(16, 2), -0.63022662401746E-04, 1E-6);
  EXPECT_NEAR(grad(17, 0), 0.10289249631693E-03, 1E-6);
  EXPECT_NEAR(grad(17, 1), -0.15747144737284E-03, 1E-6);
  EXPECT_NEAR(grad(17, 2), -0.34958185107614E-04, 1E-6);
  EXPECT_NEAR(grad(18, 0), 0.67992378875347E-03, 1E-6);
  EXPECT_NEAR(grad(18, 1), -0.20268825107959E-03, 1E-6);
  EXPECT_NEAR(grad(18, 2), 0.56995805361152E-04, 1E-6);
  EXPECT_NEAR(grad(19, 0), -0.75714494690225E-03, 1E-6);
  EXPECT_NEAR(grad(19, 1), -0.10222744775054E-05, 1E-6);
  EXPECT_NEAR(grad(19, 2), 0.73443376080339E-03, 1E-6);
  EXPECT_NEAR(grad(20, 0), 0.32739575379033E-03, 1E-6);
  EXPECT_NEAR(grad(20, 1), -0.12758497091448E-02, 1E-6);
  EXPECT_NEAR(grad(20, 2), 0.64175507129259E-03, 1E-6);
  EXPECT_NEAR(grad(21, 0), -0.52441799028590E-03, 1E-6);
  EXPECT_NEAR(grad(21, 1), -0.12853515735787E-02, 1E-6);
  EXPECT_NEAR(grad(21, 2), 0.44642213284109E-03, 1E-6);
  EXPECT_NEAR(grad(22, 0), 0.21670526023561E-03, 1E-6);
  EXPECT_NEAR(grad(22, 1), -0.28999225726610E-04, 1E-6);
  EXPECT_NEAR(grad(22, 2), 0.14676694325673E-03, 1E-6);
  EXPECT_NEAR(grad(23, 0), -0.12333261329891E-03, 1E-6);
  EXPECT_NEAR(grad(23, 1), -0.71931707092751E-04, 1E-6);
  EXPECT_NEAR(grad(23, 2), 0.44887653605974E-03, 1E-6);
  EXPECT_NEAR(grad(24, 0), -0.79657380099408E-03, 1E-6);
  EXPECT_NEAR(grad(24, 1), 0.67776107497750E-05, 1E-6);
  EXPECT_NEAR(grad(24, 2), 0.50121240439673E-03, 1E-6);
  EXPECT_NEAR(grad(25, 0), -0.65040106103026E-04, 1E-6);
  EXPECT_NEAR(grad(25, 1), 0.12650321767912E-03, 1E-6);
  EXPECT_NEAR(grad(25, 2), 0.99985078988799E-03, 1E-6);
  EXPECT_NEAR(grad(26, 0), -0.10432373823870E-03, 1E-6);
  EXPECT_NEAR(grad(26, 1), -0.55318799754943E-04, 1E-6);
  EXPECT_NEAR(grad(26, 2), 0.58744450419345E-03, 1E-6);
  EXPECT_NEAR(grad(27, 0), -0.24562200600174E-03, 1E-6);
  EXPECT_NEAR(grad(27, 1), -0.20140492475905E-03, 1E-6);
  EXPECT_NEAR(grad(27, 2), 0.47103035784947E-03, 1E-6);
  EXPECT_NEAR(grad(28, 0), -0.28488943965613E-04, 1E-6);
  EXPECT_NEAR(grad(28, 1), -0.20256025714324E-03, 1E-6);
  EXPECT_NEAR(grad(28, 2), 0.98037957845948E-03, 1E-6);
  EXPECT_NEAR(grad(29, 0), 0.11147795588742E-03, 1E-6);
  EXPECT_NEAR(grad(29, 1), 0.20795196148474E-03, 1E-6);
  EXPECT_NEAR(grad(29, 2), -0.15301660100037E-03, 1E-6);
  EXPECT_NEAR(grad(30, 0), 0.13451771804650E-03, 1E-6);
  EXPECT_NEAR(grad(30, 1), -0.21996564946691E-03, 1E-6);
  EXPECT_NEAR(grad(30, 2), -0.61757580869626E-04, 1E-6);
  EXPECT_NEAR(grad(31, 0), 0.57414127146254E-03, 1E-6);
  EXPECT_NEAR(grad(31, 1), -0.37171570269451E-03, 1E-6);
  EXPECT_NEAR(grad(31, 2), -0.11805456544775E-04, 1E-6);
  EXPECT_NEAR(grad(32, 0), 0.32530921494129E-04, 1E-6);
  EXPECT_NEAR(grad(32, 1), 0.61876148944902E-05, 1E-6);
  EXPECT_NEAR(grad(32, 2), -0.99019785973023E-04, 1E-6);
  EXPECT_NEAR(grad(33, 0), -0.17139895446928E-03, 1E-6);
  EXPECT_NEAR(grad(33, 1), -0.78032692433836E-04, 1E-6);
  EXPECT_NEAR(grad(33, 2), -0.15732827256839E-03, 1E-6);
  EXPECT_NEAR(grad(34, 0), -0.49341198602911E-03, 1E-6);
  EXPECT_NEAR(grad(34, 1), -0.61743444855760E-04, 1E-6);
  EXPECT_NEAR(grad(34, 2), 0.36648513219175E-03, 1E-6);
  EXPECT_NEAR(grad(35, 0), -0.26547658033976E-03, 1E-6);
  EXPECT_NEAR(grad(35, 1), 0.71697420789409E-04, 1E-6);
  EXPECT_NEAR(grad(35, 2), 0.64264297133084E-03, 1E-6);
  EXPECT_NEAR(grad(36, 0), -0.38465784590321E-03, 1E-6);
  EXPECT_NEAR(grad(36, 1), -0.15947120297169E-03, 1E-6);
  EXPECT_NEAR(grad(36, 2), 0.57096258492138E-03, 1E-6);
  EXPECT_NEAR(grad(37, 0), -0.54932625298560E-04, 1E-6);
  EXPECT_NEAR(grad(37, 1), 0.49249136281485E-03, 1E-6);
  EXPECT_NEAR(grad(37, 2), 0.34401752509547E-03, 1E-6);
  EXPECT_NEAR(grad(38, 0), 0.96791122176194E-05, 1E-6);
  EXPECT_NEAR(grad(38, 1), -0.25769272209720E-04, 1E-6);
  EXPECT_NEAR(grad(38, 2), -0.88815305236382E-05, 1E-6);
  EXPECT_NEAR(grad(39, 0), 0.77696687296514E-04, 1E-6);
  EXPECT_NEAR(grad(39, 1), -0.22277274905422E-03, 1E-6);
  EXPECT_NEAR(grad(39, 2), -0.24664568728692E-04, 1E-6);
  EXPECT_NEAR(grad(40, 0), 0.26974825652855E-03, 1E-6);
  EXPECT_NEAR(grad(40, 1), -0.28063870364769E-03, 1E-6);
  EXPECT_NEAR(grad(40, 2), 0.16343519077373E-04, 1E-6);
  EXPECT_NEAR(grad(41, 0), 0.28086755687923E-03, 1E-6);
  EXPECT_NEAR(grad(41, 1), 0.39351507801869E-03, 1E-6);
  EXPECT_NEAR(grad(41, 2), 0.76915390321853E-04, 1E-6);
  EXPECT_NEAR(grad(42, 0), -0.59962854513644E-03, 1E-6);
  EXPECT_NEAR(grad(42, 1), 0.28353335400055E-03, 1E-6);
  EXPECT_NEAR(grad(42, 2), 0.21538495303442E-03, 1E-6);
  EXPECT_NEAR(grad(43, 0), 0.62601743123119E-04, 1E-6);
  EXPECT_NEAR(grad(43, 1), 0.43957115723508E-03, 1E-6);
  EXPECT_NEAR(grad(43, 2), 0.26961378468819E-03, 1E-6);
  EXPECT_NEAR(grad(44, 0), 0.45470139826651E-03, 1E-6);
  EXPECT_NEAR(grad(44, 1), -0.10374796869838E-03, 1E-6);
  EXPECT_NEAR(grad(44, 2), 0.13857950034783E-03, 1E-6);
  EXPECT_NEAR(grad(45, 0), 0.14621405858143E-03, 1E-6);
  EXPECT_NEAR(grad(45, 1), 0.44049117441358E-03, 1E-6);
  EXPECT_NEAR(grad(45, 2), 0.11350077007003E-03, 1E-6);
  EXPECT_NEAR(grad(46, 0), -0.49403232081828E-03, 1E-6);
  EXPECT_NEAR(grad(46, 1), 0.53267080440346E-03, 1E-6);
  EXPECT_NEAR(grad(46, 2), 0.23135400719157E-03, 1E-6);
  EXPECT_NEAR(grad(47, 0), -0.12507698657128E-04, 1E-6);
  EXPECT_NEAR(grad(47, 1), -0.91913868078356E-04, 1E-6);
  EXPECT_NEAR(grad(47, 2), -0.15923484875254E-04, 1E-6);
  EXPECT_NEAR(grad(48, 0), -0.59920015451261E-03, 1E-6);
  EXPECT_NEAR(grad(48, 1), 0.76347893970001E-03, 1E-6);
  EXPECT_NEAR(grad(48, 2), 0.28459865567657E-03, 1E-6);
  EXPECT_NEAR(grad(49, 0), -0.54350476610919E-04, 1E-6);
  EXPECT_NEAR(grad(49, 1), 0.44517323767844E-04, 1E-6);
  EXPECT_NEAR(grad(49, 2), 0.37718403727515E-03, 1E-6);
  EXPECT_NEAR(grad(50, 0), -0.47841592399576E-03, 1E-6);
  EXPECT_NEAR(grad(50, 1), 0.23655916460177E-03, 1E-6);
  EXPECT_NEAR(grad(50, 2), 0.32110226559319E-03, 1E-6);
  EXPECT_NEAR(grad(51, 0), -0.11920946235248E-03, 1E-6);
  EXPECT_NEAR(grad(51, 1), 0.29878907176467E-04, 1E-6);
  EXPECT_NEAR(grad(51, 2), -0.96850931160190E-04, 1E-6);
  EXPECT_NEAR(grad(52, 0), 0.61330818744228E-03, 1E-6);
  EXPECT_NEAR(grad(52, 1), 0.17975767686066E-03, 1E-6);
  EXPECT_NEAR(grad(52, 2), 0.52938062027837E-04, 1E-6);
  EXPECT_NEAR(grad(53, 0), -0.44299839413577E-04, 1E-6);
  EXPECT_NEAR(grad(53, 1), -0.88609368332008E-04, 1E-6);
  EXPECT_NEAR(grad(53, 2), -0.41585309230315E-04, 1E-6);
  EXPECT_NEAR(grad(54, 0), 0.10150763871525E-02, 1E-6);
  EXPECT_NEAR(grad(54, 1), 0.18756779178569E-03, 1E-6);
  EXPECT_NEAR(grad(54, 2), 0.28681211423454E-03, 1E-6);
  EXPECT_NEAR(grad(55, 0), -0.17655992229250E-03, 1E-6);
  EXPECT_NEAR(grad(55, 1), 0.91438478076860E-04, 1E-6);
  EXPECT_NEAR(grad(55, 2), -0.48379822979946E-04, 1E-6);
  EXPECT_NEAR(grad(56, 0), 0.14763856630641E-03, 1E-6);
  EXPECT_NEAR(grad(56, 1), -0.80163757160533E-04, 1E-6);
  EXPECT_NEAR(grad(56, 2), -0.46900095211109E-04, 1E-6);
  EXPECT_NEAR(grad(57, 0), 0.95910512539158E-04, 1E-6);
  EXPECT_NEAR(grad(57, 1), -0.99726589290551E-04, 1E-6);
  EXPECT_NEAR(grad(57, 2), 0.96102877321659E-04, 1E-6);
  EXPECT_NEAR(grad(58, 0), 0.33520552531891E-03, 1E-6);
  EXPECT_NEAR(grad(58, 1), 0.29913547620154E-04, 1E-6);
  EXPECT_NEAR(grad(58, 2), -0.14124748328383E-03, 1E-6);
  EXPECT_NEAR(grad(59, 0), 0.10205051868591E-03, 1E-6);
  EXPECT_NEAR(grad(59, 1), -0.15363718400346E-04, 1E-6);
  EXPECT_NEAR(grad(59, 2), -0.73215738010412E-04, 1E-6);
  EXPECT_NEAR(grad(60, 0), 0.31206295463529E-04, 1E-6);
  EXPECT_NEAR(grad(60, 1), -0.22748496722329E-04, 1E-6);
  EXPECT_NEAR(grad(60, 2), -0.72337845666214E-04, 1E-6);
  EXPECT_NEAR(grad(61, 0), 0.18042733573728E-03, 1E-6);
  EXPECT_NEAR(grad(61, 1), -0.45779690392349E-04, 1E-6);
  EXPECT_NEAR(grad(61, 2), 0.31427285253189E-03, 1E-6);
  EXPECT_NEAR(grad(62, 0), -0.17226594573020E-04, 1E-6);
  EXPECT_NEAR(grad(62, 1), -0.45591584626299E-04, 1E-6);
  EXPECT_NEAR(grad(62, 2), 0.60752993934197E-04, 1E-6);
  EXPECT_NEAR(grad(63, 0), 0.38205654772853E-04, 1E-6);
  EXPECT_NEAR(grad(63, 1), -0.50010959823938E-04, 1E-6);
  EXPECT_NEAR(grad(63, 2), 0.10916974670623E-03, 1E-6);
  EXPECT_NEAR(grad(64, 0), -0.33748835115731E-05, 1E-6);
  EXPECT_NEAR(grad(64, 1), -0.13484587103423E-03, 1E-6);
  EXPECT_NEAR(grad(64, 2), -0.12487578225330E-03, 1E-6);
  EXPECT_NEAR(grad(65, 0), 0.79639567517741E-04, 1E-6);
  EXPECT_NEAR(grad(65, 1), -0.38404958353020E-04, 1E-6);
  EXPECT_NEAR(grad(65, 2), 0.25194631743046E-04, 1E-6);
  EXPECT_NEAR(grad(66, 0), -0.62313269442493E-04, 1E-6);
  EXPECT_NEAR(grad(66, 1), -0.16301505345603E-03, 1E-6);
  EXPECT_NEAR(grad(66, 2), 0.49945375463408E-04, 1E-6);
  EXPECT_NEAR(grad(67, 0), 0.46385309389828E-04, 1E-6);
  EXPECT_NEAR(grad(67, 1), 0.44811757577484E-04, 1E-6);
  EXPECT_NEAR(grad(67, 2), 0.11337491216874E-03, 1E-6);
  EXPECT_NEAR(grad(68, 0), 0.10020531487048E-03, 1E-6);
  EXPECT_NEAR(grad(68, 1), 0.72029779165086E-04, 1E-6);
  EXPECT_NEAR(grad(68, 2), -0.47447693661869E-04, 1E-6);
  EXPECT_NEAR(grad(69, 0), -0.71885733670794E-04, 1E-6);
  EXPECT_NEAR(grad(69, 1), 0.79197431349367E-04, 1E-6);
  EXPECT_NEAR(grad(69, 2), -0.29985666081673E-04, 1E-6);
  EXPECT_NEAR(grad(70, 0), 0.66657750351279E-03, 1E-6);
  EXPECT_NEAR(grad(70, 1), 0.80348112451879E-05, 1E-6);
  EXPECT_NEAR(grad(70, 2), 0.46509998380525E-03, 1E-6);
  EXPECT_NEAR(grad(71, 0), 0.59686857378042E-03, 1E-6);
  EXPECT_NEAR(grad(71, 1), 0.87156782141222E-04, 1E-6);
  EXPECT_NEAR(grad(71, 2), 0.28834722803016E-03, 1E-6);
  EXPECT_NEAR(grad(72, 0), 0.43457789773417E-04, 1E-6);
  EXPECT_NEAR(grad(72, 1), 0.22939963012866E-03, 1E-6);
  EXPECT_NEAR(grad(72, 2), 0.31347335281305E-03, 1E-6);
  EXPECT_NEAR(grad(73, 0), 0.53666331383611E-04, 1E-6);
  EXPECT_NEAR(grad(73, 1), 0.31076405924674E-03, 1E-6);
  EXPECT_NEAR(grad(73, 2), -0.21786753826681E-03, 1E-6);
  EXPECT_NEAR(grad(74, 0), 0.43999351148239E-03, 1E-6);
  EXPECT_NEAR(grad(74, 1), 0.11509127644925E-03, 1E-6);
  EXPECT_NEAR(grad(74, 2), -0.71534603996133E-04, 1E-6);
  EXPECT_NEAR(grad(75, 0), 0.43179822590886E-03, 1E-6);
  EXPECT_NEAR(grad(75, 1), 0.57500761210677E-04, 1E-6);
  EXPECT_NEAR(grad(75, 2), -0.16973737775986E-03, 1E-6);
  EXPECT_NEAR(grad(76, 0), 0.54492436234436E-04, 1E-6);
  EXPECT_NEAR(grad(76, 1), 0.13388535195938E-03, 1E-6);
  EXPECT_NEAR(grad(76, 2), 0.33878565294576E-03, 1E-6);
  EXPECT_NEAR(grad(77, 0), -0.34718879492883E-03, 1E-6);
  EXPECT_NEAR(grad(77, 1), 0.23069474191604E-03, 1E-6);
  EXPECT_NEAR(grad(77, 2), 0.23294024710741E-03, 1E-6);
  EXPECT_NEAR(grad(78, 0), -0.37065763710889E-03, 1E-6);
  EXPECT_NEAR(grad(78, 1), 0.32799122541849E-03, 1E-6);
  EXPECT_NEAR(grad(78, 2), 0.25589878604821E-03, 1E-6);
  EXPECT_NEAR(grad(79, 0), 0.83248427105235E-05, 1E-6);
  EXPECT_NEAR(grad(79, 1), 0.97885464966770E-04, 1E-6);
  EXPECT_NEAR(grad(79, 2), -0.64106444935496E-04, 1E-6);
  EXPECT_NEAR(grad(80, 0), -0.68012435609132E-04, 1E-6);
  EXPECT_NEAR(grad(80, 1), 0.97248355558922E-04, 1E-6);
  EXPECT_NEAR(grad(80, 2), 0.12085109824906E-03, 1E-6);
  EXPECT_NEAR(grad(81, 0), 0.10067005866307E-03, 1E-6);
  EXPECT_NEAR(grad(81, 1), 0.31985509955836E-04, 1E-6);
  EXPECT_NEAR(grad(81, 2), 0.24142441382561E-04, 1E-6);
  EXPECT_NEAR(grad(82, 0), -0.39592147417943E-03, 1E-6);
  EXPECT_NEAR(grad(82, 1), 0.51487090788829E-03, 1E-6);
  EXPECT_NEAR(grad(82, 2), -0.15697821907618E-03, 1E-6);
  EXPECT_NEAR(grad(83, 0), -0.40574768483625E-03, 1E-6);
  EXPECT_NEAR(grad(83, 1), 0.47384114765033E-03, 1E-6);
  EXPECT_NEAR(grad(83, 2), -0.46267285246846E-04, 1E-6);
  EXPECT_NEAR(grad(84, 0), 0.22385332802663E-03, 1E-6);
  EXPECT_NEAR(grad(84, 1), 0.47977995413444E-03, 1E-6);
  EXPECT_NEAR(grad(84, 2), 0.33560043398743E-04, 1E-6);
  EXPECT_NEAR(grad(85, 0), -0.41113959268474E-03, 1E-6);
  EXPECT_NEAR(grad(85, 1), -0.45939920943882E-04, 1E-6);
  EXPECT_NEAR(grad(85, 2), 0.46548649050262E-03, 1E-6);
  EXPECT_NEAR(grad(86, 0), -0.12969186143636E-03, 1E-6);
  EXPECT_NEAR(grad(86, 1), 0.37120288362057E-04, 1E-6);
  EXPECT_NEAR(grad(86, 2), 0.20356776238145E-04, 1E-6);
  EXPECT_NEAR(grad(87, 0), -0.76099512050028E-04, 1E-6);
  EXPECT_NEAR(grad(87, 1), -0.10620104574951E-03, 1E-6);
  EXPECT_NEAR(grad(87, 2), -0.70375515396661E-04, 1E-6);
  EXPECT_NEAR(grad(88, 0), -0.29991445279378E-03, 1E-6);
  EXPECT_NEAR(grad(88, 1), 0.84537364662448E-04, 1E-6);
  EXPECT_NEAR(grad(88, 2), 0.21364953246601E-03, 1E-6);
  EXPECT_NEAR(grad(89, 0), -0.10196937659892E-03, 1E-6);
  EXPECT_NEAR(grad(89, 1), -0.56884565505157E-05, 1E-6);
  EXPECT_NEAR(grad(89, 2), 0.68544773544331E-04, 1E-6);
  EXPECT_NEAR(grad(90, 0), -0.51941504567816E-04, 1E-6);
  EXPECT_NEAR(grad(90, 1), -0.38852948426937E-04, 1E-6);
  EXPECT_NEAR(grad(90, 2), 0.55508831642056E-04, 1E-6);
  EXPECT_NEAR(grad(91, 0), -0.44243048904823E-03, 1E-6);
  EXPECT_NEAR(grad(91, 1), 0.54027733043081E-03, 1E-6);
  EXPECT_NEAR(grad(91, 2), -0.38226138152842E-04, 1E-6);
  EXPECT_NEAR(grad(92, 0), -0.73391765815008E-04, 1E-6);
  EXPECT_NEAR(grad(92, 1), 0.76554117053676E-04, 1E-6);
  EXPECT_NEAR(grad(92, 2), -0.78239905949963E-05, 1E-6);
  EXPECT_NEAR(grad(93, 0), -0.85672561524956E-04, 1E-6);
  EXPECT_NEAR(grad(93, 1), 0.72740622131106E-04, 1E-6);
  EXPECT_NEAR(grad(93, 2), -0.69150038840129E-04, 1E-6);
  EXPECT_NEAR(grad(94, 0), 0.69080009754999E-03, 1E-6);
  EXPECT_NEAR(grad(94, 1), 0.22841697508313E-04, 1E-6);
  EXPECT_NEAR(grad(94, 2), -0.26794673664244E-03, 1E-6);
  EXPECT_NEAR(grad(95, 0), 0.45946479186405E-03, 1E-6);
  EXPECT_NEAR(grad(95, 1), 0.23171388498112E-03, 1E-6);
  EXPECT_NEAR(grad(95, 2), -0.12913238069659E-02, 1E-6);
  EXPECT_NEAR(grad(96, 0), 0.22278079178956E-04, 1E-6);
  EXPECT_NEAR(grad(96, 1), 0.29683680427934E-03, 1E-6);
  EXPECT_NEAR(grad(96, 2), -0.16013416646360E-02, 1E-6);
  EXPECT_NEAR(grad(97, 0), 0.41336335437682E-04, 1E-6);
  EXPECT_NEAR(grad(97, 1), -0.94003775644169E-04, 1E-6);
  EXPECT_NEAR(grad(97, 2), -0.12548312828361E-02, 1E-6);
  EXPECT_NEAR(grad(98, 0), 0.13074465597892E-03, 1E-6);
  EXPECT_NEAR(grad(98, 1), -0.85509764039333E-04, 1E-6);
  EXPECT_NEAR(grad(98, 2), -0.54758526353588E-03, 1E-6);
  EXPECT_NEAR(grad(99, 0), 0.38930609352040E-04, 1E-6);
  EXPECT_NEAR(grad(99, 1), 0.30477564547043E-03, 1E-6);
  EXPECT_NEAR(grad(99, 2), -0.11494674177346E-02, 1E-6);
  EXPECT_NEAR(grad(100, 0), 0.43022405511433E-04, 1E-6);
  EXPECT_NEAR(grad(100, 1), 0.20043850424416E-03, 1E-6);
  EXPECT_NEAR(grad(100, 2), -0.68636705558784E-03, 1E-6);
  EXPECT_NEAR(grad(101, 0), -0.46771954066317E-03, 1E-6);
  EXPECT_NEAR(grad(101, 1), 0.22179019006460E-03, 1E-6);
  EXPECT_NEAR(grad(101, 2), -0.14385576397208E-02, 1E-6);
  EXPECT_NEAR(grad(102, 0), 0.24048306696801E-03, 1E-6);
  EXPECT_NEAR(grad(102, 1), -0.36386043989416E-03, 1E-6);
  EXPECT_NEAR(grad(102, 2), -0.10931017273365E-02, 1E-6);
  EXPECT_NEAR(grad(103, 0), -0.50052743219923E-03, 1E-6);
  EXPECT_NEAR(grad(103, 1), 0.11400129721611E-03, 1E-6);
  EXPECT_NEAR(grad(103, 2), -0.12169828715277E-02, 1E-6);
  EXPECT_NEAR(grad(104, 0), 0.19425076423601E-03, 1E-6);
  EXPECT_NEAR(grad(104, 1), 0.34933249970421E-03, 1E-6);
  EXPECT_NEAR(grad(104, 2), -0.69382387119186E-03, 1E-6);
  EXPECT_NEAR(grad(105, 0), -0.11656252760373E-03, 1E-6);
  EXPECT_NEAR(grad(105, 1), 0.16226375168502E-03, 1E-6);
  EXPECT_NEAR(grad(105, 2), -0.53351140507341E-03, 1E-6);
  EXPECT_NEAR(grad(106, 0), -0.14616372557534E-03, 1E-6);
  EXPECT_NEAR(grad(106, 1), -0.17237919499979E-04, 1E-6);
  EXPECT_NEAR(grad(106, 2), -0.48309210802868E-03, 1E-6);
  EXPECT_NEAR(grad(107, 0), 0.18546506531439E-03, 1E-6);
  EXPECT_NEAR(grad(107, 1), -0.27371222604169E-03, 1E-6);
  EXPECT_NEAR(grad(107, 2), -0.57828050803819E-03, 1E-6);
  EXPECT_NEAR(grad(108, 0), 0.22531793276349E-03, 1E-6);
  EXPECT_NEAR(grad(108, 1), 0.20698949813177E-04, 1E-6);
  EXPECT_NEAR(grad(108, 2), -0.48491863713252E-03, 1E-6);
  EXPECT_NEAR(grad(109, 0), -0.74660458444787E-04, 1E-6);
  EXPECT_NEAR(grad(109, 1), -0.53619943953220E-04, 1E-6);
  EXPECT_NEAR(grad(109, 2), -0.98769664880385E-04, 1E-6);
  EXPECT_NEAR(grad(110, 0), -0.69053165797977E-04, 1E-6);
  EXPECT_NEAR(grad(110, 1), -0.81565364432975E-04, 1E-6);
  EXPECT_NEAR(grad(110, 2), -0.14032655954852E-03, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ gradient correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJGradBP86CO) {
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJ,
                                                                               systemControllerOne->getGeometry(),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 2), -0.73442092899124E-6, 1E-8);
  EXPECT_NEAR(grad(1, 2), 0.73442092899124E-6, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ gradient correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJGradBP86Jacobsen) {
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJ,
                                                                               systemControllerTwo->getGeometry(),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 0), -0.44641271320771E-03, 1E-6);
  EXPECT_NEAR(grad(0, 1), -0.34836653878647E-03, 1E-6);
  EXPECT_NEAR(grad(0, 2), 0.23794535455227E-03, 1E-6);
  EXPECT_NEAR(grad(1, 0), 0.45577469855532E-03, 1E-6);
  EXPECT_NEAR(grad(1, 1), -0.51216118328020E-03, 1E-6);
  EXPECT_NEAR(grad(1, 2), -0.22358523681505E-03, 1E-6);
  EXPECT_NEAR(grad(2, 0), -0.16787504556191E-03, 1E-6);
  EXPECT_NEAR(grad(2, 1), -0.18620322081428E-03, 1E-6);
  EXPECT_NEAR(grad(2, 2), 0.34017282987155E-03, 1E-6);
  EXPECT_NEAR(grad(3, 0), -0.60630814302006E-03, 1E-6);
  EXPECT_NEAR(grad(3, 1), -0.74541811578157E-03, 1E-6);
  EXPECT_NEAR(grad(3, 2), 0.16413329853933E-03, 1E-6);
  EXPECT_NEAR(grad(4, 0), -0.26284075145145E-03, 1E-6);
  EXPECT_NEAR(grad(4, 1), -0.89624101639072E-04, 1E-6);
  EXPECT_NEAR(grad(4, 2), 0.40298169397188E-03, 1E-6);
  EXPECT_NEAR(grad(5, 0), 0.33709084661949E-03, 1E-6);
  EXPECT_NEAR(grad(5, 1), -0.81644805432351E-03, 1E-6);
  EXPECT_NEAR(grad(5, 2), -0.18345467497986E-03, 1E-6);
  EXPECT_NEAR(grad(6, 0), 0.43843889981468E-03, 1E-6);
  EXPECT_NEAR(grad(6, 1), -0.16674409600683E-03, 1E-6);
  EXPECT_NEAR(grad(6, 2), -0.13161594754476E-03, 1E-6);
  EXPECT_NEAR(grad(7, 0), 0.24192941011567E-03, 1E-6);
  EXPECT_NEAR(grad(7, 1), -0.16832764022001E-03, 1E-6);
  EXPECT_NEAR(grad(7, 2), -0.30474423492655E-03, 1E-6);
  EXPECT_NEAR(grad(8, 0), -0.34683727112197E-03, 1E-6);
  EXPECT_NEAR(grad(8, 1), -0.71473634375525E-03, 1E-6);
  EXPECT_NEAR(grad(8, 2), 0.23349652243208E-04, 1E-6);
  EXPECT_NEAR(grad(9, 0), -0.25045716376056E-03, 1E-6);
  EXPECT_NEAR(grad(9, 1), -0.35481675979231E-03, 1E-6);
  EXPECT_NEAR(grad(9, 2), -0.47783472255487E-04, 1E-6);
  EXPECT_NEAR(grad(10, 0), -0.35809598349183E-03, 1E-6);
  EXPECT_NEAR(grad(10, 1), -0.24140359149052E-03, 1E-6);
  EXPECT_NEAR(grad(10, 2), 0.12909347585806E-03, 1E-6);
  EXPECT_NEAR(grad(11, 0), 0.24560285483310E-04, 1E-6);
  EXPECT_NEAR(grad(11, 1), -0.74514937647724E-03, 1E-6);
  EXPECT_NEAR(grad(11, 2), -0.21485206615971E-03, 1E-6);
  EXPECT_NEAR(grad(12, 0), -0.20367754473177E-03, 1E-6);
  EXPECT_NEAR(grad(12, 1), -0.26269893634963E-03, 1E-6);
  EXPECT_NEAR(grad(12, 2), -0.24559082119115E-04, 1E-6);
  EXPECT_NEAR(grad(13, 0), -0.18259474524204E-03, 1E-6);
  EXPECT_NEAR(grad(13, 1), -0.30353968247095E-03, 1E-6);
  EXPECT_NEAR(grad(13, 2), 0.14374385178206E-03, 1E-6);
  EXPECT_NEAR(grad(14, 0), 0.57658195466981E-04, 1E-6);
  EXPECT_NEAR(grad(14, 1), -0.31249288904384E-03, 1E-6);
  EXPECT_NEAR(grad(14, 2), -0.23357036167616E-03, 1E-6);
  EXPECT_NEAR(grad(15, 0), 0.56558654485810E-04, 1E-6);
  EXPECT_NEAR(grad(15, 1), -0.31072153068333E-03, 1E-6);
  EXPECT_NEAR(grad(15, 2), -0.67688361000117E-04, 1E-6);
  EXPECT_NEAR(grad(16, 0), 0.17678020442675E-03, 1E-6);
  EXPECT_NEAR(grad(16, 1), -0.40706000522237E-03, 1E-6);
  EXPECT_NEAR(grad(16, 2), 0.11029162574860E-03, 1E-6);
  EXPECT_NEAR(grad(17, 0), 0.22897539601370E-03, 1E-6);
  EXPECT_NEAR(grad(17, 1), -0.30429747189709E-03, 1E-6);
  EXPECT_NEAR(grad(17, 2), -0.14996223947871E-03, 1E-6);
  EXPECT_NEAR(grad(18, 0), 0.73903797560823E-03, 1E-6);
  EXPECT_NEAR(grad(18, 1), -0.24943995399781E-03, 1E-6);
  EXPECT_NEAR(grad(18, 2), -0.78063010908874E-04, 1E-6);
  EXPECT_NEAR(grad(19, 0), -0.61740474005408E-03, 1E-6);
  EXPECT_NEAR(grad(19, 1), -0.28919030071639E-05, 1E-6);
  EXPECT_NEAR(grad(19, 2), 0.81457294885140E-03, 1E-6);
  EXPECT_NEAR(grad(20, 0), 0.34417247784298E-03, 1E-6);
  EXPECT_NEAR(grad(20, 1), -0.61542332660393E-03, 1E-6);
  EXPECT_NEAR(grad(20, 2), 0.27938015674633E-03, 1E-6);
  EXPECT_NEAR(grad(21, 0), -0.37474625438297E-03, 1E-6);
  EXPECT_NEAR(grad(21, 1), -0.77615262525778E-03, 1E-6);
  EXPECT_NEAR(grad(21, 2), 0.25564030684309E-02, 1E-6);
  EXPECT_NEAR(grad(22, 0), 0.28156240010851E-03, 1E-6);
  EXPECT_NEAR(grad(22, 1), -0.55492643569041E-04, 1E-6);
  EXPECT_NEAR(grad(22, 2), -0.13014041559225E-03, 1E-6);
  EXPECT_NEAR(grad(23, 0), -0.13103105018074E-03, 1E-6);
  EXPECT_NEAR(grad(23, 1), 0.27799316518425E-04, 1E-6);
  EXPECT_NEAR(grad(23, 2), 0.15646036491417E-03, 1E-6);
  EXPECT_NEAR(grad(24, 0), -0.89193053184762E-03, 1E-6);
  EXPECT_NEAR(grad(24, 1), -0.15967087536234E-03, 1E-6);
  EXPECT_NEAR(grad(24, 2), 0.33601158334330E-03, 1E-6);
  EXPECT_NEAR(grad(25, 0), -0.56560909635577E-03, 1E-6);
  EXPECT_NEAR(grad(25, 1), 0.21152218527148E-03, 1E-6);
  EXPECT_NEAR(grad(25, 2), 0.97275710042093E-03, 1E-6);
  EXPECT_NEAR(grad(26, 0), -0.18379499207937E-03, 1E-6);
  EXPECT_NEAR(grad(26, 1), 0.18940050876816E-03, 1E-6);
  EXPECT_NEAR(grad(26, 2), 0.55814736958811E-03, 1E-6);
  EXPECT_NEAR(grad(27, 0), -0.28119703066395E-03, 1E-6);
  EXPECT_NEAR(grad(27, 1), -0.16966670570161E-03, 1E-6);
  EXPECT_NEAR(grad(27, 2), 0.30675614303623E-03, 1E-6);
  EXPECT_NEAR(grad(28, 0), 0.94087469680944E-04, 1E-6);
  EXPECT_NEAR(grad(28, 1), -0.16441875141781E-03, 1E-6);
  EXPECT_NEAR(grad(28, 2), 0.11479233636263E-02, 1E-6);
  EXPECT_NEAR(grad(29, 0), 0.62173835576336E-03, 1E-6);
  EXPECT_NEAR(grad(29, 1), -0.81389700838245E-04, 1E-6);
  EXPECT_NEAR(grad(29, 2), -0.22821199827511E-04, 1E-6);
  EXPECT_NEAR(grad(30, 0), 0.21843889733165E-03, 1E-6);
  EXPECT_NEAR(grad(30, 1), -0.22654742282420E-03, 1E-6);
  EXPECT_NEAR(grad(30, 2), -0.56436923651591E-04, 1E-6);
  EXPECT_NEAR(grad(31, 0), 0.68903716509116E-03, 1E-6);
  EXPECT_NEAR(grad(31, 1), -0.58973713225538E-03, 1E-6);
  EXPECT_NEAR(grad(31, 2), -0.72085095762718E-04, 1E-6);
  EXPECT_NEAR(grad(32, 0), 0.23521805385171E-03, 1E-6);
  EXPECT_NEAR(grad(32, 1), 0.91434535733175E-04, 1E-6);
  EXPECT_NEAR(grad(32, 2), -0.27059342000401E-04, 1E-6);
  EXPECT_NEAR(grad(33, 0), -0.28263706522879E-03, 1E-6);
  EXPECT_NEAR(grad(33, 1), 0.25238722703096E-04, 1E-6);
  EXPECT_NEAR(grad(33, 2), -0.17971438080203E-03, 1E-6);
  EXPECT_NEAR(grad(34, 0), -0.14349099763395E-03, 1E-6);
  EXPECT_NEAR(grad(34, 1), -0.36734470046452E-03, 1E-6);
  EXPECT_NEAR(grad(34, 2), 0.33238504618144E-03, 1E-6);
  EXPECT_NEAR(grad(35, 0), -0.39969396352155E-03, 1E-6);
  EXPECT_NEAR(grad(35, 1), 0.35905083869383E-04, 1E-6);
  EXPECT_NEAR(grad(35, 2), 0.61785762546793E-03, 1E-6);
  EXPECT_NEAR(grad(36, 0), -0.26006224750122E-03, 1E-6);
  EXPECT_NEAR(grad(36, 1), -0.21238714177715E-03, 1E-6);
  EXPECT_NEAR(grad(36, 2), 0.33696349787377E-03, 1E-6);
  EXPECT_NEAR(grad(37, 0), -0.11569947728813E-03, 1E-6);
  EXPECT_NEAR(grad(37, 1), 0.53095339335688E-03, 1E-6);
  EXPECT_NEAR(grad(37, 2), 0.33858910927632E-03, 1E-6);
  EXPECT_NEAR(grad(38, 0), 0.88463938561081E-04, 1E-6);
  EXPECT_NEAR(grad(38, 1), -0.15787530828135E-03, 1E-6);
  EXPECT_NEAR(grad(38, 2), -0.19161632052436E-04, 1E-6);
  EXPECT_NEAR(grad(39, 0), -0.21442293573086E-05, 1E-6);
  EXPECT_NEAR(grad(39, 1), -0.25193569212640E-03, 1E-6);
  EXPECT_NEAR(grad(39, 2), -0.57739268764463E-04, 1E-6);
  EXPECT_NEAR(grad(40, 0), -0.16209776098162E-03, 1E-6);
  EXPECT_NEAR(grad(40, 1), -0.34364021238956E-03, 1E-6);
  EXPECT_NEAR(grad(40, 2), -0.14763058402152E-03, 1E-6);
  EXPECT_NEAR(grad(41, 0), 0.36577755480510E-03, 1E-6);
  EXPECT_NEAR(grad(41, 1), 0.39471668527540E-03, 1E-6);
  EXPECT_NEAR(grad(41, 2), 0.99053757893348E-04, 1E-6);
  EXPECT_NEAR(grad(42, 0), -0.54886977155488E-03, 1E-6);
  EXPECT_NEAR(grad(42, 1), 0.17715312675855E-03, 1E-6);
  EXPECT_NEAR(grad(42, 2), 0.20499676161472E-03, 1E-6);
  EXPECT_NEAR(grad(43, 0), -0.88552448494570E-04, 1E-6);
  EXPECT_NEAR(grad(43, 1), 0.24207094059259E-03, 1E-6);
  EXPECT_NEAR(grad(43, 2), 0.11389353518825E-03, 1E-6);
  EXPECT_NEAR(grad(44, 0), 0.34626618506442E-03, 1E-6);
  EXPECT_NEAR(grad(44, 1), -0.16659339344682E-03, 1E-6);
  EXPECT_NEAR(grad(44, 2), 0.86740299578476E-04, 1E-6);
  EXPECT_NEAR(grad(45, 0), 0.16367720529031E-03, 1E-6);
  EXPECT_NEAR(grad(45, 1), 0.14800220164615E-03, 1E-6);
  EXPECT_NEAR(grad(45, 2), 0.72661276468407E-04, 1E-6);
  EXPECT_NEAR(grad(46, 0), -0.21447756010488E-03, 1E-6);
  EXPECT_NEAR(grad(46, 1), 0.38814549380924E-03, 1E-6);
  EXPECT_NEAR(grad(46, 2), 0.76811641717841E-03, 1E-6);
  EXPECT_NEAR(grad(47, 0), -0.43761247618718E-04, 1E-6);
  EXPECT_NEAR(grad(47, 1), 0.61632870436173E-03, 1E-6);
  EXPECT_NEAR(grad(47, 2), 0.27683860053273E-03, 1E-6);
  EXPECT_NEAR(grad(48, 0), 0.15098279301011E-04, 1E-6);
  EXPECT_NEAR(grad(48, 1), 0.93557160546381E-03, 1E-6);
  EXPECT_NEAR(grad(48, 2), -0.33050379082729E-04, 1E-6);
  EXPECT_NEAR(grad(49, 0), -0.74596918775100E-03, 1E-6);
  EXPECT_NEAR(grad(49, 1), -0.15043852117657E-03, 1E-6);
  EXPECT_NEAR(grad(49, 2), 0.13005920216928E-03, 1E-6);
  EXPECT_NEAR(grad(50, 0), -0.56163311595122E-03, 1E-6);
  EXPECT_NEAR(grad(50, 1), 0.62036182589626E-03, 1E-6);
  EXPECT_NEAR(grad(50, 2), -0.15904402086311E-03, 1E-6);
  EXPECT_NEAR(grad(51, 0), -0.52961861065490E-03, 1E-6);
  EXPECT_NEAR(grad(51, 1), 0.18542169563333E-03, 1E-6);
  EXPECT_NEAR(grad(51, 2), 0.37920803427412E-03, 1E-6);
  EXPECT_NEAR(grad(52, 0), 0.45461184156401E-03, 1E-6);
  EXPECT_NEAR(grad(52, 1), 0.55633929051552E-03, 1E-6);
  EXPECT_NEAR(grad(52, 2), -0.44506619795246E-03, 1E-6);
  EXPECT_NEAR(grad(53, 0), 0.33787323145827E-03, 1E-6);
  EXPECT_NEAR(grad(53, 1), 0.51472539709650E-03, 1E-6);
  EXPECT_NEAR(grad(53, 2), 0.17274932015987E-03, 1E-6);
  EXPECT_NEAR(grad(54, 0), 0.45860007196877E-03, 1E-6);
  EXPECT_NEAR(grad(54, 1), 0.41609792638722E-03, 1E-6);
  EXPECT_NEAR(grad(54, 2), 0.66520654768515E-03, 1E-6);
  EXPECT_NEAR(grad(55, 0), 0.26602599269625E-03, 1E-6);
  EXPECT_NEAR(grad(55, 1), -0.51313466515424E-03, 1E-6);
  EXPECT_NEAR(grad(55, 2), 0.60562121870499E-05, 1E-6);
  EXPECT_NEAR(grad(56, 0), 0.33776188964159E-03, 1E-6);
  EXPECT_NEAR(grad(56, 1), -0.77598266465638E-04, 1E-6);
  EXPECT_NEAR(grad(56, 2), 0.54284619083235E-03, 1E-6);
  EXPECT_NEAR(grad(57, 0), 0.59816288304327E-03, 1E-6);
  EXPECT_NEAR(grad(57, 1), 0.50024928219000E-04, 1E-6);
  EXPECT_NEAR(grad(57, 2), -0.23037303670062E-03, 1E-6);
  EXPECT_NEAR(grad(58, 0), 0.32585862194356E-03, 1E-6);
  EXPECT_NEAR(grad(58, 1), 0.38022784558631E-04, 1E-6);
  EXPECT_NEAR(grad(58, 2), -0.29996961476651E-03, 1E-6);
  EXPECT_NEAR(grad(59, 0), 0.32601971990793E-03, 1E-6);
  EXPECT_NEAR(grad(59, 1), 0.39677021961920E-04, 1E-6);
  EXPECT_NEAR(grad(59, 2), -0.10822338015477E-03, 1E-6);
  EXPECT_NEAR(grad(60, 0), 0.35879238987349E-03, 1E-6);
  EXPECT_NEAR(grad(60, 1), 0.13551136441295E-03, 1E-6);
  EXPECT_NEAR(grad(60, 2), -0.14024983294786E-03, 1E-6);
  EXPECT_NEAR(grad(61, 0), 0.78913900690811E-04, 1E-6);
  EXPECT_NEAR(grad(61, 1), -0.82239511390362E-04, 1E-6);
  EXPECT_NEAR(grad(61, 2), 0.43168740640081E-03, 1E-6);
  EXPECT_NEAR(grad(62, 0), 0.19101918722848E-03, 1E-6);
  EXPECT_NEAR(grad(62, 1), 0.54270232450281E-04, 1E-6);
  EXPECT_NEAR(grad(62, 2), 0.35754976768672E-03, 1E-6);
  EXPECT_NEAR(grad(63, 0), 0.19024119042316E-03, 1E-6);
  EXPECT_NEAR(grad(63, 1), -0.25233440737570E-04, 1E-6);
  EXPECT_NEAR(grad(63, 2), 0.28806297248064E-03, 1E-6);
  EXPECT_NEAR(grad(64, 0), 0.60264802794365E-04, 1E-6);
  EXPECT_NEAR(grad(64, 1), -0.36224162920476E-03, 1E-6);
  EXPECT_NEAR(grad(64, 2), 0.84571549113634E-04, 1E-6);
  EXPECT_NEAR(grad(65, 0), 0.18434142585794E-03, 1E-6);
  EXPECT_NEAR(grad(65, 1), -0.28307699026813E-03, 1E-6);
  EXPECT_NEAR(grad(65, 2), 0.17970762842745E-04, 1E-6);
  EXPECT_NEAR(grad(66, 0), 0.13334669353749E-03, 1E-6);
  EXPECT_NEAR(grad(66, 1), -0.32579829125995E-03, 1E-6);
  EXPECT_NEAR(grad(66, 2), -0.13142082595504E-03, 1E-6);
  EXPECT_NEAR(grad(67, 0), 0.26589910462975E-03, 1E-6);
  EXPECT_NEAR(grad(67, 1), 0.29825079606157E-03, 1E-6);
  EXPECT_NEAR(grad(67, 2), -0.83844867563783E-06, 1E-6);
  EXPECT_NEAR(grad(68, 0), 0.19668146944048E-03, 1E-6);
  EXPECT_NEAR(grad(68, 1), 0.26369449476365E-03, 1E-6);
  EXPECT_NEAR(grad(68, 2), 0.23007932591959E-03, 1E-6);
  EXPECT_NEAR(grad(69, 0), 0.18054354511797E-03, 1E-6);
  EXPECT_NEAR(grad(69, 1), 0.30140694652346E-03, 1E-6);
  EXPECT_NEAR(grad(69, 2), 0.98524356642812E-04, 1E-6);
  EXPECT_NEAR(grad(70, 0), 0.15185299270785E-03, 1E-6);
  EXPECT_NEAR(grad(70, 1), 0.17686093322150E-03, 1E-6);
  EXPECT_NEAR(grad(70, 2), 0.38892532177029E-03, 1E-6);
  EXPECT_NEAR(grad(71, 0), 0.13118875300895E-03, 1E-6);
  EXPECT_NEAR(grad(71, 1), 0.24816912225637E-03, 1E-6);
  EXPECT_NEAR(grad(71, 2), 0.31474219006162E-03, 1E-6);
  EXPECT_NEAR(grad(72, 0), -0.15907200387641E-04, 1E-6);
  EXPECT_NEAR(grad(72, 1), 0.15775248385612E-03, 1E-6);
  EXPECT_NEAR(grad(72, 2), 0.46390635013370E-03, 1E-6);
  EXPECT_NEAR(grad(73, 0), 0.12410706381361E-03, 1E-6);
  EXPECT_NEAR(grad(73, 1), 0.27801961361957E-03, 1E-6);
  EXPECT_NEAR(grad(73, 2), -0.37956026266356E-03, 1E-6);
  EXPECT_NEAR(grad(74, 0), 0.15536870511466E-03, 1E-6);
  EXPECT_NEAR(grad(74, 1), 0.31081958097902E-03, 1E-6);
  EXPECT_NEAR(grad(74, 2), -0.22520009108406E-03, 1E-6);
  EXPECT_NEAR(grad(75, 0), 0.22553389367182E-03, 1E-6);
  EXPECT_NEAR(grad(75, 1), 0.28403908870781E-03, 1E-6);
  EXPECT_NEAR(grad(75, 2), -0.35359755813599E-03, 1E-6);
  EXPECT_NEAR(grad(76, 0), -0.26709033471120E-04, 1E-6);
  EXPECT_NEAR(grad(76, 1), 0.33156227434239E-04, 1E-6);
  EXPECT_NEAR(grad(76, 2), 0.50128430976318E-03, 1E-6);
  EXPECT_NEAR(grad(77, 0), -0.13543120685153E-04, 1E-6);
  EXPECT_NEAR(grad(77, 1), 0.14863476663695E-03, 1E-6);
  EXPECT_NEAR(grad(77, 2), 0.39296429101396E-03, 1E-6);
  EXPECT_NEAR(grad(78, 0), -0.89411760863278E-04, 1E-6);
  EXPECT_NEAR(grad(78, 1), 0.12904042028123E-03, 1E-6);
  EXPECT_NEAR(grad(78, 2), 0.49944343279991E-03, 1E-6);
  EXPECT_NEAR(grad(79, 0), -0.91312923565871E-04, 1E-6);
  EXPECT_NEAR(grad(79, 1), 0.31995395621708E-03, 1E-6);
  EXPECT_NEAR(grad(79, 2), 0.24146123718862E-03, 1E-6);
  EXPECT_NEAR(grad(80, 0), -0.38892178776850E-04, 1E-6);
  EXPECT_NEAR(grad(80, 1), 0.42422610641330E-03, 1E-6);
  EXPECT_NEAR(grad(80, 2), 0.44633624864199E-04, 1E-6);
  EXPECT_NEAR(grad(81, 0), -0.10146706587217E-04, 1E-6);
  EXPECT_NEAR(grad(81, 1), 0.34233497328127E-03, 1E-6);
  EXPECT_NEAR(grad(81, 2), 0.14537167440825E-03, 1E-6);
  EXPECT_NEAR(grad(82, 0), 0.83923850939089E-04, 1E-6);
  EXPECT_NEAR(grad(82, 1), 0.41254213764112E-03, 1E-6);
  EXPECT_NEAR(grad(82, 2), -0.13343409534023E-03, 1E-6);
  EXPECT_NEAR(grad(83, 0), 0.63429865879225E-04, 1E-6);
  EXPECT_NEAR(grad(83, 1), 0.42134421312684E-03, 1E-6);
  EXPECT_NEAR(grad(83, 2), -0.66066235703115E-04, 1E-6);
  EXPECT_NEAR(grad(84, 0), 0.21639173723381E-03, 1E-6);
  EXPECT_NEAR(grad(84, 1), 0.44425440716771E-03, 1E-6);
  EXPECT_NEAR(grad(84, 2), -0.17153760810181E-03, 1E-6);
  EXPECT_NEAR(grad(85, 0), -0.47239903952091E-03, 1E-6);
  EXPECT_NEAR(grad(85, 1), -0.10910430897498E-03, 1E-6);
  EXPECT_NEAR(grad(85, 2), -0.29337628381223E-05, 1E-6);
  EXPECT_NEAR(grad(86, 0), -0.35234931501245E-03, 1E-6);
  EXPECT_NEAR(grad(86, 1), -0.91829422436281E-04, 1E-6);
  EXPECT_NEAR(grad(86, 2), -0.66828818369931E-04, 1E-6);
  EXPECT_NEAR(grad(87, 0), -0.33727131045948E-03, 1E-6);
  EXPECT_NEAR(grad(87, 1), -0.23679821421978E-03, 1E-6);
  EXPECT_NEAR(grad(87, 2), 0.53051597627155E-04, 1E-6);
  EXPECT_NEAR(grad(88, 0), -0.29195427463737E-03, 1E-6);
  EXPECT_NEAR(grad(88, 1), 0.94017979959145E-05, 1E-6);
  EXPECT_NEAR(grad(88, 2), 0.35418239748718E-03, 1E-6);
  EXPECT_NEAR(grad(89, 0), -0.27738181806893E-03, 1E-6);
  EXPECT_NEAR(grad(89, 1), 0.11558385118043E-03, 1E-6);
  EXPECT_NEAR(grad(89, 2), 0.19105171156565E-03, 1E-6);
  EXPECT_NEAR(grad(90, 0), -0.26597241073913E-03, 1E-6);
  EXPECT_NEAR(grad(90, 1), 0.18602452177111E-03, 1E-6);
  EXPECT_NEAR(grad(90, 2), 0.26725018437871E-03, 1E-6);
  EXPECT_NEAR(grad(91, 0), -0.29291455746557E-03, 1E-6);
  EXPECT_NEAR(grad(91, 1), 0.41280482664878E-03, 1E-6);
  EXPECT_NEAR(grad(91, 2), -0.25224026199748E-03, 1E-6);
  EXPECT_NEAR(grad(92, 0), -0.17583426496709E-03, 1E-6);
  EXPECT_NEAR(grad(92, 1), 0.39209949154896E-03, 1E-6);
  EXPECT_NEAR(grad(92, 2), -0.16327253505943E-03, 1E-6);
  EXPECT_NEAR(grad(93, 0), -0.20555563092932E-03, 1E-6);
  EXPECT_NEAR(grad(93, 1), 0.26198237489625E-03, 1E-6);
  EXPECT_NEAR(grad(93, 2), -0.17548691016592E-03, 1E-6);
  EXPECT_NEAR(grad(94, 0), 0.61658483434154E-03, 1E-6);
  EXPECT_NEAR(grad(94, 1), 0.11514427144400E-03, 1E-6);
  EXPECT_NEAR(grad(94, 2), -0.74762955235450E-03, 1E-6);
  EXPECT_NEAR(grad(95, 0), 0.34934780555561E-03, 1E-6);
  EXPECT_NEAR(grad(95, 1), 0.61832486905442E-03, 1E-6);
  EXPECT_NEAR(grad(95, 2), -0.14419018569281E-02, 1E-6);
  EXPECT_NEAR(grad(96, 0), 0.52233964664844E-03, 1E-6);
  EXPECT_NEAR(grad(96, 1), 0.83293342123291E-03, 1E-6);
  EXPECT_NEAR(grad(96, 2), -0.14946783779889E-02, 1E-6);
  EXPECT_NEAR(grad(97, 0), 0.27033204654944E-04, 1E-6);
  EXPECT_NEAR(grad(97, 1), -0.32648421936354E-04, 1E-6);
  EXPECT_NEAR(grad(97, 2), -0.13636735773680E-02, 1E-6);
  EXPECT_NEAR(grad(98, 0), 0.12670444968729E-03, 1E-6);
  EXPECT_NEAR(grad(98, 1), 0.28787781998577E-03, 1E-6);
  EXPECT_NEAR(grad(98, 2), -0.42574282950018E-03, 1E-6);
  EXPECT_NEAR(grad(99, 0), 0.26148472507445E-03, 1E-6);
  EXPECT_NEAR(grad(99, 1), 0.27654293893619E-03, 1E-6);
  EXPECT_NEAR(grad(99, 2), -0.52847207556216E-03, 1E-6);
  EXPECT_NEAR(grad(100, 0), 0.17220965745413E-03, 1E-6);
  EXPECT_NEAR(grad(100, 1), 0.30755723619206E-03, 1E-6);
  EXPECT_NEAR(grad(100, 2), -0.46063939191557E-03, 1E-6);
  EXPECT_NEAR(grad(101, 0), -0.25633475152967E-03, 1E-6);
  EXPECT_NEAR(grad(101, 1), -0.56647960951220E-04, 1E-6);
  EXPECT_NEAR(grad(101, 2), -0.14994478117036E-02, 1E-6);
  EXPECT_NEAR(grad(102, 0), -0.49308591589684E-04, 1E-6);
  EXPECT_NEAR(grad(102, 1), -0.29400572566996E-03, 1E-6);
  EXPECT_NEAR(grad(102, 2), -0.10926886700604E-02, 1E-6);
  EXPECT_NEAR(grad(103, 0), -0.31729832680110E-03, 1E-6);
  EXPECT_NEAR(grad(103, 1), -0.25500474248704E-03, 1E-6);
  EXPECT_NEAR(grad(103, 2), -0.11837421975995E-02, 1E-6);
  EXPECT_NEAR(grad(104, 0), 0.19690199139478E-03, 1E-6);
  EXPECT_NEAR(grad(104, 1), 0.36108730700931E-05, 1E-6);
  EXPECT_NEAR(grad(104, 2), -0.55156486417443E-03, 1E-6);
  EXPECT_NEAR(grad(105, 0), -0.24913241272391E-03, 1E-6);
  EXPECT_NEAR(grad(105, 1), -0.20418058626915E-03, 1E-6);
  EXPECT_NEAR(grad(105, 2), -0.62158542390943E-03, 1E-6);
  EXPECT_NEAR(grad(106, 0), -0.16343865228969E-05, 1E-6);
  EXPECT_NEAR(grad(106, 1), -0.16049067415988E-03, 1E-6);
  EXPECT_NEAR(grad(106, 2), -0.44008392762235E-03, 1E-6);
  EXPECT_NEAR(grad(107, 0), -0.25061858982336E-03, 1E-6);
  EXPECT_NEAR(grad(107, 1), -0.34655531597798E-03, 1E-6);
  EXPECT_NEAR(grad(107, 2), -0.66916344366619E-03, 1E-6);
  EXPECT_NEAR(grad(108, 0), -0.16787714792676E-03, 1E-6);
  EXPECT_NEAR(grad(108, 1), 0.56172804374167E-04, 1E-6);
  EXPECT_NEAR(grad(108, 2), -0.42315188768987E-03, 1E-6);
  EXPECT_NEAR(grad(109, 0), -0.13468544259881E-03, 1E-6);
  EXPECT_NEAR(grad(109, 1), -0.21305845759471E-03, 1E-6);
  EXPECT_NEAR(grad(109, 2), -0.15482907190840E-03, 1E-6);
  EXPECT_NEAR(grad(110, 0), -0.23573775182983E-03, 1E-6);
  EXPECT_NEAR(grad(110, 1), -0.94354164089405E-04, 1E-6);
  EXPECT_NEAR(grad(110, 2), -0.18100341338633E-03, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ gradient correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJGradBP86JacobsenGhost) {
  auto geometry = *systemControllerTwo->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJ,
                                                                               std::make_shared<Geometry>(geometry),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 0), -0.44641271320771E-03, 1E-6);
  EXPECT_NEAR(grad(0, 1), -0.34836653878647E-03, 1E-6);
  EXPECT_NEAR(grad(0, 2), 0.23794535455227E-03, 1E-6);
  EXPECT_NEAR(grad(1, 0), 0.45577469855532E-03, 1E-6);
  EXPECT_NEAR(grad(1, 1), -0.51216118328020E-03, 1E-6);
  EXPECT_NEAR(grad(1, 2), -0.22358523681505E-03, 1E-6);
  EXPECT_NEAR(grad(2, 0), -0.16787504556191E-03, 1E-6);
  EXPECT_NEAR(grad(2, 1), -0.18620322081428E-03, 1E-6);
  EXPECT_NEAR(grad(2, 2), 0.34017282987155E-03, 1E-6);
  EXPECT_NEAR(grad(3, 0), -0.60630814302006E-03, 1E-6);
  EXPECT_NEAR(grad(3, 1), -0.74541811578157E-03, 1E-6);
  EXPECT_NEAR(grad(3, 2), 0.16413329853933E-03, 1E-6);
  EXPECT_NEAR(grad(4, 0), -0.26284075145145E-03, 1E-6);
  EXPECT_NEAR(grad(4, 1), -0.89624101639072E-04, 1E-6);
  EXPECT_NEAR(grad(4, 2), 0.40298169397188E-03, 1E-6);
  EXPECT_NEAR(grad(5, 0), 0.33709084661949E-03, 1E-6);
  EXPECT_NEAR(grad(5, 1), -0.81644805432351E-03, 1E-6);
  EXPECT_NEAR(grad(5, 2), -0.18345467497986E-03, 1E-6);
  EXPECT_NEAR(grad(6, 0), 0.43843889981468E-03, 1E-6);
  EXPECT_NEAR(grad(6, 1), -0.16674409600683E-03, 1E-6);
  EXPECT_NEAR(grad(6, 2), -0.13161594754476E-03, 1E-6);
  EXPECT_NEAR(grad(7, 0), 0.24192941011567E-03, 1E-6);
  EXPECT_NEAR(grad(7, 1), -0.16832764022001E-03, 1E-6);
  EXPECT_NEAR(grad(7, 2), -0.30474423492655E-03, 1E-6);
  EXPECT_NEAR(grad(8, 0), -0.34683727112197E-03, 1E-6);
  EXPECT_NEAR(grad(8, 1), -0.71473634375525E-03, 1E-6);
  EXPECT_NEAR(grad(8, 2), 0.23349652243208E-04, 1E-6);
  EXPECT_NEAR(grad(9, 0), -0.25045716376056E-03, 1E-6);
  EXPECT_NEAR(grad(9, 1), -0.35481675979231E-03, 1E-6);
  EXPECT_NEAR(grad(9, 2), -0.47783472255487E-04, 1E-6);
  EXPECT_NEAR(grad(10, 0), -0.35809598349183E-03, 1E-6);
  EXPECT_NEAR(grad(10, 1), -0.24140359149052E-03, 1E-6);
  EXPECT_NEAR(grad(10, 2), 0.12909347585806E-03, 1E-6);
  EXPECT_NEAR(grad(11, 0), 0.24560285483310E-04, 1E-6);
  EXPECT_NEAR(grad(11, 1), -0.74514937647724E-03, 1E-6);
  EXPECT_NEAR(grad(11, 2), -0.21485206615971E-03, 1E-6);
  EXPECT_NEAR(grad(12, 0), -0.20367754473177E-03, 1E-6);
  EXPECT_NEAR(grad(12, 1), -0.26269893634963E-03, 1E-6);
  EXPECT_NEAR(grad(12, 2), -0.24559082119115E-04, 1E-6);
  EXPECT_NEAR(grad(13, 0), -0.18259474524204E-03, 1E-6);
  EXPECT_NEAR(grad(13, 1), -0.30353968247095E-03, 1E-6);
  EXPECT_NEAR(grad(13, 2), 0.14374385178206E-03, 1E-6);
  EXPECT_NEAR(grad(14, 0), 0.57658195466981E-04, 1E-6);
  EXPECT_NEAR(grad(14, 1), -0.31249288904384E-03, 1E-6);
  EXPECT_NEAR(grad(14, 2), -0.23357036167616E-03, 1E-6);
  EXPECT_NEAR(grad(15, 0), 0.56558654485810E-04, 1E-6);
  EXPECT_NEAR(grad(15, 1), -0.31072153068333E-03, 1E-6);
  EXPECT_NEAR(grad(15, 2), -0.67688361000117E-04, 1E-6);
  EXPECT_NEAR(grad(16, 0), 0.17678020442675E-03, 1E-6);
  EXPECT_NEAR(grad(16, 1), -0.40706000522237E-03, 1E-6);
  EXPECT_NEAR(grad(16, 2), 0.11029162574860E-03, 1E-6);
  EXPECT_NEAR(grad(17, 0), 0.22897539601370E-03, 1E-6);
  EXPECT_NEAR(grad(17, 1), -0.30429747189709E-03, 1E-6);
  EXPECT_NEAR(grad(17, 2), -0.14996223947871E-03, 1E-6);
  EXPECT_NEAR(grad(18, 0), 0.73903797560823E-03, 1E-6);
  EXPECT_NEAR(grad(18, 1), -0.24943995399781E-03, 1E-6);
  EXPECT_NEAR(grad(18, 2), -0.78063010908874E-04, 1E-6);
  EXPECT_NEAR(grad(19, 0), -0.61740474005408E-03, 1E-6);
  EXPECT_NEAR(grad(19, 1), -0.28919030071639E-05, 1E-6);
  EXPECT_NEAR(grad(19, 2), 0.81457294885140E-03, 1E-6);
  EXPECT_NEAR(grad(20, 0), 0.34417247784298E-03, 1E-6);
  EXPECT_NEAR(grad(20, 1), -0.61542332660393E-03, 1E-6);
  EXPECT_NEAR(grad(20, 2), 0.27938015674633E-03, 1E-6);
  EXPECT_NEAR(grad(21, 0), -0.37474625438297E-03, 1E-6);
  EXPECT_NEAR(grad(21, 1), -0.77615262525778E-03, 1E-6);
  EXPECT_NEAR(grad(21, 2), 0.25564030684309E-02, 1E-6);
  EXPECT_NEAR(grad(22, 0), 0.28156240010851E-03, 1E-6);
  EXPECT_NEAR(grad(22, 1), -0.55492643569041E-04, 1E-6);
  EXPECT_NEAR(grad(22, 2), -0.13014041559225E-03, 1E-6);
  EXPECT_NEAR(grad(23, 0), -0.13103105018074E-03, 1E-6);
  EXPECT_NEAR(grad(23, 1), 0.27799316518425E-04, 1E-6);
  EXPECT_NEAR(grad(23, 2), 0.15646036491417E-03, 1E-6);
  EXPECT_NEAR(grad(24, 0), -0.89193053184762E-03, 1E-6);
  EXPECT_NEAR(grad(24, 1), -0.15967087536234E-03, 1E-6);
  EXPECT_NEAR(grad(24, 2), 0.33601158334330E-03, 1E-6);
  EXPECT_NEAR(grad(25, 0), -0.56560909635577E-03, 1E-6);
  EXPECT_NEAR(grad(25, 1), 0.21152218527148E-03, 1E-6);
  EXPECT_NEAR(grad(25, 2), 0.97275710042093E-03, 1E-6);
  EXPECT_NEAR(grad(26, 0), -0.18379499207937E-03, 1E-6);
  EXPECT_NEAR(grad(26, 1), 0.18940050876816E-03, 1E-6);
  EXPECT_NEAR(grad(26, 2), 0.55814736958811E-03, 1E-6);
  EXPECT_NEAR(grad(27, 0), -0.28119703066395E-03, 1E-6);
  EXPECT_NEAR(grad(27, 1), -0.16966670570161E-03, 1E-6);
  EXPECT_NEAR(grad(27, 2), 0.30675614303623E-03, 1E-6);
  EXPECT_NEAR(grad(28, 0), 0.94087469680944E-04, 1E-6);
  EXPECT_NEAR(grad(28, 1), -0.16441875141781E-03, 1E-6);
  EXPECT_NEAR(grad(28, 2), 0.11479233636263E-02, 1E-6);
  EXPECT_NEAR(grad(29, 0), 0.62173835576336E-03, 1E-6);
  EXPECT_NEAR(grad(29, 1), -0.81389700838245E-04, 1E-6);
  EXPECT_NEAR(grad(29, 2), -0.22821199827511E-04, 1E-6);
  EXPECT_NEAR(grad(30, 0), 0.21843889733165E-03, 1E-6);
  EXPECT_NEAR(grad(30, 1), -0.22654742282420E-03, 1E-6);
  EXPECT_NEAR(grad(30, 2), -0.56436923651591E-04, 1E-6);
  EXPECT_NEAR(grad(31, 0), 0.68903716509116E-03, 1E-6);
  EXPECT_NEAR(grad(31, 1), -0.58973713225538E-03, 1E-6);
  EXPECT_NEAR(grad(31, 2), -0.72085095762718E-04, 1E-6);
  EXPECT_NEAR(grad(32, 0), 0.23521805385171E-03, 1E-6);
  EXPECT_NEAR(grad(32, 1), 0.91434535733175E-04, 1E-6);
  EXPECT_NEAR(grad(32, 2), -0.27059342000401E-04, 1E-6);
  EXPECT_NEAR(grad(33, 0), -0.28263706522879E-03, 1E-6);
  EXPECT_NEAR(grad(33, 1), 0.25238722703096E-04, 1E-6);
  EXPECT_NEAR(grad(33, 2), -0.17971438080203E-03, 1E-6);
  EXPECT_NEAR(grad(34, 0), -0.14349099763395E-03, 1E-6);
  EXPECT_NEAR(grad(34, 1), -0.36734470046452E-03, 1E-6);
  EXPECT_NEAR(grad(34, 2), 0.33238504618144E-03, 1E-6);
  EXPECT_NEAR(grad(35, 0), -0.39969396352155E-03, 1E-6);
  EXPECT_NEAR(grad(35, 1), 0.35905083869383E-04, 1E-6);
  EXPECT_NEAR(grad(35, 2), 0.61785762546793E-03, 1E-6);
  EXPECT_NEAR(grad(36, 0), -0.26006224750122E-03, 1E-6);
  EXPECT_NEAR(grad(36, 1), -0.21238714177715E-03, 1E-6);
  EXPECT_NEAR(grad(36, 2), 0.33696349787377E-03, 1E-6);
  EXPECT_NEAR(grad(37, 0), -0.11569947728813E-03, 1E-6);
  EXPECT_NEAR(grad(37, 1), 0.53095339335688E-03, 1E-6);
  EXPECT_NEAR(grad(37, 2), 0.33858910927632E-03, 1E-6);
  EXPECT_NEAR(grad(38, 0), 0.88463938561081E-04, 1E-6);
  EXPECT_NEAR(grad(38, 1), -0.15787530828135E-03, 1E-6);
  EXPECT_NEAR(grad(38, 2), -0.19161632052436E-04, 1E-6);
  EXPECT_NEAR(grad(39, 0), -0.21442293573086E-05, 1E-6);
  EXPECT_NEAR(grad(39, 1), -0.25193569212640E-03, 1E-6);
  EXPECT_NEAR(grad(39, 2), -0.57739268764463E-04, 1E-6);
  EXPECT_NEAR(grad(40, 0), -0.16209776098162E-03, 1E-6);
  EXPECT_NEAR(grad(40, 1), -0.34364021238956E-03, 1E-6);
  EXPECT_NEAR(grad(40, 2), -0.14763058402152E-03, 1E-6);
  EXPECT_NEAR(grad(41, 0), 0.36577755480510E-03, 1E-6);
  EXPECT_NEAR(grad(41, 1), 0.39471668527540E-03, 1E-6);
  EXPECT_NEAR(grad(41, 2), 0.99053757893348E-04, 1E-6);
  EXPECT_NEAR(grad(42, 0), -0.54886977155488E-03, 1E-6);
  EXPECT_NEAR(grad(42, 1), 0.17715312675855E-03, 1E-6);
  EXPECT_NEAR(grad(42, 2), 0.20499676161472E-03, 1E-6);
  EXPECT_NEAR(grad(43, 0), -0.88552448494570E-04, 1E-6);
  EXPECT_NEAR(grad(43, 1), 0.24207094059259E-03, 1E-6);
  EXPECT_NEAR(grad(43, 2), 0.11389353518825E-03, 1E-6);
  EXPECT_NEAR(grad(44, 0), 0.34626618506442E-03, 1E-6);
  EXPECT_NEAR(grad(44, 1), -0.16659339344682E-03, 1E-6);
  EXPECT_NEAR(grad(44, 2), 0.86740299578476E-04, 1E-6);
  EXPECT_NEAR(grad(45, 0), 0.16367720529031E-03, 1E-6);
  EXPECT_NEAR(grad(45, 1), 0.14800220164615E-03, 1E-6);
  EXPECT_NEAR(grad(45, 2), 0.72661276468407E-04, 1E-6);
  EXPECT_NEAR(grad(46, 0), -0.21447756010488E-03, 1E-6);
  EXPECT_NEAR(grad(46, 1), 0.38814549380924E-03, 1E-6);
  EXPECT_NEAR(grad(46, 2), 0.76811641717841E-03, 1E-6);
  EXPECT_NEAR(grad(47, 0), -0.43761247618718E-04, 1E-6);
  EXPECT_NEAR(grad(47, 1), 0.61632870436173E-03, 1E-6);
  EXPECT_NEAR(grad(47, 2), 0.27683860053273E-03, 1E-6);
  EXPECT_NEAR(grad(48, 0), 0.15098279301011E-04, 1E-6);
  EXPECT_NEAR(grad(48, 1), 0.93557160546381E-03, 1E-6);
  EXPECT_NEAR(grad(48, 2), -0.33050379082729E-04, 1E-6);
  EXPECT_NEAR(grad(49, 0), -0.74596918775100E-03, 1E-6);
  EXPECT_NEAR(grad(49, 1), -0.15043852117657E-03, 1E-6);
  EXPECT_NEAR(grad(49, 2), 0.13005920216928E-03, 1E-6);
  EXPECT_NEAR(grad(50, 0), -0.56163311595122E-03, 1E-6);
  EXPECT_NEAR(grad(50, 1), 0.62036182589626E-03, 1E-6);
  EXPECT_NEAR(grad(50, 2), -0.15904402086311E-03, 1E-6);
  EXPECT_NEAR(grad(51, 0), -0.52961861065490E-03, 1E-6);
  EXPECT_NEAR(grad(51, 1), 0.18542169563333E-03, 1E-6);
  EXPECT_NEAR(grad(51, 2), 0.37920803427412E-03, 1E-6);
  EXPECT_NEAR(grad(52, 0), 0.45461184156401E-03, 1E-6);
  EXPECT_NEAR(grad(52, 1), 0.55633929051552E-03, 1E-6);
  EXPECT_NEAR(grad(52, 2), -0.44506619795246E-03, 1E-6);
  EXPECT_NEAR(grad(53, 0), 0.33787323145827E-03, 1E-6);
  EXPECT_NEAR(grad(53, 1), 0.51472539709650E-03, 1E-6);
  EXPECT_NEAR(grad(53, 2), 0.17274932015987E-03, 1E-6);
  EXPECT_NEAR(grad(54, 0), 0.45860007196877E-03, 1E-6);
  EXPECT_NEAR(grad(54, 1), 0.41609792638722E-03, 1E-6);
  EXPECT_NEAR(grad(54, 2), 0.66520654768515E-03, 1E-6);
  EXPECT_NEAR(grad(55, 0), 0.26602599269625E-03, 1E-6);
  EXPECT_NEAR(grad(55, 1), -0.51313466515424E-03, 1E-6);
  EXPECT_NEAR(grad(55, 2), 0.60562121870499E-05, 1E-6);
  EXPECT_NEAR(grad(56, 0), 0.33776188964159E-03, 1E-6);
  EXPECT_NEAR(grad(56, 1), -0.77598266465638E-04, 1E-6);
  EXPECT_NEAR(grad(56, 2), 0.54284619083235E-03, 1E-6);
  EXPECT_NEAR(grad(57, 0), 0.59816288304327E-03, 1E-6);
  EXPECT_NEAR(grad(57, 1), 0.50024928219000E-04, 1E-6);
  EXPECT_NEAR(grad(57, 2), -0.23037303670062E-03, 1E-6);
  EXPECT_NEAR(grad(58, 0), 0.32585862194356E-03, 1E-6);
  EXPECT_NEAR(grad(58, 1), 0.38022784558631E-04, 1E-6);
  EXPECT_NEAR(grad(58, 2), -0.29996961476651E-03, 1E-6);
  EXPECT_NEAR(grad(59, 0), 0.32601971990793E-03, 1E-6);
  EXPECT_NEAR(grad(59, 1), 0.39677021961920E-04, 1E-6);
  EXPECT_NEAR(grad(59, 2), -0.10822338015477E-03, 1E-6);
  EXPECT_NEAR(grad(60, 0), 0.35879238987349E-03, 1E-6);
  EXPECT_NEAR(grad(60, 1), 0.13551136441295E-03, 1E-6);
  EXPECT_NEAR(grad(60, 2), -0.14024983294786E-03, 1E-6);
  EXPECT_NEAR(grad(61, 0), 0.78913900690811E-04, 1E-6);
  EXPECT_NEAR(grad(61, 1), -0.82239511390362E-04, 1E-6);
  EXPECT_NEAR(grad(61, 2), 0.43168740640081E-03, 1E-6);
  EXPECT_NEAR(grad(62, 0), 0.19101918722848E-03, 1E-6);
  EXPECT_NEAR(grad(62, 1), 0.54270232450281E-04, 1E-6);
  EXPECT_NEAR(grad(62, 2), 0.35754976768672E-03, 1E-6);
  EXPECT_NEAR(grad(63, 0), 0.19024119042316E-03, 1E-6);
  EXPECT_NEAR(grad(63, 1), -0.25233440737570E-04, 1E-6);
  EXPECT_NEAR(grad(63, 2), 0.28806297248064E-03, 1E-6);
  EXPECT_NEAR(grad(64, 0), 0.60264802794365E-04, 1E-6);
  EXPECT_NEAR(grad(64, 1), -0.36224162920476E-03, 1E-6);
  EXPECT_NEAR(grad(64, 2), 0.84571549113634E-04, 1E-6);
  EXPECT_NEAR(grad(65, 0), 0.18434142585794E-03, 1E-6);
  EXPECT_NEAR(grad(65, 1), -0.28307699026813E-03, 1E-6);
  EXPECT_NEAR(grad(65, 2), 0.17970762842745E-04, 1E-6);
  EXPECT_NEAR(grad(66, 0), 0.13334669353749E-03, 1E-6);
  EXPECT_NEAR(grad(66, 1), -0.32579829125995E-03, 1E-6);
  EXPECT_NEAR(grad(66, 2), -0.13142082595504E-03, 1E-6);
  EXPECT_NEAR(grad(67, 0), 0.26589910462975E-03, 1E-6);
  EXPECT_NEAR(grad(67, 1), 0.29825079606157E-03, 1E-6);
  EXPECT_NEAR(grad(67, 2), -0.83844867563783E-06, 1E-6);
  EXPECT_NEAR(grad(68, 0), 0.19668146944048E-03, 1E-6);
  EXPECT_NEAR(grad(68, 1), 0.26369449476365E-03, 1E-6);
  EXPECT_NEAR(grad(68, 2), 0.23007932591959E-03, 1E-6);
  EXPECT_NEAR(grad(69, 0), 0.18054354511797E-03, 1E-6);
  EXPECT_NEAR(grad(69, 1), 0.30140694652346E-03, 1E-6);
  EXPECT_NEAR(grad(69, 2), 0.98524356642812E-04, 1E-6);
  EXPECT_NEAR(grad(70, 0), 0.15185299270785E-03, 1E-6);
  EXPECT_NEAR(grad(70, 1), 0.17686093322150E-03, 1E-6);
  EXPECT_NEAR(grad(70, 2), 0.38892532177029E-03, 1E-6);
  EXPECT_NEAR(grad(71, 0), 0.13118875300895E-03, 1E-6);
  EXPECT_NEAR(grad(71, 1), 0.24816912225637E-03, 1E-6);
  EXPECT_NEAR(grad(71, 2), 0.31474219006162E-03, 1E-6);
  EXPECT_NEAR(grad(72, 0), -0.15907200387641E-04, 1E-6);
  EXPECT_NEAR(grad(72, 1), 0.15775248385612E-03, 1E-6);
  EXPECT_NEAR(grad(72, 2), 0.46390635013370E-03, 1E-6);
  EXPECT_NEAR(grad(73, 0), 0.12410706381361E-03, 1E-6);
  EXPECT_NEAR(grad(73, 1), 0.27801961361957E-03, 1E-6);
  EXPECT_NEAR(grad(73, 2), -0.37956026266356E-03, 1E-6);
  EXPECT_NEAR(grad(74, 0), 0.15536870511466E-03, 1E-6);
  EXPECT_NEAR(grad(74, 1), 0.31081958097902E-03, 1E-6);
  EXPECT_NEAR(grad(74, 2), -0.22520009108406E-03, 1E-6);
  EXPECT_NEAR(grad(75, 0), 0.22553389367182E-03, 1E-6);
  EXPECT_NEAR(grad(75, 1), 0.28403908870781E-03, 1E-6);
  EXPECT_NEAR(grad(75, 2), -0.35359755813599E-03, 1E-6);
  EXPECT_NEAR(grad(76, 0), -0.26709033471120E-04, 1E-6);
  EXPECT_NEAR(grad(76, 1), 0.33156227434239E-04, 1E-6);
  EXPECT_NEAR(grad(76, 2), 0.50128430976318E-03, 1E-6);
  EXPECT_NEAR(grad(77, 0), -0.13543120685153E-04, 1E-6);
  EXPECT_NEAR(grad(77, 1), 0.14863476663695E-03, 1E-6);
  EXPECT_NEAR(grad(77, 2), 0.39296429101396E-03, 1E-6);
  EXPECT_NEAR(grad(78, 0), -0.89411760863278E-04, 1E-6);
  EXPECT_NEAR(grad(78, 1), 0.12904042028123E-03, 1E-6);
  EXPECT_NEAR(grad(78, 2), 0.49944343279991E-03, 1E-6);
  EXPECT_NEAR(grad(79, 0), -0.91312923565871E-04, 1E-6);
  EXPECT_NEAR(grad(79, 1), 0.31995395621708E-03, 1E-6);
  EXPECT_NEAR(grad(79, 2), 0.24146123718862E-03, 1E-6);
  EXPECT_NEAR(grad(80, 0), -0.38892178776850E-04, 1E-6);
  EXPECT_NEAR(grad(80, 1), 0.42422610641330E-03, 1E-6);
  EXPECT_NEAR(grad(80, 2), 0.44633624864199E-04, 1E-6);
  EXPECT_NEAR(grad(81, 0), -0.10146706587217E-04, 1E-6);
  EXPECT_NEAR(grad(81, 1), 0.34233497328127E-03, 1E-6);
  EXPECT_NEAR(grad(81, 2), 0.14537167440825E-03, 1E-6);
  EXPECT_NEAR(grad(82, 0), 0.83923850939089E-04, 1E-6);
  EXPECT_NEAR(grad(82, 1), 0.41254213764112E-03, 1E-6);
  EXPECT_NEAR(grad(82, 2), -0.13343409534023E-03, 1E-6);
  EXPECT_NEAR(grad(83, 0), 0.63429865879225E-04, 1E-6);
  EXPECT_NEAR(grad(83, 1), 0.42134421312684E-03, 1E-6);
  EXPECT_NEAR(grad(83, 2), -0.66066235703115E-04, 1E-6);
  EXPECT_NEAR(grad(84, 0), 0.21639173723381E-03, 1E-6);
  EXPECT_NEAR(grad(84, 1), 0.44425440716771E-03, 1E-6);
  EXPECT_NEAR(grad(84, 2), -0.17153760810181E-03, 1E-6);
  EXPECT_NEAR(grad(85, 0), -0.47239903952091E-03, 1E-6);
  EXPECT_NEAR(grad(85, 1), -0.10910430897498E-03, 1E-6);
  EXPECT_NEAR(grad(85, 2), -0.29337628381223E-05, 1E-6);
  EXPECT_NEAR(grad(86, 0), -0.35234931501245E-03, 1E-6);
  EXPECT_NEAR(grad(86, 1), -0.91829422436281E-04, 1E-6);
  EXPECT_NEAR(grad(86, 2), -0.66828818369931E-04, 1E-6);
  EXPECT_NEAR(grad(87, 0), -0.33727131045948E-03, 1E-6);
  EXPECT_NEAR(grad(87, 1), -0.23679821421978E-03, 1E-6);
  EXPECT_NEAR(grad(87, 2), 0.53051597627155E-04, 1E-6);
  EXPECT_NEAR(grad(88, 0), -0.29195427463737E-03, 1E-6);
  EXPECT_NEAR(grad(88, 1), 0.94017979959145E-05, 1E-6);
  EXPECT_NEAR(grad(88, 2), 0.35418239748718E-03, 1E-6);
  EXPECT_NEAR(grad(89, 0), -0.27738181806893E-03, 1E-6);
  EXPECT_NEAR(grad(89, 1), 0.11558385118043E-03, 1E-6);
  EXPECT_NEAR(grad(89, 2), 0.19105171156565E-03, 1E-6);
  EXPECT_NEAR(grad(90, 0), -0.26597241073913E-03, 1E-6);
  EXPECT_NEAR(grad(90, 1), 0.18602452177111E-03, 1E-6);
  EXPECT_NEAR(grad(90, 2), 0.26725018437871E-03, 1E-6);
  EXPECT_NEAR(grad(91, 0), -0.29291455746557E-03, 1E-6);
  EXPECT_NEAR(grad(91, 1), 0.41280482664878E-03, 1E-6);
  EXPECT_NEAR(grad(91, 2), -0.25224026199748E-03, 1E-6);
  EXPECT_NEAR(grad(92, 0), -0.17583426496709E-03, 1E-6);
  EXPECT_NEAR(grad(92, 1), 0.39209949154896E-03, 1E-6);
  EXPECT_NEAR(grad(92, 2), -0.16327253505943E-03, 1E-6);
  EXPECT_NEAR(grad(93, 0), -0.20555563092932E-03, 1E-6);
  EXPECT_NEAR(grad(93, 1), 0.26198237489625E-03, 1E-6);
  EXPECT_NEAR(grad(93, 2), -0.17548691016592E-03, 1E-6);
  EXPECT_NEAR(grad(94, 0), 0.61658483434154E-03, 1E-6);
  EXPECT_NEAR(grad(94, 1), 0.11514427144400E-03, 1E-6);
  EXPECT_NEAR(grad(94, 2), -0.74762955235450E-03, 1E-6);
  EXPECT_NEAR(grad(95, 0), 0.34934780555561E-03, 1E-6);
  EXPECT_NEAR(grad(95, 1), 0.61832486905442E-03, 1E-6);
  EXPECT_NEAR(grad(95, 2), -0.14419018569281E-02, 1E-6);
  EXPECT_NEAR(grad(96, 0), 0.52233964664844E-03, 1E-6);
  EXPECT_NEAR(grad(96, 1), 0.83293342123291E-03, 1E-6);
  EXPECT_NEAR(grad(96, 2), -0.14946783779889E-02, 1E-6);
  EXPECT_NEAR(grad(97, 0), 0.27033204654944E-04, 1E-6);
  EXPECT_NEAR(grad(97, 1), -0.32648421936354E-04, 1E-6);
  EXPECT_NEAR(grad(97, 2), -0.13636735773680E-02, 1E-6);
  EXPECT_NEAR(grad(98, 0), 0.12670444968729E-03, 1E-6);
  EXPECT_NEAR(grad(98, 1), 0.28787781998577E-03, 1E-6);
  EXPECT_NEAR(grad(98, 2), -0.42574282950018E-03, 1E-6);
  EXPECT_NEAR(grad(99, 0), 0.26148472507445E-03, 1E-6);
  EXPECT_NEAR(grad(99, 1), 0.27654293893619E-03, 1E-6);
  EXPECT_NEAR(grad(99, 2), -0.52847207556216E-03, 1E-6);
  EXPECT_NEAR(grad(100, 0), 0.17220965745413E-03, 1E-6);
  EXPECT_NEAR(grad(100, 1), 0.30755723619206E-03, 1E-6);
  EXPECT_NEAR(grad(100, 2), -0.46063939191557E-03, 1E-6);
  EXPECT_NEAR(grad(101, 0), -0.25633475152967E-03, 1E-6);
  EXPECT_NEAR(grad(101, 1), -0.56647960951220E-04, 1E-6);
  EXPECT_NEAR(grad(101, 2), -0.14994478117036E-02, 1E-6);
  EXPECT_NEAR(grad(102, 0), -0.49308591589684E-04, 1E-6);
  EXPECT_NEAR(grad(102, 1), -0.29400572566996E-03, 1E-6);
  EXPECT_NEAR(grad(102, 2), -0.10926886700604E-02, 1E-6);
  EXPECT_NEAR(grad(103, 0), -0.31729832680110E-03, 1E-6);
  EXPECT_NEAR(grad(103, 1), -0.25500474248704E-03, 1E-6);
  EXPECT_NEAR(grad(103, 2), -0.11837421975995E-02, 1E-6);
  EXPECT_NEAR(grad(104, 0), 0.19690199139478E-03, 1E-6);
  EXPECT_NEAR(grad(104, 1), 0.36108730700931E-05, 1E-6);
  EXPECT_NEAR(grad(104, 2), -0.55156486417443E-03, 1E-6);
  EXPECT_NEAR(grad(105, 0), -0.24913241272391E-03, 1E-6);
  EXPECT_NEAR(grad(105, 1), -0.20418058626915E-03, 1E-6);
  EXPECT_NEAR(grad(105, 2), -0.62158542390943E-03, 1E-6);
  EXPECT_NEAR(grad(106, 0), -0.16343865228969E-05, 1E-6);
  EXPECT_NEAR(grad(106, 1), -0.16049067415988E-03, 1E-6);
  EXPECT_NEAR(grad(106, 2), -0.44008392762235E-03, 1E-6);
  EXPECT_NEAR(grad(107, 0), -0.25061858982336E-03, 1E-6);
  EXPECT_NEAR(grad(107, 1), -0.34655531597798E-03, 1E-6);
  EXPECT_NEAR(grad(107, 2), -0.66916344366619E-03, 1E-6);
  EXPECT_NEAR(grad(108, 0), -0.16787714792676E-03, 1E-6);
  EXPECT_NEAR(grad(108, 1), 0.56172804374167E-04, 1E-6);
  EXPECT_NEAR(grad(108, 2), -0.42315188768987E-03, 1E-6);
  EXPECT_NEAR(grad(109, 0), -0.13468544259881E-03, 1E-6);
  EXPECT_NEAR(grad(109, 1), -0.21305845759471E-03, 1E-6);
  EXPECT_NEAR(grad(109, 2), -0.15482907190840E-03, 1E-6);
  EXPECT_NEAR(grad(110, 0), -0.23573775182983E-03, 1E-6);
  EXPECT_NEAR(grad(110, 1), -0.94354164089405E-04, 1E-6);
  EXPECT_NEAR(grad(110, 2), -0.18100341338633E-03, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ-ABC gradient correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJABCGradBP86CO) {
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJABC,
                                                                               systemControllerOne->getGeometry(),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 2), -0.73442092899124E-6, 1E-8);
  EXPECT_NEAR(grad(1, 2), 0.73442092899124E-6, 1E-8);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC gradient correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJABCGradBP86Jacobsen) {
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJABC,
                                                                               systemControllerTwo->getGeometry(),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 0), -0.44419424846835E-03, 1E-6);
  EXPECT_NEAR(grad(0, 1), -0.34497327487658E-03, 1E-6);
  EXPECT_NEAR(grad(0, 2), 0.25751563332388E-03, 1E-6);
  EXPECT_NEAR(grad(1, 0), 0.44810919790884E-03, 1E-6);
  EXPECT_NEAR(grad(1, 1), -0.49797428237543E-03, 1E-6);
  EXPECT_NEAR(grad(1, 2), -0.24012114425014E-03, 1E-6);
  EXPECT_NEAR(grad(2, 0), -0.17523503799171E-03, 1E-6);
  EXPECT_NEAR(grad(2, 1), -0.17538269971019E-03, 1E-6);
  EXPECT_NEAR(grad(2, 2), 0.34753130839096E-03, 1E-6);
  EXPECT_NEAR(grad(3, 0), -0.62803624147614E-03, 1E-6);
  EXPECT_NEAR(grad(3, 1), -0.73165514136507E-03, 1E-6);
  EXPECT_NEAR(grad(3, 2), 0.15753938993305E-03, 1E-6);
  EXPECT_NEAR(grad(4, 0), -0.26356763387021E-03, 1E-6);
  EXPECT_NEAR(grad(4, 1), -0.78137545039569E-04, 1E-6);
  EXPECT_NEAR(grad(4, 2), 0.41211241149940E-03, 1E-6);
  EXPECT_NEAR(grad(5, 0), 0.34926883833109E-03, 1E-6);
  EXPECT_NEAR(grad(5, 1), -0.82490055826060E-03, 1E-6);
  EXPECT_NEAR(grad(5, 2), -0.17781855783076E-03, 1E-6);
  EXPECT_NEAR(grad(6, 0), 0.42168742012896E-03, 1E-6);
  EXPECT_NEAR(grad(6, 1), -0.15112033902461E-03, 1E-6);
  EXPECT_NEAR(grad(6, 2), -0.13966357053311E-03, 1E-6);
  EXPECT_NEAR(grad(7, 0), 0.24966360661981E-03, 1E-6);
  EXPECT_NEAR(grad(7, 1), -0.15655268010447E-03, 1E-6);
  EXPECT_NEAR(grad(7, 2), -0.31693163328506E-03, 1E-6);
  EXPECT_NEAR(grad(8, 0), -0.33886952460325E-03, 1E-6);
  EXPECT_NEAR(grad(8, 1), -0.72089553381479E-03, 1E-6);
  EXPECT_NEAR(grad(8, 2), 0.31846066021328E-04, 1E-6);
  EXPECT_NEAR(grad(9, 0), -0.25240094251490E-03, 1E-6);
  EXPECT_NEAR(grad(9, 1), -0.34577715822096E-03, 1E-6);
  EXPECT_NEAR(grad(9, 2), -0.49134085063857E-04, 1E-6);
  EXPECT_NEAR(grad(10, 0), -0.35555910668347E-03, 1E-6);
  EXPECT_NEAR(grad(10, 1), -0.23249486263795E-03, 1E-6);
  EXPECT_NEAR(grad(10, 2), 0.13156985458581E-03, 1E-6);
  EXPECT_NEAR(grad(11, 0), 0.12618654550712E-04, 1E-6);
  EXPECT_NEAR(grad(11, 1), -0.73389359942959E-03, 1E-6);
  EXPECT_NEAR(grad(11, 2), -0.22662040619494E-03, 1E-6);
  EXPECT_NEAR(grad(12, 0), -0.20518813662789E-03, 1E-6);
  EXPECT_NEAR(grad(12, 1), -0.25305081529472E-03, 1E-6);
  EXPECT_NEAR(grad(12, 2), -0.29332412721929E-04, 1E-6);
  EXPECT_NEAR(grad(13, 0), -0.18149518653715E-03, 1E-6);
  EXPECT_NEAR(grad(13, 1), -0.29166450456192E-03, 1E-6);
  EXPECT_NEAR(grad(13, 2), 0.14474045083862E-03, 1E-6);
  EXPECT_NEAR(grad(14, 0), 0.59851912379569E-04, 1E-6);
  EXPECT_NEAR(grad(14, 1), -0.30396868171866E-03, 1E-6);
  EXPECT_NEAR(grad(14, 2), -0.23338310260247E-03, 1E-6);
  EXPECT_NEAR(grad(15, 0), 0.58197898383952E-04, 1E-6);
  EXPECT_NEAR(grad(15, 1), -0.30758361044059E-03, 1E-6);
  EXPECT_NEAR(grad(15, 2), -0.64095555823846E-04, 1E-6);
  EXPECT_NEAR(grad(16, 0), 0.17270702487542E-03, 1E-6);
  EXPECT_NEAR(grad(16, 1), -0.40325736635360E-03, 1E-6);
  EXPECT_NEAR(grad(16, 2), 0.11713155473706E-03, 1E-6);
  EXPECT_NEAR(grad(17, 0), 0.22413872347098E-03, 1E-6);
  EXPECT_NEAR(grad(17, 1), -0.29544585454317E-03, 1E-6);
  EXPECT_NEAR(grad(17, 2), -0.15742229180132E-03, 1E-6);
  EXPECT_NEAR(grad(18, 0), 0.74258414458455E-03, 1E-6);
  EXPECT_NEAR(grad(18, 1), -0.23544904353269E-03, 1E-6);
  EXPECT_NEAR(grad(18, 2), -0.79629863113431E-04, 1E-6);
  EXPECT_NEAR(grad(19, 0), -0.62718318770087E-03, 1E-6);
  EXPECT_NEAR(grad(19, 1), 0.81265063752082E-05, 1E-6);
  EXPECT_NEAR(grad(19, 2), 0.82533454907008E-03, 1E-6);
  EXPECT_NEAR(grad(20, 0), 0.33646045007112E-03, 1E-6);
  EXPECT_NEAR(grad(20, 1), -0.60645359829581E-03, 1E-6);
  EXPECT_NEAR(grad(20, 2), 0.26393961962227E-03, 1E-6);
  EXPECT_NEAR(grad(21, 0), -0.38216105795221E-03, 1E-6);
  EXPECT_NEAR(grad(21, 1), -0.76192427763445E-03, 1E-6);
  EXPECT_NEAR(grad(21, 2), 0.25568559215776E-02, 1E-6);
  EXPECT_NEAR(grad(22, 0), 0.26597477353076E-03, 1E-6);
  EXPECT_NEAR(grad(22, 1), -0.37806614263501E-04, 1E-6);
  EXPECT_NEAR(grad(22, 2), -0.14069344470352E-03, 1E-6);
  EXPECT_NEAR(grad(23, 0), -0.12559297973960E-03, 1E-6);
  EXPECT_NEAR(grad(23, 1), 0.27274113327985E-04, 1E-6);
  EXPECT_NEAR(grad(23, 2), 0.15009052011999E-03, 1E-6);
  EXPECT_NEAR(grad(24, 0), -0.89435530750427E-03, 1E-6);
  EXPECT_NEAR(grad(24, 1), -0.16592489806786E-03, 1E-6);
  EXPECT_NEAR(grad(24, 2), 0.33668579975179E-03, 1E-6);
  EXPECT_NEAR(grad(25, 0), -0.56030163703870E-03, 1E-6);
  EXPECT_NEAR(grad(25, 1), 0.20168114985869E-03, 1E-6);
  EXPECT_NEAR(grad(25, 2), 0.97808373610792E-03, 1E-6);
  EXPECT_NEAR(grad(26, 0), -0.17707230556490E-03, 1E-6);
  EXPECT_NEAR(grad(26, 1), 0.18050557194038E-03, 1E-6);
  EXPECT_NEAR(grad(26, 2), 0.55571414787153E-03, 1E-6);
  EXPECT_NEAR(grad(27, 0), -0.29024705957479E-03, 1E-6);
  EXPECT_NEAR(grad(27, 1), -0.17002227258954E-03, 1E-6);
  EXPECT_NEAR(grad(27, 2), 0.30776643002116E-03, 1E-6);
  EXPECT_NEAR(grad(28, 0), 0.90380767839685E-04, 1E-6);
  EXPECT_NEAR(grad(28, 1), -0.17547077507061E-03, 1E-6);
  EXPECT_NEAR(grad(28, 2), 0.11527519013837E-02, 1E-6);
  EXPECT_NEAR(grad(29, 0), 0.60467872781642E-03, 1E-6);
  EXPECT_NEAR(grad(29, 1), -0.69815917476511E-04, 1E-6);
  EXPECT_NEAR(grad(29, 2), -0.28920962679217E-04, 1E-6);
  EXPECT_NEAR(grad(30, 0), 0.21692670351778E-03, 1E-6);
  EXPECT_NEAR(grad(30, 1), -0.23683311360878E-03, 1E-6);
  EXPECT_NEAR(grad(30, 2), -0.61253246910453E-04, 1E-6);
  EXPECT_NEAR(grad(31, 0), 0.68598143263506E-03, 1E-6);
  EXPECT_NEAR(grad(31, 1), -0.59327424588072E-03, 1E-6);
  EXPECT_NEAR(grad(31, 2), -0.78393376089589E-04, 1E-6);
  EXPECT_NEAR(grad(32, 0), 0.22121444984062E-03, 1E-6);
  EXPECT_NEAR(grad(32, 1), 0.10684997460277E-03, 1E-6);
  EXPECT_NEAR(grad(32, 2), -0.31297895549070E-04, 1E-6);
  EXPECT_NEAR(grad(33, 0), -0.27577932116675E-03, 1E-6);
  EXPECT_NEAR(grad(33, 1), 0.21286548768611E-04, 1E-6);
  EXPECT_NEAR(grad(33, 2), -0.17794077247484E-03, 1E-6);
  EXPECT_NEAR(grad(34, 0), -0.14570515846532E-03, 1E-6);
  EXPECT_NEAR(grad(34, 1), -0.37470590632276E-03, 1E-6);
  EXPECT_NEAR(grad(34, 2), 0.32853483383947E-03, 1E-6);
  EXPECT_NEAR(grad(35, 0), -0.37789482925038E-03, 1E-6);
  EXPECT_NEAR(grad(35, 1), 0.25579938175842E-04, 1E-6);
  EXPECT_NEAR(grad(35, 2), 0.62430511058046E-03, 1E-6);
  EXPECT_NEAR(grad(36, 0), -0.25240539340760E-03, 1E-6);
  EXPECT_NEAR(grad(36, 1), -0.22229593227934E-03, 1E-6);
  EXPECT_NEAR(grad(36, 2), 0.33751092471343E-03, 1E-6);
  EXPECT_NEAR(grad(37, 0), -0.99806148320629E-04, 1E-6);
  EXPECT_NEAR(grad(37, 1), 0.52318407675342E-03, 1E-6);
  EXPECT_NEAR(grad(37, 2), 0.33750626299077E-03, 1E-6);
  EXPECT_NEAR(grad(38, 0), 0.71646170857567E-04, 1E-6);
  EXPECT_NEAR(grad(38, 1), -0.14622670192495E-03, 1E-6);
  EXPECT_NEAR(grad(38, 2), -0.20930239095119E-04, 1E-6);
  EXPECT_NEAR(grad(39, 0), -0.78957269756303E-05, 1E-6);
  EXPECT_NEAR(grad(39, 1), -0.25462654976172E-03, 1E-6);
  EXPECT_NEAR(grad(39, 2), -0.57874566181241E-04, 1E-6);
  EXPECT_NEAR(grad(40, 0), -0.15943888916327E-03, 1E-6);
  EXPECT_NEAR(grad(40, 1), -0.35355599372866E-03, 1E-6);
  EXPECT_NEAR(grad(40, 2), -0.14771510349558E-03, 1E-6);
  EXPECT_NEAR(grad(41, 0), 0.34617678496598E-03, 1E-6);
  EXPECT_NEAR(grad(41, 1), 0.39580365613498E-03, 1E-6);
  EXPECT_NEAR(grad(41, 2), 0.91422688016485E-04, 1E-6);
  EXPECT_NEAR(grad(42, 0), -0.55991538722663E-03, 1E-6);
  EXPECT_NEAR(grad(42, 1), 0.18197621244918E-03, 1E-6);
  EXPECT_NEAR(grad(42, 2), 0.19826891379546E-03, 1E-6);
  EXPECT_NEAR(grad(43, 0), -0.82657213360346E-04, 1E-6);
  EXPECT_NEAR(grad(43, 1), 0.24581450328518E-03, 1E-6);
  EXPECT_NEAR(grad(43, 2), 0.10548903322119E-03, 1E-6);
  EXPECT_NEAR(grad(44, 0), 0.36352804432610E-03, 1E-6);
  EXPECT_NEAR(grad(44, 1), -0.17103399181954E-03, 1E-6);
  EXPECT_NEAR(grad(44, 2), 0.94950010168654E-04, 1E-6);
  EXPECT_NEAR(grad(45, 0), 0.16321065509730E-03, 1E-6);
  EXPECT_NEAR(grad(45, 1), 0.15806354080884E-03, 1E-6);
  EXPECT_NEAR(grad(45, 2), 0.76646347836603E-04, 1E-6);
  EXPECT_NEAR(grad(46, 0), -0.22466349605483E-03, 1E-6);
  EXPECT_NEAR(grad(46, 1), 0.38022088714113E-03, 1E-6);
  EXPECT_NEAR(grad(46, 2), 0.76795574424734E-03, 1E-6);
  EXPECT_NEAR(grad(47, 0), -0.53079207186461E-04, 1E-6);
  EXPECT_NEAR(grad(47, 1), 0.62064532904266E-03, 1E-6);
  EXPECT_NEAR(grad(47, 2), 0.27135325045856E-03, 1E-6);
  EXPECT_NEAR(grad(48, 0), 0.42294819816451E-05, 1E-6);
  EXPECT_NEAR(grad(48, 1), 0.93417455133485E-03, 1E-6);
  EXPECT_NEAR(grad(48, 2), -0.48194036745058E-04, 1E-6);
  EXPECT_NEAR(grad(49, 0), -0.73072997626114E-03, 1E-6);
  EXPECT_NEAR(grad(49, 1), -0.16790914373141E-03, 1E-6);
  EXPECT_NEAR(grad(49, 2), 0.12969542302391E-03, 1E-6);
  EXPECT_NEAR(grad(50, 0), -0.53962846497051E-03, 1E-6);
  EXPECT_NEAR(grad(50, 1), 0.61850004161158E-03, 1E-6);
  EXPECT_NEAR(grad(50, 2), -0.16487257513383E-03, 1E-6);
  EXPECT_NEAR(grad(51, 0), -0.51370679823668E-03, 1E-6);
  EXPECT_NEAR(grad(51, 1), 0.17754619315770E-03, 1E-6);
  EXPECT_NEAR(grad(51, 2), 0.39075759885388E-03, 1E-6);
  EXPECT_NEAR(grad(52, 0), 0.45908338324595E-03, 1E-6);
  EXPECT_NEAR(grad(52, 1), 0.55353094921757E-03, 1E-6);
  EXPECT_NEAR(grad(52, 2), -0.44842515537039E-03, 1E-6);
  EXPECT_NEAR(grad(53, 0), 0.34841686352288E-03, 1E-6);
  EXPECT_NEAR(grad(53, 1), 0.52118523578022E-03, 1E-6);
  EXPECT_NEAR(grad(53, 2), 0.17892604373401E-03, 1E-6);
  EXPECT_NEAR(grad(54, 0), 0.46697049469676E-03, 1E-6);
  EXPECT_NEAR(grad(54, 1), 0.41484027389534E-03, 1E-6);
  EXPECT_NEAR(grad(54, 2), 0.67521222571826E-03, 1E-6);
  EXPECT_NEAR(grad(55, 0), 0.25067499114087E-03, 1E-6);
  EXPECT_NEAR(grad(55, 1), -0.51758481789554E-03, 1E-6);
  EXPECT_NEAR(grad(55, 2), 0.83715997387041E-07, 1E-6);
  EXPECT_NEAR(grad(56, 0), 0.31860379935697E-03, 1E-6);
  EXPECT_NEAR(grad(56, 1), -0.71775181245478E-04, 1E-6);
  EXPECT_NEAR(grad(56, 2), 0.54818893716234E-03, 1E-6);
  EXPECT_NEAR(grad(57, 0), 0.58904908413006E-03, 1E-6);
  EXPECT_NEAR(grad(57, 1), 0.59035359878928E-04, 1E-6);
  EXPECT_NEAR(grad(57, 2), -0.24449075170745E-03, 1E-6);
  EXPECT_NEAR(grad(58, 0), 0.32271499852984E-03, 1E-6);
  EXPECT_NEAR(grad(58, 1), 0.42598484781490E-04, 1E-6);
  EXPECT_NEAR(grad(58, 2), -0.30169749453048E-03, 1E-6);
  EXPECT_NEAR(grad(59, 0), 0.32660165054343E-03, 1E-6);
  EXPECT_NEAR(grad(59, 1), 0.38605731082703E-04, 1E-6);
  EXPECT_NEAR(grad(59, 2), -0.10487205252228E-03, 1E-6);
  EXPECT_NEAR(grad(60, 0), 0.36088023281188E-03, 1E-6);
  EXPECT_NEAR(grad(60, 1), 0.13552760824491E-03, 1E-6);
  EXPECT_NEAR(grad(60, 2), -0.14472718762692E-03, 1E-6);
  EXPECT_NEAR(grad(61, 0), 0.74532595777986E-04, 1E-6);
  EXPECT_NEAR(grad(61, 1), -0.78743677641970E-04, 1E-6);
  EXPECT_NEAR(grad(61, 2), 0.43308838087127E-03, 1E-6);
  EXPECT_NEAR(grad(62, 0), 0.18902720151239E-03, 1E-6);
  EXPECT_NEAR(grad(62, 1), 0.53123293281693E-04, 1E-6);
  EXPECT_NEAR(grad(62, 2), 0.36198903651143E-03, 1E-6);
  EXPECT_NEAR(grad(63, 0), 0.19237491558608E-03, 1E-6);
  EXPECT_NEAR(grad(63, 1), -0.25827682670203E-04, 1E-6);
  EXPECT_NEAR(grad(63, 2), 0.28652123991647E-03, 1E-6);
  EXPECT_NEAR(grad(64, 0), 0.61361902239109E-04, 1E-6);
  EXPECT_NEAR(grad(64, 1), -0.36233159122438E-03, 1E-6);
  EXPECT_NEAR(grad(64, 2), 0.89927241296099E-04, 1E-6);
  EXPECT_NEAR(grad(65, 0), 0.18917028088016E-03, 1E-6);
  EXPECT_NEAR(grad(65, 1), -0.28227079029211E-03, 1E-6);
  EXPECT_NEAR(grad(65, 2), 0.20694866821663E-04, 1E-6);
  EXPECT_NEAR(grad(66, 0), 0.13685459876209E-03, 1E-6);
  EXPECT_NEAR(grad(66, 1), -0.32427335232816E-03, 1E-6);
  EXPECT_NEAR(grad(66, 2), -0.13300546926098E-03, 1E-6);
  EXPECT_NEAR(grad(67, 0), 0.26114841015476E-03, 1E-6);
  EXPECT_NEAR(grad(67, 1), 0.29580576192710E-03, 1E-6);
  EXPECT_NEAR(grad(67, 2), -0.90565876818288E-05, 1E-6);
  EXPECT_NEAR(grad(68, 0), 0.18866388073694E-03, 1E-6);
  EXPECT_NEAR(grad(68, 1), 0.25954541664870E-03, 1E-6);
  EXPECT_NEAR(grad(68, 2), 0.23396174434179E-03, 1E-6);
  EXPECT_NEAR(grad(69, 0), 0.17604268542655E-03, 1E-6);
  EXPECT_NEAR(grad(69, 1), 0.29497001237564E-03, 1E-6);
  EXPECT_NEAR(grad(69, 2), 0.96140210193740E-04, 1E-6);
  EXPECT_NEAR(grad(70, 0), 0.14545632653533E-03, 1E-6);
  EXPECT_NEAR(grad(70, 1), 0.16282027701304E-03, 1E-6);
  EXPECT_NEAR(grad(70, 2), 0.40023213443572E-03, 1E-6);
  EXPECT_NEAR(grad(71, 0), 0.12734994173036E-03, 1E-6);
  EXPECT_NEAR(grad(71, 1), 0.24171147927789E-03, 1E-6);
  EXPECT_NEAR(grad(71, 2), 0.31135415993650E-03, 1E-6);
  EXPECT_NEAR(grad(72, 0), -0.29837504830910E-04, 1E-6);
  EXPECT_NEAR(grad(72, 1), 0.14589381150503E-03, 1E-6);
  EXPECT_NEAR(grad(72, 2), 0.46248860864016E-03, 1E-6);
  EXPECT_NEAR(grad(73, 0), 0.11164281602966E-03, 1E-6);
  EXPECT_NEAR(grad(73, 1), 0.26733244028954E-03, 1E-6);
  EXPECT_NEAR(grad(73, 2), -0.38442977610294E-03, 1E-6);
  EXPECT_NEAR(grad(74, 0), 0.14780868927160E-03, 1E-6);
  EXPECT_NEAR(grad(74, 1), 0.30316845096713E-03, 1E-6);
  EXPECT_NEAR(grad(74, 2), -0.22411480433652E-03, 1E-6);
  EXPECT_NEAR(grad(75, 0), 0.21598953449214E-03, 1E-6);
  EXPECT_NEAR(grad(75, 1), 0.27172976063054E-03, 1E-6);
  EXPECT_NEAR(grad(75, 2), -0.36889489987574E-03, 1E-6);
  EXPECT_NEAR(grad(76, 0), -0.23011551190199E-04, 1E-6);
  EXPECT_NEAR(grad(76, 1), 0.17497475449326E-04, 1E-6);
  EXPECT_NEAR(grad(76, 2), 0.49288230676588E-03, 1E-6);
  EXPECT_NEAR(grad(77, 0), -0.10203454496426E-04, 1E-6);
  EXPECT_NEAR(grad(77, 1), 0.14242905768617E-03, 1E-6);
  EXPECT_NEAR(grad(77, 2), 0.38413180206949E-03, 1E-6);
  EXPECT_NEAR(grad(78, 0), -0.92951046630347E-04, 1E-6);
  EXPECT_NEAR(grad(78, 1), 0.11334897774169E-03, 1E-6);
  EXPECT_NEAR(grad(78, 2), 0.50394935178270E-03, 1E-6);
  EXPECT_NEAR(grad(79, 0), -0.88145265817250E-04, 1E-6);
  EXPECT_NEAR(grad(79, 1), 0.30910113573110E-03, 1E-6);
  EXPECT_NEAR(grad(79, 2), 0.24064919432476E-03, 1E-6);
  EXPECT_NEAR(grad(80, 0), -0.34686924686719E-04, 1E-6);
  EXPECT_NEAR(grad(80, 1), 0.41713008140609E-03, 1E-6);
  EXPECT_NEAR(grad(80, 2), 0.35168136596270E-04, 1E-6);
  EXPECT_NEAR(grad(81, 0), -0.71190025840919E-05, 1E-6);
  EXPECT_NEAR(grad(81, 1), 0.33388172698515E-03, 1E-6);
  EXPECT_NEAR(grad(81, 2), 0.14037391038134E-03, 1E-6);
  EXPECT_NEAR(grad(82, 0), 0.80005192253559E-04, 1E-6);
  EXPECT_NEAR(grad(82, 1), 0.40712967033319E-03, 1E-6);
  EXPECT_NEAR(grad(82, 2), -0.15037380068975E-03, 1E-6);
  EXPECT_NEAR(grad(83, 0), 0.61578385780451E-04, 1E-6);
  EXPECT_NEAR(grad(83, 1), 0.41342093506423E-03, 1E-6);
  EXPECT_NEAR(grad(83, 2), -0.68778364115291E-04, 1E-6);
  EXPECT_NEAR(grad(84, 0), 0.22008706351074E-03, 1E-6);
  EXPECT_NEAR(grad(84, 1), 0.42823921889673E-03, 1E-6);
  EXPECT_NEAR(grad(84, 2), -0.17670070301097E-03, 1E-6);
  EXPECT_NEAR(grad(85, 0), -0.47093767210147E-03, 1E-6);
  EXPECT_NEAR(grad(85, 1), -0.11045224562088E-03, 1E-6);
  EXPECT_NEAR(grad(85, 2), -0.63939801876983E-05, 1E-6);
  EXPECT_NEAR(grad(86, 0), -0.34733167111957E-03, 1E-6);
  EXPECT_NEAR(grad(86, 1), -0.89343547272491E-04, 1E-6);
  EXPECT_NEAR(grad(86, 2), -0.72215719508202E-04, 1E-6);
  EXPECT_NEAR(grad(87, 0), -0.32914022248560E-03, 1E-6);
  EXPECT_NEAR(grad(87, 1), -0.23658566934120E-03, 1E-6);
  EXPECT_NEAR(grad(87, 2), 0.49447256654468E-04, 1E-6);
  EXPECT_NEAR(grad(88, 0), -0.28410391854520E-03, 1E-6);
  EXPECT_NEAR(grad(88, 1), 0.79951680672349E-05, 1E-6);
  EXPECT_NEAR(grad(88, 2), 0.35456205757965E-03, 1E-6);
  EXPECT_NEAR(grad(89, 0), -0.27621949446983E-03, 1E-6);
  EXPECT_NEAR(grad(89, 1), 0.11502136251572E-03, 1E-6);
  EXPECT_NEAR(grad(89, 2), 0.18501567946618E-03, 1E-6);
  EXPECT_NEAR(grad(90, 0), -0.26374185825814E-03, 1E-6);
  EXPECT_NEAR(grad(90, 1), 0.18294741624755E-03, 1E-6);
  EXPECT_NEAR(grad(90, 2), 0.26743998237160E-03, 1E-6);
  EXPECT_NEAR(grad(91, 0), -0.28464477755201E-03, 1E-6);
  EXPECT_NEAR(grad(91, 1), 0.41422672051734E-03, 1E-6);
  EXPECT_NEAR(grad(91, 2), -0.25288171029992E-03, 1E-6);
  EXPECT_NEAR(grad(92, 0), -0.16659418114369E-03, 1E-6);
  EXPECT_NEAR(grad(92, 1), 0.38592109853208E-03, 1E-6);
  EXPECT_NEAR(grad(92, 2), -0.16891257434185E-03, 1E-6);
  EXPECT_NEAR(grad(93, 0), -0.20144134086315E-03, 1E-6);
  EXPECT_NEAR(grad(93, 1), 0.25390934516422E-03, 1E-6);
  EXPECT_NEAR(grad(93, 2), -0.17607477637070E-03, 1E-6);
  EXPECT_NEAR(grad(94, 0), 0.61956940549665E-03, 1E-6);
  EXPECT_NEAR(grad(94, 1), 0.12929713142598E-03, 1E-6);
  EXPECT_NEAR(grad(94, 2), -0.75110524825660E-03, 1E-6);
  EXPECT_NEAR(grad(95, 0), 0.35842823225359E-03, 1E-6);
  EXPECT_NEAR(grad(95, 1), 0.61152755961755E-03, 1E-6);
  EXPECT_NEAR(grad(95, 2), -0.14386526517954E-02, 1E-6);
  EXPECT_NEAR(grad(96, 0), 0.52686288364152E-03, 1E-6);
  EXPECT_NEAR(grad(96, 1), 0.83560382839702E-03, 1E-6);
  EXPECT_NEAR(grad(96, 2), -0.14849934797071E-02, 1E-6);
  EXPECT_NEAR(grad(97, 0), 0.45824215342805E-04, 1E-6);
  EXPECT_NEAR(grad(97, 1), -0.19573787818676E-04, 1E-6);
  EXPECT_NEAR(grad(97, 2), -0.13596429497376E-02, 1E-6);
  EXPECT_NEAR(grad(98, 0), 0.12387026472583E-03, 1E-6);
  EXPECT_NEAR(grad(98, 1), 0.28875242066941E-03, 1E-6);
  EXPECT_NEAR(grad(98, 2), -0.40897188139410E-03, 1E-6);
  EXPECT_NEAR(grad(99, 0), 0.26616086715217E-03, 1E-6);
  EXPECT_NEAR(grad(99, 1), 0.26632272603516E-03, 1E-6);
  EXPECT_NEAR(grad(99, 2), -0.52142481168618E-03, 1E-6);
  EXPECT_NEAR(grad(100, 0), 0.16809327333500E-03, 1E-6);
  EXPECT_NEAR(grad(100, 1), 0.29730991303690E-03, 1E-6);
  EXPECT_NEAR(grad(100, 2), -0.44517730312351E-03, 1E-6);
  EXPECT_NEAR(grad(101, 0), -0.24391893655121E-03, 1E-6);
  EXPECT_NEAR(grad(101, 1), -0.49919758320386E-04, 1E-6);
  EXPECT_NEAR(grad(101, 2), -0.14980324206555E-02, 1E-6);
  EXPECT_NEAR(grad(102, 0), -0.54985138737388E-04, 1E-6);
  EXPECT_NEAR(grad(102, 1), -0.29111611237035E-03, 1E-6);
  EXPECT_NEAR(grad(102, 2), -0.10893755016443E-02, 1E-6);
  EXPECT_NEAR(grad(103, 0), -0.29943334698998E-03, 1E-6);
  EXPECT_NEAR(grad(103, 1), -0.25751870955533E-03, 1E-6);
  EXPECT_NEAR(grad(103, 2), -0.11782161493584E-02, 1E-6);
  EXPECT_NEAR(grad(104, 0), 0.20893735291185E-03, 1E-6);
  EXPECT_NEAR(grad(104, 1), 0.68894745962174E-05, 1E-6);
  EXPECT_NEAR(grad(104, 2), -0.55234287566566E-03, 1E-6);
  EXPECT_NEAR(grad(105, 0), -0.24961979188683E-03, 1E-6);
  EXPECT_NEAR(grad(105, 1), -0.19789320532086E-03, 1E-6);
  EXPECT_NEAR(grad(105, 2), -0.61001233354147E-03, 1E-6);
  EXPECT_NEAR(grad(106, 0), 0.40346805616292E-05, 1E-6);
  EXPECT_NEAR(grad(106, 1), -0.15346977022997E-03, 1E-6);
  EXPECT_NEAR(grad(106, 2), -0.43010755490281E-03, 1E-6);
  EXPECT_NEAR(grad(107, 0), -0.23848344817731E-03, 1E-6);
  EXPECT_NEAR(grad(107, 1), -0.34174384810148E-03, 1E-6);
  EXPECT_NEAR(grad(107, 2), -0.65402035669312E-03, 1E-6);
  EXPECT_NEAR(grad(108, 0), -0.16202489500003E-03, 1E-6);
  EXPECT_NEAR(grad(108, 1), 0.61243359130732E-04, 1E-6);
  EXPECT_NEAR(grad(108, 2), -0.40939633341506E-03, 1E-6);
  EXPECT_NEAR(grad(109, 0), -0.13217743981632E-03, 1E-6);
  EXPECT_NEAR(grad(109, 1), -0.21326087667958E-03, 1E-6);
  EXPECT_NEAR(grad(109, 2), -0.13363237687104E-03, 1E-6);
  EXPECT_NEAR(grad(110, 0), -0.23758846996526E-03, 1E-6);
  EXPECT_NEAR(grad(110, 1), -0.92060837122905E-04, 1E-6);
  EXPECT_NEAR(grad(110, 2), -0.16064478190534E-03, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0)-ABC gradient correction with ghost atoms.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJABCGradBP86JacobsenGhost) {
  auto geometry = *systemControllerTwo->getGeometry();
  geometry += *systemControllerThree->getGeometry();
  auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJABC,
                                                                               std::make_shared<Geometry>(geometry),
                                                                               CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(grad(0, 0), -0.44419424846835E-03, 1E-6);
  EXPECT_NEAR(grad(0, 1), -0.34497327487658E-03, 1E-6);
  EXPECT_NEAR(grad(0, 2), 0.25751563332388E-03, 1E-6);
  EXPECT_NEAR(grad(1, 0), 0.44810919790884E-03, 1E-6);
  EXPECT_NEAR(grad(1, 1), -0.49797428237543E-03, 1E-6);
  EXPECT_NEAR(grad(1, 2), -0.24012114425014E-03, 1E-6);
  EXPECT_NEAR(grad(2, 0), -0.17523503799171E-03, 1E-6);
  EXPECT_NEAR(grad(2, 1), -0.17538269971019E-03, 1E-6);
  EXPECT_NEAR(grad(2, 2), 0.34753130839096E-03, 1E-6);
  EXPECT_NEAR(grad(3, 0), -0.62803624147614E-03, 1E-6);
  EXPECT_NEAR(grad(3, 1), -0.73165514136507E-03, 1E-6);
  EXPECT_NEAR(grad(3, 2), 0.15753938993305E-03, 1E-6);
  EXPECT_NEAR(grad(4, 0), -0.26356763387021E-03, 1E-6);
  EXPECT_NEAR(grad(4, 1), -0.78137545039569E-04, 1E-6);
  EXPECT_NEAR(grad(4, 2), 0.41211241149940E-03, 1E-6);
  EXPECT_NEAR(grad(5, 0), 0.34926883833109E-03, 1E-6);
  EXPECT_NEAR(grad(5, 1), -0.82490055826060E-03, 1E-6);
  EXPECT_NEAR(grad(5, 2), -0.17781855783076E-03, 1E-6);
  EXPECT_NEAR(grad(6, 0), 0.42168742012896E-03, 1E-6);
  EXPECT_NEAR(grad(6, 1), -0.15112033902461E-03, 1E-6);
  EXPECT_NEAR(grad(6, 2), -0.13966357053311E-03, 1E-6);
  EXPECT_NEAR(grad(7, 0), 0.24966360661981E-03, 1E-6);
  EXPECT_NEAR(grad(7, 1), -0.15655268010447E-03, 1E-6);
  EXPECT_NEAR(grad(7, 2), -0.31693163328506E-03, 1E-6);
  EXPECT_NEAR(grad(8, 0), -0.33886952460325E-03, 1E-6);
  EXPECT_NEAR(grad(8, 1), -0.72089553381479E-03, 1E-6);
  EXPECT_NEAR(grad(8, 2), 0.31846066021328E-04, 1E-6);
  EXPECT_NEAR(grad(9, 0), -0.25240094251490E-03, 1E-6);
  EXPECT_NEAR(grad(9, 1), -0.34577715822096E-03, 1E-6);
  EXPECT_NEAR(grad(9, 2), -0.49134085063857E-04, 1E-6);
  EXPECT_NEAR(grad(10, 0), -0.35555910668347E-03, 1E-6);
  EXPECT_NEAR(grad(10, 1), -0.23249486263795E-03, 1E-6);
  EXPECT_NEAR(grad(10, 2), 0.13156985458581E-03, 1E-6);
  EXPECT_NEAR(grad(11, 0), 0.12618654550712E-04, 1E-6);
  EXPECT_NEAR(grad(11, 1), -0.73389359942959E-03, 1E-6);
  EXPECT_NEAR(grad(11, 2), -0.22662040619494E-03, 1E-6);
  EXPECT_NEAR(grad(12, 0), -0.20518813662789E-03, 1E-6);
  EXPECT_NEAR(grad(12, 1), -0.25305081529472E-03, 1E-6);
  EXPECT_NEAR(grad(12, 2), -0.29332412721929E-04, 1E-6);
  EXPECT_NEAR(grad(13, 0), -0.18149518653715E-03, 1E-6);
  EXPECT_NEAR(grad(13, 1), -0.29166450456192E-03, 1E-6);
  EXPECT_NEAR(grad(13, 2), 0.14474045083862E-03, 1E-6);
  EXPECT_NEAR(grad(14, 0), 0.59851912379569E-04, 1E-6);
  EXPECT_NEAR(grad(14, 1), -0.30396868171866E-03, 1E-6);
  EXPECT_NEAR(grad(14, 2), -0.23338310260247E-03, 1E-6);
  EXPECT_NEAR(grad(15, 0), 0.58197898383952E-04, 1E-6);
  EXPECT_NEAR(grad(15, 1), -0.30758361044059E-03, 1E-6);
  EXPECT_NEAR(grad(15, 2), -0.64095555823846E-04, 1E-6);
  EXPECT_NEAR(grad(16, 0), 0.17270702487542E-03, 1E-6);
  EXPECT_NEAR(grad(16, 1), -0.40325736635360E-03, 1E-6);
  EXPECT_NEAR(grad(16, 2), 0.11713155473706E-03, 1E-6);
  EXPECT_NEAR(grad(17, 0), 0.22413872347098E-03, 1E-6);
  EXPECT_NEAR(grad(17, 1), -0.29544585454317E-03, 1E-6);
  EXPECT_NEAR(grad(17, 2), -0.15742229180132E-03, 1E-6);
  EXPECT_NEAR(grad(18, 0), 0.74258414458455E-03, 1E-6);
  EXPECT_NEAR(grad(18, 1), -0.23544904353269E-03, 1E-6);
  EXPECT_NEAR(grad(18, 2), -0.79629863113431E-04, 1E-6);
  EXPECT_NEAR(grad(19, 0), -0.62718318770087E-03, 1E-6);
  EXPECT_NEAR(grad(19, 1), 0.81265063752082E-05, 1E-6);
  EXPECT_NEAR(grad(19, 2), 0.82533454907008E-03, 1E-6);
  EXPECT_NEAR(grad(20, 0), 0.33646045007112E-03, 1E-6);
  EXPECT_NEAR(grad(20, 1), -0.60645359829581E-03, 1E-6);
  EXPECT_NEAR(grad(20, 2), 0.26393961962227E-03, 1E-6);
  EXPECT_NEAR(grad(21, 0), -0.38216105795221E-03, 1E-6);
  EXPECT_NEAR(grad(21, 1), -0.76192427763445E-03, 1E-6);
  EXPECT_NEAR(grad(21, 2), 0.25568559215776E-02, 1E-6);
  EXPECT_NEAR(grad(22, 0), 0.26597477353076E-03, 1E-6);
  EXPECT_NEAR(grad(22, 1), -0.37806614263501E-04, 1E-6);
  EXPECT_NEAR(grad(22, 2), -0.14069344470352E-03, 1E-6);
  EXPECT_NEAR(grad(23, 0), -0.12559297973960E-03, 1E-6);
  EXPECT_NEAR(grad(23, 1), 0.27274113327985E-04, 1E-6);
  EXPECT_NEAR(grad(23, 2), 0.15009052011999E-03, 1E-6);
  EXPECT_NEAR(grad(24, 0), -0.89435530750427E-03, 1E-6);
  EXPECT_NEAR(grad(24, 1), -0.16592489806786E-03, 1E-6);
  EXPECT_NEAR(grad(24, 2), 0.33668579975179E-03, 1E-6);
  EXPECT_NEAR(grad(25, 0), -0.56030163703870E-03, 1E-6);
  EXPECT_NEAR(grad(25, 1), 0.20168114985869E-03, 1E-6);
  EXPECT_NEAR(grad(25, 2), 0.97808373610792E-03, 1E-6);
  EXPECT_NEAR(grad(26, 0), -0.17707230556490E-03, 1E-6);
  EXPECT_NEAR(grad(26, 1), 0.18050557194038E-03, 1E-6);
  EXPECT_NEAR(grad(26, 2), 0.55571414787153E-03, 1E-6);
  EXPECT_NEAR(grad(27, 0), -0.29024705957479E-03, 1E-6);
  EXPECT_NEAR(grad(27, 1), -0.17002227258954E-03, 1E-6);
  EXPECT_NEAR(grad(27, 2), 0.30776643002116E-03, 1E-6);
  EXPECT_NEAR(grad(28, 0), 0.90380767839685E-04, 1E-6);
  EXPECT_NEAR(grad(28, 1), -0.17547077507061E-03, 1E-6);
  EXPECT_NEAR(grad(28, 2), 0.11527519013837E-02, 1E-6);
  EXPECT_NEAR(grad(29, 0), 0.60467872781642E-03, 1E-6);
  EXPECT_NEAR(grad(29, 1), -0.69815917476511E-04, 1E-6);
  EXPECT_NEAR(grad(29, 2), -0.28920962679217E-04, 1E-6);
  EXPECT_NEAR(grad(30, 0), 0.21692670351778E-03, 1E-6);
  EXPECT_NEAR(grad(30, 1), -0.23683311360878E-03, 1E-6);
  EXPECT_NEAR(grad(30, 2), -0.61253246910453E-04, 1E-6);
  EXPECT_NEAR(grad(31, 0), 0.68598143263506E-03, 1E-6);
  EXPECT_NEAR(grad(31, 1), -0.59327424588072E-03, 1E-6);
  EXPECT_NEAR(grad(31, 2), -0.78393376089589E-04, 1E-6);
  EXPECT_NEAR(grad(32, 0), 0.22121444984062E-03, 1E-6);
  EXPECT_NEAR(grad(32, 1), 0.10684997460277E-03, 1E-6);
  EXPECT_NEAR(grad(32, 2), -0.31297895549070E-04, 1E-6);
  EXPECT_NEAR(grad(33, 0), -0.27577932116675E-03, 1E-6);
  EXPECT_NEAR(grad(33, 1), 0.21286548768611E-04, 1E-6);
  EXPECT_NEAR(grad(33, 2), -0.17794077247484E-03, 1E-6);
  EXPECT_NEAR(grad(34, 0), -0.14570515846532E-03, 1E-6);
  EXPECT_NEAR(grad(34, 1), -0.37470590632276E-03, 1E-6);
  EXPECT_NEAR(grad(34, 2), 0.32853483383947E-03, 1E-6);
  EXPECT_NEAR(grad(35, 0), -0.37789482925038E-03, 1E-6);
  EXPECT_NEAR(grad(35, 1), 0.25579938175842E-04, 1E-6);
  EXPECT_NEAR(grad(35, 2), 0.62430511058046E-03, 1E-6);
  EXPECT_NEAR(grad(36, 0), -0.25240539340760E-03, 1E-6);
  EXPECT_NEAR(grad(36, 1), -0.22229593227934E-03, 1E-6);
  EXPECT_NEAR(grad(36, 2), 0.33751092471343E-03, 1E-6);
  EXPECT_NEAR(grad(37, 0), -0.99806148320629E-04, 1E-6);
  EXPECT_NEAR(grad(37, 1), 0.52318407675342E-03, 1E-6);
  EXPECT_NEAR(grad(37, 2), 0.33750626299077E-03, 1E-6);
  EXPECT_NEAR(grad(38, 0), 0.71646170857567E-04, 1E-6);
  EXPECT_NEAR(grad(38, 1), -0.14622670192495E-03, 1E-6);
  EXPECT_NEAR(grad(38, 2), -0.20930239095119E-04, 1E-6);
  EXPECT_NEAR(grad(39, 0), -0.78957269756303E-05, 1E-6);
  EXPECT_NEAR(grad(39, 1), -0.25462654976172E-03, 1E-6);
  EXPECT_NEAR(grad(39, 2), -0.57874566181241E-04, 1E-6);
  EXPECT_NEAR(grad(40, 0), -0.15943888916327E-03, 1E-6);
  EXPECT_NEAR(grad(40, 1), -0.35355599372866E-03, 1E-6);
  EXPECT_NEAR(grad(40, 2), -0.14771510349558E-03, 1E-6);
  EXPECT_NEAR(grad(41, 0), 0.34617678496598E-03, 1E-6);
  EXPECT_NEAR(grad(41, 1), 0.39580365613498E-03, 1E-6);
  EXPECT_NEAR(grad(41, 2), 0.91422688016485E-04, 1E-6);
  EXPECT_NEAR(grad(42, 0), -0.55991538722663E-03, 1E-6);
  EXPECT_NEAR(grad(42, 1), 0.18197621244918E-03, 1E-6);
  EXPECT_NEAR(grad(42, 2), 0.19826891379546E-03, 1E-6);
  EXPECT_NEAR(grad(43, 0), -0.82657213360346E-04, 1E-6);
  EXPECT_NEAR(grad(43, 1), 0.24581450328518E-03, 1E-6);
  EXPECT_NEAR(grad(43, 2), 0.10548903322119E-03, 1E-6);
  EXPECT_NEAR(grad(44, 0), 0.36352804432610E-03, 1E-6);
  EXPECT_NEAR(grad(44, 1), -0.17103399181954E-03, 1E-6);
  EXPECT_NEAR(grad(44, 2), 0.94950010168654E-04, 1E-6);
  EXPECT_NEAR(grad(45, 0), 0.16321065509730E-03, 1E-6);
  EXPECT_NEAR(grad(45, 1), 0.15806354080884E-03, 1E-6);
  EXPECT_NEAR(grad(45, 2), 0.76646347836603E-04, 1E-6);
  EXPECT_NEAR(grad(46, 0), -0.22466349605483E-03, 1E-6);
  EXPECT_NEAR(grad(46, 1), 0.38022088714113E-03, 1E-6);
  EXPECT_NEAR(grad(46, 2), 0.76795574424734E-03, 1E-6);
  EXPECT_NEAR(grad(47, 0), -0.53079207186461E-04, 1E-6);
  EXPECT_NEAR(grad(47, 1), 0.62064532904266E-03, 1E-6);
  EXPECT_NEAR(grad(47, 2), 0.27135325045856E-03, 1E-6);
  EXPECT_NEAR(grad(48, 0), 0.42294819816451E-05, 1E-6);
  EXPECT_NEAR(grad(48, 1), 0.93417455133485E-03, 1E-6);
  EXPECT_NEAR(grad(48, 2), -0.48194036745058E-04, 1E-6);
  EXPECT_NEAR(grad(49, 0), -0.73072997626114E-03, 1E-6);
  EXPECT_NEAR(grad(49, 1), -0.16790914373141E-03, 1E-6);
  EXPECT_NEAR(grad(49, 2), 0.12969542302391E-03, 1E-6);
  EXPECT_NEAR(grad(50, 0), -0.53962846497051E-03, 1E-6);
  EXPECT_NEAR(grad(50, 1), 0.61850004161158E-03, 1E-6);
  EXPECT_NEAR(grad(50, 2), -0.16487257513383E-03, 1E-6);
  EXPECT_NEAR(grad(51, 0), -0.51370679823668E-03, 1E-6);
  EXPECT_NEAR(grad(51, 1), 0.17754619315770E-03, 1E-6);
  EXPECT_NEAR(grad(51, 2), 0.39075759885388E-03, 1E-6);
  EXPECT_NEAR(grad(52, 0), 0.45908338324595E-03, 1E-6);
  EXPECT_NEAR(grad(52, 1), 0.55353094921757E-03, 1E-6);
  EXPECT_NEAR(grad(52, 2), -0.44842515537039E-03, 1E-6);
  EXPECT_NEAR(grad(53, 0), 0.34841686352288E-03, 1E-6);
  EXPECT_NEAR(grad(53, 1), 0.52118523578022E-03, 1E-6);
  EXPECT_NEAR(grad(53, 2), 0.17892604373401E-03, 1E-6);
  EXPECT_NEAR(grad(54, 0), 0.46697049469676E-03, 1E-6);
  EXPECT_NEAR(grad(54, 1), 0.41484027389534E-03, 1E-6);
  EXPECT_NEAR(grad(54, 2), 0.67521222571826E-03, 1E-6);
  EXPECT_NEAR(grad(55, 0), 0.25067499114087E-03, 1E-6);
  EXPECT_NEAR(grad(55, 1), -0.51758481789554E-03, 1E-6);
  EXPECT_NEAR(grad(55, 2), 0.83715997387041E-07, 1E-6);
  EXPECT_NEAR(grad(56, 0), 0.31860379935697E-03, 1E-6);
  EXPECT_NEAR(grad(56, 1), -0.71775181245478E-04, 1E-6);
  EXPECT_NEAR(grad(56, 2), 0.54818893716234E-03, 1E-6);
  EXPECT_NEAR(grad(57, 0), 0.58904908413006E-03, 1E-6);
  EXPECT_NEAR(grad(57, 1), 0.59035359878928E-04, 1E-6);
  EXPECT_NEAR(grad(57, 2), -0.24449075170745E-03, 1E-6);
  EXPECT_NEAR(grad(58, 0), 0.32271499852984E-03, 1E-6);
  EXPECT_NEAR(grad(58, 1), 0.42598484781490E-04, 1E-6);
  EXPECT_NEAR(grad(58, 2), -0.30169749453048E-03, 1E-6);
  EXPECT_NEAR(grad(59, 0), 0.32660165054343E-03, 1E-6);
  EXPECT_NEAR(grad(59, 1), 0.38605731082703E-04, 1E-6);
  EXPECT_NEAR(grad(59, 2), -0.10487205252228E-03, 1E-6);
  EXPECT_NEAR(grad(60, 0), 0.36088023281188E-03, 1E-6);
  EXPECT_NEAR(grad(60, 1), 0.13552760824491E-03, 1E-6);
  EXPECT_NEAR(grad(60, 2), -0.14472718762692E-03, 1E-6);
  EXPECT_NEAR(grad(61, 0), 0.74532595777986E-04, 1E-6);
  EXPECT_NEAR(grad(61, 1), -0.78743677641970E-04, 1E-6);
  EXPECT_NEAR(grad(61, 2), 0.43308838087127E-03, 1E-6);
  EXPECT_NEAR(grad(62, 0), 0.18902720151239E-03, 1E-6);
  EXPECT_NEAR(grad(62, 1), 0.53123293281693E-04, 1E-6);
  EXPECT_NEAR(grad(62, 2), 0.36198903651143E-03, 1E-6);
  EXPECT_NEAR(grad(63, 0), 0.19237491558608E-03, 1E-6);
  EXPECT_NEAR(grad(63, 1), -0.25827682670203E-04, 1E-6);
  EXPECT_NEAR(grad(63, 2), 0.28652123991647E-03, 1E-6);
  EXPECT_NEAR(grad(64, 0), 0.61361902239109E-04, 1E-6);
  EXPECT_NEAR(grad(64, 1), -0.36233159122438E-03, 1E-6);
  EXPECT_NEAR(grad(64, 2), 0.89927241296099E-04, 1E-6);
  EXPECT_NEAR(grad(65, 0), 0.18917028088016E-03, 1E-6);
  EXPECT_NEAR(grad(65, 1), -0.28227079029211E-03, 1E-6);
  EXPECT_NEAR(grad(65, 2), 0.20694866821663E-04, 1E-6);
  EXPECT_NEAR(grad(66, 0), 0.13685459876209E-03, 1E-6);
  EXPECT_NEAR(grad(66, 1), -0.32427335232816E-03, 1E-6);
  EXPECT_NEAR(grad(66, 2), -0.13300546926098E-03, 1E-6);
  EXPECT_NEAR(grad(67, 0), 0.26114841015476E-03, 1E-6);
  EXPECT_NEAR(grad(67, 1), 0.29580576192710E-03, 1E-6);
  EXPECT_NEAR(grad(67, 2), -0.90565876818288E-05, 1E-6);
  EXPECT_NEAR(grad(68, 0), 0.18866388073694E-03, 1E-6);
  EXPECT_NEAR(grad(68, 1), 0.25954541664870E-03, 1E-6);
  EXPECT_NEAR(grad(68, 2), 0.23396174434179E-03, 1E-6);
  EXPECT_NEAR(grad(69, 0), 0.17604268542655E-03, 1E-6);
  EXPECT_NEAR(grad(69, 1), 0.29497001237564E-03, 1E-6);
  EXPECT_NEAR(grad(69, 2), 0.96140210193740E-04, 1E-6);
  EXPECT_NEAR(grad(70, 0), 0.14545632653533E-03, 1E-6);
  EXPECT_NEAR(grad(70, 1), 0.16282027701304E-03, 1E-6);
  EXPECT_NEAR(grad(70, 2), 0.40023213443572E-03, 1E-6);
  EXPECT_NEAR(grad(71, 0), 0.12734994173036E-03, 1E-6);
  EXPECT_NEAR(grad(71, 1), 0.24171147927789E-03, 1E-6);
  EXPECT_NEAR(grad(71, 2), 0.31135415993650E-03, 1E-6);
  EXPECT_NEAR(grad(72, 0), -0.29837504830910E-04, 1E-6);
  EXPECT_NEAR(grad(72, 1), 0.14589381150503E-03, 1E-6);
  EXPECT_NEAR(grad(72, 2), 0.46248860864016E-03, 1E-6);
  EXPECT_NEAR(grad(73, 0), 0.11164281602966E-03, 1E-6);
  EXPECT_NEAR(grad(73, 1), 0.26733244028954E-03, 1E-6);
  EXPECT_NEAR(grad(73, 2), -0.38442977610294E-03, 1E-6);
  EXPECT_NEAR(grad(74, 0), 0.14780868927160E-03, 1E-6);
  EXPECT_NEAR(grad(74, 1), 0.30316845096713E-03, 1E-6);
  EXPECT_NEAR(grad(74, 2), -0.22411480433652E-03, 1E-6);
  EXPECT_NEAR(grad(75, 0), 0.21598953449214E-03, 1E-6);
  EXPECT_NEAR(grad(75, 1), 0.27172976063054E-03, 1E-6);
  EXPECT_NEAR(grad(75, 2), -0.36889489987574E-03, 1E-6);
  EXPECT_NEAR(grad(76, 0), -0.23011551190199E-04, 1E-6);
  EXPECT_NEAR(grad(76, 1), 0.17497475449326E-04, 1E-6);
  EXPECT_NEAR(grad(76, 2), 0.49288230676588E-03, 1E-6);
  EXPECT_NEAR(grad(77, 0), -0.10203454496426E-04, 1E-6);
  EXPECT_NEAR(grad(77, 1), 0.14242905768617E-03, 1E-6);
  EXPECT_NEAR(grad(77, 2), 0.38413180206949E-03, 1E-6);
  EXPECT_NEAR(grad(78, 0), -0.92951046630347E-04, 1E-6);
  EXPECT_NEAR(grad(78, 1), 0.11334897774169E-03, 1E-6);
  EXPECT_NEAR(grad(78, 2), 0.50394935178270E-03, 1E-6);
  EXPECT_NEAR(grad(79, 0), -0.88145265817250E-04, 1E-6);
  EXPECT_NEAR(grad(79, 1), 0.30910113573110E-03, 1E-6);
  EXPECT_NEAR(grad(79, 2), 0.24064919432476E-03, 1E-6);
  EXPECT_NEAR(grad(80, 0), -0.34686924686719E-04, 1E-6);
  EXPECT_NEAR(grad(80, 1), 0.41713008140609E-03, 1E-6);
  EXPECT_NEAR(grad(80, 2), 0.35168136596270E-04, 1E-6);
  EXPECT_NEAR(grad(81, 0), -0.71190025840919E-05, 1E-6);
  EXPECT_NEAR(grad(81, 1), 0.33388172698515E-03, 1E-6);
  EXPECT_NEAR(grad(81, 2), 0.14037391038134E-03, 1E-6);
  EXPECT_NEAR(grad(82, 0), 0.80005192253559E-04, 1E-6);
  EXPECT_NEAR(grad(82, 1), 0.40712967033319E-03, 1E-6);
  EXPECT_NEAR(grad(82, 2), -0.15037380068975E-03, 1E-6);
  EXPECT_NEAR(grad(83, 0), 0.61578385780451E-04, 1E-6);
  EXPECT_NEAR(grad(83, 1), 0.41342093506423E-03, 1E-6);
  EXPECT_NEAR(grad(83, 2), -0.68778364115291E-04, 1E-6);
  EXPECT_NEAR(grad(84, 0), 0.22008706351074E-03, 1E-6);
  EXPECT_NEAR(grad(84, 1), 0.42823921889673E-03, 1E-6);
  EXPECT_NEAR(grad(84, 2), -0.17670070301097E-03, 1E-6);
  EXPECT_NEAR(grad(85, 0), -0.47093767210147E-03, 1E-6);
  EXPECT_NEAR(grad(85, 1), -0.11045224562088E-03, 1E-6);
  EXPECT_NEAR(grad(85, 2), -0.63939801876983E-05, 1E-6);
  EXPECT_NEAR(grad(86, 0), -0.34733167111957E-03, 1E-6);
  EXPECT_NEAR(grad(86, 1), -0.89343547272491E-04, 1E-6);
  EXPECT_NEAR(grad(86, 2), -0.72215719508202E-04, 1E-6);
  EXPECT_NEAR(grad(87, 0), -0.32914022248560E-03, 1E-6);
  EXPECT_NEAR(grad(87, 1), -0.23658566934120E-03, 1E-6);
  EXPECT_NEAR(grad(87, 2), 0.49447256654468E-04, 1E-6);
  EXPECT_NEAR(grad(88, 0), -0.28410391854520E-03, 1E-6);
  EXPECT_NEAR(grad(88, 1), 0.79951680672349E-05, 1E-6);
  EXPECT_NEAR(grad(88, 2), 0.35456205757965E-03, 1E-6);
  EXPECT_NEAR(grad(89, 0), -0.27621949446983E-03, 1E-6);
  EXPECT_NEAR(grad(89, 1), 0.11502136251572E-03, 1E-6);
  EXPECT_NEAR(grad(89, 2), 0.18501567946618E-03, 1E-6);
  EXPECT_NEAR(grad(90, 0), -0.26374185825814E-03, 1E-6);
  EXPECT_NEAR(grad(90, 1), 0.18294741624755E-03, 1E-6);
  EXPECT_NEAR(grad(90, 2), 0.26743998237160E-03, 1E-6);
  EXPECT_NEAR(grad(91, 0), -0.28464477755201E-03, 1E-6);
  EXPECT_NEAR(grad(91, 1), 0.41422672051734E-03, 1E-6);
  EXPECT_NEAR(grad(91, 2), -0.25288171029992E-03, 1E-6);
  EXPECT_NEAR(grad(92, 0), -0.16659418114369E-03, 1E-6);
  EXPECT_NEAR(grad(92, 1), 0.38592109853208E-03, 1E-6);
  EXPECT_NEAR(grad(92, 2), -0.16891257434185E-03, 1E-6);
  EXPECT_NEAR(grad(93, 0), -0.20144134086315E-03, 1E-6);
  EXPECT_NEAR(grad(93, 1), 0.25390934516422E-03, 1E-6);
  EXPECT_NEAR(grad(93, 2), -0.17607477637070E-03, 1E-6);
  EXPECT_NEAR(grad(94, 0), 0.61956940549665E-03, 1E-6);
  EXPECT_NEAR(grad(94, 1), 0.12929713142598E-03, 1E-6);
  EXPECT_NEAR(grad(94, 2), -0.75110524825660E-03, 1E-6);
  EXPECT_NEAR(grad(95, 0), 0.35842823225359E-03, 1E-6);
  EXPECT_NEAR(grad(95, 1), 0.61152755961755E-03, 1E-6);
  EXPECT_NEAR(grad(95, 2), -0.14386526517954E-02, 1E-6);
  EXPECT_NEAR(grad(96, 0), 0.52686288364152E-03, 1E-6);
  EXPECT_NEAR(grad(96, 1), 0.83560382839702E-03, 1E-6);
  EXPECT_NEAR(grad(96, 2), -0.14849934797071E-02, 1E-6);
  EXPECT_NEAR(grad(97, 0), 0.45824215342805E-04, 1E-6);
  EXPECT_NEAR(grad(97, 1), -0.19573787818676E-04, 1E-6);
  EXPECT_NEAR(grad(97, 2), -0.13596429497376E-02, 1E-6);
  EXPECT_NEAR(grad(98, 0), 0.12387026472583E-03, 1E-6);
  EXPECT_NEAR(grad(98, 1), 0.28875242066941E-03, 1E-6);
  EXPECT_NEAR(grad(98, 2), -0.40897188139410E-03, 1E-6);
  EXPECT_NEAR(grad(99, 0), 0.26616086715217E-03, 1E-6);
  EXPECT_NEAR(grad(99, 1), 0.26632272603516E-03, 1E-6);
  EXPECT_NEAR(grad(99, 2), -0.52142481168618E-03, 1E-6);
  EXPECT_NEAR(grad(100, 0), 0.16809327333500E-03, 1E-6);
  EXPECT_NEAR(grad(100, 1), 0.29730991303690E-03, 1E-6);
  EXPECT_NEAR(grad(100, 2), -0.44517730312351E-03, 1E-6);
  EXPECT_NEAR(grad(101, 0), -0.24391893655121E-03, 1E-6);
  EXPECT_NEAR(grad(101, 1), -0.49919758320386E-04, 1E-6);
  EXPECT_NEAR(grad(101, 2), -0.14980324206555E-02, 1E-6);
  EXPECT_NEAR(grad(102, 0), -0.54985138737388E-04, 1E-6);
  EXPECT_NEAR(grad(102, 1), -0.29111611237035E-03, 1E-6);
  EXPECT_NEAR(grad(102, 2), -0.10893755016443E-02, 1E-6);
  EXPECT_NEAR(grad(103, 0), -0.29943334698998E-03, 1E-6);
  EXPECT_NEAR(grad(103, 1), -0.25751870955533E-03, 1E-6);
  EXPECT_NEAR(grad(103, 2), -0.11782161493584E-02, 1E-6);
  EXPECT_NEAR(grad(104, 0), 0.20893735291185E-03, 1E-6);
  EXPECT_NEAR(grad(104, 1), 0.68894745962174E-05, 1E-6);
  EXPECT_NEAR(grad(104, 2), -0.55234287566566E-03, 1E-6);
  EXPECT_NEAR(grad(105, 0), -0.24961979188683E-03, 1E-6);
  EXPECT_NEAR(grad(105, 1), -0.19789320532086E-03, 1E-6);
  EXPECT_NEAR(grad(105, 2), -0.61001233354147E-03, 1E-6);
  EXPECT_NEAR(grad(106, 0), 0.40346805616292E-05, 1E-6);
  EXPECT_NEAR(grad(106, 1), -0.15346977022997E-03, 1E-6);
  EXPECT_NEAR(grad(106, 2), -0.43010755490281E-03, 1E-6);
  EXPECT_NEAR(grad(107, 0), -0.23848344817731E-03, 1E-6);
  EXPECT_NEAR(grad(107, 1), -0.34174384810148E-03, 1E-6);
  EXPECT_NEAR(grad(107, 2), -0.65402035669312E-03, 1E-6);
  EXPECT_NEAR(grad(108, 0), -0.16202489500003E-03, 1E-6);
  EXPECT_NEAR(grad(108, 1), 0.61243359130732E-04, 1E-6);
  EXPECT_NEAR(grad(108, 2), -0.40939633341506E-03, 1E-6);
  EXPECT_NEAR(grad(109, 0), -0.13217743981632E-03, 1E-6);
  EXPECT_NEAR(grad(109, 1), -0.21326087667958E-03, 1E-6);
  EXPECT_NEAR(grad(109, 2), -0.13363237687104E-03, 1E-6);
  EXPECT_NEAR(grad(110, 0), -0.23758846996526E-03, 1E-6);
  EXPECT_NEAR(grad(110, 1), -0.92060837122905E-04, 1E-6);
  EXPECT_NEAR(grad(110, 2), -0.16064478190534E-03, 1E-6);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ-ABC hessian correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, NoneHessBP86CO) {
  auto hess = DispersionCorrectionCalculator::calcDispersionHessianCorrection(Options::DFT_DISPERSION_CORRECTIONS::NONE,
                                                                              systemControllerOne->getGeometry(),
                                                                              CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_EQ((int)hess.size(), (int)6);
  EXPECT_NEAR(hess[0](0, 1), 0.0, 1E-9);
  EXPECT_NEAR(hess[0](1, 1), 0.0, 1E-9);
  EXPECT_NEAR(hess[5](0, 1), 0.0, 1E-9);
  EXPECT_NEAR(hess[5](1, 1), 0.0, 1E-9);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ-ABC hessian correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJABCHessBP86CO) {
  auto hess = DispersionCorrectionCalculator::calcDispersionHessianCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJABC,
                                                                              systemControllerOne->getGeometry(),
                                                                              CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_EQ((int)hess.size(), (int)6);
  EXPECT_NEAR(hess[0](0, 1), -3.4488400414885987e-07, 1E-9);
  EXPECT_NEAR(hess[0](1, 1), 3.4488400414885987e-07, 1E-9);
  EXPECT_NEAR(hess[5](0, 1), 0.0, 1E-9);
  EXPECT_NEAR(hess[5](1, 1), 0.0, 1E-9);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3BJ-ABC hessian correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJHessBP86CO) {
  auto hess = DispersionCorrectionCalculator::calcDispersionHessianCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3BJ,
                                                                              systemControllerOne->getGeometry(),
                                                                              CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_EQ((int)hess.size(), (int)6);
  EXPECT_NEAR(hess[0](0, 1), -3.4488400414885987e-07, 1E-9);
  EXPECT_NEAR(hess[0](1, 1), 3.4488400414885987e-07, 1E-9);
  EXPECT_NEAR(hess[5](0, 1), 0.0, 1E-9);
  EXPECT_NEAR(hess[5](1, 1), 0.0, 1E-9);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3 hessian correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3HessBP86CO) {
  auto hess = DispersionCorrectionCalculator::calcDispersionHessianCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3,
                                                                              systemControllerOne->getGeometry(),
                                                                              CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_EQ((int)hess.size(), (int)6);
  EXPECT_NEAR(hess[0](0, 1), 1.3099683024583039e-06, 1E-9);
  EXPECT_NEAR(hess[0](1, 1), -1.3099683024583039e-06, 1E-9);
  EXPECT_NEAR(hess[5](0, 1), 0.0, 1E-9);
  EXPECT_NEAR(hess[5](1, 1), 0.0, 1E-9);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3ABC hessian correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3ABCHessBP86CO) {
  auto hess = DispersionCorrectionCalculator::calcDispersionHessianCorrection(Options::DFT_DISPERSION_CORRECTIONS::D3ABC,
                                                                              systemControllerOne->getGeometry(),
                                                                              CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_EQ((int)hess.size(), (int)6);
  EXPECT_NEAR(hess[0](0, 1), 1.3099683024583039e-06, 1E-9);
  EXPECT_NEAR(hess[0](1, 1), -1.3099683024583039e-06, 1E-9);
  EXPECT_NEAR(hess[5](0, 1), 0.0, 1E-9);
  EXPECT_NEAR(hess[5](1, 1), 0.0, 1E-9);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(BJ) interaction correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3BJInteractionBP86WaterDimer) {
  auto waterA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto waterB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  double naddDisperion = DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3BJ, waterA->getGeometry(), waterB->getGeometry(),
      CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(-0.00087863059626244278, naddDisperion, 1e-9);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) interaction correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, D3InteractionBP86WaterDimer) {
  auto waterA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto waterB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  double naddDisperion = DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::D3, waterA->getGeometry(), waterB->getGeometry(),
      CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(-0.00126192858239497, naddDisperion, 1e-9);
}

/**
 * @test
 * @brief Tests DispersionCorrectionCalculator.h/.cpp: D3(0) interaction correction.
 */
TEST_F(DispersionCorrectionCalculatorTest, NONEInteractionBP86WaterDimer) {
  auto waterA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto waterB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  double naddDisperion = DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection(
      Options::DFT_DISPERSION_CORRECTIONS::NONE, waterA->getGeometry(), waterB->getGeometry(),
      CompositeFunctionals::XCFUNCTIONALS::BP86);
  EXPECT_NEAR(0.0, naddDisperion, 1e-9);
}

} /* namespace Serenity */
