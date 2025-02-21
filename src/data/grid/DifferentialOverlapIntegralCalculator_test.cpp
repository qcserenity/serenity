/**
 * @file DifferentialOverlapIntegralCalculator_test.cpp
 *
 * @date Jan 31, 2019
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
#include "data/grid/DifferentialOverlapIntegralCalculator.h" //To be tested.
#include "basis/BasisController.h"                           //Basis function values on a grid.
#include "data/ElectronicStructure.h"                        //Run scfs.
#include "data/OrbitalController.h"                          //Orbital coefficients.
#include "data/grid/BasisFunctionOnGridControllerFactory.h"  //Needed for DOI calculation.
#include "system/SystemController.h"                         //Test systems.
#include "testsupply/SystemController__TEST_SUPPLY.h"        //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class DifferentialOverlapIntegralCalculatorTest : public ::testing::Test {
 protected:
  DifferentialOverlapIntegralCalculatorTest() {
  }
  ~DifferentialOverlapIntegralCalculatorTest() = default;
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(DifferentialOverlapIntegralCalculatorTest, dois_BasisFunctionVSBasisFunction) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::C60_MINBAS);
  unsigned int nBasisFunctions = act->getBasisController()->getNBasisFunctions();
  auto basFuncOnGridController =
      BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(), act->getGridController());
  Eigen::VectorXd c_x = Eigen::VectorXd::Zero(nBasisFunctions);
  Eigen::VectorXd c_y = Eigen::VectorXd::Zero(nBasisFunctions);
  c_x(0) = 1.0;
  c_y(0) = 1.0;
  Eigen::VectorXi dummyMap = Eigen::VectorXi::Zero(nBasisFunctions);
  dummyMap(0) = 1;
  Eigen::SparseMatrix<int> basisFunctionToXMap = dummyMap.sparseView();
  Eigen::MatrixXd dois;
  DifferentialOverlapIntegralCalculator::calculateDOI(c_x, c_y, basisFunctionToXMap, basisFunctionToXMap,
                                                      basFuncOnGridController, dois);
  EXPECT_NEAR(dois(0, 0), 2.6930610457882125, 1e-5);

  Eigen::MatrixXd dois2;
  DifferentialOverlapIntegralCalculator::calculateDOI(basFuncOnGridController, dois2);
  EXPECT_NEAR(dois(0, 0), dois2(0, 0), 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DifferentialOverlapIntegralCalculatorTest, dois_WaterOccupiedOrbitals) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  unsigned int nBasisFunctions = act->getBasisController()->getNBasisFunctions();
  auto basFuncOnGridController =
      BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(), act->getGridController());
  double e = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  (void)e;
  const auto& coef = act->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  unsigned int nOcc = act->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  Eigen::MatrixXd c_x = coef.leftCols(nOcc).eval();
  std::vector<Eigen::Triplet<int>> tripletList;
  for (unsigned int col = 0; col < nOcc; ++col) {
    for (unsigned int row = 0; row < c_x.rows(); ++row) {
      if (c_x(row, col) * c_x(row, col) > 1e-5)
        tripletList.push_back(Eigen::Triplet<int>(row, col, 1));
    }
  }
  Eigen::SparseMatrix<int> basisFunctionToXMap(nBasisFunctions, nOcc);
  basisFunctionToXMap.setFromTriplets(tripletList.begin(), tripletList.end());

  Eigen::MatrixXd dois;
  DifferentialOverlapIntegralCalculator::calculateDOI(c_x, c_x, basisFunctionToXMap, basisFunctionToXMap,
                                                      basFuncOnGridController, dois);
  /*
   * I have created this data with a working version of the code.
   * The functionality of the code was checked by calculating the overlap matrix (drop the square for x and y).
   * The results were accurate up to a precision of ~ 1e-5 to 1e-6.
   */
  EXPECT_NEAR(dois(0, 0), 4.1724123596531957, 1e-4);
  EXPECT_NEAR(dois(1, 0), 8.360601159692e-01, 1e-4);
  EXPECT_NEAR(dois(2, 0), 0.26167785576210734, 1e-4);
  EXPECT_NEAR(dois(3, 0), 0.40303050173570215, 1e-4);
  EXPECT_NEAR(dois(4, 0), 0.33241201125842823, 1e-4);
  EXPECT_NEAR(dois(1, 1), 0.28794496068802067, 1e-4);
  EXPECT_NEAR(dois(2, 2), 2.430418783570e-01, 1e-4);
  EXPECT_NEAR(dois(3, 3), 2.868143703955e-01, 1e-4);
  EXPECT_NEAR(dois(4, 4), 3.080575761869e-01, 1e-4);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
