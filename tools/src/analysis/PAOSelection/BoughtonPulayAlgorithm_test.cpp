/**
 * @file BoughtonPulayAlgorithm_test.cpp
 *
 * @date Dec 14, 2018
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
#include "analysis/PAOSelection/BoughtonPulayAlgorithm.h"             //To be tested.
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h" //Mulliken populations.
#include "basis/AtomCenteredBasisController.h"                        //Atom<-->Shell mapping.
#include "data/OrbitalController.h"                                   //Mulliken populations.
#include "integrals/OneElectronIntegralController.h"                  //getOverlapIntegrals().
#include "system/SystemController.h"                                  //Test systems.
#include "tasks/LocalizationTask.h"                                   //Orbital localization.
#include "testsupply/SystemController__TEST_SUPPLY.h"                 //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class BoughtonPulayAlgorithmTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests the function "assignAtoms" / the Boughton--Pulay algorithm.
 */
TEST_F(BoughtonPulayAlgorithmTest, test_assignAtoms) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  LocalizationTask locTask(system);
  locTask.run();
  unsigned int nOcc = system->getNOccupiedOrbitals<scfMode>();
  const auto& coefficients = system->getActiveOrbitalController<scfMode>()->getCoefficients();
  auto mullikenGrossCharges = MullikenPopulationCalculator<scfMode>::calculateAtomwiseOrbitalPopulations(
      coefficients, system->getOneElectronIntegralController()->getOverlapIntegrals(),
      system->getAtomCenteredBasisController()->getBasisIndices());
  BoughtonPulayAlgorithm bpAlgorithm(system->getOneElectronIntegralController(), system->getAtomCenteredBasisController(),
                                     std::make_shared<Eigen::MatrixXd>(mullikenGrossCharges.leftCols(nOcc).eval()),
                                     std::make_shared<Eigen::MatrixXd>(coefficients.leftCols(nOcc).eval()));
  // The reference results.
  std::vector<std::vector<unsigned int>> referenceResults;
  referenceResults.push_back({2});
  referenceResults.push_back({2, 0});
  referenceResults.push_back({2, 1});
  referenceResults.push_back({2});
  referenceResults.push_back({2});
  // Calculate and compare.
  for (unsigned int i = 0; i < nOcc; ++i) {
    auto atomIndices = bpAlgorithm.assignAtoms(mullikenGrossCharges.col(i), coefficients.col(i));
    assert(referenceResults[i].size() == atomIndices.size());
    for (unsigned int j = 0; j < atomIndices.size(); ++j) {
      EXPECT_EQ(referenceResults[i][j], atomIndices[j]);
      std::cout << atomIndices[j] << "  ";
    }
    std::cout << std::endl;
  }
}

/**
 * @test
 * @brief Tests the function "assignAtoms" / the Boughton--Pulay algorithm.
 */
TEST_F(BoughtonPulayAlgorithmTest, test_selectPAOs) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  LocalizationTask locTask(system);
  locTask.run();
  unsigned int nOcc = system->getNOccupiedOrbitals<scfMode>();
  const auto& coefficients = system->getActiveOrbitalController<scfMode>()->getCoefficients();
  auto mullikenGrossCharges = MullikenPopulationCalculator<scfMode>::calculateAtomwiseOrbitalPopulations(
      coefficients, system->getOneElectronIntegralController()->getOverlapIntegrals(),
      system->getAtomCenteredBasisController()->getBasisIndices());
  BoughtonPulayAlgorithm bpAlgorithm(system->getOneElectronIntegralController(), system->getAtomCenteredBasisController(),
                                     std::make_shared<Eigen::MatrixXd>(mullikenGrossCharges.leftCols(nOcc).eval()),
                                     std::make_shared<Eigen::MatrixXd>(coefficients.leftCols(nOcc).eval()));
  // The reference results.
  Eigen::MatrixXi referenceResults = Eigen::MatrixXi::Zero(coefficients.rows(), nOcc);
  referenceResults.block(4, 0, 14, 1) = Eigen::VectorXi::Constant(14, 1);
  referenceResults.block(2, 2, 16, 1) = Eigen::VectorXi::Constant(16, 1);
  referenceResults.block(0, 1, 2, 1) = Eigen::VectorXi::Constant(2, 1);
  referenceResults.block(4, 1, 14, 1) = Eigen::VectorXi::Constant(14, 1);
  referenceResults.block(4, 3, 14, 1) = Eigen::VectorXi::Constant(14, 1);
  referenceResults.block(4, 4, 14, 1) = Eigen::VectorXi::Constant(14, 1);
  // Calculate and compare.
  auto map = bpAlgorithm.selectPAOs();
  int diff = Eigen::MatrixXi(referenceResults - (*map)).array().abs().sum();
  EXPECT_EQ(diff, 0);
}

} /* namespace Serenity */
