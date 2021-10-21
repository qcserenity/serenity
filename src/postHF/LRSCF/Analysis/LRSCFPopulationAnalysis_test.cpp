/**
 * @file LRSCFPopulationAnalysis_test.cpp
 *
 * @date Oct 07, 2021
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
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"

/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <fstream>

namespace Serenity {

/**
 * @class LRSCFPopulationAnalysisTest
 * @brief Sets everything up for the tests of LRSCFPopulationAnalysis.h/.cpp .
 */
class LRSCFPopulationAnalysisTest : public ::testing::Test {
 protected:
  LRSCFPopulationAnalysisTest() {
  }

  virtual ~LRSCFPopulationAnalysisTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(LRSCFPopulationAnalysisTest, TransitionCharges) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, true);
  std::vector<std::shared_ptr<SystemController>> active = {sys};
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);

  lrscf.settings.densFitK = Options::DENS_FITS::RI;
  lrscf.settings.nEigen = 4;
  lrscf.settings.transitionCharges = true;
  lrscf.run();

  // Read transition charges from disk and compare.
  unsigned nAtoms = sys->getNAtoms();
  Eigen::MatrixXd charges(nAtoms, 3 + lrscf.settings.nEigen);
  std::string name = sys->getSystemPath() + sys->getSystemName() + ".transitioncharges.txt";
  std::ifstream file;
  file.open(name);
  EXPECT_TRUE(file.is_open());

  for (unsigned i = 0; i < charges.rows(); ++i) {
    for (unsigned j = 0; j < charges.cols(); ++j) {
      file >> charges(i, j);
    }
  }
  file.close();

  // Serenity Oct 2021
  Eigen::MatrixXd reference(nAtoms, 3 + lrscf.settings.nEigen);
  reference.row(0) << 2.0193613367750980e+00, -2.6664035618469617e-01, 1.9668269505105014e+00, -1.1618141936996403e-06,
      4.0173833576396437e-02, -4.2357161818948115e-02, 1.8885086677722596e-02;
  reference.row(1) << 1.1298672499137479e+00, 2.8534864481849130e-02, 1.3001315737425298e-01, -1.8911053714498404e-08,
      -7.6044139068927387e-06, -1.7696870548632299e-01, -6.0748030605035543e-05;
  reference.row(2) << 2.3974955343127147e+00, 3.7832317015007916e-01, -1.4583016503737070e+00, 1.1195908103619517e-06,
      -4.0162632654122045e-02, -4.2342517249971560e-02, -1.8916305434613132e-02;
  reference.row(3) << -1.1262767702769589e+00, -2.8534864481849130e-02, -1.2963521214932783e-01, 6.1134436887175748e-08,
      -3.5965083674226451e-06, 2.6166838455524050e-01, 9.1966787494950129e-05;

  // The sign of the transition charges is not arbitrary so it's simply squared here for comparison.
  charges = charges.cwiseProduct(charges).eval();
  reference = reference.cwiseProduct(reference).eval();
  double maxDiff = (charges - reference).cwiseAbs().maxCoeff();

  EXPECT_LE(maxDiff, 1e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFPopulationAnalysisTest, HoleParticlePopulations) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, true);
  std::vector<std::shared_ptr<SystemController>> active = {sys};
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);

  lrscf.settings.densFitK = Options::DENS_FITS::RI;
  lrscf.settings.nEigen = 4;
  lrscf.settings.transitionCharges = true;
  lrscf.run();

  // Read hole particle correlations from disk and compare.
  unsigned nAtoms = sys->getNAtoms();
  std::ifstream file;

  Eigen::MatrixXd holepartreferences(nAtoms * 4, nAtoms);
  holepartreferences.row(0) << 0.0197715, 0.0851045, 0.0197662, 0.233197;
  holepartreferences.row(1) << 0.00330251, 0.0142153, 0.00330161, 0.0389517;
  holepartreferences.row(2) << 0.0288764, 0.124295, 0.0288686, 0.340585;
  holepartreferences.row(3) << 0.00330215, 0.0142138, 0.00330126, 0.0389475;
  holepartreferences.row(4) << 0.00282781, 0.0107719, 0.00282706, 0.0187613;
  holepartreferences.row(5) << 0.0326666, 0.124436, 0.0326579, 0.216729;
  holepartreferences.row(6) << 0.0121999, 0.0464728, 0.0121967, 0.0809413;
  holepartreferences.row(7) << 0.0326684, 0.124443, 0.0326598, 0.216741;
  holepartreferences.row(8) << 0.00695813, 0.101506, 0.00695805, 0.212471;
  holepartreferences.row(9) << 0.00157968, 0.0230445, 0.00157966, 0.0482364;
  holepartreferences.row(10) << 0.0111037, 0.161982, 0.0111036, 0.339059;
  holepartreferences.row(11) << 0.0015792, 0.0230376, 0.00157918, 0.0482219;
  holepartreferences.row(12) << 0.00690119, 0.0267252, 0.00689927, 0.0500123;
  holepartreferences.row(13) << 0.00558585, 0.0216315, 0.0055843, 0.0404802;
  holepartreferences.row(14) << 0.058153, 0.2252, 0.0581368, 0.421429;
  holepartreferences.row(15) << 0.00558426, 0.0216253, 0.0055827, 0.0404686;

  for (unsigned iExc = 0; iExc < 4; iExc++) {
    Eigen::MatrixXd correlation(nAtoms, nAtoms);
    std::string name = sys->getSystemPath() + sys->getSystemName() + ".correlation" + std::to_string(iExc + 1) + ".txt";

    file.open(name);
    EXPECT_TRUE(file.is_open());
    for (unsigned i = 0; i < nAtoms; ++i) {
      for (unsigned j = 0; j < nAtoms; ++j) {
        file >> correlation(i, j);
      }
    }
    file.close();
    std::remove(name.c_str());

    double maxHolePartDiff = (correlation - holepartreferences.block(4 * iExc, 0, nAtoms, nAtoms)).cwiseAbs().maxCoeff();
    EXPECT_LE(maxHolePartDiff, 1e-6);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}
} /* namespace Serenity */
