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
  lrscf.settings.conv = 1E-6;
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

  // The sign of the transition charges is arbitrary so it's simply squared here for comparison.
  charges = charges.cwiseProduct(charges).eval();
  reference = reference.cwiseProduct(reference).eval();
  double maxDiff = (charges - reference).cwiseAbs().maxCoeff();

  EXPECT_LE(maxDiff, 1E-6);
  EXPECT_EQ(0, std::remove((sys->getSystemPath() + sys->getSystemName() + ".correlation1.txt").c_str()));
  EXPECT_EQ(0, std::remove((sys->getSystemPath() + sys->getSystemName() + ".correlation2.txt").c_str()));
  EXPECT_EQ(0, std::remove((sys->getSystemPath() + sys->getSystemName() + ".correlation3.txt").c_str()));
  EXPECT_EQ(0, std::remove((sys->getSystemPath() + sys->getSystemName() + ".correlation4.txt").c_str()));
}

TEST_F(LRSCFPopulationAnalysisTest, HoleParticlePopulations) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, true);
  std::vector<std::shared_ptr<SystemController>> active = {sys};
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);

  lrscf.settings.densFitK = Options::DENS_FITS::RI;
  lrscf.settings.nEigen = 4;
  lrscf.settings.conv = 1E-8;
  lrscf.settings.transitionCharges = true;
  lrscf.run();

  // Read hole particle correlations from disk and compare.
  unsigned nAtoms = sys->getNAtoms();
  std::ifstream file;

  // Serenity Dec 2024
  Eigen::MatrixXd holepartreferences(nAtoms * 4, nAtoms);
  holepartreferences.row(0) << 0.0202827, 0.0883255, 0.0202772, 0.241169;
  holepartreferences.row(1) << 0.00338376, 0.0147354, 0.00338285, 0.0402344;
  holepartreferences.row(2) << 0.0294375, 0.128193, 0.0294296, 0.350025;
  holepartreferences.row(3) << 0.00338342, 0.0147339, 0.0033825, 0.0402302;

  holepartreferences.row(4) << 0.00283803, 0.0108101, 0.00283727, 0.0188236;
  holepartreferences.row(5) << 0.0326934, 0.12453, 0.0326847, 0.216844;
  holepartreferences.row(6) << 0.0122162, 0.0465319, 0.012213, 0.0810258;
  holepartreferences.row(7) << 0.0326953, 0.124537, 0.0326866, 0.216856;

  holepartreferences.row(8) << 0.00737159, 0.106047, 0.00737152, 0.223194;
  holepartreferences.row(9) << 0.00165649, 0.02383, 0.00165647, 0.0501544;
  holepartreferences.row(10) << 0.0116625, 0.167776, 0.0116624, 0.353114;
  holepartreferences.row(11) << 0.001656, 0.023823, 0.00165598, 0.0501396;

  holepartreferences.row(12) << 0.00691787, 0.0267944, 0.00691594, 0.0501287;
  holepartreferences.row(13) << 0.0055954, 0.0216722, 0.00559385, 0.0405458;
  holepartreferences.row(14) << 0.0581957, 0.225404, 0.0581795, 0.421701;
  holepartreferences.row(15) << 0.0055938, 0.021666, 0.00559225, 0.0405342;

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
    EXPECT_NEAR(maxHolePartDiff, 0, 1E-6);
  }
}
} /* namespace Serenity */
