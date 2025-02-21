/**
 * @file ActiveSpaceSelectionTask_test.cpp
 *
 * @date Oct 25, 2018
 * @author: Moritz Bensberg
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
#include "tasks/ActiveSpaceSelectionTask.h"           //To be tested.
#include "data/ElectronicStructure.h"                 //getEnergy().
#include "data/OrbitalController.h"                   //Check orbital after selection.
#include "geometry/Geometry.h"                        //Check geometries after selection.
#include "settings/Settings.h"                        //Settings definition.
#include "system/SystemController.h"                  //Access to test systems.
#include "tasks/FDETask.h"                            //DOS example with FDE task.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Access to test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class ActiveSpaceSelectionTest : public ::testing::Test {
 protected:
  ActiveSpaceSelectionTest() = default;

  virtual ~ActiveSpaceSelectionTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(ActiveSpaceSelectionTest, restricted_EthaneIAOShell) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto A = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  auto B = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86);
  ActiveSpaceSelectionTask<SPIN> asTask({A, B}, {}, {});
  asTask.keepSystemPairs = true;
  asTask.settings.similarityKinEnergyThreshold = 1e-1;
  asTask.settings.similarityLocThreshold = 1e-1;
  asTask.settings.usePiBias = true;
  asTask.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell;
  asTask.run();

  auto act1 = asTask.getSystemPairs()[0].first;
  auto act2 = asTask.getSystemPairs()[1].first;
  auto env1 = asTask.getSystemPairs()[0].second;
  auto env2 = asTask.getSystemPairs()[1].second;
  auto coeffActA = act1->getActiveOrbitalController<SPIN>()->getCoefficients();
  // there should only be one occupied orbital
  auto occOrb = coeffActA.col(0);
  Eigen::VectorXd compareTo(occOrb.size());
  // This is the C-C sigma bond orbital which should be considered as active.
  compareTo << 7.876884977852e-02, -2.708802357799e-01, -8.405365955012e-02, -2.270559333949e-01, 7.225321806091e-07,
      -2.958200761817e-01, -9.723283453681e-02, -6.910491103526e-07, -1.266759020653e-01, -1.263795182561e-02,
      4.034628624195e-08, 7.553140816271e-03, -1.472047562465e-07, -3.382947885986e-03, 7.876846333397e-02,
      -2.708792898410e-01, -8.405178496326e-02, 2.270555414778e-01, 1.191516736860e-07, 2.958207390581e-01,
      9.723253610694e-02, 1.213875286177e-06, 1.266757972472e-01, -1.263803088525e-02, 1.653802189768e-07,
      7.553206793156e-03, -2.436304510273e-07, -3.382752771861e-03, 2.132082251802e-02, 3.369489848447e-02,
      -4.803680454666e-03, 1.193864133833e-03, -9.238202012026e-03, 2.131968466476e-02, 3.369254801949e-02,
      -7.782931985677e-03, 9.746564345522e-04, -6.951857512146e-03, 2.131950271331e-02, 3.369308420169e-02,
      -6.142299483056e-03, -2.168109818022e-03, -8.210958311094e-03, 2.131911294640e-02, 3.369317389085e-02,
      6.142395067744e-03, 2.168240579504e-03, 8.211040070319e-03, 2.131935479861e-02, 3.369259397654e-02,
      7.782745648485e-03, -9.745447971570e-04, 6.951822838524e-03, 2.132079030230e-02, 3.369465227585e-02,
      4.803708752552e-03, -1.193967682197e-03, 9.238472315750e-03;
  EXPECT_NEAR((occOrb - compareTo).array().abs().sum(), 0.0, 1e-6);
  unsigned int nAtoms = A->getGeometry()->getNAtoms();
  EXPECT_EQ(act1->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(env1->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(act2->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(env2->getGeometry()->getNAtoms(), nAtoms);

  auto C = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act);
  Settings cSettings = C->getSettings();
  cSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  C = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act, cSettings);

  ActiveSpaceSelectionTask<SPIN> asTask2({A, B}, {act1, C}, {env1, env2});
  asTask2.settings.similarityKinEnergyThreshold = 8e-2;
  asTask2.settings.similarityLocThreshold = 8e-2;
  asTask2.settings.load = true;
  asTask2.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell;
  asTask2.run();
  coeffActA = act1->getActiveOrbitalController<SPIN>()->getCoefficients();
  occOrb = coeffActA.col(0);
  EXPECT_NEAR((occOrb - compareTo).array().abs().sum(), 0.0, 1e-6);
  EXPECT_TRUE(C->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[0].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[0].second);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].second);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ActiveSpaceSelectionTest, restricted_EthaneIAO) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto A = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  auto B = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86);
  ActiveSpaceSelectionTask<SPIN> asTask({A, B}, {}, {});
  asTask.keepSystemPairs = true;
  asTask.settings.similarityKinEnergyThreshold = 5e-2;
  asTask.settings.similarityLocThreshold = 5e-2;
  asTask.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAO;
  asTask.run();

  auto act1 = asTask.getSystemPairs()[0].first;
  auto act2 = asTask.getSystemPairs()[1].first;
  auto env1 = asTask.getSystemPairs()[0].second;
  auto env2 = asTask.getSystemPairs()[1].second;
  auto coeffActA = act1->getActiveOrbitalController<SPIN>()->getCoefficients();
  // there should only be one occupied orbital
  auto occOrb = coeffActA.col(0);
  Eigen::VectorXd compareTo(occOrb.size());
  // This is the C-C sigma bond orbital which should be considered as active.
  compareTo << 7.876884977852e-02, -2.708802357799e-01, -8.405365955012e-02, -2.270559333949e-01, 7.225321806091e-07,
      -2.958200761817e-01, -9.723283453681e-02, -6.910491103526e-07, -1.266759020653e-01, -1.263795182561e-02,
      4.034628624195e-08, 7.553140816271e-03, -1.472047562465e-07, -3.382947885986e-03, 7.876846333397e-02,
      -2.708792898410e-01, -8.405178496326e-02, 2.270555414778e-01, 1.191516736860e-07, 2.958207390581e-01,
      9.723253610694e-02, 1.213875286177e-06, 1.266757972472e-01, -1.263803088525e-02, 1.653802189768e-07,
      7.553206793156e-03, -2.436304510273e-07, -3.382752771861e-03, 2.132082251802e-02, 3.369489848447e-02,
      -4.803680454666e-03, 1.193864133833e-03, -9.238202012026e-03, 2.131968466476e-02, 3.369254801949e-02,
      -7.782931985677e-03, 9.746564345522e-04, -6.951857512146e-03, 2.131950271331e-02, 3.369308420169e-02,
      -6.142299483056e-03, -2.168109818022e-03, -8.210958311094e-03, 2.131911294640e-02, 3.369317389085e-02,
      6.142395067744e-03, 2.168240579504e-03, 8.211040070319e-03, 2.131935479861e-02, 3.369259397654e-02,
      7.782745648485e-03, -9.745447971570e-04, 6.951822838524e-03, 2.132079030230e-02, 3.369465227585e-02,
      4.803708752552e-03, -1.193967682197e-03, 9.238472315750e-03;
  EXPECT_NEAR((occOrb - compareTo).array().abs().sum(), 0.0, 1e-6);
  unsigned int nAtoms = A->getGeometry()->getNAtoms();
  EXPECT_EQ(act1->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(env1->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(act2->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(env2->getGeometry()->getNAtoms(), nAtoms);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[0].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[0].second);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].second);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ActiveSpaceSelectionTest, restricted_EthaneMulliken) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto A = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  auto B = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86);
  ActiveSpaceSelectionTask<SPIN> asTask({A, B}, {}, {});
  asTask.keepSystemPairs = true;
  asTask.settings.similarityKinEnergyThreshold = 5e-2;
  asTask.settings.similarityLocThreshold = 5e-2;
  asTask.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN;
  asTask.run();

  auto act1 = asTask.getSystemPairs()[0].first;
  auto act2 = asTask.getSystemPairs()[1].first;
  auto env1 = asTask.getSystemPairs()[0].second;
  auto env2 = asTask.getSystemPairs()[1].second;
  auto coeffActA = act1->getActiveOrbitalController<SPIN>()->getCoefficients();
  // there should only be one occupied orbital
  auto occOrb = coeffActA.col(0);
  Eigen::VectorXd compareTo(occOrb.size());
  // This is the C-C sigma bond orbital which should be considered as active.
  compareTo << 7.876884977852e-02, -2.708802357799e-01, -8.405365955012e-02, -2.270559333949e-01, 7.225321806091e-07,
      -2.958200761817e-01, -9.723283453681e-02, -6.910491103526e-07, -1.266759020653e-01, -1.263795182561e-02,
      4.034628624195e-08, 7.553140816271e-03, -1.472047562465e-07, -3.382947885986e-03, 7.876846333397e-02,
      -2.708792898410e-01, -8.405178496326e-02, 2.270555414778e-01, 1.191516736860e-07, 2.958207390581e-01,
      9.723253610694e-02, 1.213875286177e-06, 1.266757972472e-01, -1.263803088525e-02, 1.653802189768e-07,
      7.553206793156e-03, -2.436304510273e-07, -3.382752771861e-03, 2.132082251802e-02, 3.369489848447e-02,
      -4.803680454666e-03, 1.193864133833e-03, -9.238202012026e-03, 2.131968466476e-02, 3.369254801949e-02,
      -7.782931985677e-03, 9.746564345522e-04, -6.951857512146e-03, 2.131950271331e-02, 3.369308420169e-02,
      -6.142299483056e-03, -2.168109818022e-03, -8.210958311094e-03, 2.131911294640e-02, 3.369317389085e-02,
      6.142395067744e-03, 2.168240579504e-03, 8.211040070319e-03, 2.131935479861e-02, 3.369259397654e-02,
      7.782745648485e-03, -9.745447971570e-04, 6.951822838524e-03, 2.132079030230e-02, 3.369465227585e-02,
      4.803708752552e-03, -1.193967682197e-03, 9.238472315750e-03;
  EXPECT_NEAR((occOrb - compareTo).array().abs().sum(), 0.0, 1e-6);
  unsigned int nAtoms = A->getGeometry()->getNAtoms();
  EXPECT_EQ(act1->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(env1->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(act2->getGeometry()->getNAtoms(), nAtoms);
  EXPECT_EQ(env2->getGeometry()->getNAtoms(), nAtoms);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[0].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[0].second);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].second);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ActiveSpaceSelectionTest, full_DOS_run) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto A = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  auto actA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act);
  auto envA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Env);
  auto B = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86);
  ActiveSpaceSelectionTask<SPIN> asTask({A, B}, {actA}, {envA});
  asTask.keepSystemPairs = true;
  asTask.settings.similarityKinEnergyThreshold = 5e-2;
  asTask.settings.similarityLocThreshold = 5e-2;
  asTask.settings.alignPiOrbitals = true;
  asTask.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAO;
  asTask.run();

  FDETask<SPIN> fdeTask(actA, {envA});
  fdeTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  fdeTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  fdeTask.settings.calculateEnvironmentEnergy = true;
  fdeTask.run();
  double supersystemEnergy = A->getElectronicStructure<SPIN>()->getEnergy();
  double embeddedEnergy = actA->getElectronicStructure<SPIN>()->getEnergy();
  EXPECT_NEAR(supersystemEnergy, embeddedEnergy, 5e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].second);
  auto supersystemName = actA->getSystemName() + "+" + envA->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(envA->getSystemPath() + supersystemName + "/", supersystemName);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ActiveSpaceSelectionTest, unrestricted_Ethane) {
  const auto SPIN = Options::SCF_MODES::UNRESTRICTED;
  auto A = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  auto B = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86);
  A->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  B->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  ActiveSpaceSelectionTask<SPIN> asTask({A, B}, {}, {});
  asTask.keepSystemPairs = true;
  asTask.settings.similarityKinEnergyThreshold = 5e-2;
  asTask.settings.similarityLocThreshold = 5e-2;
  asTask.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAO;
  asTask.run();
  auto activeA = asTask.getSystemPairs()[0].first;
  auto coeffActA = activeA->getActiveOrbitalController<SPIN>()->getCoefficients();
  // there should only be one occupied orbital
  auto occOrb = coeffActA.alpha.col(0);
  Eigen::VectorXd compareTo(occOrb.size());
  // This is the C-C sigma bond orbital which should be considered as active.
  compareTo << 7.876884977852e-02, -2.708802357799e-01, -8.405365955012e-02, -2.270559333949e-01, 7.225321806091e-07,
      -2.958200761817e-01, -9.723283453681e-02, -6.910491103526e-07, -1.266759020653e-01, -1.263795182561e-02,
      4.034628624195e-08, 7.553140816271e-03, -1.472047562465e-07, -3.382947885986e-03, 7.876846333397e-02,
      -2.708792898410e-01, -8.405178496326e-02, 2.270555414778e-01, 1.191516736860e-07, 2.958207390581e-01,
      9.723253610694e-02, 1.213875286177e-06, 1.266757972472e-01, -1.263803088525e-02, 1.653802189768e-07,
      7.553206793156e-03, -2.436304510273e-07, -3.382752771861e-03, 2.132082251802e-02, 3.369489848447e-02,
      -4.803680454666e-03, 1.193864133833e-03, -9.238202012026e-03, 2.131968466476e-02, 3.369254801949e-02,
      -7.782931985677e-03, 9.746564345522e-04, -6.951857512146e-03, 2.131950271331e-02, 3.369308420169e-02,
      -6.142299483056e-03, -2.168109818022e-03, -8.210958311094e-03, 2.131911294640e-02, 3.369317389085e-02,
      6.142395067744e-03, 2.168240579504e-03, 8.211040070319e-03, 2.131935479861e-02, 3.369259397654e-02,
      7.782745648485e-03, -9.745447971570e-04, 6.951822838524e-03, 2.132079030230e-02, 3.369465227585e-02,
      4.803708752552e-03, -1.193967682197e-03, 9.238472315750e-03;
  EXPECT_NEAR((occOrb - compareTo).array().abs().sum(), 0.0, 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[0].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[0].second);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].first);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(asTask.getSystemPairs()[1].second);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ActiveSpaceSelectionTest, testThrows) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto A = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, true);
  auto B = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, true);
  auto C = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  Settings settingsB = B->getSettings();
  settingsB.charge = +2;
  // Only one supersystem.
  EXPECT_THROW(ActiveSpaceSelectionTask<SPIN> asTask({A}, {}, {}), SerenityError);
  // Trying to load an non-existing electronic structure.
  ActiveSpaceSelectionTask<SPIN> asTask({A, B}, {}, {});
  asTask.settings.load = true;
  EXPECT_THROW(asTask.run(), SerenityError);
  // Wrong population algorithm.
  ActiveSpaceSelectionTask<SPIN> asTask2({A, B}, {}, {});
  asTask2.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::HIRSHFELD;
  EXPECT_THROW(asTask2.run(), SerenityError);
  ActiveSpaceSelectionTask<SPIN> asTask3({A, C}, {}, {});
  EXPECT_THROW(asTask3.run(), SerenityError);
  settingsB.charge = 0;
  settingsB.basis.label = "DEF2-TZVP";
  auto D = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, settingsB);
  // Wrong supersystem II basis label.
  ActiveSpaceSelectionTask<SPIN> asTask4({A, D}, {}, {});
  EXPECT_THROW(asTask4.run(), SerenityError);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("TestSystem_EthaneA_Def2_SVP_BP86_Act/",
                                                        "TestSystem_EthaneA_Def2_SVP_BP86_Act");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("TestSystem_EthaneA_Def2_SVP_BP86_Env/",
                                                        "TestSystem_EthaneA_Def2_SVP_BP86_Env");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("TestSystem_EthaneB_Def2_SVP_BP86_Act/",
                                                        "TestSystem_EthaneB_Def2_SVP_BP86_Act");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("TestSystem_EthaneB_Def2_SVP_BP86_Env/",
                                                        "TestSystem_EthaneB_Def2_SVP_BP86_Env");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("TestSystem_WaterMonOne_Def2_SVP_Act/",
                                                        "TestSystem_WaterMonOne_Def2_SVP_Act");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("TestSystem_WaterMonOne_Def2_SVP_Env/",
                                                        "TestSystem_WaterMonOne_Def2_SVP_Env");
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
