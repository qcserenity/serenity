/**
 * @file FCIDumpFileWriter_test.cpp
 *
 * @date Feb. 12, 2024
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
#include "io/FCIDumpFileWriter.h"                     //To be tested.
#include "data/ElectronicStructure.h"                 //Get Fock matrix.
#include "settings/Settings.h"                        //Settings/basis set change.
#include "system/SystemController.h"                  //Settings/basis set change.
#include "tasks/FCIDumpFileWriterTask.h"              //To be tested.
#include "tasks/FCIDumpFileWriterTask.h"              //Settings + task.
#include "tasks/LocalizationTask.h"                   //Virtual orbital localization.
#include "tasks/TopDownStaticEmbeddingTask.h"         //Fast embedding calculation.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test supply.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Test environment.
#include <iostream>      //Check if file exists.

namespace Serenity {

class FCIDumpFileWriterTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Test FCIDumpFileWriter for H2
 */
TEST_F(FCIDumpFileWriterTest, H2) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF);
  std::vector<std::shared_ptr<SystemController>> dummy = {};
  FCIDumpFileWriterTaskSettings settings;
  FCIDumpFileWriter<RESTRICTED> writer(system, dummy, settings);
  SpinPolarizedData<RESTRICTED, std::vector<unsigned int>> range = {0, 1};
  std::stringstream out;
  writer.writeFCIDumpFile(out, range);
  std::string line;
  std::vector<std::string> referenceLines = {"&FCI NORB= 2, NELEC= 2, MS2= 0,", "ORBSYM= 1, 1,", "ISYM= 1", "&END"};
  for (unsigned int iLine = 0; iLine < 4; ++iLine) {
    std::getline(out, line);
    ASSERT_EQ(line, referenceLines[iLine]);
  }
  std::vector<std::pair<double, std::vector<unsigned int>>> referenceIntegrals = {
      {6.5787813240e-01, {1, 1, 1, 1}},  {3.4181379567e-01, {2, 2, 1, 1}},  {4.1914410857e-02, {2, 1, 2, 1}},
      {2.8394061945e-01, {2, 2, 2, 2}},  {-1.2523523185e+00, {1, 1, 0, 0}}, {0.0, {2, 1, 0, 0}},
      {-4.7484741238e-01, {2, 2, 0, 0}},
  };
  double integral = 0.0;
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int l = 0;
  unsigned int counter = 0;
  while (std::getline(out, line) && counter < referenceIntegrals.size()) {
    std::stringstream sstream(line);
    sstream >> integral;
    sstream >> i;
    sstream >> j;
    sstream >> k;
    sstream >> l;
    ASSERT_EQ(referenceIntegrals[counter].second[0], i);
    ASSERT_EQ(referenceIntegrals[counter].second[1], j);
    ASSERT_EQ(referenceIntegrals[counter].second[2], k);
    ASSERT_EQ(referenceIntegrals[counter].second[3], l);
    ASSERT_NEAR(referenceIntegrals[counter].first, integral, 1e-10);
    counter++;
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(FCIDumpFileWriterTest, getOneParticleIntegrals) {
  auto actA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_A, true);

  auto systemSettings = actA->getSettings();
  systemSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  actA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_A, systemSettings);

  auto envB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_B, true);
  auto envC = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_C, true);
  auto envD = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_D, true);
  auto envE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_E, true);
  auto envF = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_F, true);

  TopDownStaticEmbeddingTask<RESTRICTED> topDownStaticEmbeddingTask({actA, envB, envC}, {envD, envE, envF});
  topDownStaticEmbeddingTask.settings.lcSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  topDownStaticEmbeddingTask.settings.lcSettings.embeddingSettings.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
  topDownStaticEmbeddingTask.settings.lcSettings.embeddingSettings.embeddingMode =
      Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  topDownStaticEmbeddingTask.run();

  LocalizationTask localizationTask(actA);
  localizationTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  localizationTask.settings.localizeVirtuals = true;
  localizationTask.settings.splitValenceAndCore = true;
  localizationTask.settings.useEnergyCutOff = false;

  FCIDumpFileWriterTaskSettings settings;
  settings.embedding = topDownStaticEmbeddingTask.settings.lcSettings.embeddingSettings;
  FCIDumpFileWriter<RESTRICTED> writer(actA, {envB, envC, envD, envE, envF}, settings);
  SpinPolarizedData<RESTRICTED, std::vector<unsigned int>> emptyRange = {};
  SpinPolarizedData<RESTRICTED, std::vector<unsigned int>> range = {1, 2, 3, 4, 5, 6, 7, 8};
  auto fullFock = writer.getOneParticleIntegrals(emptyRange);
  auto oneParticleIntegrals = writer.getOneParticleIntegrals(range);
  auto fockReference = actA->getElectronicStructure<RESTRICTED>()->getFockMatrix();
  auto diff = (fockReference.array() - fullFock.array()).abs().maxCoeff();
  auto diff2 = (fockReference.array() - oneParticleIntegrals.array()).abs().maxCoeff();
  EXPECT_NEAR(diff, 0.0, 1e-9);
  EXPECT_GT(diff2, 1e-6);
  EXPECT_NEAR((fockReference - oneParticleIntegrals)(0, 0), 10.470462, 1e-5);
  EXPECT_NEAR(-71.993760470656355, writer.getUncorrelatedActiveSpaceEnergy(range), 1e-7);
  EXPECT_NEAR(-385.40085681159803, writer.getTotalCoreEnergy(range), 1e-9);
  EXPECT_NEAR(writer.getTotalUncorrelatedEnergy(),
              writer.getTotalCoreEnergy(range) + writer.getUncorrelatedActiveSpaceEnergy(range), 1e-9);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(envD->getSystemPath() + "TMP_Supersystem/", "TMP_Supersystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(FCIDumpFileWriterTest, CO) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  auto settings = system->getSettings();
  settings.basis.label = "def2-svp";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.scf.energyThreshold = 1e-8;
  settings.scf.diisThreshold = 1e-8;
  settings.basis.makeSphericalBasis = true;
  system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS, settings);
  FCIDumpFileWriterTask<RESTRICTED> writerTask(system, {});
  writerTask.run();
  std::ifstream inputStream(writerTask.settings.outputFilePath);
  ASSERT_TRUE(inputStream.is_open());
  std::string line;
  unsigned int iLine = 0;
  double integral = 0.0;
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int l = 0;
  while (std::getline(inputStream, line)) {
    // Skip header.
    if (iLine < 4) {
      iLine++;
      continue;
    }
    std::stringstream sstream(line);
    sstream >> integral;
    sstream >> i;
    sstream >> j;
    sstream >> k;
    sstream >> l;
    // NaN checks
    ASSERT_TRUE(i == i);
    ASSERT_TRUE(j == j);
    ASSERT_TRUE(k == k);
    ASSERT_TRUE(l == l);
    ASSERT_NEAR(integral, integral, 1e-10);
    iLine++;
  }
  inputStream.close();
  std::remove(writerTask.settings.outputFilePath.c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */