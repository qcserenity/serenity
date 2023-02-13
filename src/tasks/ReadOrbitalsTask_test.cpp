/**
 * @file ReadOrbitalsTask_test.cpp
 *
 * @date   Dec 15, 2020
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

#include "tasks/ReadOrbitalsTask.h"                   //To be tested.
#include "data/ElectronicStructure.h"                 //GetEnergy.
#include "geometry/Atom.h"                            //Manually construct geometry.
#include "geometry/AtomTypeFactory.h"                 //Manually construct geometry.
#include "geometry/Geometry.h"                        //Manually construct geometry.
#include "potentials/bundles/PotentialBundle.h"       //getFockMatrix()
#include "settings/Settings.h"                        //Settings.
#include "system/SystemController.h"                  //Test systems.
#include "tasks/ScfTask.h"                            //Run a SCF.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
#include <string>        //SERENITY_RESOURCES

namespace Serenity {
class ReadOrbitalsTaskTest : public ::testing::Test {
 protected:
  ReadOrbitalsTaskTest() {
  }
  virtual ~ReadOrbitalsTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(ReadOrbitalsTaskTest, restricted_molpro_h2o_qzvp) {
  auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 1.8588136904 * ANGSTROM_TO_BOHR, 0.0, 0.0);
  auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.4432381548 * ANGSTROM_TO_BOHR,
                                   0.1743648456 * ANGSTROM_TO_BOHR, -0.7451535006 * ANGSTROM_TO_BOHR);
  auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.4432381548 * ANGSTROM_TO_BOHR,
                                   -0.1743648456 * ANGSTROM_TO_BOHR, 0.7451535006 * ANGSTROM_TO_BOHR);

  auto geom = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H1, H2});
  Settings settings;
  settings.name = "water_qzvp";
  settings.basis.label = "DEF2-QZVP";
  settings.basis.integralIncrementThresholdStart = 1e-12;

  auto sys = std::make_shared<SystemController>(geom, settings);
  ReadOrbitalsTask<RESTRICTED> read(sys);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/molpro_orbitals/h2o_orbitals_hf_qzvp.xml";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::MOLPRO;
  read.settings.path = pathToTestsResources;
  read.run();

  auto eStruc = sys->getElectronicStructure<RESTRICTED>();
  auto dMat = eStruc->getDensityMatrix();
  auto eCont = eStruc->getEnergyComponentController();
  auto bundle = sys->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto F = bundle->getFockMatrix(dMat, eCont);
  double energy_without_scf = eStruc->getEnergy();

  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  double energy_aftter_scf = eStruc->getEnergy();
  EXPECT_NEAR(energy_without_scf, energy_aftter_scf, 1e-7);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("./water_qzvp/", "water_qzvp");
}

TEST_F(ReadOrbitalsTaskTest, restricted_molcas_h2o_qzvp) {
  auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 0.0);
  auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.4 * ANGSTROM_TO_BOHR, 0.0 * ANGSTROM_TO_BOHR,
                                   1.4 * ANGSTROM_TO_BOHR);
  auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.4 * ANGSTROM_TO_BOHR, 0.0 * ANGSTROM_TO_BOHR,
                                   1.4 * ANGSTROM_TO_BOHR);

  auto geom = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H1, H2});
  Settings settings;
  settings.name = "water";
  settings.basis.label = "DEF2-QZVP";
  settings.basis.integralIncrementThresholdStart = 1e-12;

  auto sys = std::make_shared<SystemController>(geom, settings);
  ReadOrbitalsTask<RESTRICTED> read(sys);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/molcas_orbitals/";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::MOLCAS;
  read.settings.path = pathToTestsResources;
  read.run();

  auto eStruc = sys->getElectronicStructure<RESTRICTED>();
  auto dMat = eStruc->getDensityMatrix();
  auto eCont = eStruc->getEnergyComponentController();
  auto bundle = sys->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto F = bundle->getFockMatrix(dMat, eCont);
  double energy_without_scf = eStruc->getEnergy();

  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  double energy_aftter_scf = eStruc->getEnergy();
  EXPECT_NEAR(energy_without_scf, energy_aftter_scf, 1e-7);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("./water/", "water");
}

TEST_F(ReadOrbitalsTaskTest, restricted_molcas_h2o_qzvp_replace) {
  auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 0.0);
  auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.4 * ANGSTROM_TO_BOHR, 0.0 * ANGSTROM_TO_BOHR,
                                   1.4 * ANGSTROM_TO_BOHR);
  auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.4 * ANGSTROM_TO_BOHR, 0.0 * ANGSTROM_TO_BOHR,
                                   1.4 * ANGSTROM_TO_BOHR);

  auto geom = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H1, H2});
  Settings settings;
  settings.name = "water";
  settings.basis.label = "DEF2-QZVP";
  settings.basis.integralIncrementThresholdStart = 1e-12;

  auto sys = std::make_shared<SystemController>(geom, settings);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/molcas_orbitals/";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  ReadOrbitalsTask<RESTRICTED> replace(sys);
  replace.settings.fileFormat = Options::ORBITAL_FILE_TYPES::MOLCAS;
  replace.settings.path = pathToTestsResources;
  replace.settings.replaceInFile = true;
  replace.run();

  // Copy the replaced orbitals to read them later again.
  const std::string originalFilePath = pathToTestsResources + "/" + sys->getSystemName() + "_backup.scf.h5";
  const std::string destinationFilePath = sys->getSystemName() + ".scf.h5";
  std::ifstream src(originalFilePath, std::ios::binary);
  std::ofstream dst(destinationFilePath, std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  ReadOrbitalsTask<RESTRICTED> read(sys);
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::MOLCAS;
  read.settings.path = ".";
  read.run();

  auto eStruc = sys->getElectronicStructure<RESTRICTED>();
  auto dMat = eStruc->getDensityMatrix();
  auto eCont = eStruc->getEnergyComponentController();
  auto bundle = sys->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto F = bundle->getFockMatrix(dMat, eCont);
  double energy_without_scf = eStruc->getEnergy();

  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  double energy_aftter_scf = eStruc->getEnergy();
  EXPECT_NEAR(energy_without_scf, energy_aftter_scf, 1e-9);

  std::remove(destinationFilePath.c_str());
  std::remove(originalFilePath.c_str());
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("./water/", "water");
}

TEST_F(ReadOrbitalsTaskTest, restricted_turbomole_h2o_qzvp) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.label = "DEF2-QZVP";
  settings.basis.incrementalSteps = 1;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);
  ReadOrbitalsTask<RESTRICTED> read(sys);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/" + settings.name;
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::TURBOMOLE;
  read.settings.path = pathToTestsResources;
  read.run();

  auto eStruc = sys->getElectronicStructure<RESTRICTED>();
  auto dMat = eStruc->getDensityMatrix();
  auto eCont = eStruc->getEnergyComponentController();
  auto bundle = sys->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto F = bundle->getFockMatrix(dMat, eCont);
  double energy_without_scf = eStruc->getEnergy();

  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  double energy_aftter_scf = eStruc->getEnergy();
  EXPECT_NEAR(energy_without_scf, energy_aftter_scf, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ReadOrbitalsTaskTest, unrestricted_turbomole_h2o_qzvp) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.label = "DEF2-QZVP";
  settings.basis.incrementalSteps = 1;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);
  ReadOrbitalsTask<UNRESTRICTED> read(sys);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/" + settings.name;
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::TURBOMOLE;
  read.settings.path = pathToTestsResources;
  read.run();

  auto eStruc = sys->getElectronicStructure<UNRESTRICTED>();
  auto dMat = eStruc->getDensityMatrix();
  auto eCont = eStruc->getEnergyComponentController();
  auto bundle = sys->getPotentials<UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto F = bundle->getFockMatrix(dMat, eCont);
  double energy_without_scf = eStruc->getEnergy();

  ScfTask<UNRESTRICTED> scf(sys);
  scf.run();
  double energy_aftter_scf = eStruc->getEnergy();
  EXPECT_NEAR(energy_without_scf, energy_aftter_scf, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ReadOrbitalsTaskTest, unrestricted_turbomole_methyl_qzvp) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.label = "DEF2-QZVP";
  settings.basis.incrementalSteps = 1;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  ReadOrbitalsTask<UNRESTRICTED> read(sys);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/" + settings.name;
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::TURBOMOLE;
  read.settings.path = pathToTestsResources;
  read.run();

  auto eStruc = sys->getElectronicStructure<UNRESTRICTED>();
  auto dMat = eStruc->getDensityMatrix();
  auto eCont = eStruc->getEnergyComponentController();
  auto bundle = sys->getPotentials<UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto F = bundle->getFockMatrix(dMat, eCont);
  double energy_without_scf = eStruc->getEnergy();

  ScfTask<UNRESTRICTED> scf(sys);
  scf.run();
  double energy_aftter_scf = eStruc->getEnergy();
  EXPECT_NEAR(energy_without_scf, energy_aftter_scf, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ReadOrbitalsTaskTest, restricted_serenity_h2o_qzvp) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.label = "DEF2-QZVP";
  settings.basis.incrementalSteps = 1;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);
  ReadOrbitalsTask<RESTRICTED> read(sys);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/" + settings.name + "/def2-QZVP/";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::SERENITY;
  read.settings.path = pathToTestsResources;
  read.run();

  auto eStruc = sys->getElectronicStructure<RESTRICTED>();
  auto dMat = eStruc->getDensityMatrix();
  auto eCont = eStruc->getEnergyComponentController();
  auto bundle = sys->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto F = bundle->getFockMatrix(dMat, eCont);
  double energy_without_scf = eStruc->getEnergy();

  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  double energy_aftter_scf = eStruc->getEnergy();
  EXPECT_NEAR(energy_without_scf, energy_aftter_scf, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ReadOrbitalsTaskTest, unrestricted_serenity_h2o_qzvp) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.label = "DEF2-QZVP";
  settings.basis.incrementalSteps = 1;
  settings.scfMode = UNRESTRICTED;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);
  ReadOrbitalsTask<UNRESTRICTED> read(sys);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/" + settings.name + "/def2-QZVP/";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::SERENITY;
  read.settings.path = pathToTestsResources;
  read.settings.resetCoreOrbitals = true;
  read.run();

  auto eStruc = sys->getElectronicStructure<UNRESTRICTED>();
  auto dMat = eStruc->getDensityMatrix();
  auto eCont = eStruc->getEnergyComponentController();
  auto bundle = sys->getPotentials<UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto F = bundle->getFockMatrix(dMat, eCont);
  double energy_without_scf = eStruc->getEnergy();

  ScfTask<UNRESTRICTED> scf(sys);
  scf.run();
  double energy_aftter_scf = eStruc->getEnergy();
  EXPECT_NEAR(energy_without_scf, energy_aftter_scf, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
