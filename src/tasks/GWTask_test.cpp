/**
 * @file   GWTask_test.cpp
 *
 * @date   10.04.2020
 * @author J. Toelle
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
#include "tasks/GWTask.h"
#include "data/ElectronicStructure.h"
#include "io/HDF5.h"
#include "parameters/Constants.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/ScfTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class GWTaskTest : public ::testing::Test {
 protected:
  GWTaskTest() {
  }

  virtual ~GWTaskTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

// dRPA tests
TEST_F(GWTaskTest, RPA_R) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-76.27245320915, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform RPA calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.mbpttype = Options::MBPT::RPA;
  gw.run();
  double correlationEnergy = gw._correlationEnergy;
  // Turbomole 7.4.1 RPA/ri approximation/(gridpoints 60)
  EXPECT_NEAR(correlationEnergy, -0.30931089941, 1.0e-6);
}

TEST_F(GWTaskTest, RPA_U) {
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  // Perform SCF
  auto scf = ScfTask<UNRESTRICTED>(a);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-39.73968410648, a->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 7e-6);
  // Perform RPA calculation
  auto gw = GWTask<UNRESTRICTED>({a});
  gw.settings.mbpttype = Options::MBPT::RPA;
  gw.run();
  double correlationEnergy = gw._correlationEnergy;
  // Turbomole 7.4.1 RPA/ri approximation/(gridpoints 60)
  EXPECT_NEAR(correlationEnergy, -0.24393318612, 2.0e-6);
}

// dRPA tests
TEST_F(GWTaskTest, RPA_R_NAF) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-76.27245320915, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform RPA calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.mbpttype = Options::MBPT::RPA;
  gw.settings.nafThresh = 1e-2;
  gw.run();
  double correlationEnergy = gw._correlationEnergy;
  // Turbomole 7.4.1 RPA/ri approximation/(gridpoints 60)
  EXPECT_NEAR(correlationEnergy, -0.30931089941, 2.0e-4);
}

TEST_F(GWTaskTest, RPA_U_NAF) {
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  // Perform SCF
  auto scf = ScfTask<UNRESTRICTED>(a);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-39.73968410648, a->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 7e-6);
  // Perform RPA calculation
  auto gw = GWTask<UNRESTRICTED>({a});
  gw.settings.mbpttype = Options::MBPT::RPA;
  gw.settings.nafThresh = 1e-2;
  gw.run();
  double correlationEnergy = gw._correlationEnergy;
  // Turbomole 7.4.1 RPA/ri approximation/(gridpoints 60)
  EXPECT_NEAR(correlationEnergy, -0.24393318612, 2.0e-4);
}

// GW tests
TEST_F(GWTaskTest, GW_analytic_R) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  settings.scf.energyThreshold = 1e-8;
  settings.scf.diisThreshold = 1e-7;
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-76.27245320915, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform RPA calculation
  auto rpa = LRSCFTask<RESTRICTED>({waterMon});
  rpa.settings.func = CompositeFunctionals::XCFUNCTIONALS::HARTREE;
  rpa.settings.nEigen = 95;
  rpa.settings.analysis = false;
  rpa.run();
  // Perform analytic GW calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.gwtype = Options::GWALGORITHM::ANALYTIC;
  gw.settings.nVirt = 20;
  gw.settings.nOcc = 20;
  gw.settings.linearized = true;
  gw.settings.freq = {0.1, 0.2, 0.1};
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  // Correlation energies from Serenity with LRSCF RPA calculation
  Eigen::VectorXd correlation2(24);
  correlation2 << 17.78944, 8.56405, 0.62563, 1.41054, 1.71895, -0.48623, -0.57484, -1.47433, -1.97661, -1.94688,
      -1.12436, -2.05481, -1.54011, -1.72087, -2.14425, -1.93246, 16.12779, -3.95479, -0.98242, 8.72794, -4.55487,
      5.71577, 2.93074, -0.85810;

  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergy_spin.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlation2(i), 5.0e-6);
    }
  };
  std::string fileName = "FrequencyDependentSelfEnergy_res.txt";
  std::remove(fileName.c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
TEST_F(GWTaskTest, GW_analytic_U) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = system->getSettings();
  settings.grid.accuracy = 5;
  settings.grid.smallGridAccuracy = 3;
  system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  // Perform SCF
  auto scf = ScfTask<UNRESTRICTED>(system);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-39.73968410168, system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  // Perform RPA calculation
  auto rpa = LRSCFTask<UNRESTRICTED>({system});
  rpa.settings.func = CompositeFunctionals::XCFUNCTIONALS::HARTREE;
  rpa.settings.nEigen = 10000;
  rpa.settings.analysis = false;
  rpa.run();
  // Perform analytic GW calculation
  auto gw = GWTask<UNRESTRICTED>({system});
  gw.settings.gwtype = Options::GWALGORITHM::ANALYTIC;
  gw.settings.nVirt = 2;
  gw.settings.nOcc = 2;
  gw.run();
  // Turbomole 7.4.1
  Eigen::VectorXd correlationEnergyTurbo_eV(8);
  correlationEnergyTurbo_eV << 0.48727794, 0.48620304, 0.89457370, -0.56690728, 0.32158448, 0.32049737, -2.29496116, -0.76934675;
  auto correlationEnergy = gw._correlationEnergies;
  unsigned counter = 0;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergy_spin.size(); i++) {
      if (std::abs(correlationEnergy_spin(i)) > 1e-8) {
        EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(counter), 2.0e-5);
        counter++;
      }
    }
  };
  SystemController__TEST_SUPPLY::cleanUp();
}
TEST_F(GWTaskTest, GW_CD_R) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-76.27245320915, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 5;
  gw.settings.nOcc = 5;
  // set to zero because closer to the turbomole results
  gw.settings.eta = 0.0;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  // Turbomole 7.4.1 gw/ri approximation/analytic/eta 0.001(but seems to be done with eta =0.0)
  Eigen::VectorXd correlationEnergyTurbo_eV(10);
  correlationEnergyTurbo_eV << 17.789, 8.863, 0.625, 1.410, 1.720, -0.486, -0.575, -1.474, -1.976, -1.947;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergyTurbo_eV.size(); i++) {
      if (i == 1)
        continue;
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(i), 1.5e-3);
    }
  };
}
TEST_F(GWTaskTest, GW_CD_CD_R) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.densFitJ = Options::DENS_FITS::ACD;
  settings.basis.densFitK = Options::DENS_FITS::ACD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACD;
  settings.basis.cdThreshold = 1e-12;
  settings.basis.extendSphericalACDShells = Options::EXTEND_ACD::COMPLETE;
  settings.basis.label = "DEF2-SVP";
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Serenity without RI
  EXPECT_NEAR(-76.2723726160, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 5;
  gw.settings.nOcc = 5;
  // set to zero because closer to the turbomole results
  gw.settings.eta = 0.0;
  gw.settings.densFitCache = Options::DENS_FITS::ACD;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  // Turbomole 7.4.1 gw/ri approximation/analytic/eta 0.001(but seems to be done with eta =0.0)
  // The first value has been changed with a value obtained without RI in the full Analytic implementation in Serenity
  Eigen::VectorXd correlationEnergyTurbo_eV(10);
  correlationEnergyTurbo_eV << 17.78012, 8.863, 0.625, 1.410, 1.720, -0.486, -0.575, -1.474, -1.976, -1.947;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergyTurbo_eV.size(); i++) {
      if (i == 1)
        continue;
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(i), 3.3e-3);
    }
  };
  std::remove((waterMon->getSettings().path + "ACD-DEF2-SVP").c_str());
}

TEST_F(GWTaskTest, GW_CD_U) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = system->getSettings();
  settings.grid.accuracy = 5;
  settings.grid.smallGridAccuracy = 3;
  system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  // Perform SCF
  auto scf = ScfTask<UNRESTRICTED>(system);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-39.73968410168, system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<UNRESTRICTED>({system});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 2;
  gw.settings.nOcc = 2;
  gw.settings.eta = 0.001;
  gw.run();
  // Turbomole 7.4.1
  Eigen::VectorXd correlationEnergyTurbo_eV(8);
  correlationEnergyTurbo_eV << 0.48028999, 0.49317720, 0.89494920, -0.56708293, 0.31385060, 0.32801042, -2.29475965, -0.76964986;
  auto correlationEnergy = gw._correlationEnergies;
  unsigned counter = 0;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergy_spin.size(); i++) {
      if (std::abs(correlationEnergy_spin(i)) > 1e-8) {
        EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(counter), 1.2e-3);
        counter++;
      }
    }
  };
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(GWTaskTest, GW_CD_R_NAF) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-76.27245320915, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 2;
  gw.settings.nOcc = 2;
  // set to zero because closer to the turbomole results
  gw.settings.eta = 0.0;
  gw.run();
  auto correlationEnergy1 = gw._correlationEnergies;
  auto gw2 = GWTask<RESTRICTED>({waterMon});
  gw2.settings.gwtype = Options::GWALGORITHM::CD;
  gw2.settings.nVirt = 2;
  gw2.settings.nOcc = 2;
  gw2.settings.nafThresh = 1e-4;
  // set to zero because closer to the turbomole results
  gw2.settings.eta = 0.0;
  gw2.run();
  auto correlationEnergy2 = gw2._correlationEnergies;
  for_spin(correlationEnergy1, correlationEnergy2) {
    for (unsigned int i = 0; i < correlationEnergy1_spin.size(); i++) {
      EXPECT_NEAR(correlationEnergy1_spin(i), correlationEnergy2_spin(i), 3e-3);
    }
  };
}

TEST_F(GWTaskTest, GW_AC_R) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-76.27245320915, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.gwtype = Options::GWALGORITHM::AC;
  gw.settings.nVirt = 1;
  gw.settings.nOcc = 1;
  // set to zero because closer to the turbomole results
  gw.settings.eta = 0.0;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  // Turbomole 7.4.1 gw/ri approximation/analytic/eta 0.001(but seems to be done with eta =0.0)
  Eigen::VectorXd correlationEnergyTurbo_eV(7);
  correlationEnergyTurbo_eV << 0.0, 0.0, 0.0, 0.0, 1.720, -0.486, 0.0;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergyTurbo_eV.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(i), 1.5e-3);
    }
  };
}
TEST_F(GWTaskTest, GW_AC_U) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = system->getSettings();
  settings.grid.accuracy = 5;
  settings.grid.smallGridAccuracy = 3;
  system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  // Perform SCF
  auto scf = ScfTask<UNRESTRICTED>(system);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-39.73968410168, system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<UNRESTRICTED>({system});
  gw.settings.gwtype = Options::GWALGORITHM::AC;
  gw.settings.nVirt = 2;
  gw.settings.nOcc = 2;
  gw.settings.padePoints = 128;
  gw.settings.eta = 0.001;
  gw.settings.linearized = true;
  gw.run();
  // Turbomole 7.4.1
  Eigen::VectorXd correlationEnergyTurbo_eV(8);
  correlationEnergyTurbo_eV << 0.48028999, 0.49317720, 0.89494920, -0.56708293, 0.31385060, 0.32801042, -2.29475965, -0.76964986;
  auto correlationEnergy = gw._correlationEnergies;
  unsigned counter = 0;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergy_spin.size(); i++) {
      if (std::abs(correlationEnergy_spin(i)) > 1e-8) {
        EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(counter), 22e-3);
        counter++;
      }
    }
  };
  SystemController__TEST_SUPPLY::cleanUp();
}
TEST_F(GWTaskTest, embGW_small) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, true);

  Settings settings1 = act->getSettings();
  settings1.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings1.basis.label = "DEF2-TZVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, settings1);

  Settings settings2 = env->getSettings();
  settings2.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings2.basis.label = "DEF2-TZVP";
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, settings2);

  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.run();

  auto gw = GWTask<Options::SCF_MODES::RESTRICTED>({act}, {env});
  gw.settings.eta = 0.001;
  gw.settings.diis = true;
  gw.settings.nVirt = 1;
  gw.settings.nOcc = 1;
  gw.settings.qpiterations = 10;
  gw.settings.ConvergenceThreshold = 1e-6;
  gw.settings.integrationPoints = 128;
  gw.settings.subsystemAuxillaryBasisOnly = true;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  Eigen::VectorXd correlation2(6);
  // Serenity 19.05.2022
  correlation2 << 0.0, 0.0, 0.0, 0.0, 1.819371555985553, -0.69471338695164397;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlation2.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlation2(i), 2e-6);
    }
  };
}

TEST_F(GWTaskTest, embGW_small_AC) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, true);

  Settings settings1 = act->getSettings();
  settings1.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings1.basis.label = "DEF2-TZVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, settings1);

  Settings settings2 = env->getSettings();
  settings2.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings2.basis.label = "DEF2-TZVP";
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, settings2);

  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.run();

  auto gw = GWTask<Options::SCF_MODES::RESTRICTED>({act}, {env});
  gw.settings.eta = 0.001;
  gw.settings.diis = true;
  gw.settings.gwtype = Options::GWALGORITHM::AC;
  gw.settings.nVirt = 1;
  gw.settings.nOcc = 1;
  gw.settings.qpiterations = 10;
  gw.settings.ConvergenceThreshold = 1e-6;
  gw.settings.integrationPoints = 128;
  gw.settings.padePoints = 128;
  gw.settings.subsystemAuxillaryBasisOnly = true;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  Eigen::VectorXd correlation2(6);
  // Serenity 19.05.2022 (from CD)
  correlation2 << 0.0, 0.0, 0.0, 0.0, 1.819371555985553, -0.69471338695164397;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlation2.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlation2(i), 1e-5);
    }
  };
}

TEST_F(GWTaskTest, embGW_large) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, true);

  Settings settings1 = act->getSettings();
  settings1.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings1.basis.label = "DEF2-SVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, settings1);

  Settings settings2 = env->getSettings();
  settings2.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings2.basis.label = "DEF2-SVP";
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, settings2);

  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.run();

  auto gw = GWTask<Options::SCF_MODES::RESTRICTED>({act}, {env});
  gw.settings.eta = 0.001;
  gw.settings.evGW = true;
  gw.settings.diis = true;
  gw.settings.nVirt = 1;
  gw.settings.nOcc = 1;
  gw.settings.evGWcycles = 10;
  gw.settings.ConvergenceThreshold = 1e-6;
  gw.settings.integrationPoints = 128;
  gw.settings.gap = true;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  Eigen::VectorXd correlation2(6);
  // Serenity 19.05.2022
  correlation2 << 0.0, 0.0, 0.0, 0.0, 1.6675961049528527, -0.57892512045600975;

  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlation2.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlation2(i), 3e-6);
    }
  };
}

#ifdef SERENITY_USE_LAPLACE_MINIMAX

TEST_F(GWTaskTest, GW_CD_LT_R) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-76.27245320915, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 5;
  gw.settings.nOcc = 5;
  // set to zero because closer to the turbomole results
  gw.settings.eta = 0.0;
  gw.settings.ltconv = 1e-7;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  // Turbomole 7.4.1 gw/ri approximation/analytic/eta 0.001(but seems to be done with eta =0.0)
  Eigen::VectorXd correlationEnergyTurbo_eV(10);
  correlationEnergyTurbo_eV << 17.789, 8.863, 0.625, 1.410, 1.720, -0.486, -0.575, -1.474, -1.976, -1.947;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergyTurbo_eV.size(); i++) {
      if (i == 1)
        continue;
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(i), 1.5e-3);
    }
  };
}

TEST_F(GWTaskTest, embGW_small_LT) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, true);

  Settings settings1 = act->getSettings();
  settings1.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings1.basis.label = "DEF2-TZVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, settings1);

  Settings settings2 = env->getSettings();
  settings2.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings2.basis.label = "DEF2-TZVP";
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, settings2);

  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.run();

  auto gw = GWTask<Options::SCF_MODES::RESTRICTED>({act}, {env});
  gw.settings.eta = 0.001;
  gw.settings.diis = true;
  gw.settings.nVirt = 1;
  gw.settings.nOcc = 1;
  gw.settings.qpiterations = 10;
  gw.settings.ConvergenceThreshold = 1e-6;
  gw.settings.integrationPoints = 128;
  gw.settings.subsystemAuxillaryBasisOnly = true;
  gw.settings.ltconv = 1e-7;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  Eigen::VectorXd correlation2(6);
  // Serenity 19.05.2022 (from CD)
  correlation2 << 0.0, 0.0, 0.0, 0.0, 1.819371555985553, -0.69471338695164397;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlation2.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlation2(i), 4e-5);
    }
  };
}

TEST_F(GWTaskTest, GW_CD_R_NAF_LT) {
  auto waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  Settings settings = waterMon->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  waterMon = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(waterMon);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-76.27245320915, waterMon->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<RESTRICTED>({waterMon});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 2;
  gw.settings.nOcc = 2;
  // set to zero because closer to the turbomole results
  gw.settings.eta = 0.0;
  gw.run();
  auto correlationEnergy1 = gw._correlationEnergies;
  auto gw2 = GWTask<RESTRICTED>({waterMon});
  gw2.settings.gwtype = Options::GWALGORITHM::CD;
  gw2.settings.nVirt = 2;
  gw2.settings.nOcc = 2;
  gw2.settings.nafThresh = 1e-4;
  gw2.settings.ltconv = 1e-7;
  // set to zero because closer to the turbomole results
  gw2.settings.eta = 0.0;
  gw2.run();
  auto correlationEnergy2 = gw2._correlationEnergies;
  for_spin(correlationEnergy1, correlationEnergy2) {
    for (unsigned int i = 0; i < correlationEnergy1_spin.size(); i++) {
      EXPECT_NEAR(correlationEnergy1_spin(i), correlationEnergy2_spin(i), 3e-3);
    }
  };
}

TEST_F(GWTaskTest, GW_CD_LT_U) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = system->getSettings();
  settings.grid.accuracy = 5;
  settings.grid.smallGridAccuracy = 3;
  system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  // Perform SCF
  auto scf = ScfTask<UNRESTRICTED>(system);
  scf.run();
  // Turbomole 7.4.1 ridft/PBE/m3
  EXPECT_NEAR(-39.73968410168, system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-5);
  // Perform CD GW calculation
  auto gw = GWTask<UNRESTRICTED>({system});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 2;
  gw.settings.nOcc = 2;
  gw.settings.eta = 0.001;
  gw.settings.ltconv = 1e-7;
  gw.run();
  // Turbomole 7.4.1
  Eigen::VectorXd correlationEnergyTurbo_eV(8);
  correlationEnergyTurbo_eV << 0.48028999, 0.49317720, 0.89494920, -0.56708293, 0.31385060, 0.32801042, -2.29475965, -0.76964986;
  auto correlationEnergy = gw._correlationEnergies;
  unsigned counter = 0;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergy_spin.size(); i++) {
      if (std::abs(correlationEnergy_spin(i)) > 1e-8) {
        EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(counter), 1.2e-3);
        counter++;
      }
    }
  };
  SystemController__TEST_SUPPLY::cleanUp();
}

#endif /* SERENITY_USE_LAPLACE_MINIMAX */

} // namespace Serenity
