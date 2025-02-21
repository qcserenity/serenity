/**
 * @file   SystemController__TEST_SUPPLY.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   28. August 2015, 14:19
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
/* Include Class Header*/
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "geometry/AtomTypeFactory.h"
#include "geometry/Geometry.h"
#include "io/Filesystem.h"
#include "parameters/Constants.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

std::map<TEST_SYSTEM_CONTROLLERS, std::shared_ptr<SystemController>> SystemController__TEST_SUPPLY::_testSystemControllers;
std::map<TEST_SYSTEM_CONTROLLERS, std::shared_ptr<Geometry>> SystemController__TEST_SUPPLY::_testGeometries;

void SystemController__TEST_SUPPLY::prepare(TEST_SYSTEM_CONTROLLERS kind, bool fromScratch) {
  std::string pathToTestsResources;
  std::string basisPath;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/";
    basisPath = (std::string)env_p + "basis/";
  }
  else {
    throw SerenityError("ERROR Environment variable SERENITY_RESOURCES not set.");
  }

  switch (kind) {
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_H2_MINBAS";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_H2_MINBAS_ACTIVE";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE_LRSCF: {
      Settings settings;
      settings.basis.label = "STO-3G";
      settings.name = "TestSystem_H2_MINBAS_ACTIVE_LRSCF";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.basis.incrementalSteps = 1;
      settings.basis.integralThreshold = 1e-64;
      settings.scf.energyThreshold = 1e-10;
      settings.scf.rmsdThreshold = 1e-10;
      settings.scf.diisThreshold = 1e-9;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_H2_MINBAS_ENVIRONMENT";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT_LRSCF: {
      Settings settings;
      settings.basis.label = "STO-3G";
      settings.name = "TestSystem_H2_MINBAS_ENVIRONMENT_LRSCF";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.basis.incrementalSteps = 1;
      settings.basis.integralThreshold = 1e-64;
      settings.scf.energyThreshold = 1e-10;
      settings.scf.rmsdThreshold = 1e-10;
      settings.scf.diisThreshold = 1e-9;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LARGE_DISTANCE: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_H2_MINBAS_LARGE_DISTANCE";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::CO_MINBAS: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_CO_MINBAS";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_WATER_DISTORTED_MINBAS";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2_DEF2_TZVP";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2_DEF2_TZVP_PBE";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2_DEF2_TZVP_PW91_UNRES_0";
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.spin = 0;
      settings.charge = 0;
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PW91;
      settings.grid.accuracy = 4;
      settings.grid.smallGridAccuracy = 2;
      settings.scf.initialguess = Options::INITIAL_GUESSES::ATOM_SCF;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_BE_DEF2_TZVP_PW91_UNRES_0";
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.spin = 0;
      settings.charge = 0;
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PW91;
      settings.grid.accuracy = 4;
      settings.grid.smallGridAccuracy = 2;
      settings.scf.initialguess = Options::INITIAL_GUESSES::ATOM_SCF;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2_DEF2_TZVP_PW91_UNRES_1";
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.spin = 1;
      settings.charge = 1;
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PW91;
      settings.grid.accuracy = 4;
      settings.grid.smallGridAccuracy = 2;
      settings.scf.initialguess = Options::INITIAL_GUESSES::ATOM_SCF;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_BE_DEF2_TZVP_PW91_UNRES_1";
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.spin = 1;
      settings.charge = 1;
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PW91;
      settings.grid.accuracy = 4;
      settings.grid.smallGridAccuracy = 2;
      settings.scf.initialguess = Options::INITIAL_GUESSES::ATOM_SCF;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE_NORI: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2_DEF2_TZVP_PBE_NORI";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE_NORI), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_CAMB3LYP: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2_DEF2_TZVP_CAMB3LYP";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      settings.basis.incrementalSteps = 0;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_CAMB3LYP), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2_DEF2_TZVP_HF";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED: {
      Settings settings;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2_DEF2_TZVP_HF_UNRESTRICTED";
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_SV_P_PBE: {
      Settings settings;
      settings.basis.label = "DEF2-SV_P_";
      settings.name = "TestSystem_H2_DEF2_SV_P_PBE";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      _testSystemControllers[kind] =
          std::make_shared<SystemController>(getGeometry(TEST_SYSTEM_CONTROLLERS::H2_DEF2_SV_P_PBE), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::JACOBSEN_MINBAS: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_JACOBSEN_MINBAS";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::C60_MINBAS: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_C60_MINBAS";
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H60_Ghost_MINBAS: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_H60_Ghost_MINBAS";
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_WaterMonOne_6_31Gs";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_WaterMonOne_Def2_SVP";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_WaterMonOne_Def2_SVP_B2PLYP";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_WaterMonTwo_6_31Gs";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_WaterMonTwo_Def2_SVP";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_WaterMonTwo_6_31Gs_DFT";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.basis.basisLibPath = basisPath;
      settings.grid.smallGridAccuracy = 4;
      settings.grid.accuracy = 4;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMon1_3_DFT: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_WaterMon1_3_DFT";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.basis.basisLibPath = basisPath;
      settings.grid.smallGridAccuracy = 5;
      settings.grid.accuracy = 5;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMon2_3_DFT: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_WaterMon2_3_DFT";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.basis.basisLibPath = basisPath;
      settings.grid.smallGridAccuracy = 5;
      settings.grid.accuracy = 5;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMon3_3_DFT: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_WaterMon3_3_DFT";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.basis.basisLibPath = basisPath;
      settings.grid.smallGridAccuracy = 5;
      settings.grid.accuracy = 5;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterDim_DFT: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_WaterDim_DFT";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.basis.basisLibPath = basisPath;
      settings.grid.smallGridAccuracy = 5;
      settings.grid.accuracy = 5;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_WaterMonOne_6_31Gs_DFT";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.basis.basisLibPath = basisPath;
      settings.grid.smallGridAccuracy = 4;
      settings.grid.accuracy = 4;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Ne2_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_Ne2_6_31Gs";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Ar2_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_Ar2_6_31Gs";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Kr2_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_Kr2_6_31Gs";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::F2_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_F2_6_31Gs";
      settings.basis.basisLibPath = basisPath;
      settings.basis.incrementalSteps = 1;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Cl2_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_Cl2_6_31Gs";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Br2_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_Br2_6_31Gs";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_BP86";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::LDA;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ACTIVE_FDE";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86_Supermolecular: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86_Supermolecular";
      settings.basis.basisLibPath = basisPath;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::LDA;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_H2_def2_SVP_ACTIVE_FDE";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86";
      settings.basis.basisLibPath = basisPath;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      settings.basis.incrementalSteps = 1;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_H2_def2_SVP_ACTIVE_FDE_BP86";
      settings.basis.basisLibPath = basisPath;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_H2_def2_SVP_ENVIRONMENT_FDE_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_PBE0_UNRESTRICTED: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "DEF2-SVP";
      settings.basis.incrementalSteps = 1;
      settings.name = "TestSystem_H2_def2_SVP_ACTIVE_FDE_PBE0_UNRESTRICTED";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_PBE0_UNRESTRICTED: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "DEF2-SVP";
      settings.basis.incrementalSteps = 1;
      settings.name = "TestSystem_H2_def2_SVP_ENVIRONMENT_FDE_PBE0_UNRESTRICTED";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_B3LYP: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP";
      settings.basis.basisLibPath = basisPath;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.basis.integralIncrementThresholdStart = 1e-12;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::LDA;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      settings.basis.incrementalSteps = 1;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_B3LYP: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP";
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.basis.integralIncrementThresholdStart = 1e-12;
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He2_6_31Gs_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_He2_6_31Gs_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_He_1_6_31Gs_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_2_6_31Gs_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_He_2_6_31Gs_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_3_6_31Gs_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_He_3_6_31Gs_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_He2_def2SVP_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::HeNe_def2SVP_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_HeNe_def2SVP_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_1_def2SVP_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_He_1_def2SVP_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_2_def2SVP_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_He_2_def2SVP_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_3_def2SVP_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_He_3_def2SVP_BP86";
      settings.grid.accuracy = 7;
      settings.grid.smallGridAccuracy = 7;
      settings.grid.blockAveThreshold = 0;
      settings.grid.basFuncRadialThreshold = 0;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_MINBAS_SING: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_O2_MINBAS_SING";
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_O2_MINBAS_TRIP";
      settings.spin = 2;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP_CIS: {
      Settings settings;
      settings.name = "TestSystem_O2_MINBAS_TRIP_CIS";
      settings.load = pathToTestsResources;
      _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      if (fromScratch)
        throw SerenityError("This TestSystemController is only provided from disk: " + settings.name);
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDA: {
      Settings settings;
      settings.name = "TestSystem_O2_TRIP_DEF2_SVP_TDA";
      settings.load = pathToTestsResources;
      _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      if (fromScratch)
        throw SerenityError("This TestSystemController is only provided from disk: " + settings.name);
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDHF: {
      Settings settings;
      settings.name = "TestSystem_O2_TRIP_DEF2_SVP_TDHF";
      settings.load = pathToTestsResources;
      _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      if (fromScratch)
        throw SerenityError("This TestSystemController is only provided from disk: " + settings.name);
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDDFT: {
      Settings settings;
      settings.name = "TestSystem_O2_TRIP_DEF2_SVP_TDDFT";
      settings.load = pathToTestsResources;
      _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      if (fromScratch)
        throw SerenityError("This TestSystemController is only provided from disk: " + settings.name);
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_6_31G_PBE_TDA: {
      Settings settings;
      settings.name = "TestSystem_O2_TRIP_6_31G_PBE_TDA";
      settings.load = pathToTestsResources;
      _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      if (fromScratch)
        throw SerenityError("This TestSystemController is only provided from disk: " + settings.name);
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_6_31G_PBE_TDDFT: {
      Settings settings;
      settings.name = "TestSystem_O2_TRIP_6_31G_PBE_TDDFT";
      settings.load = pathToTestsResources;
      _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      if (fromScratch)
        throw SerenityError("This TestSystemController is only provided from disk: " + settings.name);
    } break;
    case TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs: {
      Settings settings;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_F_MINUS_6_31Gs";
      settings.charge = -1;
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_WATER_DEF2_SVP";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      settings.basis.basisLibPath = basisPath;
      settings.basis.incrementalSteps = 0;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::OH_MINBAS_PBE: {
      Settings settings;
      settings.basis.label = "STO-3G";
      settings.name = "TestSystem_OH_MINBAS_PBE";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.basis.basisLibPath = basisPath;
      settings.spin = 1;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE: {
      Settings settings;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_MethylRad_def2_SVP_PBE";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.basis.basisLibPath = basisPath;
      settings.spin = 1;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_HI_Def2_SVP_PBE";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.basis.basisLibPath = basisPath;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.makeSphericalBasis = true;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_I2_Def2_SVP_PBE";
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.basis.basisLibPath = basisPath;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.makeSphericalBasis = true;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LINEARDEPENDENT: {
      Settings settings;
      settings.basis.label = "STO-6G";
      settings.name = "TestSystem_H2_MINBAS_LINEARDEPENDENT";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_EthaneA_Def2_SVP_BP86_Act";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_EthaneA_Def2_SVP";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Env: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_EthaneA_Def2_SVP_BP86_Env";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_EthaneA_Def2_SVP_BP86";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.grid.smallGridAccuracy = 4;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_EthaneB_Def2_SVP_BP86";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.grid.smallGridAccuracy = 4;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ: {
      Settings settings;
      settings.basis.label = "AUG-CC-PVDZ";
      settings.name = "TestSystem_Formaldehyde_HF_AUG_CC_PVDZ";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ: {
      Settings settings;
      settings.basis.label = "AUG-CC-PVDZ";
      settings.name = "TestSystem_Diaziridine_HF_AUG_CC_PVDZ";
      settings.basis.basisLibPath = basisPath;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_HYBRID: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ACTIVE_FDE";
      settings.basis.basisLibPath = basisPath;
      settings.grid.accuracy = 2;
      settings.grid.smallGridAccuracy = 2;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_HYBRID: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
      settings.basis.label = "6-31GS";
      settings.name = "TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
      settings.grid.accuracy = 2;
      settings.grid.smallGridAccuracy = 2;
      settings.basis.densFitJ = Options::DENS_FITS::NONE;
      settings.basis.densFitK = Options::DENS_FITS::NONE;
      settings.basis.densFitLRK = Options::DENS_FITS::NONE;
      settings.basis.densFitCorr = Options::DENS_FITS::NONE;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2Dimer_Def2_TZVP: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_H2Dimer_Def2_TZVP";
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_Hydroxy_Def2_SVP_Act";
      settings.charge = -1;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::PHENOLATE_O_DEF2_SVP_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_Phenolate_O_Def2_SVP_BP86";
      settings.grid.accuracy = 4;
      settings.grid.smallGridAccuracy = 4;
      settings.charge = -2;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_Phenyl_Def2_SVP_Env";
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.charge = +1;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::PHENOLATE_PHENYL_DEF2_SVP_BP86: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::BP86;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_Phenolate_Phenyl_Def2_SVP_BP86";
      settings.grid.accuracy = 4;
      settings.grid.smallGridAccuracy = 4;
      settings.charge = +1;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "6-311G";
      settings.name = "TestSystem_H2_6_311G_3A";
      settings.spin = 2;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "6-311G";
      settings.name = "TestSystem_H2_6_311G_Act_3A";
      settings.spin = 1;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "6-311G";
      settings.name = "TestSystem_H2_6_311G_Env_3A";
      settings.spin = 1;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_BS_3A: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "6-311G";
      settings.name = "TestSystem_H2_6_311G_BS_3A";
      settings.spin = 0;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_BSAct_3A: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "6-311G";
      settings.name = "TestSystem_H2_6_311G_BSAct_3A";
      settings.spin = -1;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_BSEnv_3A: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "6-311G";
      settings.name = "TestSystem_H2_6_311G_BSEnv_3A";
      settings.spin = 1;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::MethylRad_Act_def2_SVP_PBE: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_MethylRad_Act_def2_SVP_PBE";
      settings.spin = 1;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::MethylRad_Env_def2_SVP_PBE: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_MethylRad_Env_def2_SVP_PBE";
      settings.spin = 1;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF: {
      // Note that a turbomole HF orbital file is available in the testresources!
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_WCCR1010_P1_Def2_SVP_HF";
      settings.charge = +1;
      settings.spin = 0;
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF: {
      Settings settings;
      settings.basis.label = "DEF2-SVP";
      settings.name = "TestSystem_Water_Dimer_def2_SVP_HF";
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Def2_TZVP_DFT: {
      // Note that a turbomole HF orbital file is available in the testresources!
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "DEF2-TZVP";
      settings.name = "TestSystem_Water_Def2-TZVP_DFT";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_Def2_TZVP_DFT: {
      // Note that a turbomole HF orbital file is available in the testresources!
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "def2-TZVP";
      settings.name = "TestSystem_He_Def2-TZVP_DFT";
      settings.basis.basisLibPath = basisPath;
      if (fromScratch) {
        _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
      }
      else {
        settings.load = pathToTestsResources;
        _testSystemControllers[kind] = std::make_shared<SystemController>(settings);
      }
    } break;
    case TEST_SYSTEM_CONTROLLERS::ETHANOL_def2_SVP_HF: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "def2-SVP";
      settings.name = "TestSystem_Ethanol_Def2-SVP_HF";
      settings.grid.smallGridAccuracy = 3;
      settings.grid.accuracy = 5;
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_A: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "DEF2-SVP";
      settings.grid.smallGridAccuracy = 4;
      settings.grid.accuracy = 4;
      settings.basis.integralIncrementThresholdStart = 1e-10;
      settings.name = "TestSystem_Water_Hexamer_Monomer_A";
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_B: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "DEF2-SVP";
      settings.grid.smallGridAccuracy = 4;
      settings.grid.accuracy = 4;
      settings.basis.integralIncrementThresholdStart = 1e-10;
      settings.name = "TestSystem_Water_Hexamer_Monomer_B";
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_C: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "DEF2-SVP";
      settings.grid.smallGridAccuracy = 4;
      settings.grid.accuracy = 4;
      settings.basis.integralIncrementThresholdStart = 1e-10;
      settings.name = "TestSystem_Water_Hexamer_Monomer_C";
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_D: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "DEF2-SVP";
      settings.grid.smallGridAccuracy = 4;
      settings.grid.accuracy = 4;
      settings.basis.integralIncrementThresholdStart = 1e-10;
      settings.name = "TestSystem_Water_Hexamer_Monomer_D";
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_E: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "DEF2-SVP";
      settings.grid.smallGridAccuracy = 4;
      settings.grid.accuracy = 4;
      settings.basis.integralIncrementThresholdStart = 1e-10;
      settings.name = "TestSystem_Water_Hexamer_Monomer_E";
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_F: {
      Settings settings;
      settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
      settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
      settings.scfMode = Options::SCF_MODES::RESTRICTED;
      settings.basis.label = "DEF2-SVP";
      settings.grid.smallGridAccuracy = 4;
      settings.grid.accuracy = 4;
      settings.basis.integralIncrementThresholdStart = 1e-10;
      settings.name = "TestSystem_Water_Hexamer_Monomer_F";
      settings.basis.basisLibPath = basisPath;
      _testSystemControllers[kind] = std::make_shared<SystemController>(getGeometry(kind), settings);
    } break;
  } // switch
}

std::shared_ptr<SystemController>
SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS kind, Settings& settings, int charge, int spin) {
  if (settings.name == "default")
    settings.name = "SomeTestSystem";
  settings.charge = charge;
  settings.spin = spin;
  return std::make_shared<SystemController>(getGeometry(kind), settings);
}

void SystemController__TEST_SUPPLY::prepareGeometry(TEST_SYSTEM_CONTROLLERS kind) {
  switch (kind) {
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE:
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE_LRSCF:
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS:
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP:
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE:
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE_NORI:
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_CAMB3LYP:
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF:
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED:
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_SV_P_PBE:
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86:
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE:
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE:
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_B3LYP:
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86:
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_BP86:
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_PBE0_UNRESTRICTED: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, -0.7);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.7);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE_LRSCF] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_MINBAS] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE_NORI] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_CAMB3LYP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_SV_P_PBE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_B3LYP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_PBE0_UNRESTRICTED] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT:
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT_LRSCF: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 1.7);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 2.4);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT_LRSCF] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LARGE_DISTANCE: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, -150.0);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 150.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LARGE_DISTANCE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::CO_MINBAS: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.0, 0.0, 0.0);
      auto O2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 1.128058 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::CO_MINBAS] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, O2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS: {
      auto O = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 0.0);
      auto H1 =
          std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.2 * ANGSTROM_TO_BOHR, -0.9 * ANGSTROM_TO_BOHR, 0.0);
      auto H2 =
          std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.9 * ANGSTROM_TO_BOHR, 0.75 * ANGSTROM_TO_BOHR, 0.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O, H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.384 * ANGSTROM_TO_BOHR, 0.0,
                                       1.184 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.384 * ANGSTROM_TO_BOHR, 0.0,
                                       1.184 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0: {
      auto Be = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Be"), 0.0, 0.0, -0.592 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{Be});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.384 * ANGSTROM_TO_BOHR, 0.0,
                                       1.184 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.384 * ANGSTROM_TO_BOHR, 0.0,
                                       1.184 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1: {
      auto Be = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Be"), 0.0, 0.0, -0.592 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{Be});
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs:
    case TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT:
    case TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP:
    case TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP:
    case TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.971647 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.938289 * ANGSTROM_TO_BOHR, 0.0,
                                       -0.252413 * ANGSTROM_TO_BOHR);
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 0.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs:
    case TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT:
    case TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.195616 * ANGSTROM_TO_BOHR,
                                       1.506265 * ANGSTROM_TO_BOHR, -0.578922 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.476701 * ANGSTROM_TO_BOHR,
                                       1.813197 * ANGSTROM_TO_BOHR, -1.380997 * ANGSTROM_TO_BOHR);
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), -1.774058 * ANGSTROM_TO_BOHR,
                                       2.235000 * ANGSTROM_TO_BOHR, -0.859006 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Ne2_6_31Gs: {
      auto Ne1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Ne"), 0.0, 0.0, 0.0);
      auto Ne2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Ne"), 0.0, 0.0, 2.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Ne2_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{Ne1, Ne2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Ar2_6_31Gs: {
      auto Ar1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Ar"), 0.0, 0.0, 0.0);
      auto Ar2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Ar"), 0.0, 0.0, 2.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Ar2_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{Ar1, Ar2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Kr2_6_31Gs: {
      auto Kr1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Kr"), 0.0, 0.0, 0.0);
      auto Kr2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Kr"), 0.0, 0.0, 2.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Kr2_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{Kr1, Kr2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::F2_6_31Gs: {
      auto F1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("F"), 0.0, 0.0, 0.0);
      auto F2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("F"), 0.0, 0.0, 2.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::F2_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{F1, F2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Cl2_6_31Gs: {
      auto Cl1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Cl"), 0.0, 0.0, 0.0);
      auto Cl2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Cl"), 0.0, 0.0, 2.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Cl2_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{Cl1, Cl2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Br2_6_31Gs: {
      auto Br1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Br"), 0.0, 0.0, 0.0);
      auto Br2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Br"), 0.0, 0.0, 2.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Br2_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{Br1, Br2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::JACOBSEN_MINBAS: {
      std::vector<std::shared_ptr<Atom>> atoms;
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.04954507139919 * ANGSTROM_TO_BOHR,
                                             0.49346459974284 * ANGSTROM_TO_BOHR, 0.19403708095030 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.68994584664293 * ANGSTROM_TO_BOHR,
                                             0.24498551943875 * ANGSTROM_TO_BOHR, -0.48728935738284 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.91938928475996 * ANGSTROM_TO_BOHR,
                                             0.26021745455615 * ANGSTROM_TO_BOHR, 1.25855031165228 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.14426501541109 * ANGSTROM_TO_BOHR,
                                             -0.38993518976336 * ANGSTROM_TO_BOHR, -0.40317206956060 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("N"), -2.28685698583691 * ANGSTROM_TO_BOHR,
                                             1.94767379863413 * ANGSTROM_TO_BOHR, 0.14878941970778 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.27951756876751 * ANGSTROM_TO_BOHR,
                                             -1.22058665005177 * ANGSTROM_TO_BOHR, -0.36550874084111 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("N"), 0.24097504669724 * ANGSTROM_TO_BOHR,
                                             1.23597051240999 * ANGSTROM_TO_BOHR, 0.07433332522267 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.79817613961669 * ANGSTROM_TO_BOHR,
                                             0.50662260467413 * ANGSTROM_TO_BOHR, -1.54994979799048 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.73139840250440 * ANGSTROM_TO_BOHR,
                                             -1.86646507169582 * ANGSTROM_TO_BOHR, -0.30459890356185 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.31426754402753 * ANGSTROM_TO_BOHR,
                                             -0.11157394102816 * ANGSTROM_TO_BOHR, -1.45104048488576 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.08740996064897 * ANGSTROM_TO_BOHR,
                                             -0.23965147661296 * ANGSTROM_TO_BOHR, 0.13329358281900 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.37237236245536 * ANGSTROM_TO_BOHR,
                                             -2.11751151850859 * ANGSTROM_TO_BOHR, -0.96740655951936 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.50036441159888 * ANGSTROM_TO_BOHR,
                                             -2.49528484810240 * ANGSTROM_TO_BOHR, -0.76674625385721 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.67692478831873 * ANGSTROM_TO_BOHR,
                                             -2.15416547391794 * ANGSTROM_TO_BOHR, 0.75385480191594 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.44921617537744 * ANGSTROM_TO_BOHR,
                                             -1.91734894037221 * ANGSTROM_TO_BOHR, -2.04497106436333 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.08572552070489 * ANGSTROM_TO_BOHR,
                                             -3.16981456558983 * ANGSTROM_TO_BOHR, -0.85948673907074 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.13020491148699 * ANGSTROM_TO_BOHR,
                                             -1.46885695236768 * ANGSTROM_TO_BOHR, 0.69354276528090 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.66954724487480 * ANGSTROM_TO_BOHR,
                                             -1.39561995215723 * ANGSTROM_TO_BOHR, -0.88425085333822 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 1.46008944970184 * ANGSTROM_TO_BOHR,
                                             0.98514467699238 * ANGSTROM_TO_BOHR, 0.39769354812489 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.45450814424032 * ANGSTROM_TO_BOHR,
                                             2.47238300834947 * ANGSTROM_TO_BOHR, 0.00318683222481 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("Mn"), -0.62740577664457 * ANGSTROM_TO_BOHR,
                                             3.00042723327634 * ANGSTROM_TO_BOHR, 0.24981138311896 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("Cl"), -0.77190534903661 * ANGSTROM_TO_BOHR,
                                             2.27052105147490 * ANGSTROM_TO_BOHR, 2.55668552389363 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.94717409770232 * ANGSTROM_TO_BOHR,
                                             3.87889063769865 * ANGSTROM_TO_BOHR, 0.76949369546203 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), -1.52528344714375 * ANGSTROM_TO_BOHR,
                                             4.56769509314170 * ANGSTROM_TO_BOHR, 0.75501905102585 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.80783982424845 * ANGSTROM_TO_BOHR,
                                             4.82121758023799 * ANGSTROM_TO_BOHR, 0.62748585937203 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.75986953640787 * ANGSTROM_TO_BOHR,
                                             3.86533208523123 * ANGSTROM_TO_BOHR, 0.15328097599859 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.27985309595695 * ANGSTROM_TO_BOHR,
                                             6.13034947382499 * ANGSTROM_TO_BOHR, 0.98003712378897 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.30569253403915 * ANGSTROM_TO_BOHR,
                                             1.81943364053180 * ANGSTROM_TO_BOHR, -0.19958098334289 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.09866363284808 * ANGSTROM_TO_BOHR,
                                             4.24099203805007 * ANGSTROM_TO_BOHR, -0.07331640048439 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.40005917295061 * ANGSTROM_TO_BOHR,
                                             1.94436165989897 * ANGSTROM_TO_BOHR, 0.90821046986723 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.83433169313173 * ANGSTROM_TO_BOHR,
                                             -0.03637112556795 * ANGSTROM_TO_BOHR, 0.29353084627507 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.10873332863783 * ANGSTROM_TO_BOHR,
                                             3.33656524353963 * ANGSTROM_TO_BOHR, 1.06263835297119 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 3.14926484875961 * ANGSTROM_TO_BOHR,
                                             4.19561876428366 * ANGSTROM_TO_BOHR, 1.54849799956785 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 3.67122070870860 * ANGSTROM_TO_BOHR,
                                             1.43860539002233 * ANGSTROM_TO_BOHR, 1.25009588091058 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.60561872832707 * ANGSTROM_TO_BOHR,
                                             6.43509333741244 * ANGSTROM_TO_BOHR, 0.70903643281861 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.54128397814805 * ANGSTROM_TO_BOHR,
                                             5.52849698547268 * ANGSTROM_TO_BOHR, 0.16428488535104 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -5.77510362299943 * ANGSTROM_TO_BOHR,
                                             3.47954238775536 * ANGSTROM_TO_BOHR, -0.44905000250600 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.35121862890659 * ANGSTROM_TO_BOHR,
                                             7.13971712191373 * ANGSTROM_TO_BOHR, 1.66334132067180 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.67789462558235 * ANGSTROM_TO_BOHR,
                                             2.25447436615883 * ANGSTROM_TO_BOHR, 1.73243457296332 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 3.83499049328653 * ANGSTROM_TO_BOHR,
                                             0.37206802842858 * ANGSTROM_TO_BOHR, 1.12313605573798 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.37783744442278 * ANGSTROM_TO_BOHR,
                                             3.62666352141399 * ANGSTROM_TO_BOHR, 1.85690263210696 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.91036151348742 * ANGSTROM_TO_BOHR,
                                             5.70045562415378 * ANGSTROM_TO_BOHR, 1.71518016308150 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -6.97770079198632 * ANGSTROM_TO_BOHR,
                                             5.98523048128377 * ANGSTROM_TO_BOHR, -0.09496820690846 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.95298360793980 * ANGSTROM_TO_BOHR,
                                             7.43475978174949 * ANGSTROM_TO_BOHR, 0.94473892775111 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 6.06482570420256 * ANGSTROM_TO_BOHR,
                                             1.74288253485095 * ANGSTROM_TO_BOHR, 2.12536886797205 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 5.15981416390186 * ANGSTROM_TO_BOHR,
                                             4.28152009200965 * ANGSTROM_TO_BOHR, 2.22423216944500 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.86964434010808 * ANGSTROM_TO_BOHR,
                                             6.53401792490724 * ANGSTROM_TO_BOHR, 3.00458213575159 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.07700897782200 * ANGSTROM_TO_BOHR,
                                             8.45837761840221 * ANGSTROM_TO_BOHR, 1.98476882685460 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.14472552625358 * ANGSTROM_TO_BOHR,
                                             7.48659578913860 * ANGSTROM_TO_BOHR, 0.76272495082705 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -7.85288616769135 * ANGSTROM_TO_BOHR,
                                             4.84662340177242 * ANGSTROM_TO_BOHR, -0.64315267966354 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -6.96939516417516 * ANGSTROM_TO_BOHR,
                                             7.13739693024225 * ANGSTROM_TO_BOHR, -1.12404261102867 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -7.60570187618051 * ANGSTROM_TO_BOHR,
                                             6.49103497803841 * ANGSTROM_TO_BOHR, 1.22364270752792 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.53942859943231 * ANGSTROM_TO_BOHR,
                                             6.33080408864722 * ANGSTROM_TO_BOHR, 0.35147003340888 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.16466266696969 * ANGSTROM_TO_BOHR,
                                             6.43091570602525 * ANGSTROM_TO_BOHR, 2.22744559397163 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 1.78113064313056 * ANGSTROM_TO_BOHR,
                                             5.93517348180064 * ANGSTROM_TO_BOHR, 2.74345085677705 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 6.19671942485128 * ANGSTROM_TO_BOHR,
                                             0.22517036797559 * ANGSTROM_TO_BOHR, 1.92600610701894 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 6.32336504220137 * ANGSTROM_TO_BOHR,
                                             2.06430527402333 * ANGSTROM_TO_BOHR, 3.61486155230123 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 7.13465911933550 * ANGSTROM_TO_BOHR,
                                             2.44584604839223 * ANGSTROM_TO_BOHR, 1.26066024918145 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 6.97972138439825 * ANGSTROM_TO_BOHR,
                                             2.22658782615567 * ANGSTROM_TO_BOHR, 0.19841063261595 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 8.13569349180792 * ANGSTROM_TO_BOHR,
                                             2.10027931172002 * ANGSTROM_TO_BOHR, 1.54390569545958 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 7.10220740280281 * ANGSTROM_TO_BOHR,
                                             3.53268524604814 * ANGSTROM_TO_BOHR, 1.39061845600487 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 5.58046178985460 * ANGSTROM_TO_BOHR,
                                             1.57256023594689 * ANGSTROM_TO_BOHR, 4.25119023000150 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 6.27313806599647 * ANGSTROM_TO_BOHR,
                                             3.14122737756939 * ANGSTROM_TO_BOHR, 3.80473758713604 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 7.31944463992090 * ANGSTROM_TO_BOHR,
                                             1.71420487889887 * ANGSTROM_TO_BOHR, 3.90912043828194 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 5.47350076827894 * ANGSTROM_TO_BOHR,
                                             -0.32542835885976 * ANGSTROM_TO_BOHR, 2.53706657152166 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 7.20027372099028 * ANGSTROM_TO_BOHR,
                                             -0.09675294866989 * ANGSTROM_TO_BOHR, 2.22286828336657 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 6.05121988503054 * ANGSTROM_TO_BOHR,
                                             -0.05673782942222 * ANGSTROM_TO_BOHR, 0.87712347748801 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 5.00924700842936 * ANGSTROM_TO_BOHR,
                                             6.33263766517312 * ANGSTROM_TO_BOHR, 1.53500724397903 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 4.47493587876587 * ANGSTROM_TO_BOHR,
                                             6.06881270407874 * ANGSTROM_TO_BOHR, 3.21386341792396 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 3.93409070743025 * ANGSTROM_TO_BOHR,
                                             7.49735523722311 * ANGSTROM_TO_BOHR, 2.32222584505330 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.85813947996631 * ANGSTROM_TO_BOHR,
                                             5.44810386681995 * ANGSTROM_TO_BOHR, 2.43598675708472 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.59611883205988 * ANGSTROM_TO_BOHR,
                                             7.01099863310485 * ANGSTROM_TO_BOHR, 2.84637518583010 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.07276768245146 * ANGSTROM_TO_BOHR,
                                             5.54167647700564 * ANGSTROM_TO_BOHR, 3.72354895445304 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 3.36199034147681 * ANGSTROM_TO_BOHR,
                                             6.20328953933432 * ANGSTROM_TO_BOHR, -0.36221504683867 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.36102275862487 * ANGSTROM_TO_BOHR,
                                             7.40462201899642 * ANGSTROM_TO_BOHR, 0.48112838840598 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.64107671634069 * ANGSTROM_TO_BOHR,
                                             5.87516839522832 * ANGSTROM_TO_BOHR, -0.06465639765120 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.72345829093801 * ANGSTROM_TO_BOHR,
                                             6.37564505895787 * ANGSTROM_TO_BOHR, 3.67281663577421 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.17497582928836 * ANGSTROM_TO_BOHR,
                                             7.22765079366384 * ANGSTROM_TO_BOHR, 3.49299588383011 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.36541900120089 * ANGSTROM_TO_BOHR,
                                             5.57845911921942 * ANGSTROM_TO_BOHR, 2.85572229502848 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.93417751366491 * ANGSTROM_TO_BOHR,
                                             8.30219801250051 * ANGSTROM_TO_BOHR, 2.64964406646150 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.42467357855893 * ANGSTROM_TO_BOHR,
                                             8.96643772449778 * ANGSTROM_TO_BOHR, 1.07761713364915 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.37660789413876 * ANGSTROM_TO_BOHR,
                                             9.12854033136797 * ANGSTROM_TO_BOHR, 2.49452668299751 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.56571588586585 * ANGSTROM_TO_BOHR,
                                             6.59930345991454 * ANGSTROM_TO_BOHR, 0.51346683103427 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.49297490713285 * ANGSTROM_TO_BOHR,
                                             8.19630170799128 * ANGSTROM_TO_BOHR, 1.28536746599081 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.48722927395676 * ANGSTROM_TO_BOHR,
                                             7.96049735492584 * ANGSTROM_TO_BOHR, -0.16382389702834 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -7.45451709864170 * ANGSTROM_TO_BOHR,
                                             4.45016084269664 * ANGSTROM_TO_BOHR, -1.58238236773522 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -8.86353654253283 * ANGSTROM_TO_BOHR,
                                             5.22204729253269 * ANGSTROM_TO_BOHR, -0.83421954841187 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -7.93162796674867 * ANGSTROM_TO_BOHR,
                                             4.02120514252091 * ANGSTROM_TO_BOHR, 0.07341081374228 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -7.62254073277562 * ANGSTROM_TO_BOHR,
                                             5.69322442060897 * ANGSTROM_TO_BOHR, 1.97391298026946 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -8.63498784699360 * ANGSTROM_TO_BOHR,
                                             6.82324002028751 * ANGSTROM_TO_BOHR, 1.04666047722293 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -7.04446420275119 * ANGSTROM_TO_BOHR,
                                             7.33487223782684 * ANGSTROM_TO_BOHR, 1.63726927115536 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -6.55437340053619 * ANGSTROM_TO_BOHR,
                                             6.79989802962478 * ANGSTROM_TO_BOHR, -2.07878381014698 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -6.36578076522199 * ANGSTROM_TO_BOHR,
                                             7.98042136264129 * ANGSTROM_TO_BOHR, -0.77301152619500 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -7.99003809968693 * ANGSTROM_TO_BOHR,
                                             7.49691091696790 * ANGSTROM_TO_BOHR, -1.29796440932058 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), -0.58117457718709 * ANGSTROM_TO_BOHR,
                                             3.14840635123511 * ANGSTROM_TO_BOHR, -1.32876522814175 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.65465162911951 * ANGSTROM_TO_BOHR,
                                             5.55158682325514 * ANGSTROM_TO_BOHR, -3.03287591760466 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.45911998087859 * ANGSTROM_TO_BOHR,
                                             5.88136721211899 * ANGSTROM_TO_BOHR, -2.52777694945093 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.12013021748694 * ANGSTROM_TO_BOHR,
                                             4.20307019079305 * ANGSTROM_TO_BOHR, -3.37124888790066 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.38457193862791 * ANGSTROM_TO_BOHR,
                                             6.34299186426688 * ANGSTROM_TO_BOHR, -3.20134425367596 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.68960292441784 * ANGSTROM_TO_BOHR,
                                             5.14717912921101 * ANGSTROM_TO_BOHR, -2.30841396123747 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.22132572604776 * ANGSTROM_TO_BOHR,
                                             6.91415426531028 * ANGSTROM_TO_BOHR, -2.29632480920601 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.28423624804405 * ANGSTROM_TO_BOHR,
                                             3.07272920174484 * ANGSTROM_TO_BOHR, -3.32336077727346 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.45944839094449 * ANGSTROM_TO_BOHR,
                                             4.01776250494701 * ANGSTROM_TO_BOHR, -3.75157175786376 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.77984832806600 * ANGSTROM_TO_BOHR,
                                             1.81217555301830 * ANGSTROM_TO_BOHR, -3.63833651476025 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.24431275101592 * ANGSTROM_TO_BOHR,
                                             3.18161412634414 * ANGSTROM_TO_BOHR, -3.03321908971706 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.11887185642819 * ANGSTROM_TO_BOHR,
                                             1.64094202344557 * ANGSTROM_TO_BOHR, -4.00025880467642 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.11246484386345 * ANGSTROM_TO_BOHR,
                                             0.95425091922413 * ANGSTROM_TO_BOHR, -3.61091891825824 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.95851314964250 * ANGSTROM_TO_BOHR,
                                             2.75334759601635 * ANGSTROM_TO_BOHR, -4.05598131498691 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -6.11541919341662 * ANGSTROM_TO_BOHR,
                                             4.88341492242030 * ANGSTROM_TO_BOHR, -3.80340449371217 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -5.49899366110141 * ANGSTROM_TO_BOHR,
                                             0.65353364283191 * ANGSTROM_TO_BOHR, -4.24642233726737 * ANGSTROM_TO_BOHR));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -7.00028167101145 * ANGSTROM_TO_BOHR,
                                             2.63804155947607 * ANGSTROM_TO_BOHR, -4.34247362306615 * ANGSTROM_TO_BOHR));
      _testGeometries[TEST_SYSTEM_CONTROLLERS::JACOBSEN_MINBAS] = std::make_shared<Geometry>(atoms);
    } break;
    case TEST_SYSTEM_CONTROLLERS::C60_MINBAS: {
      std::vector<std::shared_ptr<Atom>> atoms;
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.09836559756361, 1.10514018840352,
                                             4.90541408655266));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 5.74488397714688, 0.32279356938886,
                                             3.02400274865171));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.42047776417445, -0.59636922164672,
                                             6.00240010669277));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 5.69575109769185, -2.17164492602019,
                                             2.21822352558927));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 5.86885001146417, -2.71985447717031,
                                             -0.34878045332258));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 5.95747816709651, 2.27677039079262,
                                             1.27865169231853));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 6.13227783438843, 1.71797837329833,
                                             -1.28022646422192));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 6.08768029765233, -0.77116888893864,
                                             -2.10187938680060));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.82588275950980, 2.54530047427567,
                                             5.92529928047104));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.82876144244487, 4.26325050168200,
                                             4.84267517894085));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.05166196286655, 0.10868759853292,
                                             6.51262616257190));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 3.30033425164601, 3.53986333801338,
                                             4.33018145170227));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.44872082260061, 4.27421091325274,
                                             2.09274571036567));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.49329096012192, 5.70889099333952,
                                             3.10110357487345));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.64186752790248, 6.42434130725003,
                                             0.84231392823524));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 3.12534561174080, 5.71078071947241,
                                             0.33851294120792));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.94582060490399, -1.56523180997719,
                                             4.92072086822903));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.74101736162227, 0.87043620269912,
                                             4.33471679442119));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.59235567900821, -1.94827929711310,
                                             6.01354949087680));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.18444994596437, 2.92437953653253,
                                             4.83700600054220));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.98961918166386, 4.87552176873695,
                                             3.09978076658043));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -6.10218709513980, 0.75648571688612,
                                             2.09595824479157));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.89677386449513, 2.71216329180947,
                                             0.34588287312618));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.33208062646582, 4.76837429700234,
                                             0.84722721618074));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.92320365535341, -5.51627120861449,
                                             3.29744612008026));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.28951871895281, -5.13662522851776,
                                             2.21501099116337));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.06999230635554, -3.91964159893941,
                                             5.20380184293530));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.81188209160548, -3.15732607693334,
                                             3.02853809137064));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -6.15661120776691, -1.73833072374951,
                                             1.28224217197101));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.12190001096585, -5.69163779374628,
                                             -0.34254435708406));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.44451933137251, -4.25128853526084,
                                             -2.09488740010892));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.97538647162317, -2.26896582186380,
                                             -1.28381694387440));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.04148484096376, -3.89488618659861,
                                             3.29612331178724));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 3.19696623217716, -5.49283860456671,
                                             1.39411395903784));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.42028879156116, -3.10365785475939,
                                             5.20040033589611));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.69931520234221, -6.30882234874673,
                                             1.39581471255744));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.64579185904579, -6.41785954661423,
                                             -0.85220349512332));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.32966807622882, -4.78249055121499,
                                             -0.86089623533459));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.98815149449330, -4.87092973423404,
                                             -3.11722923800005));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.48804382068559, -5.69031498545325,
                                             -3.10796957994891));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.12062027219591, -1.10376068832651,
                                             -4.89338283029927));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.29896734961724, -3.54358609849517,
                                             -4.34611814221559));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.82437097860349, -4.25695771365949,
                                             -4.86806050011860));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.41741010862608, 0.59170159809849,
                                             -5.96901494513778));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.06355054298277, 3.87151027433481,
                                             -3.26802938340433));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.75655618543502, 2.16886702860485,
                                             -2.21601884522689));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -5.78887050230736, -0.31876845272581,
                                             -3.02784519191456));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.38830832617965, 3.08954160054673,
                                             -5.13753544666809));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.95650692902794, 5.53749283308679,
                                             -3.27558828793587));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.71155432847021, 6.33117780889876,
                                             -1.40381455331264));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.19748905628128, 5.50442262576130,
                                             -1.39965715582029));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.10272866219021, 3.90647020779320,
                                             -5.15340914618433));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 3.96816346700779, 1.56509952914788,
                                             -4.90698885845604));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.82628810395116, 3.18705146900363,
                                             -3.02784519191456));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 3.32811322579943, 5.17674411431892,
                                             -2.21847548919964));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 1.59542333455658, 1.91772242554434,
                                             -5.99074679566596));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.84728075854715, -2.56111748200792,
                                             -5.96712521900489));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 3.21359582214656, -2.93717298245216,
                                             -4.87996577475578));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 4.76336022372607, -0.88455245691178,
                                             -4.36274773218498));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.04546996036377, -0.13281940124987,
                                             -6.51968114026065));
      _testGeometries[TEST_SYSTEM_CONTROLLERS::C60_MINBAS] = std::make_shared<Geometry>(atoms);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H60_Ghost_MINBAS: {
      std::vector<std::shared_ptr<Atom>> atoms;
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 4.09836559756361, 1.10514018840352,
                                             4.90541408655266));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 5.74488397714688, 0.32279356938886,
                                             3.02400274865171));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 2.42047776417445, -0.59636922164672,
                                             6.00240010669277));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 5.69575109769185, -2.17164492602019,
                                             2.21822352558927));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 5.86885001146417, -2.71985447717031,
                                             -0.34878045332258));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 5.95747816709651, 2.27677039079262,
                                             1.27865169231853));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 6.13227783438843, 1.71797837329833,
                                             -1.28022646422192));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 6.08768029765233, -0.77116888893864,
                                             -2.10187938680060));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -0.82588275950980, 2.54530047427567,
                                             5.92529928047104));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.82876144244487, 4.26325050168200,
                                             4.84267517894085));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -0.05166196286655, 0.10868759853292,
                                             6.51262616257190));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 3.30033425164601, 3.53986333801338,
                                             4.33018145170227));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 4.44872082260061, 4.27421091325274,
                                             2.09274571036567));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -0.49329096012192, 5.70889099333952,
                                             3.10110357487345));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.64186752790248, 6.42434130725003,
                                             0.84231392823524));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 3.12534561174080, 5.71078071947241,
                                             0.33851294120792));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -3.94582060490399, -1.56523180997719,
                                             4.92072086822903));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -4.74101736162227, 0.87043620269912,
                                             4.33471679442119));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -1.59235567900821, -1.94827929711310,
                                             6.01354949087680));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -3.18444994596437, 2.92437953653253,
                                             4.83700600054220));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -2.98961918166386, 4.87552176873695,
                                             3.09978076658043));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -6.10218709513980, 0.75648571688612,
                                             2.09595824479157));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -5.89677386449513, 2.71216329180947,
                                             0.34588287312618));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -4.33208062646582, 4.76837429700234,
                                             0.84722721618074));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -0.92320365535341, -5.51627120861449,
                                             3.29744612008026));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -3.28951871895281, -5.13662522851776,
                                             2.21501099116337));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -0.06999230635554, -3.91964159893941,
                                             5.20380184293530));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -4.81188209160548, -3.15732607693334,
                                             3.02853809137064));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -6.15661120776691, -1.73833072374951,
                                             1.28224217197101));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -3.12190001096585, -5.69163779374628,
                                             -0.34254435708406));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -4.44451933137251, -4.25128853526084,
                                             -2.09488740010892));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -5.97538647162317, -2.26896582186380,
                                             -1.28381694387440));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 4.04148484096376, -3.89488618659861,
                                             3.29612331178724));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 3.19696623217716, -5.49283860456671,
                                             1.39411395903784));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 2.42028879156116, -3.10365785475939,
                                             5.20040033589611));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.69931520234221, -6.30882234874673,
                                             1.39581471255744));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -0.64579185904579, -6.41785954661423,
                                             -0.85220349512332));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 4.32966807622882, -4.78249055121499,
                                             -0.86089623533459));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 2.98815149449330, -4.87092973423404,
                                             -3.11722923800005));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.48804382068559, -5.69031498545325,
                                             -3.10796957994891));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -4.12062027219591, -1.10376068832651,
                                             -4.89338283029927));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -3.29896734961724, -3.54358609849517,
                                             -4.34611814221559));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -0.82437097860349, -4.25695771365949,
                                             -4.86806050011860));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -2.41741010862608, 0.59170159809849,
                                             -5.96901494513778));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -4.06355054298277, 3.87151027433481,
                                             -3.26802938340433));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -5.75655618543502, 2.16886702860485,
                                             -2.21601884522689));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -5.78887050230736, -0.31876845272581,
                                             -3.02784519191456));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -2.38830832617965, 3.08954160054673,
                                             -5.13753544666809));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.95650692902794, 5.53749283308679,
                                             -3.27558828793587));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -0.71155432847021, 6.33117780889876,
                                             -1.40381455331264));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -3.19748905628128, 5.50442262576130,
                                             -1.39965715582029));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.10272866219021, 3.90647020779320,
                                             -5.15340914618433));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 3.96816346700779, 1.56509952914788,
                                             -4.90698885845604));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 4.82628810395116, 3.18705146900363,
                                             -3.02784519191456));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 3.32811322579943, 5.17674411431892,
                                             -2.21847548919964));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 1.59542333455658, 1.91772242554434,
                                             -5.99074679566596));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.84728075854715, -2.56111748200792,
                                             -5.96712521900489));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 3.21359582214656, -2.93717298245216,
                                             -4.87996577475578));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 4.76336022372607, -0.88455245691178,
                                             -4.36274773218498));
      atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.04546996036377, -0.13281940124987,
                                             -6.51968114026065));
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H60_Ghost_MINBAS] = std::make_shared<Geometry>(atoms);
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 1.7);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 2.4);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_B3LYP:
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86:
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_BP86:
    case TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_PBE0_UNRESTRICTED: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 1.7);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 2.4);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_B3LYP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_PBE0_UNRESTRICTED] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86_Supermolecular: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, -0.7);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.7);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.0, 0.0, 1.7);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.0, 0.0, 2.4);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86_Supermolecular] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, H3, H4});
    } break;
    case TEST_SYSTEM_CONTROLLERS::He2_6_31Gs_BP86:
    case TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86: {
      auto He1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("He"), 0.0, 0.0, 0.0);
      auto He2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("He"), 0.0, 0.0, 1.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He2_6_31Gs_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He1, He2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He1, He2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::HeNe_def2SVP_BP86: {
      auto He1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("He"), 0.0, 0.0, 0.0);
      auto Ne2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Ne"), 0.0, 0.0, 2.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::HeNe_def2SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He1, Ne2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86:
    case TEST_SYSTEM_CONTROLLERS::He_1_def2SVP_BP86: {
      auto He = std::make_shared<Atom>(AtomTypeFactory::getAtomType("He"), 0.0, 0.0, 0.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He_1_def2SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He});
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_2_6_31Gs_BP86:
    case TEST_SYSTEM_CONTROLLERS::He_2_def2SVP_BP86: {
      auto He = std::make_shared<Atom>(AtomTypeFactory::getAtomType("He"), 0.0, 0.0, 1.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He_2_6_31Gs_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He_2_def2SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He});
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_3_6_31Gs_BP86:
    case TEST_SYSTEM_CONTROLLERS::He_3_def2SVP_BP86: {
      auto He = std::make_shared<Atom>(AtomTypeFactory::getAtomType("He"), 0.0, 0.0, 2.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He_3_6_31Gs_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He_3_def2SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He});
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_MINBAS_SING: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 0);
      auto O2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 1.21 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::O2_MINBAS_SING] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, O2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP:
    case TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP_CIS:
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDA:
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDHF:
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDDFT:
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_6_31G_PBE_TDA:
    case TEST_SYSTEM_CONTROLLERS::O2_TRIP_6_31G_PBE_TDDFT: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 0);
      auto O2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 1.21 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, O2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP_CIS] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, O2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDA] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, O2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDHF] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, O2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDDFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, O2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::O2_TRIP_6_31G_PBE_TDA] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, O2});
      _testGeometries[TEST_SYSTEM_CONTROLLERS::O2_TRIP_6_31G_PBE_TDDFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, O2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs: {
      auto F = std::make_shared<Atom>(AtomTypeFactory::getAtomType("F"), 0.0, 0.0, 0.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{F});
    } break;
    case TEST_SYSTEM_CONTROLLERS::OH_MINBAS_PBE: {
      auto O = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.0, 0.0, 0.0);
      auto H = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.97);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::OH_MINBAS_PBE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O, H});
    } break;
    case TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE: {
      auto C = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -4.55140536955085 * ANGSTROM_TO_BOHR,
                                      2.81558128447117 * ANGSTROM_TO_BOHR, -0.04650481068098 * ANGSTROM_TO_BOHR);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.53420848990172 * ANGSTROM_TO_BOHR,
                                       3.16961934366078 * ANGSTROM_TO_BOHR, 0.01972807283118 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -5.06110516479103 * ANGSTROM_TO_BOHR,
                                       2.46787734641333 * ANGSTROM_TO_BOHR, 0.83903751909563 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -5.06106097575639 * ANGSTROM_TO_BOHR,
                                       2.81228202545472 * ANGSTROM_TO_BOHR, -0.99807078124583 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C, H1, H2, H3});
    } break;
    case TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE: {
      auto I = std::make_shared<Atom>(AtomTypeFactory::getAtomType("I"), -5.79582 * ANGSTROM_TO_BOHR,
                                      3.26985 * ANGSTROM_TO_BOHR, 0.0);
      auto H = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.27766 * ANGSTROM_TO_BOHR,
                                      4.08963 * ANGSTROM_TO_BOHR, 0.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{I, H});
    } break;
    case TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE: {
      auto I1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("I"), -1.93122 * ANGSTROM_TO_BOHR,
                                       1.08036 * ANGSTROM_TO_BOHR, 0.0);
      auto I2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("I"), -0.86611 * ANGSTROM_TO_BOHR,
                                       3.63071 * ANGSTROM_TO_BOHR, 0.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{I1, I2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LINEARDEPENDENT: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, -0.370424 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.370424 * ANGSTROM_TO_BOHR);
      auto HGhost = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), 0.0, 0.0, -0.370424 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LINEARDEPENDENT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, HGhost});
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMon1_3_DFT: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.9817439 * ANGSTROM_TO_BOHR,
                                       0.875178014 * ANGSTROM_TO_BOHR, 0.0365091960 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.63911744 * ANGSTROM_TO_BOHR,
                                       0.066884843 * ANGSTROM_TO_BOHR, 0.021624441 * ANGSTROM_TO_BOHR);
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 1.60679538 * ANGSTROM_TO_BOHR,
                                       -0.066291721 * ANGSTROM_TO_BOHR, 0.087598426 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMon1_3_DFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMon2_3_DFT: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.557299588 * ANGSTROM_TO_BOHR,
                                       -0.318709677 * ANGSTROM_TO_BOHR, 0.796158231 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.48665037 * ANGSTROM_TO_BOHR,
                                       -0.493436911 * ANGSTROM_TO_BOHR, -0.764205635 * ANGSTROM_TO_BOHR);
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), -1.183706762 * ANGSTROM_TO_BOHR,
                                       0.042210901 * ANGSTROM_TO_BOHR, -0.071849209 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMon2_3_DFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterMon3_3_DFT: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.30929 * ANGSTROM_TO_BOHR,
                                       -3.46619 * ANGSTROM_TO_BOHR, 0.59972 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.41336 * ANGSTROM_TO_BOHR,
                                       -2.01591 * ANGSTROM_TO_BOHR, -0.02871 * ANGSTROM_TO_BOHR);
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 2.79197 * ANGSTROM_TO_BOHR,
                                       -2.90897 * ANGSTROM_TO_BOHR, -0.03069 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterMon3_3_DFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::WaterDim_DFT: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.9817439 * ANGSTROM_TO_BOHR,
                                       0.875178014 * ANGSTROM_TO_BOHR, 0.0365091960 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.63911744 * ANGSTROM_TO_BOHR,
                                       0.066884843 * ANGSTROM_TO_BOHR, 0.021624441 * ANGSTROM_TO_BOHR);
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 1.60679538 * ANGSTROM_TO_BOHR,
                                       -0.066291721 * ANGSTROM_TO_BOHR, 0.087598426 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.557299588 * ANGSTROM_TO_BOHR,
                                       -0.318709677 * ANGSTROM_TO_BOHR, 0.796158231 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.48665037 * ANGSTROM_TO_BOHR,
                                       -0.493436911 * ANGSTROM_TO_BOHR, -0.764205635 * ANGSTROM_TO_BOHR);
      auto O2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), -1.183706762 * ANGSTROM_TO_BOHR,
                                       0.042210901 * ANGSTROM_TO_BOHR, -0.071849209 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WaterDim_DFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, O1, H3, H4, O2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.86799 * ANGSTROM_TO_BOHR,
                                       1.86633 * ANGSTROM_TO_BOHR, 0.00000);
      auto C2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.662850 * ANGSTROM_TO_BOHR,
                                       2.791330 * ANGSTROM_TO_BOHR, 0.00000);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -4.69801 * ANGSTROM_TO_BOHR,
                                       2.32927 * ANGSTROM_TO_BOHR, 0.57419 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -3.59828 * ANGSTROM_TO_BOHR,
                                       0.89649 * ANGSTROM_TO_BOHR, 0.46861 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -4.20382 * ANGSTROM_TO_BOHR,
                                       1.68541 * ANGSTROM_TO_BOHR, -1.04280 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -2.32703 * ANGSTROM_TO_BOHR,
                                       2.97225 * ANGSTROM_TO_BOHR, 1.04280 * ANGSTROM_TO_BOHR);
      auto H5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -2.93257 * ANGSTROM_TO_BOHR,
                                       3.76118 * ANGSTROM_TO_BOHR, -0.46861 * ANGSTROM_TO_BOHR);
      auto H6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H:"), -1.83284 * ANGSTROM_TO_BOHR,
                                       2.32840 * ANGSTROM_TO_BOHR, -0.57419 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, C2, H1, H2, H3, H4, H5, H6});
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Env: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C:"), -3.86799 * ANGSTROM_TO_BOHR,
                                       1.86633 * ANGSTROM_TO_BOHR, 0.00000);
      auto C2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C:"), -2.662850 * ANGSTROM_TO_BOHR,
                                       2.791330 * ANGSTROM_TO_BOHR, 0.00000);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.69801 * ANGSTROM_TO_BOHR,
                                       2.32927 * ANGSTROM_TO_BOHR, 0.57419 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.59828 * ANGSTROM_TO_BOHR,
                                       0.89649 * ANGSTROM_TO_BOHR, 0.46861 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.20382 * ANGSTROM_TO_BOHR,
                                       1.68541 * ANGSTROM_TO_BOHR, -1.04280 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.32703 * ANGSTROM_TO_BOHR,
                                       2.97225 * ANGSTROM_TO_BOHR, 1.04280 * ANGSTROM_TO_BOHR);
      auto H5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.93257 * ANGSTROM_TO_BOHR,
                                       3.76118 * ANGSTROM_TO_BOHR, -0.46861 * ANGSTROM_TO_BOHR);
      auto H6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.83284 * ANGSTROM_TO_BOHR,
                                       2.32840 * ANGSTROM_TO_BOHR, -0.57419 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Env] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, C2, H1, H2, H3, H4, H5, H6});
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.86799 * ANGSTROM_TO_BOHR,
                                       1.86633 * ANGSTROM_TO_BOHR, 0.00000);
      auto C2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.662850 * ANGSTROM_TO_BOHR,
                                       2.791330 * ANGSTROM_TO_BOHR, 0.00000);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.69801 * ANGSTROM_TO_BOHR,
                                       2.32927 * ANGSTROM_TO_BOHR, 0.57419 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.59828 * ANGSTROM_TO_BOHR,
                                       0.89649 * ANGSTROM_TO_BOHR, 0.46861 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.20382 * ANGSTROM_TO_BOHR,
                                       1.68541 * ANGSTROM_TO_BOHR, -1.04280 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.32703 * ANGSTROM_TO_BOHR,
                                       2.97225 * ANGSTROM_TO_BOHR, 1.04280 * ANGSTROM_TO_BOHR);
      auto H5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.93257 * ANGSTROM_TO_BOHR,
                                       3.76118 * ANGSTROM_TO_BOHR, -0.46861 * ANGSTROM_TO_BOHR);
      auto H6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.83284 * ANGSTROM_TO_BOHR,
                                       2.32840 * ANGSTROM_TO_BOHR, -0.57419 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, C2, H1, H2, H3, H4, H5, H6});
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.86799 * ANGSTROM_TO_BOHR,
                                       1.86633 * ANGSTROM_TO_BOHR, 0.00000);
      auto C2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.662850 * ANGSTROM_TO_BOHR,
                                       2.791330 * ANGSTROM_TO_BOHR, 0.00000);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.69801 * ANGSTROM_TO_BOHR,
                                       2.32927 * ANGSTROM_TO_BOHR, 0.57419 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.59828 * ANGSTROM_TO_BOHR,
                                       0.89649 * ANGSTROM_TO_BOHR, 0.46861 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.20382 * ANGSTROM_TO_BOHR,
                                       1.68541 * ANGSTROM_TO_BOHR, -1.04280 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.32703 * ANGSTROM_TO_BOHR,
                                       2.97225 * ANGSTROM_TO_BOHR, 1.04280 * ANGSTROM_TO_BOHR);
      auto H5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.93257 * ANGSTROM_TO_BOHR,
                                       3.76118 * ANGSTROM_TO_BOHR, -0.46861 * ANGSTROM_TO_BOHR);
      auto H6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.83284 * ANGSTROM_TO_BOHR,
                                       2.32840 * ANGSTROM_TO_BOHR, -0.57419 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, C2, H1, H2, H3, H4, H5, H6});
    } break;
    case TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -3.86799 * ANGSTROM_TO_BOHR,
                                       1.86633 * ANGSTROM_TO_BOHR, 0.00000);
      auto C2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.207840 * ANGSTROM_TO_BOHR,
                                       3.067080 * ANGSTROM_TO_BOHR, 0.052300 * ANGSTROM_TO_BOHR);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.69801 * ANGSTROM_TO_BOHR,
                                       2.32927 * ANGSTROM_TO_BOHR, 0.57419 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.59828 * ANGSTROM_TO_BOHR,
                                       0.89649 * ANGSTROM_TO_BOHR, 0.46861 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -4.20382 * ANGSTROM_TO_BOHR,
                                       1.68541 * ANGSTROM_TO_BOHR, -1.04280 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.872020 * ANGSTROM_TO_BOHR,
                                       3.248000 * ANGSTROM_TO_BOHR, 1.095100 * ANGSTROM_TO_BOHR);
      auto H5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.477560 * ANGSTROM_TO_BOHR,
                                       4.036930 * ANGSTROM_TO_BOHR, -0.416310 * ANGSTROM_TO_BOHR);
      auto H6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.377830 * ANGSTROM_TO_BOHR,
                                       2.604150 * ANGSTROM_TO_BOHR, -0.521890 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, C2, H1, H2, H3, H4, H5, H6});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * 1.0686,
                                       ANGSTROM_TO_BOHR * -0.1411, ANGSTROM_TO_BOHR * 1.0408);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), ANGSTROM_TO_BOHR * 0.5979,
                                       ANGSTROM_TO_BOHR * 0.0151, ANGSTROM_TO_BOHR * 0.0688);
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * 1.2687,
                                       ANGSTROM_TO_BOHR * 0.2002, ANGSTROM_TO_BOHR * -0.7717);
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), ANGSTROM_TO_BOHR * -0.5960,
                                       ANGSTROM_TO_BOHR * -0.0151, ANGSTROM_TO_BOHR * -0.0686);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, C1, O1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), ANGSTROM_TO_BOHR * -0.0000000,
                                       ANGSTROM_TO_BOHR * 0.0000000, ANGSTROM_TO_BOHR * -0.6046693);
      auto N1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("N"), ANGSTROM_TO_BOHR * 0.7174854,
                                       ANGSTROM_TO_BOHR * -0.1695956, ANGSTROM_TO_BOHR * 0.6268497);
      auto N2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("N"), ANGSTROM_TO_BOHR * -0.7174854,
                                       ANGSTROM_TO_BOHR * 0.1695956, ANGSTROM_TO_BOHR * 0.6268497);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * -0.1333274,
                                       ANGSTROM_TO_BOHR * -0.9067874, ANGSTROM_TO_BOHR * -1.1877275);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * 0.1333274,
                                       ANGSTROM_TO_BOHR * 0.9067874, ANGSTROM_TO_BOHR * -1.1877275);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * -1.1850131,
                                       ANGSTROM_TO_BOHR * -0.7019020, ANGSTROM_TO_BOHR * 0.8631353);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * 1.1850131,
                                       ANGSTROM_TO_BOHR * 0.7019020, ANGSTROM_TO_BOHR * 0.8631353);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, N1, N2, H1, H2, H3, H4});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_HYBRID: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, -0.7);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.7);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_HYBRID] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_HYBRID: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 1.7);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 2.4);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_HYBRID] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2Dimer_Def2_TZVP: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, -0.7);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.7);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 1.7);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 2.4);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2Dimer_Def2_TZVP] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2, H3, H4});
    } break;
    case TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 2.300596 * ANGSTROM_TO_BOHR,
                                       -0.110893 * ANGSTROM_TO_BOHR, 0.000128 * ANGSTROM_TO_BOHR);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.680547 * ANGSTROM_TO_BOHR,
                                       0.774789 * ANGSTROM_TO_BOHR, 0.000232 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.218552 * ANGSTROM_TO_BOHR,
                                       -1.219283 * ANGSTROM_TO_BOHR, 0.000014 * ANGSTROM_TO_BOHR);
      auto C2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.168651 * ANGSTROM_TO_BOHR,
                                       -1.184791 * ANGSTROM_TO_BOHR, -0.000069 * ANGSTROM_TO_BOHR);
      auto C3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.850797 * ANGSTROM_TO_BOHR,
                                       0.028852 * ANGSTROM_TO_BOHR, -0.000109 * ANGSTROM_TO_BOHR);
      auto C4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.127187 * ANGSTROM_TO_BOHR,
                                       1.214896 * ANGSTROM_TO_BOHR, -0.000067 * ANGSTROM_TO_BOHR);
      auto C5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.263118 * ANGSTROM_TO_BOHR,
                                       1.193734 * ANGSTROM_TO_BOHR, 0.000013 * ANGSTROM_TO_BOHR);
      auto C6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.936484 * ANGSTROM_TO_BOHR,
                                       -0.025561 * ANGSTROM_TO_BOHR, 0.000056 * ANGSTROM_TO_BOHR);

      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.758987 * ANGSTROM_TO_BOHR,
                                       -2.156670 * ANGSTROM_TO_BOHR, 0.000050 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.721496 * ANGSTROM_TO_BOHR,
                                       -2.115963 * ANGSTROM_TO_BOHR, -0.000101 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.932538 * ANGSTROM_TO_BOHR,
                                       0.048567 * ANGSTROM_TO_BOHR, -0.000172 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.643317 * ANGSTROM_TO_BOHR,
                                       2.166850 * ANGSTROM_TO_BOHR, -0.000098 * ANGSTROM_TO_BOHR);
      auto H5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.823932 * ANGSTROM_TO_BOHR,
                                       2.122490 * ANGSTROM_TO_BOHR, 0.000040 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, C2, C3, C4, C5, C6, H1, H2, H3, H4, H5});
    } break;
    case TEST_SYSTEM_CONTROLLERS::PHENOLATE_O_DEF2_SVP_BP86: {
      auto O = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 2.341566 * ANGSTROM_TO_BOHR, 0.0,
                                      -0.000016 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::PHENOLATE_O_DEF2_SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O});
    } break;
    case TEST_SYSTEM_CONTROLLERS::PHENOLATE_PHENYL_DEF2_SVP_BP86: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 1.077309 * ANGSTROM_TO_BOHR, 0.0,
                                       -0.000014 * ANGSTROM_TO_BOHR);
      auto C2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.286113 * ANGSTROM_TO_BOHR,
                                       1.209101 * ANGSTROM_TO_BOHR, -0.000003 * ANGSTROM_TO_BOHR);
      auto C3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.097619 * ANGSTROM_TO_BOHR,
                                       1.196760 * ANGSTROM_TO_BOHR, 0.000012 * ANGSTROM_TO_BOHR);
      auto C4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.823396 * ANGSTROM_TO_BOHR, 0.0,
                                       0.000013 * ANGSTROM_TO_BOHR);
      auto C5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.097619 * ANGSTROM_TO_BOHR,
                                       -1.196759 * ANGSTROM_TO_BOHR, 0.000010 * ANGSTROM_TO_BOHR);
      auto C6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.286113 * ANGSTROM_TO_BOHR,
                                       -1.209102 * ANGSTROM_TO_BOHR, -0.000005 * ANGSTROM_TO_BOHR);

      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.829304 * ANGSTROM_TO_BOHR,
                                       2.149743 * ANGSTROM_TO_BOHR, -0.000005 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.634545 * ANGSTROM_TO_BOHR,
                                       2.143825 * ANGSTROM_TO_BOHR, -.000017 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.907447 * ANGSTROM_TO_BOHR, 0.0,
                                       0.000035 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.634545 * ANGSTROM_TO_BOHR,
                                       -2.143825 * ANGSTROM_TO_BOHR, -0.000013 * ANGSTROM_TO_BOHR);
      auto H5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.829303 * ANGSTROM_TO_BOHR,
                                       -2.149744 * ANGSTROM_TO_BOHR, -0.000009 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::PHENOLATE_PHENYL_DEF2_SVP_BP86] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, C2, C3, C4, C5, C6, H1, H2, H3, H4, H5});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.0);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 3.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 3.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_BS_3A: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.0);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 3.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_BSAct_3A: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 0.0);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::H2_6_311G_BSEnv_3A: {
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0, 0.0, 3.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{H1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::MethylRad_Act_def2_SVP_PBE: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.13013 * ANGSTROM_TO_BOHR,
                                       -0.38860 * ANGSTROM_TO_BOHR, 0.00000);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.08186 * ANGSTROM_TO_BOHR,
                                       -0.11833 * ANGSTROM_TO_BOHR, 0.00000);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.88833 * ANGSTROM_TO_BOHR,
                                       0.38409 * ANGSTROM_TO_BOHR, 0.00000);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.42021 * ANGSTROM_TO_BOHR,
                                       -1.43156 * ANGSTROM_TO_BOHR, 0.00000);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::MethylRad_Act_def2_SVP_PBE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, H1, H2, H3});

    } break;
    case TEST_SYSTEM_CONTROLLERS::MethylRad_Env_def2_SVP_PBE: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.13013 * ANGSTROM_TO_BOHR,
                                       -0.38860 * ANGSTROM_TO_BOHR, 1.50000 * ANGSTROM_TO_BOHR);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.08186 * ANGSTROM_TO_BOHR,
                                       -0.11833 * ANGSTROM_TO_BOHR, 1.50000 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.88833 * ANGSTROM_TO_BOHR,
                                       0.38409 * ANGSTROM_TO_BOHR, 1.50000 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.42021 * ANGSTROM_TO_BOHR,
                                       -1.43156 * ANGSTROM_TO_BOHR, 1.50000 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::MethylRad_Env_def2_SVP_PBE] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, H1, H2, H3});
    } break;
    case TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF: {
      auto Pt1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Pt"), 1.208937 * ANGSTROM_TO_BOHR,
                                        -0.205410 * ANGSTROM_TO_BOHR, 499.683201 * ANGSTROM_TO_BOHR);
      auto N2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("N"), -0.611985 * ANGSTROM_TO_BOHR,
                                       -1.191922 * ANGSTROM_TO_BOHR, 500.382538 * ANGSTROM_TO_BOHR);
      auto N3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("N"), -0.372363 * ANGSTROM_TO_BOHR,
                                       1.292367 * ANGSTROM_TO_BOHR, 499.600723 * ANGSTROM_TO_BOHR);
      auto C4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 1.776879 * ANGSTROM_TO_BOHR,
                                       0.240563 * ANGSTROM_TO_BOHR, 501.576745 * ANGSTROM_TO_BOHR);
      auto C5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.745080 * ANGSTROM_TO_BOHR,
                                       0.891072 * ANGSTROM_TO_BOHR, 498.952272 * ANGSTROM_TO_BOHR);
      auto H6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.088805 * ANGSTROM_TO_BOHR,
                                       -2.574747 * ANGSTROM_TO_BOHR, 500.321752 * ANGSTROM_TO_BOHR);
      auto C7 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 2.513482 * ANGSTROM_TO_BOHR,
                                       -1.755001 * ANGSTROM_TO_BOHR, 499.738086 * ANGSTROM_TO_BOHR);
      auto H8 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.869921 * ANGSTROM_TO_BOHR,
                                       5.964618 * ANGSTROM_TO_BOHR, 497.506418 * ANGSTROM_TO_BOHR);
      auto H9 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.261961 * ANGSTROM_TO_BOHR,
                                       -5.548917 * ANGSTROM_TO_BOHR, 499.143239 * ANGSTROM_TO_BOHR);
      auto C10 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.532296 * ANGSTROM_TO_BOHR,
                                        0.946906 * ANGSTROM_TO_BOHR, 500.016951 * ANGSTROM_TO_BOHR);
      auto C11 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -1.667337 * ANGSTROM_TO_BOHR,
                                        -0.473185 * ANGSTROM_TO_BOHR, 500.473289 * ANGSTROM_TO_BOHR);
      auto C12 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.727990 * ANGSTROM_TO_BOHR,
                                        1.822862 * ANGSTROM_TO_BOHR, 500.017209 * ANGSTROM_TO_BOHR);
      auto C13 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -2.989961 * ANGSTROM_TO_BOHR,
                                        -0.961388 * ANGSTROM_TO_BOHR, 500.930824 * ANGSTROM_TO_BOHR);
      auto C14 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.113937 * ANGSTROM_TO_BOHR,
                                        2.570757 * ANGSTROM_TO_BOHR, 499.067771 * ANGSTROM_TO_BOHR);
      auto C15 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.034663 * ANGSTROM_TO_BOHR,
                                        3.702457 * ANGSTROM_TO_BOHR, 499.867189 * ANGSTROM_TO_BOHR);
      auto C16 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.388081 * ANGSTROM_TO_BOHR,
                                        4.920172 * ANGSTROM_TO_BOHR, 499.308350 * ANGSTROM_TO_BOHR);
      auto C17 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.593977 * ANGSTROM_TO_BOHR,
                                        5.012859 * ANGSTROM_TO_BOHR, 497.942446 * ANGSTROM_TO_BOHR);
      auto C18 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.454000 * ANGSTROM_TO_BOHR,
                                        3.900511 * ANGSTROM_TO_BOHR, 497.128117 * ANGSTROM_TO_BOHR);
      auto C19 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 0.109726 * ANGSTROM_TO_BOHR,
                                        2.688153 * ANGSTROM_TO_BOHR, 497.696195 * ANGSTROM_TO_BOHR);
      auto C20 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.603252 * ANGSTROM_TO_BOHR,
                                        -2.560162 * ANGSTROM_TO_BOHR, 500.715080 * ANGSTROM_TO_BOHR);
      auto C21 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.461467 * ANGSTROM_TO_BOHR,
                                        -3.494772 * ANGSTROM_TO_BOHR, 499.688774 * ANGSTROM_TO_BOHR);
      auto C22 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.360413 * ANGSTROM_TO_BOHR,
                                        -4.846531 * ANGSTROM_TO_BOHR, 499.960309 * ANGSTROM_TO_BOHR);
      auto C23 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.383252 * ANGSTROM_TO_BOHR,
                                        -5.275235 * ANGSTROM_TO_BOHR, 501.277759 * ANGSTROM_TO_BOHR);
      auto C24 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.512104 * ANGSTROM_TO_BOHR,
                                        -4.367401 * ANGSTROM_TO_BOHR, 502.314870 * ANGSTROM_TO_BOHR);
      auto C25 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), -0.623340 * ANGSTROM_TO_BOHR,
                                        -3.015229 * ANGSTROM_TO_BOHR, 502.033033 * ANGSTROM_TO_BOHR);
      auto H26 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.781969 * ANGSTROM_TO_BOHR,
                                        0.652981 * ANGSTROM_TO_BOHR, 501.558667 * ANGSTROM_TO_BOHR);
      auto H27 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.744503 * ANGSTROM_TO_BOHR,
                                        -0.680551 * ANGSTROM_TO_BOHR, 502.154266 * ANGSTROM_TO_BOHR);
      auto H28 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.062174 * ANGSTROM_TO_BOHR,
                                        0.969801 * ANGSTROM_TO_BOHR, 501.954130 * ANGSTROM_TO_BOHR);
      auto H20 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 3.711687 * ANGSTROM_TO_BOHR,
                                        0.624472 * ANGSTROM_TO_BOHR, 499.376826 * ANGSTROM_TO_BOHR);
      auto H30 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.777291 * ANGSTROM_TO_BOHR,
                                        0.707250 * ANGSTROM_TO_BOHR, 497.872814 * ANGSTROM_TO_BOHR);
      auto H31 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.557783 * ANGSTROM_TO_BOHR,
                                        1.953004 * ANGSTROM_TO_BOHR, 499.126303 * ANGSTROM_TO_BOHR);
      auto H32 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.526934 * ANGSTROM_TO_BOHR,
                                        -4.695631 * ANGSTROM_TO_BOHR, 503.345795 * ANGSTROM_TO_BOHR);
      auto Cl33 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Cl"), -0.775162 * ANGSTROM_TO_BOHR,
                                         -1.890098 * ANGSTROM_TO_BOHR, 503.327024 * ANGSTROM_TO_BOHR);
      auto H34 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.297193 * ANGSTROM_TO_BOHR,
                                        -6.331553 * ANGSTROM_TO_BOHR, 501.498947 * ANGSTROM_TO_BOHR);
      auto Cl35 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Cl"), -0.430177 * ANGSTROM_TO_BOHR,
                                         -2.934197 * ANGSTROM_TO_BOHR, 498.061200 * ANGSTROM_TO_BOHR);
      auto Cl36 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Cl"), -0.063741 * ANGSTROM_TO_BOHR,
                                         1.286285 * ANGSTROM_TO_BOHR, 496.710009 * ANGSTROM_TO_BOHR);
      auto H37 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.612114 * ANGSTROM_TO_BOHR,
                                        3.963798 * ANGSTROM_TO_BOHR, 496.059637 * ANGSTROM_TO_BOHR);
      auto H38 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.095277 * ANGSTROM_TO_BOHR,
                                        1.955427 * ANGSTROM_TO_BOHR, 501.038381 * ANGSTROM_TO_BOHR);
      auto H39 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.505883 * ANGSTROM_TO_BOHR,
                                        2.799836 * ANGSTROM_TO_BOHR, 499.595904 * ANGSTROM_TO_BOHR);
      auto H40 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.534387 * ANGSTROM_TO_BOHR,
                                        1.363453 * ANGSTROM_TO_BOHR, 499.440989 * ANGSTROM_TO_BOHR);
      auto H41 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.976349 * ANGSTROM_TO_BOHR,
                                        -2.032586 * ANGSTROM_TO_BOHR, 501.115232 * ANGSTROM_TO_BOHR);
      auto H42 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.276780 * ANGSTROM_TO_BOHR,
                                        -0.451956 * ANGSTROM_TO_BOHR, 501.854963 * ANGSTROM_TO_BOHR);
      auto H43 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -3.756936 * ANGSTROM_TO_BOHR,
                                        -0.737796 * ANGSTROM_TO_BOHR, 500.186532 * ANGSTROM_TO_BOHR);
      auto Cl44 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("Cl"), -0.213258 * ANGSTROM_TO_BOHR,
                                         3.588070 * ANGSTROM_TO_BOHR, 501.566857 * ANGSTROM_TO_BOHR);
      auto H45 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.500297 * ANGSTROM_TO_BOHR,
                                        5.784163 * ANGSTROM_TO_BOHR, 499.949971 * ANGSTROM_TO_BOHR);
      auto H46 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 3.482835 * ANGSTROM_TO_BOHR,
                                        -1.487614 * ANGSTROM_TO_BOHR, 500.156296 * ANGSTROM_TO_BOHR);
      auto H47 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.659532 * ANGSTROM_TO_BOHR,
                                        -2.089957 * ANGSTROM_TO_BOHR, 498.706117 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{
              Pt1,  N2,  N3,   C4,   C5,  H6,  C7,  H8,  H9,  C10, C11, C12,  C13, C14, C15, C16,
              C17,  C18, C19,  C20,  C21, C22, C23, C24, C25, H26, H27, H28,  H20, H30, H31, H32,
              Cl33, H34, Cl35, Cl36, H37, H38, H39, H40, H41, H42, H43, Cl44, H45, H46, H47});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), ANGSTROM_TO_BOHR * 1.606795380,
                                       ANGSTROM_TO_BOHR * -0.066291721, ANGSTROM_TO_BOHR * 0.087598426);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * 0.639117440,
                                       ANGSTROM_TO_BOHR * 0.066884843, ANGSTROM_TO_BOHR * 0.021624441);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * 1.981743900,
                                       ANGSTROM_TO_BOHR * 0.875178014, ANGSTROM_TO_BOHR * 0.036509196);
      auto O2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), ANGSTROM_TO_BOHR * -1.183706762,
                                       ANGSTROM_TO_BOHR * 0.042210901, ANGSTROM_TO_BOHR * -0.071849209);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * -1.486650370,
                                       ANGSTROM_TO_BOHR * -0.493436911, ANGSTROM_TO_BOHR * -0.764205635);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), ANGSTROM_TO_BOHR * -1.557299588,
                                       ANGSTROM_TO_BOHR * -0.318709677, ANGSTROM_TO_BOHR * 0.796158231);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H1, H2, O2, H3, H4});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Def2_TZVP_DFT: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 0.8948052 * ANGSTROM_TO_BOHR,
                                       -0.4686013 * ANGSTROM_TO_BOHR, 0.0 * ANGSTROM_TO_BOHR);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.0257715 * ANGSTROM_TO_BOHR,
                                       -0.1879168 * ANGSTROM_TO_BOHR, 0.0 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.3778994 * ANGSTROM_TO_BOHR,
                                       0.3639341 * ANGSTROM_TO_BOHR, 0.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::Water_Def2_TZVP_DFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H1, H2});
    } break;
    case TEST_SYSTEM_CONTROLLERS::He_Def2_TZVP_DFT: {
      auto He1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("He"), -2.2469331 * ANGSTROM_TO_BOHR,
                                        0.2925840 * ANGSTROM_TO_BOHR, 0.0 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::He_Def2_TZVP_DFT] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{He1});
    } break;
    case TEST_SYSTEM_CONTROLLERS::ETHANOL_def2_SVP_HF: {
      auto C1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 27.3900000000 * ANGSTROM_TO_BOHR,
                                       0.5560000000 * ANGSTROM_TO_BOHR, 29.6700000000 * ANGSTROM_TO_BOHR);
      auto C2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("C"), 27.8400000000 * ANGSTROM_TO_BOHR,
                                       1.9910000000 * ANGSTROM_TO_BOHR, 30.0160000000 * ANGSTROM_TO_BOHR);
      auto H1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 27.5720000000 * ANGSTROM_TO_BOHR,
                                       -0.1400000000 * ANGSTROM_TO_BOHR, 30.4990000000 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 26.3150000000 * ANGSTROM_TO_BOHR,
                                       0.5060000000 * ANGSTROM_TO_BOHR, 29.4570000000 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 27.8870000000 * ANGSTROM_TO_BOHR,
                                       0.2130000000 * ANGSTROM_TO_BOHR, 28.7540000000 * ANGSTROM_TO_BOHR);
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 29.2210000000 * ANGSTROM_TO_BOHR,
                                       1.9690000000 * ANGSTROM_TO_BOHR, 30.3810000000 * ANGSTROM_TO_BOHR);
      auto H4 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 27.6760000000 * ANGSTROM_TO_BOHR,
                                       2.6600000000 * ANGSTROM_TO_BOHR, 29.1620000000 * ANGSTROM_TO_BOHR);
      auto H5 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 27.2470000000 * ANGSTROM_TO_BOHR,
                                       2.3330000000 * ANGSTROM_TO_BOHR, 30.8730000000 * ANGSTROM_TO_BOHR);
      auto H6 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 29.7410000000 * ANGSTROM_TO_BOHR,
                                       1.8070000000 * ANGSTROM_TO_BOHR, 29.5750000000 * ANGSTROM_TO_BOHR);
      _testGeometries[TEST_SYSTEM_CONTROLLERS::ETHANOL_def2_SVP_HF] =
          std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{C1, C2, H1, H2, H3, O1, H4, H5, H6});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_A: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), -0.0805817544 * ANGSTROM_TO_BOHR,
                                       0.1064693332 * ANGSTROM_TO_BOHR, -0.2594723403 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.8470675945 * ANGSTROM_TO_BOHR,
                                       -0.3167490363 * ANGSTROM_TO_BOHR, 0.1589804590 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.0989011452 * ANGSTROM_TO_BOHR,
                                       0.9050485492 * ANGSTROM_TO_BOHR, 0.2294897437 * ANGSTROM_TO_BOHR);
      _testGeometries[kind] = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H2, H3});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_B: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 1.9791350365 * ANGSTROM_TO_BOHR,
                                       -1.5279359818 * ANGSTROM_TO_BOHR, -0.9751702547 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.3919346333 * ANGSTROM_TO_BOHR,
                                       -0.9685682058 * ANGSTROM_TO_BOHR, -0.4298323393 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.5522203445 * ANGSTROM_TO_BOHR,
                                       -0.9256451130 * ANGSTROM_TO_BOHR, -1.4651178122 * ANGSTROM_TO_BOHR);
      _testGeometries[kind] = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H2, H3});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_C: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), -0.0301359408 * ANGSTROM_TO_BOHR,
                                       -3.2548441887 * ANGSTROM_TO_BOHR, -1.6486920118 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -0.7035520673 * ANGSTROM_TO_BOHR,
                                       -2.9378089905 * ANGSTROM_TO_BOHR, -1.0363643169 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 0.6157692075 * ANGSTROM_TO_BOHR,
                                       -2.5322492123 * ANGSTROM_TO_BOHR, -1.7363814116 * ANGSTROM_TO_BOHR);
      _testGeometries[kind] = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H2, H3});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_D: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), -1.9208998680 * ANGSTROM_TO_BOHR,
                                       -1.3940520287 * ANGSTROM_TO_BOHR, -1.0889282227 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -2.7490198612 * ANGSTROM_TO_BOHR,
                                       -1.7880773544 * ANGSTROM_TO_BOHR, -1.3621145487 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), -1.5633729696 * ANGSTROM_TO_BOHR,
                                       -0.9304729700 * ANGSTROM_TO_BOHR, -1.8521807194 * ANGSTROM_TO_BOHR);
      _testGeometries[kind] = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H2, H3});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_E: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 2.4306445122 * ANGSTROM_TO_BOHR,
                                       -4.1413455009 * ANGSTROM_TO_BOHR, -0.4261322320 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 1.5151513815 * ANGSTROM_TO_BOHR,
                                       -4.2161335945 * ANGSTROM_TO_BOHR, -0.7462365031 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 2.8192429543 * ANGSTROM_TO_BOHR,
                                       -3.4299900532 * ANGSTROM_TO_BOHR, -0.9407634139 * ANGSTROM_TO_BOHR);
      _testGeometries[kind] = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H2, H3});
    } break;
    case TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_F: {
      auto O1 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("O"), 4.3837327957 * ANGSTROM_TO_BOHR,
                                       -1.5445500612 * ANGSTROM_TO_BOHR, -0.4985653758 * ANGSTROM_TO_BOHR);
      auto H2 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 3.7054352760 * ANGSTROM_TO_BOHR,
                                       -2.0726823807 * ANGSTROM_TO_BOHR, -0.0492937751 * ANGSTROM_TO_BOHR);
      auto H3 = std::make_shared<Atom>(AtomTypeFactory::getAtomType("H"), 5.2253069878 * ANGSTROM_TO_BOHR,
                                       -1.8059456348 * ANGSTROM_TO_BOHR, -0.1414732039 * ANGSTROM_TO_BOHR);
      _testGeometries[kind] = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>{O1, H2, H3});
    } break;
  }
}
void SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS kind) {
  if (_testSystemControllers.find(kind) == _testSystemControllers.end())
    return;
  auto system = _testSystemControllers[kind];
  cleanUpSystemDirectory(system);
  _testSystemControllers.erase(kind);
  _testGeometries.erase(kind);
}

void SystemController__TEST_SUPPLY::cleanUp() {
  std::remove("PairIntegrals.h5");
  for (auto system : _testSystemControllers) {
    cleanUpSystemDirectory(system.second);
  }
  _testSystemControllers.clear();
  _testGeometries.clear();
}

void SystemController__TEST_SUPPLY::cleanUpSystemDirectory(std::shared_ptr<SystemController> systemController) {
  cleanUpSystemDirectory(systemController->getSystemPath(), systemController->getSystemName(),
                         systemController->getSettings().basis.label);
}

void SystemController__TEST_SUPPLY::cleanUpSystemDirectory(std::string path, std::string systemName, std::string basisLabel) {
  removeSystemFiles(path, systemName, basisLabel);
}

} /* namespace Serenity */
