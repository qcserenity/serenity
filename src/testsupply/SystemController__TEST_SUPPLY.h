/**
 * @file   SystemController__TEST_SUPPLY.h
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
#ifndef SYSTEMCONTROLLER__TEST_SUPPLY_H
#define SYSTEMCONTROLLER__TEST_SUPPLY_H
/* Include Std and External Headers */
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Serenity {
/* Forward declarations */
class Geometry;
struct Settings;
class SystemController;
/**
 * These kinds of test systems are available:
 */
enum class TEST_SYSTEM_CONTROLLERS {
  H2_MINBAS = 0,
  H2_MINBAS_ACTIVE,
  H2_MINBAS_ACTIVE_LRSCF,
  H2_MINBAS_ENVIRONMENT,
  H2_MINBAS_ENVIRONMENT_LRSCF,
  H2_MINBAS_LARGE_DISTANCE,
  H2_DEF2_TZVP,
  WATER_DISTORTED_MINBAS,
  CO_MINBAS,
  JACOBSEN_MINBAS,
  WaterMonOne_6_31Gs,
  WaterMonTwo_6_31Gs,
  WaterMon1_3_DFT,
  WaterMon2_3_DFT,
  WaterMon3_3_DFT,
  WaterDim_DFT,
  Ne2_6_31Gs,
  Ar2_6_31Gs,
  Kr2_6_31Gs,
  F2_6_31Gs,
  Cl2_6_31Gs,
  Br2_6_31Gs,
  H2_6_31Gs_BP86,
  H2_6_31Gs_ACTIVE_FDE,
  H2_6_31Gs_ENVIRONMENT_FDE,
  O2_MINBAS_SING,
  O2_MINBAS_TRIP,
  O2_MINBAS_TRIP_CIS,
  O2_TRIP_DEF2_SVP_TDA,
  O2_TRIP_DEF2_SVP_TDHF,
  O2_TRIP_DEF2_SVP_TDDFT,
  O2_TRIP_6_31G_PBE_TDA,
  O2_TRIP_6_31G_PBE_TDDFT,
  F_MINUS_6_31Gs,
  H2_DEF2_TZVP_PBE,
  H2_DEF2_SV_P_PBE,
  C60_MINBAS,
  WATER_DEF2_SVP_CAMB3LYP,
  WaterMonTwo_6_31Gs_DFT,
  WaterMonOne_6_31Gs_DFT,
  H60_Ghost_MINBAS,
  OH_MINBAS_PBE,
  MethylRad_def2_SVP_PBE,
  HI_Def2_SVP_PBE,
  I2_Def2_SVP_PBE,
  H2_MINBAS_LINEARDEPENDENT,
  EthaneA_Def2_SVP_BP86_Act,
  EthaneA_Def2_SVP_BP86_Env,
  EthaneA_Def2_SVP_BP86,
  EthaneB_Def2_SVP_BP86,
  Formaldehyde_HF_AUG_CC_PVDZ,
  Diaziridine_HF_AUG_CC_PVDZ,
  WaterMonTwo_Def2_SVP,
  WaterMonOne_Def2_SVP,
  WaterMonOne_Def2_SVP_B2PLYP,
  EthaneA_Def2_SVP,
  H2Dimer_Def2_TZVP,
  OHPhenol_Def2_SVP_Act,
  PhenylPhenol_Def2_SVP_Env,
  H2_6_31Gs_ACTIVE_FDE_HYBRID,
  H2_6_31Gs_ENVIRONMENT_FDE_HYBRID,
  H2_DEF2_TZVP_CAMB3LYP,
  H2_DEF2_TZVP_PBE_NORI,
  H2_DEF2_TZVP_HF,
  H2_DEF2_TZVP_HF_UNRESTRICTED,
  H2_6_31Gs_ACTIVE_FDE_BP86,
  H2_6_31Gs_ACTIVE_FDE_BP86_Supermolecular,
  H2_6_31Gs_ENVIRONMENT_FDE_BP86,
  He2_6_31Gs_BP86,
  H2_6_31Gs_ACTIVE_FDE_B3LYP,
  H2_6_31Gs_ENVIRONMENT_FDE_B3LYP,
  H2_def2_SVP_ACTIVE_FDE,
  H2_def2_SVP_ACTIVE_FDE_BP86,
  H2_def2_SVP_ENVIRONMENT_FDE_BP86,
  H2_def2_SVP_ACTIVE_FDE_PBE0_UNRESTRICTED,
  H2_def2_SVP_ENVIRONMENT_FDE_PBE0_UNRESTRICTED,
  PHENOLATE_O_DEF2_SVP_BP86,
  PHENOLATE_PHENYL_DEF2_SVP_BP86,
  H2_DEF2_TZVP_PW91_UNRES_0,
  H2_DEF2_TZVP_PW91_UNRES_1,
  BE_DEF2_TZVP_PW91_UNRES_0,
  BE_DEF2_TZVP_PW91_UNRES_1,
  He_1_6_31Gs_BP86,
  He_2_6_31Gs_BP86,
  He_3_6_31Gs_BP86,
  He_1_def2SVP_BP86,
  He_2_def2SVP_BP86,
  He_3_def2SVP_BP86,
  He2_def2SVP_BP86,
  HeNe_def2SVP_BP86,
  H2_6_311G_3A,
  H2_6_311G_BS_3A,
  H2_6_311G_Act_3A,
  H2_6_311G_Env_3A,
  H2_6_311G_BSAct_3A,
  H2_6_311G_BSEnv_3A,
  MethylRad_Act_def2_SVP_PBE,
  MethylRad_Env_def2_SVP_PBE,
  WCCR1010_def2_SVP_HF,
  Water_Dimer_def2_SVP_HF,
  Water_Def2_TZVP_DFT,
  He_Def2_TZVP_DFT,
  ETHANOL_def2_SVP_HF,
  Water_Hexamer_Monomer_A,
  Water_Hexamer_Monomer_B,
  Water_Hexamer_Monomer_C,
  Water_Hexamer_Monomer_D,
  Water_Hexamer_Monomer_E,
  Water_Hexamer_Monomer_F
};

/**
 * @class SystemController__TEST_SUPPLY SystemController__TEST_SUPPLY.h
 * @brief Provides SystemControllers ready to use for (internal) functionality tests.
 */
class SystemController__TEST_SUPPLY {
  SystemController__TEST_SUPPLY() = delete;

 public:
  /**
   * @brief Getter for the various test system controllers.
   * @param kind The enum flag defining the system controller.
   * @param fromScratch If true, the system is built from scratch (using the settings encoded in this .cp). Otherwise
   * the system will be loaded from disk. Note that there are many cases where the two modes don't produce the same
   * settings!
   * @return The system controller.
   */
  static std::shared_ptr<SystemController> getSystemController(TEST_SYSTEM_CONTROLLERS kind, bool fromScratch = false) {
    if (!_testSystemControllers[kind] or fromScratch)
      prepare(kind, fromScratch);
    return _testSystemControllers[kind];
  }
  /**
   * @brief Getter for a system controller with non default charge, spin and settings.
   * @param kind The enum flag.
   * @param settings The settings.
   * @param charge The charge.
   * @param spin The spin.
   * @return The system controller.
   */
  static std::shared_ptr<SystemController> getSystemController(TEST_SYSTEM_CONTROLLERS kind, Settings& settings,
                                                               int charge = 0, int spin = 0);
  /**
   * @brief Getter for the geometry associated to test system.
   * @param kind The enum flag.
   * @return The geometry.
   */
  static std::shared_ptr<Geometry> getGeometry(TEST_SYSTEM_CONTROLLERS kind) {
    if (!_testGeometries[kind])
      prepareGeometry(kind);
    return _testGeometries[kind];
  }
  /**
   * @brief This routine needs to be called if you left the system in a changed state!
   * @param The enum flag.
   */
  static void forget(TEST_SYSTEM_CONTROLLERS kind);
  /**
   * @brief Removes all files associated to the test systems.
   */
  static void cleanUp();
  /**
   * @brief Remove all system files at the given path and system name.
   * @param path          The path.
   * @param systemName    The system name-
   */
  static void cleanUpSystemDirectory(std::string path, std::string systemName, std::string basisLabel = "");
  /**
   * @brief Remove all system files for the given system controller.
   * @param systemController    The system controller.
   */
  static void cleanUpSystemDirectory(std::shared_ptr<SystemController> systemController);

 private:
  // Sets up the SystemControllers and can be called in two different modes: with fromScratch==true where the settings
  // defined in this .cpp are taken, and with fromScratch==false, where the settings and possibly some system results
  // like Fock- and Density matrix are read in from disk..
  static void prepare(TEST_SYSTEM_CONTROLLERS kind, bool fromScratch);
  // Map that stores all system controllers.
  static std::map<TEST_SYSTEM_CONTROLLERS, std::shared_ptr<SystemController>> _testSystemControllers;
  // Constructs the geometry for the given test system (with the coordinate information encoded in
  // SystemController__TEST_SUPPLY.cpp - the geometry is not read from file).
  static void prepareGeometry(TEST_SYSTEM_CONTROLLERS kind);
  // Map that stores the geometries of the test systems.
  static std::map<TEST_SYSTEM_CONTROLLERS, std::shared_ptr<Geometry>> _testGeometries;
};

} /* namespace Serenity */
#endif /* SYSTEMCONTROLLER__TEST_SUPPLY_H */
