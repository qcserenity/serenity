/**
 * @file   SystemController__TEST_SUPPLY.h
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   28. August 2015, 14:19
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
#ifndef SYSTEMCONTROLLER__TEST_SUPPLY_H
#define	SYSTEMCONTROLLER__TEST_SUPPLY_H
/* Include Serenity Internal Headers */
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Serenity {
/* Forward declarations */
class Geometry;
struct Settings;
/**
 * These kinds of test systems are available:\n
 * 
 * H2_MINBAS:\n
 * 
 * An H2 molecule in a minimal basis at the minimum geometry.\n
 * 
 * H2_MINBAS_LARGE_DISTANCE:\n
 * 
 * The same as above, but at a very large distance. The correct energy in that case is exactly -1.0
 * Hartree, but can not be obtained with all methods. It is good to test UHF, because that method
 * is capable of obtaining the correct result, in contrast to e.g. RHF.
 * 
 * WATER_DISTORTED_MINBAS:\n
 * 
 * A water molecule with unsymetrically distorted geometry, i.e. not in equilibrium/ground state
 */
enum class TEST_SYSTEM_CONTROLLERS { H2_MINBAS=0, H2_MINBAS_ACTIVE, H2_MINBAS_ENVIRONMENT, H2_MINBAS_LARGE_DISTANCE, H2_DEF2_TZVP, WATER_DISTORTED_MINBAS, CO_MINBAS, JACOBSEN_MINBAS,
                                     WaterMonOne_6_31Gs, WaterMonTwo_6_31Gs, Ne2_6_31Gs, Ar2_6_31Gs, Kr2_6_31Gs, F2_6_31Gs, Cl2_6_31Gs, Br2_6_31Gs,H2_6_31Gs_BP86,H2_6_31Gs_ACTIVE_FDE,
                                     H2_6_31Gs_ENVIRONMENT_FDE,O2_MINBAS_SING,O2_MINBAS_TRIP,F_MINUS_6_31Gs, H2_DEF2_TZVP_PBE,H2_DEF2_SV_P_PBE,C60_MINBAS,WATER_DEF2_SVP_CAMB3LYP,
                                     WaterMonTwo_6_31Gs_DFT, WaterMonOne_6_31Gs_DFT,H60_Ghost_MINBAS,OH_MINBAS_PBE,MethylRad_def2_SVP_PBE,HI_Def2_SVP_PBE,I2_Def2_SVP_PBE,H2_MINBAS_LINEARDEPENDENT,
                                     H2_6_31Gs_ACTIVE_FDE_HYBRID,H2_6_31Gs_ENVIRONMENT_FDE_HYBRID,H2_DEF2_TZVP_CAMB3LYP,H2_DEF2_TZVP_PBE_NORI,H2_DEF2_TZVP_HF,H2_DEF2_TZVP_HF_UNRESTRICTED,
                                     H2_6_31Gs_ACTIVE_FDE_BP86,H2_6_31Gs_ENVIRONMENT_FDE_BP86,He2_6_31Gs_BP86,H2_6_31Gs_ACTIVE_FDE_B3LYP,H2_6_31Gs_ENVIRONMENT_FDE_B3LYP,H2_def2_SVP_ACTIVE_FDE,
                                     H2_def2_SVP_ACTIVE_FDE_BP86, H2_def2_SVP_ENVIRONMENT_FDE_BP86,H2_def2_SVP_ACTIVE_FDE_PBE0_UNRESTRICTED,H2_def2_SVP_ENVIRONMENT_FDE_PBE0_UNRESTRICTED,
                                     PHENOLATE_O_DEF2_SVP_BP86, PHENOLATE_PHENYL_DEF2_SVP_BP86};

/**
 * @class SystemController__TEST_SUPPLY SystemController__TEST_SUPPLY.h
 * @brief Provides SystemControllers ready to use for (internal) functionality tests.
 * 
 * All functionalities of the program should be available through them.
 */
class SystemController__TEST_SUPPLY {
  SystemController__TEST_SUPPLY() = delete;
public:
  static std::shared_ptr<SystemController> getSystemController(TEST_SYSTEM_CONTROLLERS kind,
      bool fromScratch=false) {
    if (!_testSystemControllers[kind] or fromScratch) prepare(kind,fromScratch);
    return _testSystemControllers[kind];
  }
  static std::shared_ptr<SystemController> getSystemController(
      TEST_SYSTEM_CONTROLLERS kind, Settings& settings, int charge=0, int spin=0);
  static std::shared_ptr<Geometry> getGeometry(TEST_SYSTEM_CONTROLLERS kind) {
    if (!_testGeometries[kind]) prepareGeometry(kind);
    return _testGeometries[kind];
  }
  /**
   * @brief This routine needs to be called if you left the system in a changed state!
   */
  static void forget(TEST_SYSTEM_CONTROLLERS kind) {
    if(_testSystemControllers.find(kind)==_testSystemControllers.end()) return;
    auto system=_testSystemControllers[kind];
    /* TODO
     * Replace with std::filesystem when c++17 is available
     */
    std::remove((system->getSettings().path+system->getSettings().name+".settings").c_str());
    std::remove((system->getSettings().path+system->getSettings().name+".xyz").c_str());
    std::remove((system->getSettings().path+system->getSettings().name+".energies.res.h5").c_str());
    std::remove((system->getSettings().path+system->getSettings().name+".energies.unres.h5").c_str());
    std::remove((system->getSettings().path+system->getSettings().name+".orbs.res.h5").c_str());
    std::remove((system->getSettings().path+system->getSettings().name+".orbs.unres.h5").c_str());
    std::remove((system->getSettings().path+system->getSettings().name+".dmat.res.h5").c_str());
    std::remove((system->getSettings().path+system->getSettings().name+".dmat.unres.h5").c_str());
    std::remove((system->getSettings().path+system->getSettings().name+".hess.h5").c_str());
    std::remove((system->getSettings().path+system->getSettings().name).c_str());
    std::remove((system->getSettings().path).c_str());
    std::remove("WARNING");
    _testSystemControllers.erase(kind);
    _testGeometries.erase(kind);
  }
  static void cleanUp() {
    for(auto system : _testSystemControllers){
      /* TODO
       * Replace with std::filesystem when c++17 is available
       */
      std::string path = system.second->getSettings().path;
      std::string name = system.second->getSettings().name;
      std::remove((path+name+".settings").c_str());
      std::remove((path+name+".xyz").c_str());
      std::remove((path+name+".energies.res.h5").c_str());
      std::remove((path+name+".energies.unres.h5").c_str());
      std::remove((path+name+".orbs.res.h5").c_str());
      std::remove((path+name+".orbs.unres.h5").c_str());
      std::remove((path+name+".dmat.res.h5").c_str());
      std::remove((path+name+".dmat.unres.h5").c_str());
      std::remove((path+name+".basis.h5").c_str());
      std::remove((path+name+".hess.h5").c_str());
      std::remove((path+name).c_str());
      std::remove((path).c_str());
      std::remove("WARNING");
    }
    _testSystemControllers.clear();
    _testGeometries.clear();
  }
private:
  static void prepare(TEST_SYSTEM_CONTROLLERS kind,bool fromScratch);
  static std::map<TEST_SYSTEM_CONTROLLERS, std::shared_ptr<SystemController> > _testSystemControllers;
  static void prepareGeometry(TEST_SYSTEM_CONTROLLERS kind);
  static std::map<TEST_SYSTEM_CONTROLLERS, std::shared_ptr<Geometry> > _testGeometries;
};

} /* namespace Serenity */
#endif	/* SYSTEMCONTROLLER__TEST_SUPPLY_H */
