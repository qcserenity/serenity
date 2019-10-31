/**
 * @file   Settings.h
 *
 * @date   Mar 7, 2014
 * @author Thomas Dresselhaus
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
#ifndef SETTINGS_H_
#define SETTINGS_H_
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
#include "io/Filesystem.h"
#include "settings/Options.h"
#include "settings/Reflection.h"
/* Include Std and External Headers */
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>

namespace Serenity {
using namespace Serenity::Reflection;
/**
 * @class  Settings Settings.h
 * @brief  Holds one set of all switchable options for a single system.
 */
struct DFT
{
    DFT():
        functional(Options::XCFUNCTIONALS::BP86),
        densityFitting(Options::DENS_FITS::RI),
        dispersion(Options::DFT_DISPERSION_CORRECTIONS::NONE)
    {
    }
public:
    REFLECTABLE
    (
        (Options::XCFUNCTIONALS) functional,
        (Options::DENS_FITS) densityFitting,
        (Options::DFT_DISPERSION_CORRECTIONS) dispersion
    )
};

struct SCF
{
    SCF():
        initialguess(Options::INITIAL_GUESSES::ATOM_SCF),
        maxCycles(100),
        writeRestart(5),
        energyThreshold(1e-8),
        rmsdThreshold(1e-8),
        damping(Options::DAMPING_ALGORITHMS::SERIES),
        seriesDampingStart(0.7),
        seriesDampingEnd(0.6),
        seriesDampingStep(0.05),
        seriesDampingInitialSteps(2),
        staticDampingFactor(0.7),
        endDampErr(0.01),
        useLevelshift(true),
        useOffDiagLevelshift(false),
        diisFlush(1000),
        diisStartError(0.1),
        diisMaxStore(5),
        diisThreshold(1e-7),
        diisConditionNumberThreshold(1000.0),
        canOrthThreshold(1.0e-7),
        useADIIS(false)
    {
    }
public:
    REFLECTABLE
    (
        (Options::INITIAL_GUESSES) initialguess,
        (unsigned int) maxCycles,
        (unsigned int) writeRestart,
        (double) energyThreshold,
        (double) rmsdThreshold,
        (Options::DAMPING_ALGORITHMS) damping,
        (double) seriesDampingStart,
        (double) seriesDampingEnd,
        (double) seriesDampingStep,
        (double) seriesDampingInitialSteps,
        (double) staticDampingFactor,
        (double) endDampErr,
        (bool) useLevelshift,
        (bool) useOffDiagLevelshift,
        (unsigned int) diisFlush,
        (double) diisStartError,
        (unsigned int) diisMaxStore,
        (double) diisThreshold,
        (double) diisConditionNumberThreshold,
        (double) canOrthThreshold,
        (bool)   useADIIS
    )
};
struct BASIS
{
    BASIS():
     label("6-31GS"),
     auxJLabel("RI_J_WEIGEND"),
     auxCLabel(""),
     makeSphericalBasis(true),
     integralThreshold(1e-10),
     basisLibPath(""),
     firstECP(37)
    {
    }
public:
    REFLECTABLE
    (
     (std::string) label,
     (std::string) auxJLabel,
     (std::string) auxCLabel,
     (bool) makeSphericalBasis,
     (double) integralThreshold,
     (std::string) basisLibPath,
     (unsigned int) firstECP
    )
};

struct GRID
{
    GRID():
     gridType(Options::GRID_TYPES::SSF),
     radialGridType(Options::RADIAL_GRID_TYPES::AHLRICHS),
     sphericalGridType(Options::SPHERICAL_GRID_TYPES::LEBEDEV),
     blocksize(128),
     accuracy(4),
     smallGridAccuracy(2),
     blockAveThreshold(1e-11),
     basFuncRadialThreshold(1e-9),
     weightThreshold(1e-14),
     smoothing(3)
    {
    }
public:
    REFLECTABLE
    (
     (Options::GRID_TYPES) gridType,
     (Options::RADIAL_GRID_TYPES) radialGridType,
     (Options::SPHERICAL_GRID_TYPES) sphericalGridType,
     (unsigned int) blocksize,
     (unsigned int) accuracy,
     (unsigned int) smallGridAccuracy,
     (double) blockAveThreshold,
     (double) basFuncRadialThreshold,
     (double) weightThreshold,
     (unsigned int) smoothing
    )
};
struct Settings {
public:
  /**
   * @brief Default constructor
   *
   * This constructor holds all the default values for settings
   *   that are not part of any other block.
   * Each block has its defaults set in its own constructor.
   *
   * Add new (non-blocked) variables in this constructor.
   */
  Settings():
      name("default"),
      identifier(""),
      path("./"),
      charge(0),
      spin(0),
      geometry(""),
      load(""),
      scfMode(Options::SCF_MODES::RESTRICTED),
      method(Options::ELECTRONIC_STRUCTURE_THEORIES::HF)
  {
    // Set identifier
    identifier = std::to_string(std::chrono::duration_cast< std::chrono::nanoseconds >(
        std::chrono::system_clock::now().time_since_epoch()).count());
  }
  /*
   * Member Variables
   */
  REFLECTABLE(
    (std::string) name,
    (std::string) identifier,
    (std::string) path,
    (int) charge,
    (int) spin,
    (std::string) geometry,
    (std::string) load,
    (Options::SCF_MODES) scfMode,
    (Options::ELECTRONIC_STRUCTURE_THEORIES) method,
    (DFT) dft,
    (SCF) scf,
    (BASIS) basis,
    (GRID) grid
  )

  /**
   * @brief Constructor using text input.
   */
  Settings(std::ifstream& input) : Settings(){
    std::string line;
    std::string word;
    while(getline(input, line)){
      std::istringstream iss(line);
      word ="";
      iss >> word;
      // check for comments
      if (word[0] == '#') continue;
      if (word.empty()) continue;
      // blocks
      if (word[0] == '+'){
        std::string blockname = word.erase(0, 1);
        for (auto& c: blockname) c = std::toupper(c);
        if(blockname=="SYSTEM") continue;
        while(getline(input, line)){
          std::istringstream iss2(line);
          word ="";
          iss2 >> word;
          if (word[0] == '-') break;
          // check for comments
          if (word[0] == '#') continue;
          if (word.empty()) continue;
          std::string name = word;
          if (!(iss2 >> word)){
            throw SerenityError("ERROR: Value missing for keyword '"+name+"'.");
          }
          std::string value = word;
          this->set(blockname,name,value);
        }
        continue;
      }
      // end of system block
      if (!word.compare("-system")){
        /*
         * Capitalize basis label strings
         */
        for (auto& c: this->basis.label) c = std::toupper(c);
        for (auto& c: this->basis.auxCLabel) c = std::toupper(c);
        for (auto& c: this->basis.auxJLabel) c = std::toupper(c);
        // Set path string to end on "/"
        if(path.substr(path.length()-1)!="/") path=path+"/";
        if(spin!=0) scfMode=Options::SCF_MODES::UNRESTRICTED;
        return;
      }
      // unblocked options
      std::string name = word;
      if (!(iss >> word)){
        throw SerenityError("ERROR: Value missing for keyword '"+name+"'.");
      }
      std::string value = word;
      this->set("",name,value);
    }
    /*
     * Cannot be reached
     */
    throw SerenityError("Error while parsing settings file.");
  }

  /**
   * @brief Setter for options.
   *
   * @param blockname
   * @param name
   * @param value
   */
  void set(std::string blockname, std::string name, std::string value){
    /*
     * Definitions
     */
    bool check = false;

// TODO this templated lambda function version
//       of the visitor should work with C++14
//
//    auto visitor = [&check,&name,&value](auto f){
//      if(name.compare(f.name()) == 0 ){
//        std::cout << f.name() << " "<< typeid(f.get()).name() << std::endl;
//        check = true;
//      }
//    };
    /*
     * Parsing
     * (Add new blocks here aswell)
     */
    if (!blockname.compare("")){
      set_visitor visitor(name,value,check);
      visit_each(*this, visitor);
    } else if (!blockname.compare("DFT")){
      set_visitor visitor(name,value,check);
      visit_each((*this).dft, visitor);
    } else if (!blockname.compare("SCF")){
      set_visitor visitor(name,value,check);
      visit_each((*this).scf, visitor);
    } else if (!blockname.compare("BASIS")){
      set_visitor visitor(name,value,check);
      visit_each((*this).basis, visitor);
    } else if (!blockname.compare("GRID")){
      set_visitor visitor(name,value,check);
      visit_each((*this).grid, visitor);
    } else {
      throw SerenityError("ERROR: No block '"+blockname+"' known.");
    }
    if (!check){
      throw SerenityError("ERROR: No keyword '"+name+"' known in this block.");
    }
  }



  void printSettings(){
    std::string field;
    std::string value;
    // Create folder if it does not exist
    makePath(path);

    std::ofstream ofs;
    ofs.open ((*this).path+(*this).name+".settings", std::ofstream::out | std::ofstream::trunc);
    ofs << "#=======================================================" << std::endl;
    ofs << "# NOTE:" << std::endl;
    ofs << "# This file contains the entire list of settings" << std::endl;
    ofs << "#  that were stored in one instance of the Settings" << std::endl;
    ofs << "#  class." << std::endl;
    ofs << "# The settings given here have not necessarily been" << std::endl;
    ofs << "#  used, e.g. a functional might not have been needed." << std::endl;
    ofs << "#  in a HF calculations." << std::endl;
    ofs << "#=======================================================" << std::endl;
    ofs << "+system" << std::endl;
    print_visitor visitor(field,value,ofs);
    visit_each((*this),visitor);
    ofs << "+dft" << std::endl;
    visit_each((*this).dft,visitor);
    ofs << "-dft" << std::endl;
    ofs << "+scf" << std::endl;
    visit_each((*this).scf,visitor);
    ofs << "-scf" << std::endl;
    ofs << "+basis" << std::endl;
    visit_each((*this).basis,visitor);
    ofs << "-basis" << std::endl;
    ofs << "+grid" << std::endl;
    visit_each((*this).grid,visitor);
    ofs << "-grid" << std::endl;
    ofs << "-system" << std::endl;
    ofs.close();
  }



};
} /* namespace Serenity */

#endif /* SETTINGS_H_ */
