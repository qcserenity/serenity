/**
 * @file   Options.h
 *
 * @date   Jan 23, 2017
 * @author David Schnieders
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

#ifndef OPTIONS_H_
#define OPTIONS_H_

/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>


namespace Serenity {
namespace Options {

template <class T> 
void resolve(std::string& value, T& field){
  (void)value;
  (void)field;
  throw SerenityError("ERROR: No mapping to string available for this type (Options::resolve)");
}
template<> inline void resolve<int>         (std::string& value, int& field){
  if(value.empty()){
    std::ostringstream strs;
    strs << field;
    value = strs.str();
  }else{
    try{
      field = std::stoi(value);
    } catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into an integer.");
    }
  }
}
template<> inline void resolve<unsigned int>(std::string& value, unsigned int& field){
  if(value.empty()){
    std::ostringstream strs;
    strs << field;
    value = strs.str();
  }else{
    try{
      field =  std::stoi(value);
    } catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into an unsigned integer.");
    }
  }
}
template<> inline void resolve<unsigned long int>(std::string& value, unsigned long int& field){
  if(value.empty()){
    std::ostringstream strs;
    strs << field;
    value = strs.str();
  }else{
    try{
      field =  std::stoul(value);
    } catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into an unsigned integer.");
    }
  }
}
template<> inline void resolve<double>      (std::string& value, double& field){
  if(value.empty()){
    std::ostringstream strs;
    strs << field;
    value = strs.str();
  }else{
    try{
      field = std::stod(value);
    } catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into a double.");
    }
  }
}
template<> inline void resolve<bool>      (std::string& value, bool& field){
  if(value.empty()){
    if(field){
      value="true";
    }else{
      value="false";
    }
  }else{
    std::string copy = value;
    for (auto& c: copy) c = std::toupper(c);
    if (!copy.compare("TRUE") or !copy.compare("1")){
      field = true;
    } else if(!copy.compare("FALSE") or !copy.compare("0")){
      field = false;
    } else {
      throw SerenityError("ERROR: Could not convert '"+value+"' into a boolean expression.");
    }
  }
}
template<> inline void resolve<std::string> (std::string& value, std::string& field){
  if(value.empty()){
    value=field;
  }else{
    field = value;
  }
}
template<> inline void resolve<std::vector<unsigned int>> (std::string& value, std::vector<unsigned int>& field){
  if(value.empty()){
    for(unsigned int val : field){
      std::ostringstream strs;
      strs << val << ' ';
      value=strs.str();
    }
  }else{
    try{
      field.clear();
      istringstream iss(value);
      std::string word;
      while(iss >> word){
        field.push_back(std::stoi(word));
      }
    }catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into a vector of unsigned integers.");
    }
  }
}
template<> inline void resolve<std::vector<unsigned long int>> (std::string& value, std::vector<unsigned long int>& field){
  if(value.empty()){
    for(unsigned long int val : field){
      std::ostringstream strs;
      strs << val << ' ';
      value=strs.str();
    }
  }else{
    try{
      field.clear();
      istringstream iss(value);
      std::string word;
      while(iss >> word){
        field.push_back(std::stoul(word));
      }
    }catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into a vector of unsigned long integers.");
    }
  }
}
template<> inline void resolve<std::vector<int>> (std::string& value, std::vector<int>& field){
  if(value.empty()){
    for(int val : field){
      std::ostringstream strs;
      strs << val << ' ';
      value=strs.str();
    }
  }else{
    try{
      field.clear();
      istringstream iss(value);
      std::string word;
      while(iss >> word){
        field.push_back(std::stoi(word));
      }
    }catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into a vector of integers.");
    }
  }
}
template<> inline void resolve<std::vector<bool>> (std::string& value, std::vector<bool>& field){
  if(value.empty()){
    for(bool val : field){
      std::ostringstream strs;
      if(val){
        strs << "true " << ' ';
      }else{
        strs << "false " << ' ';
      }
      value=strs.str();
    }
  }else{
    field.clear();
    istringstream iss(value);
    std::string word;
    while(iss >> word){
      std::string copy = word;
      for (auto& c: copy) c = std::toupper(c);
      if (!copy.compare("TRUE") or !copy.compare("1")){
        field.push_back(true);
      } else if(!copy.compare("FALSE") or !copy.compare("0")){
        field.push_back(false);
      } else {
        throw SerenityError("ERROR: Could not convert '"+value+"' into a vector of booleans.");
      }
    }
  }
}
template<> inline void resolve<std::vector<double>> (std::string& value, std::vector<double>& field){
  if(value.empty()){
    for(double val : field){
      std::ostringstream strs;
      strs << val << ' ';
      value=strs.str();
    }
  }else{
    try{
      field.clear();
      istringstream iss(value);
      std::string word;
      while(iss >> word){
        field.push_back(std::stod(word));
      }
    }catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into a vector of doubles.");
    }
  }
}
template<> inline void resolve<std::vector<std::string>> (std::string& value, std::vector<std::string>& field){
  if(value.empty()){
    for(std::string val : field){
      std::ostringstream strs;
      strs << val << ' ';
      value=strs.str();
    }
  }else{
    try{
      field.clear();
      istringstream iss(value);
      std::string word;
      while(iss >> word){
        field.push_back(word);
      }
    }catch(...){
      throw SerenityError("ERROR: Could not convert '"+value+"' into a vector of strings.");
    }
  }
}

template<class T> 
void check(const std::map<std::string, T> m,
    std::string& key,
    T& field){
  if(key.empty()){
    for(auto const& pair : m) {
      if(pair.second==field){
        key=pair.first;
        return;
      }
    }
    throw SerenityError("ERROR: Option map missing.");
  }else{
    try{
      for (auto& c: key) c = std::toupper(c);
      field = m.at(key);
    } catch(...) {
      std::string error = "ERROR: '"+key+"' is not a valid option.\n";
      error += "Valid options are:\n";
      for ( const auto &myPair : m ) {
        error += myPair.first+"\n";
      }
      throw SerenityError(error);
    }
  }
}

/**************************************************************************************************/
/*                                    Electronic Structure                                        */
/**************************************************************************************************/
/**
 * The type of SCF calculation that is made\n
 * RESTRICTED: all electrons are paired, works only for even numbers of electrons\n
 * UNRESTRICTED: uneven numbers of electrons are allowed
 */
enum class SCF_MODES {RESTRICTED=0,UNRESTRICTED=1};
template<> inline void resolve<SCF_MODES>(std::string& value, SCF_MODES& field){
  static const std::map<std::string, SCF_MODES> m = {
      {"RESTRICTED",SCF_MODES::RESTRICTED},
      {"UNRESTRICTED",SCF_MODES::UNRESTRICTED}
  };
  check(m,value,field);
}  

enum class ELECTRONIC_STRUCTURE_THEORIES {HF=0, DFT=1};
template<> inline void resolve<ELECTRONIC_STRUCTURE_THEORIES>(std::string& value, ELECTRONIC_STRUCTURE_THEORIES& field){
  static const std::map<std::string, ELECTRONIC_STRUCTURE_THEORIES> m = {
      {"HF",ELECTRONIC_STRUCTURE_THEORIES::HF},
      {"DFT",ELECTRONIC_STRUCTURE_THEORIES::DFT}
  };
  check(m,value,field);
}  

/**
 * Readily prepared functionals (are resolved to instances of Functional) in a way a user would
 * think of a functional. These are typically composed of different parts, see
 * FunctionalClassResolver (here you also find the BASIC_FUNCTIONALS).\n
 * ATTENTION: It is VERY important that the enum numbers here match the corresponding numbers
 * in the XCFUNCTIONALS and KINFUNCTIONALS enums.\n
 * 
 */
enum class FUNCTIONALS {
  // No functional to be used
  NONE=0,
      /*
       * XC Functionals
       */
      // LDA
      SLATER=1,
      VWN3=2,
      VWN5=3,
      LDA=4,
      LDAERF=5,
      // GGA
      B97=101,
      B97_1=102,
      B97_2=103,
      OLYP=104,
      BLYP=105,
      PBE=106,
      BP86=107,
      KT1=108,
      KT2=109,
      KT3=110,
      PW91=111,
      // Meta-GGA
      // Hybrid
      PBE0=301,
      B3LYP=302,
      B3LYP_G=303,
      B3P86=304,
      B3P86_G=305,
      BPW91=306,
      HF=307,
      // Meta-Hybrid
      // Range-Separated Hybrid
      CAMB3LYP=501,
      LCBLYP=502,
      LCBLYP_047=503,
      // Double-Hybrid
      B2PLYP=601,
      B2KPLYP=602,
      B2TPLYP=603,
      B2GPPLYP=604,
      ROB2PLYP=605,
      B2PIPLYP=606,
      //    MPW2PLYP=607,
      //    MPW2KPLYP=608,
      B2PPW91=609,
      //SCS Double-Hybrid
      DSDBLYP=610,
      DUT=611,
      PUT=612,
      DSDPBEP86=613,
      //SOS Double-Hybrid
      PWPB95=614,
      // Local-Hybrid
      // Model
      SAOP=801,

      /*
       * Kinetic Functionals
       */
      // LDA
      TF=1001,
      // GGA
      TW=1101,
      PW91K=1102,
      LLP91K=1103,
      LLP91KS=1104,
      PBE2K=1105,
      PBE2KS=1106,
      PBE3K=1107,
      PBE4K=1108,
      E2000K=1109,
};
template<> inline void resolve<FUNCTIONALS>(std::string& value, FUNCTIONALS& field){
  static const std::map<std::string, FUNCTIONALS> m = {
      {"NONE",FUNCTIONALS::NONE},
      {"SLATER",FUNCTIONALS::SLATER},
      {"VWN3",FUNCTIONALS::VWN3},
      {"VWN5",FUNCTIONALS::VWN5},
      {"LDAERF",FUNCTIONALS::LDAERF},
      {"LDA",FUNCTIONALS::LDA},
      {"B97",FUNCTIONALS::B97},
      {"B97_1",FUNCTIONALS::B97_1},
      {"B97_2",FUNCTIONALS::B97_2},
      {"OLYP",FUNCTIONALS::OLYP},
      {"BLYP",FUNCTIONALS::BLYP},
      {"PBE",FUNCTIONALS::PBE},
      {"BP86",FUNCTIONALS::BP86},
      {"KT1",FUNCTIONALS::KT1},
      {"KT2",FUNCTIONALS::KT2},
      {"KT3",FUNCTIONALS::KT3},
      {"PW91",FUNCTIONALS::PW91},
      {"PBE0",FUNCTIONALS::PBE0},
      {"B3LYP",FUNCTIONALS::B3LYP},
      {"B3LYP_G",FUNCTIONALS::B3LYP_G},
      {"B3P86",FUNCTIONALS::B3P86},
      {"B3P86_G",FUNCTIONALS::B3P86_G},
      {"BPW91",FUNCTIONALS::BPW91},
      {"CAMB3LYP",FUNCTIONALS::CAMB3LYP},
      {"LCBLYP",FUNCTIONALS::LCBLYP},
      {"LCBLYP_047",FUNCTIONALS::LCBLYP_047},
      {"B2PLYP",FUNCTIONALS::B2PLYP},
      {"B2KPLYP",FUNCTIONALS::B2KPLYP},
      {"B2TPLYP",FUNCTIONALS::B2TPLYP},
      {"B2GPPLYP",FUNCTIONALS::B2GPPLYP},
      {"ROB2PLYP",FUNCTIONALS::ROB2PLYP},
      {"B2PIPLYP",FUNCTIONALS::B2PIPLYP},
      //      {"MPW2PLYP",FUNCTIONALS::MPW2PLYP},
      //      {"MPW2KPLYP",FUNCTIONALS::MPW2KPLYP},
      {"B2PPW91",FUNCTIONALS::B2PPW91},
      {"DSDBLYP",FUNCTIONALS::DSDBLYP},
      {"DUT",FUNCTIONALS::DUT},
      {"PUT",FUNCTIONALS::PUT},
      {"DSDPBEP86",FUNCTIONALS::DSDPBEP86},
      {"PWPB95",FUNCTIONALS::PWPB95},
      {"SAOP", FUNCTIONALS::SAOP},
      {"TF",FUNCTIONALS::TF},
      {"TW",FUNCTIONALS::TW},
      {"PW91K",FUNCTIONALS::PW91K},
      {"LLP91",FUNCTIONALS::LLP91K},
      {"LLP91S", FUNCTIONALS::LLP91KS},
      {"PBE2", FUNCTIONALS::PBE2K},
      {"PBE2S", FUNCTIONALS::PBE2KS},
      {"PBE3", FUNCTIONALS::PBE3K},
      {"PBE4", FUNCTIONALS::PBE4K},
      {"E2000", FUNCTIONALS::E2000K},
      {"HF", FUNCTIONALS::HF}
  };
  check(m,value,field);
}  

/**
 * Readily prepared exchange-correlation functionals (are resolved to instances of Functional)
 * in a way a user would think of a functional. These are typically composed of different parts,
 * see FunctionalClassResolver (here you also find the BASIC_FUNCTIONALS).\n
 * ATTENTION: It is VERY important that the enum numbers here match the corresponding numbers
 * in the FUNCTIONALS enum.\n
 *
 */
enum class XCFUNCTIONALS {
  // No functional to be used
  NONE=0,
      // LDA
      SLATER=1,
      VWN3=2,
      VWN5=3,
      LDA=4,
      LDAERF=5,
      // GGA
      B97=101,
      B97_1=102,
      B97_2=103,
      OLYP=104,
      BLYP=105,
      PBE=106,
      BP86=107,
      KT1=108,
      KT2=109,
      KT3=110,
      PW91=111,
      // Meta-GGA
      // Hybrid
      PBE0=301,
      B3LYP=302,
      B3LYP_G=303,
      B3P86=304,
      B3P86_G=305,
      BPW91=306,
      HF=307,
      // Meta-Hybrid
      // Range-Separated Hybrid
      CAMB3LYP=501,
      LCBLYP=502,
      LCBLYP_047=503,
      // Double-Hybrid
      B2PLYP=601,
      B2KPLYP=602,
      B2TPLYP=603,
      B2GPPLYP=604,
      ROB2PLYP=605,
      B2PIPLYP=606,
      //    MPW2PLYP=607,
      //    MPW2KPLYP=608,
      B2PPW91=609,
      //SCS Double-Hybrid
      DSDBLYP=610,
      DUT=611,
      PUT=612,
      DSDPBEP86=613,
      //SOS Double-Hybrid
      PWPB95=614,
      // Local-Hybrid
      // Model
      SAOP=801
};
template<> inline void resolve<XCFUNCTIONALS>(std::string& value, XCFUNCTIONALS& field){
  static const std::map<std::string, XCFUNCTIONALS> m = {
      {"NONE",XCFUNCTIONALS::NONE},
      {"SLATER",XCFUNCTIONALS::SLATER},
      {"VWN3",XCFUNCTIONALS::VWN3},
      {"VWN5",XCFUNCTIONALS::VWN5},
      {"LDAERF",XCFUNCTIONALS::LDAERF},
      {"LDA",XCFUNCTIONALS::LDA},
      {"B97",XCFUNCTIONALS::B97},
      {"B97_1",XCFUNCTIONALS::B97_1},
      {"B97_2",XCFUNCTIONALS::B97_2},
      {"OLYP",XCFUNCTIONALS::OLYP},
      {"BLYP",XCFUNCTIONALS::BLYP},
      {"PBE",XCFUNCTIONALS::PBE},
      {"BP86",XCFUNCTIONALS::BP86},
      {"KT1",XCFUNCTIONALS::KT1},
      {"KT2",XCFUNCTIONALS::KT2},
      {"KT3",XCFUNCTIONALS::KT3},
      {"PW91",XCFUNCTIONALS::PW91},
      {"PBE0",XCFUNCTIONALS::PBE0},
      {"B3LYP",XCFUNCTIONALS::B3LYP},
      {"B3LYP_G",XCFUNCTIONALS::B3LYP_G},
      {"B3P86",XCFUNCTIONALS::B3P86},
      {"B3P86_G",XCFUNCTIONALS::B3P86_G},
      {"BPW91",XCFUNCTIONALS::BPW91},
      {"CAMB3LYP",XCFUNCTIONALS::CAMB3LYP},
      {"LCBLYP",XCFUNCTIONALS::LCBLYP},
      {"LCBLYP_047",XCFUNCTIONALS::LCBLYP_047},
      {"B2PLYP",XCFUNCTIONALS::B2PLYP},
      {"B2KPLYP",XCFUNCTIONALS::B2KPLYP},
      {"B2TPLYP",XCFUNCTIONALS::B2TPLYP},
      {"B2GPPLYP",XCFUNCTIONALS::B2GPPLYP},
      {"ROB2PLYP",XCFUNCTIONALS::ROB2PLYP},
      {"B2PIPLYP",XCFUNCTIONALS::B2PIPLYP},
      //    {"MPW2PLYP",XCFUNCTIONALS::MPW2PLYP},
      //    {"MPW2KPLYP",XCFUNCTIONALS::MPW2KPLYP},
      {"B2PPW91",XCFUNCTIONALS::B2PPW91},
      //    {"XYG3",XCFUNCTIONALS::XYG3},
      {"DSDBLYP",XCFUNCTIONALS::DSDBLYP},
      {"DUT",XCFUNCTIONALS::DUT},
      {"PUT",XCFUNCTIONALS::PUT},
      {"DSDPBEP86",XCFUNCTIONALS::DSDPBEP86},
      {"PWPB95",XCFUNCTIONALS::PWPB95},
      {"SAOP", XCFUNCTIONALS::SAOP},
      {"HF", XCFUNCTIONALS::HF}
  };
  check(m,value,field);
}

/**
 * Readily prepared kinetic functionals (are resolved to instances of Functional)
 * in a way a user would think of a functional. These are typically composed of different parts,
 * see FunctionalClassResolver (here you also find the BASIC_FUNCTIONALS).\n
 * ATTENTION: It is VERY important that the enum numbers here match the corresponding numbers
 * in the FUNCTIONALS enum.\n
 *
 */
enum class KINFUNCTIONALS {
  // No functional to be used
  NONE=0,
      // LDA
      TF=1001,
      // GGA
      TW=1101,
      PW91K=1102,
      LLP91K=1103,
      LLP91KS=1104,
      PBE2K=1105,
      PBE2KS=1106,
      PBE3K=1107,
      PBE4K=1108,
      E2000K=1109,
};
template<> inline void resolve<KINFUNCTIONALS>(std::string& value, KINFUNCTIONALS& field){
  static const std::map<std::string, KINFUNCTIONALS> m = {
      {"NONE",KINFUNCTIONALS::NONE},
      {"TF",KINFUNCTIONALS::TF},
      {"TW",KINFUNCTIONALS::TW},
      {"PW91K",KINFUNCTIONALS::PW91K},
      {"LLP91",KINFUNCTIONALS::LLP91K},
      {"LLP91S", KINFUNCTIONALS::LLP91KS},
      {"PBE2", KINFUNCTIONALS::PBE2K},
      {"PBE2S", KINFUNCTIONALS::PBE2KS},
      {"PBE3", KINFUNCTIONALS::PBE3K},
      {"PBE4", KINFUNCTIONALS::PBE4K},
      {"E00", KINFUNCTIONALS::E2000K}
  };
  check(m,value,field);
}

/**
 * Correction levels for the DFT-D corrections by Grimme et. al.
 * This includes all corrections (energy, gradient, hessian)
 *
 * NONE:     No correction.
 * D3:       The third set of parameters (with zero damping, also called D3(0)).
 * D3ABC:    The third set of parameters (with zero damping, also called D3(0)) and 3 center correction term.
 * D3BJ:     The third set of parameters with Becke-Johnson damping.
 * D3BJiABC: The third set of parameters with Becke-Johnson damping and 3 center correction term.
 */
enum class DFT_DISPERSION_CORRECTIONS {NONE=0, D3=1, D3ABC=2, D3BJ=3, D3BJABC=4};
template<> inline void resolve<DFT_DISPERSION_CORRECTIONS>(std::string& value, DFT_DISPERSION_CORRECTIONS& field){
  static const std::map<std::string, DFT_DISPERSION_CORRECTIONS> m = {
      {"NONE",DFT_DISPERSION_CORRECTIONS::NONE},
      {"D3",DFT_DISPERSION_CORRECTIONS::D3},
      {"D3ABC",DFT_DISPERSION_CORRECTIONS::D3ABC},
      {"D3BJ",DFT_DISPERSION_CORRECTIONS::D3BJ},
      {"D3BJABC",DFT_DISPERSION_CORRECTIONS::D3BJABC}
  };
  check(m,value,field);
}
/**
 * What the DensityFunctionalController is needed for. In a SCF calculation, only the potential is needed,
 * while in a post SCF calculation one may only need higher functional derivatives of the exchange
 * correlation functional.
 *
 * SCF:   Quantities for the exchange correlation potential
 * LRSCF: Derivatives for the exchange correlation kernel
 */
enum class DFT_CONTROLLER_PURPOSES{SCF=0,LRSCF=1,ALL=2};
/**
 * How the initial guess for electronic structure calculations is determined\n
 * H_CORE:    hCore guess, i.e. ignoring electron-electron interactions,
 *            points to HCoreGuessCalculator\n
 * EHT:       Extended Hueckel Theory guess, points to ExtendedHueckel\n
 * ATOM_DENS: A guess of combined atomic densities, points to AtomicDensityGuessCalculator
 * ATOM_SCF:  Similar to ATOM_DENS, but the atomic densities themselves are calculated with a proper
 *            SCF for each atom type. This is very similar to what ADF does.
 */
enum class INITIAL_GUESSES {H_CORE=0, EHT=1, ATOM_DENS=2, ATOM_SCF=3, ATOM_SCF_INPLACE=4, SAP=5};
template<> inline void resolve<INITIAL_GUESSES>(std::string& value, INITIAL_GUESSES& field){
  static const std::map<std::string, INITIAL_GUESSES> m = {
      {"HCORE",INITIAL_GUESSES::H_CORE},
      {"EHT",INITIAL_GUESSES::EHT},
      {"ATOM_DENS",INITIAL_GUESSES::ATOM_DENS},
      {"ATOM_SCF",INITIAL_GUESSES::ATOM_SCF},
      {"ATOM_SCF_INPLACE",INITIAL_GUESSES::ATOM_SCF_INPLACE},
      {"SAP",INITIAL_GUESSES::SAP}
  };
  check(m,value,field);
}

/**
 * Which algorithm to use for the damping method.\n
 * NONE: No damping will be applied.\n
 * STATIC: Static damping using a constant factor.\n
 */
enum class DAMPING_ALGORITHMS {NONE=0, STATIC=1, SERIES=2, DYNAMIC=3};
template<> inline void resolve<DAMPING_ALGORITHMS>(std::string& value, DAMPING_ALGORITHMS& field){
  static const std::map<std::string, DAMPING_ALGORITHMS> m = {
      {"NONE",DAMPING_ALGORITHMS::NONE},
      {"STATIC",DAMPING_ALGORITHMS::STATIC},
      {"SERIES",DAMPING_ALGORITHMS::SERIES},
      {"DYNAMIC",DAMPING_ALGORITHMS::DYNAMIC}
  };
  check(m,value,field);
}

/**
 * Possible types for LRSCF problems to determine eigenvalue solving procedure:
 * TDA/CIS : AX = Xw (is Hermitian, uses Davidson)
 * TDDFT : sqrt(A-B)(A+B)sqrt(A-B) sqrt-(X+Y)
 *           = sqrt-(X+Y) w^2 (is Hermitian, uses Davidson)
 * RPA : (A+B)(X+Y) = (X-Y)w
 *       (A-B)(X-Y) = (X+Y)w (is non-Hermitian, uses a modifed OJJ)
 * The latter one includes TDHF, hybrid TDDFT, FDEc with 
 * external orthogonality and supermolecular TDDFT with local orbitals.
 */
enum class RESPONSE_PROBLEM {TDA=0,TDDFT=1,RPA=2};
template<> inline void resolve<RESPONSE_PROBLEM>(std::string& value, RESPONSE_PROBLEM& field){
  static const std::map<std::string, RESPONSE_PROBLEM> m = {
    {"TDA",RESPONSE_PROBLEM::TDA},
    {"TDDFT",RESPONSE_PROBLEM::TDDFT},
    {"RPA",RESPONSE_PROBLEM::RPA},
  };
  check(m,value,field);
}

enum class LRSCF_TYPE {ISOLATED=0, UNCOUPLED=1, COUPLED=2};
template<> inline void resolve<LRSCF_TYPE>(std::string& value, LRSCF_TYPE& field){
  static const std::map<std::string, LRSCF_TYPE> m = {
    {"ISOLATED",LRSCF_TYPE::ISOLATED},
    {"ISO",LRSCF_TYPE::ISOLATED},
    {"UNCOUPLED",LRSCF_TYPE::UNCOUPLED},
    {"FDEU",LRSCF_TYPE::UNCOUPLED},
    {"COUPLED",LRSCF_TYPE::COUPLED},
    {"FDEC",LRSCF_TYPE::COUPLED},
  };
  check(m,value,field);
}

enum class INTEGRAL_TYPE {NUMERICAL=0, ANALYTICAL=1};
template<> inline void resolve<INTEGRAL_TYPE>(std::string& value, INTEGRAL_TYPE& field){
  static const std::map<std::string, INTEGRAL_TYPE> m = {
      {"NUMERICAL",INTEGRAL_TYPE::NUMERICAL},
      {"ANALYTICAL",INTEGRAL_TYPE::ANALYTICAL}
  };
  check(m,value,field);
}

enum class GAUGE {LENGTH=0, VELOCITY=1};
template<> inline void resolve<GAUGE>(std::string& value, GAUGE& field){
  static const std::map<std::string, GAUGE> m = {
      {"LENGTH",GAUGE::LENGTH},
      {"VELOCITY",GAUGE::VELOCITY}
  };
  check(m,value,field);
}

/**
 * Multiplicity of excited state in LRSCF calculation
 * Singlet: 2*S + 1 = 1
 * Triplet: 2*S + 1 = 3
 */
enum class MULTIPLICITY {SINGLET=0,TRIPLET=1};
template<> inline void resolve<MULTIPLICITY>(std::string& value, MULTIPLICITY& field){
  static const std::map<std::string, MULTIPLICITY> m = {
      {"SINGLET",MULTIPLICITY::SINGLET},
      {"TRIPLET",MULTIPLICITY::TRIPLET}
  };
  check(m,value,field);
}
/**************************************************************************************************/
/*                                           Basis                                                */
/**************************************************************************************************/
/**
 * For what a certain basis set should be used.\n
 * DEFAULT:          The basis set in which the electronic structure is calculated and expressed\n
 * AUX_COULOMB:      The auxiliary basis for a density fitting to evaluate Coulombic electron-electron
 *                   interactions
 * HUECKEL:          The basis set used for semiempirical calculations, most probably a minimal basis
 * IAO_LOCALIZATION: The basis set in which the Intrinsic Atomic Orbitals will be expressed
 * AUX_CORREL:       The auxiliary basis for density fitting to evaluate electron-electron correlation
 *                   contributions
 */
enum class BASIS_PURPOSES {DEFAULT = 0, AUX_COULOMB=1, MINBAS=2, HUECKEL=3, IAO_LOCALIZATION=4, SCF_DENS_GUESS=5, AUX_CORREL=6};
template<> inline void resolve<BASIS_PURPOSES>(std::string& value, BASIS_PURPOSES& field){
  static const std::map<std::string, BASIS_PURPOSES> m = {
      {"DEFAULT",BASIS_PURPOSES::DEFAULT},
      {"AUX_COULOMB",BASIS_PURPOSES::AUX_COULOMB},
      {"MINBAS",BASIS_PURPOSES::MINBAS},
      {"HUECKEL",BASIS_PURPOSES::HUECKEL},
      {"IAO_LOCALIZATION",BASIS_PURPOSES::IAO_LOCALIZATION},
      {"SCF_DENS_GUESS",BASIS_PURPOSES::SCF_DENS_GUESS},
      {"AUX_CORREL",BASIS_PURPOSES::AUX_CORREL}
  };
  check(m,value,field);
}

/**
 * How a Hartree-Potential should be calculated.
 * (Old explanation:
 * In case of a split Coulomb and exchange part in the Fock update (e.g. in DFT no Exchange is
 * needed) these are the possibilities to evaluate the Coulomb part. In principle all complete
 * HARTREE_FOCK_POTENTIAL_CALCULATORS are possible, but also the famous RI approximation, which is a
 * way to reduce the scaling behaviour of DFT from O(n^4) to O(n^3).)
 */
enum class DENS_FITS{RI=0,NONE=1};
template<> inline void resolve<DENS_FITS>(std::string& value, DENS_FITS& field){
  static const std::map<std::string, DENS_FITS> m = {
      {"RI",DENS_FITS::RI},
      {"NORI",DENS_FITS::NONE},
      {"NONE",DENS_FITS::NONE}
  };
  check(m,value,field);
}
/**************************************************************************************************/
/*                                            Grid                                                */
/**************************************************************************************************/
/**
 * How the radial part of the integration grid(s) is set up\n
 * EQUI: equidistant points
 * All others are schemes named after their authors, 
 *  see the grid routines for the references.
 */
enum class RADIAL_GRID_TYPES {BECKE=0, HANDY=1, AHLRICHS=2, KNOWLES=3, EQUI=4//, MULTIEXP, EQUI
};
template<> inline void resolve<RADIAL_GRID_TYPES>(std::string& value, RADIAL_GRID_TYPES& field){
  static const std::map<std::string, RADIAL_GRID_TYPES> m = {
      {"BECKE",RADIAL_GRID_TYPES::BECKE},
      {"HANDY",RADIAL_GRID_TYPES::HANDY},
      {"AHLRICHS",RADIAL_GRID_TYPES::AHLRICHS},
      {"KNOWLES",RADIAL_GRID_TYPES::KNOWLES},
      {"EQUIDISTAND",RADIAL_GRID_TYPES::EQUI}
  };
  check(m,value,field);
}
/**
 * How the spherical part of the integration grid(s) is set up
 */
enum class SPHERICAL_GRID_TYPES {LEBEDEV=0};
template<> inline void resolve<SPHERICAL_GRID_TYPES>(std::string& value, SPHERICAL_GRID_TYPES& field){
  static const std::map<std::string, SPHERICAL_GRID_TYPES> m = {
      {"LEBEDEV",SPHERICAL_GRID_TYPES::LEBEDEV}
  };
  check(m,value,field);
}
/**
 * How the atomic cells are determined - for each atom a grid is set up, they are then added
 * upon each other and recieve an additional weight prefactor based on the partitioning\n
 * BECKE:   the fuzzy cells (see A.D. Becke, J.Chem.Phys. 88 (1988), 2547.)\n
 * VORONOI: sharp cuts between the atoms (it is contrasted to the fuzzy cells in the paper above).
 */
enum class GRID_TYPES {BECKE=0, VORONOI=1, SSF=2};
template<> inline void resolve<GRID_TYPES>(std::string& value, GRID_TYPES& field){
  static const std::map<std::string, GRID_TYPES> m = {
      {"BECKE",GRID_TYPES::BECKE},
      {"VORONOI",GRID_TYPES::VORONOI},
      {"SSF",GRID_TYPES::SSF}
  };
  check(m,value,field);
}
/**
 * For what a certain numerical grid should be used\n
 * DEFAULT: The standard integration grid to, e.g, calculate the DFT energy\n
 * SMALL:   A smaller integration grid for faster computation. This is typically used during an SCF
 *          in DFT calculations as long as the SCF is not converged\n
 * PLOT:    A grid to plot data on to for visualization purposes. Typically a cubical (equidistant)
 *          grid.
 */
enum class GRID_PURPOSES {DEFAULT = 0, SMALL=1, PLOT=2};
template<> inline void resolve<GRID_PURPOSES>(std::string& value, GRID_PURPOSES& field){
  static const std::map<std::string, GRID_PURPOSES> m = {
      {"DEFAULT",GRID_PURPOSES::DEFAULT},
      {"SMALL",GRID_PURPOSES::SMALL},
      {"PLOT",GRID_PURPOSES::PLOT}
  };
  check(m,value,field);
}
/**************************************************************************************************/
/*                                         Optimization                                           */
/**************************************************************************************************/
/**
 * Which algorithm to use for the geometry optimization.\n
 * SD: Steepest descent.\n
 * BFGS: Broyden-Fletcher-Goldfarb-Shanno.\n
 */
enum class OPTIMIZATION_ALGORITHMS {SD=0, BFGS=1};
template<> inline void resolve<OPTIMIZATION_ALGORITHMS>(std::string& value, OPTIMIZATION_ALGORITHMS& field){
  static const std::map<std::string, OPTIMIZATION_ALGORITHMS> m = {
      {"SD",OPTIMIZATION_ALGORITHMS::SD},
      {"BFGS",OPTIMIZATION_ALGORITHMS::BFGS}
  };
  check(m,value,field);
}
/**************************************************************************************************/
/*                                    Geometry optimization                                       */
/**************************************************************************************************/
/**
 * How to calculate geometry gradients.\n
 * NUMERICAL: Calculate them numerically with finite steps.\n
 * ANALYTICAL: Calculate them analytically.
 */
enum class GRADIENT_TYPES {NUMERICAL=0, ANALYTICAL=1};
template<> inline void resolve<GRADIENT_TYPES>(std::string& value, GRADIENT_TYPES& field){
  static const std::map<std::string, GRADIENT_TYPES> m = {
      {"NUMERICAL",GRADIENT_TYPES::NUMERICAL},
      {"ANALYTICAL",GRADIENT_TYPES::ANALYTICAL}
  };
  check(m,value,field);
}

/**
 * Which algorithm to use for the geometry optimization.\n
 * GROUNDSTATE: Ground state optimization.\n
 * TS:          Transition state optimization
 */
enum class GEOMETRY_OPTIMIZATION_TYPES {GROUNDSTATE=0,TS=1};
template<> inline void resolve<GEOMETRY_OPTIMIZATION_TYPES>(std::string& value, GEOMETRY_OPTIMIZATION_TYPES& field){
  static const std::map<std::string, GEOMETRY_OPTIMIZATION_TYPES> m = {
      {"GROUNDSTATE",GEOMETRY_OPTIMIZATION_TYPES::GROUNDSTATE},
      {"TS",GEOMETRY_OPTIMIZATION_TYPES::TS}
  };
  check(m,value,field);
}

/**************************************************************************************************/
/*                                    Geometry optimization                                       */
/**************************************************************************************************/
/**
 * How to calculate the Hessian.\n
 * NUMERICAL: Calculate numerically with finite steps.\n
 * ANALYTICAL: Calculate analytically.
 */
enum class HESSIAN_TYPES {NUMERICAL=0,ANALYTICAL=1};
template<> inline void resolve<HESSIAN_TYPES>(std::string& value, HESSIAN_TYPES& field){
  static const std::map<std::string, HESSIAN_TYPES> m = {
      {"NUMERICAL",HESSIAN_TYPES::NUMERICAL},
      {"ANALYTICAL",HESSIAN_TYPES::ANALYTICAL}
  };
  check(m,value,field);
}

/**************************************************************************************************/
/*                                       Localization                                             */
/**************************************************************************************************/


enum class ORBITAL_LOCALIZATION_ALGORITHMS {PIPEK_MEZEY=0, FOSTER_BOYS=1, IAO=2, IBO=3,EDMISTON_RUEDENBERG=4,NON_ORTHOGONAL=5,NONE=6};
template<> inline void resolve<ORBITAL_LOCALIZATION_ALGORITHMS>(std::string& value, ORBITAL_LOCALIZATION_ALGORITHMS& field){
  static const std::map<std::string, ORBITAL_LOCALIZATION_ALGORITHMS> m = {
      {"PM",ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY},
      {"FB",ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS},
      {"IAO",ORBITAL_LOCALIZATION_ALGORITHMS::IAO},
      {"IBO",ORBITAL_LOCALIZATION_ALGORITHMS::IBO},
      {"ER",ORBITAL_LOCALIZATION_ALGORITHMS::EDMISTON_RUEDENBERG},
      {"NO",ORBITAL_LOCALIZATION_ALGORITHMS::NON_ORTHOGONAL},
      {"NONE",ORBITAL_LOCALIZATION_ALGORITHMS::NONE}
  };
  check(m,value,field);
}

enum class KIN_EMBEDDING_MODES {
  NONE=0, NADD_FUNC=1, LEVELSHIFT=2, HUZINAGA=3, HOFFMANN=4, RECONSTRUCTION=5, FERMI_SHIFTED_HUZINAGA=6};
template<> inline void resolve<KIN_EMBEDDING_MODES>(std::string& value, KIN_EMBEDDING_MODES& field){
  static const std::map<std::string, KIN_EMBEDDING_MODES> m = {
      {"NONE",KIN_EMBEDDING_MODES::NONE},
      {"NADDFUNC",KIN_EMBEDDING_MODES::NADD_FUNC},
      {"LEVELSHIFT",KIN_EMBEDDING_MODES::LEVELSHIFT},
      {"HUZINAGA",KIN_EMBEDDING_MODES::HUZINAGA},
      {"HOFFMANN",KIN_EMBEDDING_MODES::HOFFMANN},
      {"RECONSTRUCTION",KIN_EMBEDDING_MODES::RECONSTRUCTION},
      {"FERMI_SHIFTED_HUZINAGA",KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA},
      {"FERMI",KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA}
  };
  check(m,value,field);
}
/**************************************************************************************************/
/*                                Basis Set Truncation Algorithms								  */
/**************************************************************************************************/
enum class BASIS_SET_TRUNCATION_ALGORITHMS {
  NONE=0, NET_POPULATION=1, PRIMITIVE_NET_POPULATION=2
};
template<> inline void resolve<BASIS_SET_TRUNCATION_ALGORITHMS>(std::string& value, BASIS_SET_TRUNCATION_ALGORITHMS& field) {
  static const std::map<std::string, BASIS_SET_TRUNCATION_ALGORITHMS> m = {
      {"NONE", BASIS_SET_TRUNCATION_ALGORITHMS::NONE},
      {"NETPOP", BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION},
      {"NETPOPULATION", BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION},
      {"PRIMITIVENETPOP", BASIS_SET_TRUNCATION_ALGORITHMS::PRIMITIVE_NET_POPULATION},
  };
  check(m,value,field);
}

/**************************************************************************************************/
/*                                Projection Operators          								  */
/**************************************************************************************************/
enum class PROJECTION_OPERATORS {
  LEVELSHIFT=0, HUZINAGA=1, HOFFMANN=2, NONE=3
};
template<> inline void resolve<PROJECTION_OPERATORS>(std::string& value, PROJECTION_OPERATORS& field) {
  static const std::map<std::string, PROJECTION_OPERATORS> m = {
      {"LEVELSHIFT", PROJECTION_OPERATORS::LEVELSHIFT},
      {"HUZINAGA", PROJECTION_OPERATORS::HUZINAGA},
      {"HOFFMANN", PROJECTION_OPERATORS::HOFFMANN},
      {"NONE", PROJECTION_OPERATORS::NONE},
  };
  check(m,value,field);
}

enum class CC_LEVEL {CCSD, CCSD_T};
template<> inline void resolve<CC_LEVEL>(std::string& value, CC_LEVEL& field){
  static const std::map<std::string, CC_LEVEL> m = {
      {"CCSD",CC_LEVEL::CCSD},
      {"CCSD(T)",CC_LEVEL::CCSD_T}
  };
  check(m,value,field);
}

enum class GAUGE_ORIGIN {COM, ORIGIN};
template<> inline void resolve<GAUGE_ORIGIN>(std::string& value, GAUGE_ORIGIN& field){
  static const std::map<std::string, GAUGE_ORIGIN> m = {
      {"CENTEROFMASS",GAUGE_ORIGIN::COM},
      {"ORIGIN",GAUGE_ORIGIN::ORIGIN}
  };
  check(m,value,field);
}

}
static constexpr Options::SCF_MODES RESTRICTED = Options::SCF_MODES::RESTRICTED;
static constexpr Options::SCF_MODES UNRESTRICTED = Options::SCF_MODES::UNRESTRICTED;
} /* namespace Serenity */

#endif /* OPTIONS_H_ */

