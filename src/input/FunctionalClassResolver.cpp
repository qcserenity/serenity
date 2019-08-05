/**
 * @file   FunctionalClassResolver.cpp
 *
 * @date   May 28, 2014
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
/* Include Class Header*/
#include "input/FunctionalClassResolver.h"
/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "settings/Options.h"

/* Further includes */

namespace Serenity {
using namespace Options;

Functional FunctionalClassResolver::resolveFunctional(FUNCTIONALS functional) {
  switch(functional) {
    case FUNCTIONALS::NONE:
      return Functional(
          {BASIC_FUNCTIONALS::NONE},
          {0.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::LDA:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX, BASIC_FUNCTIONALS::VWN5C},
          {1.0,1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::BLYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
          {1.0,1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::PBE:
      return Functional(
          {BASIC_FUNCTIONALS::PBEX, BASIC_FUNCTIONALS::PBEC},
          {1.0,1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::BP86:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::P86C},
          {1.0,1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::KT1:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX, BASIC_FUNCTIONALS::KTX, BASIC_FUNCTIONALS::VWN5C},
          {1.0,-0.006,1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::KT2:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX, BASIC_FUNCTIONALS::KTX, BASIC_FUNCTIONALS::VWN5C},
          {1.07173,-0.006,0.576727},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::KT3:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX, BASIC_FUNCTIONALS::KTX, BASIC_FUNCTIONALS::OPTXCORR, BASIC_FUNCTIONALS::LYPC},
          {1.092,-0.004,-0.925452,0.864409},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::LDAERF:
      return Functional(
          {BASIC_FUNCTIONALS::LDAERFX, BASIC_FUNCTIONALS::LDAERFC},
          {1.0,1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::PBE0:
      return Functional(
          {BASIC_FUNCTIONALS::PBEX, BASIC_FUNCTIONALS::PBEC},
          {0.75,1.0},
          0.25,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
    case FUNCTIONALS::B3LYP:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX, BASIC_FUNCTIONALS::BECKECORRX, BASIC_FUNCTIONALS::LYPC, BASIC_FUNCTIONALS::VWN5C},
          {0.80,0.72,0.81,0.19},
          0.20,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B3LYP_G:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX, BASIC_FUNCTIONALS::BECKECORRX, BASIC_FUNCTIONALS::LYPC, BASIC_FUNCTIONALS::VWN3C},
          {0.80,0.72,0.81,0.19},
          0.20,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B3P86:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX, BASIC_FUNCTIONALS::BECKECORRX, BASIC_FUNCTIONALS::P86CORRC, BASIC_FUNCTIONALS::VWN5C},
          {0.80,0.72,0.81,1.0},
          0.20,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B3P86_G:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX, BASIC_FUNCTIONALS::BECKECORRX, BASIC_FUNCTIONALS::P86CORRC, BASIC_FUNCTIONALS::VWN3C},
          {0.80,0.72,0.81,1.0},
          0.20,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::BPW91:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::PW91C},
          {1.0,1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B97:
      return Functional(
          {BASIC_FUNCTIONALS::B97X, BASIC_FUNCTIONALS::B97C},
          {1.0,1.0},
          0.1943,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B97_1:
      return Functional(
          {BASIC_FUNCTIONALS::B97_1X, BASIC_FUNCTIONALS::B97_1C},
          {1.0,1.0},
          0.21,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B97_2:
      return Functional(
          {BASIC_FUNCTIONALS::B97_2X, BASIC_FUNCTIONALS::B97_2C},
          {1.0,1.0},
          0.21,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::CAMB3LYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKECAMX, BASIC_FUNCTIONALS::VWN5C, BASIC_FUNCTIONALS::LYPC},
          {1.0,0.19,0.81},
          0.19,
          0.0,
          0.46,
          0.33,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::VWN5:
      return Functional(
          {BASIC_FUNCTIONALS::VWN5C},
          {1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::VWN3:
      return Functional(
          {BASIC_FUNCTIONALS::VWN3C},
          {1.0},
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::SLATER:
      return Functional(
          {BASIC_FUNCTIONALS::SLATERX},
          {1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::OLYP:
      return Functional(
          {BASIC_FUNCTIONALS::LYPC, BASIC_FUNCTIONALS::OPTX},
          {1.0, 1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::TF:
      return Functional(
          {BASIC_FUNCTIONALS::TFK},
          {1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::TW:
      return Functional(
          {BASIC_FUNCTIONALS::TW},
          {1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::LCBLYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKECAMX, BASIC_FUNCTIONALS::LYPC},
          {1.0,1.0},
          0.0,
          0.0,
          1.0,
          0.33,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::PW91K:
      return Functional(
          {BASIC_FUNCTIONALS::PW91K},
          {1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::PW91:
      return Functional(
          {BASIC_FUNCTIONALS::PW91C, BASIC_FUNCTIONALS::PW91X},
          {1.0,1.0},
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B2PLYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
          {0.47,0.73},
          0.53,
          0.27,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B2KPLYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
          {0.28,0.58},
          0.72, //a_x
          0.42, //a_c
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B2TPLYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
          {0.40,0.69},
          0.60,
          0.31,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B2GPPLYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
          {0.35,0.64},
          0.65,
          0.36,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::ROB2PLYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
          {0.41,0.72},
          0.59,
          0.28,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::B2PIPLYP:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
          {0.398,0.727},
          0.602,
          0.273,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
//    case FUNCTIONALS::MPW2PLYP:
//      return Functional(
//          {BASIC_FUNCTIONALS::MPW91X, BASIC_FUNCTIONALS::LYPC},
//          {0.45,0.75},
//          0.55,
//          0.25,
//          0.0,
//          0.0,
//          1.0,
//          1.0);
//      break;
//    case FUNCTIONALS::MPW2KPLYP:
//      return Functional(
//          {BASIC_FUNCTIONALS::MPW91X, BASIC_FUNCTIONALS::LYPC},
//          {0.28,0.58},
//          0.72,
//          0.42,
//          0.0,
//          0.0,
//          1.0,
//          1.0);
//      break;
    case FUNCTIONALS::B2PPW91:
      return Functional(
          {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::PW91C},
          {0.8,0.9},
          0.2,
          0.1,
          0.0,
          0.0,
          1.0,
          1.0);
      break;
    case FUNCTIONALS::DSDBLYP:
          return Functional(
              {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
              {0.31,0.54},
              0.69,
              1.0,
              0.0,
              0.0,
              0.37,
              0.46);
          break;
    case FUNCTIONALS::DUT:
          return Functional(
              {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
              {0.30,0.59},
              0.70,
              1.0,
              0.0,
              0.0,
              0.36,
              0.47);
          break;
    case FUNCTIONALS::PUT:
          return Functional(
              {BASIC_FUNCTIONALS::BECKEX, BASIC_FUNCTIONALS::LYPC},
              {0.32,0.63},
              0.68,
              1.0,
              0.0,
              0.0,
              0.27,
              0.46);
          break;
    case FUNCTIONALS::DSDPBEP86:
          return Functional(
              {BASIC_FUNCTIONALS::PBEX, BASIC_FUNCTIONALS::P86C},
              {0.30,0.43},
              0.70,
              1.0,
              0.0,
              0.0,
              0.25,
              0.53);
          break;
    case FUNCTIONALS::PWPB95:
          return Functional(
              //ToDo implement B95C in XCFUN and use it for PWPB95
              {BASIC_FUNCTIONALS::PW91X, BASIC_FUNCTIONALS::B97C},
              {0.50,0.731},
              0.50,
              1.0,
              0.0,
              0.0,
              0.00,
              0.269);
          break;
    case FUNCTIONALS::SAOP:
      return Functional(
          {BASIC_FUNCTIONALS::SAOP});
      break;
    case FUNCTIONALS::LLP91K:
    	return Functional(
    			{BASIC_FUNCTIONALS::LLP91K},
				{1.0},
				0.0,
				0.0,
				0.0,
				0.0,
        		1.0,
        		1.0);
    	break;
    case FUNCTIONALS::LLP91KS:
    	return Functional(
    			{BASIC_FUNCTIONALS::LLP91KS},
				{1.0},
				0.0,
				0.0,
				0.0,
				0.0,
        		1.0,
        		1.0);
    	break;
    case FUNCTIONALS::PBE2K:
    	return Functional(
    			{BASIC_FUNCTIONALS::PBE2K},
				{1.0},
				0.0,
				0.0,
				0.0,
				0.0,
        		1.0,
        		1.0);
    	break;
    case FUNCTIONALS::PBE2KS:
    	return Functional(
    			{BASIC_FUNCTIONALS::PBE2KS},
				{1.0},
				0.0,
				0.0,
				0.0,
				0.0,
        		1.0,
        		1.0);
    	break;
    case FUNCTIONALS::PBE3K:
    	return Functional(
    			{BASIC_FUNCTIONALS::PBE3K},
				{1.0},
				0.0,
				0.0,
				0.0,
				0.0,
        		1.0,
        		1.0);
    	break;
    case FUNCTIONALS::PBE4K:
    	return Functional(
    			{BASIC_FUNCTIONALS::PBE4K},
				{1.0},
				0.0,
				0.0,
				0.0,
				0.0,
        		1.0,
        		1.0);
    	break;
    case FUNCTIONALS::E2000K:
    	return Functional(
    			{BASIC_FUNCTIONALS::E2000K},
				{1.0},
				0.0,
				0.0,
				0.0,
				0.0,
        		1.0,
        		1.0);
    	break;
    default:
      assert(false && "Functional can not be resolved, missing resolve function");
      return Functional({BASIC_FUNCTIONALS::NONE},{1.0},0.0,0.0,0.0);
      break;
  }
}

Functional FunctionalClassResolver::resolveFunctional(Options::XCFUNCTIONALS functional){
  return resolveFunctional(Options::FUNCTIONALS(int(functional)));
}

Functional FunctionalClassResolver::resolveFunctional(Options::KINFUNCTIONALS functional){
  return resolveFunctional(Options::FUNCTIONALS(int(functional)));
}

FUNCTIONAL_CLASSES FunctionalClassResolver::resolveFunctionalClass(
    BASIC_FUNCTIONALS functional) {
  switch(functional) {
    case BASIC_FUNCTIONALS::NONE:
      return FUNCTIONAL_CLASSES::NONE;
      break;
    case BASIC_FUNCTIONALS::SLATERX:
    case BASIC_FUNCTIONALS::VWN5C:
    case BASIC_FUNCTIONALS::VWN3C:
    case BASIC_FUNCTIONALS::TFK:
    case BASIC_FUNCTIONALS::LDAERFC:
    case BASIC_FUNCTIONALS::LDAERFX:
    case BASIC_FUNCTIONALS::LDAERFC_JT:
    case BASIC_FUNCTIONALS::OPTX:
    case BASIC_FUNCTIONALS::OPTXCORR:
      return FUNCTIONAL_CLASSES::LDA;
      break;
    case BASIC_FUNCTIONALS::PBEC:
    case BASIC_FUNCTIONALS::PBEX:
    case BASIC_FUNCTIONALS::BECKEX:
    case BASIC_FUNCTIONALS::BECKECORRX:
    case BASIC_FUNCTIONALS::BECKESRX:
    case BASIC_FUNCTIONALS::BECKECAMX:
    case BASIC_FUNCTIONALS::LYPC:
    case BASIC_FUNCTIONALS::REVPBEX:
    case BASIC_FUNCTIONALS::B97X:
    case BASIC_FUNCTIONALS::B97C:
    case BASIC_FUNCTIONALS::B97_1X:
    case BASIC_FUNCTIONALS::B97_1C:
    case BASIC_FUNCTIONALS::B97_2X:
    case BASIC_FUNCTIONALS::B97_2C:
    case BASIC_FUNCTIONALS::P86CORRC:
    case BASIC_FUNCTIONALS::PW91C:
    case BASIC_FUNCTIONALS::PW91K:
    case BASIC_FUNCTIONALS::PW91X:
    case BASIC_FUNCTIONALS::KTX:
    case BASIC_FUNCTIONALS::VWN_PBEC:
    case BASIC_FUNCTIONALS::PW86X:
    case BASIC_FUNCTIONALS::P86C:
    case BASIC_FUNCTIONALS::TW:
    case BASIC_FUNCTIONALS::LLP91K:
    case BASIC_FUNCTIONALS::LLP91KS:
    case BASIC_FUNCTIONALS::PBE2K:
    case BASIC_FUNCTIONALS::PBE2KS:
    case BASIC_FUNCTIONALS::PBE3K:
    case BASIC_FUNCTIONALS::PBE4K:
    case BASIC_FUNCTIONALS::E2000K:

      return FUNCTIONAL_CLASSES::GGA;
      break;
    case BASIC_FUNCTIONALS::SAOP:
      return FUNCTIONAL_CLASSES::MODELL;
      break;
    case BASIC_FUNCTIONALS::M05X:
    case BASIC_FUNCTIONALS::M05X2X:
    case BASIC_FUNCTIONALS::M06X:
    case BASIC_FUNCTIONALS::M06X2X:
    case BASIC_FUNCTIONALS::M06LX:
    case BASIC_FUNCTIONALS::M06HFX:
    case BASIC_FUNCTIONALS::M05X2C:
    case BASIC_FUNCTIONALS::M05C:
    case BASIC_FUNCTIONALS::M06C:
    case BASIC_FUNCTIONALS::M06HFC:
    case BASIC_FUNCTIONALS::M06LC:
    case BASIC_FUNCTIONALS::M06X2C:
    case BASIC_FUNCTIONALS::TPSSC:
    case BASIC_FUNCTIONALS::TPSSX:
    case BASIC_FUNCTIONALS::REVTPSSC:
    case BASIC_FUNCTIONALS::REVTPSSX:
      //ToDo:
      assert(false && "Meta GGA's not implemented");
      return FUNCTIONAL_CLASSES::META_GGA;
      break;
    default:
      //ToDo: For some functionals implemented in XCFun I was too lazy to check it's correct type.
      //      If you want to use that functional, you just have to add this functional here,
      //      everything else (except meta-ggas) should already be implemented.
      throw SerenityError("FunctionalClassResolver::resolveFunctionalClass not implemented for this functional");
  }
}

} /* namespace Serenity */
