/**
 * @file   FunctionalClassResolver.h
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
#ifndef FUNCTIONALCLASSRESOLVER_H_
#define FUNCTIONALCLASSRESOLVER_H_

namespace Serenity {
/* Forward declarations */
class Functional;

namespace Options {
enum class FUNCTIONALS;
enum class XCFUNCTIONALS;
enum class KINFUNCTIONALS;
}
/**
 * Functionals available in XCFun.
 */
enum class BASIC_FUNCTIONALS {
	  NONE=0,
	  SLATERX=1,
	  PW86X=2,
	  VWN3C=3,
	  VWN5C=4,
	  PBEC=5,
	  PBEX=6,
	  PBE2K=7,
	  PBE2KS=8,
	  PBE3K=9,
	  PBE4K=10,
	  BECKEX=11,
	  BECKECORRX=12,
	  BECKESRX=13,
	  BECKECAMX=14,
	  LLP91K=15,
	  LLP91KS=16,
	  BRX=17,
	  BRC=18,
	  BRXC=19,
	  LDAERFX=20,
	  LDAERFC=21,
	  LDAERFC_JT=22,
	  LYPC=23,
	  OPTX=24,
	  OPTXCORR=25,
	  REVPBEX=26,
	  RPBEX=27,
	  SPBEC=28,
	  VWN_PBEC=29,
	  KTX=30,
	  TFK=31,
	  TW=32,
	  PW91X=33,
	  PW91K=34,
	  PW91C=35,
	  M05X=36,
	  M05X2X=37,
	  M06X=38,
	  M06X2X=39,
	  M06LX=40,
	  M06HFX=41,
	  M05X2C=42,
	  M05C=43,
	  M06C=44,
	  M06HFC=45,
	  M06LC=46,
	  M06X2C=47,
	  TPSSC=48,
	  TPSSX=49,
	  REVTPSSC=50,
	  REVTPSSX=51,
	  PZ81C=52,
	  P86C=53,
	  P86CORRC=54,
	  BTK=55,
	  VWK=56,
	  B97X=57,
	  B97C=58,
	  B97_1X=59,
	  B97_1C=60,
	  B97_2X=61,
	  B97_2C=62,
	  CSC=63,
	  APBEC=64,
	  APBEX=65,
	  ZVPBESOLC=66,
	  BLOCX=67,
	  PBEINTC=68,
	  PBEINTX=69,
	  PBELOCC=70,
	  PBESOLX=71,
	  TPSSLOCC=72,
	  ZVPBEINTC=73,
	  SAOP=74,
	  E2000K=75,
	  EXX=76
};
/**
 * Whether a functional is LDA, GGA, meta-GGA, hybrid, ...
 * Each class requires a different internal treatment. E.g. GGAs need to have the density and its
 * gradient, hybrids need to mix in Hartree-Fock exchange, ...
 */
enum class FUNCTIONAL_CLASSES {NONE, LDA, GGA, META_GGA, MODELL};//, HYBRID, META_HYBRID, RANGE_SEPARATED};
/**
 * @class FunctionalClassResolver FunctionalClassResolver.h
 * @brief Class organizing the flavors of different functionals
 */
class FunctionalClassResolver {
  FunctionalClassResolver() = delete;
public:
  static Functional resolveFunctional(Options::FUNCTIONALS functional);

  static Functional resolveFunctional(Options::XCFUNCTIONALS functional);

  static Functional resolveFunctional(Options::KINFUNCTIONALS functional);
  /**
   * @param   functional of which the type is to be determined
   * @returns the class of functional, i.e. whether it is LDA, GGA, ...
   */
  static FUNCTIONAL_CLASSES resolveFunctionalClass(BASIC_FUNCTIONALS functional);
};

} /* namespace Serenity */

#endif /* FUNCTIONALCLASSRESOLVER_H_ */
