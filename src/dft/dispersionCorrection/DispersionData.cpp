/**
 * @file DispersionData.cpp
 *
 * @date   Nov 26, 2015
 * @author Jan Unsleber
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
#include "dft/dispersionCorrection/DispersionData.h"
/* Include Serenity Internal Headers */
#include "dft/functionals/CompositeFunctionals.h"
#include "misc/SerenityError.h"
#include "settings/DFTOptions.h"
/* Include Std and External Headers */
#include <string>

namespace Serenity {

template<>
void DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(
    CompositeFunctionals::XCFUNCTIONALS functional, double& s6, double& rs6, double& s18, double& rs18, double& alp) {
  s6 = 1.0;
  alp = 14.0;
  switch (functional) {
    case CompositeFunctionals::XCFUNCTIONALS::BP86:
      rs6 = 0.3946;
      s18 = 3.2822;
      rs18 = 4.8516;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::BLYP:
      rs6 = 0.4298;
      s18 = 2.6996;
      rs18 = 4.2359;
      break;
      //           case _DF_REVPBE:
      //                rs6 =0.5238;
      //                s18 =2.3550;
      //                rs18=3.5016;
      //           break;
    case CompositeFunctionals::XCFUNCTIONALS::B97_D:
      rs6 = 0.5545;
      s18 = 2.2609;
      rs18 = 3.2297;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::PBE:
      rs6 = 0.4289;
      s18 = 0.7875;
      rs18 = 4.4407;
      break;
      //           case _DF_RPBE:
      //                rs6 =0.1820;
      //                s18 =0.8318;
      //                rs18=4.0094;
      //           break;
      //           case _DF_RPW86PBE:
      //                rs6 =0.4613;
      //                s18 =1.3845;
      //                rs18=4.5062;
      //           break;
    case CompositeFunctionals::XCFUNCTIONALS::B3LYP:
    case CompositeFunctionals::XCFUNCTIONALS::B3LYP_G:
      rs6 = 0.3981;
      s18 = 1.9889;
      rs18 = 4.4211;
      break;
      //               case FUNCTIONALS::BHANDHLYP:
      //                    rs6 =0.2793;
      //                    s18 =1.0354;
      //                    rs18=4.9615;
      //               break;
      //               case FUNCTIONALS::TPSS:
      //                    rs6 =0.4535;
      //                    s18 =1.9435;
      //                    rs18=4.4752;
      //               break;
      //           case _DF_TPSS0:
      //                rs6 =0.3768;
      //                s18 =1.2576;
      //                rs18=4.5865;
      //           break;
    case CompositeFunctionals::XCFUNCTIONALS::PBE0:
      rs6 = 0.4145;
      s18 = 1.2177;
      rs18 = 4.8593;
      break;
      //           case _DF_REVPBE38:
      //                rs6 =0.4309;
      //                s18 =1.4760;
      //                rs18=3.9446;
      //          break;
      //           case _DF_PW6B95:
      //                rs6 =0.2076;
      //                s18 =0.7257;
      //                rs18=6.3750;
      //           break;
    case CompositeFunctionals::XCFUNCTIONALS::B2PLYP:
      rs6 = 0.3065;
      s18 = 0.9147;
      rs18 = 5.0570;
      s6 = 0.64;
      break;
      //           case _DF_mPWLYP:
      //                rs6 =0.4831;
      //                s18 =2.0077;
      //                rs18=4.5323;
      //           break;
    case CompositeFunctionals::XCFUNCTIONALS::OLYP:
      rs6 = 0.5299;
      s18 = 2.6205;
      rs18 = 2.8065;
      break;
      //           case _DF_BPBE:
      //                rs6 =0.4567;
      //                s18 =4.0728;
      //                rs18=4.3908;
      //           break;
      //           case _DF_OPBE:
      //                rs6 =0.5512;
      //                s18 =3.3816;
      //                rs18=2.9444;
      //           break;
      //           case _DF_B3PW91:
      //                rs6 =0.4312;
      //                s18 =2.8524;
      //                rs18=4.4693;
      //           break;
      //           case _DF_REVPBE0:
      //                rs6 =0.4679;
      //                s18 =1.7588;
      //                rs18=3.7619;
      //           break;
      //           case _DF_TPSSh:
      //                rs6 =0.4529;
      //                s18 =2.2382;
      //                rs18=4.6550;
      //           break;
      //           case _DF_PBEH3C:
      //                rs6 =0.4860;
      //                s18 =0.0;
      //                rs18=4.5000;
      //                alp = 14.0;
      //           break;
    case CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP:
      rs6 = 0.3708;
      s18 = 2.0674;
      rs18 = 5.4743;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::B2GPPLYP:
      rs6 = 0.0000;
      s18 = 0.2597;
      rs18 = 6.3332;
      s6 = 0.560;
      break;
    // case CompositeFunctionals::XCFUNCTIONALS::PWPB95:
    // rs6 = 0.0000;
    // s18 = 0.2904;
    // rs18 = 7.3141;
    // s6 = 0.820;
    // break;
    case CompositeFunctionals::XCFUNCTIONALS::DSDBLYP:
      rs6 = 0.0000;
      s18 = 0.2130;
      rs18 = 6.0519;
      s6 = 0.5;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::DSDPBEP86:
      rs6 = 0.0000;
      s18 = 0.0;
      rs18 = 5.65;
      s6 = 0.418;
      break;
    default:
      throw SerenityError("No DFT-D3(BJ) parameters availabale for the given functional.");
      break;
  };
}

template<>
void DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3>(
    CompositeFunctionals::XCFUNCTIONALS functional, double& s6, double& rs6, double& s18, double& rs18, double& alp) {
  s6 = 1.0;
  rs18 = 1.0;
  s18 = 1.0;
  alp = 14.0;
  rs6 = 1.0;
  switch (functional) {
    case CompositeFunctionals::XCFUNCTIONALS::BLYP:
      rs6 = 1.094;
      s18 = 1.682;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::BP86:
      rs6 = 1.139;
      s18 = 1.683;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::B97_D:
      rs6 = 0.892;
      s18 = 0.909;
      break;
      //    case _DF_REVPBE:
      //            rs6=0.923;
      //            s18=1.010;
      //    break;
    case CompositeFunctionals::XCFUNCTIONALS::PBE:
      rs6 = 1.217;
      s18 = 0.722;
      break;
      //    case _DF_RPBE:
      //            rs6=0.872;
      //            s18=0.514;
      //    break;
      //        case FUNCTIONALS::TPSS:
      //                rs6=1.166;
      //                s18=1.105;
      //        break;
    case CompositeFunctionals::XCFUNCTIONALS::B3LYP:
    case CompositeFunctionals::XCFUNCTIONALS::B3LYP_G:
      rs6 = 1.261;
      s18 = 1.703;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::PBE0:
      rs6 = 1.287;
      s18 = 0.928;
      break;
      //    case _DF_PW6B95:
      //            rs6=1.523;
      //            s18=0.862;
      //    break;
      //    case _DF_TPSS0:
      //            rs6=1.252;
      //            s18=1.242;
      //    break;
    case CompositeFunctionals::XCFUNCTIONALS::B2PLYP:
      rs6 = 1.427;
      s18 = 1.022;
      s6 = 0.64;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::B2GPPLYP:
      rs6 = 1.586;
      s18 = 0.760;
      s6 = 0.56;
      break;
      //    case _DF_PWPB95:
      //          rs6=1.557;
      //          s18=0.705;
      //          s6=0.82;
      //    break;
      //    case _DF_mPWLYP:
      //          rs6=1.239;
      //          s18=1.098;
      //    break;
      //    case _DF_BPBE:
      //          rs6=1.087;
      //          s18=2.033;
      //    break;
      //    case _DF_BHANDHLYP: //cbannwarth added missing bhlyp
      //          rs6=1.370;
      //          s18=1.442;
      //    break;
      //    case _DF_TPSSh:
      //          rs6=1.223;
      //          s18=1.219;
      //    break;
      //     case _DF_REVPBE0:
      //          s18=0.792;
      //          rs6=0.949;
      //          s6=1.00;
      //     break;
      //     case _DF_REVPBE38:
      //          s18=0.862;
      //          rs6=1.021;
      //          s6=1.00;
      //     break;
      //     case _DF_RPW86PBE:
      //          s18 =0.901;
      //          rs6 =1.224;
      //          s6=1.0;
      //     break;
      //     case _DF_B3PW91:
      //          s18 =1.775;
      //          rs6 =1.176;
      //          s6=1.0;
      //     break;
      //     case _DF_M06L:
      //          rs6=1.581;
      //          s18=0.000;
      //          s6=1.0;
      //     break;
      //     case _DF_M06:
      //          rs6=1.325;
      //          s18=0.000;
      //          s6=1.0;
      //     break;
    // case CompositeFunctionals::XCFUNCTIONALS::M06_2X:
    //   rs6=1.619;
    //   s18=0.000;
    //   s6=1.0;
    //   break;
    case CompositeFunctionals::XCFUNCTIONALS::WB97X_D:
      s6 = 1.0;
      s18 = 1.0;
      rs6 = 1.281;
      rs18 = 1.094;
      break;
    case CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP:
      s6 = 1.0;
      s18 = 1.217;
      rs6 = 1.378;
      break;
    default:
      throw SerenityError("No DFT-D3(0) parameters availabale for the given functional.");
      break;
  };
  /*
   * Original comment:
  // special TZVPP optimized parameters
  // excerpt from manual
  // "Use special parameters for calculations with triple-zeta basis sets.
  // Preliminary results in the SI of the paper indicate that results are
  // slightly worse than with the default parameters and QZVP type basis
  // sets. This option should be carfully tested for future use in very
  // large computations."
         switch(Functional){
           case _DF_BLYP:
                rs6=1.243;
                s18=2.022;
           break;
           case _DF_BP86:
                rs6=1.221;
                s18=1.838;
           break;
           case _DF_REVPBE:
                rs6=0.953;
                s18=0.989;
           break;
           case _DF_PBE:
                rs6=1.277;
                s18=0.777;
           break;
           case FUNCTIONALS::TPSS:
                rs6=1.213;
                s18=1.176;
           break;
           case _DF_B3LYP:
           case _DF_B3LYP_G:
                rs6=1.314;
                s18=1.706;
           break;
           case _DF_PBE0:
                rs6=1.328;
                s18=0.926;
           break;
           case _DF_PW6B95:
                rs6=1.562;
                s18=0.821;
           break;
           case _DF_TPSS0:
                rs6=1.282;
                s18=1.250;
           break;
           case _DF_B2PLYP:
                rs6=1.551;
                s18=1.109;
                s6=0.5;
           break;
           default:
           break;
         };
   */
}
} /* namespace Serenity */
