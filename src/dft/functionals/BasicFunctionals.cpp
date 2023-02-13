/**
 * @file   BasicFunctionals.cpp
 *
 * @date   Sep 3, 2020
 * @author Jan P. Unsleber
 *
 * IMPORTANT:\n
 * This file was automatically generated please do not alter it.
 * Any required changes should be made to the generating Python script
 * which should be located close by.
 *
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
#include "dft/functionals/BasicFunctionals.h"
#include "misc/SerenityError.h"
#ifdef SERENITY_USE_LIBXC
#include <xc_funcs.h>
#endif /* SERENITY_USE_LIBXC */

namespace Serenity {
namespace BasicFunctionals {

#ifdef SERENITY_USE_LIBXC
int getLibXCAlias(const BASIC_FUNCTIONALS& functional) {
  int alias = -1;
  switch (functional) {
    case BASIC_FUNCTIONALS::X_SLATER:
      alias = XC_LDA_X;
      break;
    case BASIC_FUNCTIONALS::X_EXPONENTIAL_1D:
      alias = XC_LDA_X_1D_EXPONENTIAL;
      break;
    case BASIC_FUNCTIONALS::X_SOFT_1D:
      alias = XC_LDA_X_1D_SOFT;
      break;
    case BASIC_FUNCTIONALS::X_LDA_2D:
      alias = XC_LDA_X_2D;
      break;
    case BASIC_FUNCTIONALS::X_LDA_ERF:
      alias = XC_LDA_X_ERF;
      break;
    case BASIC_FUNCTIONALS::X_RAE:
      alias = XC_LDA_X_RAE;
      break;
    case BASIC_FUNCTIONALS::X_REL:
      alias = XC_LDA_X_REL;
      break;
    case BASIC_FUNCTIONALS::X_SLOC:
      alias = XC_LDA_X_SLOC;
      break;
    case BASIC_FUNCTIONALS::C_PMGB06:
      alias = XC_LDA_C_PMGB06;
      break;
    case BASIC_FUNCTIONALS::C_CSC_1D:
      alias = XC_LDA_C_1D_CSC;
      break;
    case BASIC_FUNCTIONALS::C_LOOS_1D:
      alias = XC_LDA_C_1D_LOOS;
      break;
    case BASIC_FUNCTIONALS::C_AMGB_2D:
      alias = XC_LDA_C_2D_AMGB;
      break;
    case BASIC_FUNCTIONALS::C_PRM_2D:
      alias = XC_LDA_C_2D_PRM;
      break;
    case BASIC_FUNCTIONALS::C_BR78:
      alias = XC_LDA_C_BR78;
      break;
    case BASIC_FUNCTIONALS::C_CHACHIYO:
      alias = XC_LDA_C_CHACHIYO;
      break;
    case BASIC_FUNCTIONALS::C_CHACHIYO_MOD:
      alias = XC_LDA_C_CHACHIYO_MOD;
      break;
    case BASIC_FUNCTIONALS::C_GK72:
      alias = XC_LDA_C_GK72;
      break;
    case BASIC_FUNCTIONALS::C_GL:
      alias = XC_LDA_C_GL;
      break;
    case BASIC_FUNCTIONALS::C_GOMBAS:
      alias = XC_LDA_C_GOMBAS;
      break;
    case BASIC_FUNCTIONALS::C_HL:
      alias = XC_LDA_C_HL;
      break;
    case BASIC_FUNCTIONALS::C_KARASIEV:
      alias = XC_LDA_C_KARASIEV;
      break;
    case BASIC_FUNCTIONALS::C_KARASIEV_MOD:
      alias = XC_LDA_C_KARASIEV_MOD;
      break;
    case BASIC_FUNCTIONALS::C_LP96:
      alias = XC_LDA_C_LP96;
      break;
    case BASIC_FUNCTIONALS::C_MCWEENY:
      alias = XC_LDA_C_MCWEENY;
      break;
    case BASIC_FUNCTIONALS::C_ML1:
      alias = XC_LDA_C_ML1;
      break;
    case BASIC_FUNCTIONALS::C_ML2:
      alias = XC_LDA_C_ML2;
      break;
    case BASIC_FUNCTIONALS::C_OB_PW:
      alias = XC_LDA_C_OB_PW;
      break;
    case BASIC_FUNCTIONALS::C_OB_PZ:
      alias = XC_LDA_C_OB_PZ;
      break;
    case BASIC_FUNCTIONALS::C_OW:
      alias = XC_LDA_C_OW;
      break;
    case BASIC_FUNCTIONALS::C_OW_LYP:
      alias = XC_LDA_C_OW_LYP;
      break;
    case BASIC_FUNCTIONALS::C_PK09:
      alias = XC_LDA_C_PK09;
      break;
    case BASIC_FUNCTIONALS::C_PW:
      alias = XC_LDA_C_PW;
      break;
    case BASIC_FUNCTIONALS::C_PW_MOD:
      alias = XC_LDA_C_PW_MOD;
      break;
    case BASIC_FUNCTIONALS::C_PW_RPA:
      alias = XC_LDA_C_PW_RPA;
      break;
    case BASIC_FUNCTIONALS::C_PZ:
      alias = XC_LDA_C_PZ;
      break;
    case BASIC_FUNCTIONALS::C_PZ_MOD:
      alias = XC_LDA_C_PZ_MOD;
      break;
    case BASIC_FUNCTIONALS::C_RC04:
      alias = XC_LDA_C_RC04;
      break;
    case BASIC_FUNCTIONALS::C_RPA:
      alias = XC_LDA_C_RPA;
      break;
    case BASIC_FUNCTIONALS::C_RPW92:
      alias = XC_LDA_C_RPW92;
      break;
    case BASIC_FUNCTIONALS::C_UPW92:
      alias = XC_LDA_C_UPW92;
      break;
    case BASIC_FUNCTIONALS::C_VBH:
      alias = XC_LDA_C_VBH;
      break;
    case BASIC_FUNCTIONALS::C_VWN:
      alias = XC_LDA_C_VWN;
      break;
    case BASIC_FUNCTIONALS::C_VWN_1:
      alias = XC_LDA_C_VWN_1;
      break;
    case BASIC_FUNCTIONALS::C_VWN_2:
      alias = XC_LDA_C_VWN_2;
      break;
    case BASIC_FUNCTIONALS::C_VWN_3:
      alias = XC_LDA_C_VWN_3;
      break;
    case BASIC_FUNCTIONALS::C_VWN_4:
      alias = XC_LDA_C_VWN_4;
      break;
    case BASIC_FUNCTIONALS::C_VWN_RPA:
      alias = XC_LDA_C_VWN_RPA;
      break;
    case BASIC_FUNCTIONALS::C_WIGNER:
      alias = XC_LDA_C_WIGNER;
      break;
    case BASIC_FUNCTIONALS::C_XALPHA:
      alias = XC_LDA_C_XALPHA;
      break;
    case BASIC_FUNCTIONALS::XC_EHWLRG_1D_1:
      alias = XC_LDA_XC_1D_EHWLRG_1;
      break;
    case BASIC_FUNCTIONALS::XC_EHWLRG_1D_2:
      alias = XC_LDA_XC_1D_EHWLRG_2;
      break;
    case BASIC_FUNCTIONALS::XC_EHWLRG_1D_3:
      alias = XC_LDA_XC_1D_EHWLRG_3;
      break;
    case BASIC_FUNCTIONALS::XC_BN05:
      alias = XC_HYB_LDA_XC_BN05;
      break;
    case BASIC_FUNCTIONALS::XC_GDSMFB:
      alias = XC_LDA_XC_GDSMFB;
      break;
    case BASIC_FUNCTIONALS::XC_KSDT:
      alias = XC_LDA_XC_KSDT;
      break;
    case BASIC_FUNCTIONALS::XC_LP_A:
      alias = XC_LDA_XC_LP_A;
      break;
    case BASIC_FUNCTIONALS::XC_LP_B:
      alias = XC_LDA_XC_LP_B;
      break;
    case BASIC_FUNCTIONALS::XC_TETER93:
      alias = XC_LDA_XC_TETER93;
      break;
    case BASIC_FUNCTIONALS::XC_TIH:
      alias = XC_LDA_XC_TIH;
      break;
    case BASIC_FUNCTIONALS::XC_ZLP:
      alias = XC_LDA_XC_ZLP;
      break;
    case BASIC_FUNCTIONALS::K_LP:
      alias = XC_LDA_K_LP;
      break;
    case BASIC_FUNCTIONALS::K_LP96_K:
      alias = XC_LDA_K_LP96;
      break;
    case BASIC_FUNCTIONALS::K_TF:
      alias = XC_LDA_K_TF;
      break;
    case BASIC_FUNCTIONALS::K_ZLP_K:
      alias = XC_LDA_K_ZLP;
      break;
    case BASIC_FUNCTIONALS::C_MGGAC:
      alias = XC_GGA_C_MGGAC;
      break;
    case BASIC_FUNCTIONALS::X_B86_2D:
      alias = XC_GGA_X_2D_B86;
      break;
    case BASIC_FUNCTIONALS::X_B86_MGC_2D:
      alias = XC_GGA_X_2D_B86_MGC;
      break;
    case BASIC_FUNCTIONALS::X_B88_2D:
      alias = XC_GGA_X_2D_B88;
      break;
    case BASIC_FUNCTIONALS::X_PBE_2D:
      alias = XC_GGA_X_2D_PBE;
      break;
    case BASIC_FUNCTIONALS::X_AIRY:
      alias = XC_GGA_X_AIRY;
      break;
    case BASIC_FUNCTIONALS::X_AK13:
      alias = XC_GGA_X_AK13;
      break;
    case BASIC_FUNCTIONALS::X_AM05:
      alias = XC_GGA_X_AM05;
      break;
    case BASIC_FUNCTIONALS::X_APBE:
      alias = XC_GGA_X_APBE;
      break;
    case BASIC_FUNCTIONALS::X_B86:
      alias = XC_GGA_X_B86;
      break;
    case BASIC_FUNCTIONALS::X_B86_MGC:
      alias = XC_GGA_X_B86_MGC;
      break;
    case BASIC_FUNCTIONALS::X_B86_R:
      alias = XC_GGA_X_B86_R;
      break;
    case BASIC_FUNCTIONALS::X_B88:
      alias = XC_GGA_X_B88;
      break;
    case BASIC_FUNCTIONALS::X_B88_6311G:
      alias = XC_GGA_X_B88_6311G;
      break;
    case BASIC_FUNCTIONALS::X_B88M:
      alias = XC_GGA_X_B88M;
      break;
    case BASIC_FUNCTIONALS::X_BAYESIAN:
      alias = XC_GGA_X_BAYESIAN;
      break;
    case BASIC_FUNCTIONALS::X_BCGP:
      alias = XC_GGA_X_BCGP;
      break;
    case BASIC_FUNCTIONALS::X_BEEFVDW:
      alias = XC_GGA_X_BEEFVDW;
      break;
    case BASIC_FUNCTIONALS::X_BPCCAC:
      alias = XC_GGA_X_BPCCAC;
      break;
    case BASIC_FUNCTIONALS::X_C09X:
      alias = XC_GGA_X_C09X;
      break;
    case BASIC_FUNCTIONALS::X_CAP:
      alias = XC_GGA_X_CAP;
      break;
    case BASIC_FUNCTIONALS::X_CHACHIYO:
      alias = XC_GGA_X_CHACHIYO;
      break;
    case BASIC_FUNCTIONALS::X_DK87_R1:
      alias = XC_GGA_X_DK87_R1;
      break;
    case BASIC_FUNCTIONALS::X_DK87_R2:
      alias = XC_GGA_X_DK87_R2;
      break;
    case BASIC_FUNCTIONALS::X_EB88:
      alias = XC_GGA_X_EB88;
      break;
    case BASIC_FUNCTIONALS::X_ECMV92:
      alias = XC_GGA_X_ECMV92;
      break;
    case BASIC_FUNCTIONALS::X_EV93:
      alias = XC_GGA_X_EV93;
      break;
    case BASIC_FUNCTIONALS::X_FD_LB94:
      alias = XC_GGA_X_FD_LB94;
      break;
    case BASIC_FUNCTIONALS::X_FD_REVLB94:
      alias = XC_GGA_X_FD_REVLB94;
      break;
    case BASIC_FUNCTIONALS::X_FT97_A:
      alias = XC_GGA_X_FT97_A;
      break;
    case BASIC_FUNCTIONALS::X_FT97_B:
      alias = XC_GGA_X_FT97_B;
      break;
    case BASIC_FUNCTIONALS::X_G96:
      alias = XC_GGA_X_G96;
      break;
    case BASIC_FUNCTIONALS::X_GAM:
      alias = XC_GGA_X_GAM;
      break;
    case BASIC_FUNCTIONALS::X_GG99:
      alias = XC_GGA_X_GG99;
      break;
    case BASIC_FUNCTIONALS::X_HCTH_A:
      alias = XC_GGA_X_HCTH_A;
      break;
    case BASIC_FUNCTIONALS::X_HERMAN:
      alias = XC_GGA_X_HERMAN;
      break;
    case BASIC_FUNCTIONALS::X_HJS_B88:
      alias = XC_GGA_X_HJS_B88;
      break;
    case BASIC_FUNCTIONALS::X_HJS_B88_V2:
      alias = XC_GGA_X_HJS_B88_V2;
      break;
    case BASIC_FUNCTIONALS::X_HJS_B97X:
      alias = XC_GGA_X_HJS_B97X;
      break;
    case BASIC_FUNCTIONALS::X_HJS_PBE:
      alias = XC_GGA_X_HJS_PBE;
      break;
    case BASIC_FUNCTIONALS::X_HJS_PBE_SOL:
      alias = XC_GGA_X_HJS_PBE_SOL;
      break;
    case BASIC_FUNCTIONALS::X_HTBS:
      alias = XC_GGA_X_HTBS;
      break;
    case BASIC_FUNCTIONALS::X_ITYH:
      alias = XC_GGA_X_ITYH;
      break;
    case BASIC_FUNCTIONALS::X_KGG99:
      alias = XC_GGA_X_KGG99;
      break;
    case BASIC_FUNCTIONALS::X_KT1:
      alias = XC_GGA_X_KT1;
      break;
    case BASIC_FUNCTIONALS::X_LAG:
      alias = XC_GGA_X_LAG;
      break;
    case BASIC_FUNCTIONALS::X_LAMBDA_CH_N:
      alias = XC_GGA_X_LAMBDA_CH_N;
      break;
    case BASIC_FUNCTIONALS::X_LAMBDA_LO_N:
      alias = XC_GGA_X_LAMBDA_LO_N;
      break;
    case BASIC_FUNCTIONALS::X_LAMBDA_OC2_N:
      alias = XC_GGA_X_LAMBDA_OC2_N;
      break;
    case BASIC_FUNCTIONALS::X_LB:
      alias = XC_GGA_X_LB;
      break;
    case BASIC_FUNCTIONALS::X_LBM:
      alias = XC_GGA_X_LBM;
      break;
    case BASIC_FUNCTIONALS::X_LG93:
      alias = XC_GGA_X_LG93;
      break;
    case BASIC_FUNCTIONALS::X_LSPBE:
      alias = XC_GGA_X_LSPBE;
      break;
    case BASIC_FUNCTIONALS::X_LSRPBE:
      alias = XC_GGA_X_LSRPBE;
      break;
    case BASIC_FUNCTIONALS::X_LV_RPW86:
      alias = XC_GGA_X_LV_RPW86;
      break;
    case BASIC_FUNCTIONALS::X_MB88:
      alias = XC_GGA_X_MB88;
      break;
    case BASIC_FUNCTIONALS::X_MPBE:
      alias = XC_GGA_X_MPBE;
      break;
    case BASIC_FUNCTIONALS::X_MPW91:
      alias = XC_GGA_X_MPW91;
      break;
    case BASIC_FUNCTIONALS::X_N12:
      alias = XC_GGA_X_N12;
      break;
    case BASIC_FUNCTIONALS::X_NCAP:
      alias = XC_GGA_X_NCAP;
      break;
    case BASIC_FUNCTIONALS::X_OL2:
      alias = XC_GGA_X_OL2;
      break;
    case BASIC_FUNCTIONALS::X_OPTB86B_VDW:
      alias = XC_GGA_X_OPTB86B_VDW;
      break;
    case BASIC_FUNCTIONALS::X_OPTB88_VDW:
      alias = XC_GGA_X_OPTB88_VDW;
      break;
    case BASIC_FUNCTIONALS::X_OPTPBE_VDW:
      alias = XC_GGA_X_OPTPBE_VDW;
      break;
    case BASIC_FUNCTIONALS::X_OPTX:
      alias = XC_GGA_X_OPTX;
      break;
    case BASIC_FUNCTIONALS::X_PBE:
      alias = XC_GGA_X_PBE;
      break;
    case BASIC_FUNCTIONALS::X_PBE_JSJR:
      alias = XC_GGA_X_PBE_JSJR;
      break;
    case BASIC_FUNCTIONALS::X_PBE_MOL:
      alias = XC_GGA_X_PBE_MOL;
      break;
    case BASIC_FUNCTIONALS::X_PBE_R:
      alias = XC_GGA_X_PBE_R;
      break;
    case BASIC_FUNCTIONALS::X_PBE_SOL:
      alias = XC_GGA_X_PBE_SOL;
      break;
    case BASIC_FUNCTIONALS::X_PBE_TCA:
      alias = XC_GGA_X_PBE_TCA;
      break;
    case BASIC_FUNCTIONALS::X_PBEA:
      alias = XC_GGA_X_PBEA;
      break;
    case BASIC_FUNCTIONALS::X_PBEFE:
      alias = XC_GGA_X_PBEFE;
      break;
    case BASIC_FUNCTIONALS::X_PBEINT:
      alias = XC_GGA_X_PBEINT;
      break;
    case BASIC_FUNCTIONALS::X_PBEK1_VDW:
      alias = XC_GGA_X_PBEK1_VDW;
      break;
    case BASIC_FUNCTIONALS::X_PBEPOW:
      alias = XC_GGA_X_PBEPOW;
      break;
    case BASIC_FUNCTIONALS::X_PBETRANS:
      alias = XC_GGA_X_PBETRANS;
      break;
    case BASIC_FUNCTIONALS::X_PW86:
      alias = XC_GGA_X_PW86;
      break;
    case BASIC_FUNCTIONALS::X_PW91:
      alias = XC_GGA_X_PW91;
      break;
    case BASIC_FUNCTIONALS::X_Q2D:
      alias = XC_GGA_X_Q2D;
      break;
    case BASIC_FUNCTIONALS::X_RGE2:
      alias = XC_GGA_X_RGE2;
      break;
    case BASIC_FUNCTIONALS::X_RPBE:
      alias = XC_GGA_X_RPBE;
      break;
    case BASIC_FUNCTIONALS::X_RPW86:
      alias = XC_GGA_X_RPW86;
      break;
    case BASIC_FUNCTIONALS::X_S12G:
      alias = XC_GGA_X_S12G;
      break;
    case BASIC_FUNCTIONALS::X_SFAT:
      alias = XC_GGA_X_SFAT;
      break;
    case BASIC_FUNCTIONALS::X_SFAT_PBE:
      alias = XC_GGA_X_SFAT_PBE;
      break;
    case BASIC_FUNCTIONALS::X_SG4:
      alias = XC_GGA_X_SG4;
      break;
    case BASIC_FUNCTIONALS::X_SOGGA:
      alias = XC_GGA_X_SOGGA;
      break;
    case BASIC_FUNCTIONALS::X_SOGGA11:
      alias = XC_GGA_X_SOGGA11;
      break;
    case BASIC_FUNCTIONALS::X_SSB:
      alias = XC_GGA_X_SSB;
      break;
    case BASIC_FUNCTIONALS::X_SSB_D:
      alias = XC_GGA_X_SSB_D;
      break;
    case BASIC_FUNCTIONALS::X_SSB_SW:
      alias = XC_GGA_X_SSB_SW;
      break;
    case BASIC_FUNCTIONALS::X_VMT84_GE:
      alias = XC_GGA_X_VMT84_GE;
      break;
    case BASIC_FUNCTIONALS::X_VMT84_PBE:
      alias = XC_GGA_X_VMT84_PBE;
      break;
    case BASIC_FUNCTIONALS::X_VMT_GE:
      alias = XC_GGA_X_VMT_GE;
      break;
    case BASIC_FUNCTIONALS::X_VMT_PBE:
      alias = XC_GGA_X_VMT_PBE;
      break;
    case BASIC_FUNCTIONALS::X_WC:
      alias = XC_GGA_X_WC;
      break;
    case BASIC_FUNCTIONALS::X_WPBEH:
      alias = XC_GGA_X_WPBEH;
      break;
    case BASIC_FUNCTIONALS::X_XPBE:
      alias = XC_GGA_X_XPBE;
      break;
    case BASIC_FUNCTIONALS::C_ACGGA:
      alias = XC_GGA_C_ACGGA;
      break;
    case BASIC_FUNCTIONALS::C_ACGGAP:
      alias = XC_GGA_C_ACGGAP;
      break;
    case BASIC_FUNCTIONALS::C_AM05:
      alias = XC_GGA_C_AM05;
      break;
    case BASIC_FUNCTIONALS::C_APBE:
      alias = XC_GGA_C_APBE;
      break;
    case BASIC_FUNCTIONALS::C_BMK:
      alias = XC_GGA_C_BMK;
      break;
    case BASIC_FUNCTIONALS::C_CHACHIYO_GGA:
      alias = XC_GGA_C_CHACHIYO;
      break;
    case BASIC_FUNCTIONALS::C_CS1:
      alias = XC_GGA_C_CS1;
      break;
    case BASIC_FUNCTIONALS::C_FT97:
      alias = XC_GGA_C_FT97;
      break;
    case BASIC_FUNCTIONALS::C_GAM:
      alias = XC_GGA_C_GAM;
      break;
    case BASIC_FUNCTIONALS::C_GAPC:
      alias = XC_GGA_C_GAPC;
      break;
    case BASIC_FUNCTIONALS::C_GAPLOC:
      alias = XC_GGA_C_GAPLOC;
      break;
    case BASIC_FUNCTIONALS::C_HCTH_A:
      alias = XC_GGA_C_HCTH_A;
      break;
    case BASIC_FUNCTIONALS::C_HYB_TAU_HCTH:
      alias = XC_GGA_C_HYB_TAU_HCTH;
      break;
    case BASIC_FUNCTIONALS::C_LM:
      alias = XC_GGA_C_LM;
      break;
    case BASIC_FUNCTIONALS::C_LYP:
      alias = XC_GGA_C_LYP;
      break;
    case BASIC_FUNCTIONALS::C_N12:
      alias = XC_GGA_C_N12;
      break;
    case BASIC_FUNCTIONALS::C_N12_SX:
      alias = XC_GGA_C_N12_SX;
      break;
    case BASIC_FUNCTIONALS::C_OP_B88:
      alias = XC_GGA_C_OP_B88;
      break;
    case BASIC_FUNCTIONALS::C_OP_G96:
      alias = XC_GGA_C_OP_G96;
      break;
    case BASIC_FUNCTIONALS::C_OP_PBE:
      alias = XC_GGA_C_OP_PBE;
      break;
    case BASIC_FUNCTIONALS::C_OP_PW91:
      alias = XC_GGA_C_OP_PW91;
      break;
    case BASIC_FUNCTIONALS::C_OP_XALPHA:
      alias = XC_GGA_C_OP_XALPHA;
      break;
    case BASIC_FUNCTIONALS::C_OPTC:
      alias = XC_GGA_C_OPTC;
      break;
    case BASIC_FUNCTIONALS::C_P86:
      alias = XC_GGA_C_P86;
      break;
    case BASIC_FUNCTIONALS::C_PBE:
      alias = XC_GGA_C_PBE;
      break;
    case BASIC_FUNCTIONALS::C_PBE_JRGX:
      alias = XC_GGA_C_PBE_JRGX;
      break;
    case BASIC_FUNCTIONALS::C_PBE_MOL:
      alias = XC_GGA_C_PBE_MOL;
      break;
    case BASIC_FUNCTIONALS::C_PBE_SOL:
      alias = XC_GGA_C_PBE_SOL;
      break;
    case BASIC_FUNCTIONALS::C_PBE_VWN:
      alias = XC_GGA_C_PBE_VWN;
      break;
    case BASIC_FUNCTIONALS::C_PBEFE:
      alias = XC_GGA_C_PBEFE;
      break;
    case BASIC_FUNCTIONALS::C_PBEINT:
      alias = XC_GGA_C_PBEINT;
      break;
    case BASIC_FUNCTIONALS::C_PBELOC:
      alias = XC_GGA_C_PBELOC;
      break;
    case BASIC_FUNCTIONALS::C_PW91:
      alias = XC_GGA_C_PW91;
      break;
    case BASIC_FUNCTIONALS::C_Q2D:
      alias = XC_GGA_C_Q2D;
      break;
    case BASIC_FUNCTIONALS::C_REGTPSS:
      alias = XC_GGA_C_REGTPSS;
      break;
    case BASIC_FUNCTIONALS::C_REVTCA:
      alias = XC_GGA_C_REVTCA;
      break;
    case BASIC_FUNCTIONALS::C_RGE2:
      alias = XC_GGA_C_RGE2;
      break;
    case BASIC_FUNCTIONALS::C_SCAN_E0:
      alias = XC_GGA_C_SCAN_E0;
      break;
    case BASIC_FUNCTIONALS::C_SG4:
      alias = XC_GGA_C_SG4;
      break;
    case BASIC_FUNCTIONALS::C_SOGGA11:
      alias = XC_GGA_C_SOGGA11;
      break;
    case BASIC_FUNCTIONALS::C_SOGGA11_X:
      alias = XC_GGA_C_SOGGA11_X;
      break;
    case BASIC_FUNCTIONALS::C_SPBE:
      alias = XC_GGA_C_SPBE;
      break;
    case BASIC_FUNCTIONALS::C_TAU_HCTH:
      alias = XC_GGA_C_TAU_HCTH;
      break;
    case BASIC_FUNCTIONALS::C_TCA:
      alias = XC_GGA_C_TCA;
      break;
    case BASIC_FUNCTIONALS::C_TM_LYP:
      alias = XC_GGA_C_TM_LYP;
      break;
    case BASIC_FUNCTIONALS::C_TM_PBE:
      alias = XC_GGA_C_TM_PBE;
      break;
    case BASIC_FUNCTIONALS::C_W94:
      alias = XC_GGA_C_W94;
      break;
    case BASIC_FUNCTIONALS::C_WI:
      alias = XC_GGA_C_WI;
      break;
    case BASIC_FUNCTIONALS::C_WI0:
      alias = XC_GGA_C_WI0;
      break;
    case BASIC_FUNCTIONALS::C_WL:
      alias = XC_GGA_C_WL;
      break;
    case BASIC_FUNCTIONALS::C_XPBE:
      alias = XC_GGA_C_XPBE;
      break;
    case BASIC_FUNCTIONALS::C_ZPBEINT:
      alias = XC_GGA_C_ZPBEINT;
      break;
    case BASIC_FUNCTIONALS::C_ZPBESOL:
      alias = XC_GGA_C_ZPBESOL;
      break;
    case BASIC_FUNCTIONALS::C_ZVPBEINT:
      alias = XC_GGA_C_ZVPBEINT;
      break;
    case BASIC_FUNCTIONALS::C_ZVPBELOC:
      alias = XC_GGA_C_ZVPBELOC;
      break;
    case BASIC_FUNCTIONALS::C_ZVPBESOL:
      alias = XC_GGA_C_ZVPBESOL;
      break;
    case BASIC_FUNCTIONALS::XC_B97_D:
      alias = XC_GGA_XC_B97_D;
      break;
    case BASIC_FUNCTIONALS::XC_B97_GGA1:
      alias = XC_GGA_XC_B97_GGA1;
      break;
    case BASIC_FUNCTIONALS::XC_BEEFVDW:
      alias = XC_GGA_XC_BEEFVDW;
      break;
    case BASIC_FUNCTIONALS::XC_EDF1:
      alias = XC_GGA_XC_EDF1;
      break;
    case BASIC_FUNCTIONALS::XC_HCTH_120:
      alias = XC_GGA_XC_HCTH_120;
      break;
    case BASIC_FUNCTIONALS::XC_HCTH_147:
      alias = XC_GGA_XC_HCTH_147;
      break;
    case BASIC_FUNCTIONALS::XC_HCTH_407:
      alias = XC_GGA_XC_HCTH_407;
      break;
    case BASIC_FUNCTIONALS::XC_HCTH_407P:
      alias = XC_GGA_XC_HCTH_407P;
      break;
    case BASIC_FUNCTIONALS::XC_HCTH_93:
      alias = XC_GGA_XC_HCTH_93;
      break;
    case BASIC_FUNCTIONALS::XC_HCTH_P14:
      alias = XC_GGA_XC_HCTH_P14;
      break;
    case BASIC_FUNCTIONALS::XC_HCTH_P76:
      alias = XC_GGA_XC_HCTH_P76;
      break;
    case BASIC_FUNCTIONALS::XC_HLE16:
      alias = XC_GGA_XC_HLE16;
      break;
    case BASIC_FUNCTIONALS::XC_KT1:
      alias = XC_GGA_XC_KT1;
      break;
    case BASIC_FUNCTIONALS::XC_KT2:
      alias = XC_GGA_XC_KT2;
      break;
    case BASIC_FUNCTIONALS::XC_KT3:
      alias = XC_GGA_XC_KT3;
      break;
    case BASIC_FUNCTIONALS::XC_LB07:
      alias = XC_HYB_GGA_XC_LB07;
      break;
    case BASIC_FUNCTIONALS::XC_MOHLYP:
      alias = XC_GGA_XC_MOHLYP;
      break;
    case BASIC_FUNCTIONALS::XC_MOHLYP2:
      alias = XC_GGA_XC_MOHLYP2;
      break;
    case BASIC_FUNCTIONALS::XC_MPWLYP1W:
      alias = XC_GGA_XC_MPWLYP1W;
      break;
    case BASIC_FUNCTIONALS::XC_NCAP:
      alias = XC_GGA_XC_NCAP;
      break;
    case BASIC_FUNCTIONALS::XC_OBLYP_D:
      alias = XC_GGA_XC_OBLYP_D;
      break;
    case BASIC_FUNCTIONALS::XC_OPBE_D:
      alias = XC_GGA_XC_OPBE_D;
      break;
    case BASIC_FUNCTIONALS::XC_OPWLYP_D:
      alias = XC_GGA_XC_OPWLYP_D;
      break;
    case BASIC_FUNCTIONALS::XC_PBE1W:
      alias = XC_GGA_XC_PBE1W;
      break;
    case BASIC_FUNCTIONALS::XC_PBELYP1W:
      alias = XC_GGA_XC_PBELYP1W;
      break;
    case BASIC_FUNCTIONALS::XC_TH1:
      alias = XC_GGA_XC_TH1;
      break;
    case BASIC_FUNCTIONALS::XC_TH2:
      alias = XC_GGA_XC_TH2;
      break;
    case BASIC_FUNCTIONALS::XC_TH3:
      alias = XC_GGA_XC_TH3;
      break;
    case BASIC_FUNCTIONALS::XC_TH4:
      alias = XC_GGA_XC_TH4;
      break;
    case BASIC_FUNCTIONALS::XC_TH_FC:
      alias = XC_GGA_XC_TH_FC;
      break;
    case BASIC_FUNCTIONALS::XC_TH_FCFO:
      alias = XC_GGA_XC_TH_FCFO;
      break;
    case BASIC_FUNCTIONALS::XC_TH_FCO:
      alias = XC_GGA_XC_TH_FCO;
      break;
    case BASIC_FUNCTIONALS::XC_TH_FL:
      alias = XC_GGA_XC_TH_FL;
      break;
    case BASIC_FUNCTIONALS::XC_VV10:
      alias = XC_GGA_XC_VV10;
      break;
    case BASIC_FUNCTIONALS::XC_XLYP:
      alias = XC_GGA_XC_XLYP;
      break;
    case BASIC_FUNCTIONALS::K_ABSP1:
      alias = XC_GGA_K_ABSP1;
      break;
    case BASIC_FUNCTIONALS::K_ABSP2:
      alias = XC_GGA_K_ABSP2;
      break;
    case BASIC_FUNCTIONALS::K_ABSP3:
      alias = XC_GGA_K_ABSP3;
      break;
    case BASIC_FUNCTIONALS::K_ABSP4:
      alias = XC_GGA_K_ABSP4;
      break;
    case BASIC_FUNCTIONALS::K_APBE:
      alias = XC_GGA_K_APBE;
      break;
    case BASIC_FUNCTIONALS::K_APBEINT:
      alias = XC_GGA_K_APBEINT;
      break;
    case BASIC_FUNCTIONALS::K_BALTIN:
      alias = XC_GGA_K_BALTIN;
      break;
    case BASIC_FUNCTIONALS::K_DK:
      alias = XC_GGA_K_DK;
      break;
    case BASIC_FUNCTIONALS::K_ERNZERHOF:
      alias = XC_GGA_K_ERNZERHOF;
      break;
    case BASIC_FUNCTIONALS::K_EXP4:
      alias = XC_GGA_K_EXP4;
      break;
    case BASIC_FUNCTIONALS::K_FR_B88:
      alias = XC_GGA_K_FR_B88;
      break;
    case BASIC_FUNCTIONALS::K_FR_PW86:
      alias = XC_GGA_K_FR_PW86;
      break;
    case BASIC_FUNCTIONALS::K_GDS08:
      alias = XC_GGA_K_GDS08;
      break;
    case BASIC_FUNCTIONALS::K_GE2:
      alias = XC_GGA_K_GE2;
      break;
    case BASIC_FUNCTIONALS::K_GHDS10:
      alias = XC_GGA_K_GHDS10;
      break;
    case BASIC_FUNCTIONALS::K_GHDS10R:
      alias = XC_GGA_K_GHDS10R;
      break;
    case BASIC_FUNCTIONALS::K_GOLDEN:
      alias = XC_GGA_K_GOLDEN;
      break;
    case BASIC_FUNCTIONALS::K_GP85:
      alias = XC_GGA_K_GP85;
      break;
    case BASIC_FUNCTIONALS::K_GR:
      alias = XC_GGA_K_GR;
      break;
    case BASIC_FUNCTIONALS::K_PW91:
      alias = XC_GGA_K_LC94;
      break;
    case BASIC_FUNCTIONALS::K_LIEB:
      alias = XC_GGA_K_LIEB;
      break;
    case BASIC_FUNCTIONALS::K_LUDENA:
      alias = XC_GGA_K_LUDENA;
      break;
    case BASIC_FUNCTIONALS::K_LLP:
      alias = XC_GGA_K_LLP;
      break;
    case BASIC_FUNCTIONALS::K_MEYER:
      alias = XC_GGA_K_MEYER;
      break;
    case BASIC_FUNCTIONALS::K_OL1:
      alias = XC_GGA_K_OL1;
      break;
    case BASIC_FUNCTIONALS::K_OL2:
      alias = XC_GGA_K_OL2;
      break;
    case BASIC_FUNCTIONALS::K_PBE3:
      alias = XC_GGA_K_PBE3;
      break;
    case BASIC_FUNCTIONALS::K_PBE4:
      alias = XC_GGA_K_PBE4;
      break;
    case BASIC_FUNCTIONALS::K_PEARSON:
      alias = XC_GGA_K_PEARSON;
      break;
    case BASIC_FUNCTIONALS::K_PERDEW:
      alias = XC_GGA_K_PERDEW;
      break;
    case BASIC_FUNCTIONALS::K_REVAPBE:
      alias = XC_GGA_K_REVAPBE;
      break;
    case BASIC_FUNCTIONALS::K_REVAPBEINT:
      alias = XC_GGA_K_REVAPBEINT;
      break;
    case BASIC_FUNCTIONALS::K_TFVW:
      alias = XC_GGA_K_TFVW;
      break;
    case BASIC_FUNCTIONALS::K_THAKKAR:
      alias = XC_GGA_K_THAKKAR;
      break;
    case BASIC_FUNCTIONALS::K_TKVLN:
      alias = XC_GGA_K_TKVLN;
      break;
    case BASIC_FUNCTIONALS::K_TW1:
      alias = XC_GGA_K_TW1;
      break;
    case BASIC_FUNCTIONALS::K_TW2:
      alias = XC_GGA_K_TW2;
      break;
    case BASIC_FUNCTIONALS::K_TW3:
      alias = XC_GGA_K_TW3;
      break;
    case BASIC_FUNCTIONALS::K_TW4:
      alias = XC_GGA_K_TW4;
      break;
    case BASIC_FUNCTIONALS::K_VJKS:
      alias = XC_GGA_K_VJKS;
      break;
    case BASIC_FUNCTIONALS::K_VSK:
      alias = XC_GGA_K_VSK;
      break;
    case BASIC_FUNCTIONALS::K_VW:
      alias = XC_GGA_K_VW;
      break;
    case BASIC_FUNCTIONALS::K_YT65:
      alias = XC_GGA_K_YT65;
      break;
    case BASIC_FUNCTIONALS::X_N12_SX:
      alias = XC_HYB_GGA_X_N12_SX;
      break;
    case BASIC_FUNCTIONALS::X_S12H:
      alias = XC_HYB_GGA_X_S12H;
      break;
    case BASIC_FUNCTIONALS::X_SOGGA11_X:
      alias = XC_HYB_GGA_X_SOGGA11_X;
      break;
    case BASIC_FUNCTIONALS::XC_APBE0:
      alias = XC_HYB_GGA_XC_APBE0;
      break;
    case BASIC_FUNCTIONALS::XC_APF:
      alias = XC_HYB_GGA_XC_APF;
      break;
    case BASIC_FUNCTIONALS::XC_B1LYP:
      alias = XC_HYB_GGA_XC_B1LYP;
      break;
    case BASIC_FUNCTIONALS::XC_B1PW91:
      alias = XC_HYB_GGA_XC_B1PW91;
      break;
    case BASIC_FUNCTIONALS::XC_B1WC:
      alias = XC_HYB_GGA_XC_B1WC;
      break;
    case BASIC_FUNCTIONALS::XC_B3LYP:
      alias = XC_HYB_GGA_XC_B3LYP;
      break;
    case BASIC_FUNCTIONALS::XC_B3LYP5:
      alias = XC_HYB_GGA_XC_B3LYP5;
      break;
    case BASIC_FUNCTIONALS::XC_B3LYP_MCM1:
      alias = XC_HYB_GGA_XC_B3LYP_MCM1;
      break;
    case BASIC_FUNCTIONALS::XC_B3LYP_MCM2:
      alias = XC_HYB_GGA_XC_B3LYP_MCM2;
      break;
    case BASIC_FUNCTIONALS::XC_B3LYPS:
      alias = XC_HYB_GGA_XC_B3LYPS;
      break;
    case BASIC_FUNCTIONALS::XC_B3P86:
      alias = XC_HYB_GGA_XC_B3P86;
      break;
    case BASIC_FUNCTIONALS::XC_B3PW91:
      alias = XC_HYB_GGA_XC_B3PW91;
      break;
    case BASIC_FUNCTIONALS::XC_B5050LYP:
      alias = XC_HYB_GGA_XC_B5050LYP;
      break;
    case BASIC_FUNCTIONALS::XC_B97:
      alias = XC_HYB_GGA_XC_B97;
      break;
    case BASIC_FUNCTIONALS::XC_B97_1:
      alias = XC_HYB_GGA_XC_B97_1;
      break;
    case BASIC_FUNCTIONALS::XC_B97_1P:
      alias = XC_HYB_GGA_XC_B97_1P;
      break;
    case BASIC_FUNCTIONALS::XC_B97_2:
      alias = XC_HYB_GGA_XC_B97_2;
      break;
    case BASIC_FUNCTIONALS::XC_B97_3:
      alias = XC_HYB_GGA_XC_B97_3;
      break;
    case BASIC_FUNCTIONALS::XC_B97_K:
      alias = XC_HYB_GGA_XC_B97_K;
      break;
    case BASIC_FUNCTIONALS::XC_BHANDH:
      alias = XC_HYB_GGA_XC_BHANDH;
      break;
    case BASIC_FUNCTIONALS::XC_BHANDHLYP:
      alias = XC_HYB_GGA_XC_BHANDHLYP;
      break;
    case BASIC_FUNCTIONALS::XC_BLYP35:
      alias = XC_HYB_GGA_XC_BLYP35;
      break;
    case BASIC_FUNCTIONALS::XC_CAM_B3LYP:
      alias = XC_HYB_GGA_XC_CAM_B3LYP;
      break;
    case BASIC_FUNCTIONALS::XC_CAM_PBEH:
      alias = XC_HYB_GGA_XC_CAM_PBEH;
      break;
    case BASIC_FUNCTIONALS::XC_CAM_QTP_00:
      alias = XC_HYB_GGA_XC_CAM_QTP_00;
      break;
    case BASIC_FUNCTIONALS::XC_CAM_QTP_01:
      alias = XC_HYB_GGA_XC_CAM_QTP_01;
      break;
    case BASIC_FUNCTIONALS::XC_CAM_QTP_02:
      alias = XC_HYB_GGA_XC_CAM_QTP_02;
      break;
    case BASIC_FUNCTIONALS::XC_CAMY_B3LYP:
      alias = XC_HYB_GGA_XC_CAMY_B3LYP;
      break;
    case BASIC_FUNCTIONALS::XC_CAMY_BLYP:
      alias = XC_HYB_GGA_XC_CAMY_BLYP;
      break;
    case BASIC_FUNCTIONALS::XC_CAMY_PBEH:
      alias = XC_HYB_GGA_XC_CAMY_PBEH;
      break;
    case BASIC_FUNCTIONALS::XC_CAP0:
      alias = XC_HYB_GGA_XC_CAP0;
      break;
    case BASIC_FUNCTIONALS::XC_EDF2:
      alias = XC_HYB_GGA_XC_EDF2;
      break;
    case BASIC_FUNCTIONALS::XC_HAPBE:
      alias = XC_HYB_GGA_XC_HAPBE;
      break;
    case BASIC_FUNCTIONALS::XC_HJS_B88:
      alias = XC_HYB_GGA_XC_HJS_B88;
      break;
    case BASIC_FUNCTIONALS::XC_HJS_B97X:
      alias = XC_HYB_GGA_XC_HJS_B97X;
      break;
    case BASIC_FUNCTIONALS::XC_HJS_PBE:
      alias = XC_HYB_GGA_XC_HJS_PBE;
      break;
    case BASIC_FUNCTIONALS::XC_HJS_PBE_SOL:
      alias = XC_HYB_GGA_XC_HJS_PBE_SOL;
      break;
    case BASIC_FUNCTIONALS::XC_HPBEINT:
      alias = XC_HYB_GGA_XC_HPBEINT;
      break;
    case BASIC_FUNCTIONALS::XC_HSE03:
      alias = XC_HYB_GGA_XC_HSE03;
      break;
    case BASIC_FUNCTIONALS::XC_HSE06:
      alias = XC_HYB_GGA_XC_HSE06;
      break;
    case BASIC_FUNCTIONALS::XC_HSE12:
      alias = XC_HYB_GGA_XC_HSE12;
      break;
    case BASIC_FUNCTIONALS::XC_HSE12S:
      alias = XC_HYB_GGA_XC_HSE12S;
      break;
    case BASIC_FUNCTIONALS::XC_HSE_SOL:
      alias = XC_HYB_GGA_XC_HSE_SOL;
      break;
    case BASIC_FUNCTIONALS::XC_KMLYP:
      alias = XC_HYB_GGA_XC_KMLYP;
      break;
    case BASIC_FUNCTIONALS::XC_LC_BLYP:
      alias = XC_HYB_GGA_XC_LC_BLYP;
      break;
    case BASIC_FUNCTIONALS::XC_LC_QTP:
      alias = XC_HYB_GGA_XC_LC_QTP;
      break;
    case BASIC_FUNCTIONALS::XC_LC_VV10:
      alias = XC_HYB_GGA_XC_LC_VV10;
      break;
    case BASIC_FUNCTIONALS::XC_LC_WPBE:
      alias = XC_HYB_GGA_XC_LC_WPBE;
      break;
    case BASIC_FUNCTIONALS::XC_LC_WPBE08_WHS:
      alias = XC_HYB_GGA_XC_LC_WPBE08_WHS;
      break;
    case BASIC_FUNCTIONALS::XC_LC_WPBE_WHS:
      alias = XC_HYB_GGA_XC_LC_WPBE_WHS;
      break;
    case BASIC_FUNCTIONALS::XC_LC_WPBEH_WHS:
      alias = XC_HYB_GGA_XC_LC_WPBEH_WHS;
      break;
    case BASIC_FUNCTIONALS::XC_LC_WPBESOL_WHS:
      alias = XC_HYB_GGA_XC_LC_WPBESOL_WHS;
      break;
    case BASIC_FUNCTIONALS::XC_LCY_BLYP:
      alias = XC_HYB_GGA_XC_LCY_BLYP;
      break;
    case BASIC_FUNCTIONALS::XC_LCY_PBE:
      alias = XC_HYB_GGA_XC_LCY_PBE;
      break;
    case BASIC_FUNCTIONALS::XC_LRC_WPBE:
      alias = XC_HYB_GGA_XC_LRC_WPBE;
      break;
    case BASIC_FUNCTIONALS::XC_LRC_WPBEH:
      alias = XC_HYB_GGA_XC_LRC_WPBEH;
      break;
    case BASIC_FUNCTIONALS::XC_MB3LYP_RC04:
      alias = XC_HYB_GGA_XC_MB3LYP_RC04;
      break;
    case BASIC_FUNCTIONALS::XC_MPW1K:
      alias = XC_HYB_GGA_XC_MPW1K;
      break;
    case BASIC_FUNCTIONALS::XC_MPW1LYP:
      alias = XC_HYB_GGA_XC_MPW1LYP;
      break;
    case BASIC_FUNCTIONALS::XC_MPW1PBE:
      alias = XC_HYB_GGA_XC_MPW1PBE;
      break;
    case BASIC_FUNCTIONALS::XC_MPW1PW:
      alias = XC_HYB_GGA_XC_MPW1PW;
      break;
    case BASIC_FUNCTIONALS::XC_MPW3LYP:
      alias = XC_HYB_GGA_XC_MPW3LYP;
      break;
    case BASIC_FUNCTIONALS::XC_MPW3PW:
      alias = XC_HYB_GGA_XC_MPW3PW;
      break;
    case BASIC_FUNCTIONALS::XC_MPWLYP1M:
      alias = XC_HYB_GGA_XC_MPWLYP1M;
      break;
    case BASIC_FUNCTIONALS::XC_O3LYP:
      alias = XC_HYB_GGA_XC_O3LYP;
      break;
    case BASIC_FUNCTIONALS::XC_PBE0_13:
      alias = XC_HYB_GGA_XC_PBE0_13;
      break;
    case BASIC_FUNCTIONALS::XC_PBE50:
      alias = XC_HYB_GGA_XC_PBE50;
      break;
    case BASIC_FUNCTIONALS::XC_PBE_MOL0:
      alias = XC_HYB_GGA_XC_PBE_MOL0;
      break;
    case BASIC_FUNCTIONALS::XC_PBE_MOLB0:
      alias = XC_HYB_GGA_XC_PBE_MOLB0;
      break;
    case BASIC_FUNCTIONALS::XC_PBE_SOL0:
      alias = XC_HYB_GGA_XC_PBE_SOL0;
      break;
    case BASIC_FUNCTIONALS::XC_PBEB0:
      alias = XC_HYB_GGA_XC_PBEB0;
      break;
    case BASIC_FUNCTIONALS::XC_PBEH:
      alias = XC_HYB_GGA_XC_PBEH;
      break;
    case BASIC_FUNCTIONALS::XC_QTP17:
      alias = XC_HYB_GGA_XC_QTP17;
      break;
    case BASIC_FUNCTIONALS::XC_RCAM_B3LYP:
      alias = XC_HYB_GGA_XC_RCAM_B3LYP;
      break;
    case BASIC_FUNCTIONALS::XC_REVB3LYP:
      alias = XC_HYB_GGA_XC_REVB3LYP;
      break;
    case BASIC_FUNCTIONALS::XC_SB98_1A:
      alias = XC_HYB_GGA_XC_SB98_1A;
      break;
    case BASIC_FUNCTIONALS::XC_SB98_1B:
      alias = XC_HYB_GGA_XC_SB98_1B;
      break;
    case BASIC_FUNCTIONALS::XC_SB98_1C:
      alias = XC_HYB_GGA_XC_SB98_1C;
      break;
    case BASIC_FUNCTIONALS::XC_SB98_2A:
      alias = XC_HYB_GGA_XC_SB98_2A;
      break;
    case BASIC_FUNCTIONALS::XC_SB98_2B:
      alias = XC_HYB_GGA_XC_SB98_2B;
      break;
    case BASIC_FUNCTIONALS::XC_SB98_2C:
      alias = XC_HYB_GGA_XC_SB98_2C;
      break;
    case BASIC_FUNCTIONALS::XC_TUNED_CAM_B3LYP:
      alias = XC_HYB_GGA_XC_TUNED_CAM_B3LYP;
      break;
    case BASIC_FUNCTIONALS::XC_WB97:
      alias = XC_HYB_GGA_XC_WB97;
      break;
    case BASIC_FUNCTIONALS::XC_WB97X:
      alias = XC_HYB_GGA_XC_WB97X;
      break;
    case BASIC_FUNCTIONALS::XC_WB97X_D:
      alias = XC_HYB_GGA_XC_WB97X_D;
      break;
    case BASIC_FUNCTIONALS::XC_WB97X_V:
      alias = XC_HYB_GGA_XC_WB97X_V;
      break;
    case BASIC_FUNCTIONALS::XC_WC04:
      alias = XC_HYB_GGA_XC_WC04;
      break;
    case BASIC_FUNCTIONALS::XC_WP04:
      alias = XC_HYB_GGA_XC_WP04;
      break;
    case BASIC_FUNCTIONALS::XC_X3LYP:
      alias = XC_HYB_GGA_XC_X3LYP;
      break;
    case BASIC_FUNCTIONALS::X_JS17_2D:
      alias = XC_MGGA_X_2D_JS17;
      break;
    case BASIC_FUNCTIONALS::X_PRHG07_2D:
      alias = XC_MGGA_X_2D_PRHG07;
      break;
    case BASIC_FUNCTIONALS::X_PRHG07_PRP10_2D:
      alias = XC_MGGA_X_2D_PRHG07_PRP10;
      break;
    case BASIC_FUNCTIONALS::X_B00:
      alias = XC_MGGA_X_B00;
      break;
    case BASIC_FUNCTIONALS::X_BJ06:
      alias = XC_MGGA_X_BJ06;
      break;
    case BASIC_FUNCTIONALS::X_BLOC:
      alias = XC_MGGA_X_BLOC;
      break;
    case BASIC_FUNCTIONALS::X_BR89:
      alias = XC_MGGA_X_BR89;
      break;
    case BASIC_FUNCTIONALS::X_BR89_1:
      alias = XC_MGGA_X_BR89_1;
      break;
    case BASIC_FUNCTIONALS::X_BR89_EXPLICIT:
      alias = XC_MGGA_X_BR89_EXPLICIT;
      break;
    case BASIC_FUNCTIONALS::X_BR89_EXPLICIT_1:
      alias = XC_MGGA_X_BR89_EXPLICIT_1;
      break;
    case BASIC_FUNCTIONALS::X_EDMGGA:
      alias = XC_MGGA_X_EDMGGA;
      break;
    case BASIC_FUNCTIONALS::X_GDME_0:
      alias = XC_MGGA_X_GDME_0;
      break;
    case BASIC_FUNCTIONALS::X_GDME_KOS:
      alias = XC_MGGA_X_GDME_KOS;
      break;
    case BASIC_FUNCTIONALS::X_GDME_NV:
      alias = XC_MGGA_X_GDME_NV;
      break;
    case BASIC_FUNCTIONALS::X_GDME_VT:
      alias = XC_MGGA_X_GDME_VT;
      break;
    case BASIC_FUNCTIONALS::X_GVT4:
      alias = XC_MGGA_X_GVT4;
      break;
    case BASIC_FUNCTIONALS::X_GX:
      alias = XC_MGGA_X_GX;
      break;
    case BASIC_FUNCTIONALS::X_LTA:
      alias = XC_MGGA_X_LTA;
      break;
    case BASIC_FUNCTIONALS::X_M06_L:
      alias = XC_MGGA_X_M06_L;
      break;
    case BASIC_FUNCTIONALS::X_M11_L:
      alias = XC_MGGA_X_M11_L;
      break;
    case BASIC_FUNCTIONALS::X_MBEEF:
      alias = XC_MGGA_X_MBEEF;
      break;
    case BASIC_FUNCTIONALS::X_MBEEFVDW:
      alias = XC_MGGA_X_MBEEFVDW;
      break;
    case BASIC_FUNCTIONALS::X_MBRXC_BG:
      alias = XC_MGGA_X_MBRXC_BG;
      break;
    case BASIC_FUNCTIONALS::X_MBRXH_BG:
      alias = XC_MGGA_X_MBRXH_BG;
      break;
    case BASIC_FUNCTIONALS::X_MGGAC:
      alias = XC_MGGA_X_MGGAC;
      break;
    case BASIC_FUNCTIONALS::X_MK00:
      alias = XC_MGGA_X_MK00;
      break;
    case BASIC_FUNCTIONALS::X_MK00B:
      alias = XC_MGGA_X_MK00B;
      break;
    case BASIC_FUNCTIONALS::X_MN12_L:
      alias = XC_MGGA_X_MN12_L;
      break;
    case BASIC_FUNCTIONALS::X_MN15_L:
      alias = XC_MGGA_X_MN15_L;
      break;
    case BASIC_FUNCTIONALS::X_MODTPSS:
      alias = XC_MGGA_X_MODTPSS;
      break;
    case BASIC_FUNCTIONALS::X_MS0:
      alias = XC_MGGA_X_MS0;
      break;
    case BASIC_FUNCTIONALS::X_MS1:
      alias = XC_MGGA_X_MS1;
      break;
    case BASIC_FUNCTIONALS::X_MS2:
      alias = XC_MGGA_X_MS2;
      break;
    case BASIC_FUNCTIONALS::X_MS2B:
      alias = XC_MGGA_X_MS2B;
      break;
    case BASIC_FUNCTIONALS::X_MS2BS:
      alias = XC_MGGA_X_MS2BS;
      break;
    case BASIC_FUNCTIONALS::X_MVS:
      alias = XC_MGGA_X_MVS;
      break;
    case BASIC_FUNCTIONALS::X_MVSB:
      alias = XC_MGGA_X_MVSB;
      break;
    case BASIC_FUNCTIONALS::X_MVSBS:
      alias = XC_MGGA_X_MVSBS;
      break;
    case BASIC_FUNCTIONALS::X_PBE_GX:
      alias = XC_MGGA_X_PBE_GX;
      break;
    case BASIC_FUNCTIONALS::X_PKZB:
      alias = XC_MGGA_X_PKZB;
      break;
    case BASIC_FUNCTIONALS::X_REGTPSS:
      alias = XC_MGGA_X_REGTPSS;
      break;
    case BASIC_FUNCTIONALS::X_REVM06_L:
      alias = XC_MGGA_X_REVM06_L;
      break;
    case BASIC_FUNCTIONALS::X_REVSCAN:
      alias = XC_MGGA_X_REVSCAN;
      break;
    case BASIC_FUNCTIONALS::X_REVSCANL:
      alias = XC_MGGA_X_REVSCANL;
      break;
    case BASIC_FUNCTIONALS::X_REVTM:
      alias = XC_MGGA_X_REVTM;
      break;
    case BASIC_FUNCTIONALS::X_REVTPSS:
      alias = XC_MGGA_X_REVTPSS;
      break;
    case BASIC_FUNCTIONALS::X_RLDA:
      alias = XC_MGGA_X_RLDA;
      break;
    case BASIC_FUNCTIONALS::X_RPP09:
      alias = XC_MGGA_X_RPP09;
      break;
    case BASIC_FUNCTIONALS::X_RSCAN:
      alias = XC_MGGA_X_RSCAN;
      break;
    case BASIC_FUNCTIONALS::X_RTPSS:
      alias = XC_MGGA_X_RTPSS;
      break;
    case BASIC_FUNCTIONALS::X_SA_TPSS:
      alias = XC_MGGA_X_SA_TPSS;
      break;
    case BASIC_FUNCTIONALS::X_SCAN:
      alias = XC_MGGA_X_SCAN;
      break;
    case BASIC_FUNCTIONALS::X_SCANL:
      alias = XC_MGGA_X_SCANL;
      break;
    case BASIC_FUNCTIONALS::X_TASK:
      alias = XC_MGGA_X_TASK;
      break;
    case BASIC_FUNCTIONALS::X_TAU_HCTH_MGGA:
      alias = XC_MGGA_X_TAU_HCTH;
      break;
    case BASIC_FUNCTIONALS::X_TB09:
      alias = XC_MGGA_X_TB09;
      break;
    case BASIC_FUNCTIONALS::X_TLDA:
      alias = XC_MGGA_X_TLDA;
      break;
    case BASIC_FUNCTIONALS::X_TM:
      alias = XC_MGGA_X_TM;
      break;
    case BASIC_FUNCTIONALS::X_TPSS:
      alias = XC_MGGA_X_TPSS;
      break;
    case BASIC_FUNCTIONALS::X_VT84:
      alias = XC_MGGA_X_VT84;
      break;
    case BASIC_FUNCTIONALS::C_B88:
      alias = XC_MGGA_C_B88;
      break;
    case BASIC_FUNCTIONALS::C_BC95:
      alias = XC_MGGA_C_BC95;
      break;
    case BASIC_FUNCTIONALS::C_CS:
      alias = XC_MGGA_C_CS;
      break;
    case BASIC_FUNCTIONALS::C_DLDF:
      alias = XC_MGGA_C_DLDF;
      break;
    case BASIC_FUNCTIONALS::C_KCIS:
      alias = XC_MGGA_C_KCIS;
      break;
    case BASIC_FUNCTIONALS::C_M05:
      alias = XC_MGGA_C_M05;
      break;
    case BASIC_FUNCTIONALS::C_M05_2X:
      alias = XC_MGGA_C_M05_2X;
      break;
    case BASIC_FUNCTIONALS::C_M06:
      alias = XC_MGGA_C_M06;
      break;
    case BASIC_FUNCTIONALS::C_M06_2X:
      alias = XC_MGGA_C_M06_2X;
      break;
    case BASIC_FUNCTIONALS::C_M06_HF:
      alias = XC_MGGA_C_M06_HF;
      break;
    case BASIC_FUNCTIONALS::C_M06_L:
      alias = XC_MGGA_C_M06_L;
      break;
    case BASIC_FUNCTIONALS::C_M08_HX:
      alias = XC_MGGA_C_M08_HX;
      break;
    case BASIC_FUNCTIONALS::C_M08_SO:
      alias = XC_MGGA_C_M08_SO;
      break;
    case BASIC_FUNCTIONALS::C_M11:
      alias = XC_MGGA_C_M11;
      break;
    case BASIC_FUNCTIONALS::C_M11_L:
      alias = XC_MGGA_C_M11_L;
      break;
    case BASIC_FUNCTIONALS::C_MN12_L:
      alias = XC_MGGA_C_MN12_L;
      break;
    case BASIC_FUNCTIONALS::C_MN12_SX:
      alias = XC_MGGA_C_MN12_SX;
      break;
    case BASIC_FUNCTIONALS::C_MN15:
      alias = XC_MGGA_C_MN15;
      break;
    case BASIC_FUNCTIONALS::C_MN15_L:
      alias = XC_MGGA_C_MN15_L;
      break;
    case BASIC_FUNCTIONALS::C_PKZB:
      alias = XC_MGGA_C_PKZB;
      break;
    case BASIC_FUNCTIONALS::C_REVM06:
      alias = XC_MGGA_C_REVM06;
      break;
    case BASIC_FUNCTIONALS::C_REVM06_L:
      alias = XC_MGGA_C_REVM06_L;
      break;
    case BASIC_FUNCTIONALS::C_REVM11:
      alias = XC_MGGA_C_REVM11;
      break;
    case BASIC_FUNCTIONALS::C_REVSCAN:
      alias = XC_MGGA_C_REVSCAN;
      break;
    case BASIC_FUNCTIONALS::C_REVSCAN_VV10:
      alias = XC_MGGA_C_REVSCAN_VV10;
      break;
    case BASIC_FUNCTIONALS::C_REVTM:
      alias = XC_MGGA_C_REVTM;
      break;
    case BASIC_FUNCTIONALS::C_REVTPSS:
      alias = XC_MGGA_C_REVTPSS;
      break;
    case BASIC_FUNCTIONALS::C_RSCAN:
      alias = XC_MGGA_C_RSCAN;
      break;
    case BASIC_FUNCTIONALS::C_SCAN:
      alias = XC_MGGA_C_SCAN;
      break;
    case BASIC_FUNCTIONALS::C_SCAN_RVV10:
      alias = XC_MGGA_C_SCAN_RVV10;
      break;
    case BASIC_FUNCTIONALS::C_SCAN_VV10:
      alias = XC_MGGA_C_SCAN_VV10;
      break;
    case BASIC_FUNCTIONALS::C_SCANL:
      alias = XC_MGGA_C_SCANL;
      break;
    case BASIC_FUNCTIONALS::C_SCANL_RVV10:
      alias = XC_MGGA_C_SCANL_RVV10;
      break;
    case BASIC_FUNCTIONALS::C_SCANL_VV10:
      alias = XC_MGGA_C_SCANL_VV10;
      break;
    case BASIC_FUNCTIONALS::C_TM:
      alias = XC_MGGA_C_TM;
      break;
    case BASIC_FUNCTIONALS::C_TPSS:
      alias = XC_MGGA_C_TPSS;
      break;
    case BASIC_FUNCTIONALS::C_TPSSLOC:
      alias = XC_MGGA_C_TPSSLOC;
      break;
    case BASIC_FUNCTIONALS::C_VSXC:
      alias = XC_MGGA_C_VSXC;
      break;
    case BASIC_FUNCTIONALS::XC_B97M_V:
      alias = XC_MGGA_XC_B97M_V;
      break;
    case BASIC_FUNCTIONALS::XC_CC06:
      alias = XC_MGGA_XC_CC06;
      break;
    case BASIC_FUNCTIONALS::XC_HLE17:
      alias = XC_MGGA_XC_HLE17;
      break;
    case BASIC_FUNCTIONALS::XC_LP90:
      alias = XC_MGGA_XC_LP90;
      break;
    case BASIC_FUNCTIONALS::XC_OTPSS_D:
      alias = XC_MGGA_XC_OTPSS_D;
      break;
    case BASIC_FUNCTIONALS::XC_TPSSLYP1W:
      alias = XC_MGGA_XC_TPSSLYP1W;
      break;
    case BASIC_FUNCTIONALS::XC_ZLP_MGGA:
      alias = XC_MGGA_XC_ZLP;
      break;
    case BASIC_FUNCTIONALS::K_PC07:
      alias = XC_MGGA_K_PC07;
      break;
    case BASIC_FUNCTIONALS::X_BMK:
      alias = XC_HYB_MGGA_X_BMK;
      break;
    case BASIC_FUNCTIONALS::X_DLDF:
      alias = XC_HYB_MGGA_X_DLDF;
      break;
    case BASIC_FUNCTIONALS::X_JS18:
      alias = XC_HYB_MGGA_X_JS18;
      break;
    case BASIC_FUNCTIONALS::X_M05:
      alias = XC_HYB_MGGA_X_M05;
      break;
    case BASIC_FUNCTIONALS::X_M05_2X:
      alias = XC_HYB_MGGA_X_M05_2X;
      break;
    case BASIC_FUNCTIONALS::X_M06:
      alias = XC_HYB_MGGA_X_M06;
      break;
    case BASIC_FUNCTIONALS::X_M06_2X:
      alias = XC_HYB_MGGA_X_M06_2X;
      break;
    case BASIC_FUNCTIONALS::X_M06_HF:
      alias = XC_HYB_MGGA_X_M06_HF;
      break;
    case BASIC_FUNCTIONALS::X_M08_HX:
      alias = XC_HYB_MGGA_X_M08_HX;
      break;
    case BASIC_FUNCTIONALS::X_M08_SO:
      alias = XC_HYB_MGGA_X_M08_SO;
      break;
    case BASIC_FUNCTIONALS::X_M11:
      alias = XC_HYB_MGGA_X_M11;
      break;
    case BASIC_FUNCTIONALS::X_MN12_SX:
      alias = XC_HYB_MGGA_X_MN12_SX;
      break;
    case BASIC_FUNCTIONALS::X_MN15:
      alias = XC_HYB_MGGA_X_MN15;
      break;
    case BASIC_FUNCTIONALS::X_MS2H:
      alias = XC_HYB_MGGA_X_MS2H;
      break;
    case BASIC_FUNCTIONALS::X_MVSH:
      alias = XC_HYB_MGGA_X_MVSH;
      break;
    case BASIC_FUNCTIONALS::X_PJS18:
      alias = XC_HYB_MGGA_X_PJS18;
      break;
    case BASIC_FUNCTIONALS::X_REVM06:
      alias = XC_HYB_MGGA_X_REVM06;
      break;
    case BASIC_FUNCTIONALS::X_REVM11:
      alias = XC_HYB_MGGA_X_REVM11;
      break;
    case BASIC_FUNCTIONALS::X_REVSCAN0:
      alias = XC_HYB_MGGA_X_REVSCAN0;
      break;
    case BASIC_FUNCTIONALS::X_SCAN0:
      alias = XC_HYB_MGGA_X_SCAN0;
      break;
    case BASIC_FUNCTIONALS::X_TAU_HCTH_HYB:
      alias = XC_HYB_MGGA_X_TAU_HCTH;
      break;
    case BASIC_FUNCTIONALS::XC_B0KCIS:
      alias = XC_HYB_MGGA_XC_B0KCIS;
      break;
    case BASIC_FUNCTIONALS::XC_B86B95:
      alias = XC_HYB_MGGA_XC_B86B95;
      break;
    case BASIC_FUNCTIONALS::XC_B88B95:
      alias = XC_HYB_MGGA_XC_B88B95;
      break;
    case BASIC_FUNCTIONALS::XC_B98:
      alias = XC_HYB_MGGA_XC_B98;
      break;
    case BASIC_FUNCTIONALS::XC_BB1K:
      alias = XC_HYB_MGGA_XC_BB1K;
      break;
    case BASIC_FUNCTIONALS::XC_EDMGGAH:
      alias = XC_HYB_MGGA_XC_EDMGGAH;
      break;
    case BASIC_FUNCTIONALS::XC_MPW1B95:
      alias = XC_HYB_MGGA_XC_MPW1B95;
      break;
    case BASIC_FUNCTIONALS::XC_MPW1KCIS:
      alias = XC_HYB_MGGA_XC_MPW1KCIS;
      break;
    case BASIC_FUNCTIONALS::XC_MPWB1K:
      alias = XC_HYB_MGGA_XC_MPWB1K;
      break;
    case BASIC_FUNCTIONALS::XC_MPWKCIS1K:
      alias = XC_HYB_MGGA_XC_MPWKCIS1K;
      break;
    case BASIC_FUNCTIONALS::XC_PBE1KCIS:
      alias = XC_HYB_MGGA_XC_PBE1KCIS;
      break;
    case BASIC_FUNCTIONALS::XC_PW6B95:
      alias = XC_HYB_MGGA_XC_PW6B95;
      break;
    case BASIC_FUNCTIONALS::XC_PW86B95:
      alias = XC_HYB_MGGA_XC_PW86B95;
      break;
    case BASIC_FUNCTIONALS::XC_PWB6K:
      alias = XC_HYB_MGGA_XC_PWB6K;
      break;
    case BASIC_FUNCTIONALS::XC_REVTPSSH:
      alias = XC_HYB_MGGA_XC_REVTPSSH;
      break;
    case BASIC_FUNCTIONALS::XC_TPSS1KCIS:
      alias = XC_HYB_MGGA_XC_TPSS1KCIS;
      break;
    case BASIC_FUNCTIONALS::XC_TPSSH:
      alias = XC_HYB_MGGA_XC_TPSSH;
      break;
    case BASIC_FUNCTIONALS::XC_WB97M_V:
      alias = XC_HYB_MGGA_XC_WB97M_V;
      break;
    case BASIC_FUNCTIONALS::XC_X1B95:
      alias = XC_HYB_MGGA_XC_X1B95;
      break;
    case BASIC_FUNCTIONALS::XC_XB1K:
      alias = XC_HYB_MGGA_XC_XB1K;
      break;
    default:
      throw SerenityError("Basic functional unknown to LibXC.");
      break;
  }
  return alias;
}
#endif /* SERENITY_USE_LIBXC */

#ifdef SERENITY_USE_XCFUN
char* getXCFunAlias(const BASIC_FUNCTIONALS& functional) {
  char* alias = (char*)"null";
  switch (functional) {
    case BASIC_FUNCTIONALS::X_SLATER:
      alias = (char*)"slaterx";
      break;
    case BASIC_FUNCTIONALS::X_LDA_ERF:
      alias = (char*)"ldaerfx";
      break;
    case BASIC_FUNCTIONALS::C_LDA_ERF:
      alias = (char*)"ldaerfc";
      break;
    case BASIC_FUNCTIONALS::C_LDA_ERF_JT:
      alias = (char*)"ldaerfc_jt";
      break;
    case BASIC_FUNCTIONALS::C_VWN:
      alias = (char*)"vwn5c";
      break;
    case BASIC_FUNCTIONALS::C_VWN_3:
      alias = (char*)"vwn3c";
      break;
    case BASIC_FUNCTIONALS::K_TF:
      alias = (char*)"tfk";
      break;
    case BASIC_FUNCTIONALS::X_APBE:
      alias = (char*)"apbex";
      break;
    case BASIC_FUNCTIONALS::X_B86:
      alias = (char*)"beckex";
      break;
    case BASIC_FUNCTIONALS::X_B88:
      alias = (char*)"beckecorrx";
      break;
    case BASIC_FUNCTIONALS::X_B88_CAM:
      alias = (char*)"beckecamx";
      break;
    case BASIC_FUNCTIONALS::X_B88_SR:
      alias = (char*)"beckesrx";
      break;
    case BASIC_FUNCTIONALS::X_KT1:
      alias = (char*)"ktx";
      break;
    case BASIC_FUNCTIONALS::X_OPTX:
      alias = (char*)"optx";
      break;
    case BASIC_FUNCTIONALS::X_PBE:
      alias = (char*)"pbex";
      break;
    case BASIC_FUNCTIONALS::X_PBE_SOL:
      alias = (char*)"pbesolx";
      break;
    case BASIC_FUNCTIONALS::X_PBEINT:
      alias = (char*)"pbeintx";
      break;
    case BASIC_FUNCTIONALS::X_PW86:
      alias = (char*)"pw86x";
      break;
    case BASIC_FUNCTIONALS::X_PW91:
      alias = (char*)"pw91x";
      break;
    case BASIC_FUNCTIONALS::X_REVPBEX:
      alias = (char*)"revpbex";
      break;
    case BASIC_FUNCTIONALS::X_RPBE:
      alias = (char*)"rpbex";
      break;
    case BASIC_FUNCTIONALS::C_APBE:
      alias = (char*)"apbec";
      break;
    case BASIC_FUNCTIONALS::C_LYP:
      alias = (char*)"lypc";
      break;
    case BASIC_FUNCTIONALS::C_OPTC:
      alias = (char*)"optxcorr";
      break;
    case BASIC_FUNCTIONALS::C_P86:
      alias = (char*)"p86c";
      break;
    case BASIC_FUNCTIONALS::C_P86CORRC:
      alias = (char*)"p86corrc";
      break;
    case BASIC_FUNCTIONALS::C_PZ81C:
      alias = (char*)"pz81c";
      break;
    case BASIC_FUNCTIONALS::C_VWN_PBEC:
      alias = (char*)"vwn_pbec";
      break;
    case BASIC_FUNCTIONALS::C_PBE:
      alias = (char*)"pbec";
      break;
    case BASIC_FUNCTIONALS::C_PBEINT:
      alias = (char*)"pbeintc";
      break;
    case BASIC_FUNCTIONALS::C_PBELOC:
      alias = (char*)"pbelocc";
      break;
    case BASIC_FUNCTIONALS::C_PW91:
      alias = (char*)"pw91c";
      break;
    case BASIC_FUNCTIONALS::C_SPBE:
      alias = (char*)"spbec";
      break;
    case BASIC_FUNCTIONALS::K_ERNZERHOF:
      alias = (char*)"E00";
      break;
    case BASIC_FUNCTIONALS::K_PW91:
      alias = (char*)"pw91k";
      break;
    case BASIC_FUNCTIONALS::K_LLP:
      alias = (char*)"llp91k";
      break;
    case BASIC_FUNCTIONALS::K_LLPS:
      alias = (char*)"llp91ks";
      break;
    case BASIC_FUNCTIONALS::K_PBE2:
      alias = (char*)"PBE2";
      break;
    case BASIC_FUNCTIONALS::K_PBE2S:
      alias = (char*)"PBE2S";
      break;
    case BASIC_FUNCTIONALS::K_PBE3:
      alias = (char*)"PBE3";
      break;
    case BASIC_FUNCTIONALS::K_PBE4:
      alias = (char*)"PBE4";
      break;
    case BASIC_FUNCTIONALS::C_B97:
      alias = (char*)"b97c";
      break;
    case BASIC_FUNCTIONALS::C_B97_1:
      alias = (char*)"b97_1c";
      break;
    case BASIC_FUNCTIONALS::C_B97_2:
      alias = (char*)"b97_2c";
      break;
    case BASIC_FUNCTIONALS::X_B97:
      alias = (char*)"b97x";
      break;
    case BASIC_FUNCTIONALS::X_B97_1:
      alias = (char*)"b97_1x";
      break;
    case BASIC_FUNCTIONALS::X_B97_2:
      alias = (char*)"b97_2x";
      break;
    case BASIC_FUNCTIONALS::X_BLOC:
      alias = (char*)"blocx";
      break;
    case BASIC_FUNCTIONALS::X_M06_L:
      alias = (char*)"m06lx";
      break;
    case BASIC_FUNCTIONALS::X_REVTPSS:
      alias = (char*)"revtpssx";
      break;
    case BASIC_FUNCTIONALS::X_TPSS:
      alias = (char*)"tpssx";
      break;
    case BASIC_FUNCTIONALS::C_M05:
      alias = (char*)"m05c";
      break;
    case BASIC_FUNCTIONALS::C_M05_2X:
      alias = (char*)"m05x2c";
      break;
    case BASIC_FUNCTIONALS::C_M06:
      alias = (char*)"m06c";
      break;
    case BASIC_FUNCTIONALS::C_M06_2X:
      alias = (char*)"m06x2c";
      break;
    case BASIC_FUNCTIONALS::C_M06_HF:
      alias = (char*)"m06hfc";
      break;
    case BASIC_FUNCTIONALS::C_M06_L:
      alias = (char*)"m06lc";
      break;
    case BASIC_FUNCTIONALS::C_REVTPSS:
      alias = (char*)"revtpssc";
      break;
    case BASIC_FUNCTIONALS::C_TPSS:
      alias = (char*)"tpssc";
      break;
    case BASIC_FUNCTIONALS::C_TPSSLOC:
      alias = (char*)"tpsslocc";
      break;
    case BASIC_FUNCTIONALS::X_M05:
      alias = (char*)"m05x";
      break;
    case BASIC_FUNCTIONALS::X_M05_2X:
      alias = (char*)"m05x2x";
      break;
    case BASIC_FUNCTIONALS::X_M06:
      alias = (char*)"m06x";
      break;
    case BASIC_FUNCTIONALS::X_M06_2X:
      alias = (char*)"m06x2x";
      break;
    case BASIC_FUNCTIONALS::X_M06_HF:
      alias = (char*)"m06hfx";
      break;
    default:
      throw SerenityError("Basic functional unknown to XCFun.");
      break;
  }
  return alias;
}
#endif /* SERENITY_USE_XCFUN */

} /* namespace BasicFunctionals */
} /* namespace Serenity */
