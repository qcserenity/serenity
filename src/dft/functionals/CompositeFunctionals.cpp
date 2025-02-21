/**
 * @file   CompositeFunctionals.cpp
 *
 * @date   Sep 3, 2020
 * @author Jan P. Unsleber
 *
 * IMPORTANT:\n
 * This file was automatically generated, please do not alter it.
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
#include "dft/functionals/CompositeFunctionals.h"
/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "dft/functionals/BasicFunctionals.h"
#include "misc/SerenityError.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Serenity {
namespace CompositeFunctionals {

Functional resolveFunctional(FUNCTIONALS functional) {
#if defined SERENITY_PREFER_XCFUN && defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC
  try {
    return resolveXCFun(functional);
  }
  catch (...) {
    try {
      return resolveLibXC(functional);
    }
    catch (...) {
      throw SerenityError("You have requested a functional that neither XCFun nor LibXC can provide.");
    }
  }
#elif defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC
  try {
    return resolveLibXC(functional);
  }
  catch (...) {
    try {
      return resolveXCFun(functional);
    }
    catch (...) {
      throw SerenityError("You have requested a functional that neither XCFun nor LibXC can provide.");
    }
  }
#elif defined SERENITY_USE_XCFUN
  return resolveXCFun(functional);
#else
  return resolveLibXC(functional);
#endif
}

Functional resolveFunctional(XCFUNCTIONALS functional) {
  return resolveFunctional(FUNCTIONALS((int)functional));
}

Functional resolveFunctional(KINFUNCTIONALS functional) {
  return resolveFunctional(FUNCTIONALS((int)functional));
}

#ifdef SERENITY_USE_LIBXC
Functional resolveLibXC(FUNCTIONALS functional) {
  using BasicFunctionals::BASIC_FUNCTIONALS;
  switch (functional) {
    case FUNCTIONALS::NONE:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::NONE}, {0.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::HF:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::NONE}, {1.0}, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::SAOP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_SAOP}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::LDA:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::C_VWN}, {1.0, 1.0},
                        0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::HARTREE:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::NONE}, {0.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::BLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_PBE, BASIC_FUNCTIONALS::C_PBE}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::BP86:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_P86}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::KT1:
      return Functional(IMPLEMENTATIONS::LIBXC,
                        {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_KT1, BASIC_FUNCTIONALS::C_VWN},
                        {1.0, -0.006, 1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::KT2:
      return Functional(IMPLEMENTATIONS::LIBXC,
                        {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_KT1, BASIC_FUNCTIONALS::C_VWN},
                        {1.07173, -0.006, 0.576727}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::KT3:
      return Functional(
          IMPLEMENTATIONS::LIBXC,
          {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_KT1, BASIC_FUNCTIONALS::C_OPTC, BASIC_FUNCTIONALS::C_LYP},
          {1.092, -0.004, -0.925452, 0.864409}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::BHLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.50, 1.0}, 0.50,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE0:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_PBE, BASIC_FUNCTIONALS::C_PBE}, {0.75, 1.0}, 0.25,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B3LYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_B3LYP5}, {1.0}, 0.20, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B3LYP_G:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_B3LYP3}, {1.0}, 0.20, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::BPW91:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_PW91}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B97:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_B97}, {1.0}, 0.1943, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B97_1:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_B97_1}, {1.0}, 0.21, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B97_2:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_B97_2}, {1.0}, 0.21, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::CAMB3LYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_CAM_B3LYP}, {1.0}, 0.19, 0.0, 0.46, 0.33, 1.0, 1.0);
    case FUNCTIONALS::B97_D:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_B97_D}, {1.0}, 0.000000, 0.0, 0.000000, 0.00, 1.0, 1.0);
    case FUNCTIONALS::WB97:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_WB97}, {1.0}, 0.000000, 0.0, 1.000000, 0.40, 1.0, 1.0);
    case FUNCTIONALS::WB97X:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_WB97X}, {1.0}, 0.157706, 0.0, 0.842294, 0.30, 1.0, 1.0);
    case FUNCTIONALS::WB97X_D:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_WB97X_D}, {1.0}, 0.222036, 0.0, 0.777964, 0.20,
                        1.0, 1.0);
    case FUNCTIONALS::WB97X_V:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_WB97X_V}, {1.0}, 0.167000, 0.0, 0.833000, 0.30,
                        1.0, 1.0);
    case FUNCTIONALS::VWN5:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::C_VWN}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::VWN3:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::C_VWN_3}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::SLATER:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_SLATER}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::OLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::C_LYP, BASIC_FUNCTIONALS::X_OPTX}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::TF:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::K_TF}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::LCBLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_LC_BLYP}, {1.0}, 0.0, 0.0, 1.0, 0.33, 1.0, 1.0);
    case FUNCTIONALS::LCBLYP_047:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_LC_BLYP}, {1.0}, 0.0, 0.0, 1.0, 0.47, 1.0, 1.0);
    case FUNCTIONALS::LCBLYP_100:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::XC_LC_BLYP}, {1.0}, 0.0, 0.0, 1.0, 1.00, 1.0, 1.0);
    case FUNCTIONALS::PW91K:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::K_PW91}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PW91:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::C_PW91, BASIC_FUNCTIONALS::X_PW91}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2PLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.47, 0.73},
                        0.53, 0.27, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2KPLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.28, 0.58},
                        0.72, 0.42, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2TPLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.40, 0.69},
                        0.60, 0.31, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2GPPLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.35, 0.64},
                        0.65, 0.36, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::ROB2PLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.41, 0.72},
                        0.59, 0.28, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2PIPLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.398, 0.727},
                        0.602, 0.273, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2PPW91:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_PW91}, {0.8, 0.9}, 0.2,
                        0.1, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::DSDBLYP:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.31, 0.54},
                        0.69, 1.0, 0.0, 0.0, 0.37, 0.46);
    case FUNCTIONALS::DUT:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.30, 0.59},
                        0.70, 1.0, 0.0, 0.0, 0.36, 0.47);
    case FUNCTIONALS::PUT:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.32, 0.63},
                        0.68, 1.0, 0.0, 0.0, 0.27, 0.46);
    case FUNCTIONALS::DSDPBEP86:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::X_PBE, BASIC_FUNCTIONALS::C_P86}, {0.30, 0.43},
                        0.70, 1.0, 0.0, 0.0, 0.25, 0.53);
    case FUNCTIONALS::LLP91K:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::K_LLP}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE3K:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::K_PBE3}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE4K:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::K_PBE4}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::E2000K:
      return Functional(IMPLEMENTATIONS::LIBXC, {BASIC_FUNCTIONALS::K_ERNZERHOF}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    default:
      throw SerenityError("Composite functional unknown to LibXC.");
      break;
  }
}

#endif /* SERENITY_USE_LIBXC */

#ifdef SERENITY_USE_XCFUN
Functional resolveXCFun(FUNCTIONALS functional) {
  using BasicFunctionals::BASIC_FUNCTIONALS;
  switch (functional) {
    case FUNCTIONALS::NONE:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::NONE}, {0.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::HF:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::NONE}, {1.0}, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::SAOP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::XC_SAOP}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::LDA:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::C_VWN}, {1.0, 1.0},
                        0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::HARTREE:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::NONE}, {0.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::BLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_PBE, BASIC_FUNCTIONALS::C_PBE}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::BP86:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_P86}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::KT1:
      return Functional(IMPLEMENTATIONS::XCFUN,
                        {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_KT1, BASIC_FUNCTIONALS::C_VWN},
                        {1.0, -0.006, 1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::KT2:
      return Functional(IMPLEMENTATIONS::XCFUN,
                        {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_KT1, BASIC_FUNCTIONALS::C_VWN},
                        {1.07173, -0.006, 0.576727}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::KT3:
      return Functional(
          IMPLEMENTATIONS::XCFUN,
          {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_KT1, BASIC_FUNCTIONALS::C_OPTC, BASIC_FUNCTIONALS::C_LYP},
          {1.092, -0.004, -0.925452, 0.864409}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::LDAERF:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_LDA_ERF, BASIC_FUNCTIONALS::C_LDA_ERF},
                        {1.0, 1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::LDAERF_JT:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_LDA_ERF, BASIC_FUNCTIONALS::C_LDA_ERF_JT},
                        {1.0, 1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::BHLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.50, 1.0}, 0.50,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE0:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_PBE, BASIC_FUNCTIONALS::C_PBE}, {0.75, 1.0}, 0.25,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B3LYP:
      return Functional(IMPLEMENTATIONS::XCFUN,
                        {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_B88_CORR, BASIC_FUNCTIONALS::C_LYP,
                         BASIC_FUNCTIONALS::C_VWN},
                        {0.80, 0.72, 0.81, 0.19}, 0.20, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B3LYP_G:
      return Functional(IMPLEMENTATIONS::XCFUN,
                        {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_B88_CORR, BASIC_FUNCTIONALS::C_LYP,
                         BASIC_FUNCTIONALS::C_VWN_3},
                        {0.80, 0.72, 0.81, 0.19}, 0.20, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B3P86:
      return Functional(IMPLEMENTATIONS::XCFUN,
                        {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_B88_CORR, BASIC_FUNCTIONALS::C_P86CORRC,
                         BASIC_FUNCTIONALS::C_VWN},
                        {0.80, 0.72, 0.81, 1.0}, 0.20, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B3P86_G:
      return Functional(IMPLEMENTATIONS::XCFUN,
                        {BASIC_FUNCTIONALS::X_SLATER, BASIC_FUNCTIONALS::X_B88_CORR, BASIC_FUNCTIONALS::C_P86CORRC,
                         BASIC_FUNCTIONALS::C_VWN_3},
                        {0.80, 0.72, 0.81, 1.0}, 0.20, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::BPW91:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_PW91}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B97:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B97, BASIC_FUNCTIONALS::C_B97}, {1.0, 1.0},
                        0.1943, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B97_1:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B97_1, BASIC_FUNCTIONALS::C_B97_1}, {1.0, 1.0},
                        0.21, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B97_2:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B97_2, BASIC_FUNCTIONALS::C_B97_2}, {1.0, 1.0},
                        0.21, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::CAMB3LYP:
      return Functional(IMPLEMENTATIONS::XCFUN,
                        {BASIC_FUNCTIONALS::X_B88_CAM, BASIC_FUNCTIONALS::C_VWN, BASIC_FUNCTIONALS::C_LYP},
                        {1.0, 0.19, 0.81}, 0.19, 0.0, 0.46, 0.33, 1.0, 1.0);
    case FUNCTIONALS::VWN5:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::C_VWN}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::VWN3:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::C_VWN_3}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::SLATER:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_SLATER}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::OLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::C_LYP, BASIC_FUNCTIONALS::X_OPTX}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::TF:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_TF}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::LCBLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88_CAM, BASIC_FUNCTIONALS::C_LYP}, {1.0, 1.0},
                        0.0, 0.0, 1.0, 0.33, 1.0, 1.0);
    case FUNCTIONALS::LCBLYP_047:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88_CAM, BASIC_FUNCTIONALS::C_LYP}, {1.0, 1.0},
                        0.0, 0.0, 1.0, 0.47, 1.0, 1.0);
    case FUNCTIONALS::LCBLYP_100:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88_CAM, BASIC_FUNCTIONALS::C_LYP}, {1.0, 1.0},
                        0.0, 0.0, 1.0, 1.00, 1.0, 1.0);
    case FUNCTIONALS::PW91K:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_PW91}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PW91:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::C_PW91, BASIC_FUNCTIONALS::X_PW91}, {1.0, 1.0}, 0.0,
                        0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2PLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.47, 0.73},
                        0.53, 0.27, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2KPLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.28, 0.58},
                        0.72, 0.42, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2TPLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.40, 0.69},
                        0.60, 0.31, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2GPPLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.35, 0.64},
                        0.65, 0.36, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::ROB2PLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.41, 0.72},
                        0.59, 0.28, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2PIPLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.398, 0.727},
                        0.602, 0.273, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::B2PPW91:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_PW91}, {0.8, 0.9}, 0.2,
                        0.1, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::DSDBLYP:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.31, 0.54},
                        0.69, 1.0, 0.0, 0.0, 0.37, 0.46);
    case FUNCTIONALS::DUT:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.30, 0.59},
                        0.70, 1.0, 0.0, 0.0, 0.36, 0.47);
    case FUNCTIONALS::PUT:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_B88, BASIC_FUNCTIONALS::C_LYP}, {0.32, 0.63},
                        0.68, 1.0, 0.0, 0.0, 0.27, 0.46);
    case FUNCTIONALS::DSDPBEP86:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::X_PBE, BASIC_FUNCTIONALS::C_P86}, {0.30, 0.43},
                        0.70, 1.0, 0.0, 0.0, 0.25, 0.53);
    case FUNCTIONALS::LLP91K:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_LLP}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::LLP91KS:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_LLPS}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE2K:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_PBE2}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE2KS:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_PBE2S}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE3K:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_PBE3}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::PBE4K:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_PBE4}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    case FUNCTIONALS::E2000K:
      return Functional(IMPLEMENTATIONS::XCFUN, {BASIC_FUNCTIONALS::K_ERNZERHOF}, {1.0}, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    default:
      throw SerenityError("Composite functional unknown to XCFun.");
      break;
  }
}

#endif /* SERENITY_USE_XCFUN */

} /* namespace CompositeFunctionals */
namespace Options {
template<>
void resolve<CompositeFunctionals::FUNCTIONALS>(std::string& value, CompositeFunctionals::FUNCTIONALS& field) {
  static const std::map<std::string, CompositeFunctionals::FUNCTIONALS> m = {
      {"NONE", CompositeFunctionals::FUNCTIONALS::NONE},
      {"SLATER", CompositeFunctionals::FUNCTIONALS::SLATER},
      {"VWN3", CompositeFunctionals::FUNCTIONALS::VWN3},
      {"VWN5", CompositeFunctionals::FUNCTIONALS::VWN5},
      {"LDAERF", CompositeFunctionals::FUNCTIONALS::LDAERF},
      {"LDA_ERF", CompositeFunctionals::FUNCTIONALS::LDAERF},
      {"LDA-ERF", CompositeFunctionals::FUNCTIONALS::LDAERF},
      {"LDAERF_JT", CompositeFunctionals::FUNCTIONALS::LDAERF_JT},
      {"LDAERFJT", CompositeFunctionals::FUNCTIONALS::LDAERF_JT},
      {"LDA_ERF_JT", CompositeFunctionals::FUNCTIONALS::LDAERF_JT},
      {"LDA-ERF-JT", CompositeFunctionals::FUNCTIONALS::LDAERF_JT},
      {"LDA", CompositeFunctionals::FUNCTIONALS::LDA},
      {"HARTREE", CompositeFunctionals::FUNCTIONALS::HARTREE},
      {"B97", CompositeFunctionals::FUNCTIONALS::B97},
      {"B97_1", CompositeFunctionals::FUNCTIONALS::B97_1},
      {"B97-1", CompositeFunctionals::FUNCTIONALS::B97_1},
      {"B97_2", CompositeFunctionals::FUNCTIONALS::B97_2},
      {"B97-2", CompositeFunctionals::FUNCTIONALS::B97_2},
      {"OLYP", CompositeFunctionals::FUNCTIONALS::OLYP},
      {"BLYP", CompositeFunctionals::FUNCTIONALS::BLYP},
      {"PBE", CompositeFunctionals::FUNCTIONALS::PBE},
      {"BP86", CompositeFunctionals::FUNCTIONALS::BP86},
      {"KT1", CompositeFunctionals::FUNCTIONALS::KT1},
      {"KT2", CompositeFunctionals::FUNCTIONALS::KT2},
      {"KT3", CompositeFunctionals::FUNCTIONALS::KT3},
      {"PW91", CompositeFunctionals::FUNCTIONALS::PW91},
      {"BHLYP", CompositeFunctionals::FUNCTIONALS::BHLYP},
      {"PBE0", CompositeFunctionals::FUNCTIONALS::PBE0},
      {"B3LYP", CompositeFunctionals::FUNCTIONALS::B3LYP},
      {"B3LYP_G", CompositeFunctionals::FUNCTIONALS::B3LYP_G},
      {"B3LYP-G", CompositeFunctionals::FUNCTIONALS::B3LYP_G},
      {"B3P86", CompositeFunctionals::FUNCTIONALS::B3P86},
      {"B3P86_G", CompositeFunctionals::FUNCTIONALS::B3P86_G},
      {"B3P86-G", CompositeFunctionals::FUNCTIONALS::B3P86_G},
      {"BPW91", CompositeFunctionals::FUNCTIONALS::BPW91},
      {"CAMB3LYP", CompositeFunctionals::FUNCTIONALS::CAMB3LYP},
      {"CAM-B3LYP", CompositeFunctionals::FUNCTIONALS::CAMB3LYP},
      {"CAM_B3LYP", CompositeFunctionals::FUNCTIONALS::CAMB3LYP},
      {"LCBLYP", CompositeFunctionals::FUNCTIONALS::LCBLYP},
      {"LCBLYP_047", CompositeFunctionals::FUNCTIONALS::LCBLYP_047},
      {"LCBLYP-047", CompositeFunctionals::FUNCTIONALS::LCBLYP_047},
      {"LCBLYP_100", CompositeFunctionals::FUNCTIONALS::LCBLYP_100},
      {"LCBLYP-100", CompositeFunctionals::FUNCTIONALS::LCBLYP_100},
      {"B2PLYP", CompositeFunctionals::FUNCTIONALS::B2PLYP},
      {"B2KPLYP", CompositeFunctionals::FUNCTIONALS::B2KPLYP},
      {"B2TPLYP", CompositeFunctionals::FUNCTIONALS::B2TPLYP},
      {"B2GPPLYP", CompositeFunctionals::FUNCTIONALS::B2GPPLYP},
      {"ROB2PLYP", CompositeFunctionals::FUNCTIONALS::ROB2PLYP},
      {"B2PIPLYP", CompositeFunctionals::FUNCTIONALS::B2PIPLYP},
      {"B2PPW91", CompositeFunctionals::FUNCTIONALS::B2PPW91},
      {"DSDBLYP", CompositeFunctionals::FUNCTIONALS::DSDBLYP},
      {"DUT", CompositeFunctionals::FUNCTIONALS::DUT},
      {"PUT", CompositeFunctionals::FUNCTIONALS::PUT},
      {"DSDPBEP86", CompositeFunctionals::FUNCTIONALS::DSDPBEP86},
      {"SAOP", CompositeFunctionals::FUNCTIONALS::SAOP},
      {"HF", CompositeFunctionals::FUNCTIONALS::HF},
      {"HARTREE-FOCK", CompositeFunctionals::FUNCTIONALS::HF},
      {"HARTREE_FOCK", CompositeFunctionals::FUNCTIONALS::HF},
      {"TF", CompositeFunctionals::FUNCTIONALS::TF},
      {"PW91K", CompositeFunctionals::FUNCTIONALS::PW91K},
      {"LLP91K", CompositeFunctionals::FUNCTIONALS::LLP91K},
      {"LLP91KS", CompositeFunctionals::FUNCTIONALS::LLP91KS},
      {"PBE2K", CompositeFunctionals::FUNCTIONALS::PBE2K},
      {"PBE2", CompositeFunctionals::FUNCTIONALS::PBE2K},
      {"PBE2KS", CompositeFunctionals::FUNCTIONALS::PBE2KS},
      {"PBE2S", CompositeFunctionals::FUNCTIONALS::PBE2KS},
      {"PBE3K", CompositeFunctionals::FUNCTIONALS::PBE3K},
      {"PBE3", CompositeFunctionals::FUNCTIONALS::PBE3K},
      {"PBE4K", CompositeFunctionals::FUNCTIONALS::PBE4K},
      {"PBE4", CompositeFunctionals::FUNCTIONALS::PBE4K},
      {"E2000K", CompositeFunctionals::FUNCTIONALS::E2000K},
      {"E2000", CompositeFunctionals::FUNCTIONALS::E2000K},
      {"E00", CompositeFunctionals::FUNCTIONALS::E2000K},
      {"B97_D", CompositeFunctionals::FUNCTIONALS::B97_D},
      {"B97-D", CompositeFunctionals::FUNCTIONALS::B97_D},
      {"WB97", CompositeFunctionals::FUNCTIONALS::WB97},
      {"WB97X", CompositeFunctionals::FUNCTIONALS::WB97X},
      {"WB97X_D", CompositeFunctionals::FUNCTIONALS::WB97X_D},
      {"WB97X-D", CompositeFunctionals::FUNCTIONALS::WB97X_D},
      {"WB97X_V", CompositeFunctionals::FUNCTIONALS::WB97X_V},
      {"WB97X-V", CompositeFunctionals::FUNCTIONALS::WB97X_V}};
  check(m, value, field);
}

template<>
void resolve<CompositeFunctionals::XCFUNCTIONALS>(std::string& value, CompositeFunctionals::XCFUNCTIONALS& field) {
  static const std::map<std::string, CompositeFunctionals::XCFUNCTIONALS> m = {
      {"NONE", CompositeFunctionals::XCFUNCTIONALS::NONE},
      {"SLATER", CompositeFunctionals::XCFUNCTIONALS::SLATER},
      {"VWN3", CompositeFunctionals::XCFUNCTIONALS::VWN3},
      {"VWN5", CompositeFunctionals::XCFUNCTIONALS::VWN5},
      {"LDAERF", CompositeFunctionals::XCFUNCTIONALS::LDAERF},
      {"LDA_ERF", CompositeFunctionals::XCFUNCTIONALS::LDAERF},
      {"LDA-ERF", CompositeFunctionals::XCFUNCTIONALS::LDAERF},
      {"LDAERF_JT", CompositeFunctionals::XCFUNCTIONALS::LDAERF_JT},
      {"LDAERFJT", CompositeFunctionals::XCFUNCTIONALS::LDAERF_JT},
      {"LDA_ERF_JT", CompositeFunctionals::XCFUNCTIONALS::LDAERF_JT},
      {"LDA-ERF-JT", CompositeFunctionals::XCFUNCTIONALS::LDAERF_JT},
      {"LDA", CompositeFunctionals::XCFUNCTIONALS::LDA},
      {"HARTREE", CompositeFunctionals::XCFUNCTIONALS::HARTREE},
      {"B97", CompositeFunctionals::XCFUNCTIONALS::B97},
      {"B97_1", CompositeFunctionals::XCFUNCTIONALS::B97_1},
      {"B97-1", CompositeFunctionals::XCFUNCTIONALS::B97_1},
      {"B97_2", CompositeFunctionals::XCFUNCTIONALS::B97_2},
      {"B97-2", CompositeFunctionals::XCFUNCTIONALS::B97_2},
      {"OLYP", CompositeFunctionals::XCFUNCTIONALS::OLYP},
      {"BLYP", CompositeFunctionals::XCFUNCTIONALS::BLYP},
      {"PBE", CompositeFunctionals::XCFUNCTIONALS::PBE},
      {"BP86", CompositeFunctionals::XCFUNCTIONALS::BP86},
      {"KT1", CompositeFunctionals::XCFUNCTIONALS::KT1},
      {"KT2", CompositeFunctionals::XCFUNCTIONALS::KT2},
      {"KT3", CompositeFunctionals::XCFUNCTIONALS::KT3},
      {"PW91", CompositeFunctionals::XCFUNCTIONALS::PW91},
      {"BHLYP", CompositeFunctionals::XCFUNCTIONALS::BHLYP},
      {"PBE0", CompositeFunctionals::XCFUNCTIONALS::PBE0},
      {"B3LYP", CompositeFunctionals::XCFUNCTIONALS::B3LYP},
      {"B3LYP_G", CompositeFunctionals::XCFUNCTIONALS::B3LYP_G},
      {"B3LYP-G", CompositeFunctionals::XCFUNCTIONALS::B3LYP_G},
      {"B3P86", CompositeFunctionals::XCFUNCTIONALS::B3P86},
      {"B3P86_G", CompositeFunctionals::XCFUNCTIONALS::B3P86_G},
      {"B3P86-G", CompositeFunctionals::XCFUNCTIONALS::B3P86_G},
      {"BPW91", CompositeFunctionals::XCFUNCTIONALS::BPW91},
      {"CAMB3LYP", CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP},
      {"CAM-B3LYP", CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP},
      {"CAM_B3LYP", CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP},
      {"LCBLYP", CompositeFunctionals::XCFUNCTIONALS::LCBLYP},
      {"LCBLYP_047", CompositeFunctionals::XCFUNCTIONALS::LCBLYP_047},
      {"LCBLYP-047", CompositeFunctionals::XCFUNCTIONALS::LCBLYP_047},
      {"LCBLYP_100", CompositeFunctionals::XCFUNCTIONALS::LCBLYP_100},
      {"LCBLYP-100", CompositeFunctionals::XCFUNCTIONALS::LCBLYP_100},
      {"B2PLYP", CompositeFunctionals::XCFUNCTIONALS::B2PLYP},
      {"B2KPLYP", CompositeFunctionals::XCFUNCTIONALS::B2KPLYP},
      {"B2TPLYP", CompositeFunctionals::XCFUNCTIONALS::B2TPLYP},
      {"B2GPPLYP", CompositeFunctionals::XCFUNCTIONALS::B2GPPLYP},
      {"ROB2PLYP", CompositeFunctionals::XCFUNCTIONALS::ROB2PLYP},
      {"B2PIPLYP", CompositeFunctionals::XCFUNCTIONALS::B2PIPLYP},
      {"B2PPW91", CompositeFunctionals::XCFUNCTIONALS::B2PPW91},
      {"DSDBLYP", CompositeFunctionals::XCFUNCTIONALS::DSDBLYP},
      {"DUT", CompositeFunctionals::XCFUNCTIONALS::DUT},
      {"PUT", CompositeFunctionals::XCFUNCTIONALS::PUT},
      {"DSDPBEP86", CompositeFunctionals::XCFUNCTIONALS::DSDPBEP86},
      {"SAOP", CompositeFunctionals::XCFUNCTIONALS::SAOP},
      {"HF", CompositeFunctionals::XCFUNCTIONALS::HF},
      {"HARTREE-FOCK", CompositeFunctionals::XCFUNCTIONALS::HF},
      {"HARTREE_FOCK", CompositeFunctionals::XCFUNCTIONALS::HF},
      {"B97_D", CompositeFunctionals::XCFUNCTIONALS::B97_D},
      {"B97-D", CompositeFunctionals::XCFUNCTIONALS::B97_D},
      {"WB97", CompositeFunctionals::XCFUNCTIONALS::WB97},
      {"WB97X", CompositeFunctionals::XCFUNCTIONALS::WB97X},
      {"WB97X_D", CompositeFunctionals::XCFUNCTIONALS::WB97X_D},
      {"WB97X-D", CompositeFunctionals::XCFUNCTIONALS::WB97X_D},
      {"WB97X_V", CompositeFunctionals::XCFUNCTIONALS::WB97X_V},
      {"WB97X-V", CompositeFunctionals::XCFUNCTIONALS::WB97X_V}};
  check(m, value, field);
}

template<>
void resolve<CompositeFunctionals::KINFUNCTIONALS>(std::string& value, CompositeFunctionals::KINFUNCTIONALS& field) {
  static const std::map<std::string, CompositeFunctionals::KINFUNCTIONALS> m = {
      {"NONE", CompositeFunctionals::KINFUNCTIONALS::NONE},
      {"TF", CompositeFunctionals::KINFUNCTIONALS::TF},
      {"PW91K", CompositeFunctionals::KINFUNCTIONALS::PW91K},
      {"LLP91K", CompositeFunctionals::KINFUNCTIONALS::LLP91K},
      {"LLP91KS", CompositeFunctionals::KINFUNCTIONALS::LLP91KS},
      {"PBE2K", CompositeFunctionals::KINFUNCTIONALS::PBE2K},
      {"PBE2", CompositeFunctionals::KINFUNCTIONALS::PBE2K},
      {"PBE2KS", CompositeFunctionals::KINFUNCTIONALS::PBE2KS},
      {"PBE2S", CompositeFunctionals::KINFUNCTIONALS::PBE2KS},
      {"PBE3K", CompositeFunctionals::KINFUNCTIONALS::PBE3K},
      {"PBE3", CompositeFunctionals::KINFUNCTIONALS::PBE3K},
      {"PBE4K", CompositeFunctionals::KINFUNCTIONALS::PBE4K},
      {"PBE4", CompositeFunctionals::KINFUNCTIONALS::PBE4K},
      {"E2000K", CompositeFunctionals::KINFUNCTIONALS::E2000K},
      {"E2000", CompositeFunctionals::KINFUNCTIONALS::E2000K},
      {"E00", CompositeFunctionals::KINFUNCTIONALS::E2000K}};
  check(m, value, field);
}
} /* namespace Options */

} /* namespace Serenity */
