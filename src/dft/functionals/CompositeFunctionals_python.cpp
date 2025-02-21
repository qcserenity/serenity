/**
 * @file   CompositeFunctionals_python.cpp
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
/* Include Serenity Internal Headers */
#include "dft/functionals/CompositeFunctionals.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

using namespace Serenity;
namespace py = pybind11;

void export_CompositeFunctionals(py::module& spy) {
  py::enum_<CompositeFunctionals::FUNCTIONALS>(spy, "FUNCTIONALS")
      .value("NONE", CompositeFunctionals::FUNCTIONALS::NONE)
      .value("SLATER", CompositeFunctionals::FUNCTIONALS::SLATER)
      .value("VWN3", CompositeFunctionals::FUNCTIONALS::VWN3)
      .value("VWN5", CompositeFunctionals::FUNCTIONALS::VWN5)
      .value("LDAERF", CompositeFunctionals::FUNCTIONALS::LDAERF)
      .value("LDAERF_JT", CompositeFunctionals::FUNCTIONALS::LDAERF_JT)
      .value("LDA", CompositeFunctionals::FUNCTIONALS::LDA)
      .value("HARTREE", CompositeFunctionals::FUNCTIONALS::HARTREE)
      .value("B97", CompositeFunctionals::FUNCTIONALS::B97)
      .value("B97_1", CompositeFunctionals::FUNCTIONALS::B97_1)
      .value("B97_2", CompositeFunctionals::FUNCTIONALS::B97_2)
      .value("OLYP", CompositeFunctionals::FUNCTIONALS::OLYP)
      .value("BLYP", CompositeFunctionals::FUNCTIONALS::BLYP)
      .value("PBE", CompositeFunctionals::FUNCTIONALS::PBE)
      .value("BP86", CompositeFunctionals::FUNCTIONALS::BP86)
      .value("KT1", CompositeFunctionals::FUNCTIONALS::KT1)
      .value("KT2", CompositeFunctionals::FUNCTIONALS::KT2)
      .value("KT3", CompositeFunctionals::FUNCTIONALS::KT3)
      .value("PW91", CompositeFunctionals::FUNCTIONALS::PW91)
      .value("BHLYP", CompositeFunctionals::FUNCTIONALS::BHLYP)
      .value("PBE0", CompositeFunctionals::FUNCTIONALS::PBE0)
      .value("B3LYP", CompositeFunctionals::FUNCTIONALS::B3LYP)
      .value("B3LYP_G", CompositeFunctionals::FUNCTIONALS::B3LYP_G)
      .value("B3P86", CompositeFunctionals::FUNCTIONALS::B3P86)
      .value("B3P86_G", CompositeFunctionals::FUNCTIONALS::B3P86_G)
      .value("BPW91", CompositeFunctionals::FUNCTIONALS::BPW91)
      .value("CAMB3LYP", CompositeFunctionals::FUNCTIONALS::CAMB3LYP)
      .value("LCBLYP", CompositeFunctionals::FUNCTIONALS::LCBLYP)
      .value("LCBLYP_047", CompositeFunctionals::FUNCTIONALS::LCBLYP_047)
      .value("LCBLYP_100", CompositeFunctionals::FUNCTIONALS::LCBLYP_100)
      .value("B2PLYP", CompositeFunctionals::FUNCTIONALS::B2PLYP)
      .value("B2KPLYP", CompositeFunctionals::FUNCTIONALS::B2KPLYP)
      .value("B2TPLYP", CompositeFunctionals::FUNCTIONALS::B2TPLYP)
      .value("B2GPPLYP", CompositeFunctionals::FUNCTIONALS::B2GPPLYP)
      .value("ROB2PLYP", CompositeFunctionals::FUNCTIONALS::ROB2PLYP)
      .value("B2PIPLYP", CompositeFunctionals::FUNCTIONALS::B2PIPLYP)
      .value("B2PPW91", CompositeFunctionals::FUNCTIONALS::B2PPW91)
      .value("DSDBLYP", CompositeFunctionals::FUNCTIONALS::DSDBLYP)
      .value("DUT", CompositeFunctionals::FUNCTIONALS::DUT)
      .value("PUT", CompositeFunctionals::FUNCTIONALS::PUT)
      .value("DSDPBEP86", CompositeFunctionals::FUNCTIONALS::DSDPBEP86)
      .value("SAOP", CompositeFunctionals::FUNCTIONALS::SAOP)
      .value("HF", CompositeFunctionals::FUNCTIONALS::HF)
      .value("TF", CompositeFunctionals::FUNCTIONALS::TF)
      .value("PW91K", CompositeFunctionals::FUNCTIONALS::PW91K)
      .value("LLP91K", CompositeFunctionals::FUNCTIONALS::LLP91K)
      .value("LLP91KS", CompositeFunctionals::FUNCTIONALS::LLP91KS)
      .value("PBE2K", CompositeFunctionals::FUNCTIONALS::PBE2K)
      .value("PBE2KS", CompositeFunctionals::FUNCTIONALS::PBE2KS)
      .value("PBE3K", CompositeFunctionals::FUNCTIONALS::PBE3K)
      .value("PBE4K", CompositeFunctionals::FUNCTIONALS::PBE4K)
      .value("E2000K", CompositeFunctionals::FUNCTIONALS::E2000K)
      .value("B97_D", CompositeFunctionals::FUNCTIONALS::B97_D)
      .value("WB97", CompositeFunctionals::FUNCTIONALS::WB97)
      .value("WB97X", CompositeFunctionals::FUNCTIONALS::WB97X)
      .value("WB97X_D", CompositeFunctionals::FUNCTIONALS::WB97X_D)
      .value("WB97X_V", CompositeFunctionals::FUNCTIONALS::WB97X_V)
      .export_values();
  py::enum_<CompositeFunctionals::XCFUNCTIONALS>(spy, "XCFUNCTIONALS")
      .value("NONE", CompositeFunctionals::XCFUNCTIONALS::NONE)
      .value("SLATER", CompositeFunctionals::XCFUNCTIONALS::SLATER)
      .value("VWN3", CompositeFunctionals::XCFUNCTIONALS::VWN3)
      .value("VWN5", CompositeFunctionals::XCFUNCTIONALS::VWN5)
      .value("LDAERF", CompositeFunctionals::XCFUNCTIONALS::LDAERF)
      .value("LDAERF_JT", CompositeFunctionals::XCFUNCTIONALS::LDAERF_JT)
      .value("LDA", CompositeFunctionals::XCFUNCTIONALS::LDA)
      .value("HARTREE", CompositeFunctionals::XCFUNCTIONALS::HARTREE)
      .value("B97", CompositeFunctionals::XCFUNCTIONALS::B97)
      .value("B97_1", CompositeFunctionals::XCFUNCTIONALS::B97_1)
      .value("B97_2", CompositeFunctionals::XCFUNCTIONALS::B97_2)
      .value("OLYP", CompositeFunctionals::XCFUNCTIONALS::OLYP)
      .value("BLYP", CompositeFunctionals::XCFUNCTIONALS::BLYP)
      .value("PBE", CompositeFunctionals::XCFUNCTIONALS::PBE)
      .value("BP86", CompositeFunctionals::XCFUNCTIONALS::BP86)
      .value("KT1", CompositeFunctionals::XCFUNCTIONALS::KT1)
      .value("KT2", CompositeFunctionals::XCFUNCTIONALS::KT2)
      .value("KT3", CompositeFunctionals::XCFUNCTIONALS::KT3)
      .value("PW91", CompositeFunctionals::XCFUNCTIONALS::PW91)
      .value("BHLYP", CompositeFunctionals::XCFUNCTIONALS::BHLYP)
      .value("PBE0", CompositeFunctionals::XCFUNCTIONALS::PBE0)
      .value("B3LYP", CompositeFunctionals::XCFUNCTIONALS::B3LYP)
      .value("B3LYP_G", CompositeFunctionals::XCFUNCTIONALS::B3LYP_G)
      .value("B3P86", CompositeFunctionals::XCFUNCTIONALS::B3P86)
      .value("B3P86_G", CompositeFunctionals::XCFUNCTIONALS::B3P86_G)
      .value("BPW91", CompositeFunctionals::XCFUNCTIONALS::BPW91)
      .value("CAMB3LYP", CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP)
      .value("LCBLYP", CompositeFunctionals::XCFUNCTIONALS::LCBLYP)
      .value("LCBLYP_047", CompositeFunctionals::XCFUNCTIONALS::LCBLYP_047)
      .value("LCBLYP_100", CompositeFunctionals::XCFUNCTIONALS::LCBLYP_100)
      .value("B2PLYP", CompositeFunctionals::XCFUNCTIONALS::B2PLYP)
      .value("B2KPLYP", CompositeFunctionals::XCFUNCTIONALS::B2KPLYP)
      .value("B2TPLYP", CompositeFunctionals::XCFUNCTIONALS::B2TPLYP)
      .value("B2GPPLYP", CompositeFunctionals::XCFUNCTIONALS::B2GPPLYP)
      .value("ROB2PLYP", CompositeFunctionals::XCFUNCTIONALS::ROB2PLYP)
      .value("B2PIPLYP", CompositeFunctionals::XCFUNCTIONALS::B2PIPLYP)
      .value("B2PPW91", CompositeFunctionals::XCFUNCTIONALS::B2PPW91)
      .value("DSDBLYP", CompositeFunctionals::XCFUNCTIONALS::DSDBLYP)
      .value("DUT", CompositeFunctionals::XCFUNCTIONALS::DUT)
      .value("PUT", CompositeFunctionals::XCFUNCTIONALS::PUT)
      .value("DSDPBEP86", CompositeFunctionals::XCFUNCTIONALS::DSDPBEP86)
      .value("SAOP", CompositeFunctionals::XCFUNCTIONALS::SAOP)
      .value("HF", CompositeFunctionals::XCFUNCTIONALS::HF)
      .value("B97_D", CompositeFunctionals::XCFUNCTIONALS::B97_D)
      .value("WB97", CompositeFunctionals::XCFUNCTIONALS::WB97)
      .value("WB97X", CompositeFunctionals::XCFUNCTIONALS::WB97X)
      .value("WB97X_D", CompositeFunctionals::XCFUNCTIONALS::WB97X_D)
      .value("WB97X_V", CompositeFunctionals::XCFUNCTIONALS::WB97X_V)
      .export_values();
  py::enum_<CompositeFunctionals::KINFUNCTIONALS>(spy, "KINFUNCTIONALS")
      .value("NONE", CompositeFunctionals::KINFUNCTIONALS::NONE)
      .value("TF", CompositeFunctionals::KINFUNCTIONALS::TF)
      .value("PW91K", CompositeFunctionals::KINFUNCTIONALS::PW91K)
      .value("LLP91K", CompositeFunctionals::KINFUNCTIONALS::LLP91K)
      .value("LLP91KS", CompositeFunctionals::KINFUNCTIONALS::LLP91KS)
      .value("PBE2K", CompositeFunctionals::KINFUNCTIONALS::PBE2K)
      .value("PBE2KS", CompositeFunctionals::KINFUNCTIONALS::PBE2KS)
      .value("PBE3K", CompositeFunctionals::KINFUNCTIONALS::PBE3K)
      .value("PBE4K", CompositeFunctionals::KINFUNCTIONALS::PBE4K)
      .value("E2000K", CompositeFunctionals::KINFUNCTIONALS::E2000K)
      .export_values();
}
