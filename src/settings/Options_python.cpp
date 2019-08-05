/**
 * @file Options_python.cpp
 *
 * @date Apr 27, 2016
 * @author: Jan Unsleber
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


/* Include Serenity Internal Headers */
#include "settings/Options.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

using namespace Serenity;
using namespace Options;
namespace py = pybind11;

void export_Options(py::module &spy) {

py::enum_<SCF_MODES>(spy,"SCF_MODES")
            .value("RESTRICTED", SCF_MODES::RESTRICTED)
            .value("UNRESTRICTED", SCF_MODES::UNRESTRICTED)
            .export_values();
py::enum_<ELECTRONIC_STRUCTURE_THEORIES>(spy,"ELECTRONIC_STRUCTURE_THEORIES")
            .value("HF", ELECTRONIC_STRUCTURE_THEORIES::HF)
            .value("DFT", ELECTRONIC_STRUCTURE_THEORIES::DFT)
            .export_values();
py::enum_<FUNCTIONALS>(spy,"FUNCTIONALS")
        .value("NONE", FUNCTIONALS::NONE)
        .value("LDA", FUNCTIONALS::LDA)
        .value("BLYP", FUNCTIONALS::BLYP)
        .value("PBE", FUNCTIONALS::PBE)
        .value("BP86", FUNCTIONALS::BP86)
        .value("KT1", FUNCTIONALS::KT1)
        .value("KT2", FUNCTIONALS::KT2)
        .value("KT3", FUNCTIONALS::KT3)
        .value("LDAERF", FUNCTIONALS::LDAERF)
        .value("PBE0", FUNCTIONALS::PBE0)
        .value("B3LYP", FUNCTIONALS::B3LYP)
        .value("B3LYP_G", FUNCTIONALS::B3LYP_G)
        .value("B3P86", FUNCTIONALS::B3P86)
        .value("B3P86_G", FUNCTIONALS::B3P86_G)
        .value("BPW91", FUNCTIONALS::BPW91)
        .value("B97", FUNCTIONALS::B97)
        .value("B97_1", FUNCTIONALS::B97_1)
        .value("B97_2", FUNCTIONALS::B97_2)
        .value("CAMB3LYP", FUNCTIONALS::CAMB3LYP)
        .value("VWN5", FUNCTIONALS::VWN5)
        .value("VWN3", FUNCTIONALS::VWN3)
        .value("SLATER", FUNCTIONALS::SLATER)
        .value("OLYP", FUNCTIONALS::OLYP)
        .value("TF", FUNCTIONALS::TF)
        .value("TW", FUNCTIONALS::TW)
        .value("LCBLYP", FUNCTIONALS::LCBLYP)
        .value("PW91K", FUNCTIONALS::PW91K)
        .value("PW91", FUNCTIONALS::PW91)
        .value("B2PLYP", FUNCTIONALS::B2PLYP)
        .value("B2KPLYP", FUNCTIONALS::B2KPLYP)
        .value("B2TPLYP", FUNCTIONALS::B2TPLYP)
        .value("B2GPPLYP", FUNCTIONALS::B2GPPLYP)
        .value("ROB2PLYP", FUNCTIONALS::ROB2PLYP)
        .value("B2PIPLYP", FUNCTIONALS::B2PIPLYP)
        .value("B2PPW91", FUNCTIONALS::B2PPW91)
        .value("DSDBLYP", FUNCTIONALS::DSDBLYP)
        .value("DUT", FUNCTIONALS::DUT)
        .value("PUT", FUNCTIONALS::PUT)
        .value("DSDPBEP86", FUNCTIONALS::DSDPBEP86)
        .value("PWPB95", FUNCTIONALS::PWPB95)
        .export_values();
py::enum_<XCFUNCTIONALS>(spy,"XCFUNCTIONALS")
        .value("NONE", XCFUNCTIONALS::NONE)
        .value("LDA", XCFUNCTIONALS::LDA)
        .value("BLYP", XCFUNCTIONALS::BLYP)
        .value("PBE", XCFUNCTIONALS::PBE)
        .value("BP86", XCFUNCTIONALS::BP86)
        .value("KT1", XCFUNCTIONALS::KT1)
        .value("KT2", XCFUNCTIONALS::KT2)
        .value("KT3", XCFUNCTIONALS::KT3)
        .value("LDAERF", XCFUNCTIONALS::LDAERF)
        .value("PBE0", XCFUNCTIONALS::PBE0)
        .value("B3LYP", XCFUNCTIONALS::B3LYP)
        .value("B3LYP_G", XCFUNCTIONALS::B3LYP_G)
        .value("B3P86", XCFUNCTIONALS::B3P86)
        .value("B3P86_G", XCFUNCTIONALS::B3P86_G)
        .value("BPW91", XCFUNCTIONALS::BPW91)
        .value("B97", XCFUNCTIONALS::B97)
        .value("B97_1", XCFUNCTIONALS::B97_1)
        .value("B97_2", XCFUNCTIONALS::B97_2)
        .value("LCBLYP", XCFUNCTIONALS::LCBLYP)
        .value("CAMB3LYP", XCFUNCTIONALS::CAMB3LYP)
        .value("VWN5", XCFUNCTIONALS::VWN5)
        .value("VWN3", XCFUNCTIONALS::VWN3)
        .value("SLATER", XCFUNCTIONALS::SLATER)
        .value("OLYP", XCFUNCTIONALS::OLYP)
        .value("PW91", XCFUNCTIONALS::PW91)
        .value("B2PLYP", XCFUNCTIONALS::B2PLYP)
        .value("B2KPLYP", XCFUNCTIONALS::B2KPLYP)
        .value("B2TPLYP", XCFUNCTIONALS::B2TPLYP)
        .value("B2GPPLYP", XCFUNCTIONALS::B2GPPLYP)
        .value("ROB2PLYP", XCFUNCTIONALS::ROB2PLYP)
        .value("B2PIPLYP", XCFUNCTIONALS::B2PIPLYP)
        .value("B2PPW91", XCFUNCTIONALS::B2PPW91)
        .value("DSDBLYP", XCFUNCTIONALS::DSDBLYP)
        .value("DUT", XCFUNCTIONALS::DUT)
        .value("PUT", XCFUNCTIONALS::PUT)
        .value("DSDPBEP86", XCFUNCTIONALS::DSDPBEP86)
        .value("PWPB95", XCFUNCTIONALS::PWPB95)
        .export_values();
py::enum_<KINFUNCTIONALS>(spy,"KINFUNCTIONALS")
        .value("NONE", KINFUNCTIONALS::NONE)
        .value("TF", KINFUNCTIONALS::TF)
        .value("TW", KINFUNCTIONALS::TW)
        .value("PW91K", KINFUNCTIONALS::PW91K)
        .export_values();
py::enum_<DFT_DISPERSION_CORRECTIONS>(spy,"DFT_DISPERSION_CORRECTIONS")
          .value("NONE", DFT_DISPERSION_CORRECTIONS::NONE)
          .value("D3", DFT_DISPERSION_CORRECTIONS::D3)
          .value("D3ABC", DFT_DISPERSION_CORRECTIONS::D3ABC)
          .value("D3BJ", DFT_DISPERSION_CORRECTIONS::D3BJ)
          .value("D3BJABC", DFT_DISPERSION_CORRECTIONS::D3BJABC)
          .export_values();
py::enum_<INITIAL_GUESSES>(spy,"INITIAL_GUESSES")
    .value("HCORE", INITIAL_GUESSES::H_CORE)
    .value("H_CORE", INITIAL_GUESSES::H_CORE)
    .value("EHT", INITIAL_GUESSES::EHT)
    .value("ATOM_DENS", INITIAL_GUESSES::ATOM_DENS)
    .value("ATOMDENS", INITIAL_GUESSES::ATOM_DENS)
    .value("ATOM_SCF", INITIAL_GUESSES::ATOM_SCF)
    .value("ATOMSCF", INITIAL_GUESSES::ATOM_SCF)
    .export_values();
py::enum_<DAMPING_ALGORITHMS>(spy,"DAMPING_ALGORITHMS")
  .value("NONE",DAMPING_ALGORITHMS::NONE)
  .value("STATIC",DAMPING_ALGORITHMS::STATIC)
  .value("SERIES",DAMPING_ALGORITHMS::SERIES)
  .value("DYNAMIC",DAMPING_ALGORITHMS::DYNAMIC)
  .export_values();
py::enum_<SIGMAVECTOR_METHODS>(spy,"SIGMAVECTOR_METHODS")
  .value("APLUSB",SIGMAVECTOR_METHODS::APLUSB)
  .value("AMINUSB",SIGMAVECTOR_METHODS::AMINUSB)
  .export_values();
py::enum_<MULTIPLICITY>(spy,"MULTIPLICITY")
  .value("SINGLET",MULTIPLICITY::SINGLET)
  .value("TRIPLET",MULTIPLICITY::TRIPLET)
  .export_values();
py::enum_<BASIS_PURPOSES>(spy,"BASIS_PURPOSES")
    .value("DEFAULT", BASIS_PURPOSES::DEFAULT)
    .value("AUX_COULOMB", BASIS_PURPOSES::AUX_COULOMB)
    .value("MINBAS", BASIS_PURPOSES::MINBAS)
    .value("HUECKEL", BASIS_PURPOSES::HUECKEL)
    .value("AUX_CORREL", BASIS_PURPOSES::AUX_CORREL)
    .export_values();
py::enum_<DENS_FITS>(spy,"DENS_FITS")
  .value("RI",DENS_FITS::RI)
  .value("NORI",DENS_FITS::NONE)
  .value("NONE",DENS_FITS::NONE)
  .export_values();
py::enum_<RADIAL_GRID_TYPES>(spy,"RADIAL_GRID_TYPES")
    .value("BECKE", RADIAL_GRID_TYPES::BECKE)
    .value("HANDY", RADIAL_GRID_TYPES::HANDY)
    .value("AHLRICHS", RADIAL_GRID_TYPES::AHLRICHS)
    .value("KNOWLES", RADIAL_GRID_TYPES::KNOWLES)
    .value("EQUI", RADIAL_GRID_TYPES::EQUI)
    .export_values()
    ;
py::enum_<SPHERICAL_GRID_TYPES>(spy,"SPHERICAL_GRID_TYPES")
    .value("LEBEDEV", SPHERICAL_GRID_TYPES::LEBEDEV)
    .export_values()
    ;
py::enum_<GRID_TYPES>(spy,"GRID_TYPES")
    .value("BECKE", GRID_TYPES::BECKE)
    .value("VORONOI", GRID_TYPES::VORONOI)
    .export_values()
    ;
py::enum_<GRID_PURPOSES>(spy,"GRID_PURPOSES")
    .value("DEFAULT", GRID_PURPOSES::DEFAULT)
    .value("SMALL", GRID_PURPOSES::SMALL)
    .value("PLOT", GRID_PURPOSES::PLOT)
    .export_values();
py::enum_<OPTIMIZATION_ALGORITHMS>(spy,"OPTIMIZATION_ALGORITHMS")
    .value("SD", OPTIMIZATION_ALGORITHMS::SD)
    .value("BFGS", OPTIMIZATION_ALGORITHMS::BFGS)
    .export_values()
    ;
py::enum_<GRADIENT_TYPES>(spy,"GRADIENT_TYPES")
    .value("NUMERICAL", GRADIENT_TYPES::NUMERICAL)
    .value("ANALYTICAL", GRADIENT_TYPES::ANALYTICAL)
    .export_values()
    ;
py::enum_<GEOMETRY_OPTIMIZATION_TYPES>(spy,"GEOMETRY_OPTIMIZATION_TYPES")
    .value("GROUNDSTATE", GEOMETRY_OPTIMIZATION_TYPES::GROUNDSTATE)
    .value("TS", GEOMETRY_OPTIMIZATION_TYPES::TS)
    .export_values()
    ;
py::enum_<HESSIAN_TYPES>(spy,"HESSIAN_TYPES")
    .value("NUMERICAL", HESSIAN_TYPES::NUMERICAL)
    .value("ANALYTICAL", HESSIAN_TYPES::ANALYTICAL)
    .export_values()
    ;
py::enum_<ORBITAL_LOCALIZATION_ALGORITHMS>(spy,"ORBITAL_LOCALIZATION_ALGORITHMS")
    .value("PM", ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY)
    .value("FB", ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS)
    .value("IAO", ORBITAL_LOCALIZATION_ALGORITHMS::IAO)
    .value("IBO", ORBITAL_LOCALIZATION_ALGORITHMS::IBO)
    .value("ER", ORBITAL_LOCALIZATION_ALGORITHMS::EDMISTON_RUEDENBERG)
    .value("NO", ORBITAL_LOCALIZATION_ALGORITHMS::NON_ORTHOGONAL)
    .export_values()
    ;
py::enum_<KIN_EMBEDDING_MODES>(spy,"KIN_EMBEDDING_MODES")
    .value("NONE", KIN_EMBEDDING_MODES::NONE)
    .value("NADDFUNC", KIN_EMBEDDING_MODES::NADD_FUNC)
    .value("LEVELSHIFT", KIN_EMBEDDING_MODES::LEVELSHIFT)
    .value("HUZINAGA", KIN_EMBEDDING_MODES::HUZINAGA)
    .value("RECONSTRUCTION", KIN_EMBEDDING_MODES::RECONSTRUCTION)
    .export_values()
    ;
py::enum_<CC_LEVEL>(spy,"CC_LEVEL")
    .value("CCSD",CC_LEVEL::CCSD)
    .value("CCSD_T",CC_LEVEL::CCSD_T)
    .export_values()
    ;
py::enum_<REFERENCE_POINT>(spy,"REFERENCE_POINT")
    .value("CENTEROFMASS",REFERENCE_POINT::COM)
    .value("ORIGIN",REFERENCE_POINT::ORIGIN)
    .export_values()
    ;
}



