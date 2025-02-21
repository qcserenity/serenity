/**
 * @file Options_python.cpp
 *
 * @date Apr 27, 2016
 * @author: Jan Unsleber
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
#include "settings/BasisOptions.h"
#include "settings/CorrelatedMethodsOptions.h"
#include "settings/DFTOptions.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/EmbeddingOptions.h"
#include "settings/GeometryOptions.h"
#include "settings/GridOptions.h"
#include "settings/LRSCFOptions.h"
#include "settings/LocalizationOptions.h"
#include "settings/MBPTOptions.h"
#include "settings/MiscOptions.h"
#include "settings/SCFOptions.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

using namespace Serenity::Options;
namespace py = pybind11;

void export_Options(py::module& spy) {
  py::enum_<SCF_MODES>(spy, "SCF_MODES")
      .value("RESTRICTED", SCF_MODES::RESTRICTED)
      .value("UNRESTRICTED", SCF_MODES::UNRESTRICTED)
      .export_values();
  py::enum_<ELECTRONIC_STRUCTURE_THEORIES>(spy, "ELECTRONIC_STRUCTURE_THEORIES")
      .value("HF", ELECTRONIC_STRUCTURE_THEORIES::HF)
      .value("DFT", ELECTRONIC_STRUCTURE_THEORIES::DFT)
      .export_values();
  py::enum_<DFT_DISPERSION_CORRECTIONS>(spy, "DFT_DISPERSION_CORRECTIONS")
      .value("NONE", DFT_DISPERSION_CORRECTIONS::NONE)
      .value("D3", DFT_DISPERSION_CORRECTIONS::D3)
      .value("D3ABC", DFT_DISPERSION_CORRECTIONS::D3ABC)
      .value("D3BJ", DFT_DISPERSION_CORRECTIONS::D3BJ)
      .value("D3BJABC", DFT_DISPERSION_CORRECTIONS::D3BJABC)
      .export_values();
  py::enum_<INITIAL_GUESSES>(spy, "INITIAL_GUESSES")
      .value("HCORE", INITIAL_GUESSES::H_CORE)
      .value("H_CORE", INITIAL_GUESSES::H_CORE)
      .value("EHT", INITIAL_GUESSES::EHT)
      .value("ATOM_SCF", INITIAL_GUESSES::ATOM_SCF)
      .value("ATOMSCF", INITIAL_GUESSES::ATOM_SCF)
      .value("ATOM_SCF_INPLACE", INITIAL_GUESSES::ATOM_SCF_INPLACE)
      .value("ATOMSCF_INPLACE", INITIAL_GUESSES::ATOM_SCF_INPLACE)
      .value("SAP", INITIAL_GUESSES::SAP)
      .export_values();
  py::enum_<DAMPING_ALGORITHMS>(spy, "DAMPING_ALGORITHMS")
      .value("NONE", DAMPING_ALGORITHMS::NONE)
      .value("STATIC", DAMPING_ALGORITHMS::STATIC)
      .value("SERIES", DAMPING_ALGORITHMS::SERIES)
      .value("DYNAMIC", DAMPING_ALGORITHMS::DYNAMIC)
      .export_values();
  py::enum_<BASIS_PURPOSES>(spy, "BASIS_PURPOSES")
      .value("DEFAULT", BASIS_PURPOSES::DEFAULT)
      .value("AUX_COULOMB", BASIS_PURPOSES::AUX_COULOMB)
      .value("MINBAS", BASIS_PURPOSES::MINBAS)
      .value("HUECKEL", BASIS_PURPOSES::HUECKEL)
      .value("IAO_LOCALIZATION", BASIS_PURPOSES::IAO_LOCALIZATION)
      .value("SCF_DENS_GUESS", BASIS_PURPOSES::SCF_DENS_GUESS)
      .value("AUX_CORREL", BASIS_PURPOSES::AUX_CORREL)
      .value("ATOMIC_CHOLESKY", BASIS_PURPOSES::ATOMIC_CHOLESKY)
      .value("ATOMIC_COMPACT_CHOLESKY", BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY)
      .value("ATOMIC_CHOLESKY_ERF", BASIS_PURPOSES::ERF_ATOMIC_CHOLESKY)
      .value("ATOMIC_COMPACT_CHOLESKY_ERF", BASIS_PURPOSES::ERF_ATOMIC_COMPACT_CHOLESKY)
      .value("AUX_JK", BASIS_PURPOSES::AUX_JK)
      .export_values();
  py::enum_<AUX_BASIS_PURPOSES>(spy, "AUX_BASIS_PURPOSES")
      .value("COULOMB", AUX_BASIS_PURPOSES::COULOMB)
      .value("EXCHANGE", AUX_BASIS_PURPOSES::EXCHANGE)
      .value("LREXCHANGE", AUX_BASIS_PURPOSES::LREXCHANGE)
      .value("CORRELATION", AUX_BASIS_PURPOSES::CORRELATION)
      .export_values();
  py::enum_<DENS_FITS>(spy, "DENS_FITS")
      .value("RI", DENS_FITS::RI)
      .value("NORI", DENS_FITS::NONE)
      .value("NONE", DENS_FITS::NONE)
      .value("CD", DENS_FITS::CD)
      .value("ACD", DENS_FITS::ACD)
      .value("ACCD", DENS_FITS::ACCD)
      .export_values();
  py::enum_<EXTEND_ACD>(spy, "EXTEND_ACD")
      .value("NONE", EXTEND_ACD::NONE)
      .value("SIMPLE", EXTEND_ACD::SIMPLE)
      .value("FIRST", EXTEND_ACD::FIRST)
      .value("COMPLETE", EXTEND_ACD::COMPLETE)
      .export_values();
  py::enum_<RADIAL_GRID_TYPES>(spy, "RADIAL_GRID_TYPES")
      .value("BECKE", RADIAL_GRID_TYPES::BECKE)
      .value("HANDY", RADIAL_GRID_TYPES::HANDY)
      .value("AHLRICHS", RADIAL_GRID_TYPES::AHLRICHS)
      .value("KNOWLES", RADIAL_GRID_TYPES::KNOWLES)
      .value("EQUI", RADIAL_GRID_TYPES::EQUI)
      .value("EQUIDISTAND", RADIAL_GRID_TYPES::EQUI)
      .export_values();
  py::enum_<SPHERICAL_GRID_TYPES>(spy, "SPHERICAL_GRID_TYPES").value("LEBEDEV", SPHERICAL_GRID_TYPES::LEBEDEV).export_values();
  py::enum_<GRID_TYPES>(spy, "GRID_TYPES")
      .value("BECKE", GRID_TYPES::BECKE)
      .value("VORONOI", GRID_TYPES::VORONOI)
      .value("SSF", GRID_TYPES::SSF)
      .export_values();
  py::enum_<GRID_PURPOSES>(spy, "GRID_PURPOSES")
      .value("DEFAULT", GRID_PURPOSES::DEFAULT)
      .value("SMALL", GRID_PURPOSES::SMALL)
      .value("PLOT", GRID_PURPOSES::PLOT)
      .export_values();
  py::enum_<OPTIMIZATION_ALGORITHMS>(spy, "OPTIMIZATION_ALGORITHMS")
      .value("SQNM", OPTIMIZATION_ALGORITHMS::SQNM)
      .value("BFGS", OPTIMIZATION_ALGORITHMS::BFGS)
      .export_values();
  py::enum_<GRADIENT_TYPES>(spy, "GRADIENT_TYPES")
      .value("NUMERICAL", GRADIENT_TYPES::NUMERICAL)
      .value("ANALYTICAL", GRADIENT_TYPES::ANALYTICAL)
      .export_values();
  py::enum_<GEOMETRY_OPTIMIZATION_TYPES>(spy, "GEOMETRY_OPTIMIZATION_TYPES")
      .value("GROUNDSTATE", GEOMETRY_OPTIMIZATION_TYPES::GROUNDSTATE)
      .value("TS", GEOMETRY_OPTIMIZATION_TYPES::TS)
      .export_values();
  py::enum_<HESSIAN_TYPES>(spy, "HESSIAN_TYPES")
      .value("NUMERICAL", HESSIAN_TYPES::NUMERICAL)
      .value("ANALYTICAL", HESSIAN_TYPES::ANALYTICAL)
      .export_values();
  py::enum_<ORBITAL_LOCALIZATION_ALGORITHMS>(spy, "ORBITAL_LOCALIZATION_ALGORITHMS")
      .value("PM", ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY)
      .value("FB", ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS)
      .value("IAO", ORBITAL_LOCALIZATION_ALGORITHMS::IAO)
      .value("IBO", ORBITAL_LOCALIZATION_ALGORITHMS::IBO)
      .value("ER", ORBITAL_LOCALIZATION_ALGORITHMS::EDMISTON_RUEDENBERG)
      .value("NO", ORBITAL_LOCALIZATION_ALGORITHMS::NON_ORTHOGONAL)
      .value("ALIGN", ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN)
      .value("NONE", ORBITAL_LOCALIZATION_ALGORITHMS::NONE)
      .export_values();
  py::enum_<POPULATION_ANALYSIS_ALGORITHMS>(spy, "POPULATION_ANALYSIS_ALGORITHMS")
      .value("MUL", POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN)
      .value("HIRSHFELD", POPULATION_ANALYSIS_ALGORITHMS::HIRSHFELD)
      .value("IAO", POPULATION_ANALYSIS_ALGORITHMS::IAO)
      .value("IAOSHELL", POPULATION_ANALYSIS_ALGORITHMS::IAOShell)
      .value("BECKE", POPULATION_ANALYSIS_ALGORITHMS::BECKE)
      .export_values();
  py::enum_<DOS_SETTINGS>(spy, "DOS_SETTINGS")
      .value("LOOSE", DOS_SETTINGS::LOOSE)
      .value("NORMAL", DOS_SETTINGS::NORMAL)
      .value("TIGHT", DOS_SETTINGS::TIGHT)
      .value("VERYTIGHT", DOS_SETTINGS::VERY_TIGHT)
      .value("VERY_TIGHT", DOS_SETTINGS::VERY_TIGHT)
      .value("EXTREME", DOS_SETTINGS::EXTREME)
      .value("SPREAD", DOS_SETTINGS::SPREAD)
      .export_values();
  py::enum_<KIN_EMBEDDING_MODES>(spy, "KIN_EMBEDDING_MODES")
      .value("NONE", KIN_EMBEDDING_MODES::NONE)
      .value("NADDFUNC", KIN_EMBEDDING_MODES::NADD_FUNC)
      .value("LEVELSHIFT", KIN_EMBEDDING_MODES::LEVELSHIFT)
      .value("HUZINAGA", KIN_EMBEDDING_MODES::HUZINAGA)
      .value("HOFFMANN", KIN_EMBEDDING_MODES::HOFFMANN)
      .value("RECONSTRUCTION", KIN_EMBEDDING_MODES::RECONSTRUCTION)
      .value("FERMI_SHIFTED_HUZINAGA", KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA)
      .value("FERMI", KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA)
      .value("LOEWDIN", KIN_EMBEDDING_MODES::LOEWDIN)
      .value("ALMO", KIN_EMBEDDING_MODES::ALMO)
      .export_values();
  py::enum_<EMBEDDING_SCHEME>(spy, "EMBEDDING_SCHEME")
      .value("NONE", EMBEDDING_SCHEME::NONE)
      .value("ISOLATED", EMBEDDING_SCHEME::ISOLATED)
      .value("FDE", EMBEDDING_SCHEME::FDE)
      .value("FAT", EMBEDDING_SCHEME::FAT)
      .export_values();
  py::enum_<CC_LEVEL>(spy, "CC_LEVEL")
      .value("CCSD", CC_LEVEL::CCSD)
      .value("CCSD_T", CC_LEVEL::CCSD_T)
      .value("CCSD(T)", CC_LEVEL::CCSD_T)
      .value("DLPNO-CCSD", CC_LEVEL::DLPNO_CCSD)
      .value("DLPNO-CCSD(T0)", CC_LEVEL::DLPNO_CCSD_T0)
      .export_values();
  py::enum_<MP2_TYPES>(spy, "MP2_TYPES")
      .value("AO", MP2_TYPES::AO)
      .value("DF", MP2_TYPES::DF)
      .value("LOCAL", MP2_TYPES::LOCAL)
      .value("LT", MP2_TYPES::LT)
      .export_values();
  py::enum_<PNO_SETTINGS>(spy, "PNO_SETTINGS")
      .value("LOOSE", PNO_SETTINGS::LOOSE)
      .value("NORMAL", PNO_SETTINGS::NORMAL)
      .value("TIGHT", PNO_SETTINGS::TIGHT)
      .export_values();
  py::enum_<PNO_METHOD>(spy, "PNO_METHOD")
      .value("DLPNO-MP2", PNO_METHOD::DLPNO_MP2)
      .value("LMP2", PNO_METHOD::DLPNO_MP2)
      .value("DLPNO-CCSD", PNO_METHOD::DLPNO_CCSD)
      .value("DLPNO-CCSD(T0)", PNO_METHOD::DLPNO_CCSD_T0)
      .value("SC-MP2", PNO_METHOD::SC_MP2)
      .value("NONE", PNO_METHOD::NONE)
      .value("HF", PNO_METHOD::NONE)
      .export_values();
  py::enum_<SYSTEM_SPLITTING_ALGORITHM>(spy, "SYSTEM_SPLITTING_ALGORITHM")
      .value("ENFORCE_CHARGES", SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES)
      .value("ENFORCECHARGES", SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES)
      .value("BEST_MATCH", SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH)
      .value("BESTMATCH", SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH)
      .value("POPULATION_THRESHOLD", SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD)
      .value("POPULATION", SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD)
      .value("SPADE", SYSTEM_SPLITTING_ALGORITHM::SPADE)
      .value("SPADE_ENFORCE_CHARGES", SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES)
      .value("SPADE_ENFORCECHARGES", SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES)
      .value("SPADEENFORCECHARGES", SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES)
      .export_values();
  py::enum_<BASIS_SET_TRUNCATION_ALGORITHMS>(spy, "BASIS_SET_TRUNCATION_ALGORITHMS")
      .value("NONE", BASIS_SET_TRUNCATION_ALGORITHMS::NONE)
      .value("NETPOP", BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION)
      .value("NETPOPULATION", BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION)
      .value("PRIMITIVENETPOP", BASIS_SET_TRUNCATION_ALGORITHMS::PRIMITIVE_NET_POPULATION)
      .export_values();
  py::enum_<GAUGE_ORIGIN>(spy, "GAUGE_ORIGIN")
      .value("CENTEROFMASS", GAUGE_ORIGIN::COM)
      .value("ORIGIN", GAUGE_ORIGIN::ORIGIN)
      .export_values();
  py::enum_<ORBITAL_FILE_TYPES>(spy, "ORBITAL_FILE_TYPES")
      .value("SERENITY", ORBITAL_FILE_TYPES::SERENITY)
      .value("TURBOMOLE", ORBITAL_FILE_TYPES::TURBOMOLE)
      .value("MOLPRO", ORBITAL_FILE_TYPES::MOLPRO)
      .value("MOLCAS", ORBITAL_FILE_TYPES::MOLCAS)
      .value("MOLDEN", ORBITAL_FILE_TYPES::MOLDEN)
      .export_values();
  py::enum_<LRSCF_TYPE>(spy, "LRSCF_TYPE")
      .value("ISOLATED", LRSCF_TYPE::ISOLATED)
      .value("ISO", LRSCF_TYPE::ISOLATED)
      .value("UNCOUPLED", LRSCF_TYPE::UNCOUPLED)
      .value("FDEU", LRSCF_TYPE::UNCOUPLED)
      .value("COUPLED", LRSCF_TYPE::COUPLED)
      .value("FDEC", LRSCF_TYPE::COUPLED)
      .export_values();
  py::enum_<INTEGRAL_TYPE>(spy, "INTEGRAL_TYPE")
      .value("NUMERICAL", INTEGRAL_TYPE::NUMERICAL)
      .value("ANALYTICAL", INTEGRAL_TYPE::ANALYTICAL)
      .export_values();
  py::enum_<GAUGE>(spy, "GAUGE").value("LENGTH", GAUGE::LENGTH).value("VELOCITY", GAUGE::VELOCITY).export_values();
  py::enum_<LR_METHOD>(spy, "LR_METHOD")
      .value("CIS", LR_METHOD::TDA)
      .value("TDA", LR_METHOD::TDA)
      .value("RPA", LR_METHOD::TDDFT)
      .value("TDDFT", LR_METHOD::TDDFT)
      .value("CC2", LR_METHOD::CC2)
      .value("CISDINF", LR_METHOD::CISDINF)
      .value("CISD", LR_METHOD::CISD)
      .value("ADC2", LR_METHOD::ADC2)
      .export_values();
  py::enum_<STABILITY_ANALYSIS>(spy, "STABILITY_ANALYSIS")
      .value("NONE", STABILITY_ANALYSIS::NONE)
      .value("REAL", STABILITY_ANALYSIS::REAL)
      .value("NONREAL", STABILITY_ANALYSIS::NONREAL)
      .value("SPINFLIP", STABILITY_ANALYSIS::SPINFLIP)
      .export_values();
  py::enum_<Serenity::CompositeFunctionals::IMPLEMENTATIONS>(spy, "XC_IMPLEMENTATIONS")
      .value("XCFUN", Serenity::CompositeFunctionals::IMPLEMENTATIONS::XCFUN)
      .value("LIBXC", Serenity::CompositeFunctionals::IMPLEMENTATIONS::LIBXC)
      .export_values();
  py::enum_<MBPT>(spy, "MBPT").value("GW", MBPT::GW).value("RPA", MBPT::RPA).export_values();
  py::enum_<GWALGORITHM>(spy, "GWALGORITHM")
      .value("CD", GWALGORITHM::CD)
      .value("AC", GWALGORITHM::AC)
      .value("ANALYTIC", GWALGORITHM::ANALYTIC)
      .export_values();
}
