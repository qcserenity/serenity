/**
 * @file Settings_python.cpp
 *
 * @date Apr 25, 2016
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
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

// Settings deepcopySettings(const Settings settings){return Settings(settings);};

void export_Settings(py::module& spy) {
  py::class_<DFT>(spy, "DFTSettings", "DFT settings, all set to their default values.")
      .def_readwrite("functional", &DFT::functional)
      .def_readwrite("dispersion", &DFT::dispersion);

  py::class_<CUSTOMFUNCTIONAL>(spy, "CustomFunctionalSettings", "Custom functional settings, all set to their default values.")
      .def_readwrite("impl", &CUSTOMFUNCTIONAL::impl)
      .def_readwrite("basicFunctionals", &CUSTOMFUNCTIONAL::basicFunctionals)
      .def_readwrite("mixingFactors", &CUSTOMFUNCTIONAL::mixingFactors)
      .def_readwrite("hfExchangeRatio", &CUSTOMFUNCTIONAL::hfExchangeRatio)
      .def_readwrite("hfCorrelRatio", &CUSTOMFUNCTIONAL::hfCorrelRatio)
      .def_readwrite("lrExchangeRatio", &CUSTOMFUNCTIONAL::lrExchangeRatio)
      .def_readwrite("mu", &CUSTOMFUNCTIONAL::mu)
      .def_readwrite("ssScaling", &CUSTOMFUNCTIONAL::ssScaling)
      .def_readwrite("osScaling", &CUSTOMFUNCTIONAL::osScaling);

  py::class_<SCF>(spy, "SCFSettings", "SCF settings, all set to their default values.")
      .def_readwrite("initialguess", &SCF::initialguess)
      .def_readwrite("maxCycles", &SCF::maxCycles)
      .def_readwrite("writeRestart", &SCF::writeRestart)
      .def_readwrite("energyThreshold", &SCF::energyThreshold)
      .def_readwrite("rmsdThreshold", &SCF::rmsdThreshold)
      .def_readwrite("damping", &SCF::damping)
      .def_readwrite("seriesDampingStart", &SCF::seriesDampingStart)
      .def_readwrite("seriesDampingEnd", &SCF::seriesDampingEnd)
      .def_readwrite("seriesDampingStep", &SCF::seriesDampingStep)
      .def_readwrite("seriesDampingInitialSteps", &SCF::seriesDampingInitialSteps)
      .def_readwrite("staticDampingFactor", &SCF::staticDampingFactor)
      .def_readwrite("endDampErr", &SCF::endDampErr)
      .def_readwrite("useLevelshift", &SCF::useLevelshift)
      .def_readwrite("useOffDiagLevelshift", &SCF::useOffDiagLevelshift)
      .def_readwrite("minimumLevelshift", &SCF::minimumLevelshift)
      .def_readwrite("diisFlush", &SCF::diisFlush)
      .def_readwrite("diisStartError", &SCF::diisStartError)
      .def_readwrite("diisMaxStore", &SCF::diisMaxStore)
      .def_readwrite("diisThreshold", &SCF::diisThreshold)
      .def_readwrite("canOrthThreshold", &SCF::canOrthThreshold)
      .def_readwrite("useADIIS", &SCF::useADIIS)
      .def_readwrite("allowNotConverged", &SCF::allowNotConverged)
      .def_readwrite("rohf", &SCF::rohf)
      .def_readwrite("suhfLambda", &SCF::suhfLambda)
      .def_readwrite("degeneracyThreshold", &SCF::degeneracyThreshold);

  py::class_<BASIS>(spy, "BasisSettings", "Basis settings, all set to their default values.")
      .def_readwrite("label", &BASIS::label)
      .def_readwrite("auxJLabel", &BASIS::auxJLabel)
      .def_readwrite("auxJKLabel", &BASIS::auxJKLabel)
      .def_readwrite("auxCLabel", &BASIS::auxCLabel)
      .def_readwrite("makeSphericalBasis", &BASIS::makeSphericalBasis)
      .def_readwrite("integralThreshold", &BASIS::integralThreshold)
      .def_readwrite("integralIncrementThresholdStart", &BASIS::integralIncrementThresholdStart)
      .def_readwrite("integralIncrementThresholdEnd", &BASIS::integralIncrementThresholdEnd)
      .def_readwrite("incrementalSteps", &BASIS::incrementalSteps)
      .def_readwrite("basisLibPath", &BASIS::basisLibPath)
      .def_readwrite("firstECP", &BASIS::firstECP)
      .def_readwrite("densFitJ", &BASIS::densFitJ)
      .def_readwrite("densFitK", &BASIS::densFitK)
      .def_readwrite("densFitLRK", &BASIS::densFitLRK)
      .def_readwrite("densFitCorr", &BASIS::densFitCorr)
      .def_readwrite("cdThreshold", &BASIS::cdThreshold)
      .def_readwrite("extendSphericalACDShells", &BASIS::extendSphericalACDShells)
      .def_readwrite("intCondition", &BASIS::intCondition)
      .def_readwrite("secondCD", &BASIS::secondCD)
      .def_readwrite("cdOffset", &BASIS::cdOffset);

  py::class_<GRID>(spy, "GridSettings", "Grid settings, all set to their default values.")
      .def_readwrite("gridType", &GRID::gridType)
      .def_readwrite("radialGridType", &GRID::radialGridType)
      .def_readwrite("sphericalGridType", &GRID::sphericalGridType)
      .def_readwrite("blocksize", &GRID::blocksize)
      .def_readwrite("accuracy", &GRID::accuracy)
      .def_readwrite("smallGridAccuracy", &GRID::smallGridAccuracy)
      .def_readwrite("blockAveThreshold", &GRID::blockAveThreshold)
      .def_readwrite("basFuncRadialThreshold", &GRID::basFuncRadialThreshold)
      .def_readwrite("weightThreshold", &GRID::weightThreshold)
      .def_readwrite("smoothing", &GRID::smoothing)
      .def_readwrite("gridPointSorting", &GRID::gridPointSorting);

  py::class_<EFIELD>(spy, "EFIELD", "Electric Field settings, all set to their default values.")
      .def_readwrite("use", &EFIELD::use)
      .def_readwrite("analytical", &EFIELD::analytical)
      .def_readwrite("pos1", &EFIELD::pos1)
      .def_readwrite("pos2", &EFIELD::pos2)
      .def_readwrite("distance", &EFIELD::distance)
      .def_readwrite("nRings", &EFIELD::nRings)
      .def_readwrite("radius", &EFIELD::radius)
      .def_readwrite("fieldStrength", &EFIELD::fieldStrength)
      .def_readwrite("nameOutput", &EFIELD::nameOutput);

  py::class_<Settings>(spy, "Settings")
      .def(py::init<>())
      .def_readwrite("name", &Settings::name)
      .def_readwrite("identifier", &Settings::identifier)
      .def_readwrite("path", &Settings::path)
      .def_readwrite("charge", &Settings::charge)
      .def_readwrite("ignoreCharge", &Settings::ignoreCharge)
      .def_readwrite("spin", &Settings::spin)
      .def_readwrite("geometry", &Settings::geometry)
      .def_readwrite("load", &Settings::load)
      .def_readwrite("scfMode", &Settings::scfMode)
      .def_readwrite("method", &Settings::method)
      .def_readwrite("dft", &Settings::dft)
      .def_readwrite("scf", &Settings::scf)
      .def_readwrite("basis", &Settings::basis)
      .def_readwrite("grid", &Settings::grid)
      .def_readwrite("efield", &Settings::efield)
      .def_readwrite("pcm", &Settings::pcm)
      .def_readwrite("extCharges", &Settings::extCharges)
      .def_readwrite("customFunc", &Settings::customFunc);
}