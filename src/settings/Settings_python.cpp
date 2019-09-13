/**
 * @file Settings_python.cpp
 *
 * @date Apr 25, 2016
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
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

using namespace Serenity;
using namespace Options;
namespace py = pybind11;


//Settings deepcopySettings(const Settings settings){return Settings(settings);};

void export_Settings(py::module &spy) {

  py::class_<DFT>(spy,"DFTSettings",
      "DFT settings all set to their default values.")
              .def_readwrite("functional",&DFT::functional)
              .def_readwrite("densityFitting",&DFT::densityFitting)
              .def_readwrite("dispersion",&DFT::dispersion);

  py::class_<SCF>(spy,"SCFSettings",
      "SCF settings all set to their default values.")
      .def_readwrite("initialguess",&SCF::initialguess)
      .def_readwrite("maxCycles",&SCF::maxCycles)
      .def_readwrite("writeRestart",&SCF::writeRestart)
      .def_readwrite("energyThreshold",&SCF::energyThreshold)
      .def_readwrite("rmsdThreshold",&SCF::rmsdThreshold)
      .def_readwrite("damping",&SCF::damping)
      .def_readwrite("seriesDampingStart",&SCF::seriesDampingStart)
      .def_readwrite("seriesDampingEnd",&SCF::seriesDampingEnd)
      .def_readwrite("seriesDampingStep",&SCF::seriesDampingStep)
      .def_readwrite("seriesDampingInitialSteps",&SCF::seriesDampingInitialSteps)
      .def_readwrite("staticDampingFactor",&SCF::staticDampingFactor)
      .def_readwrite("endDampErr",&SCF::endDampErr)
      .def_readwrite("useLevelshift",&SCF::useLevelshift)
      .def_readwrite("useOffDiagLevelshift",&SCF::useOffDiagLevelshift)
      .def_readwrite("diisFlush",&SCF::diisFlush)
      .def_readwrite("diisStartError",&SCF::diisStartError)
      .def_readwrite("diisMaxStore",&SCF::diisMaxStore)
      .def_readwrite("diisThreshold",&SCF::diisThreshold)
      .def_readwrite("diisConditionNumberThreshol",&SCF::diisConditionNumberThreshold)
      .def_readwrite("useADIIS",&SCF::useADIIS);

  py::class_<BASIS>(spy,"BasisSettings",
      "Basis settings all set to their default values.")
      .def_readwrite("label",&BASIS::label)
      .def_readwrite("auxCLabel",&BASIS::auxCLabel)
      .def_readwrite("auxJLabel",&BASIS::auxJLabel)
      .def_readwrite("makeSphericalBasis",&BASIS::makeSphericalBasis)
      .def_readwrite("integralThreshold",&BASIS::integralThreshold)
      .def_readwrite("basisLibPath",&BASIS::basisLibPath);

  py::class_<GRID>(spy,"GridSettings",
      "Grid settings all set to their default values.")
        .def_readwrite("gridType",&GRID::gridType)
        .def_readwrite("radialGridType",&GRID::radialGridType)
        .def_readwrite("sphericalGridType",&GRID::sphericalGridType)
        .def_readwrite("blocksize",&GRID::blocksize)
        .def_readwrite("accuracy",&GRID::accuracy)
        .def_readwrite("smallGridAccuracy",&GRID::smallGridAccuracy)
        .def_readwrite("blockAveThreshold",&GRID::blockAveThreshold)
        .def_readwrite("basFuncRadialThreshold",&GRID::basFuncRadialThreshold)
        .def_readwrite("weightThreshold",&GRID::weightThreshold)
        .def_readwrite("smoothing",&GRID::smoothing);

  py::class_<Settings>(spy, "Settings")
        .def(py::init<>())
        .def_readwrite("name",&Settings::name)
        .def_readwrite("identifier",&Settings::identifier)
        .def_readwrite("path",&Settings::path)
        .def_readwrite("charge",&Settings::charge)
        .def_readwrite("spin",&Settings::spin)
        .def_readwrite("geometry",&Settings::geometry)
        .def_readwrite("load",&Settings::load)
        .def_readwrite("scfMode",&Settings::scfMode)
        .def_readwrite("method",&Settings::method)
        .def_readwrite("dft",&Settings::dft)
        .def_readwrite("scf",&Settings::scf)
        .def_readwrite("basis",&Settings::basis)
        .def_readwrite("grid",&Settings::grid);

}


