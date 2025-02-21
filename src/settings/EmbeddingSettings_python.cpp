/**
 * @file EmbeddingSettings_python.cpp
 *
 * @date 18 Aug 2019
 * @author Moritz Bensberg
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
#include "settings/EmbeddingSettings.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_EmbeddingSettings(py::module& spy) {
  py::class_<EmbeddingSettings>(spy, "EmbeddingSettings", "Embedding settings all set to their default values.")
      .def_readwrite("levelShiftParameter", &EmbeddingSettings::levelShiftParameter)
      .def_readwrite("naddXCFunc", &EmbeddingSettings::naddXCFunc)
      .def_readwrite("naddKinFunc", &EmbeddingSettings::naddKinFunc)
      .def_readwrite("longRangeNaddKinFunc", &EmbeddingSettings::longRangeNaddKinFunc)
      .def_readwrite("embeddingMode", &EmbeddingSettings::embeddingMode)
      .def_readwrite("dispersion", &EmbeddingSettings::dispersion)
      .def_readwrite("smoothFactor", &EmbeddingSettings::smoothFactor)
      .def_readwrite("singValThreshold", &EmbeddingSettings::singValThreshold)
      .def_readwrite("potentialBasis", &EmbeddingSettings::potentialBasis)
      .def_readwrite("singValThreshold", &EmbeddingSettings::singValThreshold)
      .def_readwrite("lbDamping", &EmbeddingSettings::lbDamping)
      .def_readwrite("lbCycles", &EmbeddingSettings::lbCycles)
      .def_readwrite("carterCycles", &EmbeddingSettings::carterCycles)
      .def_readwrite("borderAtomThreshold", &EmbeddingSettings::borderAtomThreshold)
      .def_readwrite("basisFunctionRatio", &EmbeddingSettings::basisFunctionRatio)
      .def_readwrite("truncateProjector", &EmbeddingSettings::truncateProjector)
      .def_readwrite("projecTruncThresh", &EmbeddingSettings::projecTruncThresh)
      .def_readwrite("fermiShift", &EmbeddingSettings::fermiShift)
      .def_readwrite("calculateMP2Correction", &EmbeddingSettings::calculateMP2Correction)
      .def_readwrite("customNaddXCFunc", &EmbeddingSettings::customNaddXCFunc)
      .def_readwrite("customNaddKinFunc", &EmbeddingSettings::customNaddKinFunc)
      .def_readwrite("customLongRangeNaddKinFunc", &EmbeddingSettings::customLongRangeNaddKinFunc)
      .def_readwrite("loewdinOrder", &EmbeddingSettings::loewdinOrder)
      .def_readwrite("loewdinWeights", &EmbeddingSettings::loewdinWeights)
      .def_readwrite("embeddingModeList", &EmbeddingSettings::embeddingModeList)
      .def_readwrite("naddXCFuncList", &EmbeddingSettings::naddXCFuncList)
      .def_readwrite("naddKinFuncList", &EmbeddingSettings::naddKinFuncList);
}
