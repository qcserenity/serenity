/**
 * @file serenipy.cpp
 *
 * @date: Apr 25, 2016
 * @author: Jan Unsleber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#include <omp.h>
#include <pybind11/iostream.h> // Include to allow ostream redirect
#include <pybind11/pybind11.h>
#include <iostream>
#include <memory>
namespace py = pybind11;

/*
 * Forward declare everything
 */
void export_BasicFunctionals(py::module& spy);
void export_BasisController(py::module& spy);
void export_CompositeFunctionals(py::module& spy);
void export_CoupledClusterTask(py::module& spy);
void export_DensityMatrixController(py::module& spy);
void export_DispersionCorrectionTask(py::module& spy);
void export_ElectronicStructure(py::module& spy);
void export_ElectronTransferTask(py::module& spy);
void export_EmbeddingSettings(py::module& spy);
void export_GeneralTaskSettings(py::module& spy);
void export_EnergyContributions(py::module& spy);
void export_ExportCavityTask(py::module& spy);
void export_FDETask(py::module& spy);
void export_FreezeAndThawTask(py::module& spy);
void export_FXDTask(py::module& spy);
void export_GWTask(py::module& spy);
void export_Geometry(py::module& spy);
void export_GeometryOptimizationTask(py::module& spy);
void export_GeneralizedDOSTask(py::module& spy);
void export_GradientTask(py::module& spy);
void export_GridController(py::module& spy);
void export_ImportCavityTask(py::module& spy);
void export_Libint(py::module& spy);
void export_Looper(py::module& spy);
void export_LocalizationTask(py::module& spy);
void export_LRSCFTask(py::module& spy);
void export_MP2Task(py::module& spy);
void export_MultipoleMomentTask(py::module& spy);
void export_Options(py::module& spy);
void export_PlotTask(py::module& spy);
void export_ProjectionBasedEmbTask(py::module& spy);
void export_OrbitalsIOTask(py::module& spy);
void export_ScfTask(py::module& spy);
void export_Settings(py::module& spy);
void export_SystemController(py::module& spy);
void export_SystemAdditionTask(py::module& spy);
void export_SystemSplittingTask(py::module& spy);
void export_Timings(py::module& spy);
void export_VirtualOrbitalSpaceSelectionTask(py::module& spy);

PYBIND11_MODULE(serenipy, spy) {
  setvbuf(stdout, NULL, _IONBF, BUFSIZ);

  // Initialize the stream redirector
  // static StreamRedirector streamRedirector;

  spy.doc() = "Serenipy - Python Bindings for Serenity ";

  export_BasicFunctionals(spy);
  export_BasisController(spy);
  export_CompositeFunctionals(spy);
  export_DensityMatrixController(spy);
  export_DispersionCorrectionTask(spy);
  export_ElectronicStructure(spy);
  export_EmbeddingSettings(spy);
  export_EnergyContributions(spy);
  export_GeneralTaskSettings(spy);
  export_Geometry(spy);
  export_GridController(spy);
  export_Libint(spy);
  export_Options(spy);
  export_Settings(spy);
  export_SystemController(spy);
  export_Timings(spy);

  export_CoupledClusterTask(spy);
  export_ElectronTransferTask(spy);
  export_ExportCavityTask(spy);
  export_FDETask(spy);
  export_FreezeAndThawTask(spy);
  export_FXDTask(spy);
  export_GWTask(spy);
  export_GeometryOptimizationTask(spy);
  export_GeneralizedDOSTask(spy);
  export_ImportCavityTask(spy);
  export_Looper(spy);
  export_LocalizationTask(spy);
  export_LRSCFTask(spy);
  export_MP2Task(spy);
  export_GradientTask(spy);
  export_MultipoleMomentTask(spy);
  export_PlotTask(spy);
  export_ProjectionBasedEmbTask(spy);
  export_ScfTask(spy);
  export_SystemAdditionTask(spy);
  export_SystemSplittingTask(spy);
  export_OrbitalsIOTask(spy);
  export_VirtualOrbitalSpaceSelectionTask(spy);
}
