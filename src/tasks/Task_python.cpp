/**
 * @file Task_python.cpp
 *
 * @author Moritz Bensberg
 * @date Oct 23, 2019
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
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

void export_GeneralTaskSettings(py::module& spy) {
  py::class_<GeneralTaskSettings>(spy, "GeneralTaskSettings", "GeneralTaskSettings settings all set to their default values.")
      .def_readwrite("printLevel", &GeneralTaskSettings::printLevel);
}