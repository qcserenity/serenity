#!/usr/bin/env python3
# @file   __init__.py
#
# @date   Feb 19, 2025
# @author Anton Rikus
# @copyright \n
# This file is part of the program Serenity.\n\n
# Serenity is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.\n\n
# Serenity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.\n\n
# You should have received a copy of the GNU Lesser General
# Public License along with Serenity.
# If not, see <http://www.gnu.org/licenses/>.\n

import os
from pathlib import Path

# Automatically set SERENITY_RESOURCES
package_dir = Path(__file__).resolve().parent.parent
resources_dir = package_dir / "share" / "serenity" / "data"
if resources_dir.exists():
    os.environ["SERENITY_RESOURCES"] = str(resources_dir) + "/"

laplace_dir = package_dir / "share" / "laplace-minimax"
if laplace_dir.exists():
    os.environ["LAPLACE_ROOT"] = str(laplace_dir)

from .couple import couple_excitations
from .cubeUtils import addCubes, subtractCubes
from .fragment import fragment
from .redirectOutput import redirectOutputToFile
from .spectrum import parseSpectrum

__all__ = ["fragment", "addCubes", "subtractCubes", "redirectOutputToFile", "parseSpectrum", "couple_excitations"]

try:
    from .holeParticleCorrelation import plotHoleParticleCorrelation
    __all__.append("plotHoleParticleCorrelation")
except ImportError:
    pass

try:
    from .showCube import showCube
    __all__.append("showCube")
except ImportError:
    pass

try:
    from .spectrum import plotSpectrum
    __all__.append("plotSpectrum")
except ImportError:
    pass

try:
    from .twoDimheatMap import plot_2Dheatmap
    __all__.append("plot_2Dheatmap")
except ImportError:
    pass