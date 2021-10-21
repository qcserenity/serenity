#!/usr/bin/python3
# @file   addCubes.py
#
# @date   Oct 11, 2021
# @author Patrick Eschenbach
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

import os, sys
serenityDir = ""
if (os.getenv("SERENITY_HOME") != None):
    serenityDir = os.getenv("SERENITY_HOME")
else:
    print("SERENITY_HOME variable not set, you need to source the serenity.sh in the Serenity home directory!")
    sys.exit()
toolsDir = os.path.join(serenityDir,"tools")
sys.path.insert(0,toolsDir)
import CubeUtils as cu

if (len(sys.argv) > 3):
    outfile = sys.argv[3]
else:
    outfile = "sum.cube"
cu.addCubes(sys.argv[1],sys.argv[2], outfile)