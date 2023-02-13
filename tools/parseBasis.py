# @file   parseBasis.py
#
# @date   Mar 27, 2022
# @author Niklas Niemeyer
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

import re
import sys

###################################################################################################
###  Purpose of this script:                                                                    ###
###  This script's purpose is to "translate" basis set files downloaded from                    ###
###  basissetexchange.org in TURBOMOLE format so that Serenity can read them.                   ###
###                                                                                             ###
###  Requirements for the output basis file:                                                    ###
###   1. TURBOMOLE format.                                                                      ###  
###   2. Three spaces between element symbols and the basis name.                               ###
###   3. Name of basis file must be identical to how the basis is referenced within the file.   ###
###   4. If scientific notation is present, use "E" and never "D".                              ###
###                                                                                             ###
###   Usage: parseBasis.py downloaded_basis_file                                                ###
###################################################################################################

file_name = sys.argv[1]
data = open(file_name, "r").readlines()
basis_name = ""

# Find out basis name.
for i in range(len(data)):
    if (data[max(0,i-1)].strip() == "*" and data[min(i+1,len(data) - 1)].strip() == "*"):
      basis_name = data[i].split()[1]

print(" ---------------------------------------------------------------------------------------------------")
print("   Found basis: " + basis_name.upper())
basis_name = basis_name.upper()
if ("RIFIT" in  basis_name):
  basis_name = re.sub("RIFIT", "RI-C", basis_name)
  print(" Renaming to ..", basis_name)
output = open(basis_name, "w")

for i in range(len(data)):
    data[i] = re.sub("D-", "E-", data[i])
    # I don't know why but D+ has to be escaped explicitly ..
    data[i] = re.sub("D\+", "E+", data[i])
    if (data[max(0,i-1)].strip() == "*" and data[min(i+1,len(data) - 1)].strip() == "*"):
      data[i] = re.sub(data[i].split()[1], "  " + basis_name, data[i])
    output.write(data[i])

print(" ---------------------------------------------------------------------------------------------------")
print("   Move the file " + basis_name + " to the ~/serenity/data/basis folder")
print("   and reference it as " + basis_name + " in the input file.")
print(" ---------------------------------------------------------------------------------------------------")
