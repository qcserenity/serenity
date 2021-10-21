#!/usr/bin/python3
# @file   spectrum.py
#
# @date   Jun 14, 2021
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


import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if (len(sys.argv) == 1):
  exit("Specify a Serenity output file.")
if (len(sys.argv) == 4):
  exit("Only two command-line arguments allowed: file name and 'abs' for absorption spectra and 'cd' for CD spectra.")

outfile = open(sys.argv[1], "r")
content = outfile.readlines()

neigen = 0
exspectrum = np.zeros((neigen, 5))

i = 0
for line in content:
  if ("summary" in line):
    method = line.split()[0].upper()
  if ("absorption spectrum (dipole-length)" in line):
    j = 5
    finding = True
    while(finding):
      if (content[i + j][0:5] == "-----"):
        finding = False
        break
      j += 1
    n_eigen = j - 5
    exspectrum.resize((n_eigen, 5))

    exenergies = [float(l.split()[1]) for l in content[i+5:i+5+n_eigen]]
    osc_l = [float(l.split()[4]) for l in content[i+5:i+5+n_eigen]]
    for iEigen in range(len(osc_l)):
      exspectrum[iEigen, 0] = exenergies[iEigen]
      exspectrum[iEigen, 1] = osc_l[iEigen]
  if ("absorption spectrum (dipole-velocity)" in line):
    osc_v = [float(l.split()[4]) for l in content[i+5:i+5+n_eigen]]
    for iEigen in range(len(osc_v)):
      exspectrum[iEigen, 2] = osc_v[iEigen]
  if ("cd spectrum (dipole-length)" in line):
    rot_l = [float(l.split()[4]) for l in content[i+5:i+5+n_eigen]]
    for iEigen in range(len(rot_l)):
      exspectrum[iEigen, 3] = rot_l[iEigen]
  if ("cd spectrum (dipole-velocity)" in line):
    rot_v = [float(l.split()[4]) for l in content[i+5:i+5+n_eigen]]
    for iEigen in range(len(rot_v)):
      exspectrum[iEigen, 4] = rot_v[iEigen]
  i += 1

freqStart = exspectrum[0,0] - 1
freqEnd = exspectrum[-1,0] + 1
nPoints = 1000

cd = 0
if(len(sys.argv) == 2 or sys.argv[2] == "abs"):
  cd = 0
  method += " Absorption Spectrum"
elif(sys.argv[2] == "cd"):
  cd = 2
  method += " CD Spectrum"
else:
  exit("Your second command-line argument should be 'abs' or 'cd'.")

def gaussBand(freq, ev, fl):
    stddev = 0.15
    return fl * np.exp(-((freq - ev)/stddev)**2)

frequencies = np.linspace(freqStart, freqEnd, nPoints)
absorptionsLength = np.linspace(0, 0, nPoints)
absorptionsVelocity = np.linspace(0, 0, nPoints)

for iPoint in range(nPoints):
  for jExc in range(0, n_eigen):
    absorptionsLength[iPoint] += gaussBand(frequencies[iPoint], exspectrum[jExc, 0], exspectrum[jExc, 1 + cd])
    absorptionsVelocity[iPoint] += gaussBand(frequencies[iPoint], exspectrum[jExc, 0], exspectrum[jExc, 2 + cd])

# Normalize.
absorptionsLength /= max(absorptionsLength)
absorptionsVelocity /= max(absorptionsVelocity)

plt.title(method)
plt.xlabel("Frequency / eV")
plt.ylabel("Absorption")
plt.plot(frequencies, absorptionsLength, label = "dipole-length")
plt.plot(frequencies, absorptionsVelocity, label = "dipole-velocity")
plt.legend(frameon=0)
plt.savefig("spectrum.pdf", bbox_inches = "tight");
