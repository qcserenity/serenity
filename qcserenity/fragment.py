#!/usr/bin/python3
# @file   fragment.py
#
# @date   Jun 14, 2021
# @author Moritz Bensberg
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
import numpy as np
from collections import Counter


# This class handles a xyz-file without needing ASE, which may not be installed by default.
class Geometry(object):
    """A geometry defined by nuclear symbols and coordinates"""

    def __init__(self, symbols, coordinates):
        super(Geometry, self).__init__()
        self.symbols = symbols
        self.coordinates = coordinates

    @classmethod
    def fromFile(cls, path):
        f = open(path, "r")
        lines = f.readlines()
        symbols = []
        nAtoms = int(lines[0].split()[0])
        coordinates = np.zeros((nAtoms, 3))
        for iAtom in range(nAtoms):
            line = lines[iAtom + 2]
            words = line.split()
            symbols.append(words[0].upper())
            coordinates[iAtom, 0] = float(words[1])
            coordinates[iAtom, 1] = float(words[2])
            coordinates[iAtom, 2] = float(words[3])
        return cls(symbols, coordinates)

    def toXYZFile(self, name, comment=""):
        nAtoms = len(self.symbols)
        xyzFileString = str(nAtoms) + "\n"
        xyzFileString += comment + "\n"
        for iAtom in range(nAtoms):
            xyzFileString += "{:<5}".format(self.symbols[iAtom])
            xyzFileString += "{:<15}".format(str(self.coordinates[iAtom, 0])) + " "
            xyzFileString += "{:<15}".format(str(self.coordinates[iAtom, 1])) + " "
            xyzFileString += "{:<15}".format(str(self.coordinates[iAtom, 2])) + "\n"
        with open(name, "w") as xyzFile:
            xyzFile.write(xyzFileString)

    def getNAtoms(self):
        return len(self.symbols)

    def getSumFormula(self):
        formula = ""
        count = Counter(self.symbols)
        for element in count:
            formula += element + str(count[element])
        return formula

    def print(self):
        for iAtom in range(len(self.symbols)):
            sym = self.symbols[iAtom]
            x = self.coordinates[iAtom, 0]
            y = self.coordinates[iAtom, 1]
            z = self.coordinates[iAtom, 2]
            print(f'  {sym:5}  {x:10f}  {y:10f}  {z:10f}')
        return

    def calculateDistanceMatrix(self):
        nAtoms = len(self.symbols)
        distances = np.zeros((nAtoms, nAtoms))
        for i in range(nAtoms):
            for j in range(nAtoms):
                if (j >= i): break;
                dist = np.linalg.norm(np.subtract(self.coordinates[i], self.coordinates[j]))
                distances[i, j] = dist
                distances[j, i] = dist
        return distances


# Construct the next atom cluster starting with the atom with index "seed"
def getNextCluster(seed, assigned, nAtoms, distances, radii, factor, sym):
    iFrag = [seed]
    assigned[seed] = True
    for iAtom in iFrag:
        iRad = radii[sym[iAtom]]
        for jAtom in range(nAtoms):
            if assigned[jAtom]:
                continue
            dist = distances[iAtom, jAtom]
            cutOff = factor * (iRad + radii[sym[jAtom]])
            if dist < cutOff:
                iFrag.append(jAtom)
                assigned[jAtom] = True
    return iFrag


# Get the next atom not already assigned.
def getNextSeed(assigned):
    for i in range(len(assigned)):
        if not assigned[i]:
            return i


# Check if NOT all atoms were assigned yet.
def checkIfNotAllAssigned(assigned):
    return False in assigned

def fragment(inputFile, factor=0.5):

    # The UFF-radii according to:  J. Am. Chem. Soc. 1992, 114, 25
    allRadii = [
        1.4430, 1.81, 1.2255, 1.3725, 2.0515, 1.9255, 1.83, 1.75, 1.682, 1.6215, 1.4915, 1.5105, 2.2495, 2.1475, 2.0735,
        2.0175, 1.9735, 1.934, 1.9060, 1.6995, 1.6475, 1.5875, 1.5720, 1.5115, 1.4805, 1.4560, 1.4360, 1.4170, 1.7475,
        1.3815, 2.1915, 2.14, 2.115, 2.1025, 2.0945, 2.0705, 2.0570, 1.8205, 1.6725, 1.5620, 1.5825, 1.526, 1.499, 1.4815,
        1.4645, 1.4495, 1.5740, 1.4240, 2.2315, 2.1960, 2.2100, 2.2350, 2.25, 2.2020, 2.2585, 1.8515, 1.7610, 1.7780,
        1.8030, 1.7875, 1.7735, 1.7600, 1.7465, 1.6840, 1.7255, 1.7140, 1.7045, 1.6955, 1.6870, 1.6775, 1.8200, 1.5705,
        1.850, 1.5345, 1.4770, 1.5600, 1.4200, 1.3770, 1.6465, 1.3525, 2.1735, 2.1485, 2.1850, 2.3545, 2.3750, 2.3825,
        2.4500, 1.8385, 1.7390, 1.6980, 1.7120, 1.6975, 1.7120, 1.7120, 1.6905, 1.6630, 1.6695, 1.6565, 1.6495, 1.6430,
        1.6370, 1.6240, 1.6180, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
    ]
    # All atom symbols.
    allSymbols = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
        "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
        "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
        "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
        "Ts", "Og"
    ]
    # Combine the atom symbols and the radii to a dictionary.
    radii = dict([(allSymbols[i].upper(), allRadii[i]) for i in range(len(allSymbols))])

    # Read the input file.
    totalGeometry = Geometry.fromFile(inputFile)
    nAtoms = totalGeometry.getNAtoms()
    coords = totalGeometry.coordinates
    sym = totalGeometry.symbols
    # Container for the cluster/molecule atom indices.
    # This will become a list of lists. Each list contains the
    # indices of the atom cluster.
    fragmentsIndices = []
    # Calculate all inter-atom distances.
    distances = totalGeometry.calculateDistanceMatrix()
    # Keep track on which atoms were already assigned.
    assigned = [False for i in range(nAtoms)]

    # Construct atom clusters.
    while checkIfNotAllAssigned(assigned):
        seed = getNextSeed(assigned)
        fragmentsIndices.append(getNextCluster(seed, assigned, nAtoms, distances, radii, factor, sym))

    # Construct new geometries and write them to new files.
    iFrag = 0
    for frag in fragmentsIndices:
        iFrag += 1
        nFragAtoms = len(frag)
        fragCoords = np.zeros((nFragAtoms, 3))
        fragSymbols = []
        counter = 0
        for iAtom in frag:
            fragSymbols.append(sym[iAtom])
            fragCoords[counter] = coords[iAtom]
            counter += 1
        newGeom = Geometry(fragSymbols, fragCoords)
        fileName = str(iFrag) + "_" + newGeom.getSumFormula() + ".xyz"
        newGeom.toXYZFile(fileName, "automatically fragmented")

if __name__ == "__main__":
    # Check for an input file.
    if len(sys.argv) < 2:
        raise Exception("No input file given. Exiting.")
    inputFile = sys.argv[1]
    try:
        open(inputFile, "r")
    except IOError:
        raise Exception("The input file does not exist or is corrupted!")
    factor = 0.5
    # Check if a user-specified radii-scaling is given (default is 0.5).
    if len(sys.argv) > 2:
        factor = float(sys.argv[2])
    fragment(inputFile, factor)
