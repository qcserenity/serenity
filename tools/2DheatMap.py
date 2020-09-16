#!/usr/bin/python3
#@file   2DheatMap.py
#
#@date   Feb 4, 2020
#@author Anja Massolle
#@copyright \n
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

import argparse
import sys
from distutils import util
from ase.io import read
from ase import Atom
from ase.data import covalent_radii as radii
from ase.data.colors import jmol_colors
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.interpolate import griddata

parser = argparse.ArgumentParser(description='Plot heat map from dat and corresponding xyz file.')
parser.add_argument('-dat', nargs=1, metavar='datFile', required=True,
                   help='dat file containing x, y and the value of the porperty (*_XYPLANE.dat).')
parser.add_argument('-xyz', metavar='xyzFile', nargs=1, required=True,
                   help='xyz file of the molecule in the plane (*_MOLECULE_ROTATED_TO_XYPLANE.xyz).')
parser.add_argument('--N', metavar='nPoints', type=int, default=500, nargs=1, required=False,
                   help='Number of points in x and y direction (only for the generation of the grid if the heatmap is interpolated).')
parser.add_argument('--title', metavar='title', type=str, default=None, required=False,
                   help='Title of the plot.')
parser.add_argument('--cut', metavar='cutoff', type=float, default=1.0, required=False,
                   help='Cut off radius for the displayed atoms in z direction.')
parser.add_argument('--min', metavar='vmin', type=float, default=None, required=False,
                   help='Minimum value for the property.')
parser.add_argument('--max', metavar='vmax', type=float, default=None, required=False,
                   help='Maximum value for the property.')
parser.add_argument('--latex', metavar='bool', type=lambda x:bool(util.strtobool(x)), default=False, \
                    required=False,help='Use the LaTeX font engine.')
parser.add_argument('--diverging', metavar='bool', type=lambda x:bool(util.strtobool(x)), default=False,\
                    required=False, help='Use a diverging colormap.')
parser.add_argument('--interpolate', metavar='bool', type=lambda x:bool(util.strtobool(x)), default=True,\
                    required=False,help='Interpolate the z values.')
args = parser.parse_args()

# Use LaTeX engine for text
if args.latex:
    from matplotlib import rcParams
    plt.rc('text', usetex=True)

# Load xyz file
molecule = read(args.xyz[0],format="xyz")
# Load data from .dat
dat = np.loadtxt(args.dat[0],skiprows=1)
X = dat[:,0].round(8)
Y = dat[:,1].round(8)
Z = dat[:,2]
if(len(dat[0])>3):
    sys.exit("The dat file has more than 3 columns, please use the _XYPLANE.dat file.")

# dimension of the picture
extent = [X.min(), X.max(), Y.min(), Y.max()]
fig, ax = plt.subplots()
if args.title != None:
    plt.title(args.title)
if args.latex:
    plt.ylabel('$y$ / \AA')
    plt.xlabel('$x$ / \AA')
else:
    plt.ylabel('$y$ / Å')
    plt.xlabel('$x$ / Å')
# Colors for the colorbar
colors = [
        (255/255,   0/255,   0/255, 1.0), # Red
        (255/255, 132/255,  17/255, 1.0), # Orange
        (255/255, 246/255,  35/255, 1.0), # Yellow
        ( 23/255, 230/255,   5/255, 1.0), # Green
        (107/255, 255/255, 237/255, 1.0), # Cyan
        (114/255, 135/255, 255/255, 1.0), # Blue
        (233/255, 134/255, 240/255, 1.0), # Violett
        (255/255, 255/255, 255/255, 1.0)  # White
        ]

cm_r = LinearSegmentedColormap.from_list(
        'cm_r', colors[::-1], N=200)

if (args.diverging):
    colormap = plt.cm.seismic
else:
    colormap = cm_r

if (args.interpolate):
    # create x-y points to be used in heatmap
    xi = np.linspace(X.min(),X.max(),args.N)
    yi = np.linspace(Y.min(),Y.max(),args.N)
    # Z is a matrix of x-y values
    z = griddata((X, Y), Z, (xi[None,:], yi[:,None]),method='cubic')
    plt.imshow(z, extent=extent,origin='lower', vmin=args.min, vmax=args.max, cmap=colormap)
else:
    unique, counts = np.unique(X, return_counts=True)
    lenX = counts[0]
    unique, counts = np.unique(Y, return_counts=True)
    lenY = counts[0]
    x=X.reshape(lenY,lenX)
    y=Y.reshape(lenY,lenX)
    z=Z.reshape(lenY,lenX)
    plt.pcolormesh(x,y,z,vmin=args.min, vmax=args.max, cmap=colormap)
plt.colorbar()
for atom in molecule:
    if abs(atom.z) <= args.cut:
        distancePercent = float(abs(atom.z) / args.cut)
        if atom.z < 0:
            color = jmol_colors[atom.number]*(1-distancePercent)
        else:
            color = jmol_colors[atom.number] + (1-jmol_colors[atom.number])* distancePercent
        # The atom number of dummy atoms is 0
        if atom.number == 0:
            color='k'
        radius = radii[atom.number]*0.4
        circle = Circle((atom.x, atom.y), radius, facecolor=color,\
                        edgecolor='k', linewidth=0.5, alpha=0.75*(1-distancePercent))
        ax.add_patch(circle)
        if atom.number != 0:
            label = ax.annotate(atom.symbol, xy=(atom.x, atom.y), ha="center", va="center")
plt.savefig("heatmap.pdf", dpi=600, format="pdf")
