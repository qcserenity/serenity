#!/usr/bin/env python3
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
from ase.data import covalent_radii as radii
from ase.data.colors import jmol_colors
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.interpolate import griddata

def plot_2Dheatmap(dat_file, xyz_file, N=500, title=None, cut=1.0, vmin=None, vmax=None, latex=False, diverging=False, interpolate=True):
    """
    Plot a 2D heatmap from a .dat file and corresponding .xyz file.

    Parameters:
    dat_file (str): Path to the .dat file containing x, y, and the value of the property.
    xyz_file (str): Path to the .xyz file of the molecule in the plane.
    N (int): Number of points in x and y direction for grid generation if interpolated. Default is 500.
    title (str): Title of the plot. Default is None.
    cut (float): Cut off radius for the displayed atoms in z direction. Default is 1.0.
    vmin (float): Minimum value for the property. Default is None.
    vmax (float): Maximum value for the property. Default is None.
    latex (bool): Use the LaTeX font engine. Default is False.
    diverging (bool): Use a diverging colormap. Default is False.
    interpolate (bool): Interpolate the z values. Default is True.

    Returns:
    None
    """
    # Use LaTeX engine for text
    if latex:
        from matplotlib import rcParams
        plt.rc('text', usetex=True)

    # Load xyz file
    molecule = read(xyz_file, format="xyz")
    # Load data from .dat
    dat = np.loadtxt(dat_file, skiprows=1)
    X = dat[:, 0].round(8)
    Y = dat[:, 1].round(8)
    Z = dat[:, 2]
    if len(dat[0]) > 3:
        sys.exit("The dat file has more than 3 columns, please use the _XYPLANE.dat file.")

    # Dimension of the picture
    extent = [X.min(), X.max(), Y.min(), Y.max()]
    fig, ax = plt.subplots()
    if title:
        plt.title(title)
    if latex:
        plt.ylabel(r'$y$ / \AA')
        plt.xlabel(r'$x$ / \AA')
    else:
        plt.ylabel('$y$ / Å')
        plt.xlabel('$x$ / Å')

    # Colors for the colorbar
    colors = [
        (255/255, 0/255, 0/255, 1.0),  # Red
        (255/255, 132/255, 17/255, 1.0),  # Orange
        (255/255, 246/255, 35/255, 1.0),  # Yellow
        (23/255, 230/255, 5/255, 1.0),  # Green
        (107/255, 255/255, 237/255, 1.0),  # Cyan
        (114/255, 135/255, 255/255, 1.0),  # Blue
        (233/255, 134/255, 240/255, 1.0),  # Violet
        (255/255, 255/255, 255/255, 1.0)  # White
    ]

    cm_r = LinearSegmentedColormap.from_list('cm_r', colors[::-1], N=200)

    colormap = plt.cm.seismic if diverging else cm_r

    if interpolate:
        # Create x-y points to be used in heatmap
        xi = np.linspace(X.min(), X.max(), N)
        yi = np.linspace(Y.min(), Y.max(), N)
        # Z is a matrix of x-y values
        z = griddata((X, Y), Z, (xi[None, :], yi[:, None]), method='cubic')
        plt.imshow(z, extent=extent, origin='lower', vmin=vmin, vmax=vmax, cmap=colormap)
    else:
        unique, counts = np.unique(X, return_counts=True)
        lenX = counts[0]
        unique, counts = np.unique(Y, return_counts=True)
        lenY = counts[0]
        x = X.reshape(lenY, lenX)
        y = Y.reshape(lenY, lenX)
        z = Z.reshape(lenY, lenX)
        plt.pcolormesh(x, y, z, vmin=vmin, vmax=vmax, cmap=colormap)

    plt.colorbar()
    for atom in molecule:
        if abs(atom.z) <= cut:
            distancePercent = float(abs(atom.z) / cut)
            if atom.z < 0:
                color = jmol_colors[atom.number] * (1 - distancePercent)
            else:
                color = jmol_colors[atom.number] + (1 - jmol_colors[atom.number]) * distancePercent
            # The atom number of dummy atoms is 0
            if atom.number == 0:
                color = 'k'
            radius = radii[atom.number] * 0.4
            circle = Circle((atom.x, atom.y), radius, facecolor=color,
                            edgecolor='k', linewidth=0.5, alpha=0.75 * (1 - distancePercent))
            ax.add_patch(circle)
            if atom.number != 0:
                ax.annotate(atom.symbol, xy=(atom.x, atom.y), ha="center", va="center")

    plt.savefig("heatmap.pdf", dpi=600, format="pdf")
    print("Plotted the heatmap to heatmap.pdf")

def main():
    parser = argparse.ArgumentParser(description='Plot heat map from dat and corresponding xyz file.')
    parser.add_argument('-dat', required=True, help='dat file containing x, y and the value of the property (*_XYPLANE.dat).')
    parser.add_argument('-xyz', required=True, help='xyz file of the molecule in the plane (*_MOLECULE_ROTATED_TO_XYPLANE.xyz).')
    parser.add_argument('--N', type=int, default=500, help='Number of points in x and y direction (only for the generation of the grid if the heatmap is interpolated).')
    parser.add_argument('--title', type=str, default=None, help='Title of the plot.')
    parser.add_argument('--cut', type=float, default=1.0, help='Cut off radius for the displayed atoms in z direction.')
    parser.add_argument('--min', type=float, default=None, help='Minimum value for the property.')
    parser.add_argument('--max', type=float, default=None, help='Maximum value for the property.')
    parser.add_argument('--latex', type=lambda x: bool(util.strtobool(x)), default=False, help='Use the LaTeX font engine.')
    parser.add_argument('--diverging', type=lambda x: bool(util.strtobool(x)), default=False, help='Use a diverging colormap.')
    parser.add_argument('--interpolate', type=lambda x: bool(util.strtobool(x)), default=True, help='Interpolate the z values.')
    args = parser.parse_args()

    plot_2Dheatmap(args.dat, args.xyz, args.N, args.title, args.cut, args.min, args.max, args.latex, args.diverging, args.interpolate)

if __name__ == "__main__":
    main()