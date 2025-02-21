#!/usr/bin/env python3
# @file   holeParticleCorrelation.py
#
# @date   Oct 12, 2021
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

import numpy as np
import matplotlib.pyplot as plt
import argparse

def plotHoleParticleCorrelation(filename):
    """
    Plot the hole-particle correlation from a given file.
    This script works with the files: <sysname> + ".correlation" + <iExc> + ".txt"
    which are generated from an LRSCFTask with the "transitionCharges"-keyword set to True.

    Parameters:
    filename (str): The path to the file containing the correlation data.

    Returns:
    None
    """
    correlationplot = np.loadtxt(filename)

    fig, ax1 = plt.subplots()

    im = ax1.imshow(correlationplot, "Greys", vmin=0, vmax=.8)

    ax1.tick_params(top=False, bottom=True,
                    labeltop=False, labelbottom=True)

    nAtoms = len(correlationplot)
    ax1.set_xticks(np.arange(nAtoms))
    ax1.set_yticks(np.arange(nAtoms))
    ax1.set_xticklabels(range(1, nAtoms), fontsize=9)
    ax1.set_yticklabels(range(nAtoms, 0, -1), fontsize=9)
    ax1.set_xlabel("Hole")
    ax1.set_ylabel("Particle")
    ax1.set_title("Electron correlation")

    ax1.figure.colorbar(im)

    fig.tight_layout()
    plt.savefig(filename + ".pdf", format="pdf", bbox_inches='tight')
    print("Plotted the hole-particle correlation plot to " + filename + ".pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the hole-particle correlation from a given file.")
    parser.add_argument("filename", type=str, help="The path to the file containing the correlation data.")
    args = parser.parse_args()

    plotHoleParticleCorrelation(args.filename)