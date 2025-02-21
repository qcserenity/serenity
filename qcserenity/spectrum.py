#!/usr/bin/env python3
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

import numpy as np
import argparse

def gaussBand(freq, ev, fl):
    """
    Calculate the value that a Gaussian-formed peak at ev has at position freq.

    Parameters:
    freq (float): The frequency at which to calculate the Gaussian band.
    ev (float): The excitation energy.
    fl (float): The oscillator strength or rotator strength.

    Returns:
    float: The Gaussian band value at the given frequency.
    """
    stddev = 0.15
    return fl * np.exp(-((freq - ev)/stddev)**2)
    
def parseSpectrum(outputfile):
    """
    Parse the spectrum data from an output file.

    Parameters:
    outputfile (str): The path to the output file containing the spectrum data.

    Returns:
    np.ndarray: A 2D numpy array containing the parsed spectrum data.
                The array has shape (n_eigen, 5), where n_eigen is the number of eigenvalues.
                Columns represent:
                - Column 0: Excitation energies
                - Column 1: Oscillator strengths (dipole-length)
                - Column 2: Oscillator strengths (dipole-velocity)
                - Column 3: Rotator strengths (dipole-length)
                - Column 4: Rotator strengths (dipole-velocity)
    """
    
    outfile = open(outputfile, "r")
    content = outfile.readlines()
    
    neigen = 0
    exspectrum = np.zeros((neigen, 5))
    
    i = 0
    for line in content:
      if ("Summary" in line):
        method = line.split()[0].upper()
      if ("Absorption Spectrum (dipole-length)" in line):
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
      if ("Absorption Spectrum (dipole-velocity)" in line):
        osc_v = [float(l.split()[4]) for l in content[i+5:i+5+n_eigen]]
        for iEigen in range(len(osc_v)):
          exspectrum[iEigen, 2] = osc_v[iEigen]
      if ("CD Spectrum (dipole-length)" in line):
        rot_l = [float(l.split()[4]) for l in content[i+5:i+5+n_eigen]]
        for iEigen in range(len(rot_l)):
          exspectrum[iEigen, 3] = rot_l[iEigen]
      if ("CD Spectrum (dipole-velocity)" in line):
        rot_v = [float(l.split()[4]) for l in content[i+5:i+5+n_eigen]]
        for iEigen in range(len(rot_v)):
          exspectrum[iEigen, 4] = rot_v[iEigen]
      i += 1
    
    return exspectrum

def plotSpectrum(exspectrum, savefile="spectrum.pdf", vlabel="dipole-velocity", abs=True, unit="eV", style="plt.rc", vlines=False, xlowlim=100, xupplim=500, **kwargs):
    """
    Plot the spectrum data.

    Parameters:
    exspectrum (np.ndarray): A 2D numpy array containing the spectrum data.
                             The array should have shape (n_eigen, 5), where n_eigen is the number of eigenvalues.
                             Columns represent:
                             - Column 0: Excitation energies
                             - Column 1: Oscillator strengths (dipole-length)
                             - Column 2: Oscillator strengths (dipole-velocity)
                             - Column 3: Rotator strengths (dipole-length)
                             - Column 4: Rotator strengths (dipole-velocity)
    savefile (str): The path to save the plot as a PDF file. Default is "spectrum.pdf".
    vlabel (str): The label for the velocity plot. Default is "dipole-velocity".
    abs (bool): If True, plot the absorption spectrum. If False, plot the CD spectrum. Default is True.
    unit (str): The unit for the x-axis. Can be "eV" or "nm". Default is "eV".
    style (str): The matplotlib style to use for the plot. Default is "plt.rc".
    vlines (bool): If True, plot vertical lines at the excitation energies. Default is False.
    xlowlim (float): The lower limit for the x-axis. Default is 100.
    xupplim (float): The upper limit for the x-axis. Default is 500.
    **kwargs: Additional keyword arguments to pass to the plt.plot() function.

    Returns:
    None
    """
    import matplotlib.pyplot as plt
    plt.style.use(style)
    freqStart = exspectrum[0,0] - 1
    freqEnd = exspectrum[-1,0] + 1
    nPoints = 1000
    cd = 0
    method = "CC2"
    if(abs):
      cd = 0
      method += " Absorption Spectrum"
    else:
      cd = 2
      method += " CD Spectrum"
    
    frequencies = np.linspace(freqStart, freqEnd, nPoints)
    absorptionsLength = np.linspace(0, 0, nPoints)
    absorptionsVelocity = np.linspace(0, 0, nPoints)
    
    for iPoint in range(nPoints):
      for jExc in range(0, np.shape(exspectrum)[0]):
        absorptionsLength[iPoint] += gaussBand(frequencies[iPoint], exspectrum[jExc, 0], exspectrum[jExc, 1 + cd])
        absorptionsVelocity[iPoint] += gaussBand(frequencies[iPoint], exspectrum[jExc, 0], exspectrum[jExc, 2 + cd])
    
    # Normalize.
    absorptionsLength /= max(absorptionsLength)
    absorptionsVelocity /= max(absorptionsVelocity)
    if unit == "nm":
        wavelengths = 1240 / frequencies  # Convert eV to nm
        plt.xlabel("Wavelength / nm")
        plt.plot(wavelengths[::-1], absorptionsVelocity[::-1], label=vlabel, **kwargs)  # Reverse for ascending order of wavelength
        plt.xlim(xlowlim, xupplim)
        if vlines:
          for i, wl in enumerate(wavelengths):
             if 100 <= wl <= 500:  # Plot bars only within visible range
                 plt.vlines(wl, 0, 0.8 * absorptionsVelocity[i], linestyle="dashed", linewidth=0.8)
    else:
        plt.xlabel("Frequency / eV")
        plt.plot(frequencies, absorptionsVelocity, label=vlabel, **kwargs)
        if xlowlim != 100 and xupplim != 500:
          plt.xlim(xlowlim, xupplim)
        if vlines:
          for i, freq in enumerate(exspectrum[:, 0]):
              plt.vlines(freq, 0, 1, linestyle="dashed", linewidth=0.8)

    
    plt.title(method)
    plt.ylabel("Absorption")
    #plt.plot(frequencies, absorptionsLength, label = "dipole-length")
    #plt.plot(frequencies, absorptionsVelocity, label=vlabel)
    plt.legend(frameon=0)
    plt.savefig(savefile, bbox_inches = "tight")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse and plot spectrum data.")
    parser.add_argument("outputfile", type=str, help="The path to the output file containing the spectrum data.")
    parser.add_argument("--savefile", type=str, default="spectrum.pdf", help="The path to save the plot as a PDF file. Default is 'spectrum.pdf'.")
    parser.add_argument("--vlabel", type=str, default="dipole-velocity", help="The label for the velocity plot. Default is 'dipole-velocity'.")
    parser.add_argument("--abs", type=bool, default=True, help="If True, plot the absorption spectrum. If False, plot the CD spectrum. Default is True.")
    parser.add_argument("--unit", type=str, choices=["eV", "nm"], default="eV", help="The unit for the x-axis. Can be 'eV' or 'nm'. Default is 'eV'.")
    parser.add_argument("--style", type=str, default="plt.rc", help="The matplotlib style to use for the plot. Default is 'plt.rc'.")
    parser.add_argument("--vlines", action="store_true", help="If True, plot vertical lines at the excitation energies. Default is False.")
    parser.add_argument("--xlowlim", type=float, default=100, help="The lower limit for the x-axis. Default is 100.")
    parser.add_argument("--xupplim", type=float, default=500, help="The upper limit for the x-axis. Default is 500.")
    args = parser.parse_args()

    exspectrum = parseSpectrum(args.outputfile)
    plotSpectrum(exspectrum, savefile=args.savefile, vlabel=args.vlabel, abs=args.abs, unit=args.unit, style=args.style, vlines=args.vlines, xlowlim=args.xlowlim, xupplim=args.xupplim)