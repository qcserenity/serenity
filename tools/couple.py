#@file   coupleTDA.py
#
#@date   Sep 28, 2021
#@author Niklas Niemeyer
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

#################################################################################################
###                                                                                           ###
### This script accepts Löwdin transition charges of an arbitrary                             ###
### number of subsystems and excitations computed beforehand and                              ###
### couples those excitations based on the transition charges.                                ###
### Analytical couplings from subsystem TDA can also be included.                             ###
###                                                                                           ###
### Usage:                                                                                    ###
### -------                                                                                   ###
### python couple.py <N_subsystems> <sys1> <sys2> .. <sysN> <N_analyticalCouplings> <12> <23> ###
###                                                                                           ###
### Here, N_subsystems are coupled named sys1, sys2, ... Further, N_analyticalCouplings       ###
### analytical couplings between the subsystems 1 and 2 as well as 2 and 3 are included.      ###
### These are read from disk just like the transition                                         ###
### charges. IMPORTANT: In all (FDEu) tasks set the transitionCharges keyword to true and     ###
### in all (FDEc) tasks for analytical couplings set the saveResponseMatrix keyword to true.  ###
### Otherwise, this script will fail because it won't be able to find everything it needs.    ###
###                                                                                           ###
#################################################################################################
import sys
import numpy as np

# Constants.
AU_TO_CGS = 64604.8164
HARTREE_TO_EV = 27.21138602
HARTREE_TO_NM = 45.56337117

# Simple distance function.
def distance(atom1, atom2):
  x1, y1, z1 = atom1[0], atom1[1], atom1[2]
  x2, y2, z2 = atom2[0], atom2[1], atom2[2]
  return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

nSubsystems = int(sys.argv[1])
names = sys.argv[2:2 + nSubsystems]
print("\n Number of systems      : ", nSubsystems)
print(" ----------------------------")
charges = []
spectra = []
nEigen = np.zeros((nSubsystems), dtype = int)
nAtoms = np.zeros((nSubsystems), dtype = int)

for I in range(nSubsystems):
  nameI = names[I]
  pathI = nameI + "/"
  charges.append(np.loadtxt(pathI + nameI + ".transitioncharges.txt", dtype = float))
  spectra.append(np.loadtxt(pathI + nameI + ".exspectrum.txt", dtype = float))
  if (len(charges[I].shape) == 1):
    # Only one atom.
    nEigen[I] = charges[I].shape[0] - 3
    nAtoms[I] = 1
    tmp = charges[I]
    charges[I] = np.zeros((1, nEigen[I] + 3), dtype = float)
    charges[I][0] = tmp
  else:
    # Several atoms.
    nEigen[I] = charges[I].shape[1] - 3
    nAtoms[I] = charges[I].shape[0]

  # Print some information about this subsystem.
  print("\n  System                : ", names[I])
  print("  Number of atoms       : ", nAtoms[I])
  print("  Number of excitations : ", nEigen[I])

# Prepare coupling matrix.
H = np.zeros((nEigen.sum(), nEigen.sum()))

# 1. Fill diagonal with uncoupled excitation energies.
for I in range(nSubsystems):
  for iEigen in range(nEigen[I]):
    H[nEigen[0:I].sum() + iEigen][nEigen[0:I].sum() + iEigen] = spectra[I][0][iEigen]

# 2. Insert requested blocks with loaded sTDA couplings.
couplings = []
if (len(sys.argv) > 2 + nSubsystems):
  nCouplings = int(sys.argv[2 + nSubsystems])
  print("\n Including %2i analytical coupling(s) from sTDA:"%(nCouplings))
  print(" ------------------------------------------------")
  for iCoupling in range(nCouplings):
    coupling = int(sys.argv[3 + nSubsystems + iCoupling])
    I = coupling // 10 - 1
    J = coupling % 10 - 1
    couplings.append((I, J))
    couplings.append((J, I))
    print("  %3i.   %-16s <---> %16s"%(iCoupling+1, names[I], names[J]))

    # Load block.
    nameI = names[I]
    pathI = nameI + "/"
    nameJ = names[J]
    pathJ = nameJ + "/"
    IJ = np.loadtxt(pathI + nameJ + ".TDACoupling.txt", dtype = float)
    JI = np.loadtxt(pathJ + nameI + ".TDACoupling.txt", dtype = float)
    if (np.amax(IJ - JI.transpose()) > 1e-6):
      print("Warning: analytical sTDA coupling blocks may not be symmetric.")
    
    # Insert into coupling matrix.
    for iEigen in range(nEigen[I]):
      for jEigen in range(nEigen[J]):
        H[nEigen[0:I].sum() + iEigen][nEigen[0:J].sum() + jEigen] = IJ[iEigen, jEigen]
        H[nEigen[0:J].sum() + jEigen][nEigen[0:I].sum() + iEigen] = JI[jEigen, iEigen]
else:
  print("\n No analytical couplings from sTDA will be included.")


# 3. Insert missing blocks with transition charge couplings.
for I in range(nSubsystems):
  for J in range(I + 1, nSubsystems):
    # Skip if already present.
    if ((I, J) in couplings):
      continue
    # Evaluate this coupling block.
    for iEigen in range(nEigen[I]):
      for jEigen in range(nEigen[J]):
        # Evaluate sum for this exciton pair.
        for iAtom in range(nAtoms[I]):
          for jAtom in range(nAtoms[J]):
            V = charges[I][iAtom][3+iEigen] * charges[J][jAtom][3+jEigen] / distance(charges[I][iAtom], charges[J][jAtom])
            H[nEigen[0:I].sum() + iEigen][nEigen[0:J].sum() + jEigen] += V
            H[nEigen[0:J].sum() + jEigen][nEigen[0:I].sum() + iEigen] += V

# Solve eigenvalue problem.
eigenvalues, eigenvectors = np.linalg.eigh(H)

# Calculate coupled transition moments.
len = np.zeros((3, nEigen.sum()))
vel = np.zeros((3, nEigen.sum()))
mag = np.zeros((3, nEigen.sum()))

iStart = 0
for I in range(nSubsystems):
  nEigenI = spectra[I].shape[1]
  len += spectra[I][1: 4] @ eigenvectors[iStart:iStart + nEigenI]
  vel += spectra[I][4: 7] @ eigenvectors[iStart:iStart + nEigenI]
  mag += spectra[I][7:10] @ eigenvectors[iStart:iStart + nEigenI]
  iStart += nEigenI

vel = vel @ np.diag(np.reciprocal(eigenvalues))

# Calculate transition strengths
S_ll, S_lv, S_vv, S_lm, S_vm = [], [], [], [], [], []

for iEigen in range(nEigen.sum()):
  # Here we do not take care of complex algebra so a -1 needs to be included 
  # where only one operator is imaginary.
  S_ll.append(np.outer(len[:, iEigen], len[:, iEigen]))
  S_lv.append(np.outer(-1.0 * len[:, iEigen], vel[:, iEigen]))
  S_vv.append(np.outer(vel[:, iEigen], vel[:, iEigen]))
  S_lm.append(np.outer(-1.0 * AU_TO_CGS * len[:, iEigen], mag[:, iEigen]))
  S_vm.append(np.outer(AU_TO_CGS * vel[:, iEigen], mag[:, iEigen]))

# Calculate n -> subsystem map.
nToSubsystem = np.zeros((nEigen.sum()), dtype = int)
nToExcitation = np.zeros((nEigen.sum()), dtype = int)
n = 0
for I in range(nSubsystems):
  for iEigen in range(nEigen[I]):
    nToSubsystem[n] = I
    nToExcitation[n] = iEigen
    n += 1

if (nEigen.sum() < 15):
  print("\n  Coupling Matrix / eV:")
  print(" -----------------------")
  for iEigen in range(nEigen.sum()):
    for jEigen in range(nEigen.sum()):
      print("%10.4f"%(H[iEigen, jEigen] * HARTREE_TO_EV), end="")
    print("")

print("---------------------------------------------------------------------------------------")
print("                               Dominant Contributions                                  ")
print("---------------------------------------------------------------------------------------")
print(" state       energy      wavelength       sys      excitation       contribution       ")
print("              (eV)          (nm)                                      100*|c|^2        ")
print("---------------------------------------------------------------------------------------")
for iEigen in range(nEigen.sum()):
  contributions = eigenvectors[:, iEigen] * eigenvectors[:, iEigen]
  first = True
  for c in range(nEigen.sum()):
    # Print only when contribution dominant.
    if (contributions[c] > 0.1):
      # Print contribution.
      if (first):
        print(" %3i %15.5f %12.1f %10i %14i %18.2f"%(
          iEigen + 1, eigenvalues[iEigen] * HARTREE_TO_EV, HARTREE_TO_NM / eigenvalues[iEigen],
          nToSubsystem[c] + 1, nToExcitation[c] + 1, 100 * contributions[c]))
        first = False
      else:
        print(" %32s %10i %14i %18.2f"%(
          "", nToSubsystem[c] + 1, nToExcitation[c] + 1, 100 * contributions[c]))
print("---------------------------------------------------------------------------------------")
print("                          Absorption Spectrum (dipole-length)                          ")
print("---------------------------------------------------------------------------------------")
print(" state       energy      wavelength        fosc          Sxx        Syy        Szz     ")
print("              (eV)          (nm)           (au)                    (au)                ")
print("---------------------------------------------------------------------------------------")
for iEigen in range(nEigen.sum()):
  print(" %3i %15.5f %12.1f %15.6f %12.5f %10.5f %10.5f"%(
    iEigen + 1, eigenvalues[iEigen] * HARTREE_TO_EV, HARTREE_TO_NM / eigenvalues[iEigen], 
    2.0/3.0 * eigenvalues[iEigen] * S_ll[iEigen].trace(), S_ll[iEigen][0, 0], S_ll[iEigen][1, 1], S_ll[iEigen][2, 2]))
print("---------------------------------------------------------------------------------------")
print("                         Absorption Spectrum (dipole-velocity)                         ")
print("---------------------------------------------------------------------------------------")
print(" state       energy      wavelength        fosc          Sxx        Syy        Szz     ")
print("              (eV)          (nm)           (au)                    (au)                ")
print("---------------------------------------------------------------------------------------")
for iEigen in range(nEigen.sum()):
  print(" %3i %15.5f %12.1f %15.6f %12.5f %10.5f %10.5f"%(
    iEigen + 1, eigenvalues[iEigen] * HARTREE_TO_EV, HARTREE_TO_NM / eigenvalues[iEigen], 
    2.0/3.0 * eigenvalues[iEigen] * S_vv[iEigen].trace(), S_vv[iEigen][0, 0], S_vv[iEigen][1, 1], S_vv[iEigen][2, 2]))
print("---------------------------------------------------------------------------------------")
print("                              CD Spectrum (dipole-length)                              ")
print("---------------------------------------------------------------------------------------")
print(" state       energy      wavelength         R            Sxx        Syy        Szz     ")
print("              (eV)          (nm)        (1e-40cgs)               (1e-40cgs)            ")
print("---------------------------------------------------------------------------------------")
for iEigen in range(nEigen.sum()):
  print(" %3i %15.5f %12.1f %15.4f %12.5f %10.5f %10.5f"%(
    iEigen + 1, eigenvalues[iEigen] * HARTREE_TO_EV, HARTREE_TO_NM / eigenvalues[iEigen], 
    S_lm[iEigen].trace(), S_lm[iEigen][0, 0], S_lm[iEigen][1, 1], S_lm[iEigen][2, 2]))
print("---------------------------------------------------------------------------------------")
print("                             CD Spectrum (dipole-velocity)                             ")
print("---------------------------------------------------------------------------------------")
print(" state       energy      wavelength         R            Sxx        Syy        Szz     ")
print("              (eV)          (nm)        (1e-40cgs)               (1e-40cgs)            ")
print("---------------------------------------------------------------------------------------")
for iEigen in range(nEigen.sum()):
  print(" %3i %15.5f %12.1f %15.4f %12.5f %10.5f %10.5f"%(
    iEigen + 1, eigenvalues[iEigen] * HARTREE_TO_EV, HARTREE_TO_NM / eigenvalues[iEigen], 
    S_vm[iEigen].trace(), S_vm[iEigen][0, 0], S_vm[iEigen][1, 1], S_vm[iEigen][2, 2]))
print("---------------------------------------------------------------------------------------")
