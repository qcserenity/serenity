#@file   couple.py
#
#@date   Sep 29, 2021
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

##############################################################################################
##                                                                                          ##
## This script accepts Löwdin transition charges of an arbitrary number of subsystems and   ##
## excitations computed beforehand and couples those excitations based on the transition    ##
## charges. Analytical couplings (FDEc) and simple dipole--dipole couplings based on the    ##
## transition dipole moments can also be included.                                          ##
##                                                                                          ##
## Usage:                                                                                   ##
## -------                                                                                  ##
## python couple.py <N_subsystems> <sys1> <sys2> .. <sysN>                                  ##
##                  <N_FDEcCouplings> <sys1#sys2> <sys2#sys3> ..                            ##
##                  <N_TCCouplings> <sys1#sys3> <sys2#sys4> ..                              ##
##                  <N_DDCouplings> <sys1#sys3> <sys2#sys4> ..                              ##
##                                                                                          ##
## This script must be executed in the same folder all system folders are located in.       ##
##                                                                                          ##
## IMPORTANT:                                                                               ##
## In all (FDEu) tasks set the transitionCharges keyword to true.                           ##
## in all (FDEc) tasks for FDEc couplings set partialresponseConstruction keyword to true.  ##
## Otherwise, this script will fail because it won't be able to find everything it needs.   ##
##                                                                                          ##
##############################################################################################

import sys
import numpy as np
import time
import argparse
import multiprocessing
from multiprocessing import Pool


# Constants.
AU_TO_CGS = 64604.8164
HARTREE_TO_EV = 27.21138602
HARTREE_TO_NM = 45.56337117

# Simple distance function.
def distance(atom1, atom2):
  x1, y1, z1 = atom1[0], atom1[1], atom1[2]
  x2, y2, z2 = atom2[0], atom2[1], atom2[2]
  return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def couple_excitations(args):
  start_time = time.time()
  all_time = time.time()
  nSubsystems = args.n_subsystems
  nFDECouplings = 0
  nTCCouplings = 0
  names = args.systems
  name_to_index = {}
  print("\n Number of systems      : ", nSubsystems)
  print(" ----------------------------")
  charges = []
  spectra = []
  nEigen = np.zeros((nSubsystems), dtype = int)
  nAtoms = np.zeros((nSubsystems), dtype = int)

  for I in range(nSubsystems):
    nameI = names[I]
    name_to_index[nameI] = I
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


  # Print information.
  print("\n Number of roots to be determined : ", nEigen.sum())

  n_threads = multiprocessing.cpu_count()
  print(" Number of threads used : ", n_threads)

  # Prepare coupling matrix.
  H = np.zeros((nEigen.sum(), nEigen.sum()))

  # 1. Fill diagonal with uncoupled excitation energies.
  for I in range(nSubsystems):
    for iEigen in range(nEigen[I]):
      H[nEigen[0:I].sum() + iEigen][nEigen[0:I].sum() + iEigen] = spectra[I][0][iEigen]

  ###########################
  ### START FDEc COUPLING ###
  ###########################

  # 2. Insert requested blocks with loaded sTDA couplings.
  couplings = []
  if args.n_fdecouplings > 0:
    nFDECouplings = args.n_fdecouplings
    print("\n Including %5i TDA couplings from Serenity :"%(nFDECouplings))
    print(" ------------------------------------------------")
    for iCoupling in range(nFDECouplings):
      coupling = args.fdecouplings[iCoupling].split("#")
      I = name_to_index[coupling[0]]
      J = name_to_index[coupling[1]]
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
    print("\n No TDA couplings from FDEc will be included.")

  print("\n Time for loading (min):  %7.3f\n"%((time.time() - start_time) / 60))

  #########################
  ### START TC COUPLING ###
  #########################

  start_time = time.time()

  # To parallelize, must define function to return coupling for one exciton pair.
  def getTCCoupling(ij):
    # Coupling variable.
    V = 0
    # Identify subsystem indices and transition indices
    i = ij // nEigen.sum()
    j = ij % nEigen.sum()
    I, J = 0, 0
    iEigen, jEigen = 0, 0
    n = 0
    for sys in range(nSubsystems):
      for exc in range(nEigen[sys]):
        if (n == i):
          I = sys
          iEigen = exc
        if (n == j):
          J = sys
          jEigen = exc
        n += 1

    # Evaluate sum for this exciton pair.
    V = 0.0
    for iAtom in range(nAtoms[I]):
      for jAtom in range(nAtoms[J]):
        V += charges[I][iAtom][3+iEigen] * charges[J][jAtom][3+jEigen] / distance(charges[I][iAtom], charges[J][jAtom])
    
    return V

  # Determine composite index of coupling to be calculated.
  TC_coupling_incides = []
  if args.n_tccouplings > 0:
    nTCCouplings = args.n_tccouplings
    print("\n Including %5i transition-charge couplings :"%(nTCCouplings))
    print(" ------------------------------------------------")
    for iCoupling in range(nTCCouplings):
      coupling = args.tccouplings[iCoupling].split("#")
      I = name_to_index[coupling[0]]
      J = name_to_index[coupling[1]]
      
      print("  %3i.   %-16s <---> %16s"%(iCoupling+1, names[I], names[J]))
      # Skip if already present.
      if ((I, J) in couplings):
        continue
      # Evaluate this coupling block.
      for iEigen in range(nEigen[I]):
        for jEigen in range(nEigen[J]):
          i = nEigen[0:I].sum() + iEigen
          j = nEigen[0:J].sum() + jEigen
          TC_coupling_incides.append(i * nEigen.sum() + j)
  else:
    print(" No transition-charge couplings will be included.")

  print("\n Number of TC couplings to be calculated : ", len(TC_coupling_incides))

  # Calculate couplings.
  TC_coupling_values = []
  with Pool(processes=n_threads) as pool:
      TC_coupling_values = pool.map(getTCCoupling, TC_coupling_incides)

  # Distribute couplings.
  for iCoupl in range(len(TC_coupling_incides)):
    V = TC_coupling_values[iCoupl]
    ij = TC_coupling_incides[iCoupl]
    i = ij // nEigen.sum()
    j = ij % nEigen.sum()
    
    H[i][j] = V
    H[j][i] = V

  print(" Time for TC couplings (min):  %7.3f\n"%((time.time() - start_time) / 60))

  #########################
  ### START DD COUPLING ###
  #########################
  start_time = time.time()

  # To parallelize, must define function to return coupling for one exciton pair.
  def getDDCoupling(ij):
    # Coupling variable.
    V = 0
    # Identify subsystem indices and transition indices
    i = ij // nEigen.sum()
    j = ij % nEigen.sum()
    I, J = 0, 0
    iEigen, jEigen = 0, 0
    n = 0
    for sys in range(nSubsystems):
      for exc in range(nEigen[sys]):
        if (n == i):
          I = sys
          iEigen = exc
        if (n == j):
          J = sys
          jEigen = exc
        n += 1
    
    # Get centers of both monomers.
    R_I = np.array([0.0,0.0,0.0])
    for iAtom in range(nAtoms[I]):
      R_I += np.array(charges[I][iAtom][0:3])

    R_J = np.array([0.0,0.0,0.0])
    for jAtom in range(nAtoms[J]):
      R_J += np.array(charges[J][jAtom][0:3])

    R_I /= nAtoms[I]
    R_J /= nAtoms[J]

    # Get distance between monomers.
    R = R_J - R_I

    mu_i = np.array(spectra[I][1:4,iEigen])
    mu_j = np.array(spectra[J][1:4,jEigen])

    # Evaluate dipole-dipole coupling for this exciton pair.
    V = np.dot(mu_i,mu_j) / np.linalg.norm(R)**3 - 3 * np.dot(mu_i,R) * np.dot(mu_j,R) / np.linalg.norm(R)**5
    
    return V

  # Determine composite index of coupling to be calculated.
  DD_coupling_incides = []
  if args.n_ddcouplings > 0:
    nDDCouplings = args.n_ddcouplings
    print("\n Including %5i dipole-dipole couplings :"%(nDDCouplings))
    print(" ------------------------------------------------")
    for iCoupling in range(nDDCouplings):
      coupling = args.ddcouplings[iCoupling].split("#")
      I = name_to_index[coupling[0]]
      J = name_to_index[coupling[1]]
      
      print("  %3i.   %-16s <---> %16s"%(iCoupling+1, names[I], names[J]))
      # Skip if already present.
      if ((I, J) in couplings):
        continue
      # Evaluate this coupling block.
      for iEigen in range(nEigen[I]):
        for jEigen in range(nEigen[J]):
          i = nEigen[0:I].sum() + iEigen
          j = nEigen[0:J].sum() + jEigen
          DD_coupling_incides.append(i * nEigen.sum() + j)
  else:
    print(" No dipole-dipole couplings will be included.")

  print("\n Number of dipole-dipole couplings to be calculated : ", len(DD_coupling_incides))

  # Calculate couplings.
  DD_coupling_values = []
  with Pool(processes=n_threads) as pool:
      DD_coupling_values = pool.map(getDDCoupling, DD_coupling_incides)

  # Distribute couplings.
  for iCoupl in range(len(DD_coupling_incides)):
    V = DD_coupling_values[iCoupl]
    ij = DD_coupling_incides[iCoupl]
    i = ij // nEigen.sum()
    j = ij % nEigen.sum()
    
    H[i][j] = V
    H[j][i] = V


  print(" Time for DD couplings (min):  %7.3f\n"%((time.time() - start_time) / 60))
  start_time = time.time()


  # Solve eigenvalue problem.
  print(" ********************************************* ")
  print(" ****  Hamilton matrix assembly done  ******** ")
  print(" ********************************************* ")

  print("\n Diagonalizing Hamilton matrix ...")
  eigenvalues, eigenvectors = np.linalg.eigh(H)

  print(" Time H Diagonalization (min):  %7.3f\n"%((time.time() - start_time) / 60))
  start_time = time.time()

  print(" Calculating transition moments ...")
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
  S_ll, S_lv, S_vv, S_lm, S_lm_mod, S_vm = [], [], [], [], [], []

  for iEigen in range(nEigen.sum()):
    # Here we do not take care of complex algebra so a -1 needs to be included 
    # where only one operator is imaginary.
    S_ll.append(np.outer(len[:, iEigen], len[:, iEigen]))
    S_lv.append(np.outer(-1.0 * len[:, iEigen], vel[:, iEigen]))
    S_vv.append(np.outer(vel[:, iEigen], vel[:, iEigen]))
    S_lm.append(np.outer(-1.0 * AU_TO_CGS * len[:, iEigen], mag[:, iEigen]))
    S_vm.append(np.outer(AU_TO_CGS * vel[:, iEigen], mag[:, iEigen]))

    U, S, Vt = np.linalg.svd(S_lv[iEigen])
    S_lm_mod.append(U.transpose() @ S_lm[iEigen] @ Vt.transpose())

  # Calculate n -> subsystem map.
  nToSubsystem = np.zeros((nEigen.sum()), dtype = int)
  nToExcitation = np.zeros((nEigen.sum()), dtype = int)
  n = 0
  for I in range(nSubsystems):
    for iEigen in range(nEigen[I]):
      nToSubsystem[n] = I
      nToExcitation[n] = iEigen
      n += 1

  print(" Time for transition moments (min):  %7.3f\n"%((time.time() - start_time) / 60))

  if (nEigen.sum() < 1024):
    print("  Coupling Matrix / eV:")
    print(" -----------------------")
    for iEigen in range(nEigen.sum()):
      for jEigen in range(nEigen.sum()):
        print("%12.3e"%(H[iEigen, jEigen] * HARTREE_TO_EV), end="")
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
      if (contributions[c] > 0.01):
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
  print("                            CD Spectrum (mod. dipole-length)                           ")
  print("---------------------------------------------------------------------------------------")
  print(" state       energy      wavelength         R            Sxx        Syy        Szz     ")
  print("              (eV)          (nm)        (1e-40cgs)               (1e-40cgs)            ")
  print("---------------------------------------------------------------------------------------")
  for iEigen in range(nEigen.sum()):
    print(" %3i %15.5f %12.1f %15.4f %12.5f %10.5f %10.5f"%(
      iEigen + 1, eigenvalues[iEigen] * HARTREE_TO_EV, HARTREE_TO_NM / eigenvalues[iEigen], 
      S_lm_mod[iEigen].trace(), S_lm_mod[iEigen][0, 0], S_lm_mod[iEigen][1, 1], S_lm_mod[iEigen][2, 2]))
  print("---------------------------------------------------------------------------------------")

  print(" Total elapsed time for coupling (min):  %7.3f\n"%((time.time() - all_time) / 60))
  print(" All done. Have a nice day!")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Couple excitations based on transition charges.")
    parser.add_argument("N_subsystems", type=int, help="Number of subsystems")
    parser.add_argument("subsystems", nargs='+', help="List of subsystem names")
    parser.add_argument("N_FDEcCouplings", type=int, help="Number of FDEc couplings")
    parser.add_argument("FDEcCouplings", nargs='+', help="List of FDEc couplings in the format sys1#sys2")
    parser.add_argument("N_TCCouplings", type=int, help="Number of TC couplings")
    parser.add_argument("TCCouplings", nargs='+', help="List of TC couplings in the format sys1#sys3")
    parser.add_argument("N_DDCouplings", type=int, help="Number of DD couplings")
    parser.add_argument("DDCouplings", nargs='+', help="List of DD couplings in the format sys1#sys3")
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Print parsed arguments for debugging
    print(f"N_subsystems: {args.N_subsystems}")
    print(f"N_FDEcCouplings: {args.N_FDEcCouplings}")
    print(f"N_TCCouplings: {args.N_TCCouplings}")
    print(f"N_DDCouplings: {args.N_DDCouplings}")
    
    couple_excitations(args)

if __name__ == "__main__":
    main()