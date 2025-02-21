#!/usr/bin/python3
# @file   cubeUtils.py
#
# @date   Oct 11, 2021
# @author Patrick Eschenbach, Albrecht Goetz
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

### GLOBAL PARAMETERS ###

import numpy as np
import argparse

# Bohr to Angstrom conversion
bohr_radius = 0.52917721092

# Element numbers and symbols
elements = {
   'H'  :  1 ,
   'He' :  2 ,
   'Li' :  3 ,
   'Be' :  4 ,
   'B'  :  5 ,
   'C'  :  6 ,
   'N'  :  7 ,
   'O'  :  8 ,
   'F'  :  9 ,
   'Ne' : 10 ,
   'Na' : 11 ,
   'Mg' : 12 ,
   'Al' : 13 ,
   'Si' : 14 ,
   'P'  : 15 ,
   'S'  : 16 ,
   'Cl' : 17 ,
   'Ar' : 18 ,
   'K'  : 19 ,
   'Ca' : 20 ,
   'Sc' : 21 ,
   'Ti' : 22 ,
   'V'  : 23 ,
   'Cr' : 24 ,
   'Mn' : 25 ,
   'Fe' : 26 ,
   'Co' : 27 ,
   'Ni' : 28 ,
   'Cu' : 29 ,
   'Zn' : 30 ,
   'Ga' : 31 ,
   'Ge' : 32 ,
   'As' : 33 ,
   'Se' : 34 ,
   'Br' : 35 ,
   'Kr' : 36 ,
   'Rb' : 37 ,
   'Sr' : 38 ,
   'Y'  : 39 ,
   'Zr' : 40 ,
   'Nb' : 41 ,
   'Mo' : 42 ,
   'Tc' : 43 ,
   'Ru' : 44 ,
   'Rh' : 45 ,
   'Pd' : 46 ,
   'Ag' : 47 ,
   'Cd' : 48 ,
   'In' : 49 ,
   'Sn' : 50 ,
   'Sb' : 51 ,
   'Te' : 52 ,
   'I'  : 53 ,
   'Xe' : 54 ,
   'Cs' : 55 ,
   'Ba' : 56 ,
   'La' : 57 ,
   'Ce' : 58 ,
   'Pr' : 59 ,
   'Nd' : 60 ,
   'Pm' : 61 ,
   'Sm' : 62 ,
   'Eu' : 63 ,
   'Gd' : 64 ,
   'Tb' : 65 ,
   'Dy' : 66 ,
   'Ho' : 67 ,
   'Er' : 68 ,
   'Tm' : 69 ,
   'Yb' : 70 ,
   'Lu' : 71 ,
   'Hf' : 72 ,
   'Ta' : 73 ,
   'W'  : 74 ,
   'Re' : 75 ,
   'Os' : 76 ,
   'Ir' : 77 ,
   'Pt' : 78 ,
   'Au' : 79 ,
   'Hg' : 80 ,
   'Tl' : 81 ,
   'Pb' : 82 ,
   'Bi' : 83 ,
   'Po' : 84 ,
   'At' : 85 ,
   'Rn' : 86 ,
   'Fr' : 87 ,
   'Ra' : 88 ,
   'Ac' : 89 ,
   'Th' : 90 ,
   'Pa' : 91 ,
   'U'  : 92 ,
   'Np' : 93 ,
   'Pu' : 94 ,
   'Am' : 95 ,
   'Cm' : 96 ,
   'Bk' : 97 ,
   'Cf' : 98 ,
   'Es' : 99 ,
   'Fm' : 100 ,
   'Md' : 101 ,
   'No' : 102 ,
   'Lr' : 103 ,
   'Rf' : 104 ,
   'Db' : 105 ,
   'Sg' : 106 ,
   'Bh' : 107 ,
   'Hs' : 108 ,
   'Mt' : 109 ,
   'Ds' : 110 ,
   'Rg' : 111 ,
   'Cn' : 112 
}

### CLASSES ###

class molecule:
   """ Class to represent a molecule. """
   __name__ = 'molecule'

   def __init__(self, coords, symbols):
      """ Constructor for molecule. """

      if len(coords) != len(symbols):
         raise RuntimeError('Number of atomic coordinates and symbols not equal.')

      else:
         self.coords = coords
         self.symbols = symbols


class grid:
   """ Class to represent a grid. """
   __name__ = 'grid'

   def __init__(self, coords, weights):
      """ Constructor for grid. """

      self.coords = coords
      self.weights = weights

   def get_coords(self):

      return self.coords

   def get_weights(self):

      return self.weights

   def get_unique(self):

      xyz = [[a[i] for a in self.get_coords()] for i in range(3)]
      xyz_unique = [sorted(list(set(xyz[i]))) for i in range(3)]

      return xyz_unique

   def sort_regular(self, style='cub'):
      """ Return new grid sorted according to Cube or PLT style.
          Works only for full cube grids. """

      newcoords = []
      uxyz = self.get_unique()

      # Check if this is a regular grid
      if len(uxyz[0]) * len(uxyz[1]) * len(uxyz[2]) != len(self.coords):
         raise RuntimeError('Only regular cubic grids can be sorted with sort_regular.')

      if style == 'cub':
         for x in uxyz[0] :
            for y in uxyz[1] :
               for z in uxyz[2] :
                  newcoords.append([x, y, z])
 
      elif style == 'plt':
         for z in uxyz[2] :
            for y in uxyz[1] :
               for x in uxyz[0] :
                  newcoords.append([x, y, z])

      else:
         raise RuntimeError('Unknown order type: '+style)

      # Regular cube grid can only have identical weights
      return grid(np.array(newcoords), self.weights)

   def sort_irregular(self, style='cub'):
      """ Returns new grid sorted according to Cube or PLT style.
          Works for all grids, but is rather inefficient. """

      if style == 'cub':
         ind = np.lexsort((self.coords[:,2], self.coords[:,1], self.coords[:,0]))

      elif style == 'plt':
         ind = np.lexsort((self.coords[:,0], self.coords[:,1], self.coords[:,2]))

      else:
         raise RuntimeError('Unknown order type: '+style)

      return grid(self.coords[ind], self.weights[ind])


class property:
   """ Class to represent a molecular property on a grid. """
   __name__ = 'property'

   def __init__(self, grid, values, type=None):

      self.grid = grid
      self.values = values
      self.type = type

   def __add__(self, other):

      return property(self.grid, np.add(self.values, other.values))

   def __sub__(self, other):

      return property(self.grid, np.subtract(self.values, other.values))

   def get_grid(self):

      return self.grid

   def get_values(self):

      return self.values

   def integrate(self, abs=False):

      if abs:
         return np.sum(np.absolute(self.values))
      else:
         return np.sum(self.values)

   def sort_irregular(self, style='cub'):

      if style == 'cub':
         ind = np.lexsort((self.grid.coords[:,2], self.grid.coords[:,1], self.grid.coords[:,0]))

      elif style == 'plt':
         ind = np.lexsort((self.grid.coords[:,0], self.grid.coords[:,1], self.grid.coords[:,2]))

      else:
         raise RuntimeError('Unknown order type: '+style)

      return property(grid(self.grid.coords[ind], self.grid.weights[ind]), self.values[ind], type)


class cubefile:
   """ Class to represent an existing Gaussian Cube file. """

   def __init__(self, fname):

      # General attributes
      self.fname = fname

      f = open(fname, 'r')
      self.comment1 = f.readline()
      self.comment2 = f.readline()

      # Number of atoms
      line = f.readline().split()
      self.natom = int(line[0])

      # Grid properties
      self.origin = [float(x) for x in line[1:4]]
      self.vox = []
      self.incr = []
      self.grid = []
      for i in range(3):
         line = f.readline().split()
         self.vox.append(int(line[0]))
         # Make sure we have a cubic grid
         for j in range(3):
            if i != j and float(line[j+1]) != 0:
               raise RuntimeError('Currently, only cubic grids are supported.')
         self.incr.append(float(line[i+1]))

      self.nvox = self.vox[0]*self.vox[1]*self.vox[2]

      # Molecular geometry
      self.coords = np.empty((0,3), float)
      self.symbols = np.empty((0,0), int)
      for i in range(self.natom):
         line = f.readline().split()
         self.symbols = np.append(self.symbols, int(line[0]))
         self.coords = np.vstack((self.coords, np.array([float(line[2]), float(line[3]), float(line[4])])))
      # Read property values
      self.values = []
      for line in f:
         self.values.extend(map(float, line.split()))
      self.values = np.array(self.values)
      f.close()

   def get_molecule(self):
      """ Return the molecule read from the cube file. """

      mol = molecule(self.coords, self.symbols)

      return mol

   def get_grid(self):
      """ Construct grid object from cube file. """

      coords = []
      for x in range(self.vox[0]) :
         for y in range(self.vox[1]) :
            for z in range(self.vox[2]) :
               coords.append([self.origin[0] + x*self.incr[0], self.origin[1] + y*self.incr[1], self.origin[2] + z*self.incr[2]])

      coords = np.array(coords)
      weights = np.ones(len(coords))
      self.grid = grid(coords, weights)
      return grid(coords, weights)

   def set_grid(self, newgrid):
      self.grid = newgrid
 
   def get_values(self):
      """ Extract property values from cube file. """

      return self.values

   def get_property(self):
      """ Extract property object from cube file. """

      return property(self.get_grid(), self.get_values())

class xyzfile:
   """ Class to represent an existing XYZ file. """
   __name__ = 'xyzfile'

   def __init__(self, fname):

      f = open(fname, 'r')

      self.natoms = int(f.readline())
      self.comment = f.readline()

      self.symbols = []
      self.coords = []
      self.atoms = []

      for i in range(self.natoms):
         line = f.readline().split()
         try:
            self.symbols.append(elements[line[0]])
            self.atoms.append(line[0])
         except:
            self.symbols.append(0)

         self.coords.append([x for x in map(float, line[1:4])])

      self.symbols = np.array(self.symbols)
      self.coords = np.array(self.coords)

      f.close()


   def get_molecule(self):

       mol = molecule(self.coords, self.symbols)

       return mol
      
   def get_coordinates(self):
      return self.coords

   def set_coordinates(self,newcoords):
      self.coords=newcoords



class xyzvfile:
   """ Class to represent an existing XYZV file. """
   __name__ = 'xyzvfile'

   def __init__(self, fname, style='pyadf'):

      # General attributes
      self.fname = fname

      f = open(fname, 'r')

      # Try to guess number of comment lines
      self.comment = ''
      if style == 'pyadf':
         nclines = 2

      elif style == 'tm':
         nclines = 8

      else:
         nclines = 0

      for i in range(nclines):
         self.comment += f.readline()

      # Read grid coordinates and property values
      self.gcoords = []
      self.values = []

      for line in f:
         l = line.split()
         # Hopefully ignore any additional non-relevant lines
         if len(l) == 4 and l[0][0] != '#':
            self.gcoords.append(map(float, l[0:3]))
            self.values.append(float(l[3]))

      self.gcoords = np.array(self.gcoords)
      self.values = np.array(self.values)

      f.close()

   def get_grid(self):
      """ Construct grid object from XYZV file. """

      weights = np.ones(len(self.gcoords))

      return grid(self.gcoords, weights)

   def get_values(self):
      """ Extract property values from XYZV file. """

      return self.values
 
   def get_property(self):
      """ Extract property object from XYZV file. """

      return property(self.get_grid(), self.get_values())

### FUNCTIONS ###

def write_cubefile(fname, molecule, prop, origin=None, vox=None, incr=None, comment=None):

   f = open(fname, 'w')

   grid = prop.get_grid()

   # Check if grid and values match
   nvox = len(prop.get_grid().get_coords())
   if nvox != len(prop.get_values()):
      raise RuntimeError('No match between number of grid points and values.')

   # Try to determine cube grid automatically
   if vox == None or incr == None:

      xyz_unique = grid.get_unique()
      vox = [len(xyz_unique[i]) for i in range(3)]
      origin = [xyz_unique[0][0], xyz_unique[1][0], xyz_unique[2][0]]

      incr =[]
      for dim in xyz_unique:
         spacing = dim[1] - dim[0]
         for el in range(1,len(dim)):
            if abs(dim[el] - dim[el-1] - spacing) > 1e-7:
               print(str(dim[el] - dim[el-1])+' is not equal to increment '+str(spacing)+'.')
               raise RuntimeError('Given grid not equally spaced.')
         else:
            incr.append(spacing)

    
   if comment == None:
      comment = 'Cube file generated by cubeUtils\n'+'Number of points: '+str(nvox)+'\n'

   if comment.count('\n') != 2:
      raise RuntimeError('Comment string must contain exactly two newlines.')

   cubestr = comment
   cubestr += "%5d%12.6f%12.6f%12.6f\n" % \
              (len(molecule.coords), origin[0], origin[1], origin[2])
   cubestr += "%5d%12.6f%12.6f%12.6f\n" % \
              (vox[0], incr[0], 0.0, 0.0)
   cubestr += "%5d%12.6f%12.6f%12.6f\n" % \
              (vox[1], 0.0, incr[1], 0.0)
   cubestr += "%5d%12.6f%12.6f%12.6f\n" % \
              (vox[2], 0.0, 0.0, incr[2])

   for symbol, coord in zip(molecule.symbols, molecule.coords):
      cubestr += "%5d%12.6f%12.6f%12.6f%12.6f\n" % \
                 (symbol, 0.0, coord[0], coord[1], coord[2])
  
   for i in range(nvox):
      if (i+1) % 6 == 0:
         cubestr += "%12.8e\n" % prop.values[i]
      else:
         cubestr += "%12.8e  " % prop.values[i]

   f.write(cubestr)  
   f.close()

def blowup_cube(smallprop, largergrid):
   """ Transfer properties on small grid to large one,
       filling up with zeros where necessary. This is
       useful for converting data on e.g. an isosurface
       to a full cube representation.
       Returns a property object for the full grid. """

   lv = np.zeros(len(largergrid.coords))

   lindex = 0
   
   sp_sorted = smallprop.sort_irregular()
   lg_sorted = largergrid.sort_irregular()

   for i,coord1 in enumerate(sp_sorted.grid.coords):
      for j,coord2 in enumerate(lg_sorted.coords[lindex:], lindex):
         for c in range(3):
            print(i,j,c)
            if abs(coord1[c] - coord2[c]) > 1e-7:
               print("break")
               break
         else:
            lv[j] = smallprop.values[i]
            lindex = j
            break

      else:
         raise RuntimeError('Small grid is not a subset of large grid.')

   return property(lg_sorted, lv)

def subtractCubes(fname1,fname2,outfile):
   f1 = cubefile(fname1)
   f2 = cubefile(fname2)
   diff = f1.get_property()-f2.get_property()
   write_cubefile(outfile, f1.get_molecule(), diff, origin=f1.origin, vox=f1.vox, incr=f1.incr)

def addCubes(fname1,fname2,outfile):
   f1 = cubefile(fname1)
   f2 = cubefile(fname2)
   diff = f1.get_property()+f2.get_property()
   write_cubefile(outfile, f1.get_molecule(), diff, origin=f1.origin, vox=f1.vox, incr=f1.incr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cube file utilities")
    subparsers = parser.add_subparsers(dest="command")

    # Subparser for addCubes
    parser_add = subparsers.add_parser("add", help="Add two cube files")
    parser_add.add_argument("fname1", type=str, help="First cube file")
    parser_add.add_argument("fname2", type=str, help="Second cube file")
    parser_add.add_argument("outfile", type=str, nargs="?", default="sum.cube", help="Output cube file")

    # Subparser for subtractCubes
    parser_subtract = subparsers.add_parser("subtract", help="Subtract two cube files")
    parser_subtract.add_argument("fname1", type=str, help="First cube file")
    parser_subtract.add_argument("fname2", type=str, help="Second cube file")
    parser_subtract.add_argument("outfile", type=str, nargs="?", default="diff.cube", help="Output cube file")

    args = parser.parse_args()

    if args.command == "add":
        addCubes(args.fname1, args.fname2, args.outfile)
    elif args.command == "subtract":
        subtractCubes(args.fname1, args.fname2, args.outfile)
    else:
        parser.print_help()