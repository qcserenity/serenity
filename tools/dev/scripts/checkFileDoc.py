#!/usr/bin/python3
#@file   checkFileDoc.py
#
#@date   Jun 26, 2017
#@author Jan Unsleber
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

import sys,os

try:
  src = sys.argv[1]
  os.path.isfile(os.path.join(src,"/serenity.cpp"))
except:
  print("Path to src/ folder needed as first argument to this script.")
  sys.exit(1)

docerrors = []
codeerrors = []

whitelist = [
"./serenity.cpp",
"./python/serenipy.cpp",
"./io/Eigen3HDF5.h",
"./io/Eigen3HDF5_test.cpp",
"./grid/construction/sphere_lebedev_rule.cpp",
"./grid/construction/sphere_lebedev_rule.h"
]

cwd = os.getcwd()
os.chdir(src)
for root, dirnames, filenames in os.walk('.'):
  for filename in filenames:
    fname = os.path.join(root, filename)
    if fname in whitelist: continue
    if not (fname.split('.')[-1] in ['cpp','h']): continue
    firstlinedoxy = False
    atfile = False
    atauthor = False
    atdate = False
    atcopyright = False
    copyright1 = False
    copyright2 = False
    includeguards1 = fname.split('.')[-1]=='cpp'
    includeguards2 = fname.split('.')[-1]=='cpp'
    includeguards3 = fname.split('.')[-1]=='cpp'
    namespace = ('_python' in fname)
    guard = filename.split('.')[0]
    with open(fname,'r') as f:
      for l,line in enumerate(f):
        if (l==0) and ("/**" in line):
          firstlinedoxy = True
        if "@file" in line:
          if fname.split('/')[-1] == line.split()[-1]: atfile = True
        if "@author" in line: atauthor = True
        if "@date" in line: atdate = True
        if "@copyright" in line: atcopyright = True
        if "GNU Lesser General Public License" in line: copyright1 = True
        if "WITHOUT ANY WARRANTY" in line: copyright2 = True
        if ('#ifndef' in line) and (guard.upper() in line): includeguards1 = True
        if ('#define' in line) and (guard.upper() in line): includeguards2 = True
        if ('#endif' in line): includeguards3 = True
        if ('namespace' in line) and ('Serenity' in line) and not ('using' in line): namespace = True
    if not copyright2*copyright1*atcopyright*atdate*atauthor*atfile*firstlinedoxy:
      docerrors.append(fname)
    if not includeguards1*includeguards2*includeguards3*namespace:
      codeerrors.append(fname)

os.chdir(cwd)

fail = False
if len(docerrors) >0:
  print("The file documentation on the following files was missing, inclomplete or wrong:")
  for e in docerrors:
    print(os.path.join(src,e[2:]))
  fail= True
  print()

if len(codeerrors) >0:
  print("There was a namespace or include guard wrong or missing in the following files:")
  for e in codeerrors:
    print(os.path.join(src,e[2:]))
  fail= True

if fail: sys.exit(1);
