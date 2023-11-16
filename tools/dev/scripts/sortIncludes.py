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

whitelist = [
"./serenity.cpp",
"./memory/MemoryManager.cpp",
"./python/serenipy.cpp",
"./integrals/wrappers/Libint.h",
"./basis/Shell.h",
"./grid/construction/sphere_lebedev_rule.cpp",
"./grid/construction/sphere_lebedev_rule.h"
]

# cd to src/
cwd = os.getcwd()
os.chdir(src)

# get all possible internal headers
paths = {}
headers = []
for root, dirnames, filenames in os.walk('.'):
  for filename in filenames:
    if filename.split('.')[-1] == 'h':
      headers.append(filename)
      paths[filename] = os.path.join(root,filename)[2:]

for root, dirnames, filenames in os.walk('.'):
  for filename in filenames:
    fname = os.path.join(root, filename)
    if fname in whitelist: continue
    if not (fname.split('.')[-1] in ['cpp','h']): continue
   
    # read all includes
    includes = []
    firstincl = False
    lastincl = -1
    with open(fname,'r') as f:
      lines = f.readlines()
    for l,line in enumerate(lines):
      if "#include" in line: 
        includes.append(line.split()[1])
        if not firstincl:
          if '/* Include' in lines[l-1]:
            firstincl = l-1
          else:
            firstincl = l
        lastincl = l

    # clean and sort
    classheader = None
    internal = []
    external = []
    for incl in includes:
      if incl[1:-1].split('/')[-1].split('.')[0] == filename.split('.')[0]:
        classheader = incl[1:-1].split('/')[-1]
      elif incl[1:-1].split('/')[-1] in headers: 
        internal.append(incl[1:-1].split('/')[-1])
      else:
        external.append(incl[1:-1])
     
    internal.sort() 
    external.sort()

    # print to file
    with open(fname,'w') as f:
      for l,line in enumerate(lines):
        if l==firstincl:
          if classheader:
            f.write("/* Include Class Header*/\n")
            f.write('#include "{:s}"\n'.format(paths[classheader]))
          if len(internal) > 0:
            f.write("/* Include Serenity Internal Headers */\n")
            for i in internal:
              f.write('#include "{:s}"\n'.format(paths[i]))
          if len(external) > 0:
            f.write("/* Include Std and External Headers */\n")
            for e in external:
              f.write('#include <{:s}>\n'.format(e))
        if l>= firstincl and l<=lastincl: continue
        f.write(line)


   


