#!/usr/bin/python3
#@file   updateFileList.py
#
#@date   May 29, 2020
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

cwd = os.getcwd()
os.chdir(src)
test_files = []
header_files = []
cpp_files = []
python_files = []
for root, dirnames, filenames in os.walk('.'):
  for filename in filenames:
    fname = os.path.join(*(root.split(os.path.sep)[1:]),filename)
    if '_test.cpp' == fname[-9:] or \
            '__TEST_SUPPLY.h' == fname[-15:] or \
            '__TEST_SUPPLY.cpp' == fname[-17:]:
        test_files.append(fname)
    elif '_python.cpp' == fname[-11:] or\
            'serenipy.cpp' == fname[-12:]:
        python_files.append(fname)
    elif '.cpp' == fname[-4:] and not 'serenity.cpp' == fname[-12:]:
        cpp_files.append(fname)
    elif '.h' == fname[-2:]:
        header_files.append(fname)

test_files.sort()
header_files.sort()
cpp_files.sort()
python_files.sort()

with open('Files.cmake', 'w') as f:
    f.write('set(SERENITY_CPPS\n')
    for i in cpp_files:
        f.write('  {}\n'.format(os.path.join('${PROJECT_SOURCE_DIR}','src',i)))
    f.write(')\n\n')
    f.write('set(SERENITY_HEADERS\n')
    for i in header_files:
        f.write('  {}\n'.format(os.path.join('${PROJECT_SOURCE_DIR}','src',i)))
    f.write(')\n\n')
    f.write('set(SERENITY_PYTHON_FILES\n')
    for i in python_files:
        f.write('  {}\n'.format(os.path.join('${PROJECT_SOURCE_DIR}','src',i)))
    f.write(')\n\n')
    f.write('set(SERENITY_TEST_FILES\n')
    for i in test_files:
        f.write('  {}\n'.format(os.path.join('${PROJECT_SOURCE_DIR}','src',i)))
    f.write(')\n\n')

os.chdir(cwd)
