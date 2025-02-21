# @file   redirectOutput.py
#
# @date   Feb 19, 2025
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

import os
import sys

class redirectOutputToFile(object):
    '''
    A context manager that redirects stdout to a file for its scope, usage:

    with redirectOutputToFile('output.txt'):
        os.system('ls -l')
    '''
    
    def __init__(self, filepath, *args, **kw):
        sys.stdout.flush()
        self._origstdout = sys.stdout
        self._oldstdout_fno = os.dup(sys.stdout.fileno())
        self._file = open(filepath, 'w')  # Open your file for writing

    def __enter__(self):
        self._newstdout = os.dup(1)
        os.dup2(self._file.fileno(), 1)  # Redirect stdout to the file
        sys.stdout = os.fdopen(self._newstdout, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._origstdout
        sys.stdout.flush()
        os.dup2(self._oldstdout_fno, 1)  # Restore original stdout
        self._file.close()  # Ensure the file is closed