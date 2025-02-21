# @file test_import.py
# 
# @date: Jun 20, 2018
# @author: Jan Unsleber
# @copyright \n
#  This file is part of the program Serenity.\n\n
#  Serenity is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.\n\n
#  Serenity is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.\n\n
#  You should have received a copy of the GNU Lesser General
#  Public License along with Serenity.
#  If not, see <http://www.gnu.org/licenses/>.\n

import unittest

class Import(unittest.TestCase):

    def test_import(self):
        import qcserenity.serenipy as spy
        self.assertIsNotNone(spy)
