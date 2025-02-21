# @file test_tasks.py
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
import qcserenity.serenipy as spy
import shutil
import os

class TestTasks(unittest.TestCase):

    def test_SCFTask_HF_restricted(self):
        sett = spy.Settings()
        sett.geometry = os.path.join(os.environ['SERENITY_RESOURCES'],'xyzfiles','water.xyz')
        sett.name = "tmp"
        sett.method = spy.ELECTRONIC_STRUCTURE_THEORIES.HF
        sett.basis.label = "STO-3G"
        sys = spy.System(sett)
        task = spy.ScfTask_R(sys)
        task.run()
        self.assertAlmostEqual(sys.getEnergy(), -74.96056774259112, places=7, msg=None, delta=None)
        shutil.rmtree('tmp')

    def test_SCFTask_HF_unrestricted(self):
        sett = spy.Settings()
        sett.geometry = os.path.join(os.environ['SERENITY_RESOURCES'],'xyzfiles','water.xyz')
        sett.name = "tmp"
        sett.method = spy.ELECTRONIC_STRUCTURE_THEORIES.HF
        sett.basis.label = "STO-3G"
        sys = spy.System(sett)
        task = spy.ScfTask_U(sys)
        task.run()
        self.assertAlmostEqual(sys.getEnergy(), -74.96056774259112, places=7, msg=None, delta=None)
        shutil.rmtree('tmp')

    def test_SCFTask_DFT_restricted(self):
        sett = spy.Settings()
        sett.geometry = os.path.join(os.environ['SERENITY_RESOURCES'],'xyzfiles','water.xyz')
        sett.name = "tmp"
        sett.method = spy.ELECTRONIC_STRUCTURE_THEORIES.DFT
        sett.basis.label = "STO-3G"
        sett.dft.functional = spy.PBE
        sett.customFunc.basicFunctionals = [spy.BASIC_FUNCTIONALS.C_PBE, spy.BASIC_FUNCTIONALS.X_PBE]
        sett.customFunc.mixingFactors = [1.0, 1.0]
        sett.customFunc.impl = spy.XC_IMPLEMENTATIONS.LIBXC
        sys = spy.System(sett)
        task = spy.ScfTask_R(sys)
        task.run()
        self.assertAlmostEqual(sys.getEnergy(), -75.2218044178, places=7, msg=None, delta=None)
        shutil.rmtree('tmp')

    def test_SCFTask_DFT_unrestricted(self):
        sett = spy.Settings()
        sett.geometry = os.path.join(os.environ['SERENITY_RESOURCES'],'xyzfiles','water.xyz')
        sett.name = "tmp"
        sett.method = spy.ELECTRONIC_STRUCTURE_THEORIES.DFT
        sett.basis.label = "STO-3G"
        sett.dft.functional = spy.PBE
        sys = spy.System(sett)
        task = spy.ScfTask_U(sys)
        task.run()
        self.assertAlmostEqual(sys.getEnergy(), -75.22180876683022, places=7, msg=None, delta=None)
        shutil.rmtree('tmp')
    
    def test_LRSCFTask(self):
        sett = spy.Settings()
        sett.geometry = os.path.join(os.environ['SERENITY_RESOURCES'],'xyzfiles','water.xyz')
        sett.name = "tmp"
        sett.method = spy.ELECTRONIC_STRUCTURE_THEORIES.HF
        sett.basis.label = "def2-svp"
        sys = spy.System(sett)
        task = spy.LRSCFTask_R([sys], [])
        task.settings.nEigen = 1
        task.settings.ltconv = 1e-7
        task.settings.method = spy.LR_METHOD.CC2
        task.run()
        self.assertAlmostEqual(task.getTransitions()[0][0], 0.3175885, places=6, msg=None, delta=None)
        shutil.rmtree('tmp')
