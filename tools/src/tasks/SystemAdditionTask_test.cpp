/**
 * @file SystemAdditionTask_test.cpp
 *
 * @date 26 Jan 2020
 * @author Moritz Bensberg
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "tasks/SystemAdditionTask.h"                 //To be tested.
#include "basis/BasisController.h"                    //GetNBasisFunctions.
#include "data/ElectronicStructure.h"                 //GetDensityMatrix.
#include "data/OrbitalController.h"                   //GetCoefficients.
#include "geometry/Geometry.h"                        //Get number of atoms.
#include "integrals/OneElectronIntegralController.h"  //Overlap matrix.
#include "settings/Settings.h"                        //Settings
#include "system/SystemController.h"                  //GetElectronicStructure.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class SystemAdditionTaskTest : public ::testing::Test {
 protected:
  SystemAdditionTaskTest() {
  }
  virtual ~SystemAdditionTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(SystemAdditionTaskTest, mixedBasis) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  Settings settingsSys2 = sys2->getSettings();
  settingsSys2.basis.label = "DEF2-SVP";
  settingsSys2.name = "TMP_SUPERSYSTEM";
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settingsSys2);
  SystemAdditionTask<RESTRICTED> add(supersystem, {sys1, sys2});
  add.run();

  DensityMatrix<RESTRICTED> superDens = supersystem->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  MatrixInBasis<RESTRICTED> superOverlap = supersystem->getOneElectronIntegralController()->getOverlapIntegrals();
  // The orbitals of the water molecules will be projected to the new basis. Check the
  // number of electrons that "survive" the projection.
  double nElec = (superDens.array() * superOverlap.array()).sum();
  EXPECT_NEAR(20, nElec, 0.01);
  EXPECT_EQ(20, supersystem->getNElectrons<RESTRICTED>());

  std::string supersystemName = supersystem->getSystemName();
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".settings").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".xyz").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".energies.res").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".orbs.res.h5").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".dmat.res.h5").c_str());
  std::remove((supersystem->getSystemPath()).c_str());
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(SystemAdditionTaskTest, electronicStrucure) {
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  Settings settingsSys2 = sys2->getSettings();
  settingsSys2.name = "TMP_SUPERSYSTEM";
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settingsSys2);
  SystemAdditionTask<RESTRICTED> add(supersystem, {sys1, sys2});
  add.run();
  // The occupied orbitals of the water monomers will be sorted to the new basis set.
  // The coefficients have to be identical. Virtual orbitals will be set to zero.
  CoefficientMatrix<RESTRICTED> superCoeffs = supersystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  CoefficientMatrix<RESTRICTED> coeffsSys1 = sys1->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  CoefficientMatrix<RESTRICTED> coeffsSys2 = sys2->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  unsigned int nBasisFuncWaterMono = sys1->getBasisController()->getNBasisFunctions();

  double diff = (superCoeffs.topLeftCorner(nBasisFuncWaterMono, 5) - coeffsSys1.topLeftCorner(nBasisFuncWaterMono, 5))
                    .array()
                    .abs()
                    .sum();
  EXPECT_NEAR(0.0, diff, 1e-12);
  diff = (superCoeffs.block(nBasisFuncWaterMono, 5, nBasisFuncWaterMono, 5) - coeffsSys2.topLeftCorner(nBasisFuncWaterMono, 5))
             .array()
             .abs()
             .sum();
  EXPECT_NEAR(0.0, diff, 1e-12);
  diff = (superCoeffs.rightCols(2 * nBasisFuncWaterMono - 10)).array().abs().sum();
  EXPECT_NEAR(0.0, diff, 1e-12);

  std::string supersystemName = supersystem->getSystemName();
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".settings").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".xyz").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".energies.res").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".orbs.res.h5").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".dmat.res.h5").c_str());
  std::remove((supersystem->getSystemPath()).c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(SystemAdditionTaskTest, onlyAtoms) {
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  Settings settingsSys2 = sys2->getSettings();
  settingsSys2.basis.label = "DEF2-SVP";
  settingsSys2.name = "TMP_SUPERSYSTEM";
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settingsSys2);
  SystemAdditionTask<RESTRICTED> add(supersystem, {sys1, sys2});
  add.settings.addOccupiedOrbitals = false;
  add.run();

  EXPECT_EQ(6, supersystem->getGeometry()->getNAtoms());

  std::string supersystemName = supersystem->getSystemName();
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".settings").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".xyz").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".energies.res").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".orbs.res.h5").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".dmat.res.h5").c_str());
  std::remove((supersystem->getSystemPath()).c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
TEST_F(SystemAdditionTaskTest, checkGeometryAndCharge) {
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  Settings settingsSys2 = sys2->getSettings();
  settingsSys2.name = "TMP_SUPERSYSTEM";
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settingsSys2);
  SystemAdditionTask<RESTRICTED> add(supersystem, {sys1, sys2});
  add.run();

  SystemAdditionTask<RESTRICTED> add2(supersystem, {sys1, sys2});
  add2.settings.checkSuperGeom = true;
  add2.settings.checkSuperCharge = true;
  EXPECT_NO_THROW(add2.run());

  std::string supersystemName = supersystem->getSystemName();
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".settings").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".xyz").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".energies.res").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".orbs.res.h5").c_str());
  std::remove((supersystem->getSystemPath() + "/" + supersystemName + ".dmat.res.h5").c_str());
  std::remove((supersystem->getSystemPath()).c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
} /*namespace Serenity*/
