/**
 * @file ElectronicStructureCopyTask_test.cpp
 *
 * @author Moritz Bensberg
 * @date Feb 13, 2020
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
#include "tasks/ElectronicStructureCopyTask.h"        //To be tested.
#include "data/ElectronicStructure.h"                 //Get density matrix controller.
#include "data/OrbitalController.h"                   //Access to the coefficients matrix.
#include "data/matrices/DensityMatrixController.h"    //Density matrix construction and check for number of electrons.
#include "integrals/OneElectronIntegralController.h"  //Overlap integrals.
#include "settings/Settings.h"                        //Get system with altered settings.
#include "system/SystemController.h"                  //Test systems.
#include "tasks/ScfTask.h"                            //Run SCF.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test resources.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework

namespace Serenity {

class ElectronicStructureCopyTaskTest : public ::testing::Test {
 protected:
  ElectronicStructureCopyTaskTest() {
  }

  virtual ~ElectronicStructureCopyTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(ElectronicStructureCopyTaskTest, restricted) {
  auto source = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto target = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  ScfTask<RESTRICTED> scf(source);
  scf.run();

  ElectronicStructureCopyTask<RESTRICTED> copyTask(source, {target});
  copyTask.settings.atomFrameIndices = {0, 1, 2};
  copyTask.run();

  DensityMatrix<RESTRICTED> targetDens = target->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  MatrixInBasis<RESTRICTED> targetOverlap = target->getOneElectronIntegralController()->getOverlapIntegrals();
  double nElec = (targetDens.array() * targetOverlap.array()).sum();
  EXPECT_NEAR(10, nElec, 1e-5);
  CoefficientMatrix<RESTRICTED> targetCoefficients = target->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  Eigen::MatrixXd moOverlapMatrix = targetCoefficients.transpose() * targetOverlap * targetCoefficients;
  double diffToIdentity =
      (moOverlapMatrix - Eigen::MatrixXd::Identity(moOverlapMatrix.rows(), moOverlapMatrix.cols())).array().abs().sum();
  EXPECT_NEAR(0.0, diffToIdentity, 1e-4);

  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ElectronicStructureCopyTaskTest, restricted_linear) {
  auto source = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
  auto target = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, true);
  ScfTask<RESTRICTED> scf(source);
  scf.run();

  ElectronicStructureCopyTask<RESTRICTED> copyTask(source, {target});
  copyTask.settings.orthogonalize = true;
  copyTask.run();

  DensityMatrix<RESTRICTED> targetDens = target->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  MatrixInBasis<RESTRICTED> targetOverlap = target->getOneElectronIntegralController()->getOverlapIntegrals();
  double nElec = (targetDens.array() * targetOverlap.array()).sum();
  EXPECT_NEAR(2, nElec, 1e-9);
  CoefficientMatrix<RESTRICTED> targetCoefficients = target->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  Eigen::MatrixXd moOverlapMatrix = targetCoefficients.transpose() * targetOverlap * targetCoefficients;
  double diffToIdentity =
      (moOverlapMatrix - Eigen::MatrixXd::Identity(moOverlapMatrix.rows(), moOverlapMatrix.cols())).array().abs().sum();
  EXPECT_NEAR(0.0, diffToIdentity, 1e-9);

  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ElectronicStructureCopyTaskTest, restricted_orthogonalization) {
  auto source = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto target = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  ScfTask<RESTRICTED> scf(source);
  scf.run();

  ElectronicStructureCopyTask<RESTRICTED> copyTask(source, {target});
  copyTask.settings.orthogonalize = true;
  copyTask.run();

  DensityMatrix<RESTRICTED> targetDens = target->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  MatrixInBasis<RESTRICTED> targetOverlap = target->getOneElectronIntegralController()->getOverlapIntegrals();
  double nElec = (targetDens.array() * targetOverlap.array()).sum();
  EXPECT_NEAR(10, nElec, 1e-5);
  CoefficientMatrix<RESTRICTED> targetCoefficients = target->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  Eigen::MatrixXd moOverlapMatrix = targetCoefficients.transpose() * targetOverlap * targetCoefficients;
  double diffToIdentity =
      (moOverlapMatrix - Eigen::MatrixXd::Identity(moOverlapMatrix.rows(), moOverlapMatrix.cols())).array().abs().sum();
  EXPECT_NEAR(0.0, diffToIdentity, 1e-9);

  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ElectronicStructureCopyTaskTest, restrictedBasisSetChange) {
  auto source = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto target = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  Settings settings = target->getSettings();
  settings.basis.label = "DEF2-SVP";
  target = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, settings);
  ScfTask<RESTRICTED> scf(source);
  scf.run();

  ElectronicStructureCopyTask<RESTRICTED> copyTask(source, {target});
  copyTask.run();

  DensityMatrix<RESTRICTED> targetDens = target->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  MatrixInBasis<RESTRICTED> targetOverlap = target->getOneElectronIntegralController()->getOverlapIntegrals();
  double nElec = (targetDens.array() * targetOverlap.array()).sum();
  EXPECT_NEAR(10, nElec, 5e-3);

  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ElectronicStructureCopyTaskTest, unrestricted) {
  auto source = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto target = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  source->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  target->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  ScfTask<UNRESTRICTED> scf(source);
  scf.run();

  ElectronicStructureCopyTask<UNRESTRICTED> copyTask(source, {target});
  copyTask.run();

  DensityMatrix<UNRESTRICTED> targetDens = target->getElectronicStructure<UNRESTRICTED>()->getDensityMatrix();
  MatrixInBasis<RESTRICTED> targetOverlap = target->getOneElectronIntegralController()->getOverlapIntegrals();
  double nElec = (targetDens.alpha.array() * targetOverlap.array()).sum() * 2.0;
  EXPECT_NEAR(10, nElec, 1e-5);
  CoefficientMatrix<UNRESTRICTED> targetCoefficients = target->getActiveOrbitalController<UNRESTRICTED>()->getCoefficients();
  Eigen::MatrixXd moOverlapMatrix = targetCoefficients.alpha.transpose() * targetOverlap * targetCoefficients.alpha;
  double diffToIdentity =
      (moOverlapMatrix - Eigen::MatrixXd::Identity(moOverlapMatrix.rows(), moOverlapMatrix.cols())).array().abs().sum();
  EXPECT_NEAR(0.0, diffToIdentity, 1e-4);

  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ElectronicStructureCopyTaskTest, restrictedMultiCopy) {
  auto source = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto target1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  auto settings = target1->getSettings();
  auto target2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  auto target3 =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP, settings);
  ScfTask<RESTRICTED> scf(source);
  scf.run();

  ElectronicStructureCopyTask<RESTRICTED> copyTask(source, {target1, target2, target3});
  copyTask.run();

  DensityMatrix<RESTRICTED> targetDens = target2->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  MatrixInBasis<RESTRICTED> targetOverlap = target2->getOneElectronIntegralController()->getOverlapIntegrals();
  double nElec = (targetDens.array() * targetOverlap.array()).sum();
  EXPECT_NEAR(10, nElec, 1e-5);
  CoefficientMatrix<RESTRICTED> targetCoefficients = target2->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  Eigen::MatrixXd moOverlapMatrix = targetCoefficients.transpose() * targetOverlap * targetCoefficients;
  double diffToIdentity =
      (moOverlapMatrix - Eigen::MatrixXd::Identity(moOverlapMatrix.rows(), moOverlapMatrix.cols())).array().abs().sum();
  EXPECT_NEAR(0.0, diffToIdentity, 1e-4);

  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
