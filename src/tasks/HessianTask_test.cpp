/**
 * @file HessianTask_test.cpp
 *
 * @date Feb 23, 2017
 * @author Kevin Klahr
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
#include "tasks/HessianTask.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "io/HDF5.h"
#include "math/Matrix.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class HessianTaskTest
 * @brief Sets everything up for the tests of HessianTask.h/.cpp .
 *
 */
class HessianTaskTest : public ::testing::Test {
 protected:
  HessianTaskTest() {
  }

  virtual ~HessianTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test HessianTaskTest
 * @brief Tests the HF hessian
 */
TEST_F(HessianTaskTest, h2Hf_res) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.scf.energyThreshold = 1e-6;
  settings.scf.rmsdThreshold = 1e-6;
  settings.basis.label = "DEF2-SVP";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  std::vector<std::shared_ptr<SystemController>> activesystems;
  activesystems.push_back(systemController);
  auto hessTask = HessianTask<Options::SCF_MODES::RESTRICTED>(activesystems);
  hessTask.settings.printToFile = false;
  hessTask.run();
  auto path = systemController->getSystemPath();

  HDF5::Filepath name(path + settings.name + ".hess.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "hessian");
  HDF5::dataset_exists(file, "eigenvectors");
  HDF5::dataset_exists(file, "eigenvalues");
  HDF5::dataset_exists(file, "frequencies");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", settings.identifier);

  Matrix<double> hessian;
  HDF5::load(file, "hessian", hessian);

  Eigen::VectorXd eigenValues;
  HDF5::load(file, "eigenvalues", eigenValues);

  Eigen::VectorXd frequencies;
  HDF5::load(file, "frequencies", frequencies);

  file.close();

  EXPECT_NEAR(hessian(0, 0), -0.0034, 1e-3);
  EXPECT_NEAR(hessian(0, 1), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 2), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 3), 0.0034, 1e-3);
  EXPECT_NEAR(hessian(0, 4), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 5), 0.0, 1e-3);
  EXPECT_NEAR(hessian(1, 1), -0.0034, 1e-3);
  EXPECT_NEAR(hessian(2, 2), 0.4140, 1e-3);
  EXPECT_NEAR(hessian(3, 3), -0.0034, 1e-3);

  EXPECT_NEAR(frequencies[0], 4677.733, 1e-0);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test HessianTaskTest
 * @brief Tests the HF hessian for an open-shell system
 */
TEST_F(HessianTaskTest, h2oHf_unres) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.scf.energyThreshold = 1e-10;
  settings.scf.rmsdThreshold = 1e-10;
  settings.basis.label = "DEF2-SVP";
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.charge = 1;
  settings.spin = 1;
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  systemController->setSpin(1);
  systemController->setCharge(1);
  std::vector<std::shared_ptr<SystemController>> activesystems;
  activesystems.push_back(systemController);
  auto hessTask = HessianTask<Options::SCF_MODES::UNRESTRICTED>(activesystems);
  hessTask.settings.numHessStepSize = 1e-3;
  hessTask.settings.printToFile = false;
  hessTask.run();
  auto path = systemController->getSystemPath();

  HDF5::Filepath name(path + settings.name + ".hess.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "rawHessian");
  HDF5::dataset_exists(file, "eigenvectors");
  HDF5::dataset_exists(file, "eigenvalues");
  HDF5::dataset_exists(file, "frequencies");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", settings.identifier);

  Eigen::MatrixXd hessian;
  HDF5::load(file, "rawHessian", hessian);

  Eigen::VectorXd eigenValues;
  HDF5::load(file, "eigenvalues", eigenValues);

  Eigen::VectorXd frequencies;
  HDF5::load(file, "frequencies", frequencies);

  file.close();

  // ORCA 6.0.1 18.02.2025
  Eigen::MatrixXd orca_ref = Eigen::MatrixXd::Zero(9, 9);
  orca_ref << 0.03936, 0.00000, -0.01398, 0.00369, 0.00000, 0.04833, -0.04305, 0.00000, -0.03436, 0.00000, -0.00693,
      0.00000, 0.00000, 0.00516, 0.00000, 0.00000, 0.00177, 0.00000, -0.01398, 0.00000, 0.53809, 0.00449, 0.00000,
      -0.01052, 0.00948, 0.00000, -0.52757, 0.00369, 0.00000, 0.00449, 0.49742, 0.00000, -0.13720, -0.50111, 0.00000,
      0.13271, 0.00000, 0.00516, 0.00000, 0.00000, -0.00693, 0.00000, 0.00000, 0.00177, 0.00000, 0.04833, 0.00000,
      -0.01052, -0.13720, 0.00000, 0.08003, 0.08886, 0.00000, -0.06950, -0.04305, 0.00000, 0.00948, -0.50111, 0.00000,
      0.08886, 0.54416, 0.00000, -0.09835, 0.00000, 0.00177, 0.00000, 0.00000, 0.00177, 0.00000, 0.00000, -0.00354,
      0.00000, -0.03436, 0.00000, -0.52757, 0.13271, 0.00000, -0.06950, -0.09835, 0.00000, 0.59707;

  EXPECT_NEAR((hessian - orca_ref).cwiseAbs().maxCoeff(), 0, 1e-3);
  EXPECT_NEAR(frequencies[0], 1559.941, 1e-0);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test HessianTaskTest
 * @brief Tests the HF hessian
 */
TEST_F(HessianTaskTest, h2Hf_unres) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.scf.energyThreshold = 1e-6;
  settings.scf.rmsdThreshold = 1e-6;
  settings.basis.label = "DEF2-SVP";
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  std::vector<std::shared_ptr<SystemController>> activesystems;
  activesystems.push_back(systemController);
  auto hessTask = HessianTask<Options::SCF_MODES::UNRESTRICTED>(activesystems);
  hessTask.settings.printToFile = false;
  hessTask.run();
  auto path = systemController->getSystemPath();

  HDF5::Filepath name(path + settings.name + ".hess.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "hessian");
  HDF5::dataset_exists(file, "eigenvectors");
  HDF5::dataset_exists(file, "eigenvalues");
  HDF5::dataset_exists(file, "frequencies");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", settings.identifier);

  Matrix<double> hessian;
  HDF5::load(file, "hessian", hessian);

  Eigen::VectorXd eigenValues;
  HDF5::load(file, "eigenvalues", eigenValues);

  Eigen::VectorXd frequencies;
  HDF5::load(file, "frequencies", frequencies);

  file.close();

  EXPECT_NEAR(hessian(0, 0), -0.0034, 1e-3);
  EXPECT_NEAR(hessian(0, 1), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 2), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 3), 0.0034, 1e-3);
  EXPECT_NEAR(hessian(0, 4), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 5), 0.0, 1e-3);
  EXPECT_NEAR(hessian(1, 1), -0.0034, 1e-3);
  EXPECT_NEAR(hessian(2, 2), 0.4140, 1e-3);
  EXPECT_NEAR(hessian(3, 3), -0.0034, 1e-3);

  EXPECT_NEAR(frequencies[0], 4677.733, 1e-0);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test HessianTaskTest
 * @brief Tests the FaT hessian
 */
TEST_F(HessianTaskTest, h2FaT_res) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  settings.scf.energyThreshold = 1e-6;
  settings.scf.rmsdThreshold = 1e-6;
  settings.basis.label = "STO-3G";
  auto activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  settings.name = "HessianTaskTest_h2FaT_env";
  auto activeSystem2 =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, settings);
  auto hessTask = HessianTask<Options::SCF_MODES::RESTRICTED>({activeSystem, activeSystem2});
  hessTask.settings.printToFile = false;
  hessTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  hessTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  hessTask.run();
  auto path = activeSystem->getSystemPath();

  HDF5::Filepath name(path + "SomeTestSystem.hess.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "hessian");
  HDF5::dataset_exists(file, "eigenvectors");
  HDF5::dataset_exists(file, "eigenvalues");
  HDF5::dataset_exists(file, "frequencies");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", settings.identifier);

  Matrix<double> hessian;
  HDF5::load(file, "hessian", hessian);

  Eigen::VectorXd eigenValues;
  HDF5::load(file, "eigenvalues", eigenValues);

  Eigen::VectorXd frequencies;
  HDF5::load(file, "frequencies", frequencies);

  file.close();

  EXPECT_NEAR(hessian(0, 0), -0.0048920487988684666, 1e-3);
  EXPECT_NEAR(hessian(0, 1), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 2), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 3), -0.0073323618703177917, 1e-3);
  EXPECT_NEAR(hessian(0, 4), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 5), 0.0, 1e-3);
  EXPECT_NEAR(hessian(1, 1), -0.0048920490632923452, 1e-3);
  EXPECT_NEAR(hessian(2, 2), 0.49347353116360243, 1e-3);
  EXPECT_NEAR(hessian(3, 3), -0.57724820590568904, 1e-3);

  // For some reason, even GGA funcationals use the ERIPotential functions for
  // the gradient. These were changed at 14.09.2020 to use the RI approximation
  // if the corresponding fitting setting is set. Thus, these frequencies changed
  // slightly. Before, the full four-center integrals were used!
  EXPECT_NEAR(frequencies[0], -10409.978397631383, 1e-1);
  EXPECT_NEAR(frequencies[5], 8669.932458869529, 1e-1);

  EXPECT_EQ(0, std::remove(((std::string) "WARNING").c_str()));
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(activeSystem);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(activeSystem2);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test HessianTaskTest
 * @brief Tests the FaT hessian
 */
TEST_F(HessianTaskTest, h2FaT_unres) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  settings.scf.energyThreshold = 1e-6;
  settings.scf.rmsdThreshold = 1e-6;
  settings.basis.label = "STO-3G";
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  auto activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  settings.name = "HessianTaskTest_h2FaT_env";
  auto activeSystem2 =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, settings);
  auto hessTask = HessianTask<Options::SCF_MODES::UNRESTRICTED>({activeSystem, activeSystem2});
  hessTask.settings.printToFile = false;
  hessTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  hessTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  hessTask.run();
  auto path = activeSystem->getSystemPath();

  HDF5::Filepath name(path + "SomeTestSystem.hess.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "hessian");
  HDF5::dataset_exists(file, "eigenvectors");
  HDF5::dataset_exists(file, "eigenvalues");
  HDF5::dataset_exists(file, "frequencies");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", settings.identifier);

  Matrix<double> hessian;
  HDF5::load(file, "hessian", hessian);

  Eigen::VectorXd eigenValues;
  HDF5::load(file, "eigenvalues", eigenValues);

  Eigen::VectorXd frequencies;
  HDF5::load(file, "frequencies", frequencies);

  file.close();

  EXPECT_NEAR(hessian(0, 0), -0.0048920487988684666, 1e-3);
  EXPECT_NEAR(hessian(0, 1), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 2), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 3), -0.0073323618703177917, 1e-3);
  EXPECT_NEAR(hessian(0, 4), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 5), 0.0, 1e-3);
  EXPECT_NEAR(hessian(1, 1), -0.0048920490632923452, 1e-3);
  EXPECT_NEAR(hessian(2, 2), 0.49347353116360243, 1e-3);
  EXPECT_NEAR(hessian(3, 3), -0.57724820590568904, 1e-3);

  // For some reason, even GGA funcationals use the ERIPotential functions for
  // the gradient. These were changed at 14.09.2020 to use the RI approximation
  // if the corresponding fitting setting is set. Thus, these frequencies changed
  // slightly. Before, the full four-center integrals were used!
  EXPECT_NEAR(frequencies[0], -10409.978397631383, 1e-1);
  EXPECT_NEAR(frequencies[5], 8669.932458869529, 1e-1);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(activeSystem);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(activeSystem2);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test HessianTaskTest
 * @brief Tests the active system FaT hessian
 */
TEST_F(HessianTaskTest, h2actFaT_res) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  settings.scf.energyThreshold = 1e-6;
  settings.scf.rmsdThreshold = 1e-6;
  settings.basis.label = "STO-3G";
  auto activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  settings.name = "HessianTaskTest_h2actFaT_env";
  auto environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, settings);
  auto hessTask = HessianTask<Options::SCF_MODES::RESTRICTED>({activeSystem}, {environmentSystem});
  hessTask.settings.printToFile = false;
  hessTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  hessTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  hessTask.run();
  auto path = activeSystem->getSystemPath();

  HDF5::Filepath name(path + "SomeTestSystem.hess.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "hessian");
  HDF5::dataset_exists(file, "eigenvectors");
  HDF5::dataset_exists(file, "eigenvalues");
  HDF5::dataset_exists(file, "frequencies");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", settings.identifier);

  Matrix<double> hessian;
  HDF5::load(file, "hessian", hessian);

  Eigen::VectorXd eigenValues;
  HDF5::load(file, "eigenvalues", eigenValues);

  Eigen::VectorXd frequencies;
  HDF5::load(file, "frequencies", frequencies);

  file.close();

  EXPECT_NEAR(hessian(0, 0), -0.0046669480843530304, 1e-3);
  EXPECT_NEAR(hessian(0, 1), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 2), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 3), -0.0067836026842090806, 1e-3);
  EXPECT_NEAR(hessian(0, 4), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 5), 0.0, 1e-3);
  EXPECT_NEAR(hessian(1, 1), -0.0046669480869835243, 1e-3);
  EXPECT_NEAR(hessian(2, 2), 0.49433100937944013, 1e-3);
  EXPECT_NEAR(hessian(3, 3), -0.58142462724513666, 1e-3);

  EXPECT_NEAR(frequencies[0], 6919.5972220852491, 1e-0);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(activeSystem);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(environmentSystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test HessianTaskTest
 * @brief Tests the active system FaT hessian
 */
TEST_F(HessianTaskTest, h2actFaT_unres) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  settings.scf.energyThreshold = 1e-6;
  settings.scf.rmsdThreshold = 1e-6;
  settings.basis.label = "STO-3G";
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  auto activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  settings.name = "HessianTaskTest_h2actFaT_env";
  auto environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, settings);
  auto hessTask = HessianTask<Options::SCF_MODES::UNRESTRICTED>({activeSystem}, {environmentSystem});
  hessTask.settings.printToFile = false;
  hessTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  hessTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  hessTask.run();
  auto path = activeSystem->getSystemPath();

  HDF5::Filepath name(path + "SomeTestSystem.hess.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "hessian");
  HDF5::dataset_exists(file, "eigenvectors");
  HDF5::dataset_exists(file, "eigenvalues");
  HDF5::dataset_exists(file, "frequencies");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", settings.identifier);

  Matrix<double> hessian;
  HDF5::load(file, "hessian", hessian);

  Eigen::VectorXd eigenValues;
  HDF5::load(file, "eigenvalues", eigenValues);

  Eigen::VectorXd frequencies;
  HDF5::load(file, "frequencies", frequencies);

  file.close();

  EXPECT_NEAR(hessian(0, 0), -0.0046669480843530304, 1e-3);
  EXPECT_NEAR(hessian(0, 1), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 2), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 3), -0.0067836026842090806, 1e-3);
  EXPECT_NEAR(hessian(0, 4), 0.0, 1e-3);
  EXPECT_NEAR(hessian(0, 5), 0.0, 1e-3);
  EXPECT_NEAR(hessian(1, 1), -0.0046669480869835243, 1e-3);
  EXPECT_NEAR(hessian(2, 2), 0.49433100937944013, 1e-3);
  EXPECT_NEAR(hessian(3, 3), -0.58142462724513666, 1e-3);

  EXPECT_NEAR(frequencies[0], 6919.5972220852491, 1e-0);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(activeSystem);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(environmentSystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

}; /* NameSpace Serenity */
