/**
 * @file ExportCavityTask_test.cpp
 *
 * @date   Nov 25, 2024
 * @author Lukas Paetow
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
#include "tasks/ExportCavityTask.h"
#include "data/ElectronicStructure.h"
#include "geometry/MolecularSurfaceController.h"
#include "io/HDF5.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <fstream>
#include <string>

namespace Serenity {

class ExportCavityTaskTest : public ::testing::Test {
 protected:
  ExportCavityTaskTest() {
  }

  virtual ~ExportCavityTaskTest() = default;

  inline bool fileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    if (f.good()) {
      f.close();
      return true;
    }
    else {
      f.close();
      return false;
    }
  }

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(ExportCavityTaskTest, h2) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  Settings settings = systemController->getSettings();
  settings.pcm.use = true;
  settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  settings.pcm.solverType = Options::PCM_SOLVER_TYPES::CPCM;
  settings.pcm.cavityFormation = true;
  settings.pcm.saveCharges = true;
  systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto es = systemController->getElectronicStructure<RESTRICTED>();

  ExportCavityTask task(systemController);
  task.run();

  std::string filepath = systemController->getSystemPath() + "CavityData.h5";
  HDF5::H5File file(filepath, H5F_ACC_RDONLY);
  Eigen::Matrix3Xd gridPoints;
  Eigen::VectorXd weights;
  Eigen::Matrix3Xd normalVectors;
  std::string label;
  Eigen::VectorXi pointToSphereMappingVector;
  Eigen::VectorXi firstIndices;
  Eigen::VectorXi lastIndices;
  std::string sphereString = "nSpheres";
  unsigned int nSpheres = HDF5::load_attribute<unsigned int>(file, sphereString);
  Eigen::VectorXd radii(nSpheres);
  Eigen::VectorXi angularMomenta(nSpheres);
  Eigen::MatrixXd coordMatrix(3, nSpheres);
  Eigen::VectorXi sphereTypeVector(nSpheres);

  HDF5::load_scalar_attribute(file, "label", label);
  EXPECT_EQ(label, "DELLEY");
  HDF5::load(file, "gridPoints", gridPoints);
  EXPECT_NEAR(2.74884, gridPoints(0, 0), 1e-4);
  HDF5::load(file, "weights", weights);
  EXPECT_NEAR(1.20575, weights(0, 0), 1e-4);
  HDF5::load(file, "normalVectors", normalVectors);
  if (normalVectors(0, 0) < 0.5)
    EXPECT_NEAR(0.00000, normalVectors(0, 0), 1e-4);
  else
    EXPECT_NEAR(1.00000, normalVectors(0, 0), 1e-4);
  HDF5::load(file, "pointToSphereMappingvector", pointToSphereMappingVector);
  EXPECT_TRUE((pointToSphereMappingVector[0] == 0.0) || (pointToSphereMappingVector[0] == 1.0));
  HDF5::load(file, "radii", radii);
  EXPECT_NEAR(2.72121, radii[0], 1e-4);
  HDF5::load(file, "angularMomenta", angularMomenta);
  EXPECT_EQ(angularMomenta[0], 4);
  HDF5::load(file, "coordMatrix", coordMatrix);
  EXPECT_EQ(coordMatrix(0, 0), 0);
  HDF5::load(file, "sphereTypeVector", sphereTypeVector);
  EXPECT_EQ(sphereTypeVector[0], 0);
  HDF5::load(file, "firstIndices", firstIndices);
  EXPECT_TRUE((firstIndices[0] == 0) || (firstIndices[0] == 37));
  HDF5::load(file, "lastIndices", lastIndices);
  EXPECT_TRUE((lastIndices[0] == 74) || (lastIndices[0] == 37));
  file.close();

  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/CavityData.h5").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/VDWCavityData.h5").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/PCMChargesData.h5").c_str()));
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
