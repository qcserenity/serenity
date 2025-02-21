/**
 * @file   ImportCavityTask.cpp
 *
 * @date   Nov 08, 2024
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

/* Include Class Header*/
#include "tasks/ImportCavityTask.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "geometry/MolecularSurface.h"
#include "geometry/MolecularSurfaceController.h"
#include "io/HDF5.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

ImportCavityTask::ImportCavityTask(std::shared_ptr<SystemController> system) : _systemController(system) {
}

void ImportCavityTask::run() {
  if (settings.cavityPath == "")
    throw SerenityError("No path given for the molecular surface.");

  std::string filename1 = "CavityData.h5";

  std::string filepath = settings.cavityPath + filename1;
  if (settings.fdecavity) {
    printf("\n  Setting molecular surface (supersystem).\n");
    importCavityFunction(filepath, MOLECULAR_SURFACE_TYPES::FDE);
  }
  else {
    printf("\n  Setting molecular surface (active system).\n");
    importCavityFunction(filepath, MOLECULAR_SURFACE_TYPES::ACTIVE);
  }

  if (settings.vdwcavityPath != "") {
    filename1 = "VDWCavityData.h5";
    filepath = settings.vdwcavityPath + filename1;
    printf("\n  Setting auxiliary molecular surface (VDW).\n");
    importCavityFunction(filepath, MOLECULAR_SURFACE_TYPES::ACTIVE_VDW);
  }
}

void ImportCavityTask::importCavityFunction(std::string filepath, MOLECULAR_SURFACE_TYPES surfaceType) {
  HDF5::Filepath name(filepath);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);

  auto pcmsettings = _systemController->getSettings().pcm;
  pcmsettings.loadedPCM = true;
  pcmsettings.cavityPath = settings.cavityPath;

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
  Eigen::MatrixXd coordMatrix(3, nSpheres); // 3 rows for each element of Vector3d, one column per vector
  Eigen::VectorXi sphereTypeVector(nSpheres);
  std::string sphereTypeString = ""; // might remove this line

  HDF5::attribute_exists(file, "label");
  HDF5::load_scalar_attribute(file, "label", label);
  HDF5::dataset_exists(file, "gridPoints");
  HDF5::load(file, "gridPoints", gridPoints);
  HDF5::dataset_exists(file, "weights");
  HDF5::load(file, "weights", weights);
  HDF5::dataset_exists(file, "normalVectors");
  HDF5::load(file, "normalVectors", normalVectors);
  HDF5::dataset_exists(file, "pointToSphereMappingvector");
  HDF5::load(file, "pointToSphereMappingvector", pointToSphereMappingVector);
  HDF5::dataset_exists(file, "radii");
  HDF5::load(file, "radii", radii);
  HDF5::dataset_exists(file, "angularMomenta");
  HDF5::load(file, "angularMomenta", angularMomenta);
  HDF5::dataset_exists(file, "coordMatrix");
  HDF5::load(file, "coordMatrix", coordMatrix);
  HDF5::dataset_exists(file, "sphereTypeVector");
  HDF5::load(file, "sphereTypeVector", sphereTypeVector);
  HDF5::dataset_exists(file, "firstIndices");
  HDF5::load(file, "firstIndices", firstIndices);
  HDF5::dataset_exists(file, "lastIndices");
  HDF5::load(file, "lastIndices", lastIndices);
  file.close();

  std::vector<std::pair<unsigned int, unsigned int>> sphereIndices;
  for (int i = 0; i < firstIndices.size(); ++i) {
    sphereIndices.push_back(std::pair<unsigned int, unsigned int>(firstIndices(i), lastIndices(i)));
  }

  // construct new molecularsurface
  std::unique_ptr<Eigen::Matrix3Xd> gridPointsptr = std::make_unique<Eigen::Matrix3Xd>(gridPoints);
  std::unique_ptr<Eigen::Matrix3Xd> normalVectorsptr = std::make_unique<Eigen::Matrix3Xd>(normalVectors);

  std::unique_ptr<Eigen::VectorXd> weightsptr = std::make_unique<Eigen::VectorXd>(weights);
  std::vector<unsigned int> pointWiseSphereIndices2(pointToSphereMappingVector.data(),
                                                    pointToSphereMappingVector.data() + pointToSphereMappingVector.size());

  auto r_s = 0.0; // molecularsurfaceController->getCavityFormationProbeRadius(); // private. just get
                  // it later from the solvent compound automatically

  std::vector<Sphere> sphereVector;
  for (unsigned int i = 0; i < nSpheres; i++) {
    double x = coordMatrix(i, 0);
    double y = coordMatrix(i, 1);
    double z = coordMatrix(i, 2);
    Point newCenter(x, y, z);
    auto sphereType = static_cast<SphereType>(sphereTypeVector(i));
    sphereVector.push_back(Sphere(newCenter, radii(i), sphereType, angularMomenta(i)));
  }

  auto newsurface =
      std::make_unique<MolecularSurface>(std::move(gridPointsptr), std::move(weightsptr), std::move(normalVectorsptr),
                                         label, sphereIndices, pointWiseSphereIndices2, r_s, sphereVector);

  std::shared_ptr<Geometry> geom = std::make_shared<Geometry>(*(_systemController->getGeometry()));

  auto newMolecularSurfaceController = std::make_shared<MolecularSurfaceController>(geom, pcmsettings);
  newMolecularSurfaceController->setSurface(std::move(newsurface));

  if (surfaceType == MOLECULAR_SURFACE_TYPES::FDE) {
    _systemController->setMolecularSurface(newMolecularSurfaceController, MOLECULAR_SURFACE_TYPES::FDE);
  }
  else if (surfaceType == MOLECULAR_SURFACE_TYPES::ACTIVE) {
    _systemController->setMolecularSurface(newMolecularSurfaceController, MOLECULAR_SURFACE_TYPES::ACTIVE);
  }
  else if (surfaceType == MOLECULAR_SURFACE_TYPES::ACTIVE_VDW) {
    _systemController->setMolecularSurface(newMolecularSurfaceController, MOLECULAR_SURFACE_TYPES::ACTIVE_VDW);
  }

  newMolecularSurfaceController->printInfo();
}

} // namespace Serenity
