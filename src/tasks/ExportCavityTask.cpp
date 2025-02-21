/**
 * @file   ExportCavityTask.cpp
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
#include "tasks/ExportCavityTask.h"
/* Include Serenity Internal Headers */
#include "geometry/MolecularSurface.h"
#include "geometry/MolecularSurfaceController.h"
#include "io/HDF5.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

ExportCavityTask::ExportCavityTask(std::shared_ptr<SystemController> system) : _systemController(system) {
}

void ExportCavityTask::run() {
  std::shared_ptr<MolecularSurfaceController> molecularSurfaceController;
  if (!settings.fdecavity) {
    molecularSurfaceController = _systemController->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE);
    exportCavityFunction(molecularSurfaceController, MOLECULAR_SURFACE_TYPES::ACTIVE);
  }
  else {
    molecularSurfaceController = _systemController->getMolecularSurface(MOLECULAR_SURFACE_TYPES::FDE);
    exportCavityFunction(molecularSurfaceController, MOLECULAR_SURFACE_TYPES::FDE);
  }

  if (_systemController->getSettings().pcm.cavityFormation) {
    molecularSurfaceController = _systemController->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE_VDW);
    exportCavityFunction(molecularSurfaceController, MOLECULAR_SURFACE_TYPES::ACTIVE_VDW);
  }
}

void ExportCavityTask::exportCavityFunction(std::shared_ptr<MolecularSurfaceController> molecularSurfaceController,
                                            MOLECULAR_SURFACE_TYPES surfaceType) {
  const auto& molecularSurface = molecularSurfaceController->getMolecularSurface();
  auto gridPoints = molecularSurfaceController->getGridPoints();
  auto weights = molecularSurfaceController->getWeights();
  auto normalVectors = molecularSurfaceController->getNormalVectors();
  auto label = molecularSurface.getLabel();
  auto pointToSphereMapping = molecularSurfaceController->getPointToSphereMapping(); // pointwise sphere indices.
  Eigen::VectorXi pointToSphereMappingvector(pointToSphereMapping.size());
  for (size_t i = 0; i < pointToSphereMapping.size(); ++i) {
    pointToSphereMappingvector[i] = static_cast<int>(pointToSphereMapping[i]);
  }
  auto gridIndices = molecularSurface.getGridIndicesOfAtoms();
  Eigen::VectorXi firstIndices(gridIndices.size());
  Eigen::VectorXi lastIndices(gridIndices.size());
  for (size_t i = 0; i < gridIndices.size(); ++i) {
    firstIndices(i) = gridIndices[i].first;
    lastIndices(i) = gridIndices[i].second;
  }

  // spheres
  const auto& spheres = molecularSurface.getSpheres();
  unsigned int nSpheres = spheres.size();
  Eigen::VectorXd radii(nSpheres);
  Eigen::VectorXi angularMomenta(nSpheres);
  Eigen::MatrixXd coordMatrix(3, nSpheres); // 3 rows for each element of Vector3d, one column per vector
  Eigen::VectorXi sphereTypeVector(nSpheres);
  std::string sphereTypeString = "";

  unsigned int iterator = 0;
  for (auto itr : spheres) {
    const double& radius = itr.getRadius();
    unsigned int angularMomentum = itr.getAngularMomentum();
    coordMatrix.col(iterator) = itr.getCenterCoords();
    sphereTypeVector(iterator) = static_cast<int>(itr.getSphereType());
    radii(iterator) = radius;
    angularMomenta(iterator) = angularMomentum;
    iterator++;
  }

  std::string filename1 = "";
  if (surfaceType == MOLECULAR_SURFACE_TYPES::ACTIVE || surfaceType == MOLECULAR_SURFACE_TYPES::FDE)
    filename1 = _systemController->getSystemPath() + "CavityData.h5";
  else if (surfaceType == MOLECULAR_SURFACE_TYPES::ACTIVE_VDW)
    filename1 = _systemController->getSystemPath() + "VDWCavityData.h5";
  HDF5::H5File file(filename1, H5F_ACC_TRUNC);
  HDF5::save_scalar_attribute(file, "label", label);
  HDF5::save(file, "gridPoints", gridPoints);
  HDF5::save(file, "weights", weights);
  HDF5::save(file, "normalVectors", normalVectors);
  HDF5::save(file, "firstIndices", firstIndices);
  HDF5::save(file, "lastIndices", lastIndices);
  HDF5::save(file, "pointToSphereMappingvector", pointToSphereMappingvector);
  HDF5::save(file, "radii", radii);
  HDF5::save(file, "angularMomenta", angularMomenta);
  HDF5::save(file, "coordMatrix", coordMatrix);
  HDF5::save(file, "sphereTypeVector", sphereTypeVector);
  std::string nSpheresString = "nSpheres";
  HDF5::save_attribute<unsigned int>(file, nSpheresString, nSpheres);
  if (surfaceType == MOLECULAR_SURFACE_TYPES::ACTIVE || surfaceType == MOLECULAR_SURFACE_TYPES::FDE)
    printf("\n  Saved all quantities of the molecular surface.\n");
  else if (surfaceType == MOLECULAR_SURFACE_TYPES::ACTIVE_VDW)
    printf("\n  Saved all quantities of the auxiliary VDW molecular surface.\n");
  file.close();
}

} // namespace Serenity
