/**
 * @file MolecularSurfaceGridFileWriter.cpp
 *
 * @author Moritz Bensberg
 * @date May 19, 2020
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
#include "io/MolecularSurfaceGridFileWriter.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"                       //getAtomType/AtomType
#include "geometry/Geometry.h"                   //getAtoms
#include "geometry/MolecularSurfaceController.h" //Underlying gridController
#include "parameters/Constants.h"                //Convert to angstrom
#include "settings/Settings.h"                   //PCM Settings

namespace Serenity {

MolecularSurfaceGridFileWriter::MolecularSurfaceGridFileWriter(const Settings& settings,
                                                               std::shared_ptr<MolecularSurfaceController> surfaceController)
  : GeneralGridFileWriter(settings, surfaceController), _surfaceGridController(surfaceController) {
}

void MolecularSurfaceGridFileWriter::writeData(std::string filename, const Eigen::VectorXd& data) {
  std::string fullFileName = filename + ".grid.data.xyz";
  std::FILE* file = std::fopen(fullFileName.data(), "a");
  const Eigen::MatrixXd coordinates = BOHR_TO_ANGSTROM * _surfaceGridController->getGridPoints();
  auto atomToPointMapping = _surfaceGridController->getSphereToPointMapping();
  unsigned int nSpheres = atomToPointMapping.size();
  auto geom = _surfaceGridController->getGeometry();
  const auto& atoms = geom->getAtoms();
  for (unsigned int iSphere = 0; iSphere < nSpheres; ++iSphere) {
    std::string symbol = this->_settings.pcm.cavity == Options::PCM_CAVITY_TYPES::DELLEY
                             ? atoms[iSphere]->getAtomType()->getElementSymbol()
                             : "P ";
    for (unsigned int i = atomToPointMapping[iSphere].first; i < atomToPointMapping[iSphere].second; ++i) {
      fprintf(file, "%5s %10g %10g %10g %10g\n", symbol.c_str(), coordinates(0, i), coordinates(1, i),
              coordinates(2, i), data[i]);
    }
  }
  fclose(file);
}

} /* namespace Serenity */
