/**
 * @file BasisSetTruncationTask.cpp
 *
 * @date Sep 24, 2018
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

/* Include Class Header*/
#include "tasks/BasisSetTruncationTask.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h" //Shell-wise truncation.
#include "basis/CustomBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h" //Overlap integrals
#include "misc/BasisSetTruncationAlgorithms.h"
#include "misc/SystemSplittingTools.h"
#include "misc/WarningTracker.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
BasisSetTruncationTask<SCFMode>::BasisSetTruncationTask(std::shared_ptr<SystemController> system) : _system(system) {
}

template<Options::SCF_MODES SCFMode>
void BasisSetTruncationTask<SCFMode>::run() {
  assert(_system->hasElectronicStructure<SCFMode>() && "Error: The system needs to have an electronic structure" &&
         " for any type of basis set truncation performed with the basisSetTruncationTask. Check your input.");
  // Search for dummy atoms and give warning if non are found.
  const auto activeAtoms = _system->getGeometry()->getAtoms();
  bool hasDummy = false;
  for (const auto& atom : activeAtoms) {
    if (atom->isDummy()) {
      hasDummy = true;
      break;
    } // of atom->isDummy()
  }   // for atom
  if (!hasDummy) {
    WarningTracker::printWarning(
        (std::string) "The basis set truncation task was called for a system with no dummy atoms!\n" +
            (std::string) " Only shells centered on dummy atoms can be truncated!",
        iOOptions.printSCFCycleInfo);
  }
  // We need to explicitly copy the basis controller, since it may be reset during basis
  // set truncation. Since we would like to project the density matrix below. We need to
  // keep the basis set information.
  Basis basis = _system->getBasisController()->getBasis();
  auto oldBasis = std::make_shared<CustomBasisController>(basis, "nonTruncated");
  // active density matrix
  DensityMatrix<SCFMode> newActiveDensityMatrix(oldBasis);
  auto oldActiveDensityMatrix = _system->getElectronicStructure<SCFMode>()->getDensityMatrix();
  const auto nCoreOrbitals = _system->template getActiveOrbitalController<SCFMode>()->getNCoreOrbitals();
  // Explcitly copy the 'old' density matrix.
  for_spin(newActiveDensityMatrix, oldActiveDensityMatrix) {
    newActiveDensityMatrix_spin.setZero();
    newActiveDensityMatrix_spin += Eigen::MatrixXd(oldActiveDensityMatrix_spin);
  };
  // truncate the basis set
  BasisSetTruncationAlgorithms truncAlgorithm(_system);
  truncAlgorithm.truncateBasis(settings.truncAlgorithm, settings.truncationFactor,
                               std::make_shared<DensityMatrix<Options::SCF_MODES::RESTRICTED>>(oldActiveDensityMatrix.total()),
                               settings.netThreshold);

  // project former density matrix into the new basis set.
  DensityMatrix<SCFMode> projectedDens = SystemSplittingTools<SCFMode>::projectMatrixIntoNewBasis(
      newActiveDensityMatrix, _system->getBasisController(),
      std::make_shared<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>(
          _system->getOneElectronIntegralController()->getOverlapIntegrals()));

  // set new electronic structure
  auto newES = std::make_shared<ElectronicStructure<SCFMode>>(_system->getOneElectronIntegralController(),
                                                              _system->getNOccupiedOrbitals<SCFMode>(), nCoreOrbitals);
  newES->getDensityMatrixController()->setDensityMatrix(projectedDens);
  _system->setElectronicStructure<SCFMode>(newES);
}
template class BasisSetTruncationTask<Options::SCF_MODES::RESTRICTED>;
template class BasisSetTruncationTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
