/**
 * @file PCMPotential.cpp
 *
 * @author Moritz Bensberg
 * @date Apr 28, 2020
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
#include "potentials/PCMPotential.h"
/* Include Serenity Internal Headers */
#include "data/grid/ElectrostaticPotentialOnGridController.h" //getDensityMatrixController()
#include "data/matrices/DensityMatrixController.h"            //getDensityMatrix()
#include "geometry/Geometry.h"                                //getNAtoms()
#include "integrals/wrappers/Libint.h"                        //compute1eInts()
#include "misc/SerenityError.h"                               //throw SerenityError
#include "misc/Timing.h"                                      //Timings.
#include "settings/PCMSettings.h"                             //Member .use.
#include "solvation/ContinuumModel.h"                         //getPCMCharges() and getActiveEnergy().

namespace Serenity {

template<Options::SCF_MODES SCFMode>
PCMPotential<SCFMode>::PCMPotential(const PCMSettings& pcmSettings, std::shared_ptr<BasisController> basisController,
                                    std::shared_ptr<const Geometry> geometry,
                                    std::shared_ptr<MolecularSurfaceController> molecularSurface,
                                    std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>> activePotential,
                                    std::vector<std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>> environmentPotentials)
  : Potential<SCFMode>(basisController), _geometry(geometry), _activePotential(activePotential) {
  if (pcmSettings.use) {
    if (!molecularSurface || !activePotential)
      throw SerenityError("ERROR: PCM potential construction for a non-existing system.");
    _continuumModel =
        std::make_shared<ContinuumModel<SCFMode>>(pcmSettings, molecularSurface, activePotential, environmentPotentials);
    activePotential->getDensityMatrixController()->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
    for (auto envPotential : environmentPotentials)
      envPotential->getDensityMatrixController()->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  }
  else {
    _continuumModel = nullptr;
  }
  basisController->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  auto& libint = Libint::getInstance();
  libint.keepEngines(libint2::Operator::nuclear, 0, 2);
}

template<Options::SCF_MODES SCFMode>
PCMPotential<SCFMode>::~PCMPotential() {
  auto& libint = Libint::getInstance();
  libint.freeEngines(libint2::Operator::nuclear, 0, 2);
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd PCMPotential<SCFMode>::integrateToFockMatrix(const GridPotential<RESTRICTED>& pcmCharges) {
  Timings::takeTime(" Tech. - PCM Charge Integration");
  Eigen::MatrixXd result = _activePotential->integrateFockMatrix(pcmCharges);
  Timings::timeTaken(" Tech. - PCM Charge Integration");
  return result;
}

template<Options::SCF_MODES SCFMode>
double PCMPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  (void)P;
  if (!_continuumModel)
    return 0.0;
  return _continuumModel->getActivePCMEnergy();
}

template<Options::SCF_MODES SCFMode>
double PCMPotential<SCFMode>::getTotalEnergy() {
  if (!_continuumModel)
    return 0.0;
  return _continuumModel->getTotalPCMEnergy();
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& PCMPotential<SCFMode>::getMatrix() {
  if (!_potential) {
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    if (_continuumModel) {
      Eigen::MatrixXd pcmFockContribution = integrateToFockMatrix(_continuumModel->getPCMCharges());
      for_spin(pot) {
        pot_spin.setZero();
        pot_spin += pcmFockContribution;
      };
    }
    else {
      for_spin(pot) {
        pot_spin.setZero();
      };
    }
  }
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd PCMPotential<SCFMode>::getGeomGradients() {
  Eigen::MatrixXd gradient(_geometry->getNAtoms(), 3);
  gradient.setZero();
  if (!_continuumModel)
    return gradient;
  throw SerenityError("Geometrical gradients for PCM are not implemented yet.");
  return gradient;
}

template class PCMPotential<Options::SCF_MODES::RESTRICTED>;
template class PCMPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
