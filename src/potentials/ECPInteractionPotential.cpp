/**
 * @file ECPInteractionPotential.cpp
 *
 * @date Jun 11, 2018
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
#include "potentials/ECPInteractionPotential.h"
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "geometry/Geometry.h"
#include "potentials/EffectiveCorePotential.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ECPInteractionPotential<SCFMode>::ECPInteractionPotential(std::shared_ptr<SystemController> actSystem,
                                                          std::vector<std::shared_ptr<Atom>> actAtoms,
                                                          std::vector<std::shared_ptr<Atom>> envAtoms,
                                                          std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensities,
                                                          std::shared_ptr<BasisController> basis)
  : Potential<SCFMode>(basis), _actSystem(actSystem), _actAtoms(actAtoms), _envAtoms(envAtoms), _envDensities(envDensities) {
  // Build potential for interaction between active density and env ECPs.
  _envActDensECP = std::make_shared<EffectiveCorePotential<SCFMode>>(actSystem, _envAtoms, this->_basis);
  // Build potential for interaction between env density and act ECPs.
  for (const auto& envDensity : _envDensities) {
    auto newECP = std::make_shared<EffectiveCorePotential<SCFMode>>(actSystem, _actAtoms,
                                                                    envDensity->getDensityMatrix().getBasisController());
    _actEnvDensECPs.push_back(newECP);
  }
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& ECPInteractionPotential<SCFMode>::getMatrix() {
  return _envActDensECP->getMatrix();
}

template<Options::SCF_MODES SCFMode>
double ECPInteractionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  double energy = 0.0;
  // Energy contribution from the interaction between active density and environment ECPs.
  energy += _envActDensECP->getEnergy(P);
  // Energy contribution from the interaction between environment density and active ECPs.
  unsigned int counter = 0;
  for (const auto& ecp : _actEnvDensECPs) {
    // Calculate the energy.
    energy += ecp->getEnergy(_envDensities[counter]->getDensityMatrix());
    ++counter;
  }
  return energy;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd ECPInteractionPotential<SCFMode>::getGeomGradients() {
  auto actSystem = _actSystem.lock();
  Eigen::MatrixXd gradient(actSystem->getGeometry()->getAtoms().size(), 3);
  gradient.setZero();
  gradient += _envActDensECP->getGeomGradients();
  for (const auto& ecp : _actEnvDensECPs)
    gradient += ecp->getGeomGradients();
  return gradient;
}

template class ECPInteractionPotential<Options::SCF_MODES::RESTRICTED>;
template class ECPInteractionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
