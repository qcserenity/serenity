/**
 * @file EffectiveCorePotential.cpp
 *
 * @date May 24, 2018
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
#include "potentials/EffectiveCorePotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libecpint.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
EffectiveCorePotential<SCFMode>::EffectiveCorePotential(std::shared_ptr<SystemController> system,
                                                        std::vector<std::shared_ptr<Atom>> atoms,
                                                        std::shared_ptr<BasisController> basis)
  : Potential<SCFMode>(basis), _system(system), _atoms(atoms) {
  _notZero = false;
  for (const auto& atom : _atoms) {
    if (atom->usesECP())
      _notZero = true;
  }
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& EffectiveCorePotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -     1e-Int Pot.");
  if (!_potential) {
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& f = *_potential;
    // Perform the Libecpint call only if there are atoms with an ECP.
    if (_notZero) {
      FockMatrix<Options::SCF_MODES::RESTRICTED> ecp(this->_basis);
      ecp = Libecpint::computeECPIntegrals(this->_basis, _atoms);
      for_spin(f) {
        f_spin += ecp;
      };
    }
  }
  Timings::timeTaken("Active System -     1e-Int Pot.");
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
double EffectiveCorePotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (!_potential)
    this->getMatrix();
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };
  return energy;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd EffectiveCorePotential<SCFMode>::getGeomGradients() {
  auto system = _system.lock();
  if (_notZero) {
    auto densMatrix = system->template getElectronicStructure<SCFMode>()->getDensityMatrix().total();
    auto acb = std::dynamic_pointer_cast<AtomCenteredBasisController>(this->_basis);
    return Libecpint::computeECPGradientContribution(acb, _atoms, densMatrix);
  }
  else {
    return Eigen::MatrixXd::Zero(system->getGeometry()->getAtoms().size(), 3);
  }
}

template class EffectiveCorePotential<Options::SCF_MODES::RESTRICTED>;
template class EffectiveCorePotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
