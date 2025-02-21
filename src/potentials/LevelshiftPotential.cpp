/**
 * @file LevelshiftPotential.cpp
 *
 * @date Nov 24, 2016
 * @author Jan Unsleber
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
#include "potentials/LevelshiftPotential.h"
/* Include Serenity Internal Headers */
#include "integrals/wrappers/Libint.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
LevelshiftPotential<SCFMode>::LevelshiftPotential(const std::shared_ptr<BasisController> actBasis,
                                                  std::shared_ptr<DensityMatrixController<SCFMode>> envDensityMatrixController,
                                                  const double levelShiftParameter)
  : Potential<SCFMode>(actBasis),
    _envDMatController(envDensityMatrixController),
    _levelShiftParameter(levelShiftParameter),
    _potential(nullptr) {
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  _envDMatController->getDensityMatrix().getBasisController()->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
};

template<Options::SCF_MODES SCFMode>
double LevelshiftPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (!_potential)
    this->getMatrix();
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += pot_spin.cwiseProduct(P_spin.transpose()).sum();
  };
  return energy;
};

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& LevelshiftPotential<SCFMode>::getMatrix() {
  if (!_potential) {
    auto& libint = Libint::getInstance();
    auto overlap = libint.compute1eInts(LIBINT_OPERATOR::overlap, this->_basis,
                                        _envDMatController->getDensityMatrix().getBasisController());
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& F = *_potential;
    auto envMat = _envDMatController->getDensityMatrix();
    // A factor of one half for restricted to account for the factor of two in the density matrix.
    double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 0.5 : 1.0;
    for_spin(F, envMat) {
      F_spin.setZero();
      F_spin = overlap.transpose() * envMat_spin * overlap;
      F_spin *= scfFactor * _levelShiftParameter;
    };
  }
  return *_potential;
};

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd LevelshiftPotential<SCFMode>::getGeomGradients() {
  throw SerenityError("Gradients from the LevelshiftPotential are not sensible!");

  return Eigen::MatrixXd::Zero(0, 0);
}

template class LevelshiftPotential<Options::SCF_MODES::RESTRICTED>;
template class LevelshiftPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
