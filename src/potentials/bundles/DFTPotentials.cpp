/**
 * @file DFTPotentials.cpp
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
#include "potentials/bundles/DFTPotentials.h"
/* Include Serenity Internal Headers */
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h"
#include "potentials/ERIPotential.h"
#include "potentials/FuncPotential.h"
#include "potentials/HCorePotential.h"
#include "potentials/LRXPotential.h"
#include "potentials/Potential.h"

namespace Serenity {
template<Options::SCF_MODES SCFMode>
DFTPotentials<SCFMode>::DFTPotentials(std::shared_ptr<HCorePotential<SCFMode>> hcore,
                                      std::shared_ptr<Potential<SCFMode>> J, std::shared_ptr<FuncPotential<SCFMode>> Vxc,
                                      std::shared_ptr<Potential<SCFMode>> pcm, std::shared_ptr<const Geometry> geom,
                                      std::shared_ptr<DensityMatrixController<SCFMode>> dMatController,
                                      const double prescreeningThreshold)
  : _h(hcore), _J(J), _Vxc(Vxc), _pcm(pcm), _geom(geom), _dMatController(dMatController), _prescreeningThreshold(prescreeningThreshold) {
  assert(_h);
  assert(_J);
  assert(_Vxc);
  assert(_pcm);
  assert(_geom);
  assert(_dMatController);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> DFTPotentials<SCFMode>::getFockMatrix(const DensityMatrix<SCFMode>& P,
                                                          std::shared_ptr<EnergyComponentController> energies) {
  // one elctron part
  const FockMatrix<SCFMode>& h = _h->getMatrix();
  double eh = _h->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::ONE_ELECTRON_ENERGY, eh);
  // XC part
  auto& V = _Vxc->getMatrix();
  double eV = _Vxc->getEnergy(P);
  // Coulomb part
  const auto& J = _J->getMatrix();
  double eJ = _J->getEnergy(P);
  // PCM Part
  const auto& pcm = _pcm->getMatrix();
  double ePCM = _pcm->getEnergy(P);
  if (!_Vxc->getFunctional().isHybrid()) {
    energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXACT_EXCHANGE, 0.0);
    energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::ELECTRON_ELECTRON_COULOMB, eJ);
  }
  else {
    auto ptr = std::dynamic_pointer_cast<ERIPotential<SCFMode>>(_J);
    double ex = ptr->getXEnergy(P);
    energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXACT_EXCHANGE, ex);
    energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::ELECTRON_ELECTRON_COULOMB, eJ - ex);
  }
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXCHANGE_CORRELATION_NO_HF, eV);
  // NN part
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::NUCLEUS_NUCLEUS_COULOMB, _geom->getCoreCoreRepulsion());
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::SOLVATION_ENERGY, ePCM);
  return h + J + V + pcm;
}

template<Options::SCF_MODES SCFMode>
void DFTPotentials<SCFMode>::replaceFunctionalPotential(std::shared_ptr<FuncPotential<SCFMode>> newVxc) {
  _Vxc = newVxc;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd DFTPotentials<SCFMode>::getGradients() {
  if (_Vxc->getFunctional().isDoubleHybrid())
    throw SerenityError("no gradients for double hybrid functionals implemented yet!");
  auto gradients = _h->getGeomGradients();
  gradients += _J->getGeomGradients();
  gradients += _Vxc->getGeomGradients();
  gradients += _pcm->getGeomGradients();
  return gradients;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd DFTPotentials<SCFMode>::getPointChargeGradients() {
  return _h->getPointChargeGradients();
}

template class DFTPotentials<Options::SCF_MODES::RESTRICTED>;
template class DFTPotentials<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
