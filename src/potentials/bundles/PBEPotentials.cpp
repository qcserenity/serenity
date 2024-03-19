/**
 * @file PBEPotentials.cpp
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
#include "potentials/bundles/PBEPotentials.h"
/* Include Serenity Internal Headers */
#include "energies/EnergyContributions.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
PBEPotentials<SCFMode>::PBEPotentials(std::shared_ptr<PotentialBundle<SCFMode>> activeSystemPot,
                                      std::shared_ptr<PotentialBundle<SCFMode>> esiPot,
                                      std::shared_ptr<Potential<SCFMode>> naddXC, std::shared_ptr<Potential<SCFMode>> projection,
                                      std::shared_ptr<Potential<SCFMode>> ecp, std::shared_ptr<Potential<SCFMode>> pcm)
  : _activeSystemPot(activeSystemPot), _esiPot(esiPot), _naddXC(naddXC), _projection(projection), _ecp(ecp), _pcm(pcm) {
  assert(_activeSystemPot);
  assert(_esiPot);
  assert(_naddXC);
  assert(_projection);
  assert(_ecp);
  assert(pcm);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> PBEPotentials<SCFMode>::getFockMatrix(const DensityMatrix<SCFMode>& P,
                                                          std::shared_ptr<EnergyComponentController> energies) {
  // nadd kin (sets ElectronicStructure in case of reconstruction)
  auto F = _projection->getMatrix();
  // active system
  F += _activeSystemPot->getFockMatrix(P, energies);
  // esi
  F += _esiPot->getFockMatrix(P, energies);
  // nadd XC
  F += _naddXC->getMatrix();
  double eXC = _naddXC->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_XC, eXC);
  double eKin = _projection->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::PBE_LINEAR_CORRECTION, eKin);
  // ECP
  F += _ecp->getMatrix();
  double eECP = _ecp->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_ECP_ENV_ENERGY, eECP);
  F += _pcm->getMatrix();
  double ePCM = _pcm->getEnergy(P);
  // Overwrite solvation energy set by the active system potential.
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::SOLVATION_ENERGY, ePCM);

  return F;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd PBEPotentials<SCFMode>::getGradients() {
  auto gradients = _activeSystemPot->getGradients();
  gradients += _esiPot->getGradients();
  gradients += _naddXC->getGeomGradients();
  gradients += _projection->getGeomGradients();
  gradients += _ecp->getGeomGradients();
  gradients += _pcm->getGeomGradients();
  return gradients;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd PBEPotentials<SCFMode>::getPointChargeGradients() {
  return _activeSystemPot->getPointChargeGradients();
}

template class PBEPotentials<Options::SCF_MODES::RESTRICTED>;
template class PBEPotentials<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
