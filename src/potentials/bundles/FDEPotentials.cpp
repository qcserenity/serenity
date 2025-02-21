/**
 * @file FDEPotentials.cpp
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
#include "potentials/bundles/FDEPotentials.h"
/* Include Serenity Internal Headers */
#include "energies/EnergyContributions.h"
#include "potentials/NAddFuncPotential.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FDEPotentials<SCFMode>::FDEPotentials(std::shared_ptr<PotentialBundle<SCFMode>> activeSystemPot,
                                      std::shared_ptr<PotentialBundle<SCFMode>> esiPot,
                                      std::vector<std::shared_ptr<NAddFuncPotential<SCFMode>>> naddXC,
                                      std::vector<std::shared_ptr<Potential<SCFMode>>> naddKin,
                                      std::shared_ptr<Potential<SCFMode>> ecp, std::shared_ptr<Potential<SCFMode>> pcm)
  : _activeSystemPot(activeSystemPot), _esiPot(esiPot), _naddXC(naddXC), _naddKin(naddKin), _ecp(ecp), _pcm(pcm) {
  assert(_activeSystemPot);
  assert(_esiPot);
  assert(_ecp);
  assert(_pcm);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> FDEPotentials<SCFMode>::getFockMatrix(const DensityMatrix<SCFMode>& P,
                                                          std::shared_ptr<EnergyComponentController> energies) {
  // Energies:
  double eKin = 0.0;
  double eXC = 0.0;
  // esi
  auto F = _esiPot->getFockMatrix(P, energies);
  // nadd XC
  for (auto naddxc : _naddXC) {
    F += naddxc->getMatrix();
    eXC += naddxc->getEnergy(P);
  }
  // nadd kin
  for (auto naddkin : _naddKin) {
    F += naddkin->getMatrix();
    eKin += naddkin->getEnergy(P);
  }
  // active system
  F += _activeSystemPot->getFockMatrix(P, energies);
  // ecp
  F += _ecp->getMatrix();
  // pcm
  F += _pcm->getMatrix();

  double eECP = _ecp->getEnergy(P);
  double ePCM = _pcm->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_XC, eXC);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_KINETIC, eKin);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_ECP_ENV_ENERGY, eECP);
  // Overwrite the PCM energy contribution set by the active system potential.
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::SOLVATION_ENERGY, ePCM);
  return F;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd FDEPotentials<SCFMode>::getGradients() {
  auto gradients = _activeSystemPot->getGradients();
  gradients += _esiPot->getGradients();
  // nadd kin
  for (auto naddkin : _naddKin) {
    gradients += naddkin->getGeomGradients();
  }
  // nadd XC
  for (auto naddxc : _naddXC) {
    gradients += naddxc->getGeomGradients();
  }
  gradients += _ecp->getGeomGradients();
  gradients += _pcm->getGeomGradients();
  return gradients;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<PotentialBundle<SCFMode>> FDEPotentials<SCFMode>::getActiveSystemPotentials() {
  return _activeSystemPot;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<PotentialBundle<SCFMode>> FDEPotentials<SCFMode>::getESIPotentials() {
  return _esiPot;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<NAddFuncPotential<SCFMode>>> FDEPotentials<SCFMode>::getNaddXCPotentials() {
  return _naddXC;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<Potential<SCFMode>>> FDEPotentials<SCFMode>::getNaddKinPotentials() {
  return _naddKin;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Potential<SCFMode>> FDEPotentials<SCFMode>::getECPInteractionPotential() {
  return _ecp;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Potential<SCFMode>> FDEPotentials<SCFMode>::getPCMPotential() {
  return _pcm;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd FDEPotentials<SCFMode>::getPointChargeGradients() {
  return _activeSystemPot->getPointChargeGradients();
}

template class FDEPotentials<Options::SCF_MODES::RESTRICTED>;
template class FDEPotentials<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
