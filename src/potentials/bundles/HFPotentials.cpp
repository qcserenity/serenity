/**
 * @file HFPotentials.cpp
 *
 * @date Nov 24, 2016
 * @author Jan Unsleeber
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
#include "potentials/bundles/HFPotentials.h"
/* Include Serenity Internal Headers */
#include "energies/EnergyContributions.h"
#include "potentials/HCorePotential.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HFPotentials<SCFMode>::HFPotentials(std::shared_ptr<HCorePotential<SCFMode>> hcore, std::shared_ptr<Potential<SCFMode>> g,
                                    std::shared_ptr<Potential<SCFMode>> pcm, std::shared_ptr<const Geometry> geom)
  : _h(hcore), _g(g), _pcm(pcm), _geom(geom) {
  assert(_h);
  assert(_g);
  assert(_pcm);
  assert(_geom);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> HFPotentials<SCFMode>::getFockMatrix(const DensityMatrix<SCFMode>& P,
                                                         std::shared_ptr<EnergyComponentController> energies) {
  // one elctron part
  const FockMatrix<SCFMode>& h = _h->getMatrix();
  double eh = _h->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::ONE_ELECTRON_ENERGY, eh);

  // two electron part
  const FockMatrix<SCFMode>& g = _g->getMatrix();
  double eg = _g->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::HF_TWO_ELECTRON_ENERGY, eg);

  // PCM part
  const FockMatrix<SCFMode>& pcm = _pcm->getMatrix();
  double ePCM = _pcm->getEnergy(P);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::SOLVATION_ENERGY, ePCM);

  // nn part
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::NUCLEUS_NUCLEUS_COULOMB, _geom->getCoreCoreRepulsion());
  return g + h + pcm;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd HFPotentials<SCFMode>::getGradients() {
  auto gradients = _h->getGeomGradients();
  gradients += _g->getGeomGradients();
  gradients += _pcm->getGeomGradients();
  return gradients;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd HFPotentials<SCFMode>::getPointChargeGradients() {
  return _h->getPointChargeGradients();
}

template class HFPotentials<Options::SCF_MODES::RESTRICTED>;
template class HFPotentials<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
