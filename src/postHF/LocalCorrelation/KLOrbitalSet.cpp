/**
 * @file KLOrbitalSet.cpp
 *
 * @author Moritz Bensberg
 * @date Dec 10, 2019
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
#include "postHF/LocalCorrelation/KLOrbitalSet.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalPair.h"                                      //Definition of an OrbitalPair.
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h" //Overlap matrices between domains.

namespace Serenity {

KLOrbitalSet::KLOrbitalSet(std::shared_ptr<OrbitalPair> ijPair, std::shared_ptr<OrbitalPair> klPair)
  : _ijPair(ijPair), _klPair(klPair) {
}

const Eigen::MatrixXd& KLOrbitalSet::getS_ij_kl() {
  if (!_s_ij_kl) {
    _s_ij_kl = _domainSController.lock()->getS(_ijPair.lock(), _klPair.lock());
  }
  return *_s_ij_kl;
}
const Eigen::MatrixXd& KLOrbitalSet::getS_ij_k() {
  if (!_s_ij_k) {
    _s_ij_k = _domainSController.lock()->getS(_ijPair.lock(), _klPair.lock()->singles_i);
  }
  return *_s_ij_k;
}
const Eigen::MatrixXd& KLOrbitalSet::getS_ij_l() {
  if (!_s_ij_l) {
    _s_ij_l = _domainSController.lock()->getS(_ijPair.lock(), _klPair.lock()->singles_j);
  }
  return *_s_ij_l;
}
void KLOrbitalSet::setOverlapMatrixController(std::shared_ptr<DomainOverlapMatrixController> domainSController) {
  _domainSController = domainSController;
}

void KLOrbitalSet::cleanUp() {
  ki_la.resize(0);
  kj_la.resize(0);
  li_ka.resize(0);
  lj_ka.resize(0);
}

} /* namespace Serenity */
