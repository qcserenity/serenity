/**
 * @file LTMP2.cpp
 *
 * @date Dec 30, 2022
 * @author Niklas Niemeyer
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
#include "postHF/MPn/LTMP2.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "io/FormattedOutputStream.h"
#include "math/LaplaceMinimaxWrapper.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
LTMP2<SCFMode>::LTMP2(std::shared_ptr<SystemController> systemController, const double osScaling, const double ltConv)
  : _systemController(systemController), _Jia(nullptr), _oss(osScaling), _ltConv(ltConv) {
}

template<Options::SCF_MODES SCFMode>
double LTMP2<SCFMode>::calculateCorrection() {
  printSmallCaption("LT-SOS-MP2 Calculation");
  OutputControl::nOut << " OS-Scaling: " << _oss << std::endl;

  RIIntegrals<SCFMode> ri(_systemController);
  _Jia = ri.getJiaPtr();
  return this->calculateEnergy();
}

template<>
double LTMP2<UNRESTRICTED>::calculateEnergy() {
  auto no = _systemController->getNOccupiedOrbitals<UNRESTRICTED>();
  auto nv = _systemController->getNVirtualOrbitals<UNRESTRICTED>();

  unsigned noa = no.alpha, nob = no.beta;
  unsigned nva = nv.alpha, nvb = nv.beta;

  double energy = 0.0;

  auto e = _systemController->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  double LUMO = std::min(e.alpha(noa), e.beta(nob));
  double HOMO = std::max(e.alpha(noa - 1), e.beta(nob - 1));
  double HUMO = std::max(e.alpha(e.alpha.size() - 1), e.beta(e.beta.size() - 1));
  double LOMO = std::min(e.alpha(0), e.beta(0));

  Eigen::VectorXd roots(0);
  Eigen::VectorXd weights(0);
  getMinimaxRoots(roots, weights, 2 * (LUMO - HOMO), 2 * (HUMO - LOMO), _ltConv);

  for (unsigned m = 0; m < roots.size(); ++m) {
    Eigen::VectorXd param_a(nva * noa);
    for (size_t i = 0, ia = 0; i < noa; ++i) {
      for (size_t a = noa; a < noa + nva; ++a, ++ia) {
        param_a(ia) = std::exp(-(e.alpha(a) - e.alpha(i)) * roots(m));
      }
    }
    Eigen::VectorXd param_b(nvb * nob);
    for (size_t j = 0, jb = 0; j < nob; ++j) {
      for (size_t b = nob; b < nob + nvb; ++b, ++jb) {
        param_b(jb) = std::exp(-(e.beta(b) - e.beta(j)) * roots(m));
      }
    }
    Eigen::MatrixXd QP_a = _Jia->alpha.transpose() * param_a.asDiagonal() * _Jia->alpha;
    Eigen::MatrixXd QP_b = _Jia->beta.transpose() * param_b.asDiagonal() * _Jia->beta;

    energy -= weights(m) * QP_a.cwiseProduct(QP_b).sum();
  }

  return _oss * energy;
}

template<>
double LTMP2<RESTRICTED>::calculateEnergy() {
  unsigned no = _systemController->getNOccupiedOrbitals<RESTRICTED>();
  unsigned nv = _systemController->getNVirtualOrbitals<RESTRICTED>();
  auto e = _systemController->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  double energy = 0.0;

  double LUMO = e(no);
  double HOMO = e(no - 1);
  double HUMO = e(e.size() - 1);
  double LOMO = e(0);

  Eigen::VectorXd roots(0);
  Eigen::VectorXd weights(0);
  getMinimaxRoots(roots, weights, 2 * (LUMO - HOMO), 2 * (HUMO - LOMO), _ltConv);

  for (unsigned m = 0; m < roots.size(); ++m) {
    Eigen::VectorXd param(nv * no);
    for (size_t i = 0, ia = 0; i < no; ++i) {
      for (size_t a = no; a < no + nv; ++a, ++ia) {
        param(ia) = std::exp(-(e(a) - e(i)) * roots(m));
      }
    }
    Eigen::MatrixXd QP = (*_Jia).transpose() * param.asDiagonal() * (*_Jia);
    energy -= weights(m) * QP.squaredNorm();
  }

  return _oss * energy;
}

template class LTMP2<Options::SCF_MODES::RESTRICTED>;
template class LTMP2<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
