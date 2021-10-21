/**
 * @file RIMP2.cpp
 *
 * @date Aug 14, 2017
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
#include "postHF/MPn/RIMP2.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
RIMP2<SCFMode>::RIMP2(std::shared_ptr<SystemController> systemController, const double ssScaling, const double osScaling)
  : _systemController(systemController), _Jia(nullptr), _sss(ssScaling), _oss(osScaling) {
}

template<Options::SCF_MODES SCFMode>
double RIMP2<SCFMode>::calculateCorrection() {
  printSmallCaption("RI-MP2 Calculation");

  RIIntegrals<SCFMode> ri(_systemController);
  _Jia = ri.getJiaPtr();
  return this->calculateEnergy();
}

template<>
double RIMP2<UNRESTRICTED>::calculateEnergy() {
  auto no = _systemController->getNOccupiedOrbitals<UNRESTRICTED>();
  auto nv = _systemController->getNVirtualOrbitals<UNRESTRICTED>();

  long noa = no.alpha, nob = no.beta;
  long nva = nv.alpha, nvb = nv.beta;

  auto e = _systemController->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  double ssEnergy = 0.0, osEnergy = 0.0;

  // alpha
  Eigen::MatrixXd eVirta = Eigen::MatrixXd::Zero(nva, nva);
  eVirta.colwise() += e.alpha.tail(nva);
  eVirta.rowwise() += e.alpha.tail(nva).transpose();
  // beta
  Eigen::MatrixXd eVirtb = Eigen::MatrixXd::Zero(nvb, nvb);
  eVirtb.colwise() += e.beta.tail(nvb);
  eVirtb.rowwise() += e.beta.tail(nvb).transpose();
  // mixed spin
  Eigen::MatrixXd eVirtab = Eigen::MatrixXd::Zero(nva, nvb);
  eVirtab.colwise() += e.alpha.tail(nva);
  eVirtab.rowwise() += e.beta.tail(nvb).transpose();

  Eigen::setNbThreads(1);
#pragma omp parallel for reduction(+ : ssEnergy) schedule(dynamic)
  for (long ij = 0; ij < noa * noa; ++ij) {
    long i = ij / noa, j = ij % noa;
    if (i > j) {
      continue;
    }
    Eigen::MatrixXd iajb = _Jia->alpha.middleRows(i * nva, nva) * _Jia->alpha.middleRows(j * nva, nva).transpose();
    Eigen::MatrixXd amp = iajb.cwiseQuotient(Eigen::MatrixXd::Constant(nva, nva, e.alpha(i) + e.alpha(j)) - eVirta);
    ssEnergy += (i == j ? 1.0 : 2.0) * (amp - amp.transpose()).cwiseProduct(iajb).sum();
  }
#pragma omp parallel for reduction(+ : ssEnergy) schedule(dynamic)
  for (long ij = 0; ij < nob * nob; ++ij) {
    long i = ij / nob, j = ij % nob;
    if (i > j) {
      continue;
    }
    Eigen::MatrixXd iajb = _Jia->beta.middleRows(i * nvb, nvb) * _Jia->beta.middleRows(j * nvb, nvb).transpose();
    Eigen::MatrixXd amp = iajb.cwiseQuotient(Eigen::MatrixXd::Constant(nvb, nvb, e.beta(i) + e.beta(j)) - eVirtb);
    ssEnergy += (i == j ? 1.0 : 2.0) * (amp - amp.transpose()).cwiseProduct(iajb).sum();
  }
#pragma omp parallel for reduction(+ : osEnergy) schedule(dynamic)
  for (long ij = 0; ij < noa * nob; ++ij) {
    long i_a = ij / nob, j_b = ij % nob;
    Eigen::MatrixXd iajb = _Jia->alpha.middleRows(i_a * nva, nva) * _Jia->beta.middleRows(j_b * nvb, nvb).transpose();
    Eigen::MatrixXd amp = iajb.cwiseQuotient(Eigen::MatrixXd::Constant(nva, nvb, e.alpha(i_a) + e.beta(j_b)) - eVirtab);
    osEnergy += amp.cwiseProduct(iajb).sum();
  }
  Eigen::setNbThreads(0);

  return 0.5 * _sss * ssEnergy + _oss * osEnergy;
}

template<>
double RIMP2<RESTRICTED>::calculateEnergy() {
  long no = _systemController->getNOccupiedOrbitals<RESTRICTED>();
  long nv = _systemController->getNVirtualOrbitals<RESTRICTED>();
  auto e = _systemController->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  double energy = 0.0;
  Eigen::MatrixXd eVirt = Eigen::MatrixXd::Zero(nv, nv);
  eVirt.colwise() += e.tail(nv);
  eVirt.rowwise() += e.tail(nv).transpose();

  Eigen::setNbThreads(1);
#pragma omp parallel for reduction(+ : energy) schedule(dynamic)
  for (long ij = 0; ij < no * no; ++ij) {
    long i = ij / no, j = ij % no;
    if (i > j) {
      continue;
    }
    Eigen::MatrixXd iajb = _Jia->middleRows(i * nv, nv) * _Jia->middleRows(j * nv, nv).transpose();
    Eigen::MatrixXd amp = iajb.cwiseQuotient(Eigen::MatrixXd::Constant(nv, nv, e(i) + e(j)) - eVirt);
    energy += (i == j ? 1.0 : 2.0) * ((_sss + _oss) * amp - _sss * amp.transpose()).cwiseProduct(iajb).sum();
  }
  Eigen::setNbThreads(0);

  return energy;
}

template class RIMP2<Options::SCF_MODES::RESTRICTED>;
template class RIMP2<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
