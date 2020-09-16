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
#include "data/matrices/DensityMatrixController.h"
#include "integrals/RI_J_IntegralController.h"
#include "math/linearAlgebra/MatrixFunctions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
RIMP2<SCFMode>::RIMP2(std::shared_ptr<SystemController> systemController, const double ssScaling, const double osScaling)
  : _systemController(systemController), _ssScaling(ssScaling), _osScaling(osScaling) {
  assert(_systemController);
}

template<Options::SCF_MODES SCFMode>
double RIMP2<SCFMode>::calculateCorrection() {
  printSmallCaption("RI-MP2 Calculation");

  auto es = _systemController->getElectronicStructure<SCFMode>();
  auto basisController = _systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  auto auxBasisController = _systemController->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
  CoefficientMatrix<SCFMode> coeff = es->getMolecularOrbitals()->getCoefficients();

  if (_systemController->getSettings().basis.auxCLabel == "") {
    printf("\n  Auxiliary Basis Set         : %15s\n\n", (_systemController->getSettings().basis.label + "-RI-C").c_str());
  }
  else {
    printf("\n  Auxiliary Basis Set         : %15s\n\n", (_systemController->getSettings().basis.auxCLabel).c_str());
  }
  auto no = _systemController->getNOccupiedOrbitals<SCFMode>();
  auto nv = _systemController->getNVirtualOrbitals<SCFMode>();
  unsigned nb = basisController->getNBasisFunctions();
  unsigned nx = auxBasisController->getNBasisFunctions();

  double gbFactor = (double)sizeof(double) / (1024 * 1024 * 1024);

  unsigned iSpin = 0;
  for_spin(no, nv) {
    if (iSpin == 0) {
      printf("  Occupied Orbitals (alpha)   : %15i\n", no_spin);
      printf("  Virtual Orbitals  (alpha)   : %15i\n", nv_spin);
      printf("  Memory (ip|Q)     (alpha)   : %11.3f GiB\n\n", nx * no_spin * nb * gbFactor);
    }
    else {
      printf("  Occupied Orbitals (beta)    : %15i\n", no_spin);
      printf("  Virtual Orbitals  (beta)    : %15i\n", nv_spin);
      printf("  Memory (ip|Q)     (beta)    : %11.3f GiB\n\n", nx * no_spin * nb * gbFactor);
    }
    ++iSpin;
  };
  printf("  Basis Functions             : %15i\n", nb);
  printf("  Auxiliary Functions         : %15i\n\n", nx);

  auto riints = std::make_shared<RI_J_IntegralController>(basisController, auxBasisController);
  auto looper = std::make_shared<TwoElecThreeCenterIntLooper>(libint2::Operator::coulomb, 0, basisController,
                                                              auxBasisController, 1E-10);

  for_spin(coeff, _BiaQ, no, nv) {
    _BiaQ_spin = Eigen::MatrixXd::Zero(no_spin * nb, nx);
    auto biaptr = _BiaQ_spin.data();
    auto coeffptr = coeff_spin.data();
    Eigen::MatrixXd coeffVirtT = coeff_spin.rightCols(nv_spin).transpose();

    auto distribute = [&](const unsigned mu, const unsigned nu, const unsigned Q, const double integral, const unsigned) {
      unsigned long long offset = Q * no_spin * nb;
      for (unsigned i = 0; i < no_spin; ++i) {
        biaptr[offset + i * nb + mu] += coeffptr[i * nb + nu] * integral;
        if (mu != nu) {
          biaptr[offset + i * nb + nu] += coeffptr[i * nb + mu] * integral;
        }
      }
    };

    looper->loopNoDerivative(distribute);

    for (unsigned Q = 0; Q < nx; ++Q) {
      _BiaQ_spin.block(0, Q, no_spin * nv_spin, 1) =
          coeffVirtT * Eigen::Map<Eigen::MatrixXd>(biaptr + Q * nb * no_spin, nb, no_spin);
    }

    // throw away useless info in Bia
    _BiaQ_spin.conservativeResize(no_spin * nv_spin, nx);

    // avoid temporary object for everything at once
    for (unsigned i = 0; i < no_spin; ++i) {
      _BiaQ_spin.middleRows(i * nv_spin, nv_spin) *= riints->getInverseMSqrt();
    }
  };

  riints = nullptr;
  looper = nullptr;

  printf("  Done setting up (ia|Q) integrals.\n\n");
  return this->calculateEnergy();
}

template<>
double RIMP2<UNRESTRICTED>::calculateEnergy() {
  auto no = _systemController->getNOccupiedOrbitals<UNRESTRICTED>();
  auto nv = _systemController->getNVirtualOrbitals<UNRESTRICTED>();

  unsigned noa = no.alpha, nob = no.beta;
  unsigned nva = nv.alpha, nvb = nv.beta;

  auto e = _systemController->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  double ssEnergy = 0.0, osEnergy = 0.0;

  // alpha
  Eigen::MatrixXd eVirta = Eigen::MatrixXd::Zero(nva, nva);
  eVirta.colwise() += e.alpha.tail(nva);
  eVirta.rowwise() += e.alpha.tail(nva).transpose();

  for (unsigned i = 0; i < noa; ++i) {
    for (unsigned j = i; j < noa; ++j) {
      Eigen::MatrixXd iajb = _BiaQ.alpha.middleRows(i * nva, nva) * _BiaQ.alpha.middleRows(j * nva, nva).transpose();
      Eigen::MatrixXd amp = iajb.cwiseQuotient(Eigen::MatrixXd::Constant(nva, nva, e.alpha(i) + e.alpha(j)) - eVirta);
      ssEnergy += (i == j ? 1.0 : 2.0) * (amp - amp.transpose()).cwiseProduct(iajb).sum();
    }
  }

  // beta
  Eigen::MatrixXd eVirtb = Eigen::MatrixXd::Zero(nvb, nvb);
  eVirtb.colwise() += e.beta.tail(nvb);
  eVirtb.rowwise() += e.beta.tail(nvb).transpose();

  for (unsigned i = 0; i < nob; ++i) {
    for (unsigned j = i; j < nob; ++j) {
      Eigen::MatrixXd iajb = _BiaQ.beta.middleRows(i * nvb, nvb) * _BiaQ.beta.middleRows(j * nvb, nvb).transpose();
      Eigen::MatrixXd amp = iajb.cwiseQuotient(Eigen::MatrixXd::Constant(nvb, nvb, e.beta(i) + e.beta(j)) - eVirtb);
      ssEnergy += (i == j ? 1.0 : 2.0) * (amp - amp.transpose()).cwiseProduct(iajb).sum();
    }
  }

  // mixed spin
  Eigen::MatrixXd eVirtab = Eigen::MatrixXd::Zero(nva, nvb);
  eVirtab.colwise() += e.alpha.tail(nva);
  eVirtab.rowwise() += e.beta.tail(nvb).transpose();

  for (unsigned i = 0; i < noa; ++i) {
    for (unsigned j = 0; j < nob; ++j) {
      Eigen::MatrixXd iajb = _BiaQ.alpha.middleRows(i * nva, nva) * _BiaQ.beta.middleRows(j * nvb, nvb).transpose();
      Eigen::MatrixXd amp = iajb.cwiseQuotient(Eigen::MatrixXd::Constant(nva, nvb, e.alpha(i) + e.beta(j)) - eVirtab);
      osEnergy += amp.cwiseProduct(iajb).sum();
    }
  }

  return 0.5 * _ssScaling * ssEnergy + _osScaling * osEnergy;
}

template<>
double RIMP2<RESTRICTED>::calculateEnergy() {
  unsigned no = _systemController->getNOccupiedOrbitals<RESTRICTED>();
  unsigned nv = _systemController->getNVirtualOrbitals<RESTRICTED>();
  auto e = _systemController->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  double ssEnergy = 0.0, osEnergy = 0.0;
  Eigen::MatrixXd eVirt = Eigen::MatrixXd::Zero(nv, nv);
  eVirt.colwise() += e.tail(nv);
  eVirt.rowwise() += e.tail(nv).transpose();

  for (unsigned i = 0; i < no; ++i) {
    for (unsigned j = i; j < no; ++j) {
      Eigen::MatrixXd iajb = _BiaQ.middleRows(i * nv, nv) * _BiaQ.middleRows(j * nv, nv).transpose();
      Eigen::MatrixXd amp = iajb.cwiseQuotient(Eigen::MatrixXd::Constant(nv, nv, e(i) + e(j)) - eVirt);
      ssEnergy += (i == j ? 1.0 : 2.0) * (amp - amp.transpose()).cwiseProduct(iajb).sum();
      osEnergy += (i == j ? 1.0 : 2.0) * amp.cwiseProduct(iajb).sum();
    }
  }

  return _ssScaling * ssEnergy + _osScaling * osEnergy;
}

template class RIMP2<Options::SCF_MODES::RESTRICTED>;
template class RIMP2<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
