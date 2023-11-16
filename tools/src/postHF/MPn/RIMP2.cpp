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
DensityMatrix<UNRESTRICTED> RIMP2<UNRESTRICTED>::calculateDensityCorrection() {
  if (_Jia == 0) {
    RIIntegrals<UNRESTRICTED> ri(_systemController);
    _Jia = ri.getJiaPtr();
  }

  takeTime("Unrelaxed MP2 density");

  auto no = _systemController->getNOccupiedOrbitals<UNRESTRICTED>();
  auto nv = _systemController->getNVirtualOrbitals<UNRESTRICTED>();
  unsigned noa = no.alpha, nob = no.beta;
  unsigned nva = nv.alpha, nvb = nv.beta;
  auto e = _systemController->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  // alpha density matrix

  // alpha eigenvalues
  Eigen::MatrixXd eVirta = Eigen::MatrixXd::Zero(nva, nva);
  eVirta.colwise() += e.alpha.tail(nva);
  eVirta.rowwise() += e.alpha.tail(nva).transpose();
  Eigen::MatrixXd eOcca = Eigen::MatrixXd::Zero(noa, noa);
  eOcca.colwise() += e.alpha.head(noa);
  eOcca.rowwise() += e.alpha.head(noa).transpose();
  // mixed spin eigenvalues
  Eigen::MatrixXd eVirtab = Eigen::MatrixXd::Zero(nva, nvb);
  eVirtab.colwise() += e.alpha.tail(nva);
  eVirtab.rowwise() += e.beta.tail(nvb).transpose();
  Eigen::MatrixXd eOccab = Eigen::MatrixXd::Zero(noa, nob);
  eOccab.colwise() += e.alpha.head(noa);
  eOccab.rowwise() += e.beta.head(nob).transpose();

  // occupied block
  Eigen::MatrixXd dOcca = Eigen::MatrixXd::Zero(noa, noa);
  for (unsigned i = 0; i < noa; ++i) {
    for (unsigned j = i; j < noa; ++j) {
      for (unsigned k = 0; k < noa; ++k) {
        Eigen::MatrixXd iakb = _Jia->alpha.middleRows(i * nva, nva) * _Jia->alpha.middleRows(k * nva, nva).transpose();
        Eigen::MatrixXd jakb = _Jia->alpha.middleRows(j * nva, nva) * _Jia->alpha.middleRows(k * nva, nva).transpose();
        Eigen::MatrixXd ampik = iakb.cwiseQuotient(Eigen::MatrixXd::Constant(nva, nva, e.alpha(i) + e.alpha(k)) - eVirta);
        Eigen::MatrixXd ampjk = jakb.cwiseQuotient(Eigen::MatrixXd::Constant(nva, nva, e.alpha(j) + e.alpha(k)) - eVirta);
        dOcca(i, j) += (-0.5) * (ampik - ampik.transpose()).cwiseProduct(ampjk - ampjk.transpose()).sum();
      }
      for (unsigned k = 0; k < nob; ++k) {
        Eigen::MatrixXd iakb = _Jia->alpha.middleRows(i * nva, nva) * _Jia->beta.middleRows(k * nvb, nvb).transpose();
        Eigen::MatrixXd jakb = _Jia->alpha.middleRows(j * nva, nva) * _Jia->beta.middleRows(k * nvb, nvb).transpose();
        Eigen::MatrixXd ampik = iakb.cwiseQuotient(Eigen::MatrixXd::Constant(nva, nvb, e.alpha(i) + e.beta(k)) - eVirtab);
        Eigen::MatrixXd ampjk = jakb.cwiseQuotient(Eigen::MatrixXd::Constant(nva, nvb, e.alpha(j) + e.beta(k)) - eVirtab);
        dOcca(i, j) += (-1.0) * ampik.cwiseProduct(ampjk).sum();
      }
      dOcca(j, i) = dOcca(i, j);
    }
  }

  // virtual block
  Eigen::MatrixXd dVirta = Eigen::MatrixXd::Zero(nva, nva);
  Eigen::MatrixXd iaQ = Eigen::MatrixXd::Zero(noa, _Jia->alpha.cols());
  Eigen::MatrixXd ibQ = Eigen::MatrixXd::Zero(noa, _Jia->alpha.cols());
  Eigen::MatrixXd jcQ = Eigen::MatrixXd::Zero(noa, _Jia->alpha.cols());
  Eigen::MatrixXd jcQ_beta = Eigen::MatrixXd::Zero(nob, _Jia->beta.cols());
  for (unsigned a = 0; a < nva; ++a) {
    for (unsigned b = a; b < nva; ++b) {
      for (unsigned c = 0; c < nva; ++c) {
        for (unsigned i = 0; i < noa; ++i) {
          iaQ.row(i) = _Jia->alpha.row(a + i * nva);
          ibQ.row(i) = _Jia->alpha.row(b + i * nva);
          jcQ.row(i) = _Jia->alpha.row(c + i * nva);
        }
        Eigen::MatrixXd iajc = iaQ * jcQ.transpose();
        Eigen::MatrixXd ibjc = ibQ * jcQ.transpose();
        Eigen::MatrixXd ampac =
            iajc.cwiseQuotient(eOcca - Eigen::MatrixXd::Constant(noa, noa, e.alpha(noa + a) + e.alpha(noa + c)));
        Eigen::MatrixXd ampbc =
            ibjc.cwiseQuotient(eOcca - Eigen::MatrixXd::Constant(noa, noa, e.alpha(noa + b) + e.alpha(noa + c)));
        dVirta(a, b) += 0.5 * (ampac - ampac.transpose()).cwiseProduct(ampbc - ampbc.transpose()).sum();
      }
      for (unsigned c = 0; c < nvb; ++c) {
        for (unsigned i = 0; i < noa; ++i) {
          iaQ.row(i) = _Jia->alpha.row(a + i * nva);
          ibQ.row(i) = _Jia->alpha.row(b + i * nva);
        }
        for (unsigned i = 0; i < nob; ++i) {
          jcQ_beta.row(i) = _Jia->beta.row(c + i * nvb);
        }
        Eigen::MatrixXd iajc = iaQ * jcQ_beta.transpose();
        Eigen::MatrixXd ibjc = ibQ * jcQ_beta.transpose();
        Eigen::MatrixXd ampac =
            iajc.cwiseQuotient(eOccab - Eigen::MatrixXd::Constant(noa, nob, e.alpha(noa + a) + e.beta(nob + c)));
        Eigen::MatrixXd ampbc =
            ibjc.cwiseQuotient(eOccab - Eigen::MatrixXd::Constant(noa, nob, e.alpha(noa + b) + e.beta(nob + c)));
        dVirta(a, b) += ampac.cwiseProduct(ampbc).sum();
      }
      dVirta(b, a) = dVirta(a, b);
    }
  }

  // beta density matrix

  // beta eigenvalues
  Eigen::MatrixXd eVirtb = Eigen::MatrixXd::Zero(nvb, nvb);
  eVirtb.colwise() += e.beta.tail(nvb);
  eVirtb.rowwise() += e.beta.tail(nvb).transpose();
  Eigen::MatrixXd eOccb = Eigen::MatrixXd::Zero(nob, nob);
  eOccb.colwise() += e.beta.head(nob);
  eOccb.rowwise() += e.beta.head(nob).transpose();

  // occupied block
  Eigen::MatrixXd dOccb = Eigen::MatrixXd::Zero(nob, nob);
  for (unsigned i = 0; i < nob; ++i) {
    for (unsigned j = i; j < nob; ++j) {
      for (unsigned k = 0; k < nob; ++k) {
        Eigen::MatrixXd iakb = _Jia->beta.middleRows(i * nvb, nvb) * _Jia->beta.middleRows(k * nvb, nvb).transpose();
        Eigen::MatrixXd jakb = _Jia->beta.middleRows(j * nvb, nvb) * _Jia->beta.middleRows(k * nvb, nvb).transpose();
        Eigen::MatrixXd ampik = iakb.cwiseQuotient(Eigen::MatrixXd::Constant(nvb, nvb, e.beta(i) + e.beta(k)) - eVirtb);
        Eigen::MatrixXd ampjk = jakb.cwiseQuotient(Eigen::MatrixXd::Constant(nvb, nvb, e.beta(j) + e.beta(k)) - eVirtb);
        dOccb(i, j) += (-0.5) * (ampik - ampik.transpose()).cwiseProduct(ampjk - ampjk.transpose()).sum();
      }
      for (unsigned k = 0; k < noa; ++k) {
        Eigen::MatrixXd iakb = _Jia->beta.middleRows(i * nvb, nvb) * _Jia->alpha.middleRows(k * nva, nva).transpose();
        Eigen::MatrixXd jakb = _Jia->beta.middleRows(j * nvb, nvb) * _Jia->alpha.middleRows(k * nva, nva).transpose();
        Eigen::MatrixXd ampik =
            iakb.cwiseQuotient(Eigen::MatrixXd::Constant(nvb, nva, e.beta(i) + e.alpha(k)) - eVirtab.transpose());
        Eigen::MatrixXd ampjk =
            jakb.cwiseQuotient(Eigen::MatrixXd::Constant(nvb, nva, e.beta(j) + e.alpha(k)) - eVirtab.transpose());
        dOccb(i, j) += (-1.0) * ampik.cwiseProduct(ampjk).sum();
      }
      dOccb(j, i) = dOccb(i, j);
    }
  }

  // virtual block
  Eigen::MatrixXd dVirtb = Eigen::MatrixXd::Zero(nvb, nvb);
  Eigen::MatrixXd iaQb = Eigen::MatrixXd::Zero(nob, _Jia->beta.cols());
  Eigen::MatrixXd ibQb = Eigen::MatrixXd::Zero(nob, _Jia->beta.cols());
  Eigen::MatrixXd jcQb = Eigen::MatrixXd::Zero(nob, _Jia->beta.cols());
  Eigen::MatrixXd jcQ_alpha = Eigen::MatrixXd::Zero(noa, _Jia->alpha.cols());
  for (unsigned a = 0; a < nvb; ++a) {
    for (unsigned b = a; b < nvb; ++b) {
      for (unsigned c = 0; c < nvb; ++c) {
        for (unsigned i = 0; i < nob; ++i) {
          iaQb.row(i) = _Jia->beta.row(a + i * nvb);
          ibQb.row(i) = _Jia->beta.row(b + i * nvb);
          jcQb.row(i) = _Jia->beta.row(c + i * nvb);
        }
        Eigen::MatrixXd iajc = iaQb * jcQb.transpose();
        Eigen::MatrixXd ibjc = ibQb * jcQb.transpose();
        Eigen::MatrixXd ampac =
            iajc.cwiseQuotient(eOccb - Eigen::MatrixXd::Constant(nob, nob, e.beta(nob + a) + e.beta(nob + c)));
        Eigen::MatrixXd ampbc =
            ibjc.cwiseQuotient(eOccb - Eigen::MatrixXd::Constant(nob, nob, e.beta(nob + b) + e.beta(nob + c)));
        dVirtb(a, b) += 0.5 * (ampac - ampac.transpose()).cwiseProduct(ampbc - ampbc.transpose()).sum();
      }
      for (unsigned c = 0; c < nva; ++c) {
        for (unsigned i = 0; i < nob; ++i) {
          iaQb.row(i) = _Jia->beta.row(a + i * nvb);
          ibQb.row(i) = _Jia->beta.row(b + i * nvb);
        }
        for (unsigned i = 0; i < noa; ++i) {
          jcQ_alpha.row(i) = _Jia->alpha.row(c + i * nva);
        }
        Eigen::MatrixXd iajc = iaQb * jcQ_alpha.transpose();
        Eigen::MatrixXd ibjc = ibQb * jcQ_alpha.transpose();
        Eigen::MatrixXd ampac =
            iajc.cwiseQuotient(eOccab.transpose() - Eigen::MatrixXd::Constant(nob, noa, e.beta(nob + a) + e.alpha(noa + c)));
        Eigen::MatrixXd ampbc =
            ibjc.cwiseQuotient(eOccab.transpose() - Eigen::MatrixXd::Constant(nob, noa, e.beta(nob + b) + e.alpha(noa + c)));
        dVirtb(a, b) += ampac.cwiseProduct(ampbc).sum();
      }
      dVirtb(b, a) = dVirtb(a, b);
    }
  }

  // transformation from MO to AO

  DensityMatrix<UNRESTRICTED> densityCorrection(_systemController->getBasisController());
  auto coefficients = _systemController->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();

  // alpha
  Eigen::MatrixXd densityCorrectionMOAOa = Eigen::MatrixXd::Zero(noa + nva, noa + nva);
  densityCorrectionMOAOa.topRows(noa) = dOcca * coefficients.alpha.transpose().topRows(noa);
  densityCorrectionMOAOa.bottomRows(nva) = dVirta * coefficients.alpha.transpose().bottomRows(nva);
  densityCorrection.alpha += coefficients.alpha * densityCorrectionMOAOa;
  // beta
  Eigen::MatrixXd densityCorrectionMOAOb = Eigen::MatrixXd::Zero(nob + nvb, nob + nvb);
  densityCorrectionMOAOb.topRows(nob) = dOccb * coefficients.beta.transpose().topRows(nob);
  densityCorrectionMOAOb.bottomRows(nvb) = dVirtb * coefficients.beta.transpose().bottomRows(nvb);
  densityCorrection.beta += coefficients.beta * densityCorrectionMOAOb;

  timeTaken(0, "Unrelaxed MP2 density");

  return densityCorrection;
}

template<>
DensityMatrix<RESTRICTED> RIMP2<RESTRICTED>::calculateDensityCorrection() {
  if (_Jia == 0) {
    RIIntegrals<RESTRICTED> ri(_systemController);
    _Jia = ri.getJiaPtr();
  }

  takeTime("Unrelaxed MP2 density");

  unsigned no = _systemController->getNOccupiedOrbitals<RESTRICTED>();
  unsigned nv = _systemController->getNVirtualOrbitals<RESTRICTED>();
  auto e = _systemController->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  Eigen::MatrixXd eVirt = Eigen::MatrixXd::Zero(nv, nv);
  eVirt.colwise() += e.tail(nv);
  eVirt.rowwise() += e.tail(nv).transpose();
  Eigen::MatrixXd eOcc = Eigen::MatrixXd::Zero(no, no);
  eOcc.colwise() += e.head(no);
  eOcc.rowwise() += e.head(no).transpose();

  // occupied block
  Eigen::MatrixXd dOcc = Eigen::MatrixXd::Zero(no, no);
  for (unsigned i = 0; i < no; ++i) {
    for (unsigned j = i; j < no; ++j) {
      for (unsigned k = 0; k < no; ++k) {
        Eigen::MatrixXd iakb = _Jia->middleRows(i * nv, nv) * _Jia->middleRows(k * nv, nv).transpose();
        Eigen::MatrixXd jakb = _Jia->middleRows(j * nv, nv) * _Jia->middleRows(k * nv, nv).transpose();
        Eigen::MatrixXd ampik = iakb.cwiseQuotient(Eigen::MatrixXd::Constant(nv, nv, e(i) + e(k)) - eVirt);
        Eigen::MatrixXd ampjk = jakb.cwiseQuotient(Eigen::MatrixXd::Constant(nv, nv, e(j) + e(k)) - eVirt);
        dOcc(i, j) +=
            (-1.0) *
            ((ampik - ampik.transpose()).cwiseProduct(ampjk - ampjk.transpose()) + 2 * ampik.cwiseProduct(ampjk)).sum();
      }
      dOcc(j, i) = dOcc(i, j);
    }
  }

  // virtual block
  Eigen::MatrixXd dVirt = Eigen::MatrixXd::Zero(nv, nv);
  Eigen::MatrixXd iaQ = Eigen::MatrixXd::Zero(no, _Jia->cols());
  Eigen::MatrixXd ibQ = Eigen::MatrixXd::Zero(no, _Jia->cols());
  Eigen::MatrixXd jcQ = Eigen::MatrixXd::Zero(no, _Jia->cols());
  for (unsigned a = 0; a < nv; ++a) {
    for (unsigned b = a; b < nv; ++b) {
      for (unsigned c = 0; c < nv; ++c) {
        for (unsigned i = 0; i < no; ++i) {
          iaQ.row(i) = _Jia->row(a + i * nv);
          ibQ.row(i) = _Jia->row(b + i * nv);
          jcQ.row(i) = _Jia->row(c + i * nv);
        }
        Eigen::MatrixXd iajc = iaQ * jcQ.transpose();
        Eigen::MatrixXd ibjc = ibQ * jcQ.transpose();
        Eigen::MatrixXd ampac = iajc.cwiseQuotient(eOcc - Eigen::MatrixXd::Constant(no, no, e(no + a) + e(no + c)));
        Eigen::MatrixXd ampbc = ibjc.cwiseQuotient(eOcc - Eigen::MatrixXd::Constant(no, no, e(no + b) + e(no + c)));
        dVirt(a, b) +=
            ((ampac - ampac.transpose()).cwiseProduct(ampbc - ampbc.transpose()) + 2 * ampac.cwiseProduct(ampbc)).sum();
      }
      dVirt(b, a) = dVirt(a, b);
    }
  }

  // transformation from MO to AO
  DensityMatrix<RESTRICTED> densityCorrection(_systemController->getBasisController());
  auto coefficients = _systemController->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  Eigen::MatrixXd densityCorrectionMOAO = Eigen::MatrixXd::Zero(no + nv, no + nv);
  densityCorrectionMOAO.topRows(no) = dOcc * coefficients.transpose().topRows(no);
  densityCorrectionMOAO.bottomRows(nv) = dVirt * coefficients.transpose().bottomRows(nv);
  densityCorrection += coefficients * densityCorrectionMOAO;

  timeTaken(0, "Unrelaxed MP2 density");

  return densityCorrection;
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
