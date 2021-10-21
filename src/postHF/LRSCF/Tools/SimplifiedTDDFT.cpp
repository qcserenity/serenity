/**
 * @file SimplifiedTDDFT.cpp
 *
 * @date Oct 01, 2021
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
#include "postHF/LRSCF/Tools/SimplifiedTDDFT.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "dft/Functional.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/linearAlgebra/MatrixFunctions.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SimplifiedTDDFT<SCFMode>::SimplifiedTDDFT(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf)
  : _lrscf(lrscf) {
  this->setupSimplifiedTDDFT();
}

template<Options::SCF_MODES SCFMode>
void SimplifiedTDDFT<SCFMode>::setupSimplifiedTDDFT() {
  printBigCaption("Simplified TDDFT");

  unsigned nSub = _lrscf.size();
  for (unsigned I = 0; I < nSub; ++I) {
    auto sysI = _lrscf[I]->getSys();
    for (unsigned J = I; J < nSub; ++J) {
      auto sysJ = _lrscf[J]->getSys();

      // First things first
      auto nAtoms = sysI->getNAtoms();
      auto geom = sysI->getGeometry()->getAtoms();
      auto no = _lrscf[I]->getNOccupied();
      auto nv = _lrscf[I]->getNVirtual();

      double a1 = 1.42, a2 = 0.48;
      double b1 = 0.20, b2 = 1.83;
      double alpha, beta;

      CompositeFunctionals::XCFUNCTIONALS funcEnum;
      if (I == J) {
        if (_lrscf[I]->getSys()->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
          funcEnum = sysI->getSettings().dft.functional;
        }
        else {
          funcEnum = CompositeFunctionals::XCFUNCTIONALS::HF;
        }
      }
      else {
        funcEnum = _lrscf[I]->getLRSCFSettings().embedding.naddXCFunc;
      }

      // Source: sTDA github Oct 2021 (manual)
      Functional func = CompositeFunctionals::resolveFunctional(funcEnum);
      if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP) {
        _hfExchangeRatio = 0.38;
        alpha = 0.90;
        beta = 1.86;
      }
      else if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::LCBLYP) {
        _hfExchangeRatio = 0.53;
        alpha = 4.50;
        beta = 8.00;
      }
      else if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::WB97) {
        _hfExchangeRatio = 0.61;
        alpha = 4.41;
        beta = 8.00;
      }
      else if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::WB97X) {
        _hfExchangeRatio = 0.56;
        alpha = 4.58;
        beta = 8.00;
      }
      else if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::WB97X_D) {
        _hfExchangeRatio = 0.51;
        alpha = 4.51;
        beta = 8.00;
      }
      else {
        _hfExchangeRatio = func.getHfExchangeRatio();
        alpha = a1 + a2 * _hfExchangeRatio;
        beta = b1 + b2 * _hfExchangeRatio;
      }

      std::string funcString;
      Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(funcString, funcEnum);
      printf("    Functional : %-20s\n", funcString.c_str());
      printf("    ax         : %-10.4f\n", _hfExchangeRatio);
      printf("    alpha      : %-10.4f\n", alpha);
      printf("    beta       : %-10.4f\n\n", beta);

      Eigen::MatrixXd gammaK(nAtoms, nAtoms);
      Eigen::MatrixXd gammaJ(nAtoms, nAtoms);
      for (unsigned iAtomA = 0; iAtomA < nAtoms; ++iAtomA) {
        for (unsigned iAtomB = 0; iAtomB < nAtoms; ++iAtomB) {
          double hardA = geom[iAtomA]->getChemicalHardness();
          double hardB = geom[iAtomB]->getChemicalHardness();
          double avgHard = (hardA + hardB) / 2.0;
          double atomDist = distance(*geom[iAtomA], *geom[iAtomB]);

          gammaK(iAtomA, iAtomB) = std::pow(1.0 / (std::pow(atomDist, alpha) + std::pow(avgHard, -alpha)), 1.0 / alpha);
          gammaJ(iAtomA, iAtomB) =
              std::pow(1.0 / (std::pow(atomDist, beta) + std::pow(_hfExchangeRatio * avgHard, -beta)), 1.0 / beta);
        }
      }

      // Integrals will be calculated as (pq|rs) = <Q_pq|G|Q_rs>.
      // To save half the multiplications, factorize G = sqrt(G)sqrt(G) and
      // the integral above thus as (pq|rs) = <T_pq|T_rs> where T = Q * sqrt(G).
      Eigen::MatrixXd sqrtK = mSqrt_Sym(gammaK);
      Eigen::MatrixXd sqrtJ = mSqrt_Sym(gammaJ);
      Eigen::MatrixXd sqrtS = mSqrt_Sym(sysI->getOneElectronIntegralController()->getOverlapIntegrals());

      auto& C = _lrscf[I]->getCoefficients();
      auto indices = sysI->getAtomCenteredBasisController()->getBasisIndices();
      for_spin(no, nv, C, _Jij, _Jai, _Jab) {
        Eigen::MatrixXd L = sqrtS * C_spin;
        Eigen::MatrixXd Qij(no_spin * no_spin, nAtoms);
        Eigen::MatrixXd Qai(nv_spin * no_spin, nAtoms);
        Eigen::MatrixXd Qab(nv_spin * nv_spin, nAtoms);
        for (unsigned iAtom = 0; iAtom < nAtoms; ++iAtom) {
          unsigned start = indices[iAtom].first;
          unsigned stopAfter = indices[iAtom].second - indices[iAtom].first;
          Eigen::Ref<Eigen::MatrixXd> Li = L.block(start, 0, stopAfter, no_spin);
          Eigen::Ref<Eigen::MatrixXd> La = L.block(start, no_spin, stopAfter, nv_spin);
          Eigen::MatrixXd Lij = Li.transpose() * Li;
          Eigen::MatrixXd Lai = La.transpose() * Li;
          Eigen::MatrixXd Lab = La.transpose() * La;
          Qij.col(iAtom) = Eigen::Map<Eigen::VectorXd>(Lij.data(), Lij.size());
          Qai.col(iAtom) = Eigen::Map<Eigen::VectorXd>(Lai.data(), Lai.size());
          Qab.col(iAtom) = Eigen::Map<Eigen::VectorXd>(Lab.data(), Lab.size());
        }
        _Jai_spin = Qai * sqrtK;
        _Jij_spin = Qij * sqrtJ;
        _Jab_spin = Qab * sqrtJ;
      };
    }
  }
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::MatrixXd>& SimplifiedTDDFT<SCFMode>::getJij() {
  return _Jij;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::MatrixXd>& SimplifiedTDDFT<SCFMode>::getJai() {
  return _Jai;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::MatrixXd>& SimplifiedTDDFT<SCFMode>::getJab() {
  return _Jab;
}

template<Options::SCF_MODES SCFMode>
double SimplifiedTDDFT<SCFMode>::getHFExchangeRatio() {
  return _hfExchangeRatio;
}

template class SimplifiedTDDFT<Options::SCF_MODES::RESTRICTED>;
template class SimplifiedTDDFT<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
