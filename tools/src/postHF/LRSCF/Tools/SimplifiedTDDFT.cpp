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
  Timings::takeTime("LRSCF -    Simpl. TDDFT Prep.");

  unsigned nSub = _lrscf.size();

  _hfExchangeRatio.resize(nSub);

  _gammaK.resize(nSub);
  _gammaJ.resize(nSub);

  _Jij.resize(nSub);
  _Jai.resize(nSub);
  _Jab.resize(nSub);

  for (unsigned I = 0; I < nSub; ++I) {
    _gammaK[I].resize(nSub);
    _gammaJ[I].resize(nSub);

    auto sysI = _lrscf[I]->getSys();
    auto nAtomsI = sysI->getNAtoms();
    auto geomI = sysI->getGeometry()->getAtoms();

    for (unsigned J = 0; J < nSub; ++J) {
      auto sysJ = _lrscf[J]->getSys();
      auto nAtomsJ = sysJ->getNAtoms();
      auto geomJ = sysJ->getGeometry()->getAtoms();

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
        if (I == J) {
          _hfExchangeRatio[I] = 0.38;
        }
        alpha = 0.90;
        beta = 1.86;
      }
      else if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::LCBLYP) {
        if (I == J) {
          _hfExchangeRatio[I] = 0.53;
        }
        alpha = 4.50;
        beta = 8.00;
      }
      else if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::WB97) {
        if (I == J) {
          _hfExchangeRatio[I] = 0.61;
        }
        alpha = 4.41;
        beta = 8.00;
      }
      else if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::WB97X) {
        if (I == J) {
          _hfExchangeRatio[I] = 0.56;
        }
        alpha = 4.58;
        beta = 8.00;
      }
      else if (funcEnum == CompositeFunctionals::XCFUNCTIONALS::WB97X_D) {
        if (I == J) {
          _hfExchangeRatio[I] = 0.51;
        }
        alpha = 4.51;
        beta = 8.00;
      }
      else {
        if (I == J) {
          _hfExchangeRatio[I] = func.getHfExchangeRatio();
        }
        alpha = a1 + a2 * _hfExchangeRatio[I];
        beta = b1 + b2 * _hfExchangeRatio[I];
      }

      std::string funcString;
      Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(funcString, funcEnum);
      if (I == J && _lrscf[I]->getLRSCFSettings().grimme) {
        printf("    Intra-subsystem :  %2i\n", I + 1);
        printf("   - - - - - - - - - - - - - - - - - - -\n");
        printf("    ax         : %-10.4f\n", _hfExchangeRatio[I]);
        printf("   - - - - - - - - - - - - - - - - - - -\n");
        printf("    Functional : %-20s\n", funcString.c_str());
        printf("    alpha      : %-10.4f\n", alpha);
        printf("    beta       : %-10.4f\n", beta);
        printf("   -------------------------------------\n");
      }

      auto& gammaK = _gammaK[I][J];
      auto& gammaJ = _gammaJ[I][J];
      gammaK.resize(nAtomsI, nAtomsJ);
      gammaJ.resize(nAtomsI, nAtomsJ);
      for (unsigned iAtomA = 0; iAtomA < nAtomsI; ++iAtomA) {
        for (unsigned iAtomB = 0; iAtomB < nAtomsJ; ++iAtomB) {
          double avgHard = (geomI[iAtomA]->getChemicalHardness() + geomJ[iAtomB]->getChemicalHardness()) / 2.0;
          double dist = distance(*geomI[iAtomA], *geomJ[iAtomB]);

          gammaK(iAtomA, iAtomB) =
              (I != J) ? (1.0 / dist) : (std::pow(1.0 / (std::pow(dist, alpha) + std::pow(avgHard, -alpha)), 1.0 / alpha));
          gammaJ(iAtomA, iAtomB) =
              std::pow(1.0 / (std::pow(dist, beta) + std::pow(_hfExchangeRatio[I] * avgHard, -beta)), 1.0 / beta);
        }
      }
    }
  }

  // Calculate charges (resembling three-center integrals).
  for (unsigned I = 0; I < nSub; ++I) {
    Eigen::MatrixXd sqrtS = mSqrt_Sym(_lrscf[I]->getSys()->getOneElectronIntegralController()->getOverlapIntegrals());

    auto& C = _lrscf[I]->getCoefficients();
    auto indices = _lrscf[I]->getSys()->getAtomCenteredBasisController()->getBasisIndices();
    auto nAtoms = _lrscf[I]->getSys()->getNAtoms();

    auto no = _lrscf[I]->getNOccupied();
    auto nv = _lrscf[I]->getNVirtual();

    auto& Jij = _Jij[I];
    auto& Jai = _Jai[I];
    auto& Jab = _Jab[I];
    for_spin(no, nv, C, Jij, Jai, Jab) {
      Eigen::MatrixXd L = sqrtS * C_spin;
      Jij_spin.resize(no_spin * no_spin, nAtoms);
      Jai_spin.resize(nv_spin * no_spin, nAtoms);
      Jab_spin.resize(nv_spin * nv_spin, nAtoms);
      for (unsigned iAtom = 0; iAtom < nAtoms; ++iAtom) {
        unsigned start = indices[iAtom].first;
        unsigned stopAfter = indices[iAtom].second - indices[iAtom].first;
        Eigen::Ref<Eigen::MatrixXd> Li = L.block(start, 0, stopAfter, no_spin);
        Eigen::Ref<Eigen::MatrixXd> La = L.block(start, no_spin, stopAfter, nv_spin);
        Eigen::MatrixXd Lij = Li.transpose() * Li;
        Eigen::MatrixXd Lai = La.transpose() * Li;
        Eigen::MatrixXd Lab = La.transpose() * La;
        Jij_spin.col(iAtom) = Eigen::Map<Eigen::VectorXd>(Lij.data(), Lij.size());
        Jai_spin.col(iAtom) = Eigen::Map<Eigen::VectorXd>(Lai.data(), Lai.size());
        Jab_spin.col(iAtom) = Eigen::Map<Eigen::VectorXd>(Lab.data(), Lab.size());
      }
    };
  }
  Timings::timeTaken("LRSCF -    Simpl. TDDFT Prep.");
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>& SimplifiedTDDFT<SCFMode>::getJij() {
  return _Jij;
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>& SimplifiedTDDFT<SCFMode>::getJia() {
  return _Jai;
}

template<Options::SCF_MODES SCFMode>
std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>& SimplifiedTDDFT<SCFMode>::getJab() {
  return _Jab;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::vector<Eigen::MatrixXd>>& SimplifiedTDDFT<SCFMode>::getGammaK() {
  return _gammaK;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::vector<Eigen::MatrixXd>>& SimplifiedTDDFT<SCFMode>::getGammaJ() {
  return _gammaJ;
}

template<Options::SCF_MODES SCFMode>
std::vector<double> SimplifiedTDDFT<SCFMode>::getHFExchangeRatio() {
  return _hfExchangeRatio;
}

template class SimplifiedTDDFT<Options::SCF_MODES::RESTRICTED>;
template class SimplifiedTDDFT<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
