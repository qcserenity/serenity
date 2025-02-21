/**
 * @file GrimmeSigmavector.cpp
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
#include "postHF/LRSCF/Sigmavectors/GrimmeSigmavector.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Tools/SimplifiedTDDFT.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GrimmeSigmavector<SCFMode>::GrimmeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                              std::vector<Eigen::MatrixXd> b, std::vector<int> pm,
                                              std::shared_ptr<SimplifiedTDDFT<SCFMode>> simplifiedTDDFT)
  : Sigmavector<SCFMode>(lrscf, b), _pm(pm), _simplifiedTDDFT(simplifiedTDDFT) {
}

template<Options::SCF_MODES SCFMode>
GrimmeSigmavector<SCFMode>::~GrimmeSigmavector() = default;

template<>
void GrimmeSigmavector<Options::SCF_MODES::RESTRICTED>::calcSigma() {
  Timings::takeTime("LRSCF -   Sigmavector: Grimme");

  unsigned nSub = _lrscf.size();

  long iaStart = 0;
  for (unsigned I = 0; I < nSub; ++I) {
    unsigned noI = _lrscf[I]->getNOccupied();
    unsigned nvI = _lrscf[I]->getNVirtual();

    auto& Jij_I = _simplifiedTDDFT->getJij()[I];
    auto& Jia_I = _simplifiedTDDFT->getJia()[I];
    auto& Jab_I = _simplifiedTDDFT->getJab()[I];

    for (unsigned J = 0; J < nSub; ++J) {
      unsigned noJ = _lrscf[J]->getNOccupied();
      unsigned nvJ = _lrscf[J]->getNVirtual();

      if (this->skipThisInteraction(I, J)) {
        continue;
      }

      // Only calculate intra-subsystem blocks if the grimme keyword is true.
      if (I == J && !this->_lrscf[I]->getLRSCFSettings().grimme) {
        continue;
      }

      double distance =
          this->_lrscf[I]->getSys()->getGeometry()->getMinimumDistance(*this->_lrscf[J]->getSys()->getGeometry());
      if (I != J && distance < this->_lrscf[I]->getLRSCFSettings().approxCoulomb[1]) {
        continue;
      }

      auto& Jia_J = _simplifiedTDDFT->getJia()[J];

      auto& gammaK_IJ = _simplifiedTDDFT->getGammaK()[I][J];
      auto& gammaJ_IJ = _simplifiedTDDFT->getGammaJ()[I][J];

      for (unsigned iSet = 0; iSet < _b.size(); ++iSet) {
        _sigma[iSet].middleRows(iaStart, nvI * noI).noalias() +=
            (2.0 + _pm[iSet] * 2.0) * Jia_I * gammaK_IJ * (Jia_J.transpose() * _b[iSet][J]).eval();
      }

      if (I == J) {
        for (unsigned iSet = 0; iSet < _b.size(); ++iSet) {
          Eigen::MatrixXd tJij = Jij_I * gammaJ_IJ;
          Eigen::MatrixXd tJai = _pm[iSet] * _simplifiedTDDFT->getHFExchangeRatio()[I] * Jia_I * gammaK_IJ;
          std::vector<Eigen::MatrixXd> sigmaParallel(omp_get_max_threads(),
                                                     Eigen::MatrixXd::Zero(nvI * noI, _sigma[iSet].cols()));
          Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic) collapse(2)
          for (unsigned iGuess = 0; iGuess < _b[iSet][I].cols(); ++iGuess) {
            for (unsigned iAtom = 0; iAtom < Jia_I.cols(); ++iAtom) {
              Eigen::Map<Eigen::MatrixXd> sigma(sigmaParallel[omp_get_thread_num()].col(iGuess).data(), nvI, noI);
              Eigen::Map<Eigen::MatrixXd> guess(_b[iSet][I].col(iGuess).data(), nvJ, noJ);

              Eigen::Map<Eigen::MatrixXd> ji(tJij.col(iAtom).data(), noI, noI);
              Eigen::Map<Eigen::MatrixXd> ab(Jab_I.col(iAtom).data(), nvI, nvI);
              sigma.noalias() -= (ab * guess) * ji;

              if (_pm[iSet]) {
                Eigen::Map<Eigen::MatrixXd> ia(tJai.col(iAtom).data(), nvI, noI);
                Eigen::Map<Eigen::MatrixXd> ai(Jia_I.col(iAtom).data(), nvI, noI);
                sigma.noalias() -= ai * (guess.transpose() * ia);
              }
            }
          }
          Eigen::setNbThreads(0);

          for (int iThread = 0; iThread < omp_get_max_threads(); ++iThread) {
            _sigma[iSet].middleRows(iaStart, nvI * noI) += sigmaParallel[iThread];
          }
        }
      }
    }
    iaStart += nvI * noI;
  }

  _hasBeenCalculated = true;
  Timings::timeTaken("LRSCF -   Sigmavector: Grimme");
}

template<>
void GrimmeSigmavector<Options::SCF_MODES::UNRESTRICTED>::calcSigma() {
  Timings::takeTime("LRSCF -   Sigmavector: Grimme");

  unsigned nSub = _lrscf.size();

  long iaStart = 0;
  for (unsigned I = 0; I < nSub; ++I) {
    auto noI = _lrscf[I]->getNOccupied();
    auto nvI = _lrscf[I]->getNVirtual();

    unsigned alphaI = nvI.alpha * noI.alpha;
    unsigned betaI = nvI.beta * noI.beta;

    auto& Jij_I = _simplifiedTDDFT->getJij()[I];
    auto& Jia_I = _simplifiedTDDFT->getJia()[I];
    auto& Jab_I = _simplifiedTDDFT->getJab()[I];

    for (unsigned J = 0; J < nSub; ++J) {
      auto noJ = _lrscf[J]->getNOccupied();
      auto nvJ = _lrscf[J]->getNVirtual();

      unsigned alphaJ = nvJ.alpha * noJ.alpha;
      unsigned betaJ = nvJ.beta * noJ.beta;

      if (this->skipThisInteraction(I, J)) {
        continue;
      }

      // Only calculate intra-subsystem blocks if the grimme keyword is true.
      if (I == J && !this->_lrscf[I]->getLRSCFSettings().grimme) {
        continue;
      }

      // Calculate Coulomb interactions to be calculated via Grimme integrals.
      double distance =
          this->_lrscf[I]->getSys()->getGeometry()->getMinimumDistance(*this->_lrscf[J]->getSys()->getGeometry());
      if (I != J && distance < this->_lrscf[I]->getLRSCFSettings().approxCoulomb[1]) {
        continue;
      }

      auto& Jia_J = _simplifiedTDDFT->getJia()[J];

      auto& gammaK = _simplifiedTDDFT->getGammaK()[I][J];
      auto& gammaJ = _simplifiedTDDFT->getGammaJ()[I][J];

      for (unsigned iSet = 0; iSet < _b.size(); ++iSet) {
        Eigen::MatrixXd X = (Jia_J.alpha.transpose() * _b[iSet][J].middleRows(0, alphaJ)).eval();
        X.noalias() += (Jia_J.beta.transpose() * _b[iSet][J].middleRows(alphaJ, betaJ)).eval();

        _sigma[iSet].middleRows(iaStart, alphaI).noalias() += (1.0 + _pm[iSet] * 1.0) * Jia_I.alpha * gammaK * X;
        _sigma[iSet].middleRows(iaStart + alphaI, betaI).noalias() += (1.0 + _pm[iSet] * 1.0) * Jia_I.beta * gammaK * X;
      }

      for (unsigned iSet = 0; iSet < _b.size(); ++iSet) {
        if (I == J) {
          Eigen::MatrixXd tJija = Jij_I.alpha * gammaJ;
          Eigen::MatrixXd tJijb = Jij_I.beta * gammaJ;
          Eigen::MatrixXd tJaia = _pm[iSet] * _simplifiedTDDFT->getHFExchangeRatio()[I] * Jia_I.alpha * gammaK;
          Eigen::MatrixXd tJaib = _pm[iSet] * _simplifiedTDDFT->getHFExchangeRatio()[I] * Jia_I.beta * gammaK;

          std::vector<Eigen::MatrixXd> sigmaParallel(
              omp_get_max_threads(), Eigen::MatrixXd::Zero(nvI.alpha * noI.alpha + nvI.beta * noI.beta, _sigma[iSet].cols()));
          Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic) collapse(2)
          for (unsigned iGuess = 0; iGuess < _b[iSet][I].cols(); ++iGuess) {
            for (unsigned iAtom = 0; iAtom < Jia_I.alpha.cols(); ++iAtom) {
              Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaParallel[omp_get_thread_num()].col(iGuess).data() + iaStart,
                                                  nvI.alpha, noI.alpha);
              Eigen::Map<Eigen::MatrixXd> guess_a(_b[iSet][I].col(iGuess).data(), nvI.alpha, noI.alpha);

              Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaParallel[omp_get_thread_num()].col(iGuess).data() + iaStart + alphaI,
                                                  nvI.beta, noI.beta);
              Eigen::Map<Eigen::MatrixXd> guess_b(_b[iSet][I].col(iGuess).data() + alphaI, nvI.beta, noI.beta);

              Eigen::Map<Eigen::MatrixXd> ji_a(tJija.col(iAtom).data(), noI.alpha, noI.alpha);
              Eigen::Map<Eigen::MatrixXd> ia_a(tJaia.col(iAtom).data(), nvI.alpha, noI.alpha);
              Eigen::Map<Eigen::MatrixXd> ai_a(Jia_I.alpha.col(iAtom).data(), nvI.alpha, noI.alpha);
              Eigen::Map<Eigen::MatrixXd> ab_a(Jab_I.alpha.col(iAtom).data(), nvI.alpha, nvI.alpha);

              Eigen::Map<Eigen::MatrixXd> ji_b(tJijb.col(iAtom).data(), noI.beta, noI.beta);
              Eigen::Map<Eigen::MatrixXd> ia_b(tJaib.col(iAtom).data(), nvI.beta, noI.beta);
              Eigen::Map<Eigen::MatrixXd> ai_b(Jia_I.beta.col(iAtom).data(), nvI.beta, noI.beta);
              Eigen::Map<Eigen::MatrixXd> ab_b(Jab_I.beta.col(iAtom).data(), nvI.beta, nvI.beta);

              sigma_a.noalias() -= (ab_a * guess_a) * ji_a;
              sigma_b.noalias() -= (ab_b * guess_b) * ji_b;
              if (_pm[iSet]) {
                sigma_a.noalias() -= ai_a * (guess_a.transpose() * ia_a);
                sigma_b.noalias() -= ai_b * (guess_b.transpose() * ia_b);
              }
            }
          }
          Eigen::setNbThreads(0);

          for (int iThread = 0; iThread < omp_get_max_threads(); ++iThread) {
            _sigma[iSet].middleRows(iaStart, nvI.alpha * noI.alpha + nvI.beta * noI.beta) += sigmaParallel[iThread];
          }
        }
      }
    }
    iaStart += alphaI + betaI;
  }

  _hasBeenCalculated = true;
  Timings::timeTaken("LRSCF -   Sigmavector: Grimme");
}

template class GrimmeSigmavector<Options::SCF_MODES::RESTRICTED>;
template class GrimmeSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
