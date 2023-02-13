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
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Tools/SimplifiedTDDFT.h"
#include "settings/Settings.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GrimmeSigmavector<SCFMode>::GrimmeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf)
  : _lrscf(lrscf), _simplifiedTDDFT(std::make_unique<SimplifiedTDDFT<SCFMode>>(lrscf)) {
}

template<Options::SCF_MODES SCFMode>
GrimmeSigmavector<SCFMode>::~GrimmeSigmavector() = default;

template<>
std::vector<Eigen::MatrixXd> GrimmeSigmavector<RESTRICTED>::getSigmavectors(std::vector<Eigen::MatrixXd>& guessVectors,
                                                                            std::vector<int> pm) {
  Timings::takeTime("LRSCF -   Sigmavector: Grimme");
  std::vector<Eigen::MatrixXd> sigmaVectors = guessVectors;
  for (auto& sigma : sigmaVectors) {
    sigma.setZero();
  }

  unsigned nSub = _lrscf.size();
  auto ax = _simplifiedTDDFT->getHFExchangeRatio();

  long iaStart = 0;
  for (unsigned I = 0; I < nSub; ++I) {
    unsigned noI = _lrscf[I]->getNOccupied();
    unsigned nvI = _lrscf[I]->getNVirtual();

    auto& Jij = _simplifiedTDDFT->getJij()[I];
    auto& Jai = _simplifiedTDDFT->getJai()[I];
    auto& Jab = _simplifiedTDDFT->getJab()[I];

    long jbStart = 0;
    for (unsigned J = 0; J < nSub; ++J) {
      unsigned noJ = _lrscf[J]->getNOccupied();
      unsigned nvJ = _lrscf[J]->getNVirtual();

      // Check if this block is to be skipped.
      bool exploitSymmetry = (_lrscf.size() > 1 && !_lrscf[I]->getLRSCFSettings().fullFDEc);
      exploitSymmetry = exploitSymmetry && _lrscf[I]->getLRSCFSettings().partialResponseConstruction;
      if (exploitSymmetry) {
        if (_lrscf[I]->getLRSCFSettings().method == Options::LR_METHOD::TDA) {
          if (I == J || I > J) {
            jbStart += nvJ * noJ;
            continue;
          }
        }
        else if (_lrscf[I]->getLRSCFSettings().method == Options::LR_METHOD::TDDFT) {
          if (I > J) {
            jbStart += nvJ * noJ;
            continue;
          }
        }
      }

      auto& Jia = _simplifiedTDDFT->getJai()[J];

      auto& gammaK = _simplifiedTDDFT->getGammaK()[I][J];
      auto& gammaJ = _simplifiedTDDFT->getGammaJ()[I][J];

      for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
        sigmaVectors[iSet].middleRows(iaStart, nvI * noI).noalias() +=
            (2.0 + pm[iSet] * 2.0) * Jai * gammaK *
            (Jia.transpose() * guessVectors[iSet].middleRows(jbStart, nvJ * noJ)).eval();
        if (I == J) {
          Eigen::MatrixXd tJij = Jij * gammaJ;
          Eigen::MatrixXd tJai = Jai * gammaK;
          for (unsigned iGuess = 0; iGuess < guessVectors[iSet].cols(); ++iGuess) {
            Eigen::Map<Eigen::MatrixXd> sigma(sigmaVectors[iSet].col(iGuess).data() + iaStart, nvI, noI);
            Eigen::Map<Eigen::MatrixXd> guess(guessVectors[iSet].col(iGuess).data() + jbStart, nvJ, noJ);
            for (unsigned iAtom = 0; iAtom < Jai.cols(); ++iAtom) {
              Eigen::Map<Eigen::MatrixXd> ji(tJij.col(iAtom).data(), noI, noI);
              Eigen::Map<Eigen::MatrixXd> ia(tJai.col(iAtom).data(), nvI, noI);
              Eigen::Map<Eigen::MatrixXd> ai(Jai.col(iAtom).data(), nvI, noI);
              Eigen::Map<Eigen::MatrixXd> ab(Jab.col(iAtom).data(), nvI, nvI);
              sigma.noalias() -= ab * guess * ji;
              if (pm[iSet]) {
                sigma.noalias() -= pm[iSet] * ax[I] * ai * guess.transpose() * ia;
              }
            }
          }
        }
      }
      jbStart += nvJ * noJ;
    }
    iaStart += nvI * noI;
  }

  Timings::timeTaken("LRSCF -   Sigmavector: Grimme");
  return sigmaVectors;
} /* this->getSigmavectors() restricted */

template<>
std::vector<Eigen::MatrixXd> GrimmeSigmavector<UNRESTRICTED>::getSigmavectors(std::vector<Eigen::MatrixXd>& guessVectors,
                                                                              std::vector<int> pm) {
  Timings::takeTime("LRSCF -   Sigmavector: Grimme");
  std::vector<Eigen::MatrixXd> sigmaVectors = guessVectors;
  for (auto& sigma : sigmaVectors) {
    sigma.setZero();
  }

  unsigned nSub = _lrscf.size();
  auto ax = _simplifiedTDDFT->getHFExchangeRatio();

  long iaStart = 0;
  for (unsigned I = 0; I < nSub; ++I) {
    auto noI = _lrscf[I]->getNOccupied();
    auto nvI = _lrscf[I]->getNVirtual();

    unsigned alphaI = nvI.alpha * noI.alpha;
    unsigned betaI = nvI.beta * noI.beta;

    auto& Jij = _simplifiedTDDFT->getJij()[I];
    auto& Jai = _simplifiedTDDFT->getJai()[I];
    auto& Jab = _simplifiedTDDFT->getJab()[I];

    long jbStart = 0;
    for (unsigned J = 0; J < nSub; ++J) {
      auto noJ = _lrscf[J]->getNOccupied();
      auto nvJ = _lrscf[J]->getNVirtual();

      unsigned alphaJ = nvJ.alpha * noJ.alpha;
      unsigned betaJ = nvJ.beta * noJ.beta;

      // Check if this block is to be skipped.
      bool exploitSymmetry = (_lrscf.size() > 1 && !_lrscf[I]->getLRSCFSettings().fullFDEc);
      exploitSymmetry = exploitSymmetry && _lrscf[I]->getLRSCFSettings().partialResponseConstruction;
      if (exploitSymmetry) {
        if (_lrscf[I]->getLRSCFSettings().method == Options::LR_METHOD::TDA) {
          if (I == J || I > J) {
            jbStart += alphaJ + betaJ;
            continue;
          }
        }
        else if (_lrscf[I]->getLRSCFSettings().method == Options::LR_METHOD::TDDFT) {
          if (I > J) {
            jbStart += alphaJ + betaJ;
            continue;
          }
        }
      }

      auto& Jia = _simplifiedTDDFT->getJai()[J];

      auto& gammaK = _simplifiedTDDFT->getGammaK()[I][J];
      auto& gammaJ = _simplifiedTDDFT->getGammaJ()[I][J];

      for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
        Eigen::MatrixXd X = (Jia.alpha.transpose() * guessVectors[iSet].middleRows(jbStart, alphaJ)).eval();
        X.noalias() += (Jia.beta.transpose() * guessVectors[iSet].middleRows(jbStart + alphaJ, betaJ)).eval();

        sigmaVectors[iSet].middleRows(iaStart, alphaI) += (1.0 + pm[iSet] * 1.0) * Jai.alpha * gammaK * X;
        sigmaVectors[iSet].middleRows(iaStart + alphaI, betaI) += (1.0 + pm[iSet] * 1.0) * Jai.beta * gammaK * X;
        if (I == J) {
          Eigen::MatrixXd tJija = Jij.alpha * gammaJ;
          Eigen::MatrixXd tJijb = Jij.beta * gammaJ;
          Eigen::MatrixXd tJaia = Jai.alpha * gammaK;
          Eigen::MatrixXd tJaib = Jai.beta * gammaK;

          for (unsigned iGuess = 0; iGuess < guessVectors[iSet].cols(); ++iGuess) {
            Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVectors[iSet].col(iGuess).data() + iaStart, nvI.alpha, noI.alpha);
            Eigen::Map<Eigen::MatrixXd> guess_a(guessVectors[iSet].col(iGuess).data() + jbStart, nvJ.alpha, noJ.alpha);

            Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVectors[iSet].col(iGuess).data() + iaStart + alphaI, nvI.beta, noI.beta);
            Eigen::Map<Eigen::MatrixXd> guess_b(guessVectors[iSet].col(iGuess).data() + jbStart + alphaJ, nvJ.beta, noJ.beta);

            for (unsigned iAtom = 0; iAtom < Jai.alpha.cols(); ++iAtom) {
              Eigen::Map<Eigen::MatrixXd> ji_a(tJija.col(iAtom).data(), noI.alpha, noI.alpha);
              Eigen::Map<Eigen::MatrixXd> ia_a(tJaia.col(iAtom).data(), nvI.alpha, noI.alpha);
              Eigen::Map<Eigen::MatrixXd> ai_a(Jai.alpha.col(iAtom).data(), nvI.alpha, noI.alpha);
              Eigen::Map<Eigen::MatrixXd> ab_a(Jab.alpha.col(iAtom).data(), nvI.alpha, nvI.alpha);

              Eigen::Map<Eigen::MatrixXd> ji_b(tJijb.col(iAtom).data(), noI.beta, noI.beta);
              Eigen::Map<Eigen::MatrixXd> ia_b(tJaib.col(iAtom).data(), nvI.beta, noI.beta);
              Eigen::Map<Eigen::MatrixXd> ai_b(Jai.beta.col(iAtom).data(), nvI.beta, noI.beta);
              Eigen::Map<Eigen::MatrixXd> ab_b(Jab.beta.col(iAtom).data(), nvI.beta, nvI.beta);

              sigma_a.noalias() -= ab_a * guess_a * ji_a;
              sigma_b.noalias() -= ab_b * guess_b * ji_b;
              if (pm[iSet]) {
                sigma_a.noalias() -= pm[iSet] * ax[I] * ai_a * guess_a.transpose() * ia_a;
                sigma_b.noalias() -= pm[iSet] * ax[I] * ai_b * guess_b.transpose() * ia_b;
              }
            }
          }
        }
      }
      jbStart += alphaJ + betaJ;
    }
    iaStart += alphaI + betaI;
  }

  Timings::timeTaken("LRSCF -   Sigmavector: Grimme");
  return sigmaVectors;
} /* this->getSigmavectors() unrestricted */

template class GrimmeSigmavector<Options::SCF_MODES::RESTRICTED>;
template class GrimmeSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
