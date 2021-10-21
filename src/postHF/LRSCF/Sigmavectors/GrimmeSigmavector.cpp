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
  Timings::takeTime("LRSCF -    Grimme Sigmavector");
  std::vector<Eigen::MatrixXd> sigmaVectors = guessVectors;
  unsigned nSub = _lrscf.size();

  auto& Jij = _simplifiedTDDFT->getJij();
  auto& Jai = _simplifiedTDDFT->getJai();
  auto& Jab = _simplifiedTDDFT->getJab();

  double ax = _simplifiedTDDFT->getHFExchangeRatio();
  long iaStart = 0;
  for (unsigned I = 0; I < nSub; ++I) {
    unsigned noI = _lrscf[I]->getNOccupied();
    unsigned nvI = _lrscf[I]->getNVirtual();

    long jbStart = 0;
    for (unsigned J = 0; J < nSub; ++J) {
      unsigned noJ = _lrscf[J]->getNOccupied();
      unsigned nvJ = _lrscf[J]->getNVirtual();

      for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
        sigmaVectors[iSet] = (2.0 + pm[iSet] * 2.0) * Jai * Jai.transpose() * guessVectors[iSet];
        for (unsigned iGuess = 0; iGuess < guessVectors[iSet].cols(); ++iGuess) {
          Eigen::Map<Eigen::MatrixXd> sigma(sigmaVectors[iSet].col(iGuess).data(), nvI, noI);
          Eigen::Map<Eigen::MatrixXd> guess(guessVectors[iSet].col(iGuess).data(), nvJ, noJ);
          for (unsigned iAtom = 0; iAtom < Jai.cols(); ++iAtom) {
            Eigen::Map<Eigen::MatrixXd> ji(Jij.col(iAtom).data(), noI, noI);
            Eigen::Map<Eigen::MatrixXd> ai(Jai.col(iAtom).data(), nvI, noI);
            Eigen::Map<Eigen::MatrixXd> ab(Jab.col(iAtom).data(), nvI, nvI);
            sigma.noalias() -= ab * guess * ji;
            if (pm[iSet]) {
              sigma.noalias() += pm[iSet] * ax * ai * guess.transpose() * ai;
            }
          }
        }
      }
      jbStart += nvI * noI;
    }
    iaStart += nvI * noI;
  }

  Timings::timeTaken("LRSCF -    Grimme Sigmavector");
  return sigmaVectors;
} /* this->getSigmavectors() restricted */

template<>
std::vector<Eigen::MatrixXd> GrimmeSigmavector<UNRESTRICTED>::getSigmavectors(std::vector<Eigen::MatrixXd>& guessVectors,
                                                                              std::vector<int> pm) {
  Timings::takeTime("LRSCF -    Grimme Sigmavector");
  std::vector<Eigen::MatrixXd> sigmaVectors = guessVectors;
  unsigned nSub = _lrscf.size();

  auto& Jij = _simplifiedTDDFT->getJij();
  auto& Jai = _simplifiedTDDFT->getJai();
  auto& Jab = _simplifiedTDDFT->getJab();

  double ax = _simplifiedTDDFT->getHFExchangeRatio();

  long iaStart = 0;
  for (unsigned I = 0; I < nSub; ++I) {
    auto noI = _lrscf[I]->getNOccupied();
    auto nvI = _lrscf[I]->getNVirtual();
    unsigned alphaI = nvI.alpha * noI.alpha;
    unsigned betaI = nvI.beta * noI.beta;

    long jbStart = 0;
    for (unsigned J = 0; J < nSub; ++J) {
      auto noJ = _lrscf[J]->getNOccupied();
      auto nvJ = _lrscf[J]->getNVirtual();
      unsigned alphaJ = nvJ.alpha * noJ.alpha;
      unsigned betaJ = nvJ.beta * noJ.beta;

      for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
        Eigen::MatrixXd X = Jai.alpha.transpose() * guessVectors[iSet].topRows(alphaJ);
        X.noalias() += Jai.beta.transpose() * guessVectors[iSet].bottomRows(betaJ);
        sigmaVectors[iSet].topRows(alphaI) = (1.0 + pm[iSet] * 1.0) * Jai.alpha * X;
        sigmaVectors[iSet].bottomRows(betaI) = (1.0 + pm[iSet] * 1.0) * Jai.beta * X;
        for (unsigned iGuess = 0; iGuess < guessVectors[iSet].cols(); ++iGuess) {
          Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVectors[iSet].col(iGuess).data(), nvI.alpha, noI.alpha);
          Eigen::Map<Eigen::MatrixXd> guess_a(guessVectors[iSet].col(iGuess).data(), nvJ.alpha, noJ.alpha);

          Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVectors[iSet].col(iGuess).data() + alphaI, nvI.beta, noI.beta);
          Eigen::Map<Eigen::MatrixXd> guess_b(guessVectors[iSet].col(iGuess).data() + alphaJ, nvJ.beta, noJ.beta);

          for (unsigned iAtom = 0; iAtom < Jai.alpha.cols(); ++iAtom) {
            Eigen::Map<Eigen::MatrixXd> ji_a(Jij.alpha.col(iAtom).data(), noI.alpha, noI.alpha);
            Eigen::Map<Eigen::MatrixXd> ai_a(Jai.alpha.col(iAtom).data(), nvI.alpha, noI.alpha);
            Eigen::Map<Eigen::MatrixXd> ab_a(Jab.alpha.col(iAtom).data(), nvI.alpha, nvI.alpha);

            Eigen::Map<Eigen::MatrixXd> ji_b(Jij.beta.col(iAtom).data(), noI.beta, noI.beta);
            Eigen::Map<Eigen::MatrixXd> ai_b(Jai.beta.col(iAtom).data(), nvI.beta, noI.beta);
            Eigen::Map<Eigen::MatrixXd> ab_b(Jab.beta.col(iAtom).data(), nvI.beta, nvI.beta);

            sigma_a.noalias() -= ab_a * guess_a * ji_a;
            sigma_b.noalias() -= ab_b * guess_b * ji_b;
            if (pm[iSet]) {
              sigma_a.noalias() += pm[iSet] * ax * ai_a * guess_a.transpose() * ai_a;
              sigma_b.noalias() += pm[iSet] * ax * ai_b * guess_b.transpose() * ai_b;
            }
          }
        }
      }
      jbStart += alphaJ + betaJ;
    }
    iaStart += alphaI + betaI;
  }

  Timings::timeTaken("LRSCF -    Grimme Sigmavector");
  return sigmaVectors;
} /* this->getSigmavectors() unrestricted */

template class GrimmeSigmavector<Options::SCF_MODES::RESTRICTED>;
template class GrimmeSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
