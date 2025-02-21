/**
 * @file Sigmavector.cpp
 *
 * @date Dec 06, 2018
 * @author Michael Boeckers
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
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "data/OrbitalPair.h"
#include "data/PAOController.h"
#include "io/FormattedOutputStream.h"
#include "memory/MemoryManager.h"
#include "misc/Timing.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/LRSCFOptions.h"
#include "settings/Settings.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
Sigmavector<SCFMode>::Sigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b)
  : _lrscf(lrscf), _nSub(lrscf.size()), _nGuess(b[0].cols()), _nSet(b.size()), _hasBeenCalculated(false) {
  // Construct guess vectors for each subsystem for easier access (b holds guess vectors of total system,
  // i.e. the full subsystem response problem).
  _b.resize(_nSet);
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    _b[iSet].resize(_nSub);
    unsigned iCount = 0;
    for (unsigned I = 0; I < _nSub; ++I) {
      // Get position of guess vector for subsystem I
      auto nOccI = _lrscf[I]->getNOccupied();
      auto nVirtI = _lrscf[I]->getNVirtual();
      unsigned nDimI = 0;
      for_spin(nOccI, nVirtI) {
        nDimI += nOccI_spin * nVirtI_spin;
        iCount += nOccI_spin * nVirtI_spin;
      };
      unsigned iStartI = iCount - nDimI;
      _b[iSet][I] = b[iSet].block(iStartI, 0, nDimI, _nGuess);
    }
  }
  // Set dimensions for sigma vector and null
  _sigma.resize(b.size());
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    _sigma[iSet].resize(b[iSet].rows(), _nGuess);
    _sigma[iSet].setZero();
  }
}

template<Options::SCF_MODES SCFMode>
Sigmavector<SCFMode>::Sigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf)
  : _lrscf(lrscf), _nSub(lrscf.size()), _nGuess(1), _nSet(1), _hasBeenCalculated(false) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> Sigmavector<SCFMode>::calcP(unsigned J) {
  Timings::takeTime("LRSCF -  AO-MO Transformation");
  auto P_J = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(_nSet);
  auto C = _lrscf[J]->getCoefficients();
  auto no = _lrscf[J]->getNOccupied();
  auto nv = _lrscf[J]->getNVirtual();
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
      (*P_J)[iSet].emplace_back(_lrscf[J]->getBasisController());
      auto& P = (*P_J)[iSet][iGuess];
      unsigned iStartSpin = 0;
      for_spin(P, C, no, nv) {
        Eigen::Map<Eigen::MatrixXd> b(_b[iSet][J].col(iGuess).data() + iStartSpin, nv_spin, no_spin);
        P_spin = C_spin.middleCols(0, no_spin) * b.transpose() * C_spin.middleCols(no_spin, nv_spin).transpose();
        iStartSpin += nv_spin * no_spin;
      };
    }
  }
  Timings::timeTaken("LRSCF -  AO-MO Transformation");
  return P_J;
}

template<Options::SCF_MODES SCFMode>
void Sigmavector<SCFMode>::addToSigma(unsigned I, unsigned iStartI,
                                      std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> F_I) {
  Timings::takeTime("LRSCF -  AO-MO Transformation");
  auto C = _lrscf[I]->getCoefficients();
  auto no = _lrscf[I]->getNOccupied();
  auto nv = _lrscf[I]->getNVirtual();
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
      auto& F = (*F_I)[iSet][iGuess];
      unsigned iStartSpin = 0;
      for_spin(F, C, no, nv) {
        F_spin = C_spin.middleCols(no_spin, nv_spin).transpose() * F_spin.transpose() * C_spin.middleCols(0, no_spin);
        _sigma[iSet].block(iStartI + iStartSpin, iGuess, nv_spin * no_spin, 1) +=
            Eigen::Map<Eigen::VectorXd>(F_spin.data(), F_spin.size());
        F_spin.resize(0, 0);
        iStartSpin += nv_spin * no_spin;
      };
    }
  }
  F_I.reset();
  Timings::timeTaken("LRSCF -  AO-MO Transformation");
}

template<Options::SCF_MODES SCFMode>
void Sigmavector<SCFMode>::calcSigma() {
  // Set sigma vectors for each set to zero.
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    _sigma[iSet].setZero();
  }
  // Calculate matrix vector product M_IJ*b_J for each set.
  unsigned iCount = 0;
  for (unsigned I = 0; I < _nSub; ++I) {
    auto nOccI = _lrscf[I]->getNOccupied();
    auto nVirtI = _lrscf[I]->getNVirtual();
    unsigned nDimI = 0;
    for_spin(nOccI, nVirtI) {
      nDimI += nOccI_spin * nVirtI_spin;
      iCount += nOccI_spin * nVirtI_spin;
    };
    unsigned iStartI = iCount - nDimI;

    for (unsigned J = 0; J < _nSub; ++J) {
      if (this->skipThisInteraction(I, J)) {
        continue;
      }

      auto P_J = this->calcP(J);

      // Determine prescreening threshold.
      _prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                        this->_lrscf[J]->getSysSettings().basis.integralThreshold);

      if (_prescreeningThreshold == 0) {
        double prescreenTresholdI = this->_lrscf[I]->getBasisController()->getPrescreeningThreshold();
        double prescreenTresholdJ = this->_lrscf[J]->getBasisController()->getPrescreeningThreshold();
        _prescreeningThreshold = std::min(prescreenTresholdI, prescreenTresholdJ);
      }

      // Make prescreening threshold adaptive.
      double minGuessVectorNorm = std::numeric_limits<double>::infinity();
      if (_lrscf[J]->getLRSCFSettings().adaptivePrescreening) {
        for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
          minGuessVectorNorm = std::min(_b[iSet][J].colwise().norm().minCoeff(), minGuessVectorNorm);
        }
        // Start with a prescreening threshold of 1e-9 for standard settings.
        double startScreening = _lrscf[I]->getLRSCFSettings().conv * 0.0001;
        // Make it smaller based on the minimum guess vector norm (so that every needed integral is calculated).
        double adaptive = startScreening * minGuessVectorNorm;
        _prescreeningThreshold = std::max(adaptive, _prescreeningThreshold);
      }

      // Restrict number of Fock matrix copies for integral-direct algorithms.
      size_t nb_I = _lrscf[I]->getBasisController()->getNBasisFunctions();
      size_t memDemand = nb_I * nb_I * _nGuess * _nSet * sizeof(double);
      size_t available = 0.8 * MemoryManager::getInstance()->getAvailableSystemMemory();

      unsigned nFock = available / memDemand;
      _nThreads = std::min(nFock, _maxThreads);
      if (_nThreads < _maxThreads) {
        OutputControl::nOut << "   Maximum Threads            : " << _maxThreads << std::endl;
        OutputControl::nOut << "   Used Threads               : " << _nThreads << std::endl;
        OutputControl::nOut << "   Number of Fock matrices    : " << _nGuess * _nSet << std::endl;
        OutputControl::nOut << "   Memory demand for max (GB) : "
                            << (nb_I * nb_I * _nGuess * _nSet * sizeof(double) * _maxThreads) / 1e9 << std::endl;
        OutputControl::nOut << "   Actual memory demand  (GB) : "
                            << (nb_I * nb_I * _nGuess * _nSet * sizeof(double) * _nThreads) / 1e9 << std::endl;
      }
      omp_set_num_threads(_nThreads);
      auto F_IJ = this->calcF(I, J, std::move(P_J));
      omp_set_num_threads(_maxThreads);

      if (F_IJ) {
        this->addToSigma(I, iStartI, std::move(F_IJ));
      }
    }
  }
  _hasBeenCalculated = true;
}

template<Options::SCF_MODES SCFMode>
std::vector<MatrixInBasis<SCFMode>> Sigmavector<SCFMode>::getPerturbedFockMatrix() {
  // Calculate matrix vector product M_IJ*b_J for each set.
  unsigned iCount = 0;
  if (_nSub > 1)
    throw SerenityError("More than one subsystem not supported!");

  auto nOccI = _lrscf[0]->getNOccupied();
  auto nVirtI = _lrscf[0]->getNVirtual();
  unsigned nDimI = 0;
  for_spin(nOccI, nVirtI) {
    nDimI += nOccI_spin * nVirtI_spin;
    iCount += nOccI_spin * nVirtI_spin;
  };
  auto P_J = calcP(0);

  // Simple significance screening.
  double maxDens = 0.0;
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
      auto& P = (*P_J)[iSet][iGuess];
      for_spin(P) {
        maxDens = std::max(maxDens, P_spin.array().abs().maxCoeff());
      };
    }
  }

  // Restrict number of Fock matrix copies for integral-direct algorithms.
  size_t nb_I = _lrscf[0]->getBasisController()->getNBasisFunctions();
  size_t memDemand = nb_I * nb_I * _nGuess * _nSet * sizeof(double);
  size_t available = 0.8 * MemoryManager::getInstance()->getAvailableSystemMemory();

  unsigned nFock = available / memDemand;
  _nThreads = std::min(nFock, _maxThreads);
  if (_nThreads < _maxThreads) {
    OutputControl::nOut << "   Maximum Threads            : " << _maxThreads << std::endl;
    OutputControl::nOut << "   Used Threads               : " << _nThreads << std::endl;
    OutputControl::nOut << "   Number of Fock matrices    : " << _nGuess * _nSet << std::endl;
    OutputControl::nOut << "   Memory demand for max (GB) : "
                        << (nb_I * nb_I * _nGuess * _nSet * sizeof(double) * _maxThreads) / 1e9 << std::endl;
    OutputControl::nOut << "   Actual memory demand  (GB) : "
                        << (nb_I * nb_I * _nGuess * _nSet * sizeof(double) * _nThreads) / 1e9 << std::endl;
  }
  omp_set_num_threads(_nThreads);
  auto F_IJ = this->calcF(0, 0, std::move(P_J));
  omp_set_num_threads(_maxThreads);

  auto perturbedFockMatrix = (*F_IJ)[0];
  if ((*F_IJ).size() > 1) {
    for (unsigned iState = 0; iState < perturbedFockMatrix.size(); iState++) {
      perturbedFockMatrix[iState] += (*F_IJ)[1][iState];
    }
  }

  return perturbedFockMatrix;
}

template<Options::SCF_MODES SCFMode>
void Sigmavector<SCFMode>::setShellWiseMaxDens(unsigned J, std::vector<std::vector<MatrixInBasis<SCFMode>>>& densityMatrices) {
  Timings::takeTime("LRSCF - Shell-Wise Abs. Calc.");
  const auto& basis = this->_lrscf[J]->getBasisController()->getBasis();
  unsigned nShells = basis.size();
  _maxDensMat = Eigen::MatrixXd::Zero(nShells, nShells);

  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      auto blockMax = densityMatrices[iSet][iGuess].shellWiseAbsMax().total();
      for (unsigned iShell = 0; iShell < nShells; ++iShell) {
        for (unsigned jShell = 0; jShell < nShells; ++jShell) {
          _maxDensMat(iShell, jShell) = std::max(_maxDensMat(iShell, jShell), blockMax(iShell, jShell));
        }
      }
    }
  }
  _maxDens = _maxDensMat.lpNorm<Eigen::Infinity>();
  Timings::timeTaken("LRSCF - Shell-Wise Abs. Calc.");
}

template<Options::SCF_MODES SCFMode>
bool Sigmavector<SCFMode>::skipThisInteraction(unsigned I, unsigned J) {
  // Do not calculate off-diagonal blocks whatsoever if requested.
  if (_lrscf[I]->getLRSCFSettings().noCoupling && I != J) {
    return true;
  }

  bool exploitSymmetry = (_lrscf.size() > 1 && !_lrscf[I]->getLRSCFSettings().fullFDEc);
  exploitSymmetry = exploitSymmetry && _lrscf[I]->getLRSCFSettings().partialResponseConstruction;

  if (exploitSymmetry) {
    if (_lrscf[I]->getLRSCFSettings().method == Options::LR_METHOD::TDA) {
      if (I == J || I > J) {
        return true;
      }
    }
    else if (_lrscf[I]->getLRSCFSettings().method == Options::LR_METHOD::TDDFT) {
      if (I > J) {
        return true;
      }
    }
  }

  // Do not calculate diagonal blocks for CC2.
  bool isCC2 = _lrscf[I]->getLRSCFSettings().method != Options::LR_METHOD::TDA &&
               _lrscf[I]->getLRSCFSettings().method != Options::LR_METHOD::TDDFT;
  if (isCC2 && I == J) {
    return true;
  }

  return false;
}

template class Sigmavector<Options::SCF_MODES::RESTRICTED>;
template class Sigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
