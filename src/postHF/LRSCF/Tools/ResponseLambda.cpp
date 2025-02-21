/**
 * @file ResponseLambda.cpp
 * @date Nov 13, 2020
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
#include "postHF/LRSCF/Tools/ResponseLambda.h"
/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "postHF/LRSCF/Kernel/Kernel.h"
#include "postHF/LRSCF/RICC2/CC2Controller.h"
#include "postHF/LRSCF/Sigmavectors/CoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/EOSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/ExchangeSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/FockSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/GrimmeSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/KernelSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RI/MRICoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RI/RICoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RI/RIExchangeSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
#include "postHF/LRSCF/Tools/SimplifiedTDDFT.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ResponseLambda<SCFMode>::~ResponseLambda() = default;

template<Options::SCF_MODES SCFMode>
ResponseLambda<SCFMode>::ResponseLambda(std::vector<std::shared_ptr<SystemController>> act,
                                        std::vector<std::shared_ptr<SystemController>> env,
                                        std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                        const Eigen::VectorXd& diagonal, LRSCFTaskSettings& settings)
  : _act(act), _env(env), _lrscf(lrscf), _diagonal(diagonal), _settings(settings) {
  // Kernel to be used
  for (const auto& sys : _act) {
    if (sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      if (sys->getSettings().dft.functional != CompositeFunctionals::XCFUNCTIONALS::HF) {
        if (sys->getSettings().dft.functional != CompositeFunctionals::XCFUNCTIONALS::NONE) {
          if (_settings.func != CompositeFunctionals::XCFUNCTIONALS::HF) {
            _usesKernel = true;
          }
        }
      }
    }
  }
  for (const auto& sys : _env) {
    if (sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      if (sys->getSettings().dft.functional != CompositeFunctionals::XCFUNCTIONALS::HF) {
        if (sys->getSettings().dft.functional != CompositeFunctionals::XCFUNCTIONALS::NONE) {
          _usesKernel = true;
        }
      }
    }
  }
  if (settings.embedding.naddXCFunc != CompositeFunctionals::XCFUNCTIONALS::NONE) {
    if (settings.embedding.naddXCFunc != CompositeFunctionals::XCFUNCTIONALS::HF) {
      _usesKernel = true;
    }
  }
  if (settings.embedding.naddKinFunc != CompositeFunctionals::KINFUNCTIONALS::NONE) {
    if (settings.embedding.naddXCFunc != CompositeFunctionals::XCFUNCTIONALS::HF) {
      _usesKernel = true;
    }
  }

  // If requested by the input, skip the calculation of the XC kernel.
  if (_settings.noKernel) {
    _usesKernel = false;
  }

  // Exchange to be used
  for (unsigned I = 0; I < _lrscf.size(); ++I) {
    if (_lrscf[I]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
      _usesExchange = true;
    }
    else {
      auto func_enum = this->_lrscf[I]->getLRSCFSettings().func;
      auto func = this->_lrscf[I]->getLRSCFSettings().customFunc.basicFunctionals.size()
                      ? Functional(this->_lrscf[I]->getLRSCFSettings().customFunc)
                      : resolveFunctional(func_enum);

      // Resolve functional.
      if (func_enum == CompositeFunctionals::XCFUNCTIONALS::NONE) {
        func = this->_lrscf[I]->getSysSettings().customFunc.basicFunctionals.size()
                   ? Functional(this->_lrscf[I]->getSysSettings().customFunc)
                   : resolveFunctional(_lrscf[I]->getSysSettings().dft.functional);
      }
      _usesExchange = func.isHybrid() ? func.isHybrid() : _usesExchange;
      _usesLRExchange = func.isRSHybrid() ? func.isRSHybrid() : _usesLRExchange;
      _usesDoubleHybrid = func.isDoubleHybrid() ? func.isDoubleHybrid() : _usesDoubleHybrid;
    }
  }
  auto naddXCfunc = _settings.embedding.customNaddXCFunc.basicFunctionals.size()
                        ? Functional(_settings.embedding.customNaddXCFunc)
                        : resolveFunctional(_settings.embedding.naddXCFunc);
  _usesExchange = naddXCfunc.isHybrid() ? naddXCfunc.isHybrid() : _usesExchange;
  _usesLRExchange = naddXCfunc.isRSHybrid() ? naddXCfunc.isRSHybrid() : _usesLRExchange;
  for (auto naddXCfunc : _settings.embedding.naddXCFuncList) {
    auto func = resolveFunctional(naddXCfunc);
    _usesExchange = func.isHybrid() ? func.isHybrid() : _usesExchange;
    _usesLRExchange = func.isRSHybrid() ? func.isRSHybrid() : _usesLRExchange;
  }

  // Simplified TDDFT to be used
  if (_settings.grimme || _settings.approxCoulomb[1] != std::numeric_limits<double>::infinity()) {
    _simplifiedTDDFT = std::make_unique<SimplifiedTDDFT<SCFMode>>(_lrscf);
    if (_settings.grimme) {
      _usesKernel = false;
      _usesExchange = false;
      _usesLRExchange = false;
      _usesDoubleHybrid = false;
    }
  }

  if (_usesLRExchange && _settings.rpaScreening) {
    throw SerenityError("RS-Hybrid Functional not supported for BSE calculations!");
  }

  if (_usesDoubleHybrid) {
    printBigCaption("Double-Hybrid TDDFT");
    Functional func = _lrscf[0]->getSysSettings().customFunc.basicFunctionals.size()
                          ? Functional(_lrscf[0]->getSysSettings().customFunc)
                          : resolveFunctional(_lrscf[0]->getSysSettings().dft.functional);
    _doubleHybridCorrRatio = func.getHfCorrelRatio();
    // Set spin scaling factors from functional.
    _settings.sss = func.getssScaling();
    _settings.oss = func.getosScaling();
    printf("  WF corr. ratio        : %7.3f\n", _doubleHybridCorrRatio);
    printf("  Same-spin scaling     : %7.3f\n", _settings.sss);
    printf("  Opposite-spin scaling : %7.3f\n\n", _settings.oss);
    if (_usesDoubleHybrid && _lrscf.size() > 1) {
      throw SerenityError("Cannot use double-hybrids in a coupled calculation!");
    }
  }

  // EO to be used
  if (_settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
    _usesEO = true;
  }
  if (_settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::HUZINAGA) {
    _usesEO = true;
  }
  if (_settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::HOFFMANN) {
    _usesEO = true;
  }
  if (_settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA) {
    _usesEO = true;
  }
  for (auto eo : _settings.embedding.embeddingModeList) {
    if (eo == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
      _usesEO = true;
    }
    if (eo == Options::KIN_EMBEDDING_MODES::HUZINAGA) {
      _usesEO = true;
    }
    if (eo == Options::KIN_EMBEDDING_MODES::HOFFMANN) {
      _usesEO = true;
    }
    if (eo == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA) {
      _usesEO = true;
    }
  }

  _densFitJ = _settings.densFitJ != Options::DENS_FITS::NONE;
  _densFitK = _settings.densFitK != Options::DENS_FITS::NONE;
  _densFitLRK = _settings.densFitLRK != Options::DENS_FITS::NONE;

  // Setup RI Integral cache for ADC(2) or CC2.
  bool isCC2 = !(settings.method == Options::LR_METHOD::TDA || settings.method == Options::LR_METHOD::TDDFT);
  if (isCC2 || _settings.rpaScreening || _usesDoubleHybrid) {
    if (!_usesDoubleHybrid) {
      _densFitK = true;
      _settings.densFitK = _settings.densFitCache;
    }
    for (auto& lrscf : _lrscf) {
      lrscf->initializeRIIntegrals(LIBINT_OPERATOR::coulomb, 0.0, true);
      if (_settings.aocache) {
        lrscf->getRIIntegrals(LIBINT_OPERATOR::coulomb)->cacheAOIntegrals();
      }
    }
  }
  if (isCC2 || _usesDoubleHybrid) {
    if (_settings.sss != 1.0 || _settings.oss != 1.0) {
      printBigCaption("Custom spin-component scaling");
      printf(" Same-spin     : %-6.3f\n", _settings.sss);
      printf(" Opposite-spin : %-6.3f\n\n", _settings.oss);
    }
    if (_settings.ltconv != 0) {
      printBigCaption("Laplace-Transform algorithms");
      printf(" Threshold     : %-6.1e\n", _settings.ltconv);
      printf(" Opposite-spin : %-6.3f\n\n", _settings.oss);
    }
    printSubSectionTitle("Ground State");
    for (auto& lrscf : _lrscf) {
      // In case of double-hybrids, we fake the method to be CIS(D) to get the
      // reduced sigma vector build which does not contain the CIS contributions.
      if (_usesDoubleHybrid) {
        if (!settings.frequencies.empty()) {
          throw SerenityError("LRSCFTask: Response functions are not properly defined for double-hybrid functionals!");
        }
        auto oldmethod = _settings.method;
        _settings.method = Options::LR_METHOD::CISD;
        lrscf->initializeCC2Controller();
        _settings.method = oldmethod;
      }
      else {
        lrscf->initializeCC2Controller();
      }
    }
  }

  if (_settings.rpaScreening) {
    _usesKernel = false;
    _usesExchange = true;
    Timings::takeTime("LRSCF -         RPA Screening");
    _lrscf[0]->calculateScreening(_diagonal);
    Timings::timeTaken("LRSCF -         RPA Screening");
  }

  if (_usesKernel) {
    printBigCaption("Kernel");
    Timings::takeTime("LRSCF -  Kernel on Grid Eval.");
    this->setupKernel(Options::GRID_PURPOSES::SMALL);
    printf("  Grid:    %1i/%1i\n\n", _settings.grid.smallGridAccuracy, _settings.grid.accuracy);
    printf(" .. done.\n\n");
    Timings::timeTaken("LRSCF -  Kernel on Grid Eval.");
  }

  if (!_settings.grimme) {
    printBigCaption("Density Fitting");
    std::string fitting;
    Options::resolve<Options::DENS_FITS>(fitting, _settings.densFitJ);
    printf("  Coulomb      :  %-10s\n", fitting.c_str());
    if (_usesExchange) {
      std::string fitting;
      Options::resolve<Options::DENS_FITS>(fitting, _settings.densFitK);
      printf("  Exchange     :  %-10s\n", fitting.c_str());
    }
    if (_usesLRExchange) {
      std::string fitting;
      Options::resolve<Options::DENS_FITS>(fitting, _settings.densFitLRK);
      printf("  LR-Exchange  :  %-10s\n", fitting.c_str());
    }
    printf("\n");

    if (_settings.adaptivePrescreening) {
      printf("  Using adaptive prescreening thresholds.\n\n");
    }
  }
}

template<Options::SCF_MODES SCFMode>
void ResponseLambda<SCFMode>::setupTDDFTLambdas() {
  // Sigmacalculator: (A+B)*b stored in [0,2,4,..], (A-B)*b stored in [1,3,5,..]
  _RPA = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::make_unique<std::vector<Eigen::MatrixXd>>(
        guessVectors.size(), Eigen::MatrixXd::Zero(guessVectors[0].rows(), guessVectors[0].cols()));

    std::vector<int> pm(guessVectors.size());
    std::vector<Eigen::MatrixXd> guessAPB(guessVectors.size() / 2);
    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      // (A+B)
      if (iSet % 2 == 0) {
        pm[iSet] = 1;
        guessAPB[iSet / 2] = guessVectors[iSet];
      }
      // (A-B)
      else {
        pm[iSet] = -1;
      }
    }

    D = std::make_unique<FockSigmavector<SCFMode>>(_lrscf, guessVectors);

    if (!_settings.grimme) {
      if (_densFitJ) {
        J = std::make_unique<RICoulombSigmavector<SCFMode>>(_lrscf, guessAPB);
      }
      else {
        J = std::make_unique<CoulombSigmavector<SCFMode>>(_lrscf, guessAPB);
      }
      if (_settings.approxCoulomb[0] != std::numeric_limits<double>::infinity()) {
        DFJ = std::make_unique<MRICoulombSigmavector<SCFMode>>(_lrscf, guessAPB);
      }
    }

    if (_usesExchange || _usesLRExchange) {
      if (_densFitK || _densFitLRK) {
        DFK = std::make_unique<RIExchangeSigmavector<SCFMode>>(_lrscf, guessVectors, pm, _densFitK, _densFitLRK);
      }
      if (!_densFitK || !_densFitLRK) {
        K = std::make_unique<ExchangeSigmavector<SCFMode>>(_lrscf, guessVectors, pm, _densFitK, _densFitLRK);
      }
    }

    if (_usesKernel) {
      F = std::make_unique<KernelSigmavector<SCFMode>>(_lrscf, guessAPB, _kernel, _ukernel);
    }

    if (_usesEO) {
      EO = std::make_unique<EOSigmavector<SCFMode>>(_lrscf, guessVectors, _settings.embedding);
    }

    if (_simplifiedTDDFT) {
      G = std::make_unique<GrimmeSigmavector<SCFMode>>(_lrscf, guessVectors, pm, _simplifiedTDDFT);
    }

    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      (*sigma)[iSet] += D->getSigma()[iSet];
      if (iSet % 2 == 0) {
        if (!(SCFMode == RESTRICTED && _settings.triplet)) {
          if (J) {
            (*sigma)[iSet] += 2.0 * _scfFactor * J->getSigma()[iSet / 2];
          }
          if (DFJ) {
            (*sigma)[iSet] += 2.0 * _scfFactor * DFJ->getSigma()[iSet / 2];
          }
        }
        if (_usesKernel) {
          double tripletFactor = (_settings.triplet ? 0.5 : 1.0);
          (*sigma)[iSet] += 2.0 * tripletFactor * _scfFactor * F->getSigma()[iSet / 2];
        }
      }
      if (DFK) {
        (*sigma)[iSet] -= DFK->getSigma()[iSet];
      }
      if (K) {
        (*sigma)[iSet] -= K->getSigma()[iSet];
      }
      if (EO) {
        (*sigma)[iSet] += EO->getSigma()[iSet];
      }
      if (G) {
        (*sigma)[iSet] += G->getSigma()[iSet];
      }
    }

    return sigma;
  }; /* RPA */

  // Sigmacalculator: A*b
  // For stability analyses, this lambda is also used for the (A+B), (A-B), and SF-TDDFT matrices.
  _TDA = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::make_unique<std::vector<Eigen::MatrixXd>>(
        guessVectors.size(), Eigen::MatrixXd::Zero(guessVectors[0].rows(), guessVectors[0].cols()));

    D = std::make_unique<FockSigmavector<SCFMode>>(_lrscf, guessVectors);

    if (!_settings.grimme) {
      if (_densFitJ) {
        J = std::make_unique<RICoulombSigmavector<SCFMode>>(_lrscf, guessVectors);
      }
      else {
        J = std::make_unique<CoulombSigmavector<SCFMode>>(_lrscf, guessVectors);
      }
      if (_settings.approxCoulomb[0] != std::numeric_limits<double>::infinity()) {
        DFJ = std::make_unique<MRICoulombSigmavector<SCFMode>>(_lrscf, guessVectors);
      }
    }

    int exchange_sign = 0;
    if (_settings.scfstab == Options::STABILITY_ANALYSIS::REAL) {
      exchange_sign = 1;
    }
    else if (_settings.scfstab == Options::STABILITY_ANALYSIS::NONREAL) {
      exchange_sign = -1;
    }
    else if (_settings.scfstab == Options::STABILITY_ANALYSIS::SPINFLIP) {
      exchange_sign = 0;
    }

    std::vector<int> pm(guessVectors.size(), exchange_sign);

    if (_usesExchange || _usesLRExchange) {
      if (_densFitK || _densFitLRK) {
        DFK = std::make_unique<RIExchangeSigmavector<SCFMode>>(_lrscf, guessVectors, pm, _densFitK, _densFitLRK);
      }
      if (!_densFitK || !_densFitLRK) {
        K = std::make_unique<ExchangeSigmavector<SCFMode>>(_lrscf, guessVectors, pm, _densFitK, _densFitLRK);
      }
    }

    if (_usesKernel) {
      F = std::make_unique<KernelSigmavector<SCFMode>>(_lrscf, guessVectors, _kernel, _ukernel);
    }

    if (_usesEO) {
      EO = std::make_unique<EOSigmavector<SCFMode>>(_lrscf, guessVectors, _settings.embedding);
    }

    if (_simplifiedTDDFT) {
      G = std::make_unique<GrimmeSigmavector<SCFMode>>(_lrscf, guessVectors, pm, _simplifiedTDDFT);
    }

    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      (*sigma)[iSet] += D->getSigma()[iSet];
    }

    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      if (_settings.scfstab == Options::STABILITY_ANALYSIS::NONE || _settings.scfstab == Options::STABILITY_ANALYSIS::REAL) {
        double aplusbFactor = (_settings.scfstab == Options::STABILITY_ANALYSIS::REAL ? 2.0 : 1.0);
        if (!(SCFMode == RESTRICTED && _settings.triplet)) {
          if (J) {
            (*sigma)[iSet] += aplusbFactor * _scfFactor * J->getSigma()[iSet];
          }
          if (DFJ) {
            (*sigma)[iSet] += aplusbFactor * _scfFactor * DFJ->getSigma()[iSet];
          }
        }
        if (_usesKernel) {
          double tripletFactor = (_settings.triplet ? 0.5 : 1.0);
          (*sigma)[iSet] += tripletFactor * aplusbFactor * _scfFactor * F->getSigma()[iSet];
        }
      }

      if (DFK) {
        (*sigma)[iSet] -= DFK->getSigma()[iSet];
      }
      if (K) {
        (*sigma)[iSet] -= K->getSigma()[iSet];
      }
      if (EO) {
        (*sigma)[iSet] += EO->getSigma()[iSet];
      }
      if (G) {
        (*sigma)[iSet] += G->getSigma()[iSet];
      }
    }

    return sigma;
  }; /* TDA */
}

template<Options::SCF_MODES SCFMode>
void ResponseLambda<SCFMode>::setupCC2Lambdas() {
  // lambda function for right transformation CC2/CIS(Dinf)/ADC(2)
  _rightCC2 = [&](Eigen::MatrixXd& guessVectors, Eigen::VectorXd guessValues) {
    auto sigma = std::make_unique<Eigen::MatrixXd>(guessVectors.rows(), guessVectors.cols());

    for (unsigned iGuess = 0; iGuess < guessVectors.cols(); ++iGuess) {
      long iStart = 0;
      for (auto& lrscf : _lrscf) {
        auto no = lrscf->getNOccupied();
        auto nv = lrscf->getNVirtual();
        long nDimI = 0;
        for_spin(nv, no) {
          nDimI += nv_spin * no_spin;
        };
        Eigen::Ref<Eigen::VectorXd> sigmaI = (*sigma).col(iGuess).segment(iStart, nDimI);
        Eigen::Ref<Eigen::VectorXd> guessI = guessVectors.col(iGuess).segment(iStart, nDimI);
        sigmaI = lrscf->getCC2Controller()->getRightCC2Sigma(guessI, guessValues(iGuess));
        iStart += nDimI;
      }
    }

    // CC2 coupling.
    if (_lrscf.size() > 1) {
      std::vector<Eigen::MatrixXd> guess = {guessVectors};
      J = std::make_unique<RICoulombSigmavector<SCFMode>>(_lrscf, guess);
      if (!(SCFMode == RESTRICTED && _settings.triplet)) {
        (*sigma) += _scfFactor * J->getSigma()[0];
      }
      if (_usesExchange || _usesLRExchange) {
        K = std::make_unique<RIExchangeSigmavector<SCFMode>>(_lrscf, guess, std::vector<int>(1, 0), true, true);
        (*sigma) -= K->getSigma()[0];
      }
      if (_usesKernel) {
        F = std::make_unique<KernelSigmavector<SCFMode>>(_lrscf, guess, _kernel, _ukernel);
        double tripletFactor = (_settings.triplet ? 0.5 : 1.0);
        (*sigma) += tripletFactor * _scfFactor * F->getSigma()[0];
      }
      if (_usesEO) {
        EO = std::make_unique<EOSigmavector<SCFMode>>(_lrscf, guess, _settings.embedding);
        (*sigma) += EO->getSigma()[0];
      }
    }

    return sigma;
  };

  // lambda function for left transformation CC2/CIS(Dinf)/ADC(2)
  _leftCC2 = [&](Eigen::MatrixXd& guessVectors, Eigen::VectorXd guessValues) {
    auto sigma = std::make_unique<Eigen::MatrixXd>(guessVectors.rows(), guessVectors.cols());

    for (unsigned iGuess = 0; iGuess < guessVectors.cols(); ++iGuess) {
      long iStart = 0;
      for (auto& lrscf : _lrscf) {
        auto no = lrscf->getNOccupied();
        auto nv = lrscf->getNVirtual();
        long nDimI = 0;
        for_spin(nv, no) {
          nDimI += nv_spin * no_spin;
        };
        Eigen::Ref<Eigen::VectorXd> sigmaI = (*sigma).col(iGuess).segment(iStart, nDimI);
        Eigen::Ref<Eigen::VectorXd> guessI = guessVectors.col(iGuess).segment(iStart, nDimI);
        sigmaI = lrscf->getCC2Controller()->getLeftCC2Sigma(guessI, guessValues(iGuess));
        iStart += nDimI;
      }
    }

    // CC2 coupling.
    if (_lrscf.size() > 1) {
      std::vector<Eigen::MatrixXd> guess = {guessVectors};
      J = std::make_unique<RICoulombSigmavector<SCFMode>>(_lrscf, guess);
      if (!(SCFMode == RESTRICTED && _settings.triplet)) {
        (*sigma) += _scfFactor * J->getSigma()[0];
      }
      if (_usesExchange || _usesLRExchange) {
        K = std::make_unique<RIExchangeSigmavector<SCFMode>>(_lrscf, guess, std::vector<int>(1, 0), true, true);
        (*sigma) -= K->getSigma()[0];
      }
      if (_usesKernel) {
        F = std::make_unique<KernelSigmavector<SCFMode>>(_lrscf, guess, _kernel, _ukernel);
        double tripletFactor = (_settings.triplet ? 0.5 : 1.0);
        (*sigma) += tripletFactor * _scfFactor * F->getSigma()[0];
      }
      if (_usesEO) {
        EO = std::make_unique<EOSigmavector<SCFMode>>(_lrscf, guess, _settings.embedding);
        (*sigma) += EO->getSigma()[0];
      }
    }

    return sigma;
  };
}

template<Options::SCF_MODES SCFMode>
void ResponseLambda<SCFMode>::setupKernel(Options::GRID_PURPOSES gridAccuracy) {
  _kernel = std::make_shared<Kernel<SCFMode>>(_act, _env, _settings, _lrscf.size(), gridAccuracy);
  if (SCFMode == RESTRICTED && _settings.triplet) {
    printf("  Setting up triplet kernel!\n\n");
    _ukernel = std::make_shared<Kernel<UNRESTRICTED>>(_act, _env, _settings, _lrscf.size(), gridAccuracy);
  }
}

template<Options::SCF_MODES SCFMode>
SigmaCalculator ResponseLambda<SCFMode>::getTDASigma() {
  return _TDA;
}

template<Options::SCF_MODES SCFMode>
SigmaCalculator ResponseLambda<SCFMode>::getRPASigma() {
  return _RPA;
}

template<Options::SCF_MODES SCFMode>
NonlinearSigmaCalculator ResponseLambda<SCFMode>::getRightCC2Sigma() {
  return _rightCC2;
}

template<Options::SCF_MODES SCFMode>
NonlinearSigmaCalculator ResponseLambda<SCFMode>::getLeftCC2Sigma() {
  return _leftCC2;
}

template<Options::SCF_MODES SCFMode>
bool ResponseLambda<SCFMode>::usesKernel() {
  return _usesKernel;
}

template<Options::SCF_MODES SCFMode>
bool ResponseLambda<SCFMode>::usesDoubleHybrid() {
  return _usesDoubleHybrid;
}

template<Options::SCF_MODES SCFMode>
double ResponseLambda<SCFMode>::getDHRatio() {
  return _doubleHybridCorrRatio;
}

template<Options::SCF_MODES SCFMode>
bool ResponseLambda<SCFMode>::usesExchange() {
  return _usesExchange;
}

template<Options::SCF_MODES SCFMode>
bool ResponseLambda<SCFMode>::usesLRExchange() {
  return _usesLRExchange;
}

template<Options::SCF_MODES SCFMode>
bool ResponseLambda<SCFMode>::usesEO() {
  return _usesEO;
}

template class ResponseLambda<Options::SCF_MODES::RESTRICTED>;
template class ResponseLambda<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
