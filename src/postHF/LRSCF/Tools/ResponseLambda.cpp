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
#include "postHF/LRSCF/Sigmavectors/CoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/EOSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/ExchangeSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/FockSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/GrimmeSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/KernelSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RI/RICoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RI/RIExchangeSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RICC2/XWFController.h"
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
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
      _usesKernel = true;
    }
  }
  for (const auto& sys : _env) {
    if (sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      _usesKernel = true;
    }
  }

  // Simplified TDDFT to be used
  if (_settings.grimme) {
    _grimme = std::make_unique<GrimmeSigmavector<SCFMode>>(_lrscf);
    _usesKernel = false;
  }

  // Exchange to be used
  for (unsigned I = 0; I < _lrscf.size(); ++I) {
    if (_lrscf[I]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
      _usesExchange = true;
    }
    else {
      auto func_enum = this->_lrscf[I]->getLRSCFSettings().func;
      auto func = resolveFunctional(func_enum);

      // Resolve functional.
      if (func_enum == CompositeFunctionals::XCFUNCTIONALS::NONE) {
        func = resolveFunctional(_lrscf[I]->getSysSettings().dft.functional);
      }
      _usesExchange = func.isHybrid() ? func.isHybrid() : _usesExchange;
      _usesLRExchange = func.isRSHybrid() ? func.isRSHybrid() : _usesLRExchange;
      _usesDoubleHybrid = func.isDoubleHybrid() ? func.isDoubleHybrid() : _usesDoubleHybrid;
    }
  }
  auto naddXCfunc = resolveFunctional(_settings.embedding.naddXCFunc);
  _usesExchange = naddXCfunc.isHybrid() ? naddXCfunc.isHybrid() : _usesExchange;
  _usesLRExchange = naddXCfunc.isRSHybrid() ? naddXCfunc.isRSHybrid() : _usesLRExchange;
  for (auto naddXCfunc : _settings.embedding.naddXCFuncList) {
    auto func = resolveFunctional(naddXCfunc);
    _usesExchange = func.isHybrid() ? func.isHybrid() : _usesExchange;
    _usesLRExchange = func.isRSHybrid() ? func.isRSHybrid() : _usesLRExchange;
  }

  if (_usesLRExchange && _settings.rpaScreening) {
    throw SerenityError("RS-Hybrid Functional not supported for BSE calculations!");
  }

  if (_usesDoubleHybrid) {
    printBigCaption("Double-Hybrid TDDFT");
    _doubleHybridCorrRatio = resolveFunctional(_lrscf[0]->getSysSettings().dft.functional).getHfCorrelRatio();
    // Set spin scaling factors from functional.
    _settings.sss = resolveFunctional(_lrscf[0]->getSysSettings().dft.functional).getssScaling();
    _settings.oss = resolveFunctional(_lrscf[0]->getSysSettings().dft.functional).getosScaling();
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
  }

  _densFitJ = _settings.densFitJ != Options::DENS_FITS::NONE;
  _densFitK = _settings.densFitK != Options::DENS_FITS::NONE;
  _densFitLRK = _settings.densFitLRK != Options::DENS_FITS::NONE;

  // Setup RI Integral cache for ADC(2) or CC2.
  bool isXWF = !(settings.method == Options::LR_METHOD::TDA || settings.method == Options::LR_METHOD::TDDFT);
  if (isXWF || _settings.rpaScreening || _usesDoubleHybrid) {
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
  if (isXWF || _usesDoubleHybrid) {
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
    printSubSectionTitle("Ground-state Calculation");
    for (auto& lrscf : _lrscf) {
      // In case of double-hybrids, we fake the method to be CIS(D) to get the
      // reduced sigma vector build which does not contain the CIS contributions.
      if (_usesDoubleHybrid) {
        auto oldmethod = _settings.method;
        _settings.method = Options::LR_METHOD::CISD;
        lrscf->initializeXWFController();
        _settings.method = oldmethod;
      }
      else {
        lrscf->initializeXWFController();
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

  if (!_settings.grimme && !isXWF) {
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
  }
  if (_settings.adaptivePrescreening) {
    printf("  Using adaptive prescreening thresholds.\n\n");
  }
}

template<Options::SCF_MODES SCFMode>
void ResponseLambda<SCFMode>::setupTDDFTLambdas() {
  // Sigmacalculator: (A+B)*b stored in [0,2,4,..], (A-B)*b stored in [1,3,5,..]
  _RPA = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::make_unique<std::vector<Eigen::MatrixXd>>(
        guessVectors.size(), Eigen::MatrixXd::Zero(guessVectors[0].rows(), guessVectors[0].cols()));

    std::vector<int> xc(guessVectors.size());
    std::vector<Eigen::MatrixXd> guessAPB(guessVectors.size() / 2);
    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      // (A+B)
      if (iSet % 2 == 0) {
        xc[iSet] = 1;
        guessAPB[iSet / 2] = guessVectors[iSet];
      }
      // (A-B)
      else {
        xc[iSet] = -1;
      }
    }

    D = std::make_unique<FockSigmavector<SCFMode>>(_lrscf, guessVectors);

    if (_densFitJ) {
      J = std::make_unique<RICoulombSigmavector<SCFMode>>(_lrscf, guessAPB);
    }
    else {
      J = std::make_unique<CoulombSigmavector<SCFMode>>(_lrscf, guessAPB);
    }

    if (_usesExchange || _usesLRExchange) {
      if (_densFitK || _densFitLRK) {
        DFK = std::make_unique<RIExchangeSigmavector<SCFMode>>(_lrscf, guessVectors, xc, _densFitK, _densFitLRK);
      }
      if (!_densFitK || !_densFitLRK) {
        K = std::make_unique<ExchangeSigmavector<SCFMode>>(_lrscf, guessVectors, xc, _densFitK, _densFitLRK);
      }
    }

    if (_usesKernel) {
      F = std::make_unique<KernelSigmavector<SCFMode>>(_lrscf, guessAPB, _kernel, _ukernel);
    }

    if (_usesEO) {
      EO = std::make_unique<EOSigmavector<SCFMode>>(_lrscf, guessVectors, _settings.embedding.levelShiftParameter,
                                                    _settings.embedding.embeddingMode);
    }

    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      (*sigma)[iSet] += D->getSigma()[iSet];
      if (iSet % 2 == 0) {
        if (!(SCFMode == RESTRICTED && _settings.triplet)) {
          (*sigma)[iSet] += 2.0 * _scfFactor * J->getSigma()[iSet / 2];
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
    }

    return sigma;
  }; /* RPA */

  // Sigmacalculator: sqrt(A-B)*(A+B)*sqrt(A-B)*b
  _TDDFT = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::make_unique<std::vector<Eigen::MatrixXd>>(
        guessVectors.size(), Eigen::MatrixXd::Zero(guessVectors[0].rows(), guessVectors[0].cols()));

    std::vector<Eigen::MatrixXd> transGuess(guessVectors.size());

    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      transGuess[iSet] = _diagonal.cwiseSqrt().asDiagonal() * guessVectors[iSet];
    }

    D = std::make_unique<FockSigmavector<SCFMode>>(_lrscf, transGuess);
    if (_densFitJ) {
      J = std::make_unique<RICoulombSigmavector<SCFMode>>(_lrscf, transGuess);
    }
    else {
      J = std::make_unique<CoulombSigmavector<SCFMode>>(_lrscf, transGuess);
    }

    F = std::make_unique<KernelSigmavector<SCFMode>>(_lrscf, transGuess, _kernel, _ukernel);

    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      (*sigma)[iSet] += D->getSigma()[iSet];
      if (!(SCFMode == RESTRICTED && _settings.triplet)) {
        (*sigma)[iSet] += 2 * _scfFactor * J->getSigma()[iSet];
      }
      if (_usesKernel) {
        double tripletFactor = (_settings.triplet ? 0.5 : 1.0);
        (*sigma)[iSet] += 2 * _scfFactor * tripletFactor * F->getSigma()[iSet];
      }
      (*sigma)[iSet] = (_diagonal.cwiseSqrt().asDiagonal() * (*sigma)[iSet]).eval();
    }

    return sigma;
  }; /* TDDFT */

  // Sigmacalculator: A*b
  _TDA = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::make_unique<std::vector<Eigen::MatrixXd>>(
        guessVectors.size(), Eigen::MatrixXd::Zero(guessVectors[0].rows(), guessVectors[0].cols()));

    D = std::make_unique<FockSigmavector<SCFMode>>(_lrscf, guessVectors);

    if (_densFitJ) {
      J = std::make_unique<RICoulombSigmavector<SCFMode>>(_lrscf, guessVectors);
    }
    else {
      J = std::make_unique<CoulombSigmavector<SCFMode>>(_lrscf, guessVectors);
    }

    if (_usesExchange || _usesLRExchange) {
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
      EO = std::make_unique<EOSigmavector<SCFMode>>(_lrscf, guessVectors, _settings.embedding.levelShiftParameter,
                                                    _settings.embedding.embeddingMode);
    }

    if (_grimme) {
      (*sigma) = _grimme->getSigmavectors(guessVectors, std::vector<int>(guessVectors.size(), 0));
    }

    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      (*sigma)[iSet] += D->getSigma()[iSet];
    }

    for (unsigned iSet = 0; iSet < guessVectors.size() && !_settings.grimme; ++iSet) {
      if (_settings.scfstab == Options::STABILITY_ANALYSIS::NONE || _settings.scfstab == Options::STABILITY_ANALYSIS::REAL) {
        double aplusbFactor = (_settings.scfstab == Options::STABILITY_ANALYSIS::REAL ? 2.0 : 1.0);
        if (!(SCFMode == RESTRICTED && _settings.triplet)) {
          if (J) {
            (*sigma)[iSet] += aplusbFactor * _scfFactor * J->getSigma()[iSet];
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
    }

    return sigma;
  }; /* TDA */

  if (_settings.grimme) {
    _TDA = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
      std::vector<int> pm(guessVectors.size(), 0);
      auto sigma = std::make_unique<std::vector<Eigen::MatrixXd>>(_grimme->getSigmavectors(guessVectors, pm));
      D = std::make_unique<FockSigmavector<SCFMode>>(_lrscf, guessVectors);
      for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
        (*sigma)[iSet] += D->getSigma()[iSet];
      }
      return sigma;
    };

    _RPA = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
      std::vector<int> pm(guessVectors.size());
      for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
        pm[iSet] = (iSet % 2 == 0) ? 1 : -1;
      }
      auto sigma = std::make_unique<std::vector<Eigen::MatrixXd>>(_grimme->getSigmavectors(guessVectors, pm));
      D = std::make_unique<FockSigmavector<SCFMode>>(_lrscf, guessVectors);
      for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
        (*sigma)[iSet] += D->getSigma()[iSet];
      }
      return sigma;
    };
  }
}

template<Options::SCF_MODES SCFMode>
void ResponseLambda<SCFMode>::setupCC2Lambdas() {
  // lambda function for right transformation CC2/CIS(Dinf)/ADC(2)
  _rightXWF = [&](Eigen::Ref<Eigen::MatrixXd> guessVectors, Eigen::VectorXd guessValues) {
    auto sigma = std::make_unique<Eigen::MatrixXd>(guessVectors.rows(), guessVectors.cols());
    for (unsigned iCol = 0; iCol < guessVectors.cols(); ++iCol) {
      sigma->col(iCol) = _lrscf[0]->getXWFController()->getRightXWFSigma(guessVectors.col(iCol), guessValues(iCol));
    }
    return sigma;
  };

  // lambda function for left transformation CC2/CIS(Dinf)/ADC(2)
  _leftXWF = [&](Eigen::Ref<Eigen::MatrixXd> guessVectors, Eigen::VectorXd guessValues) {
    auto sigma = std::make_unique<Eigen::MatrixXd>(guessVectors.rows(), guessVectors.cols());
    for (unsigned iCol = 0; iCol < guessVectors.cols(); ++iCol) {
      sigma->col(iCol) = _lrscf[0]->getXWFController()->getLeftXWFSigma(guessVectors.col(iCol), guessValues(iCol));
    }
    return sigma;
  };
}

template<Options::SCF_MODES SCFMode>
SigmaCalculator ResponseLambda<SCFMode>::setEigensolverMode(bool fockDiagonal) {
  // If TDA: ignore user input and solve symmetric problem.
  if (_settings.method != Options::LR_METHOD::TDDFT) {
    _settings.algorithm = Options::RESPONSE_ALGORITHM::SYMMETRIC;
    return _TDA;
  }
  else {
    // If not TDA and RPA requested, solve RPA.
    if (_settings.algorithm == Options::RESPONSE_ALGORITHM::SYMPLECTIC) {
      return _RPA;
    }
    // If not explicitly requested, solve symmetrized problem if possible and RPA problem otherwise.
    else {
      if (!fockDiagonal || _usesEO || _usesExchange || _usesLRExchange || _settings.partialResponseConstruction ||
          _settings.grimme) {
        // (A-B) is not diagonal, have to solve symplectic problem
        _settings.algorithm = Options::RESPONSE_ALGORITHM::SYMPLECTIC;
        return _RPA;
      }
      else {
        // (A-B) is diagonal, might also solve the symmetrized problem
        _settings.algorithm = Options::RESPONSE_ALGORITHM::SYMMETRIZED;
        return _TDDFT;
      }
    }
  }
}

template<Options::SCF_MODES SCFMode>
void ResponseLambda<SCFMode>::setupKernel(Options::GRID_PURPOSES gridFineness) {
  _kernel = std::make_shared<Kernel<SCFMode>>(_act, _env, _settings, gridFineness);
  if (SCFMode == RESTRICTED && _settings.triplet) {
    printf("  Setting up triplet kernel!\n\n");
    _ukernel = std::make_shared<Kernel<UNRESTRICTED>>(_act, _env, _settings, gridFineness);
  }
}

template<Options::SCF_MODES SCFMode>
SigmaCalculator ResponseLambda<SCFMode>::getTDASigma() {
  return _TDA;
}

template<Options::SCF_MODES SCFMode>
SigmaCalculator ResponseLambda<SCFMode>::getTDDFTSigma() {
  return _TDDFT;
}

template<Options::SCF_MODES SCFMode>
SigmaCalculator ResponseLambda<SCFMode>::getRPASigma() {
  return _RPA;
}

template<Options::SCF_MODES SCFMode>
XWFSigmaCalculator ResponseLambda<SCFMode>::getRightXWFSigma() {
  return _rightXWF;
}

template<Options::SCF_MODES SCFMode>
XWFSigmaCalculator ResponseLambda<SCFMode>::getLeftXWFSigma() {
  return _leftXWF;
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
