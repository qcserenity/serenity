/**
 * @file RIIntegrals.cpp
 *
 * @date Feb 14, 2020
 * @author Niklas Niemeyer
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "postHF/LRSCF/Tools/RIIntegrals.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "geometry/Geometry.h"
#include "integrals/RI_J_IntegralController.h"
#include "integrals/looper/TwoElecThreeCenterCalculator.h"
#include "integrals/looper/TwoElecThreeCenterIntLooper.h"
#include "misc/WarningTracker.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

/* Include Std and External Headers */

namespace Serenity {

template<Options::SCF_MODES SCFMode>
RIIntegrals<SCFMode>::RIIntegrals(std::shared_ptr<SystemController> sys, LIBINT_OPERATOR op, double mu, bool calcJia,
                                  unsigned pStart, unsigned pEnd, double nafThresh, std::shared_ptr<Geometry> geo,
                                  Options::DENS_FITS densFitCache) {
  LRSCFTaskSettings settings;
  settings.nafThresh = nafThresh;
  settings.densFitCache = densFitCache;
  (*this) =
      RIIntegrals<SCFMode>(std::make_shared<LRSCFController<SCFMode>>(sys, settings), op, mu, calcJia, pStart, pEnd, geo);
}

template<Options::SCF_MODES SCFMode>
RIIntegrals<SCFMode>::RIIntegrals(std::shared_ptr<LRSCFController<SCFMode>> lrscf, LIBINT_OPERATOR op, double mu,
                                  bool calcJia, unsigned pStart, unsigned pEnd, std::shared_ptr<Geometry> geo)
  : _lrscf(lrscf),
    _op(op),
    _mu(mu),
    _no(lrscf->getNOccupied()),
    _nv(lrscf->getNVirtual()),
    _pStart(pStart),
    _pEnd(pEnd),
    _calcJia(calcJia),
    _calcJpq(_pStart != 0 || _pEnd != 0),
    _basContr(lrscf->getBasisController(Options::BASIS_PURPOSES::DEFAULT)),
    _geo(geo),
    _fullyMOCached(false) {
  this->printInfo();
  this->calculateIntegrals();
}

template<Options::SCF_MODES SCFMode>
void RIIntegrals<SCFMode>::printInfo() {
  std::string cap;
  if (_op == LIBINT_OPERATOR::coulomb) {
    cap = "RI Integrals - Coulomb";
  }
  else if (_op == LIBINT_OPERATOR::erf_coulomb) {
    cap = "RI Integrals - erf-Coulomb";
  }
  else {
    throw SerenityError("Operator for RI integrals not yet supported.");
  }
  printBigCaption(cap);

  auto densFit = _lrscf.lock()->getLRSCFSettings().densFitCache;
  if (densFit == Options::DENS_FITS::NONE) {
    WarningTracker::printWarning("You have chosen NONE as density fitting for RIIntegrals. Will use RI instead", true);
    densFit = Options::DENS_FITS::RI;
  }
  if (densFit == Options::DENS_FITS::CD) {
    WarningTracker::printWarning("You have chosen CD as density fitting for RIIntegrals. Will use ACD instead", true);
    densFit = Options::DENS_FITS::ACD;
  }

  // New custom auxbasis controller dependent on the given geometry.
  if (_geo) {
    _auxBasContr = AtomCenteredBasisControllerFactory::produce(
        _geo, _lrscf.lock()->getSysSettings().basis.basisLibPath,
        _lrscf.lock()->getSysSettings().basis.makeSphericalBasis, false, 999999999,
        _lrscf.lock()->getSys()->getAuxBasisController(Options::AUX_BASIS_PURPOSES::CORRELATION, densFit)->getBasisString());
  }
  else {
    _auxBasContr = _lrscf.lock()->getSys()->getAuxBasisController(Options::AUX_BASIS_PURPOSES::CORRELATION, densFit);
  }

  printf("  Auxiliary Basis Set         : %15s\n\n", _auxBasContr->getBasisString().c_str());

  _nb = _basContr->getNBasisFunctions();
  _nxb = _auxBasContr->getNBasisFunctions();
  _nx = _nxb;

  // Info about orbitals and basis functions.
  bool isAlpha = true;
  for_spin(_no, _nv) {
    if (isAlpha) {
      printf("  Occupied orbitals (alpha)   : %15i\n", _no_spin);
      printf("  Virtual orbitals  (alpha)   : %15i\n", _nv_spin);
    }
    else {
      printf("  Occupied orbitals (beta)    : %15i\n", _no_spin);
      printf("  Virtual orbitals  (beta)    : %15i\n", _nv_spin);
    }
    isAlpha = false;
  };
  printf("  Basis functions             : %15lu\n", _nb);
  printf("  Auxiliary basis functions   : %15lu\n\n", _nxb);

  // Memory demand.
  double memDemand = 0.0;
  for_spin(_no, _nv) {
    memDemand += _nxb * _no_spin * _no_spin;
    if (_calcJia) {
      memDemand += _nxb * _no_spin * _nv_spin;
    }
    if (_calcJpq) {
      memDemand += (_pEnd - _pStart) * (_nv_spin + _no_spin) * _no_spin;
    }
  };
  memDemand *= sizeof(double);
  auto mem = MemoryManager::getInstance();
  size_t freeMem = 0.8 * MemoryManager::getInstance()->getAvailableSystemMemory();
  double gbFactor = 1e-9;
  printf("  Available memory            : %12.3f GB\n", gbFactor * freeMem);
  printf("  Minimal memory demand       : %12.3f GB\n\n", gbFactor * memDemand);
  if (freeMem < memDemand) {
    throw SerenityError("Sorry, RI-integral caching needs more memory than your machine provides.");
  }
}

template<Options::SCF_MODES SCFMode>
void RIIntegrals<SCFMode>::calculateIntegrals() {
  Timings::takeTime("RIIntegrals -        MO Integrals");

  RI_J_IntegralController riints(_basContr, _auxBasContr, nullptr, _op, _mu);
  _integrals = std::make_shared<TwoElecThreeCenterCalculator>(_op, _mu, _basContr, _basContr, _auxBasContr,
                                                              _basContr->getPrescreeningThreshold());
  _auxTrafo = std::make_shared<Eigen::MatrixXd>(riints.getInverseMSqrt());

  _metric = std::make_shared<Eigen::MatrixXd>(riints.getMetric());

  _Jij = std::make_shared<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>();
  _Jia = std::make_shared<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>();
  _Jpq = std::make_shared<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>();

  auto& Jij = *_Jij;
  auto& Jia = *_Jia;
  auto& Jpq = *_Jpq;

  auto C = _lrscf.lock()->getCoefficients();

  for_spin(Jij, Jia, Jpq, _no, _nv) {
    Jij_spin = Eigen::MatrixXd((long)_no_spin * _no_spin, _nxb);
    if (_calcJia) {
      Jia_spin = Eigen::MatrixXd((long)_nv_spin * _no_spin, _nxb);
    }
    if (_calcJpq) {
      Jpq_spin = Eigen::MatrixXd((long)(_no_spin + _nv_spin) * (_pEnd - _pStart), _nxb);
    }
  };

  auto distribute = [&](Eigen::Map<Eigen::MatrixXd> AO, size_t P, unsigned) {
    for_spin(C, Jij, Jia, Jpq, _no, _nv) {
      Eigen::Map<Eigen::MatrixXd> oo(Jij_spin.col(P).data(), _no_spin, _no_spin);
      oo = C_spin.middleCols(0, _no_spin).transpose() * AO * C_spin.middleCols(0, _no_spin);

      if (_calcJia) {
        Eigen::Map<Eigen::MatrixXd> vo(Jia_spin.col(P).data(), _nv_spin, _no_spin);
        vo = C_spin.middleCols(_no_spin, _nv_spin).transpose() * AO * C_spin.middleCols(0, _no_spin);
      }
      if (_calcJpq) {
        Eigen::Map<Eigen::MatrixXd> pq(Jpq_spin.col(P).data(), (_no_spin + _nv_spin), (_pEnd - _pStart));
        pq = C_spin.middleCols(0, (_no_spin + _nv_spin)).transpose() * AO * C_spin.middleCols(_pStart, (_pEnd - _pStart));
      }
    };
  };
  _integrals->loop(distribute);

  for_spin(Jij, Jia, Jpq, _no, _nv) {
    // P -> Q transformation.
    Jij_spin *= (*_auxTrafo);

    if (_calcJia) {
      // Avoid large temporary object.
      for (size_t i = 0; i < _no_spin; ++i) {
        Jia_spin.middleRows(i * _nv_spin, _nv_spin) *= (*_auxTrafo);
      }
    }

    if (_calcJpq) {
      Jpq_spin *= (*_auxTrafo);
    }
  };

  // NAF approximation.
  if (_lrscf.lock()->getLRSCFSettings().nafThresh != 0) {
    this->applyNAFApproximation();
  }

  _fullyMOCached = true;
  Timings::timeTaken("RIIntegrals -        MO Integrals");
}

template<Options::SCF_MODES SCFMode>
void RIIntegrals<SCFMode>::cacheAOIntegrals() {
  Timings::takeTime("RIIntegrals -          AO Caching");
  _integrals->cacheIntegrals();
  Timings::timeTaken("RIIntegrals -          AO Caching");
}

template<Options::SCF_MODES SCFMode>
void RIIntegrals<SCFMode>::applyNAFApproximation() {
  auto& settings = _lrscf.lock()->getLRSCFSettings();

  Timings::takeTime("RIIntegrals -   NAF Approximation");
  printf("\n    NAF (%-5.1e):\n\n", settings.nafThresh);
  double _scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;

  Eigen::MatrixXd naf = Eigen::MatrixXd::Zero(_nxb, _nxb);
  auto& Jij = *_Jij;
  auto& Jia = *_Jia;
  auto& Jpq = *_Jpq;
  for_spin(Jij, Jia) {
    naf.noalias() += _scfFactor * Jij_spin.transpose() * Jij_spin;
    if (_calcJia) {
      naf.noalias() += _scfFactor * 2.0 * Jia_spin.transpose() * Jia_spin;
    }
  };

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(naf);
  unsigned toKeep = 0;
  double thresh = settings.nafThresh;
  if (!_Jia) {
    thresh *= 1e-4;
    printf("    Modifying NAF threshold since only (ij|Q) integrals are available.\n");
    printf("    New NAF threshold: %-5.1e\n", thresh);
  }
  for (size_t Q = 0; Q < _nxb; ++Q) {
    if (eigensolver.eigenvalues()(Q) > thresh) {
      toKeep = _nxb - Q;
      break;
    }
  }

  printf("    Retaining %4i/%4i aux. functions\n", toKeep, (int)_nxb);
  printf("    -> removed %2.0f %%.\n\n", (1 - (double)toKeep / _nxb) * 100);

  Eigen::MatrixXd nafTrafo = eigensolver.eigenvectors().rightCols(toKeep);

  // The full auxiliary transformation now is P -> Q -> ~Q.
  (*_auxTrafo) *= (nafTrafo);

  // Set transformed aux basis dimension (now {~Q} due to NAF).
  _nx = _auxTrafo->cols();

  // Apply NAF transformation Q -> ~Q.
  for_spin(Jij, Jia, Jpq, _no, _nv) {
    Jij_spin *= nafTrafo;
    if (_calcJia) {
      for (size_t i = 0; i < _no_spin; ++i) {
        Jia_spin.block(i * _nv_spin, 0, _nv_spin, _nx) = (Jia_spin.block(i * _nv_spin, 0, _nv_spin, _nxb) * nafTrafo).eval();
      }
      Jia_spin.conservativeResize(Eigen::NoChange, _nx);
    }
    if (_calcJpq) {
      Jpq_spin *= nafTrafo;
    }
  };
  Timings::timeTaken("RIIntegrals -   NAF Approximation");
}

template<Options::SCF_MODES SCFMode>
void RIIntegrals<SCFMode>::clearMOCache() {
  if (_Jia) {
    auto Jia = *_Jia;
    for_spin(Jia) {
      Jia_spin.resize(0, 0);
    };
  }
  if (_Jij) {
    auto Jij = *_Jij;
    for_spin(Jij) {
      Jij_spin.resize(0, 0);
    };
  }
  if (_Jpq) {
    auto Jpq = *_Jpq;
    for_spin(Jpq) {
      Jpq_spin.resize(0, 0);
    };
  }
  _fullyMOCached = false;
}

template<Options::SCF_MODES SCFMode>
void RIIntegrals<SCFMode>::clearAOCache() {
  _integrals->clearCache();
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, unsigned> RIIntegrals<SCFMode>::getNOccupied() {
  return _no;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, unsigned> RIIntegrals<SCFMode>::getNVirtual() {
  return _nv;
}

template<Options::SCF_MODES SCFMode>
unsigned RIIntegrals<SCFMode>::getNBasisFunctions() {
  return _nb;
}

template<Options::SCF_MODES SCFMode>
unsigned RIIntegrals<SCFMode>::getNTransformedAuxBasisFunctions() {
  return _nx;
}

template<Options::SCF_MODES SCFMode>
unsigned RIIntegrals<SCFMode>::getNAuxBasisFunctions() {
  return _nxb;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> RIIntegrals<SCFMode>::getJijPtr() {
  return _Jij;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> RIIntegrals<SCFMode>::getJiaPtr() {
  return _Jia;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> RIIntegrals<SCFMode>::getJpqPtr() {
  return _Jpq;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Eigen::MatrixXd> RIIntegrals<SCFMode>::getAuxTrafoPtr() {
  return _auxTrafo;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Eigen::MatrixXd> RIIntegrals<SCFMode>::getAuxMetricPtr() {
  return _metric;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<TwoElecThreeCenterCalculator> RIIntegrals<SCFMode>::getIntegralPtr() {
  return _integrals;
}

template<Options::SCF_MODES SCFMode>
bool RIIntegrals<SCFMode>::isFullyMOCached() {
  return _fullyMOCached;
}

template<Options::SCF_MODES SCFMode>
int RIIntegrals<SCFMode>::getPStart() {
  return _pStart;
}

template<Options::SCF_MODES SCFMode>
int RIIntegrals<SCFMode>::getPEnd() {
  return _pEnd;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Geometry> RIIntegrals<SCFMode>::getGeo() {
  return _geo;
}

template<Options::SCF_MODES SCFMode>
void RIIntegrals<SCFMode>::setGeo(std::shared_ptr<Geometry> geo) {
  _geo = geo;
}

template class RIIntegrals<Options::SCF_MODES::RESTRICTED>;
template class RIIntegrals<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
