/**
 * @file   ERIPotential.cpp
 *
 * @date   Nov 24, 2016
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
#include "potentials/ERIPotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/CDIntegralController.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "misc/Timing.h"
#include "potentials/CDExchangePotential.h"
#include "potentials/CDHFPotential.h"
#include "potentials/CoulombPotential.h"
#include "potentials/ExchangePotential.h"
#include "potentials/HFPotential.h"
#include "potentials/LRXPotential.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <algorithm>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ERIPotential<SCFMode>::ERIPotential(std::shared_ptr<SystemController> systemController,
                                    std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double xRatio,
                                    const double prescreeningThreshold, double prescreeningIncrementStart,
                                    double prescreeningIncrementEnd, unsigned int incrementSteps,
                                    bool externalSplitting, double lrxRatio, double mu, bool clear4CenterCache)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _systemController(systemController),
    _xRatio(xRatio),
    _lrxRatio(lrxRatio),
    _mu(mu),
    _dMatController(dMat),
    _fullpotential(nullptr),
    _fullXpotential(nullptr),
    _outOfDate(true) {
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);

  _fullpotential = std::make_unique<FockMatrix<SCFMode>>(FockMatrix<SCFMode>(this->_basis));
  auto& temp = *_fullpotential;
  for_spin(temp) {
    temp_spin.setZero();
  };
  _fullXpotential = std::make_unique<FockMatrix<SCFMode>>(FockMatrix<SCFMode>(this->_basis));
  auto& temp2 = *_fullXpotential;
  for_spin(temp2) {
    temp2_spin.setZero();
  };
  switch (systemController->getSettings().basis.densityFitting) {
    case Options::DENS_FITS::NONE:
    case Options::DENS_FITS::RI:
      if (systemController->getSettings().basis.densityFitting == Options::DENS_FITS::RI &&
          systemController->getSettings().method != Options::ELECTRONIC_STRUCTURE_THEORIES::HF && externalSplitting) {
        _coulomb = std::make_shared<CoulombPotential<SCFMode>>(
            systemController, dMat,
            RI_J_IntegralControllerFactory::getInstance().produce(
                systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                systemController->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)),
            prescreeningThreshold, prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps);
        // This is for the case that RI should be used with exact exchange
        if (_xRatio != 0.0)
          _exchange = std::make_shared<ExchangePotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                                   prescreeningIncrementStart, prescreeningIncrementEnd,
                                                                   incrementSteps, clear4CenterCache);
      }
      else {
        _hf = std::make_shared<HFPotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                     prescreeningIncrementStart, prescreeningIncrementEnd,
                                                     incrementSteps, clear4CenterCache);
      }
      if (_lrxRatio != 0.0) {
        _lrexchange = std::make_shared<LRXPotential<SCFMode>>(systemController, _dMatController, _lrxRatio,
                                                              prescreeningThreshold, prescreeningIncrementStart,
                                                              prescreeningIncrementEnd, incrementSteps, _mu);
      }
      break;
    case Options::DENS_FITS::CD:
      _hf = std::make_shared<CDHFPotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                     prescreeningIncrementStart, prescreeningIncrementEnd);
      if (_lrxRatio != 0.0) {
        WarningTracker::printWarning(
            "No LRExchange implemented for full CD. Using unfitted Potential. Try ACD/ACCD for accelerated LRExchange.", true);
        _lrexchange = std::make_shared<LRXPotential<SCFMode>>(systemController, _dMatController, _lrxRatio,
                                                              prescreeningThreshold, prescreeningIncrementStart,
                                                              prescreeningIncrementEnd, incrementSteps, _mu);
      }
      break;
    case Options::DENS_FITS::ACD:
      if (_xRatio != 0.0) {
        if (systemController->getCDIntegralController()->getACDVectors(
                systemController->getBasisController(),
                systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_CHOLESKY))) {
          // Integrals sum_{P} (munu|P)(P|Q)^{-1/2} are stored in memory
          _hf = std::make_shared<CDHFPotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                         prescreeningIncrementStart, prescreeningIncrementEnd);
        }
        else {
          _coulomb = std::make_shared<CoulombPotential<SCFMode>>(
              systemController, dMat,
              RI_J_IntegralControllerFactory::getInstance().produce(
                  systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                  systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_CHOLESKY)),
              prescreeningThreshold, prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps);
          _exchange = std::make_shared<CDExchangePotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold);
        }
      }
      else {
        _coulomb = std::make_shared<CoulombPotential<SCFMode>>(
            systemController, dMat,
            RI_J_IntegralControllerFactory::getInstance().produce(
                systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_CHOLESKY)),
            prescreeningThreshold, prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps);
      }
      if (_lrxRatio != 0.0) {
        _lrexchange = std::make_shared<CDExchangePotential<SCFMode>>(
            systemController, dMat, _lrxRatio, prescreeningThreshold, LIBINT_OPERATOR::erf_coulomb, _mu);
      }
      break;
    case Options::DENS_FITS::ACCD:
      if (_xRatio != 0.0) {
        if (systemController->getCDIntegralController()->getACDVectors(
                systemController->getBasisController(),
                systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY))) {
          // Integrals sum_{P} (munu|P)(P|Q)^{-1/2} are stored in memory
          _hf = std::make_shared<CDHFPotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                         prescreeningIncrementStart, prescreeningIncrementEnd);
        }
        else {
          _coulomb = std::make_shared<CoulombPotential<SCFMode>>(
              systemController, dMat,
              RI_J_IntegralControllerFactory::getInstance().produce(
                  systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                  systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY)),
              prescreeningThreshold, prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps);
          _exchange = std::make_shared<CDExchangePotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold);
        }
      }
      else {
        _coulomb = std::make_shared<CoulombPotential<SCFMode>>(
            systemController, dMat,
            RI_J_IntegralControllerFactory::getInstance().produce(
                systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY)),
            prescreeningThreshold, prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps);
      }
      if (_lrxRatio != 0.0) {
        _lrexchange = std::make_shared<CDExchangePotential<SCFMode>>(
            systemController, dMat, _lrxRatio, prescreeningThreshold, LIBINT_OPERATOR::erf_coulomb, _mu);
      }
      break;
  }
  _screening = prescreeningIncrementStart;
};

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& ERIPotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -   Coul./XC Pot.");
  if (_outOfDate) {
    _fullpotential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_fullpotential;
    for_spin(pot) {
      pot_spin.setZero();
    };
    if (_hf) {
      *_fullpotential += _hf->getMatrix();
    }
    else {
      assert(_coulomb);
      *_fullpotential += _coulomb->getMatrix();
      if (_exchange) {
        *_fullpotential += _exchange->getMatrix();
      }
    }
    if (_lrexchange)
      *_fullpotential += _lrexchange->getMatrix();
    _outOfDate = false;
  }
  Timings::timeTaken("Active System -   Coul./XC Pot.");
  return *_fullpotential;
}

template<Options::SCF_MODES SCFMode>
double ERIPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (_outOfDate)
    this->getMatrix();
  Timings::takeTime("Active System -   Coul./XC Pot.");
  auto& pot = *_fullpotential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += 0.5 * pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("Active System -   Coul./XC Pot.");
  return energy;
};

template<Options::SCF_MODES SCFMode>
double ERIPotential<SCFMode>::getXEnergy(const DensityMatrix<SCFMode>& P) {
  double energy = 0.0;
  if (_hf) {
    energy += _hf->getXEnergy(P);
  }
  else {
    if (_exchange)
      energy += _exchange->getEnergy(P);
  }
  if (_lrexchange)
    energy += _lrexchange->getEnergy(P);
  return energy;
};

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd ERIPotential<SCFMode>::getGeomGradients() {
  auto systemController = _systemController.lock();
  auto atoms = systemController->getAtoms();
  unsigned int nAtoms = atoms.size();

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(nAtoms, 3);
  if (_hf) {
    grad += _hf->getGeomGradients();
  }
  else {
    grad += _coulomb->getGeomGradients();
    if (_exchange) {
      grad += _exchange->getGeomGradients();
    }
  }
  if (_lrexchange)
    grad += _lrexchange->getGeomGradients();
  return grad;
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& ERIPotential<SCFMode>::getXPotential() {
  if (!_fullXpotential) {
    _fullXpotential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& exchange = *_fullXpotential;
    if (_hf) {
      exchange = _hf->getXPotential();
    }
    else {
      if (_exchange)
        exchange = _exchange->getMatrix();
    }
    if (_lrexchange)
      exchange += _lrexchange->getMatrix();
  }
  return *_fullXpotential;
}

template class ERIPotential<Options::SCF_MODES::RESTRICTED>;
template class ERIPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
