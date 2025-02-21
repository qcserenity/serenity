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
#include "potentials/CDHFPotential.h"
#include "potentials/ExchangePotential.h"
#include "potentials/HFPotential.h"
#include "potentials/LRXPotential.h"
#include "potentials/RICoulombPotential.h"
#include "potentials/RIExchangePotential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <algorithm>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ERIPotential<SCFMode>::ERIPotential(std::shared_ptr<SystemController> systemController,
                                    std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double xRatio,
                                    const double prescreeningThreshold, double prescreeningIncrementStart,
                                    double prescreeningIncrementEnd, unsigned int incrementSteps, double lrxRatio,
                                    double mu, bool clear4CenterCache)
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

  auto densFitJ = systemController->getSettings().basis.densFitJ;
  auto densFitK = systemController->getSettings().basis.densFitK;
  auto densFitLRK = systemController->getSettings().basis.densFitLRK;

  /*
   * First: Coulomb and HF exchange potential.
   */

  // HF: Ignore keywords and do not use and density fitting for either the Coulomb or the exchange potential.
  if (systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    _hf = std::make_shared<HFPotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                 prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps,
                                                 clear4CenterCache);
  }
  // DFT: Respect all keywords.
  else {
    if (densFitJ == Options::DENS_FITS::NONE and densFitK == Options::DENS_FITS::NONE) {
      _hf = std::make_shared<HFPotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                   prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps,
                                                   clear4CenterCache);
    }
    else if (densFitJ == Options::DENS_FITS::CD and densFitK == Options::DENS_FITS::CD) {
      _hf = std::make_shared<CDHFPotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                     prescreeningIncrementStart, prescreeningIncrementEnd);
    }
    else {
      // Coulomb.
      switch (densFitJ) {
        case Options::DENS_FITS::NONE:
          _coulomb = std::make_shared<HFPotential<SCFMode>>(systemController, dMat, 0.0, prescreeningThreshold,
                                                            prescreeningIncrementStart, prescreeningIncrementEnd,
                                                            incrementSteps, clear4CenterCache);
          break;
        case Options::DENS_FITS::CD:
          _coulomb = std::make_shared<CDHFPotential<SCFMode>>(systemController, dMat, 0.0, prescreeningThreshold,
                                                              prescreeningIncrementStart, prescreeningIncrementEnd);
          break;
        case Options::DENS_FITS::RI:
        case Options::DENS_FITS::ACD:
        case Options::DENS_FITS::ACCD:
          _coulomb = std::make_shared<RICoulombPotential<SCFMode>>(
              systemController, dMat,
              RI_J_IntegralControllerFactory::getInstance().produce(
                  systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                  systemController->getAuxBasisController(Options::AUX_BASIS_PURPOSES::COULOMB, densFitJ)),
              prescreeningThreshold, prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps);
          break;
      }

      // HF Exchange.
      if (_xRatio != 0.0) {
        switch (densFitK) {
          case Options::DENS_FITS::NONE:
            _exchange = std::make_shared<ExchangePotential<SCFMode>>(systemController, dMat, _xRatio, prescreeningThreshold,
                                                                     prescreeningIncrementStart, prescreeningIncrementEnd,
                                                                     incrementSteps, clear4CenterCache);
            break;
          case Options::DENS_FITS::CD:
            throw SerenityError(
                "Isolated full CD for the Exchange not implemented. Set densFitJ CD for complete CD treatment.");
            break;
          case Options::DENS_FITS::RI:
          case Options::DENS_FITS::ACD:
          case Options::DENS_FITS::ACCD:
            _exchange = std::make_shared<RIExchangePotential<SCFMode>>(
                systemController, dMat,
                RI_J_IntegralControllerFactory::getInstance().produce(
                    systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                    systemController->getAuxBasisController(Options::AUX_BASIS_PURPOSES::EXCHANGE, densFitK)),
                _xRatio, prescreeningThreshold);
            break;
        }
      }
    }
  }

  /*
   * Second: Long-range exchange potential.
   */
  if (_lrxRatio != 0.0) {
    switch (densFitLRK) {
      case Options::DENS_FITS::NONE:
        _lrexchange = std::make_shared<LRXPotential<SCFMode>>(systemController, _dMatController, _lrxRatio,
                                                              prescreeningThreshold, prescreeningIncrementStart,
                                                              prescreeningIncrementEnd, incrementSteps, _mu);
        break;
      case Options::DENS_FITS::CD:
        throw SerenityError("LR-Exchange for full CD not implemented. Chose a different setting for densFitLRK");
        break;
      case Options::DENS_FITS::RI:
      case Options::DENS_FITS::ACD:
      case Options::DENS_FITS::ACCD:
        _lrexchange = std::make_shared<RIExchangePotential<SCFMode>>(
            systemController, dMat,
            RI_J_IntegralControllerFactory::getInstance().produce(
                systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                systemController->getAuxBasisController(Options::AUX_BASIS_PURPOSES::LREXCHANGE, densFitLRK),
                LIBINT_OPERATOR::erf_coulomb, _mu),
            _lrxRatio, prescreeningThreshold, LIBINT_OPERATOR::erf_coulomb, _mu);
        break;
    }
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
