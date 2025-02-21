/**
 * @file   BasisController.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   30. Juli 2015, 17:19
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
#include "basis/BasisController.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "misc/SerenityError.h"

/* Include Std and External Headers */

namespace Serenity {

BasisController::BasisController(const std::string& basisString)
  : _basis(nullptr), _shellPairList(nullptr), _RIPrescreeningFactors(nullptr), _basisString(basisString) {
}

void BasisController::produceBasis() {
  assert(!_basis);
  _basis = produceBasisFunctionVector();
  assert(_basis);
  for (auto& basFunc : *_basis) {
    basFunc->addSensitiveObject(this->_self);
  }
  /*
   * For a nice handling of the basis some further information is very helpful which is
   * precalculated and stored here. In particular, some fiddling around with indices is necessary,
   * because a BasisFunction is in fact a complete shell of basis functions, i.e. the Basis contains
   * less entries than the actual number of basis functions. The quantities calculated here make
   * it possible to always convert an index in the vector of shells to an actual index for the
   * basis functions (which e.g. determines the size of matrices).
   */
  const unsigned int nBasRed = _basis->size();
  _extendedIndex.resize(nBasRed);
  _extendedIndexCart.resize(nBasRed);
  _extendedIndexSpherical.resize(nBasRed);
  _nBasisFunctions = 0;
  _nBasisFunctionsCart = 0;
  _nBasisFunctionsSpherical = 0;
  _maxAngularMomentum = 0;
  /*
   * Determine _nBasisFunctions and store the relationship
   * between extended->reduced indices.
   */
  for (unsigned int i = 0; i < nBasRed; ++i) {
    const unsigned int angMom = (*_basis)[i]->getAngularMomentum();
    _extendedIndex[i] = _nBasisFunctions;
    _extendedIndexCart[i] = _nBasisFunctionsCart;
    _extendedIndexSpherical[i] = _nBasisFunctionsSpherical;
    _nBasisFunctions += (*_basis)[i]->getNContracted();
    _nBasisFunctionsCart += N_SHELL_CART[angMom];
    _nBasisFunctionsSpherical += N_SHELL_SPH[angMom];
    /*
     * Check for maxAM
     */
    if (angMom > _maxAngularMomentum) {
      _maxAngularMomentum = angMom;
    }
  }
  /*
   * Precalculate the backrelation reduced->extended indices
   */
  _reducedIndex.reserve(_nBasisFunctionsCart);
  for (unsigned int i = 0, iRed = 0; iRed < nBasRed; ++iRed) {
    for (unsigned int j = 0; j < (*_basis)[iRed]->getNContracted(); ++j) {
      _reducedIndex[i + j] = iRed;
    }
    i += (*_basis)[iRed]->getNContracted();
  }
  _pureCartesian = (_nBasisFunctions == _nBasisFunctionsCart);
  _pureSpherical = (_nBasisFunctions == _nBasisFunctionsSpherical);
  _maxPrim = 1;
  for (const auto& shell : *_basis)
    _maxPrim = std::max(_maxPrim, shell->getNPrimitives());
  postConstruction();
}

void BasisController::calculateRIPrescreeningFactors() {
  _RIPrescreeningFactors = std::make_shared<std::vector<ShellPairData>>();
  // intialize libint
  auto& libint = Libint::getInstance();
  libint.initialize_plain(LIBINT_OPERATOR::coulomb, 2, std::numeric_limits<double>::epsilon() / 1e4, 10,
                          this->getMaxNumberOfPrimitives());
  Eigen::MatrixXd integrals;
  const unsigned int nShells = (*_basis).size();
  bool normAux = !(this->isAtomicCholesky());
  for (unsigned int i = 0; i < nShells; ++i) {
    const auto& shellI = *(*_basis)[i];
    // calculate integrals
    if (libint.compute(LIBINT_OPERATOR::coulomb, 0, shellI, shellI, integrals, normAux)) {
      (*_RIPrescreeningFactors).push_back(ShellPairData(i, i, sqrt(integrals.maxCoeff())));
    } /* if (prescreen) */
  }   /* i/shellI */
  // finalize libint
  libint.finalize(LIBINT_OPERATOR::coulomb, 0, 2);
}

void BasisController::createShellPairData() {
  unsigned nShells = (*_basis).size();
  _schwarzParams = Eigen::MatrixXd::Zero(nShells, nShells);
  _shellPairList = std::make_shared<std::vector<ShellPairData>>();
  // intialize libint
  auto& libint = Libint::getInstance();
  libint.initialize_plain(LIBINT_OPERATOR::coulomb, 4, std::numeric_limits<double>::epsilon() / 1e4, 10,
                          this->getMaxNumberOfPrimitives());

  // loops over shells
  Eigen::MatrixXd integrals;
  for (unsigned int i = 0; i < nShells; ++i) {
    const auto& shellI = *(*_basis)[i];
    for (unsigned int j = 0; j <= i; ++j) {
      const auto& shellJ = *(*_basis)[j];
      // calculate integrals
      if (libint.compute(LIBINT_OPERATOR::coulomb, 0, shellI, shellJ, shellI, shellJ, integrals)) {
        double schwarz = std::sqrt(integrals.maxCoeff());
        (*_shellPairList).push_back(ShellPairData(i, j, schwarz));
        _schwarzParams(i, j) = schwarz;
        _schwarzParams(j, i) = schwarz;
      } /* if (prescreen) */
    }   /* j/shellJ */
  }     /* i/shellI */
  // finalize libint
  libint.finalize(LIBINT_OPERATOR::coulomb, 0, 4);
  // sort the list
  std::sort((*_shellPairList).begin(), (*_shellPairList).end());
  std::reverse((*_shellPairList).begin(), (*_shellPairList).end());
}

const SparseMap& BasisController::getFunctionToShellMap() {
  if (!_functionToShellMap) {
    std::vector<Eigen::Triplet<int>> triplets;
    unsigned int rowCounter = 0;
    for (unsigned int iShell = 0; iShell < _basis->size(); ++iShell) {
      const auto& shell = (*_basis)[iShell];
      const unsigned int nCont = shell->getNContracted();
      for (unsigned int i = 0; i < nCont; ++i) {
        triplets.push_back(Eigen::Triplet<int>(rowCounter, iShell));
        ++rowCounter;
      } // for i
    }   // for iShell
    _functionToShellMap = std::make_shared<SparseMap>(getNBasisFunctions(), _basis->size());
    _functionToShellMap->setFromTriplets(triplets.begin(), triplets.end());
  } // if !_functionToShellMap
  return *_functionToShellMap;
}

unsigned int BasisController::getMaxNumberOfPrimitives() {
  return _maxPrim;
}

void BasisController::notify() {
  this->notifyObjects();
  this->_shellPairList = nullptr;
  this->_RIPrescreeningFactors = nullptr;
}

size_t BasisController::getReducedNBasisFunctions() {
  if (!_basis)
    produceBasis();
  return _basis->size();
}

Eigen::VectorXd BasisController::shellWiseAbsMax(const Eigen::VectorXd& target) {
  unsigned ns = _basis->size();
  auto basis = *_basis;
  Eigen::VectorXd result = Eigen::VectorXd::Zero(ns);
  for (unsigned i = 0; i < ns; ++i) {
    unsigned ni = basis[i]->size();
    unsigned iStart = this->extendedIndex(i);
    result(i) = target.segment(iStart, ni).array().abs().maxCoeff();
  } // for i
  return result;
}

double BasisController::getPrescreeningThreshold() {
  if (!_basis) {
    produceBasis();
  }
  return 1e-8 / (3.0 * _nBasisFunctionsCart);
}

const Eigen::MatrixXd& BasisController::getSchwarzParams(LIBINT_OPERATOR op, double mu) {
  if (op == LIBINT_OPERATOR::coulomb) {
    return _schwarzParams;
  }
  else if (op == LIBINT_OPERATOR::erf_coulomb) {
    // Represent mu as long. I know that this is kind of hacky and
    // I don't like it either, but it turned out that passing a
    // functional enum here makes things unnecessarily complicated.
    long int muInt = (long int)(mu * 1e7 + 0.5);
    if (_schwarzParamsErf.count(muInt) > 0) {
      return _schwarzParamsErf[muInt];
    }
    else {
      unsigned nShells = (*_basis).size();
      Eigen::MatrixXd schwarzParams = Eigen::MatrixXd::Zero(nShells, nShells);
      auto& libint = Libint::getInstance();
      libint.initialize(LIBINT_OPERATOR::erf_coulomb, 0, 4, std::vector<std::shared_ptr<Atom>>(0), mu,
                        std::numeric_limits<double>::epsilon() / 1e4, 10, this->getMaxNumberOfPrimitives());
      Eigen::MatrixXd integrals;
      for (auto& shellPair : (*_shellPairList)) {
        const auto& shellI = *(*_basis)[shellPair.bf1];
        const auto& shellJ = *(*_basis)[shellPair.bf2];
        if (libint.compute(LIBINT_OPERATOR::erf_coulomb, 0, shellI, shellJ, shellI, shellJ, integrals)) {
          double schwarz = std::sqrt(integrals.maxCoeff());
          schwarzParams(shellPair.bf1, shellPair.bf2) = schwarz;
          schwarzParams(shellPair.bf2, shellPair.bf1) = schwarz;
        } /* if (prescreen) */
      }
      libint.finalize(LIBINT_OPERATOR::erf_coulomb, 0, 4);
      _schwarzParamsErf[muInt] = schwarzParams;
      return _schwarzParamsErf[muInt];
    }
  }
  else {
    throw SerenityError("No Schwarz-prescreening parameters available for this operator.");
  }
}

bool BasisController::isAtomicCholesky() {
  return (this->_basisString.substr(0, 4) == "ACD-" or this->_basisString.substr(0, 5) == "ACCD-") ? true : false;
}

std::shared_ptr<std::vector<ShellPairData>> BasisController::getDOIPrescreeningFactors() {
  if (!_doiPrescreeningFactors) {
    unsigned nShells = (*_basis).size();
    _doiPrescreeningFactors = std::make_shared<std::vector<ShellPairData>>();
    // intialize libint
    auto& libint = Libint::getInstance();
    libint.initialize_plain(LIBINT_OPERATOR::delta, 4, std::numeric_limits<double>::epsilon() / 1e4, 10,
                            this->getMaxNumberOfPrimitives());

    // loops over shells
    Eigen::MatrixXd integrals;
    for (unsigned int i = 0; i < nShells; ++i) {
      const auto& shellI = *(*_basis)[i];
      for (unsigned int j = 0; j <= i; ++j) {
        const auto& shellJ = *(*_basis)[j];
        // calculate integrals
        if (libint.compute(LIBINT_OPERATOR::delta, 0, shellI, shellJ, shellI, shellJ, integrals)) {
          double integral = std::sqrt(integrals.maxCoeff());
          if (integral > 1e-12) {
            (*_doiPrescreeningFactors).push_back(ShellPairData(i, j, integral));
          }
        } /* if (prescreen) */
      }   /* j/shellJ */
    }     /* i/shellI */
    // finalize libint
    libint.finalize(LIBINT_OPERATOR::delta, 0, 4);
    // sort the list
    std::sort((*_doiPrescreeningFactors).begin(), (*_doiPrescreeningFactors).end());
    std::reverse((*_doiPrescreeningFactors).begin(), (*_doiPrescreeningFactors).end());
  }
  return _doiPrescreeningFactors;
}

BasisController::~BasisController() = default;

} /* namespace Serenity */
