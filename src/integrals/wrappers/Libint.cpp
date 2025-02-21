/**
 * @file Libint.cpp
 * @author: Thomas Dresselhaus, Jan Unsleber
 *
 * @date 28. Juli 2013, 14:12
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
#include "integrals/wrappers/Libint.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h" //Loop shells.
#include "basis/BasisController.h"
#include "basis/Shell.h"
#include "geometry/Atom.h"
#include "integrals/Normalization.h"
#include "math/IntegerMaths.h"
#include "math/Matrix.h" //Engine matrix.
#include "misc/Timing.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <libint2.hpp>
#pragma GCC diagnostic pop
#include <omp.h>
#include <algorithm>
#include <array>

namespace Serenity {

/**
 * @brief Comparison between IntegralTypes in oder to make them usable as keys in std::maps.
 * @param lhs Left IntegralType.
 * @param rhs Right IntegralType.
 * @return Returns true if the left IntegralType is 'smaller' than the rhs IntegralType.
 */
bool operator<(const IntegralType& lhs, const IntegralType& rhs) {
  if (static_cast<int>(lhs.op) != static_cast<int>(rhs.op)) {
    return static_cast<int>(lhs.op) < static_cast<int>(rhs.op);
  }
  else {
    if (lhs.deriv != rhs.deriv) {
      return lhs.deriv < rhs.deriv;
    }
    else {
      return lhs.nCenter < rhs.nCenter;
    }
  }
}

static_assert(AM_MAX <= LIBINT2_MAX_AM_eri,
              "Libint is compiled to a lower angular momentum than may occur in this program (4-center ERIs)");
static_assert(AM_MAX <= LIBINT2_MAX_AM_2eri,
              "Libint is compiled to a lower angular momentum than may occur in this program (2-center ERIs)");
static_assert(AM_MAX <= LIBINT2_MAX_AM_3eri,
              "Libint is compiled to a lower angular momentum than may occur in this program (3-center ERIs)");

Libint::Libint() {
  _lock.lock();
  std::vector<Matrix<std::vector<std::unique_ptr<libint2::Engine>>>> mats;
  unsigned int nprocs = omp_get_max_threads();
  for (unsigned int k = 0; k < 2; ++k) {
    for (unsigned int j = 2; j < 5; ++j) {
      IntegralType it;
      it = IntegralType(LIBINT_OPERATOR::overlap, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::kinetic, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::nuclear, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::emultipole1, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::emultipole2, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::delta, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::coulomb, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::cgtg, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::cgtg_x_coulomb, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::delcgtg2, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::erf_coulomb, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(LIBINT_OPERATOR::erfc_coulomb, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
    }
  }

  for (unsigned int k = 0; k < 3; ++k) {
    for (unsigned int j = 2; j < 5; ++j) {
      _keep[IntegralType(LIBINT_OPERATOR::overlap, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::kinetic, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::nuclear, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::emultipole1, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::emultipole2, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::delta, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::coulomb, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::cgtg, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::cgtg_x_coulomb, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::delcgtg2, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::erf_coulomb, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::erfc_coulomb, k, j)] = 0;
    }
  }

  _lock.unlock();
}

long double Libint::getFinalPrecision(double precision, double maxD, double nPrim, unsigned int nCenters) {
  // If we assume that the error introduced by each primitive function is adds up, we have to consider
  // the total number of primitive functions occurring in a contracted integral.
  double nPrimPowN = 1.0;
  for (unsigned int center = 0; center < nCenters; ++center) {
    nPrimPowN /= (double)nPrim;
  }
  // Just to be really sure that this is safe we divide by 1e2. This may be quite conservative.
  return std::min(precision * nPrimPowN, precision / maxD) / 1e2;
}

void Libint::initialize(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter,
                        const std::vector<std::pair<double, std::array<double, 3>>>& pointCharges, double mu,
                        double precision, double maxD, unsigned int maxNPrim) {
  const long double finalPrecision = getFinalPrecision(precision, maxD, maxNPrim, nCenter);
  libint2::Operator libintOp = resolveLibintOperator(op);
  if (deriv > 1)
    throw SerenityError("Libint2 is only configured for 1st derivatives.");
  Timings::takeTime("Tech. - Libint Initializations");
  if (!libint2::initialized())
    libint2::initialize();
  IntegralType it(op, deriv, nCenter);
  auto& evector = _engines.at(it);
  /*
   * Create unique engine for operator, derivative, number of centers and OMP thread
   */
  if (evector[0]) {
    if (op == LIBINT_OPERATOR::nuclear) {
      if (pointCharges.size() > 0) {
#pragma omp parallel for schedule(static, 1)
        for (unsigned int i = 0; i < (unsigned int)evector.size(); ++i) {
          evector[i]->set_params(pointCharges);
        }
      }
    }
    else if (op == LIBINT_OPERATOR::erf_coulomb || op == LIBINT_OPERATOR::erfc_coulomb) {
#pragma omp parallel for schedule(static, 1)
      for (unsigned int i = 0; i < (unsigned int)evector.size(); ++i) {
        evector[i]->set_params(mu);
      }
    }
    Timings::timeTaken("Tech. - Libint Initializations");
    return;
  }

  // Set braket
  libint2::BraKet braket = libint2::BraKet::xx_xx;
  if (nCenter == 3)
    braket = libint2::BraKet::xs_xx;
  if (nCenter == 2)
    braket = libint2::BraKet::xs_xs;

  const int am = (deriv > 0) ? AM_MAX - 1 : AM_MAX;

  if (op == LIBINT_OPERATOR::coulomb) {
    // Exception 1: if coulomb, the correct number of centers needs to be set here,
    evector[0].reset(new libint2::Engine(libintOp, N_PRIM_MAX, am, deriv, finalPrecision,
                                         libint2::operator_traits<libint2::Operator::coulomb>::default_params(), braket));
  }
  else if (op == LIBINT_OPERATOR::nuclear) {
    // Exception 2: if nuclear, the nuclear charges (or other point charges) need to be added
    evector[0].reset(new libint2::Engine(libintOp, N_PRIM_MAX, am, deriv, finalPrecision));
    if (pointCharges.size() > 0) {
      evector[0]->set_params(pointCharges);
    }
  }
  else if (op == LIBINT_OPERATOR::erf_coulomb) {
    // Exception 3: if erf_coulomb, the range separation parameter has to be included
    evector[0].reset(new libint2::Engine(libintOp, N_PRIM_MAX, am, deriv, finalPrecision,
                                         libint2::operator_traits<libint2::Operator::erf_coulomb>::default_params(), braket));
    evector[0]->set_params(mu);
  }
  else if (op == LIBINT_OPERATOR::erfc_coulomb) {
    // Exception 3: if erfc_coulomb, the range separation parameter has to be included
    evector[0].reset(new libint2::Engine(libintOp, N_PRIM_MAX, am, deriv, finalPrecision,
                                         libint2::operator_traits<libint2::Operator::erfc_coulomb>::default_params(), braket));
    evector[0]->set_params(mu);
  }
  else {
    // Default
    evector[0].reset(new libint2::Engine(libintOp, N_PRIM_MAX, am, deriv, finalPrecision));
  }
#pragma omp parallel for schedule(static, 1)
  for (unsigned int i = 1; i < (unsigned int)evector.size(); ++i) {
    evector[i] = std::make_unique<libint2::Engine>(*evector[0]);
  }
  Timings::timeTaken("Tech. - Libint Initializations");
}

void Libint::initialize(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter,
                        const std::vector<std::shared_ptr<Atom>>& atoms, double mu, double precision, double maxD,
                        unsigned int maxNPrim) {
  if (deriv > 1)
    throw SerenityError("Libint2 is only configured for 1st derivatives.");
  std::vector<std::pair<double, std::array<double, 3>>> q;
  for (unsigned int j = 0; j < atoms.size(); ++j) {
    q.push_back(
        {static_cast<double>(atoms[j]->getEffectiveCharge()), {{atoms[j]->getX(), atoms[j]->getY(), atoms[j]->getZ()}}});
  }
  this->initialize(op, deriv, nCenter, q, mu, precision, maxNPrim, maxD);
}

void Libint::initialize(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter,
                        const Point multipoleOrigin, double precision, double maxD, unsigned int maxNPrim) {
  const long double finalPrecision = getFinalPrecision(precision, maxD, maxNPrim, nCenter);
  if (deriv > 1)
    throw SerenityError("Libint2 is only configured for 1st derivatives.");
  assert((op == LIBINT_OPERATOR::emultipole1 || op == LIBINT_OPERATOR::emultipole2) &&
         "If you put a multipole origin in the Libint initializer you should also use a multipole operator.");

  Timings::takeTime("Tech. - Libint Initializations");
  if (!libint2::initialized())
    libint2::initialize();
  IntegralType it(op, deriv, nCenter);
  auto& evector = _engines.at(it);

  const int am = (deriv > 0) ? AM_MAX - 1 : AM_MAX;

  // Libint expects a std::array as parameter
  std::array<double, 3> origin = {multipoleOrigin.getX(), multipoleOrigin.getY(), multipoleOrigin.getZ()};
  if (evector[0]) {
#pragma omp parallel for schedule(static, 1)
    for (unsigned int i = 0; i < (unsigned int)evector.size(); ++i) {
      evector[i]->set_params(origin);
    }
    Timings::timeTaken("Tech. - Libint Initializations");
    return;
  }

  evector[0].reset(new libint2::Engine(resolveLibintOperator(op), N_PRIM_MAX, am, deriv, finalPrecision));
  evector[0]->set_params(origin);
#pragma omp parallel for schedule(static, 1)
  for (unsigned int i = 1; i < (unsigned int)evector.size(); ++i) {
    evector[i] = std::make_unique<libint2::Engine>(*evector[0]);
  }
  Timings::timeTaken("Tech. - Libint Initializations");
}

void Libint::finalize(LIBINT_OPERATOR op, unsigned int deriv, unsigned int nCenter) {
  IntegralType it(op, deriv, nCenter);
  if (_keep[it])
    return;
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    _engines.at(it)[i].reset(nullptr);
  }
  for (auto const& engineMat : _engines) {
    for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
      if (engineMat.second[i] != nullptr)
        return;
    }
  }
  if (libint2::initialized())
    libint2::finalize();
}

void Libint::keepEngines(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter) {
  IntegralType it(op, deriv, nCenter);
  _keep[it] += 1;
}

void Libint::freeEngines(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter) {
  IntegralType it(op, deriv, nCenter);
  if (_keep[it] > 0)
    _keep[it] -= 1;
  if (_keep[it] == 0)
    this->finalize(op, deriv, nCenter);
}

Libint::~Libint() {
  if (!libint2::initialized())
    libint2::finalize();
}

void Libint::clearAllEngines() {
  for (auto& engineMat : _engines) {
    for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
      engineMat.second[i].reset(nullptr);
    }
  }
  libint2::finalize();
  for (unsigned int k = 0; k < 3; ++k) {
    for (unsigned int j = 2; j < 5; ++j) {
      _keep[IntegralType(LIBINT_OPERATOR::overlap, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::kinetic, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::nuclear, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::emultipole1, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::emultipole2, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::delta, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::coulomb, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::cgtg, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::cgtg_x_coulomb, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::delcgtg2, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::erf_coulomb, k, j)] = 0;
      _keep[IntegralType(LIBINT_OPERATOR::erfc_coulomb, k, j)] = 0;
    }
  }
}

bool Libint::compute(LIBINT_OPERATOR op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     Eigen::MatrixXd& ints, bool normAux) {
  return this->compute(resolveLibintOperator(op), deriv, a, b, ints, normAux);
}

bool Libint::compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     Eigen::MatrixXd& ints, bool normAux) {
  IntegralType it(resolveLibintOperator(op), deriv, 2);
  (void)normAux;
  const int threadID = omp_get_thread_num();

  // sanity check
  assert(_engines.at(it)[threadID]);
  assert(a.contr[0].pure == b.contr[0].pure);

  // get the result dimensions
  const auto nInts = a.size() * b.size();
  const auto nIntTypes = _engines.at(it)[threadID]->nshellsets();

  // buffer for the raw libint2 data
  const auto& buf_vec = _engines.at(it)[threadID]->results();

  // compute the integrals
  _engines.at(it)[threadID]->compute(a, b);

  if (buf_vec[0] == nullptr) {
    // integrals screened out
    ints.resize(0, 0);
    return false;
  }
  else {
    // store and normalize integral data
    ints.resize(nInts, nIntTypes);
    for (unsigned int i = 0; i < nIntTypes; ++i) {
      for (unsigned int j = 0; j < nInts; ++j) {
        ints(j, i) = buf_vec[i][j];
      }
      if (!a.contr[0].pure) {
        if (normAux)
          Normalization::normalizeShell(ints.col(i), a.contr[0].l, b.contr[0].l);
      }
    }
    return true;
  }
}

bool Libint::compute(LIBINT_OPERATOR op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     const libint2::Shell& c, Eigen::MatrixXd& ints, bool normAux) {
  return this->compute(resolveLibintOperator(op), deriv, a, b, c, ints, normAux);
}

bool Libint::compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     const libint2::Shell& c, Eigen::MatrixXd& ints, bool normAux) {
  (void)normAux;
  IntegralType it(resolveLibintOperator(op), deriv, 3);

  const int threadID = omp_get_thread_num();

  // sanity check
  assert(_engines.at(it)[threadID]);
  assert(a.contr[0].pure == b.contr[0].pure);
  assert(a.contr[0].pure == c.contr[0].pure);

  // get the result dimensions
  const auto nInts = a.size() * b.size() * c.size();
  const auto nIntTypes = _engines.at(it)[threadID]->nshellsets();
  // buffer for the raw libint2 data
  const auto& buf_vec = _engines.at(it)[threadID]->results();

  // compute the integrals
  _engines.at(it)[threadID]->compute(a, b, c);

  if (buf_vec[0] == nullptr) {
    // integrals screened out
    ints.resize(0, 0);
    return false;
  }
  else {
    // store and normalize integral data
    ints.resize(nInts, nIntTypes);
    for (unsigned int i = 0; i < nIntTypes; ++i) {
      for (unsigned int j = 0; j < nInts; ++j) {
        ints(j, i) = buf_vec[i][j];
      }
      if (!a.contr[0].pure) {
        if (normAux) {
          Normalization::normalizeShell(ints.col(i), a.contr[0].l, b.contr[0].l, c.contr[0].l);
        }
        else {
          Normalization::normalizeShellNoAux(ints.col(i), a.contr[0].l, b.contr[0].l, c.contr[0].l);
        }
      }
    }
    return true;
  }
}

bool Libint::compute(LIBINT_OPERATOR op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     const libint2::Shell& c, const libint2::Shell& d, Eigen::MatrixXd& ints) {
  return this->compute(resolveLibintOperator(op), deriv, a, b, c, d, ints);
}

bool Libint::compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     const libint2::Shell& c, const libint2::Shell& d, Eigen::MatrixXd& ints) {
  IntegralType it(resolveLibintOperator(op), deriv, 4);

  const int threadID = omp_get_thread_num();

  // sanity check
  assert(_engines.at(it)[threadID]);
  assert(a.contr[0].pure == b.contr[0].pure);
  assert(a.contr[0].pure == c.contr[0].pure);
  assert(a.contr[0].pure == d.contr[0].pure);

  // get the result dimensions
  const auto nInts = a.size() * b.size() * c.size() * d.size();
  const auto nIntTypes = _engines.at(it)[threadID]->nshellsets();

  // buffer for the raw libint2 data
  const auto& buf_vec = _engines.at(it)[threadID]->results();

  // compute the integrals
  _engines.at(it)[threadID]->compute(a, b, c, d);

  if (buf_vec[0] == nullptr) {
    // integrals screened out
    ints.resize(0, 0);
    return false;
  }
  else {
    // store and normalize integral data
    ints.resize(nInts, nIntTypes);
    for (unsigned int i = 0; i < nIntTypes; ++i) {
      for (unsigned int j = 0; j < nInts; ++j) {
        ints(j, i) = buf_vec[i][j];
      }
      if (!a.contr[0].pure) {
        Normalization::normalizeShell(ints.col(i), a.contr[0].l, b.contr[0].l, c.contr[0].l, d.contr[0].l);
      }
    }
    return true;
  }
}

Eigen::MatrixXd Libint::compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis,
                                      const std::vector<std::shared_ptr<Atom>>& atoms, double precision, double maxD,
                                      unsigned int maxNPrim) {
  bool normAux = !(basis->isAtomicCholesky());
  libint2::Operator libintOp = resolveLibintOperator(op);
  const auto& bas = basis->getBasis();
  const unsigned int nBFs = basis->getNBasisFunctions();
  Eigen::MatrixXd results(nBFs, nBFs);
  results.setZero();
  this->initialize(op, 0, 2, atoms, 0.0, precision, maxD, maxNPrim);
#pragma omp parallel
  {
    Eigen::MatrixXd ints;
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < bas.size(); ++i) {
      unsigned int offI = basis->extendedIndex(i);
      const unsigned int nI = bas[i]->getNContracted();
      for (unsigned int j = 0; j <= i; ++j) {
        unsigned int offJ = basis->extendedIndex(j);
        const unsigned int nJ = bas[j]->getNContracted();
        if (this->compute(libintOp, 0, *bas[i], *bas[j], ints, normAux)) {
          Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(), nJ, nI);
          results.block(offJ, offI, nJ, nI) = tmp;
          if (i != j)
            results.block(offI, offJ, nI, nJ) = tmp.transpose();
        }
      }
    }
  } /* END OpenMP parallel */
  this->finalize(op, 0, 2);
  return results;
}

Eigen::MatrixXd Libint::compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis,
                                      const std::vector<std::pair<double, std::array<double, 3>>> pointCharges,
                                      double precision, double maxD, unsigned int maxNPrim,
                                      std::shared_ptr<std::vector<ShellPairData>> shellPairData) {
  if (!shellPairData) {
    shellPairData = std::make_shared<std::vector<ShellPairData>>();
    const unsigned int nShells = basis->getReducedNBasisFunctions();
    for (unsigned int iShell = 0; iShell < nShells; ++iShell) {
      for (unsigned int jShell = 0; jShell <= iShell; ++jShell) {
        shellPairData->emplace_back(iShell, jShell, 1.0);
      }
    }
  }

  bool normAux = !(basis->isAtomicCholesky());
  libint2::Operator libintOp = resolveLibintOperator(op);
  const auto& bas = basis->getBasis();
  const unsigned int nBFs = basis->getNBasisFunctions();
  Eigen::MatrixXd results(nBFs, nBFs);
  results.setZero();
  this->initialize(op, 0, 2, pointCharges, 0.0, precision, maxD, maxNPrim);
  const auto& pairData = *shellPairData;
#pragma omp parallel
  {
    Eigen::MatrixXd ints;
#pragma omp for schedule(dynamic, 1)
    for (unsigned int iShellPair = 0; iShellPair < pairData.size(); ++iShellPair) {
      const unsigned int i = pairData[iShellPair].bf1;
      const unsigned int j = pairData[iShellPair].bf2;
      const unsigned int offI = basis->extendedIndex(i);
      const unsigned int offJ = basis->extendedIndex(j);
      const unsigned int nI = bas[i]->getNContracted();
      const unsigned int nJ = bas[j]->getNContracted();
      if (this->compute(libintOp, 0, *bas[i], *bas[j], ints, normAux)) {
        Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(), nJ, nI);
        results.block(offJ, offI, nJ, nI) = tmp;
        if (i != j)
          results.block(offI, offJ, nI, nJ) = tmp.transpose();
      }
    }
  } /* END OpenMP parallel */
  this->finalize(op, 0, 2);
  return results;
}

Eigen::MatrixXd Libint::compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis1,
                                      std::shared_ptr<BasisController> basis2, const std::vector<std::shared_ptr<Atom>>& atoms,
                                      double mu, double precision, double maxD, unsigned int maxNPrim) {
  bool normAux = !(basis1->isAtomicCholesky());
  libint2::Operator libintOp = resolveLibintOperator(op);
  const auto& bas1 = basis1->getBasis();
  const auto& bas2 = basis2->getBasis();
  const unsigned int nBF1s = basis1->getNBasisFunctions();
  const unsigned int nBF2s = basis2->getNBasisFunctions();
  // TODO: The basis ordering here does not make any sense: Why first basis2 and then basis1?
  //      It may be more intuitive to clean this up.
  Eigen::MatrixXd results(nBF2s, nBF1s);
  results.setZero();
  this->initialize(op, 0, 2, atoms, mu, precision, maxD, maxNPrim);
#pragma omp parallel
  {
    Eigen::MatrixXd ints;
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < bas1.size(); ++i) {
      unsigned int offI = basis1->extendedIndex(i);
      const unsigned int nI = bas1[i]->getNContracted();
      for (unsigned int j = 0; j < bas2.size(); ++j) {
        unsigned int offJ = basis2->extendedIndex(j);
        const unsigned int nJ = bas2[j]->getNContracted();
        if (this->compute(libintOp, 0, *bas1[i], *bas2[j], ints, normAux)) {
          Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(), nJ, nI);
          results.block(offJ, offI, nJ, nI) = tmp;
        }
      }
    }
  } /* END OpenMP parallel */
  this->finalize(op, 0, 2);
  return results;
}

Eigen::MatrixXd Libint::compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis1,
                                      std::shared_ptr<BasisController> basis2,
                                      const std::vector<std::pair<double, std::array<double, 3>>> pointCharges,
                                      double precision, double maxD, unsigned int maxNPrim) {
  bool normAux = !(basis1->isAtomicCholesky());
  libint2::Operator libintOp = resolveLibintOperator(op);
  const auto& bas1 = basis1->getBasis();
  const auto& bas2 = basis2->getBasis();
  const unsigned int nBF1s = basis1->getNBasisFunctions();
  const unsigned int nBF2s = basis2->getNBasisFunctions();
  Eigen::MatrixXd results(nBF1s, nBF2s);
  results.setZero();
  this->initialize(op, 0, 2, pointCharges, 0.0, precision, maxD, maxNPrim);
#pragma omp parallel
  {
    Eigen::MatrixXd ints;
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < bas1.size(); ++i) {
      unsigned int offI = basis1->extendedIndex(i);
      const unsigned int nI = bas1[i]->getNContracted();
      for (unsigned int j = 0; j < bas2.size(); ++j) {
        unsigned int offJ = basis2->extendedIndex(j);
        const unsigned int nJ = bas2[j]->getNContracted();
        if (this->compute(libintOp, 0, *bas1[i], *bas2[j], ints, normAux)) {
          Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(), nJ, nI);
          results.block(offJ, offI, nJ, nI) = tmp;
        }
      }
    }
  } /* END OpenMP parallel */
  this->finalize(op, 0, 2);
  return results;
}

void Libint::initialize_plain(LIBINT_OPERATOR op, const unsigned int nCenter, double precision, double maxD,
                              unsigned int maxNPrim) {
  this->initialize(op, 0, nCenter, std::vector<std::shared_ptr<Atom>>(0), 0.0, precision, maxD, maxNPrim);
}

libint2::Operator Libint::resolveLibintOperator(LIBINT_OPERATOR op) {
  switch (op) {
    case LIBINT_OPERATOR::overlap:
      return libint2::Operator::overlap;
    case LIBINT_OPERATOR::kinetic:
      return libint2::Operator::kinetic;
    case LIBINT_OPERATOR::nuclear:
      return libint2::Operator::nuclear;
    case LIBINT_OPERATOR::emultipole1:
      return libint2::Operator::emultipole1;
    case LIBINT_OPERATOR::emultipole2:
      return libint2::Operator::emultipole2;
    case LIBINT_OPERATOR::delta:
      return libint2::Operator::delta;
    case LIBINT_OPERATOR::coulomb:
      return libint2::Operator::coulomb;
    case LIBINT_OPERATOR::cgtg:
      return libint2::Operator::cgtg;
    case LIBINT_OPERATOR::cgtg_x_coulomb:
      return libint2::Operator::cgtg_x_coulomb;
    case LIBINT_OPERATOR::delcgtg2:
      return libint2::Operator::delcgtg2;
    case LIBINT_OPERATOR::erf_coulomb:
      return libint2::Operator::erf_coulomb;
    case LIBINT_OPERATOR::erfc_coulomb:
      return libint2::Operator::erfc_coulomb;
      /*
       * No default in order to generate a compiler warning
       * in case of an incomplete switch.
       */
  }
  return libint2::Operator::overlap;
}

LIBINT_OPERATOR Libint::resolveLibintOperator(libint2::Operator op) {
  switch (op) {
    case libint2::Operator::overlap:
      return LIBINT_OPERATOR::overlap;
    case libint2::Operator::kinetic:
      return LIBINT_OPERATOR::kinetic;
    case libint2::Operator::nuclear:
      return LIBINT_OPERATOR::nuclear;
    case libint2::Operator::emultipole1:
      return LIBINT_OPERATOR::emultipole1;
    case libint2::Operator::emultipole2:
      return LIBINT_OPERATOR::emultipole2;
    case libint2::Operator::delta:
      return LIBINT_OPERATOR::delta;
    case libint2::Operator::coulomb:
      return LIBINT_OPERATOR::coulomb;
    case libint2::Operator::cgtg:
      return LIBINT_OPERATOR::cgtg;
    case libint2::Operator::cgtg_x_coulomb:
      return LIBINT_OPERATOR::cgtg_x_coulomb;
    case libint2::Operator::delcgtg2:
      return LIBINT_OPERATOR::delcgtg2;
    case libint2::Operator::erf_coulomb:
      return LIBINT_OPERATOR::erf_coulomb;
    case libint2::Operator::erfc_coulomb:
      return LIBINT_OPERATOR::erfc_coulomb;
    default:
      throw SerenityError("ERROR: Operator unknown to Serenity.");
  }
  return LIBINT_OPERATOR::overlap;
}

std::vector<std::unique_ptr<libint2::Engine>>& Libint::getFourCenterEngines(LIBINT_OPERATOR op) {
  return _engines.at(IntegralType(op, 0, 4));
}

unsigned int Libint::getNPrimMax() {
  return N_PRIM_MAX;
}
void Libint::initialize(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter,
                        const std::vector<std::pair<double, Point>>& pointCharges, double mu, double precision,
                        double maxD, unsigned int maxNPrim) {
  std::vector<std::pair<double, std::array<double, 3>>> convertedCharges;
  for (const auto& charge : pointCharges) {
    convertedCharges.emplace_back(charge.first,
                                  std::array<double, 3>{charge.second.getX(), charge.second.getY(), charge.second.getZ()});
  }
  this->initialize(op, deriv, nCenter, convertedCharges, mu, precision, maxD, maxNPrim);
}
Eigen::MatrixXd Libint::compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis,
                                      const std::vector<std::pair<double, Point>>& pointCharges, double precision, double maxD,
                                      unsigned int maxNPrim, std::shared_ptr<std::vector<ShellPairData>> shellPairData) {
  std::vector<std::pair<double, std::array<double, 3>>> convertedCharges;
  for (const auto& charge : pointCharges) {
    convertedCharges.emplace_back(charge.first,
                                  std::array<double, 3>{charge.second.getX(), charge.second.getY(), charge.second.getZ()});
  }
  return compute1eInts(op, basis, convertedCharges, precision, maxD, maxNPrim, shellPairData);
}

} /* namespace Serenity */
