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
#include "basis/Shell.h"
#include "geometry/Atom.h"
#include "integrals/Normalization.h"
#include "math/IntegerMaths.h"
#include "misc/Timing.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <omp.h>
#include <algorithm>
#include <array>

namespace Serenity {
using namespace std;

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
      it = IntegralType(libint2::Operator::overlap, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::kinetic, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::nuclear, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::emultipole1, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::emultipole2, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::delta, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::coulomb, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::cgtg, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::cgtg_x_coulomb, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::delcgtg2, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
      it = IntegralType(libint2::Operator::erf_coulomb, k, j);
      _engines[it] = std::vector<std::unique_ptr<libint2::Engine>>(nprocs);
    }
  }

  for (unsigned int k = 0; k < 3; ++k) {
    for (unsigned int j = 2; j < 5; ++j) {
      _keep[IntegralType(libint2::Operator::overlap, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::kinetic, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::nuclear, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::emultipole1, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::emultipole2, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::delta, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::coulomb, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::cgtg, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::cgtg_x_coulomb, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::delcgtg2, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::erf_coulomb, k, j)] = 0;
    }
  }

  _lock.unlock();
}

void Libint::initialize(libint2::Operator op, const unsigned int deriv, const unsigned int nCenter,
                        const std::vector<std::pair<double, std::array<double, 3>>> pointCharges, double mu) {
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
    if (op == libint2::Operator::nuclear) {
      if (pointCharges.size() > 0) {
#pragma omp parallel for schedule(static, 1)
        for (unsigned int i = 0; i < (unsigned int)evector.size(); ++i) {
          evector[i]->set_params(pointCharges);
        }
      }
    }
    else if (op == libint2::Operator::erf_coulomb || op == libint2::Operator::erfc_coulomb) {
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

  if (op == libint2::Operator::coulomb) {
    // Exception 1: if coulomb, the correct number of centers needs to be set here,
    evector[0].reset(new libint2::Engine(op, N_PRIM_MAX, am, deriv, std::numeric_limits<double>::epsilon(),
                                         libint2::operator_traits<libint2::Operator::coulomb>::default_params(), braket));
  }
  else if (op == libint2::Operator::nuclear) {
    // Exception 2: if nuclear, the nuclear charges (or other point charges) need to be added
    evector[0].reset(new libint2::Engine(op, N_PRIM_MAX, am, deriv, std::numeric_limits<double>::epsilon()));
    if (pointCharges.size() > 0) {
      evector[0]->set_params(pointCharges);
    }
  }
  else if (op == libint2::Operator::erf_coulomb) {
    // Exception 3: if erf_coulomb, the range separation parameter has to be included
    evector[0].reset(new libint2::Engine(op, N_PRIM_MAX, am, deriv, std::numeric_limits<double>::epsilon(),
                                         libint2::operator_traits<libint2::Operator::erf_coulomb>::default_params(), braket));
    evector[0]->set_params(mu);
  }
  else if (op == libint2::Operator::erfc_coulomb) {
    // Exception 3: if erfc_coulomb, the range separation parameter has to be included
    evector[0].reset(new libint2::Engine(op, N_PRIM_MAX, am, deriv, std::numeric_limits<double>::epsilon(),
                                         libint2::operator_traits<libint2::Operator::erfc_coulomb>::default_params(), braket));
    evector[0]->set_params(mu);
  }
  else {
    // Default
    evector[0].reset(new libint2::Engine(op, N_PRIM_MAX, am, deriv, std::numeric_limits<double>::epsilon()));
  }
#pragma omp parallel for schedule(static, 1)
  for (unsigned int i = 1; i < (unsigned int)evector.size(); ++i) {
    evector[i] = std::make_unique<libint2::Engine>(*evector[0]);
  }
  Timings::timeTaken("Tech. - Libint Initializations");
}

void Libint::initialize(libint2::Operator op, const unsigned int deriv, const unsigned int nCenter,
                        const std::vector<std::shared_ptr<Atom>>& atoms, double mu) {
  if (deriv > 1)
    throw SerenityError("Libint2 is only configured for 1st derivatives.");
  std::vector<std::pair<double, std::array<double, 3>>> q;
  for (unsigned int j = 0; j < atoms.size(); ++j) {
    q.push_back(
        {static_cast<double>(atoms[j]->getEffectiveCharge()), {{atoms[j]->getX(), atoms[j]->getY(), atoms[j]->getZ()}}});
  }
  this->initialize(op, deriv, nCenter, q, mu);
}

void Libint::initialize(libint2::Operator op, const unsigned int deriv, const unsigned int nCenter, const Point multipoleOrigin) {
  if (deriv > 1)
    throw SerenityError("Libint2 is only configured for 1st derivatives.");
  assert((op == libint2::Operator::emultipole1 || op == libint2::Operator::emultipole2) &&
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

  evector[0].reset(new libint2::Engine(op, N_PRIM_MAX, am, deriv, std::numeric_limits<double>::epsilon()));
  evector[0]->set_params(origin);
#pragma omp parallel for schedule(static, 1)
  for (unsigned int i = 1; i < (unsigned int)evector.size(); ++i) {
    evector[i] = std::make_unique<libint2::Engine>(*evector[0]);
  }
  Timings::timeTaken("Tech. - Libint Initializations");
}

void Libint::finalize(libint2::Operator op, unsigned int deriv, unsigned int nCenter) {
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

void Libint::keepEngines(libint2::Operator op, const unsigned int deriv, const unsigned int nCenter) {
  IntegralType it(op, deriv, nCenter);
  _keep[it] += 1;
}

void Libint::freeEngines(libint2::Operator op, const unsigned int deriv, const unsigned int nCenter) {
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
      _keep[IntegralType(libint2::Operator::overlap, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::kinetic, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::nuclear, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::emultipole1, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::emultipole2, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::delta, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::coulomb, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::cgtg, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::cgtg_x_coulomb, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::delcgtg2, k, j)] = 0;
      _keep[IntegralType(libint2::Operator::erf_coulomb, k, j)] = 0;
    }
  }
}

bool Libint::compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     Eigen::MatrixXd& ints) {
  IntegralType it(op, deriv, 2);
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
        Normalization::normalizeShell(ints.col(i), a.contr[0].l, b.contr[0].l);
      }
    }
    return true;
  }
}

bool Libint::compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     const libint2::Shell& c, Eigen::MatrixXd& ints) {
  IntegralType it(op, deriv, 3);

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
        Normalization::normalizeShell(ints.col(i), a.contr[0].l, b.contr[0].l, c.contr[0].l);
      }
    }
    return true;
  }
}

bool Libint::compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
                     const libint2::Shell& c, const libint2::Shell& d, Eigen::MatrixXd& ints) {
  IntegralType it(op, deriv, 4);

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

Eigen::MatrixXd Libint::compute1eInts(libint2::Operator op, std::shared_ptr<BasisController> basis,
                                      const std::vector<std::shared_ptr<Atom>>& atoms) {
  const auto& bas = basis->getBasis();
  const unsigned int nBFs = basis->getNBasisFunctions();
  Eigen::MatrixXd results(nBFs, nBFs);
  results.setZero();
  this->initialize(op, 0, 2, atoms);
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
        if (this->compute(op, 0, *bas[i], *bas[j], ints)) {
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

Eigen::MatrixXd Libint::compute1eInts(libint2::Operator op, std::shared_ptr<BasisController> basis,
                                      const std::vector<std::pair<double, std::array<double, 3>>> pointCharges) {
  const auto& bas = basis->getBasis();
  const unsigned int nBFs = basis->getNBasisFunctions();
  Eigen::MatrixXd results(nBFs, nBFs);
  results.setZero();
  this->initialize(op, 0, 2, pointCharges, 0.0);
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
        if (this->compute(op, 0, *bas[i], *bas[j], ints)) {
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

Eigen::MatrixXd Libint::compute1eInts(libint2::Operator op, std::shared_ptr<BasisController> basis1,
                                      std::shared_ptr<BasisController> basis2,
                                      const std::vector<std::shared_ptr<Atom>>& atoms) {
  const auto& bas1 = basis1->getBasis();
  const auto& bas2 = basis2->getBasis();
  const unsigned int nBF1s = basis1->getNBasisFunctions();
  const unsigned int nBF2s = basis2->getNBasisFunctions();
  Eigen::MatrixXd results(nBF2s, nBF1s);
  results.setZero();
  this->initialize(op, 0, 2, atoms);
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
        if (this->compute(op, 0, *bas1[i], *bas2[j], ints)) {
          Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(), nJ, nI);
          results.block(offJ, offI, nJ, nI) = tmp;
        }
      }
    }
  } /* END OpenMP parallel */
  this->finalize(op, 0, 2);
  return results;
}

Eigen::MatrixXd Libint::compute1eInts(libint2::Operator op, std::shared_ptr<BasisController> basis1,
                                      std::shared_ptr<BasisController> basis2,
                                      const std::vector<std::pair<double, std::array<double, 3>>> pointCharges) {
  const auto& bas1 = basis1->getBasis();
  const auto& bas2 = basis2->getBasis();
  const unsigned int nBF1s = basis1->getNBasisFunctions();
  const unsigned int nBF2s = basis2->getNBasisFunctions();
  Eigen::MatrixXd results(nBF1s, nBF2s);
  results.setZero();
  this->initialize(op, 0, 2, pointCharges, 0.0);
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
        if (this->compute(op, 0, *bas1[i], *bas2[j], ints)) {
          Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(), nJ, nI);
          results.block(offJ, offI, nJ, nI) = tmp;
        }
      }
    }
  } /* END OpenMP parallel */
  this->finalize(op, 0, 2);
  return results;
}

std::vector<std::unique_ptr<libint2::Engine>>& Libint::getFourCenterEngines(libint2::Operator op) {
  return _engines.at(IntegralType(op, 0, 4));
}

} /* namespace Serenity */
