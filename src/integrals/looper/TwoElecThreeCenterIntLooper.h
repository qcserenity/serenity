/**
 * @file TwoElecThreeCenterIntLooper.h
 *
 * @date Oct 28, 2016
 * @author: Jan Unsleber
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
#ifndef TWOELECTHREECENTERINTLOOPER_H_
#define TWOELECTHREECENTERINTLOOPER_H_
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "integrals/wrappers/Libint.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/**
 * @class TwoElecThreeCenterIntLooper TwoElecThreeCenterIntLooper.h
 * @brief A looper for 2e- 3 center integrals between two basis sets.
 */
class TwoElecThreeCenterIntLooper {
 public:
  /**
   * @brief Constructor.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param basis The basis.
   * @param auxbasis The auxiliary basis.
   * @param kRange The range for the K shells.
   *        The first shell index is 0, the second index should end
   *        at the number of shells at maximum.
   *        kRange.first and kRange.second are the boundaries for a loop.
   *        @code
   *         for (unsigned int K=kRange.first;K < kRange.second;++K){
   *         }
   *        @endcode
   * @param mu Range separation parameter if operator is erf_coulomb.
   */
  TwoElecThreeCenterIntLooper(LIBINT_OPERATOR op, const unsigned deriv, std::shared_ptr<BasisController> basis,
                              std::shared_ptr<BasisController> auxbasis, double prescreeningThreshold,
                              std::pair<unsigned int, unsigned int> kRange = {0, 0}, double mu = 0.0)
    : _op(op),
      _deriv(deriv),
      _basis(basis),
      _auxbasis(auxbasis),
      _prescreeningThreshold(prescreeningThreshold),
      _kRange(kRange),
      _mu(mu) {
  }
  /**
   * @brief Default destructor.
   */
  virtual ~TwoElecThreeCenterIntLooper() = default;
  /**
   * @brief Loops over 2e-3-center integrals.
   *
   * This function expects a lambda function as input.
   * The function will be called for each integral
   * (be aware of symmetry!).
   * The function has to expect four arguments:
   * @code
   *   void std::function<void (unsigned& i,
   *   unsigned& j,
   *   unsigned& K,
   *   Eigen::VectorXd& intValues)>
   * @endcode
   * Here i and j are the indices of the 'normal' basis while
   * K is the index of the auxiliary basis.
   * \f[ \langle K | ij \rangle \f]
   * All indices are unfolded (running over all contracted functions
   * in all shells.)
   * The routine uses the symmetry within i and j thus all function calls
   * should expect i<=j and account for the 'missing' calls.
   *
   * The function can then do anything with the given integrals,
   * e.g. sum up:
   * @code
   *  TwoElecThreeCenterIntLooper looper(op, deriv, basis, auxbasis);
   *  Eigen::VectorXd sum = 0.0;
   *  auto const sumInts = [&sum]
   *                       (unsigned&  i,
   *                        unsigned&  j,
   *                        unsigned&  K,
   *                        Eigen::VectorXd& intValues,
   *                        unsigned threadId) {
   *  sum[threadId] += intValues[0];
   *  if (i!=j) sum[threadId] += intValues[0];
   *  };
   *  looper.loop(sumInts);
   *  double threadSum = sum.sum();
   * @endcode
   *
   * For interaction terms i and j should be centered on one system and K
   * on the other, for the complete interaction also the switched version is
   * usually needed.
   *
   * @param distribute The function to use each integral, see above for extended description.
   * @param prescreen A function to use for more detailed integral prescreening.
   * @param maxD      Absolute maximum coefficient to be contracted with the integrals. Used to
   *                  adjust the precision in case of absurdly large coefficients.
   */
  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loop(Func distribute, PrescreenFunc prescreen, double maxD = 1) {
    // intialize libint
    auto& libint = Libint::getInstance();
    libint.initialize(_op, _deriv, 3, std::vector<std::shared_ptr<Atom>>(0), _mu, std::numeric_limits<double>::epsilon(),
                      maxD, std::max(_basis->getMaxNumberOfPrimitives(), _auxbasis->getMaxNumberOfPrimitives()));
    bool normAux = !(_auxbasis->isAtomicCholesky());
    auto& basis = _basis->getBasis();
    auto& auxbasis = _auxbasis->getBasis();
    const auto& shellPairs = _basis->getShellPairData();
    const auto& riPrescreeningFactors = _auxbasis->getRIPrescreeningFactors();
    if (_kRange.first == _auxbasis->getNBasisFunctions())
      return;
    if (_kRange.second == 0)
      _kRange.second = _auxbasis->getNBasisFunctions();

#ifdef _OPENMP
    std::vector<Eigen::MatrixXd> integrals(omp_get_max_threads());
#else
    std::vector<Eigen::MatrixXd> integrals(1);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int kIndex = _auxbasis->reducedIndex(_kRange.second - 1);
         kIndex >= (int)_auxbasis->reducedIndex(_kRange.first); --kIndex) {
      auto& k = (*riPrescreeningFactors)[kIndex];

#ifdef _OPENMP
      unsigned threadId = omp_get_thread_num();
#else
      unsigned threadId = 0;
#endif
      unsigned shellK = k.bf1;
      const auto& auxbasK = *auxbasis[shellK];
      unsigned nK = auxbasis[shellK]->getNContracted();
      // loops over shells
      for (auto& p : *shellPairs) {
        if (p.factor * k.factor < _prescreeningThreshold)
          break;

        const bool swap = (basis[p.bf1]->getAngularMomentum() < basis[p.bf2]->getAngularMomentum());

        unsigned i = (swap) ? p.bf2 : p.bf1;
        unsigned j = (swap) ? p.bf1 : p.bf2;

        const auto& basI = *basis[i];
        const auto& basJ = *basis[j];

        unsigned nI = basis[i]->getNContracted();
        unsigned nJ = basis[j]->getNContracted();
        /*
         * Optional advanced prescreening
         */
        if (prescreen(i, j, shellK, p.factor * k.factor)) {
          continue;
        }
        // calculate integrals
        if (libint.compute(_op, _deriv, auxbasK, basI, basJ, integrals[threadId], normAux)) {
          // angular momentum reordering
          if (swap) {
            if (integrals[threadId].cols() == 9) {
              integrals[threadId].col(3).swap(integrals[threadId].col(6));
              integrals[threadId].col(4).swap(integrals[threadId].col(7));
              integrals[threadId].col(5).swap(integrals[threadId].col(8));
            }
            else if (integrals[threadId].cols() > 9) {
              std::cout << "2nd Derivatives and higher not yet supported!" << std::endl;
              assert(false);
            }
          }
          // unpack and run
          for (unsigned int K = 0; K < nK; ++K) {
            unsigned kk = _auxbasis->extendedIndex(shellK) + K;
            if (kk >= _kRange.second)
              continue;
            if (kk < _kRange.first)
              continue;
            for (unsigned int I = 0; I < nI; ++I) {
              unsigned ii = _basis->extendedIndex(i) + I;
              for (unsigned int J = 0; J < nJ; ++J) {
                unsigned jj = _basis->extendedIndex(j) + J;
                if (swap) {
                  if (ii > jj)
                    continue;
                  Eigen::VectorXd set(integrals[threadId].row(K * nJ * nI + I * nJ + J));
                  distribute(jj, ii, kk, set, threadId);
                }
                else {
                  if (jj > ii)
                    continue;
                  Eigen::VectorXd set(integrals[threadId].row(K * nJ * nI + I * nJ + J));
                  distribute(ii, jj, kk, set, threadId);
                }
              } /* bfs of j -> J */
            }   /* bfs of i -> I */
          }     /* bfs of k -> K */
        }       /* if (prescreen) */
      }         /* p/shellpair */
    }           /* K/shellK */
    // finalize libint
    libint.finalize(_op, _deriv, 3);
  }

  /**
   * @brief Loops over 2e-3-center integrals.
   *        Cannot handle derivatives, however, is much faster at the cost of generality.
   *
   * @param distribute The function to use each integral, see above for extended description.
   * @param prescreen A function to use for more detailed integral prescreening.
   * @param maxD      Absolute maximum coefficient to be contracted with the integrals. Used to
   *                  adjust the precision in case of absurdly large coefficients.
   */
  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loopNoDerivative(Func distribute, PrescreenFunc prescreen, double maxD = 1) {
    bool normAux = !(_auxbasis->isAtomicCholesky());
    // intialize libint
    auto& libint = Libint::getInstance();
    libint.initialize(_op, 0, 3, std::vector<std::shared_ptr<Atom>>(0), _mu, std::numeric_limits<double>::epsilon(),
                      maxD, std::max(_basis->getMaxNumberOfPrimitives(), _auxbasis->getMaxNumberOfPrimitives()));
    auto& basis = _basis->getBasis();
    auto& auxbasis = _auxbasis->getBasis();
    const auto& shellPairs = _basis->getShellPairData();
    const auto& riPrescreeningFactors = _auxbasis->getRIPrescreeningFactors();
    if (_kRange.first == _auxbasis->getNBasisFunctions())
      return;
    if (_kRange.second == 0)
      _kRange.second = _auxbasis->getNBasisFunctions();

#ifdef _OPENMP
    std::vector<Eigen::MatrixXd> integrals(omp_get_max_threads());
#else
    std::vector<Eigen::MatrixXd> integrals(1);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int kIndex = _auxbasis->reducedIndex(_kRange.second - 1);
         kIndex >= (int)_auxbasis->reducedIndex(_kRange.first); --kIndex) {
      auto& k = (*riPrescreeningFactors)[kIndex];

#ifdef _OPENMP
      unsigned threadId = omp_get_thread_num();
#else
      unsigned threadId = 0;
#endif
      unsigned shellK = k.bf1;
      const auto& auxbasK = *auxbasis[shellK];
      unsigned nK = auxbasis[shellK]->getNContracted();
      // loops over shells
      for (auto& p : *shellPairs) {
        if (p.factor * k.factor < _prescreeningThreshold)
          break;

        const bool swap = (basis[p.bf1]->getAngularMomentum() < basis[p.bf2]->getAngularMomentum());

        unsigned i = (swap) ? p.bf2 : p.bf1;
        unsigned j = (swap) ? p.bf1 : p.bf2;

        const auto& basI = *basis[i];
        const auto& basJ = *basis[j];

        unsigned nI = basis[i]->getNContracted();
        unsigned nJ = basis[j]->getNContracted();
        /*
         * Optional advanced prescreening
         */
        if (prescreen(i, j, shellK, p.factor * k.factor)) {
          continue;
        }
        // calculate integrals
        if (libint.compute(_op, 0, auxbasK, basI, basJ, integrals[threadId], normAux)) {
          const double* intptr = integrals[threadId].data();
          for (unsigned int K = 0; K < nK; ++K) {
            unsigned kk = _auxbasis->extendedIndex(shellK) + K;
            if (kk >= _kRange.second || kk < _kRange.first) {
              intptr += nI * nJ;
              continue;
            }
            for (unsigned int I = 0; I < nI; ++I) {
              unsigned ii = _basis->extendedIndex(i) + I;
              for (unsigned int J = 0; J < nJ; ++J, ++intptr) {
                unsigned jj = _basis->extendedIndex(j) + J;
                if (swap) {
                  if (ii > jj) {
                    continue;
                  }
                  distribute(jj, ii, kk, *intptr, threadId);
                }
                else {
                  if (jj > ii) {
                    continue;
                  }
                  distribute(ii, jj, kk, *intptr, threadId);
                }
              } /* bfs of j -> J */
            }   /* bfs of i -> I */
          }     /* bfs of k -> K */
        }       /* if (prescreen) */
      }         /* p/shellpair */
    }           /* K/shellK */
    // finalize libint
    libint.finalize(_op, 0, 3);
  }

  /**
   * @brief Loops over 2e-3-center integrals.
   *        No prescreening is being done here.
   *
   * @param distribute The function to use each integral, see above for extended description.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loop(Func distribute, double maxD = 1) {
    auto prescreen = [](const unsigned, const unsigned, const unsigned, const double) { return false; };
    loop(distribute, prescreen, maxD);
  }

  /**
   * @brief Loops over 2e-3-center integrals.
   *        Cannot handle derivatives, however, is much faster at the cost of generality.
   *        No prescreening is being done here.
   *
   * @param distribute The function to use each integral, see above for extended description.
   * @param maxD  Maximum coefficient to be contracted with the integrals. Used to adjust
   *             the integral precision in case of absurdly large coefficients.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loopNoDerivative(Func distribute, double maxD = 1) {
    auto prescreen = [](const unsigned, const unsigned, const unsigned, const double) { return false; };
    loopNoDerivative(distribute, prescreen, maxD);
  }

 private:
  /// @brief The kernel/operator as libint enum.
  LIBINT_OPERATOR _op;
  /// @brief The derivative level.
  unsigned _deriv;
  /// @brief The basis.
  std::shared_ptr<BasisController> _basis;
  /// @brief The auxiliary basis.
  std::shared_ptr<BasisController> _auxbasis;
  /// @brief A threshold for the Schwartz conditions.
  double _prescreeningThreshold;
  /// @brief The range for the K shells.
  std::pair<unsigned int, unsigned int> _kRange;
  /// @brief Range separation parameter mu
  double _mu;
};
} /* namespace Serenity */
#endif /* TWOELECTHREECENTERINTLOOPER_H_ */
