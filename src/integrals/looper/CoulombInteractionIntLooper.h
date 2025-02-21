/**
 * @file CoulombInteractionIntLooper.h
 *
 * @date Oct 29, 2016
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

#ifndef COULOMBINTERACTIONINTLOOPER_H_
#define COULOMBINTERACTIONINTLOOPER_H_

/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "integrals/wrappers/Libint.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
/* Include Std and External Headers */
#include <libint2/engine.h>
#pragma GCC diagnostic pop

namespace Serenity {
/**
 * @class CoulombInteractionIntLooper CoulombInteractionIntLooper.h
 * @brief A looper for 2e- 4 center integrals between two basis sets.
 */
class CoulombInteractionIntLooper {
 public:
  /**
   * @brief Constructor.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param basisOne The basis of one system.
   * @param basisTwo The basis of the second interacting system.
   */
  CoulombInteractionIntLooper(LIBINT_OPERATOR op, const unsigned deriv, std::shared_ptr<BasisController> basisOne,
                              std::shared_ptr<BasisController> basisTwo, double prescreeningThreshold)
    : _op(op), _deriv(deriv), _basisOne(basisOne), _basisTwo(basisTwo), _prescreeningThreshold(prescreeningThreshold) {
  }

  /**
   * @brief Default destructor.
   */
  virtual ~CoulombInteractionIntLooper() = default;

  /**
   * @brief Loops over 2e- 4-center interaction integrals.
   *
   * This function expects a lambda function as input.
   * The function will be called for each integral
   * (be aware of symmetry!).
   * The function has to expect four arguments:
   * @code
   *   void std::function<void (const unsigned int& i,
   *   const unsigned int& j,
   *   const unsigned int&  a,
   *   const unsigned int&  b,
   *   Eigen::VectorXd& intValues
   *   const unsigned int&  threadID)>
   * @endcode
   * Here i and j are the indices of the first basis while
   * a and b are the indces of the second basis.
   * \f[ \langle ij | ab \rangle \f]
   * All indices are unfolded (running over all contracted functions
   * in all shells.)
   * The routine uses the symmetry within i and j thus all function calls
   * should expect i<=j and account for the 'missing' calls.
   * The same holds for a and b.
   *
   * The function can than do anything with the given integrals,
   * e.g. sum up:
   * @code
   *  TwoElecFourCenterInteractionIntLooper looper(op, deriv, basis1, basis2);
   *  Eigen::VectorXd sum(nThreads);
   *  sum.setZero();
   *  auto const sumInts = [&sum]
   *                       (const unsigned int&  i,
   *                        const unsigned int&  j,
   *                        const unsigned int&  a,
   *                        const unsigned int&  b,
   *                        Eigen::VectorXd& intValues,
   *                        const unsigned int&  threadID) {
   *  sum[threadID] += intValues(0);
   *  if (i!=j) sum[threadID] += intValues[0];
   *  if (a!=b) sum[threadID] += intValues[0];
   *  if (a!=b && i!=j) sum[threadID] += intValues[0];
   *  };
   *  looper.loop(sumInts);
   *  // sum over thread contributions
   * @endcode
   *
   *
   * @param loopEvalFunction The function to use each integral, see above for extended description.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loop(Func loopEvalFunction) {
    loop(loopEvalFunction, [](unsigned int, unsigned int, unsigned int, unsigned int, double) { return false; });
  }

  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loop(Func loopEvalFunction, PrescreenFunc prescreenFunc) {
    // intialize libint

    auto& libint = Libint::getInstance();
    libint.keepEngines(_op, _deriv, 4);
    libint.initialize(_op, _deriv, 4);

    const auto& shellPairsB1 = _basisOne->getShellPairData();
    const auto& basis1 = _basisOne->getBasis();
    const auto& shellPairsB2 = _basisTwo->getShellPairData();
    const auto& basis2 = _basisTwo->getBasis();

    /*
     * Thread safety: one buffer for each thread
     */
#ifdef _OPENMP
    std::vector<Eigen::MatrixXd> ints(omp_get_max_threads());
#else
    std::vector<Eigen::MatrixXd> ints(1);
#endif

#pragma omp parallel for schedule(static, 1)
    // loops over shells
    for (int pIndex = shellPairsB1->size() - 1; pIndex >= 0; --pIndex) {
      auto& p = (*shellPairsB1)[pIndex];

#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif

      const unsigned int i = p.bf1;
      const unsigned int j = p.bf2;
      const auto& basI = *basis1[i];
      const auto& basJ = *basis1[j];
      const unsigned int firstI = _basisOne->extendedIndex(i);
      const unsigned int firstJ = _basisOne->extendedIndex(j);

      for (auto& q : *shellPairsB2) {
        /*
         * Simple Prescreening, break out of loops early
         * sqrt[(ij|ij)*(rs|rs)] < threshold
         */
        if (p.factor * q.factor < _prescreeningThreshold)
          break;
        const unsigned int a = q.bf1;
        const unsigned int b = q.bf2;
        const auto& basA = *basis2[a];
        const auto& basB = *basis2[b];
        const unsigned int firstA = _basisTwo->extendedIndex(a);
        const unsigned int firstB = _basisTwo->extendedIndex(b);
        /*
         * Optional advanced prescreening
         */
        if (prescreenFunc(i, j, a, b, p.factor * q.factor))
          continue;
        if (firstJ > firstI || firstB > firstA)
          continue;

        // calculate integrals
        if (libint.compute(_op, _deriv, basI, basJ, basA, basB, ints[threadId])) {
          // unpack and run
          unsigned int counter = 0;
          for (unsigned int I = 0; I < basis1[i]->getNContracted(); ++I) {
            const unsigned int ii = _basisOne->extendedIndex(i) + I;
            for (unsigned int J = 0; J < basis1[j]->getNContracted(); ++J) {
              const unsigned int jj = _basisOne->extendedIndex(j) + J;
              for (unsigned int A = 0; A < basis2[a]->getNContracted(); ++A) {
                const unsigned int aa = _basisTwo->extendedIndex(a) + A;
                for (unsigned int B = 0; B < basis2[b]->getNContracted(); ++B, ++counter) {
                  const unsigned int bb = _basisTwo->extendedIndex(b) + B;
                  if (jj > ii || bb > aa)
                    continue;

                  Eigen::VectorXd set(ints[threadId].row(counter));
                  loopEvalFunction(ii, jj, aa, bb, set, threadId);
                } /* primitives of b -> B  */
              }   /* primitives of a -> A  */
            }     /* primitives of j -> J  */
          }       /* primitives of i -> I  */

        } /* if (prescreen) */

      } /* q/shellPairsTwo */
    }   /* p/shellPairsTwo */
    // finalize libint
    libint.freeEngines(_op, _deriv, 4);
  }

  /**
   * @brief Loops over 2e- 4-center interaction integrals.
   *
   * This function expects a lambda function as input.
   * The function will be called for each integral
   * (be aware of symmetry!).
   * The function has to expect four arguments:
   * @code
   *   void std::function<void (const unsigned int& i,
   *   const unsigned int& j,
   *   const unsigned int&  a,
   *   const unsigned int&  b,
   *   double intValues
   *   const unsigned int&  threadID)>
   * @endcode
   * Here i and j are the indices of the first basis while
   * a and b are the indces of the second basis.
   * \f[ \langle ij | ab \rangle \f]
   * All indices are unfolded (running over all contracted functions
   * in all shells.)
   * The routine uses the symmetry within i and j thus all function calls
   * should expect i<=j and account for the 'missing' calls.
   * The same holds for a and b.
   *
   * The function can than do anything with the given integrals,
   * e.g. sum up:
   * @code
   *  TwoElecFourCenterInteractionIntLooper looper(op, deriv, basis1, basis2);
   *  Eigen::VectorXd sum(nThreads);
   *  sum.setZero();
   *  auto const sumInts = [&sum]
   *                       (const unsigned int&  i,
   *                        const unsigned int&  j,
   *                        const unsigned int&  a,
   *                        const unsigned int&  b,
   *                        double intValues,
   *                        const unsigned int&  threadID) {
   *  sum[threadID] += intValues;
   *  if (i!=j) sum[threadID] += intValues;
   *  if (a!=b) sum[threadID] += intValues;
   *  if (a!=b && i!=j) sum[threadID] += intValues;
   *  };
   *  looper.loop(sumInts);
   *  // sum over thread contributions
   * @endcode
   *
   *
   * @param loopEvalFunction The function to use each integral, see above for extended description.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loopNoDerivative(Func distribute) {
    loopNoDerivative(distribute, [](unsigned int, unsigned int, unsigned int, unsigned int, double) { return false; });
  }

  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loopNoDerivative(Func distribute, PrescreenFunc prescreenFunc) {
    // intialize libint

    auto& libint = Libint::getInstance();
    libint.keepEngines(_op, 0, 4);
    libint.initialize(_op, 0, 4);
    auto& engines = libint.getFourCenterEngines(_op);

    const auto& shellPairsB1 = _basisOne->getShellPairData();
    const auto& basis1 = _basisOne->getBasis();
    const auto& shellPairsB2 = _basisTwo->getShellPairData();
    const auto& basis2 = _basisTwo->getBasis();

    /*
     * Thread safety: one buffer for each thread
     */
#ifdef _OPENMP
    std::vector<Eigen::MatrixXd> ints(omp_get_max_threads());
#else
    std::vector<Eigen::MatrixXd> ints(1);
#endif

#pragma omp parallel for schedule(static, 1)
    // loops over shells
    for (int pIndex = shellPairsB1->size() - 1; pIndex >= 0; --pIndex) {
      auto& p = (*shellPairsB1)[pIndex];

#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif

      const unsigned int i = p.bf1;
      const unsigned int j = p.bf2;
      const auto& basI = *basis1[i];
      const auto& basJ = *basis1[j];
      const unsigned int nI = basis1[i]->getNContracted();
      const unsigned int nJ = basis1[j]->getNContracted();
      const unsigned int firstI = _basisOne->extendedIndex(i);
      const unsigned int firstJ = _basisOne->extendedIndex(j);

      for (auto& q : *shellPairsB2) {
        /*
         * Simple Prescreening, break out of loops early
         * sqrt[(ij|ij)*(rs|rs)] < threshold
         */
        if (p.factor * q.factor < _prescreeningThreshold)
          break;
        const unsigned int a = q.bf1;
        const unsigned int b = q.bf2;
        const auto& basA = *basis2[a];
        const auto& basB = *basis2[b];
        const unsigned int firstA = _basisTwo->extendedIndex(a);
        const unsigned int firstB = _basisTwo->extendedIndex(b);
        const unsigned int nA = basis2[a]->getNContracted();
        const unsigned int nB = basis2[b]->getNContracted();
        /*
         * Optional advanced prescreening
         */
        if (prescreenFunc(i, j, a, b, p.factor * q.factor))
          continue;
        if (firstJ > firstI || firstB > firstA)
          continue;

        const auto& rawints = engines[threadId]->results();
        engines[threadId]->compute(basI, basJ, basA, basB);
        if (rawints[0] == nullptr) {
          continue;
        }

        const double* intptr = &(rawints[0][0]);
        // unpack and run
        for (unsigned int I = 0; I < nI; ++I) {
          const unsigned int ii = firstI + I;
          for (unsigned int J = 0; J < nJ; ++J) {
            const unsigned int jj = firstJ + J;
            for (unsigned int A = 0; A < nA; ++A) {
              const unsigned int aa = firstA + A;
              for (unsigned int B = 0; B < nB; ++B, ++intptr) {
                const unsigned int bb = firstB + B;
                if (jj > ii || bb > aa)
                  continue;

                distribute(ii, jj, aa, bb, (*intptr), threadId);
              } /* primitives of b -> B  */
            }   /* primitives of a -> A  */
          }     /* primitives of j -> J  */
        }       /* primitives of i -> I  */

      } /* q/shellPairsTwo */
    }   /* p/shellPairsTwo */
    // finalize libint
    libint.freeEngines(_op, 0, 4);
  }

 private:
  /// @brief The kernel/operator as libint enum.
  LIBINT_OPERATOR _op;
  /// @brief The derivative level.
  const unsigned int _deriv;
  /// @brief The basis.
  std::shared_ptr<BasisController> _basisOne;
  /// @brief The auxiliary basis.
  std::shared_ptr<BasisController> _basisTwo;
  /// @brief A threshold for the Schwartz conditions.
  double _prescreeningThreshold;
};

} /* namespace Serenity */

#endif /* COULOMBINTERACTIONINTLOOPER_H_ */
