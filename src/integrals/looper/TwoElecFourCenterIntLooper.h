/**
 * @file   TwoElecFourCenterIntLooper.h
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   29. Oktober 2016, 10:53
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
#ifndef TWOELECFOURCENTERINTLOOPER_H
#define TWOELECFOURCENTERINTLOOPER_H
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "basis/Shell.h"
#include "integrals/IntegralCachingController.h"
#include "integrals/Normalization.h"
#include "integrals/wrappers/Libint.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <libint2/engine.h>
#pragma GCC diagnostic pop
#include <omp.h>
#include <Eigen/Dense>
#include <limits>
#include <memory>

namespace Serenity {
/* Forward Declarations */

/**
 * @class TwoElecFourCenterIntLooper TwoElecFourCenterIntLooper.h
 * @brief A looper for prescreened 2e- 4 center integrals.
 */
class TwoElecFourCenterIntLooper {
 public:
  /**
   * @brief Constructor
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param basisController The basis to be used.
   * @param prescreeningThreshold A threshold for the Schwarz conditions.
   * @param mu Parameter for range sepraration. To be used with erf_coulomb operator.
   */
  TwoElecFourCenterIntLooper(LIBINT_OPERATOR op, const unsigned deriv, std::shared_ptr<BasisController> basisController,
                             double prescreeningThreshold, double mu = 0.0)
    : _op(op), _deriv(deriv), _basisController(basisController), _prescreeningThreshold(prescreeningThreshold), _mu(mu) {
  }
  /**
   * @brief Default destructor.
   */
  virtual ~TwoElecFourCenterIntLooper() = default;

  /**
   * @brief Loops over prescreened 2e- 4-center integrals.
   *
   * This function expects a lambda function as input.
   * The function will be called for each integral
   * (be aware of symmetry!).
   * The function has to expect four arguments:
   * @code
   *   void std::function<void (const unsigned int& i,
   *                            const unsigned int&  j,
   *                            const unsigned int&  k,
   *                            const unsigned int&  l,
   *                            Eigen::VectorXd& intValues)>
   * @endcode
   * All indices are unfolded (running over all contracted functions
   * in all shells.)
   * The routine uses the symmetry within i,j,k,l to calculate the integral
   * values.
   * Thus for i!=j!=k!=l eight cases will have to be handled in one function.
   * call.
   *
   * The function can then do anything with the given integrals,
   * e.g. sum up:
   * @code
   *  TwoElecFourCenterIntLooper looper(op, deriv, basisi, 1E-10);
   *  Eigen::VectorXd sum(nThreads);
   *  sum.setZero();
   *  auto const sumInts = [&sum]
   *                       (const unsigned int&  i,
   *                        const unsigned int&  j,
   *                        const unsigned int&  k,
   *                        const unsigned int&  l,
   *                        Eigen::VectorXd& intValues,
   *                        const unsigned int&  threadID) {
   *  sum[threadID] += intValues(0);
   *  };
   *  looper.loop(sumInts);
   * @endcode
   *
   * For interaction terms between two basis sets see
   * the CoulombInteractionIntLooper and also the TwoElecThreeCenterIntLooper
   * for Rij approximated versions.
   *
   * MAKE SURE YOUR FUNCTION IS THREADSAFE!
   *
   * @param distributionFunction The function to use each integral, see above for extended description.
   *                             MAKE SURE YOUR FUNCTION IS THREADSAFE!
   * @param maxD                 Maximum coefficient to be contracted with the integrals. Used to adjust
   *                             the integral precision in case of absurdly large coefficients.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loop(Func distributionFunction, double maxD = 1, bool scaleIntegrals = false) {
    loop(
        distributionFunction, [](unsigned int, unsigned int, unsigned int, unsigned int, double) { return false; },
        maxD, scaleIntegrals);
  }
  /**
   * @brief Loops over all four-center two-electron Coulomb integrals and calls a function for each.
   * @param distributionFunction This function is called for each (significant) integral. CAUTION:
   *                             MAKE SURE YOUR FUNCTION IS THREADSAFE!
   * @param prescreenFunc a function for a more detailed prescreening, e.g. using the
   *                      density matrix. Also this function must, of course, be THREADSAFE!
   * @param maxD          Maximum coefficient to be contracted with the integrals. Used to adjust
   *                      the integral precision in case of absurdly large coefficients.
   */
  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loop(Func distributionFunction, PrescreenFunc prescreenFunc,
                                                  double maxD = 1, bool scaleIntegrals = false) {
    /*
     * Initialization
     */
    takeTime("init");
    const auto& shellPairs = _basisController->getShellPairData();
    const auto& basis = _basisController->getBasis();
    const unsigned int nBasisFunctions = _basisController->getNBasisFunctions();
    auto& libint = Libint::getInstance();
    libint.initialize(_op, _deriv, 4, std::vector<std::shared_ptr<Atom>>(0), _mu,
                      std::numeric_limits<double>::epsilon(), maxD, _basisController->getMaxNumberOfPrimitives());
    libint2::Operator libintOp = Libint::resolveLibintOperator(_op);
    /*
     * Thread safety: one buffer for each thread
     */
#ifdef _OPENMP
    std::vector<Eigen::MatrixXd> ints(omp_get_max_threads());
#else
    std::vector<Eigen::MatrixXd> ints(1);
#endif
    timeTaken(3, "init");
    /*
     * Begin loops
     */
    takeTime("calc");
/*
 * The static scheduling is important to make the results reproducible. It imposes a fixed
 * distribution of the loop to the different threads. A dynamic schedule would lead to a
 * different distribution to threads in each run. This would lead to numerical differences.
 */
#pragma omp parallel for schedule(static, 1)
    /* (Also copied)
     * We go this loop backwards, because then the execution time gets smaller the further
     * we advance in the loop. This can help in parallelizing efficiently.
     */
    for (int pIndex = shellPairs->size() - 1; pIndex >= 0; --pIndex) {
      auto& p = (*shellPairs)[pIndex];
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      const unsigned int i = p.bf1;
      const unsigned int j = p.bf2;
      const auto& basI = *basis[i];
      const auto& basJ = *basis[j];
      const unsigned int nI = basis[i]->getNContracted();
      const unsigned int nJ = basis[j]->getNContracted();
      const unsigned int firstI = _basisController->extendedIndex(i);
      const unsigned int firstJ = _basisController->extendedIndex(j);
      for (auto& q : *shellPairs) {
        /*
         * Simple Prescreening, break out of loops early
         */
        if (p.factor * q.factor < _prescreeningThreshold)
          break;
        const unsigned int k = q.bf1;
        const unsigned int l = q.bf2;
        const unsigned int firstK = _basisController->extendedIndex(k);
        const unsigned int firstL = _basisController->extendedIndex(l);
        const unsigned int nK = basis[k]->getNContracted();
        const unsigned int nL = basis[l]->getNContracted();
        /*
         * Optional advanced prescreening
         */
        if (prescreenFunc(i, j, k, l, p.factor * q.factor))
          continue;
        /*
         * Symmetry; only calculate integrals for which the combined index kl <= combined index ij.
         */
        if (firstK * nBasisFunctions + firstL > (firstI + nI) * nBasisFunctions + firstJ + nJ)
          continue;
        const auto& basK = *basis[k];
        const auto& basL = *basis[l];
        /*
         * Calculation of the integrals
         */
        bool significant = libint.compute(libintOp, _deriv, basI, basJ, basK, basL, ints[threadId]);
        if (!significant)
          continue;

        /*
         * Two-electron-four-center integrals have an eightfold permutational symmetry, i.e.
         * (ab|cd) = (ab|dc) = (ba|cd) = (ba|dc) = (cd|ab) = (cd|ba) = (dc|ab) = (dc|ba).
         * Due to the index restrictions above, we only loop over permutationally unique
         * quartets of shells. If, e.g., (a) == (b), that means that we have to scale the
         * integrals in this quartet by a factor of a half to account for the fact that we
         * will perform twice as many operations as actually necessary in any distribute
         * functions (which intrinsically assume that we have to do any contractions eight
         * times by default). Similar arguments hold for the other index combinations.
         *
         * Please note that the integrals are therefore coming out of this looper scaled,
         * so this is only a valid approach when the distribute-function assumes that any operation has to be performed
         * for all eight symmetric variants, e.g. when they are contracted on-the-fly in integral- direct approaches.
         * Use the boolean of this function to prevent the scaling if you need the integrals as they actually are.
         */
        /*
         * Unpack the integrals for the shell quadruple
         */
        double perm = 1.0;
        for (unsigned int ii = 0; ii < nI; ++ii) {
          const unsigned int iii = ii + firstI;
          for (unsigned int jj = 0; jj < nJ; ++jj) {
            const unsigned int jjj = jj + firstJ;
            if (jjj > iii)
              continue;
            for (unsigned int kk = 0; kk < nK; ++kk) {
              const unsigned int kkk = kk + firstK;
              for (unsigned int ll = 0; ll < nL; ++ll) {
                const unsigned int lll = ll + firstL;
                if (lll > kkk)
                  continue;
                if ((kkk)*nBasisFunctions + lll > (iii)*nBasisFunctions + jjj)
                  continue;
                perm = 1.0;
                if (scaleIntegrals) {
                  perm *= (iii == jjj) ? 0.5 : 1.0;
                  perm *= (kkk == lll) ? 0.5 : 1.0;
                  perm *= (iii == kkk) ? (jjj == lll ? 0.5 : 1.0) : 1.0;
                }
                const Eigen::VectorXd integral = perm * ints[threadId].row(ii * nJ * nK * nL + jj * nK * nL + kk * nL + ll);
                /*
                 * Call function to use the integral; once for each permutational symmetry,
                 * but only if a new index combination is actually achieved.
                 */
                distributionFunction(iii, jjj, kkk, lll, integral, threadId);
              }
            }
          }
        }
      }
    }
    libint.finalize(_op, _deriv, 4);
    timeTaken(3, "calc");
  }

  /**
   * @brief Loops over all four-center two-electron Coulomb integrals and calls a function for each.
   * @param distributionFunction This function is called for each (significant) integral. CAUTION:
   *                             MAKE SURE YOUR FUNCTION IS THREADSAFE!
   * @param maxD  Maximum coefficient to be contracted with the integrals. Used to adjust
   *             the integral precision in case of absurdly large coefficients.
   * @param scaleIntegrals If true, the integrals are scaled according to their permutational symmetry and reduncant
   *                       operations done during Fock matrix symmetrization as necessary for standard Fock matrix
   * builds. This function cannot handle derivatives of with respect to the nuclear coordinates. However, the integral
   * contraction is much more efficient due to double passing instead of Eigen::VectorXd.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loopNoDerivative(Func distribute, double maxD = 1, bool scaleIntegrals = false) {
    loopNoDerivative(
        distribute, [](unsigned int, unsigned int, unsigned int, unsigned int, double) { return false; }, maxD, nullptr,
        scaleIntegrals);
  }
  /**
   * @brief Loops over all four-center two-electron Coulomb integrals and calls a function for each.
   * @param distributionFunction This function is called for each (significant) integral. CAUTION:
   *                             MAKE SURE YOUR FUNCTION IS THREADSAFE!
   * @param prescreenFunc a function for a more detailed prescreening, e.g. using the
   *                      density matrix. Also this function must, of course, be THREADSAFE!
   * @param maxD  Maximum coefficient to be contracted with the integrals. Used to adjust
   *             the integral precision in case of absurdly large coefficients.
   * @param intcache The 4-center integral cache if available.
   * @param scaleIntegrals If true, the integrals are scaled according to their permutational symmetry and reduncant
   *                       operations done during Fock matrix symmetrization as necessary for standard Fock matrix
   * builds.
   *
   * This function cannot handle derivatives of with respect to the nuclear coordinates. However,
   * the integral contraction is much more efficient due to double passing instead of Eigen::VectorXd.
   */
  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loopNoDerivative(Func distribute, PrescreenFunc prescreen, double maxD = 1,
                                                              std::shared_ptr<IntegralCachingController> intcache = nullptr,
                                                              bool scaleIntegrals = false) {
    const auto& basis = _basisController->getBasis();
    const auto& shellPairs = _basisController->getShellPairData();
    unsigned ns = basis.size();
    auto schwarzParams = _basisController->getSchwarzParams(_op, _mu).data();

    auto& libint = Libint::getInstance();
    libint.initialize(_op, 0, 4, std::vector<std::shared_ptr<Atom>>(0), _mu, std::numeric_limits<double>::epsilon(),
                      maxD, _basisController->getMaxNumberOfPrimitives());
    auto& engines = libint.getFourCenterEngines(_op);

    if (intcache) {
      intcache->creatMemManagingVec();
      intcache->setPrescreeningThreshold(_prescreeningThreshold);
    }
#pragma omp parallel for schedule(dynamic)
    for (unsigned ab = 0; ab < shellPairs->size(); ++ab) {
      unsigned threadId = omp_get_thread_num();
      auto& abs = (*shellPairs)[ab];
      unsigned a = abs.bf1;
      unsigned b = abs.bf2;
      auto& aa = *basis[a];
      auto& bb = *basis[b];
      unsigned na = aa.getNContracted();
      unsigned nb = bb.getNContracted();
      unsigned aStart = _basisController->extendedIndex(a);
      unsigned bStart = _basisController->extendedIndex(b);
      unsigned long cdcounter = 0;
      for (unsigned c = 0; c <= a; ++c) {
        auto& cc = *basis[c];
        unsigned nc = cc.getNContracted();
        unsigned cStart = _basisController->extendedIndex(c);
        unsigned dmax = (a == c) ? b : c;
        for (unsigned d = 0; d <= dmax; ++d) {
          auto& dd = *basis[d];
          unsigned nd = dd.getNContracted();
          unsigned dStart = _basisController->extendedIndex(d);

          // crude schwarz and integral screening
          double schwarz = schwarzParams[a * ns + b] * schwarzParams[c * ns + d];
          if (schwarz < _prescreeningThreshold) {
            continue;
          }
          // selective integral caching
          const double* intptr;
          if (intcache) {
            if (intcache->timeCondition(aa, bb, cc, dd)) {
              intptr = intcache->getIntegral(ab, cdcounter);
              cdcounter += 1;
              if (!intptr) {
                if (intcache->checkMem()) {
                  const auto& rawints = engines[threadId]->results();
                  engines[threadId]->compute(aa, bb, cc, dd);
                  if (rawints[0] == nullptr) {
                    continue;
                  }
                  intptr = &(rawints[0][0]);
                  intcache->cacheIntegral(ab, intptr, na * nb * nc * nd, threadId);
                }
                else {
                  if (prescreen(a, b, c, d, schwarz)) {
                    continue;
                  }
                  const auto& rawints = engines[threadId]->results();
                  engines[threadId]->compute(aa, bb, cc, dd);
                  if (rawints[0] == nullptr) {
                    continue;
                  }
                  intptr = &(rawints[0][0]);
                }
              }
            }
            else {
              if (prescreen(a, b, c, d, schwarz)) {
                continue;
              }
              const auto& rawints = engines[threadId]->results();
              engines[threadId]->compute(aa, bb, cc, dd);
              if (rawints[0] == nullptr) {
                continue;
              }
              intptr = &(rawints[0][0]);
            }
          }
          else {
            if (prescreen(a, b, c, d, schwarz)) {
              continue;
            }
            const auto& rawints = engines[threadId]->results();
            engines[threadId]->compute(aa, bb, cc, dd);
            if (rawints[0] == nullptr) {
              continue;
            }
            intptr = &(rawints[0][0]);
          }

          double perm = 1.0;
          /*
           * Two-electron-four-center integrals have an eightfold permutational symmetry, i.e.
           * (ab|cd) = (ab|dc) = (ba|cd) = (ba|dc) = (cd|ab) = (cd|ba) = (dc|ab) = (dc|ba).
           * Due to the index restrictions above, we only loop over permutationally unique
           * quartets of shells. If, e.g., (a) == (b), that means that we have to scale the
           * integrals in this quartet by a factor of a half to account for the fact that we
           * will perform twice as many operations as actually necessary in any distribute
           * functions (which intrinsically assume that we have to do any contractions eight
           * times by default). Similar arguments hold for the other index combinations.
           *
           * Please note that the integrals are therefore coming out of this looper scaled,
           * so this is only a valid approach when they are contracted on-the-fly in integral-
           * direct approaches. Use the boolean of this function to prevent the scaling
           * if you need the integrals as they actually are.
           */
          if (scaleIntegrals) {
            perm *= (a == b) ? 0.5 : 1.0;
            perm *= (c == d) ? 0.5 : 1.0;
            perm *= (a == c) ? (b == d ? 0.5 : 1.0) : 1.0;
          }

          // Handle normalization in case of cartesian functions
          // Then distribute
          // Assumes that either all shells are cartesian or all shells are spherical
          if (aa.contr[0].pure) {
            for (unsigned aaa = aStart; aaa < aStart + na; ++aaa) {
              for (unsigned bbb = bStart; bbb < bStart + nb; ++bbb) {
                for (unsigned ccc = cStart; ccc < cStart + nc; ++ccc) {
                  for (unsigned ddd = dStart; ddd < dStart + nd; ++ddd, ++intptr) {
                    distribute(aaa, bbb, ccc, ddd, perm * (*intptr), threadId);
                  } /* bfs d */
                }   /* bfs c */
              }     /* bfs b */
            }       /* bfs a */
          }
          else {
            const auto nInts = na * nb * nc * nd;
            Eigen::VectorXd local(Eigen::Map<const Eigen::VectorXd>(intptr, nInts));
            Normalization::normalizeShell(local, aa.contr[0].l, bb.contr[0].l, cc.contr[0].l, dd.contr[0].l);
            const double* intptr = local.data();
            for (unsigned aaa = aStart; aaa < aStart + na; ++aaa) {
              for (unsigned bbb = bStart; bbb < bStart + nb; ++bbb) {
                for (unsigned ccc = cStart; ccc < cStart + nc; ++ccc) {
                  for (unsigned ddd = dStart; ddd < dStart + nd; ++ddd, ++intptr) {
                    distribute(aaa, bbb, ccc, ddd, perm * (*intptr), threadId);
                  } /* bfs d */
                }   /* bfs c */
              }     /* bfs b */
            }       /* bfs a */
          }         /* cartesisan vs. spherical */
        }           /* shell d */
      }             /* shell c */
    }               /* shell a/b */

    libint.finalize(_op, 0, 4);
  }

  /**
   * @brief Loops over all diagonal elements of the four-center two-electron Coulomb integrals and
   *        calls a function for each.
   * @param distributionFunction This function is called for each (significant) integral. CAUTION:
   *                             MAKE SURE YOUR FUNCTION IS THREADSAFE! For more details look at the loop()
   *                             function.
   */
  template<class Func>
  void loopDiagonal(Func distributionFunction) {
    loopDiagonal(distributionFunction, [](unsigned int, unsigned int, unsigned int, unsigned int, unsigned int,
                                          unsigned int, unsigned int, unsigned int, double) { return false; });
  }
  /**
   * @brief Loops over all diagonal elements of the four-center two-electron Coulomb integrals and
   *        calls a function for each.
   * @param distributionFunction This function is called for each (significant) integral. CAUTION:
   *                             MAKE SURE YOUR FUNCTION IS THREADSAFE!
   * @param prescreenFunc a function for a more detailed prescreening, e.g. using the
   *                      density matrix. Also this function must, of course, be THREADSAFE!
   */
  template<class Func, class PrescreenFunc>
  void loopDiagonal(Func distributionFunction, PrescreenFunc prescreenFunc) {
    /*
     * Initialization
     */
    takeTime("init");
    const auto& shellPairs = _basisController->getShellPairData();
    const auto& basis = _basisController->getBasis();
    auto& libint = Libint::getInstance();
    libint.finalize(_op, _deriv, 4);
    libint.initialize(_op, _deriv, 4, std::vector<std::shared_ptr<Atom>>(0), _mu);
    /*
     * Thread safety: one buffer for each thread
     */
#ifdef _OPENMP
    std::vector<Eigen::MatrixXd> ints(omp_get_max_threads());
#else
    std::vector<Eigen::MatrixXd> ints(1);
#endif
    timeTaken(3, "init");
    takeTime("calc");
    /*
     * TODO the parallelization here at the outermost loop is very inefficient. However, parallelizing
     *      the next innermost loop is also not good, because breaking out of it is then not possible
     *      any more.
     * (Copied from an older routine, I guess it should actually be the same here.)
     * The static scheduling is important to make the results reproducible. It imposes a fixed
     * distribution of the loop to the different threads. A dynamic schedule would lead to a
     * different distribution to threads in each run. This would lead to numerical differences.
     */
#pragma omp parallel for schedule(dynamic)
    /* (Also copied)
     * We go this loop backwards, because then the execution time gets smaller the further
     * we advance in the loop. This can help in parallelizing efficiently.
     */
    for (int pIndex = shellPairs->size() - 1; pIndex >= 0; --pIndex) {
      auto& p = (*shellPairs)[pIndex];
      if (p.factor * p.factor < _prescreeningThreshold)
        continue;
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      const unsigned int i = p.bf1;
      const unsigned int j = p.bf2;
      const auto& basI = *basis[i];
      const auto& basJ = *basis[j];
      const unsigned int nI = basis[i]->getNContracted();
      const unsigned int nJ = basis[j]->getNContracted();
      const unsigned int firstI = _basisController->extendedIndex(i);
      const unsigned int firstJ = _basisController->extendedIndex(j);
      if (prescreenFunc(firstI, firstJ, firstI, firstJ, nI, nJ, nI, nJ, p.factor * p.factor))
        continue;
      bool significant = libint.compute(_op, _deriv, basI, basJ, basI, basJ, ints[threadId]);
      if (!significant)
        continue;
      for (unsigned int ii = 0; ii < nI; ++ii) {
        const unsigned int iii = ii + firstI;
        for (unsigned int jj = 0; jj < nJ; ++jj) {
          const unsigned int jjj = jj + firstJ;
          if (jjj > iii)
            continue;
          const Eigen::VectorXd integral = ints[threadId].row(ii * nJ * nI * nJ + jj * nI * nJ + ii * nJ + jj);
          distributionFunction(iii, jjj, integral, threadId);
        }
      }
    }
    libint.finalize(_op, _deriv, 4);
    timeTaken(3, "calc");
  }

 private:
  /// @brief The kernel/operator as libint enum.
  LIBINT_OPERATOR _op;
  /// @brief The derivative level.
  const unsigned int _deriv;
  /// @brief The basis.
  std::shared_ptr<BasisController> _basisController;
  /// @brief A threshold for the Schwarz conditions.
  double _prescreeningThreshold;
  /// @brief Range-separation parameter
  double _mu;
};

} /* namespace Serenity */
#endif /* TWOELECFOURCENTERINTLOOPER_H */
