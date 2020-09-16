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
#include "basis/BasisController.h"
#include "basis/Shell.h"
#include "integrals/wrappers/Libint.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <omp.h>
#include <Eigen/Dense>
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
  TwoElecFourCenterIntLooper(libint2::Operator op, const unsigned deriv, std::shared_ptr<BasisController> basisController,
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
   * values but it will call the function for each of the symmetric ones.
   * Thus for i!=j!=k!=l eight calls will be done.
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
   */
  template<class Func>
  __attribute__((always_inline)) inline void loop(Func distributionFunction) {
    loop(distributionFunction, [](unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int,
                                  unsigned int, unsigned int, double) { return false; });
  }
  /**
   * @brief Loops over all four-center two-electron Coulomb integrals and calls a function for each.
   * @param distributionFunction This function is called for each (significant) integral. CAUTION:
   *                             MAKE SURE YOUR FUNCTION IS THREADSAFE!
   * @param prescreenFunc a function for a more detailed prescreening, e.g. using the
   *                      density matrix. Also this function must, of course, be THREADSAFE!
   */
  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loop(Func distributionFunction, PrescreenFunc prescreenFunc) {
    /*
     * Initialization
     */
    takeTime("init");
    const auto& shellPairs = _basisController->getShellPairData();
    const auto& basis = _basisController->getBasis();
    const unsigned int nBasisFunctions = _basisController->getNBasisFunctions();
    auto& libint = Libint::getInstance();
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
        if (prescreenFunc(firstI, firstJ, firstK, firstL, nI, nJ, nK, nL, p.factor * q.factor))
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
        bool significant = libint.compute(_op, _deriv, basI, basJ, basK, basL, ints[threadId]);
        if (!significant)
          continue;
        /*
         * Unpack the integrals for the shell quadruple
         */
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
                const Eigen::VectorXd integral(ints[threadId].row(ii * nJ * nK * nL + jj * nK * nL + kk * nL + ll));
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
   * This function cannot handle derivatives of with respect to the nuclear coordinates. However,
   * the integral contraction is much more efficient due to double passing instead of Eigen::VectorXd.
   *
   * It cannot handle normalization, either.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loopNoDerivative(Func distribute) {
    loopNoDerivative(distribute, [](unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int,
                                    unsigned int, unsigned int, double) { return false; });
  }
  /**
   * @brief Loops over all four-center two-electron Coulomb integrals and calls a function for each.
   * @param distributionFunction This function is called for each (significant) integral. CAUTION:
   *                             MAKE SURE YOUR FUNCTION IS THREADSAFE!
   * @param prescreenFunc a function for a more detailed prescreening, e.g. using the
   *                      density matrix. Also this function must, of course, be THREADSAFE!
   *
   * This function cannot handle derivatives of with respect to the nuclear coordinates. However,
   * the integral contraction is much more efficient due to double passing instead of Eigen::VectorXd.
   *
   * It cannot handle normalization, either.
   */
  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loopNoDerivative(Func distribute, PrescreenFunc prescreen) {
    const auto& basis = _basisController->getBasis();
    const auto& shellPairs = _basisController->getShellPairData();
    unsigned ns = basis.size();
    unsigned nbf = _basisController->getNBasisFunctions();
    auto schwartzParams = _basisController->getSchwartzParams().data();

    auto& libint = Libint::getInstance();
    libint.initialize(_op, 0, 4, std::vector<std::shared_ptr<Atom>>(0), _mu);
    auto& engines = libint.getFourCenterEngines(_op);

#pragma omp parallel for schedule(dynamic)
    for (unsigned ab = 0; ab < shellPairs->size(); ++ab) {
      unsigned threadId = omp_get_thread_num();
      auto& abs = (*shellPairs)[ab];
      unsigned a = abs.bf1;
      unsigned b = abs.bf2;
      auto& aa = *basis[a];
      auto& bb = *basis[b];
      unsigned na = aa.size();
      unsigned nb = bb.size();
      unsigned aStart = _basisController->extendedIndex(a);
      unsigned bStart = _basisController->extendedIndex(b);
      for (unsigned c = 0; c <= a; ++c) {
        auto& cc = *basis[c];
        unsigned nc = cc.size();
        unsigned cStart = _basisController->extendedIndex(c);
        for (unsigned d = 0; d <= c; ++d) {
          auto& dd = *basis[d];
          unsigned nd = dd.size();
          unsigned dStart = _basisController->extendedIndex(d);
          unsigned ncnd = nc * nd;

          // crude schwartz and integral screening
          double schwartz = schwartzParams[a * ns + b] * schwartzParams[c * ns + d];
          if (schwartz < _prescreeningThreshold || prescreen(aStart, bStart, cStart, dStart, na, nb, nc, nd, schwartz)) {
            continue;
          }

          const auto& rawints = engines[threadId]->results();
          engines[threadId]->compute(aa, bb, cc, dd);
          if (rawints[0] == nullptr) {
            continue;
          }

          const double* intptr = &(rawints[0][0]);
          for (unsigned aaa = aStart; aaa < aStart + na; ++aaa) {
            for (unsigned bbb = bStart; bbb < bStart + nb; ++bbb) {
              if (bbb > aaa) {
                intptr += ncnd;
                continue;
              }
              for (unsigned ccc = cStart; ccc < cStart + nc; ++ccc) {
                for (unsigned ddd = dStart; ddd < dStart + nd; ++ddd, ++intptr) {
                  if (ccc * nbf + ddd > aaa * nbf + bbb || ddd > ccc) {
                    continue;
                  }
                  distribute(aaa, bbb, ccc, ddd, *intptr, threadId);
                } /* bfs d */
              }   /* bfs c */
            }     /* bfs b */
          }       /* bfs a */
        }         /* shell d */
      }           /* shell c */
    }             /* shell a/b */

    libint.finalize(_op, 0, 4);
  }

 private:
  /// @brief The kernel/operator as libint enum.
  libint2::Operator _op;
  /// @brief The derivative level.
  const unsigned int _deriv;
  /// @brief The basis.
  std::shared_ptr<BasisController> _basisController;
  /// @brief A threshold for the Schwartz conditions.
  double _prescreeningThreshold;
  /// @brief Range separation parameter mu
  double _mu;
};

} /* namespace Serenity */
#endif /* TWOELECFOURCENTERINTLOOPER_H */
