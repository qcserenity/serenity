/**
 * @file   IntegralCachingController.h
 *
 * @date   Sep 24, 2021
 * @author Moritz Bensberg, Johannes Scheffler
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

#ifndef INTEGRALCACHINGCONTROLLER_H_
#define INTEGRALCACHINGCONTROLLER_H_

/* Include Serenity Internal Headers */
#include "basis/Shell.h"
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class BasisController;
class MemoryManager;
class Basis;

/**
 * @class IntegralCachingController IntegralCachingController.h
 * @brief This class handles the caching of 4-center integrals defined in a given basis set.
 * The caching works as follows:\n
 * If integrals from a shell quadruple qualify for caching (vide infra) and memory is
 * available for caching they are stored. The storing is organized for each outer
 * shell pair index in the 4-center loop individually. During the loop over all
 * 4-center integrals we keep track of the number of qualified integral sets for each
 * outer shell pair. Hence, if the we need an integral set, we check if it qualifies
 * for caching and then check if the counter for the number of qualified integral sets
 * is smaller than the size of the integral list stored for the given outer-shell pair.\n\n
 *
 * Integrals of type \f$ (a b|cd) \f$ are cached if \f$ t(A,B,C,D) \geq \tau \f$, where
 * \f$ \tau \f$ is a predefined selection threshold and the function \f$ t(A,B,C,D) \f$
 * is defined for the shells \f$ A,B,C \f$ and \f$ D \f$ as
 * \f$ t(A,B,C,D) = n_A n_B n_C n_D t_l \f$.\n
 * Here, \f$ n_A \f$ etc. is the contraction level of the basis function.
 * The factor \f$ t_l \f$ is given with the angular momenta of the shells
 * \f$ l_A, l_B, ... \f$ as \n
 * \f$ t_l = \begin{cases}
 * 6 , \text{if } l_A + l_B + l_C + l_D = 0 \\
 * 3 , \text{if } l_A + l_B + l_C + l_D = 1 \\
 * 2 , \text{if } l_A + l_B + l_C + l_D = 2 \\
 * 1 , \text{else}.
 * \end{cases}
 * \f$
 */
class IntegralCachingController : public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @brief Constructor.
   * @param basisController  The basis controller.
   * @param intCondition     The integral condition threshold.
   */
  IntegralCachingController(std::shared_ptr<BasisController> basisController, int intCondition);

  /**
   * @brief Destructor. Frees the memory used for the cache if necessary.
   */
  virtual ~IntegralCachingController();
  /**
   * @brief Add an integral set to the cache as calculated by libint.
   * @param ij        The outermost shell-pair index.
   * @param intptr    The pointer to the integral array.
   * @param size      The number of integrals in the set.
   * @param threadID  The thread id.
   */
  void cacheIntegral(unsigned ij, const double* intptr, size_t size, unsigned threadID);
  /**
   * @brief Getter for the integral array. Note that the integral caching depends on the loop
   *        structure used for the initial integral calculation.
   * @param ij          The outermost shell-pair index.
   * @param cdcounter   The innermost shell-pair counter.
   * @return The pointer to the integral set. Returns nullptr if the integral set was not cached.
   */
  const double* getIntegral(unsigned int ij, size_t cdcounter);
  /**
   * @brief Calculate the time condition \f$ t(aa,bb,cc,dd) \f$
   * @param aa The shell aa.
   * @param bb The shell bb.
   * @param cc The shell cc.
   * @param dd the shell dd.
   * @return Returns true if \f$ t(aa,bb,cc,dd) > \tau \f$
   */
  bool timeCondition(const Shell& aa, const Shell& bb, const Shell& cc, const Shell& dd);
  /**
   * @brief Reinitialize the amount of memory available for each thread.
   */
  void creatMemManagingVec();
  /**
   * @brief Check if memory is available for caching.
   * @return True if memory is available.
   */
  bool checkMem();
  /**
   * @brief Remove all integrals from memory.
   */
  void clearCache();
  /**
   * @brief If the basis controller associated to this object is changed,
   *        the cache has to be cleared..
   */
  void notify();
  /**
   * @brief If the prescreening is changed, the cache needs to be cleared. Note that
   * this threshold corresponds to the non-adjusted threshold during the SCF iterations,
   * i.e. integralThreshold in Settings.h!
   * @param threshold The potentially updated prescreening threshold.
   */
  void setPrescreeningThreshold(double threshold);

 private:
  bool _memAvailable = true;
  std::unique_ptr<std::vector<std::vector<std::vector<double>>>> _cacheVector;
  std::shared_ptr<BasisController> _basisController;
  int _intCondition;
  const unsigned _nThreads;
  std::shared_ptr<MemoryManager> _memManager;
  std::vector<double> _memoryPerThread;
  double _prescreeningThreshold = std::numeric_limits<double>::infinity();
};
} // namespace Serenity
#endif /* INTEGRALCHACHINGCONTROLLER_H_ */
