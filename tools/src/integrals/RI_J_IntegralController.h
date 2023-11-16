/**
 * @file RI_J_IntegralController.h
 *
 * @date Mar 8, 2016
 * @author Kevin Klahr
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

#ifndef BASICS_INTEGRALS_RI_J_INTEGRALCONTROLLER_H_
#define BASICS_INTEGRALS_RI_J_INTEGRALCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "integrals/looper/ABTwoElecThreeCenterIntLooper.h"
#include "integrals/looper/TwoElecThreeCenterIntLooper.h"
#include "integrals/wrappers/Libint.h"
#include "math/Matrix.h"
#include "memory/MemoryManager.h"
#include "misc/SerenityError.h"
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <memory>
#include <mutex>
#include <vector>

namespace Serenity {

/* Forward declarations */
class Basis;
class BasisController;
class MemoryManager;
/**
 * @class RI_J_IntegralController RI_J_IntegralController.h
 * @brief A controller that handles the integrals for the RI-J approximation.
 * All of the 2 center ERIs as well as all or parts of the three center ERIs
 * are calculated once when needed and stored here for further use.
 *
 * Implementation according to:
 * [1] Weigend, F.; Kattannek, M.; Ahlrichs, R.; J.Chem.Phys (2009), 130, 164106 (eq. 4)
 *
 * See also:
 * [2] Neese, F.; J.Comput.Chem. (2003), 24, 1740 (Scheme 3)
 *
 */

class RI_J_IntegralController : public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @brief Constructor.
   * @param basisControllerA The first basis, for which integrals are to be managed.
   * @param auxBasisController The basis, the auxiliary basis.
   * @param basisControllerB A second basis, if the second index belongs to a different basis than first.
   * @param op The libint operator for these integrals.
   * @param mu Range-separation parameter if the operator is erf_coulomb.
   */
  RI_J_IntegralController(std::shared_ptr<BasisController> basisControllerA, std::shared_ptr<BasisController> auxBasisController,
                          std::shared_ptr<BasisController> basisControllerB = nullptr,
                          LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb, double mu = 0.0);
  /**
   * @brief Default destructor.
   */
  virtual ~RI_J_IntegralController() = default;

  /**
   * @brief Loops over 2e- 3-center integrals.
   *
   * This function expects a lambda function as input.
   * The function will be called for each integral
   * (be aware of symmetry!).
   * The function has to expect four arguments:
   * @code
   *   void std::function<void (const unsigned int& i,
   *   const unsigned int& j,
   *   const unsigned int& K,
   *   const double& intValue,
   *   const unsigned int threadId )>
   * @endcode
   * Here i and j are the indices of the 'normal' basis while
   * K is the index of the auxiliary basis.
   * All indices are unfolded (running over all contracted functions
   * in all shells.)
   * The routine uses the symmetry within i and j thus all function calls
   * should expect i<=j and account for the 'missing' calls.
   *
   * The function can than do anything with the given integrals,
   * e.g. sum up:
   * @code
   *  double sum = 0.0;
   *  auto const sumInts = [&sum]
   *                       (const unsigned int&  i,
   *                        const unsigned int&  j,
   *                        const unsigned int&  K,
   *                        const double& intValue,
   *                        const unsigned int threadId) {
   *  sum += intValues(0);
   *  if (i!=j) sum += intValues[0];
   *  };
   *  rijController.loopOver3CInts(sumInts);
   * @endcode
   *
   * The difference between this function and the looper is that this version
   * allows the RI_J_Controller to cache some of the values.
   * The loop then runs over cached values and afterwards all the other values
   * will be recalculated.
   *
   *
   * @param distribute The function to use each integral, see above for extended description.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loopOver3CInts(Func distribute) {
    auto prescreen = [](const unsigned, const unsigned, const unsigned, const double) { return false; };
    loopOver3CInts(distribute, prescreen);
  }

  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loopOver3CInts(Func distribute, PrescreenFunc prescreen) {
    // Determine the mode the method is run in.
    const bool twoBasisMode = _basisControllerB != nullptr;

    // Cache integrals if none are cached. We found that keeping three-center integrals
    // in memory is only really beneficial if the number of (SCF and not auxiliary)
    // basis functions is below approximately 1500.
    if (!_cache && !twoBasisMode && _basisControllerA->getNBasisFunctions() < 1500) {
      cache3CInts();
    }
    if (!_cache) {
      if (!twoBasisMode) {
        TwoElecThreeCenterIntLooper looper(_op, 0, _basisControllerA, _auxBasisController,
                                           _basisControllerA->getPrescreeningThreshold(),
                                           std::pair<unsigned int, unsigned int>(0, 0), _mu);
        looper.loopNoDerivative(distribute, prescreen);
      }
      else {
        double prescreeningThresholdA = _basisControllerA->getPrescreeningThreshold();
        double prescreeningThresholdB = _basisControllerB->getPrescreeningThreshold();
        double prescreeningThreshold = std::min(prescreeningThresholdA, prescreeningThresholdB);
        ABTwoElecThreeCenterIntLooper looper(_op, 0, _basisControllerA, _basisControllerB, _auxBasisController,
                                             prescreeningThreshold, std::pair<unsigned int, unsigned int>(0, 0), _mu);
        looper.loopNoDerivative(distribute);
      }
    }
    else {
      // if there was a cache, use it,
      auto data = _cache->data();
      const unsigned int colsize = _cache->rows();
      const unsigned int nBFs_A = _basisControllerA->getNBasisFunctions();
#pragma omp parallel for schedule(dynamic)
      for (unsigned int i = 0; i < nBFs_A; ++i) {
        const unsigned int jEnd = (twoBasisMode) ? _basisControllerB->getNBasisFunctions() - 1 : i;
        for (unsigned int j = 0; j <= jEnd; ++j) {
          const unsigned int threadId = omp_get_thread_num();
          const unsigned long long col =
              (!twoBasisMode) ? (unsigned long long)((i * (i + 1) / 2) + j) * colsize
                              : (unsigned long long)((i * (_basisControllerB->getNBasisFunctions()) + j) * colsize);
          for (unsigned int K = 0; K < colsize; ++K) {
            distribute(i, j, K, data[(unsigned long long)col + K], threadId);
          }
        }
      }
      // then run the rest.
      if (_cache->rows() != _auxBasisController->getNBasisFunctions()) {
        if (!twoBasisMode) {
          TwoElecThreeCenterIntLooper looper(
              _op, 0, _basisControllerA, _auxBasisController, _basisControllerA->getPrescreeningThreshold(),
              std::pair<unsigned int, unsigned int>(_cache->rows(), _auxBasisController->getNBasisFunctions()), _mu);
          looper.loopNoDerivative(distribute, prescreen);
        }
        else {
          throw SerenityError("Three center integral caching not possible for two regular basis controllers.");
        }
      }
    }
  };

  /**
   * @brief Empties the cache.
   */
  void clearCache() {
    _cache.reset(nullptr);
  }

 private:
  /**
   * @brief Checks if 3 center ints can be cached.
   *        Note: expects the cache to be a nullptr.
   */
  void cache3CInts();

  /// @brief The 2e3e integral cache.
  std::unique_ptr<Eigen::MatrixXd> _cache;

 public:
  /**
   * @brief Getter for the Cholesky decomposition of the Coulomb metric.
   * @return The Cholesky decomposition of the Coulomb metric.
   */
  const Eigen::LLT<Eigen::MatrixXd>& getLLTMetric();
  /**
   * @brief Getter for the Coulomb metric.
   * @return The Coulomb metric.
   */
  const Eigen::MatrixXd& getMetric();
  /**
   * @brief Getter for the (pseudo) inverse of the Coulomb metric.
   * @return The pseudo inverse of the Coulomb metric.
   */
  const Eigen::MatrixXd& getInverseM();
  /**
   * @brief Getter for the pseudo inverse square root of the Coulomb metric.
   * @return The pseudo inverse square root of the Coulomb metric.
   */
  const Eigen::MatrixXd& getInverseMSqrt();
  /**
   * @brief Getter for the basis controller.
   * @return The basis controller.
   */
  const std::shared_ptr<BasisController> getBasisController();
  /**
   * @brief Getter for the basis controller B.
   * @return The basis controller B.
   */
  const std::shared_ptr<BasisController> getBasisControllerB();
  /**
   * @brief Getter for the aux. basis controller.
   * @return The aux. basis controller.
   */
  const std::shared_ptr<BasisController> getAuxBasisController();
  /**
   * @brief Calculate the Coulomb metric.
   */
  void calculate2CenterIntegrals();
  /**
   * @brief Initialize the controller.
   */
  void initialize();
  /**
   * @brief Reset the controller.
   */
  void notify() override final {
    _M = nullptr;
    _lltM = nullptr;
    _inverseM.resize(0, 0);
    _inverseMSqrt.resize(0, 0);
    _cache = nullptr;
  }

 private:
  const std::shared_ptr<BasisController> _basisControllerA;
  std::shared_ptr<BasisController> _basisControllerB;
  const std::shared_ptr<BasisController> _auxBasisController;
  /// @brief The kernel/operator as libint enum.
  LIBINT_OPERATOR _op;
  /// @brief Range separation parameter mu
  double _mu;
  /*
   * The number of basis functions inside the 'normal' basis, i.e. the basis in which e.g.
   * the density and fock matrix are expressed.
   */
  const unsigned int _nBasisFunctions;
  /*
   * The number of basis functions forming the auxiliary basis
   */
  const unsigned int _nAuxFunctions;
  const unsigned int _nAuxFunctionsRed;

  std::vector<double> normFactors;
  /*
   * The Coulomb metric of the auxiliary basis
   */
  std::shared_ptr<Eigen::MatrixXd> _M;
  // The inverted metric.
  Eigen::MatrixXd _inverseM;
  // The inverted square root of the metric.
  Eigen::MatrixXd _inverseMSqrt;
  // The LLT decomposition of the metric.
  std::shared_ptr<Eigen::LLT<Eigen::MatrixXd>> _lltM;

  // The libint instance.
  const std::shared_ptr<Libint> _libint = Libint::getSharedPtr();
  // The memory manager.
  std::shared_ptr<MemoryManager> _memManager;
};

} /* namespace Serenity */
#endif /* BASICS_INTEGRALS_RI_J_INTEGRALCONTROLLER_H_ */
