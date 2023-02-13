/**
 * @file TwoElecThreeCenterCalculator.h
 *
 * @date Jul 12, 2021
 * @author Niklas Niemeyer
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
#ifndef TWOELECTHREECENTERCALCULATOR_H_
#define TWOELECTHREECENTERCALCULATOR_H_
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "integrals/wrappers/Libint.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/**
 * @class TwoElecThreeCenterCalculator TwoElecThreeCenterCalculator.h
 * @brief A looper for 2e- 3 center integrals between two basis sets.
 */
class TwoElecThreeCenterCalculator {
 public:
  /**
   * @brief Constructor.
   * @param op The kernel/operator as libint enum.
   * @param mu Range separation parameter if operator is erf_coulomb.
   * @param basisA The basis A.
   * @param basisB The basis B.
   * @param auxbasis The auxiliary basis.
   * @param prescreeningThreshold Internal prescreening threshold.
   * @param maxDens Maximum density matrix entry for libint initializations.
   */
  TwoElecThreeCenterCalculator(LIBINT_OPERATOR op, double mu, std::shared_ptr<BasisController> basisA,
                               std::shared_ptr<BasisController> basisB, std::shared_ptr<BasisController> auxbasis,
                               double prescreeningThreshold, double maxDens = 10);

  /**
   * @brief Destructor.
   * Finalizes libint engines used for this calculator.
   */
  virtual ~TwoElecThreeCenterCalculator();

  /**
   * @brief Calculates all (mu nu|p) integrals for p in the aux basis shell with the index P.
   *        Note: Please be careful since this function requires libint to be initialized correctly beforehand.
   * @param P Shell index of the aux basis.
   * @param iThread The thread ID.
   * @return A matrix with nu*mu rows and nBF(shell P) columns for easy mapping.
   */
  Eigen::Ref<Eigen::MatrixXd> calculateIntegrals(unsigned P, unsigned iThread);

  /**
   * @brief Sets up shell pairs for the regular bases, i.e. the mu and nu indices in (mu nu|p). Note that
   * mu and nu might be from different basis sets and this class will then setup those mixed shell pairs upon
   * construction
   */
  void setupShellPairs();

  /**
   * @brief Loops over all sets of integrals (mu nu|P).
   *
   * @param distribute The distribute function.
   */
  __attribute__((always_inline)) inline void
  loop(std::function<void(Eigen::Map<Eigen::MatrixXd> AO, size_t P, unsigned iThread)> distribute) {
    auto& auxbasis = _auxbasis->getBasis();

#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < _auxbasis->getReducedNBasisFunctions(); ++iShell) {
      unsigned long P_in_iShell = auxbasis[iShell]->getNContracted();
      unsigned iThread = omp_get_thread_num();
      auto integrals = (iShell < _nss) ? _cache.middleCols(_offsets[iShell], P_in_iShell)
                                       : this->calculateIntegrals(iShell, iThread);

      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned long P_all = _auxbasis->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> AO(integrals.col(P).data(), _nb_B, _nb_A);
        distribute(AO, P_all, iThread);
      }
    }
  }

  /**
   * @brief Caches integrals of the type (mu nu|P).
   */
  void cacheIntegrals();

  /**
   * @brief Clear integral cache.
   */
  void clearCache();

 private:
  // Libint engine.
  const std::shared_ptr<Libint> _libint = Libint::getSharedPtr();
  /// @brief The kernel/operator as libint enum.
  LIBINT_OPERATOR _op;
  /// @brief The first basis.
  std::shared_ptr<BasisController> _basisControllerA;
  /// @brief The second basis.
  std::shared_ptr<BasisController> _basisControllerB;
  /// @brief Bool indicating whether( _basisControllerA == _basisControllerB).
  bool _twoBasisMode;
  /// @brief Shell-pair data (one basis set or two)
  std::shared_ptr<std::vector<ShellPairData>> _shellPairData;
  /// @brief The auxiliary basis.
  std::shared_ptr<BasisController> _auxbasis;
  /// @brief A threshold for the Schwarz conditions.
  double _prescreeningThreshold;
  /// @brief Range-separation parameter mu.
  double _mu;
  /// @brief The integrals.
  std::vector<std::unique_ptr<Eigen::MatrixXd>> _integrals;
  /// @brief Number of basis functions in the first basis.
  unsigned _nb_A;
  /// @brief Number of basis functions in the second basis.
  unsigned _nb_B;
  /// @brief Cached integrals.
  Eigen::MatrixXd _cache;
  /// @brief Stored shells.
  unsigned _nss;
  /// @brief Basis function offsets for all shell indices.
  std::vector<size_t> _offsets;
};
} /* namespace Serenity */
#endif /* TWOELECTHREECENTERCALCULATOR_H_ */
