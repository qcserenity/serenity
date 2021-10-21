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
#include "integrals/wrappers/Libint.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

class BasisController;
class ShellPairData;

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
  Eigen::MatrixXd& calculateIntegrals(unsigned P, unsigned iThread);

  /**
   * @brief Sets up shell pairs for the regular bases, i.e. the mu and nu indices in (mu nu|p). Note that
   * mu and nu might be from different basis sets and this class will then setup those mixed shell pairs upon
   * construction
   */
  void setupShellPairs();

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
};
} /* namespace Serenity */
#endif /* TWOELECTHREECENTERCALCULATOR_H_ */
