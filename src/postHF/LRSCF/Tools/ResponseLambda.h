/**
 * @file ResponseLambda.h
 *
 * @date Nov 13, 2020
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

#ifndef LRSCF_RESPONSELAMBDA
#define LRSCF_RESPONSELAMBDA

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Tools/SigmaCalculator.h"
#include "settings/GridOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

struct LRSCFTaskSettings;

template<Options::SCF_MODES SCFMode>
class SimplifiedTDDFT;

template<Options::SCF_MODES SCFMode>
class Sigmavector;

template<Options::SCF_MODES SCFMode>
class Kernel;

/**
 * @class ResponseLambda.
 * @brief Sets up all needed sigma vector calculators.
 */
template<Options::SCF_MODES SCFMode>
class ResponseLambda {
 public:
  /**
   * @brief Constructor.
   * @param act The active systems.
   * @param env The environment systems.
   * @param lrscf The LRSCFController.
   * @param diagonal The orbital-energy differences.
   * @param settings The LRSCFTask settings.
   */
  ResponseLambda(std::vector<std::shared_ptr<SystemController>> act, std::vector<std::shared_ptr<SystemController>> env,
                 std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, const Eigen::VectorXd& diagonal,
                 LRSCFTaskSettings& settings);

  /**
   * @brief Destructor.
   */
  virtual ~ResponseLambda();

  /**
   * @brief Sets up the sigma vector calculator lambda functions for TDA/TDDFT/BSE.
   */
  void setupTDDFTLambdas();

  /**
   * @brief Returns the TDA sigma vector calculator.
   * @return The TDA sigma vector calculator.
   */
  SigmaCalculator getTDASigma();

  /**
   * @brief Returns the TDDFT sigma vector calculator.
   * @return The TDDFT sigma vector calculator.
   */
  SigmaCalculator getRPASigma();

  /**
   * @brief Sets up the sigma vector calculator lambda functions for CC2/CIS(Dinf)/ADC(2).
   */
  void setupCC2Lambdas();

  /**
   * @brief Returns the right Jacobian sigma vector calculator.
   * @return The right Jacobian sigma vector calculator.
   */
  NonlinearSigmaCalculator getRightCC2Sigma();

  /**
   * @brief Returns the left Jacobian sigma vector calculator.
   * @return The left Jacobian sigma vector calculator.
   */
  NonlinearSigmaCalculator getLeftCC2Sigma();

  /**
   * @brief Returns whether a double hybrid is used.
   * @return True if a double hybrid is used.
   */
  bool usesDoubleHybrid();

  /**
   * @brief Returns MP2 correlation ratio of the active systems DH functional. Only supports isolated and FDEu
   * calculations.
   * @return The MP2 correlation ratio of the active systems DH functional.
   */
  double getDHRatio();

  /**
   * @brief Returns whether the kernel is used or not.
   * @return True if the kernel is used.
   */
  bool usesKernel();

  /**
   * @brief Returns whether exact exchange used or not.
   * @return True if exact exchange is used.
   */
  bool usesExchange();

  /**
   * @brief Returns whether exact LR exchange used or not.
   * @return True if LR exchange is used.
   */
  bool usesLRExchange();

  /**
   * @brief Returns whether external orthogonality is used or not.
   * @return True if external orthogonality is used.
   */
  bool usesEO();

  /**
   * @brief Setup kernel with given grid purpose.
   * @param gridAccuracy The grid accuracy.
   */
  void setupKernel(Options::GRID_PURPOSES gridAccuracy = Options::GRID_PURPOSES::DEFAULT);

 private:
  ///@brief The active systems.
  std::vector<std::shared_ptr<SystemController>> _act;

  ///@brief The environment systems.
  std::vector<std::shared_ptr<SystemController>> _env;

  ///@brief The LRSCFController.
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;

  ///@brief The orbital-energy differences.
  const Eigen::VectorXd& _diagonal;

  ///@brief The LRSCFTask settings.
  LRSCFTaskSettings& _settings;

  ///@brief Using double-hybrid functional?
  bool _usesDoubleHybrid = false;

  ///@brief Using double-hybrid functional?
  double _doubleHybridCorrRatio = 0.0;

  ///@brief Using kernel?
  bool _usesKernel = false;

  ///@brief Using exact exchange?
  bool _usesExchange = false;

  ///@brief Using external orthogonality?
  bool _usesEO = false;

  ///@brief Using long-range exact exchange?
  bool _usesLRExchange = false;

  ///@brief Using RI-J, -K, -erfK?
  bool _densFitJ, _densFitK, _densFitLRK;

  ///@brief SCF factor.
  double _scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;

  ///@brief Fock, Coulomb, exchange, RI exchange, kernel and EO sigma vector.
  std::unique_ptr<Sigmavector<SCFMode>> D, J, DFJ, K, DFK, F, EO, G;

  ///@brief The kernel.
  std::shared_ptr<Kernel<SCFMode>> _kernel = nullptr;
  std::shared_ptr<Kernel<UNRESTRICTED>> _ukernel = nullptr;

  ///@brief The TDA sigma vector calculator.
  SigmaCalculator _TDA;

  ///@brief The TDDFT sigma vector calculator.
  SigmaCalculator _RPA;

  ///@brief The right Jacobian sigma vector calculator.
  NonlinearSigmaCalculator _rightCC2;

  ///@brief The left Jacobian sigma vector calculator.
  NonlinearSigmaCalculator _leftCC2;

  ///@brief The simplified TDDFT object.
  std::shared_ptr<SimplifiedTDDFT<SCFMode>> _simplifiedTDDFT;
};

} /* namespace Serenity */

#endif /* LRSCF_RESPONSELAMBDA */
