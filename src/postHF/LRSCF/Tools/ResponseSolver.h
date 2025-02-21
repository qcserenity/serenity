/**
 * @file ResponseSolver.h
 *
 * @date May 30, 2019
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

#ifndef LRSCF_RESPONSESOLVER
#define LRSCF_RESPONSESOLVER

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Tools/IterativeSolver.h"
#include "postHF/LRSCF/Tools/SigmaCalculator.h"

namespace Serenity {

namespace Options {
enum class GAUGE;
}
/**
 * @class ResponseSolver
 * @brief A class for the iterative solution of the response problem (in the form of a linear system).\n
 *
 * For undamped response problems (frequency is purely real) it utilizes a modified OJJ-Algorithm\n
 * as proposed by Filip Furche in his PhD thesis.\n
 * (F. Furche, Göttingen 2002. Fak. f. Chemie, Diss. v. 1.2.2002., PhD thesis, 2002.)\n
 *
 * For damped response problems (frequency is complex) it utilizes a modified OJJ-Algorithm adapted\n
 * for this purpose (see Niklas' Master Thesis for further reference).
 *
 * AR: Caution! This class respects only the first three columns of the right-hand side matrix ppmq. See the seed()
 * function
 */
class ResponseSolver : public IterativeSolver {
 public:
  /**
   * @brief Constructor.
   * @param diagonal Orbital-energy difference of the TDDFT problem.
   * @param convergenceCriterion If the norm of a residual vector of a root falls beneath this
   *        threshold, the root will be considered converged.
   * @param maxIterations If not converged after this number of iterations, abort.
   * @param maxSubspaceDimension Will perform a subspace collapse if the number of guess vectors
   *        exceeds this threshold.
   * @param frequencies The frequencies to be solved for.
   * @param damping The implicitly imaginary damping factor (in eV).
   * @param ppmq The right-hand side of the linear system.
   * @param sigmaCalculator A lambda to conveniently form response matrix -- guess vector products.\n
   *        Takes a set of guessvectors as an argument and returns a pointer to the sigmavectors.
   * @param initialGuess The initial guess space might also be passed to the response solver.
   * @param writeToDisk A lambda function to store temporary iteration data to disk.
   */
  ResponseSolver(
      const Eigen::VectorXd& diagonal, double convergenceCriterion, unsigned maxIterations,
      unsigned maxSubspaceDimension, std::vector<double>& frequencies, double damping, std::vector<Eigen::MatrixXd> ppmq,
      SigmaCalculator sigmaCalculator, std::shared_ptr<std::vector<Eigen::MatrixXd>> initialGuess = nullptr,
      std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk = [](std::vector<Eigen::MatrixXd>&,
                                                                                            Eigen::VectorXd&) {});

  /**
   * @brief Default destructor
   */
  virtual ~ResponseSolver() = default;

  /**
   * @brief Performs one iteration of the response solver.
   */
  void iterate() override final;

  /**
   * @brief No postprocessing needed for response solver (inherited from IterativeSolver).
   */
  void postProcessing() override final{};

 private:
  /**
   * @brief Initializes the response solver by printing some info, getting the needed
   *        matrices right, calculating the seed and the corresponding sigma vectors.
   */
  void initialize() override final;

  /**
   * @brief Generates the initial guess spaces for each frequency and right-hand side.
   *        Uses a block-diagonal approximation for the coefficient matrix of the linear problem.
   */
  void seed();

  /**
   * @brief Preconditiones the residual vectors with the inverse of the approximate
   *        coefficient matrix of the linear problem.
   *        The same block-diagonal approximation as for the seed() function is used.
   */
  void precondition();

  ///@brief The frequencies to be solved for.
  std::vector<double> _frequencies;
  ///@brief A vector determining if a frequency is already converged (so it can be skipped).
  Eigen::VectorXi _freqConverged;
  ///@brief The damping factor.
  double _damping;
  ///@brief Simple bool for if statements.
  bool _damped;
  ///@brief The number of frequencies.
  unsigned int _nFreqs;
  ///@brief The right-hand sides.
  std::vector<Eigen::MatrixXd> _ppmq;

}; /* class ResponseSolver */
} /* namespace Serenity */

#endif /* LRSCF_RESPONSESOLVER */
