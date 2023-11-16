/**
 * @file EigenvalueSolver.h
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

#ifndef LRSCF_EIGENVALUESOLVER
#define LRSCF_EIGENVALUESOLVER

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Tools/IterativeSolver.h"
#include "settings/LRSCFOptions.h"

namespace Serenity {

namespace Options {
enum class LR_METHOD;
}

/**
 * @class EigenvalueSolver
 * @brief A class for the iterative solution of the response eigenvalue problem.\n
 *
 * This class includes solution techniques for Hermitian and symplectic eigenvalue\n
 * problems, as required by (Hybrid-) TDDFT.\n
 *
 * For Hermitian problems, the standard Davidson-Liu algorithm will be used and for\n
 * symplectic problems a modified OJJ-Algorithm as proposed by Filip Furche in his\n
 * PhD thesis.\n
 * (F. Furche, Göttingen 2002. Fak. f. Chemie, Diss. v. 1.2.2002., PhD thesis, 2002.)\n
 *
 * Note that the Hermitian solver utilizes non-orthonormal guess spaces for\n
 * convergence acceleration. (J. Chem. Theory Comput. 2016 12 (7), 3003-3007)
 */
class EigenvalueSolver : public IterativeSolver {
 public:
  /**
   * @brief Constructor.
   * @param nDimension Dimension of the TDDFT problem.
   * @param nEigen Number of the lowest eigenvalues to be obtained.
   * @param diagonal Orbital-energy difference of the TDDFT problem.
   * @param convergenceCriterion If the norm of a residual vector of a root falls beneath this
   *        threshold, the root will be considered converged.
   * @param maxIterations If not converged after this number of iterations, abort.
   * @param maxSubspaceDimension Will perform a subspace collapse if the number of guess vectors
   *        exceeds this threshold.
   * @param initialSubspace Multiplied with nEigen to obtain the size of the initial subspace
   *        (usually two or three is a common choice).
   * @param algorithm Variable to tell the eigenvalue solver which kind of TDDFT problem is to be solved
   *        to apply the according solution technique.
   * @param sigmaCalculator A lambda to conveniently form response matrix -- guess vector products.\n
   *        Takes a set of guessvectors as an argument and returns a pointer to the sigmavectors.
   * @param initialGuess The initial guess space might also be passed to the eigenvalue solver.
   * @param writeToDisk A lambda function to store temporary iteration data to disk.
   */
  EigenvalueSolver(
      unsigned nDimension, unsigned nEigen, const Eigen::VectorXd& diagonal, double convergenceCriterion,
      unsigned maxIterations, unsigned maxSubspaceDimension, unsigned initialSubspace, Options::LR_METHOD algorithm,
      std::function<std::unique_ptr<std::vector<Eigen::MatrixXd>>(std::vector<Eigen::MatrixXd>& guessVectors)> sigmaCalculator,
      std::shared_ptr<std::vector<Eigen::MatrixXd>> initialGuess = nullptr,
      std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk = [](std::vector<Eigen::MatrixXd>&,
                                                                                            Eigen::VectorXd&) {});

  /**
   * @brief Default destructor
   */
  virtual ~EigenvalueSolver() = default;

  /**
   * @brief Performs one iteration of the eigenvalue solver.
   */
  void iterate() override final;

  /**
   * @brief Obtains the left eigenvectors if the symmetric TDDFT problem was solved
   *        and/or normlizes the eigenvectors. Prints the converged eigenvalues in au.
   */
  void postProcessing() override;

 private:
  /**
   * @brief Initializes the eigenvalue solver by printing some info, getting the needed
   *        matrices right, calculating the seed and the corresponding sigma vectors.
   */
  void initialize() override final;

  ///@brief Number of initial unit-guess vectors.
  unsigned _initialSubspace;

  ///@brief Response type of the eigenvalue problem (symmetric, symmetrized or symplectic).
  Options::LR_METHOD _algorithm;
}; /* class EigenvalueSolver */
} /* namespace Serenity */

#endif /* LRSCF_EIGENVALUESOLVER */
