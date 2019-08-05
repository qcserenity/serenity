/**
 * @file DavidsonSolver.h
 *
 * @date Aug 31, 2016
 * @author Michael Boeckers
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef BASICS_MATH_LINEARALGEBRA_DAVIDSONSOLVER_H_
#define BASICS_MATH_LINEARALGEBRA_DAVIDSONSOLVER_H_

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Solver/IterativeEigenvalueSolver.h"
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <functional>


namespace Serenity {

class DavidsonSolver: public IterativeEigenvalueSolver {
public:
  /**
   *
   * @param nDimension
   * @param nEigen
   * @param diagonalElements
   * @param convergenceCriterion
   * @param nMaxIterations
   * @param output
   * @param sigmaVectorCalculator  A function to calculate sigma vectors. See e.g. LRSCFSigmaVectorCalculator
   * @param nMaxSubspaceDimension  The maximum dimension of the subspace. Will perform a subspace collapse
   *                               if dimension of the subspace becomes larger than this threshold.
   * @brief                        Finds the nEigen lowest eigenpairs of a symmetric diagonal
   *                               dominant matrix using the Davidson-Liu algorithm.
   * @note                         This is an actual implementation of the IterativeEigenvalueSolver class.
   *                               Since this class is used for hermitian, as well as non-hermitian
   *                               eigenvalue solvers, it works with std::vector's of matrices for
   *                               left and right eigenvectors. Though we here only
   *                               have one set of eigenvectors, we have to work with std::vector's of matrices
   *                               in order to be able to use the IterativeEigenvalueSolver class.
   */
  DavidsonSolver(
      unsigned int nDimension,
      unsigned int nEigen,
      Eigen::VectorXd& diagonalElements,
      double convergenceCriterion,
      unsigned int nMaxIterations,
      std::function<Eigen::MatrixXd(
          Eigen::MatrixXd& guessVectors)> sigmaVectorCalculator,
      unsigned int nMaxSubspaceDimension,
      std::string loadPath = "",
      std::string filePath = "",
      std::string id = "");

  virtual ~DavidsonSolver() = default;

private:
  /*
   *
   * Input arguments not covered by IterativeEigenvalueSolver class
   *
   */

  /*
   * A function to calculate sigma vectors. The input and output arguments are std::vector<Eigen::MatrixXd >
   * of guess vectors or sigma vectors to be consistent with the IterativeEigenvalueSolver class and
   * the LRSCF SigmaVectorCalculator. Note that the length of this vector is one since only one set of guess
   * vectors is used here.
   */
  std::function<Eigen::MatrixXd(
            Eigen::MatrixXd& guessVectors)> _sigmaVectorCalculator;

  /*
   * The maximum dimension of the subspace. Will perform subspace collapse if reached.
   */
  const unsigned int _nMaxSubspaceDimension;

  /*
   *
   * Functions not covered by IterativeEigenvalueSolver class
   *
   */

  /*
   * Performs one iteration of the Davidson algorithm
   */
  void iterate() override final;

  /*
   * Post-processing function, mandatory for inheritance but not needed here.
   */
  void postProcessing() override final;

  /*
   * Computes the correction vectors using diagonal preconditioning.
   */
  void computeCorrectionVectors(
      Eigen::MatrixXd& correctionVectors,
      Eigen::MatrixXd& residualVectors);

  /*
   * Restart the Davidson algorithm using the current approximation for the true
   * eigenvectors as new starting vectors if the dimension of the subspace is larger
   * than nMaxSubspaceDimension or if no correction vector has a Schmidt norm above 1.0e-3.
   * This will of course decrease the speed of convergence.
   */
  void subspaceCollapse(
      Eigen::MatrixXd& correctionVectors,
      Eigen::MatrixXd& guessVectors,
      Eigen::MatrixXd& sigmaVectors,
      Eigen::MatrixXd& eigenvectors);

  /*
   * Update guess vectors and discard correction vectors with small Schmidt norms below 1.0e-3.
   */
  void updateGuessVectors(
      Eigen::MatrixXd& correctionVectors,
      Eigen::MatrixXd& guessVectors);

  /*
   * Set dimensions, print header and initialize sigma vectors
   */
  void initialize();


};

} /* namespace Serenity */

#endif /* BASICS_MATH_LINEARALGEBRA_DAVIDSONSOLVER_H_ */
