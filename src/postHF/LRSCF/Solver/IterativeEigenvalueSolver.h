/**
 * @file IterativeEigenvalueSolver.h
 *
 * @date Oct 17, 2016
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

#ifndef BASICS_MATH_LINEARALGEBRA_ITERATIVEEIGENVALUESOLVER_H_
#define BASICS_MATH_LINEARALGEBRA_ITERATIVEEIGENVALUESOLVER_H_

/* Include Serenity Internal Headers */
#include "io/FormattedOutput.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <iomanip>
#include <vector>


namespace Serenity {


/**
 * @brief Class holding objects for iterative eigenvalue solvers like the Davidson algorithm.
 *        This class is written in a very general fashion, all complex mathematics should be
 *        handled by actual implementations of iterative eigenvalue solvers. It is used for
 *        implementations of iterative eigenvalue solvers for hermitian as well as non-hermitian
 *        eigenvalue problems.
 */
class IterativeEigenvalueSolver {
public:
  /**
   *
   * @param nDimension              The dimension of the eigenvalue problem
   * @param nEigen                  The number of eigenpairs to be determined
   * @param diagonalElements        The diagonal elements of the matrix to be diagonalized
   * @param convergenceCriterion    The convergence threshold
   * @param nMaxIterations          The maximum number of iterations
   */
  IterativeEigenvalueSolver(
      unsigned int nDimension,
      unsigned int nEigen,
      Eigen::VectorXd& diagonalElements,
      double convergenceCriterion,
      unsigned int nMaxIterations,
      std::string loadPath = "",
      std::string filePath = "",
      std::string id = "");

  virtual ~IterativeEigenvalueSolver() = default;

  /**
   * @return Returns the desired eigenvectors.
   */
  std::vector<Eigen::MatrixXd >& getEigenvectors();

  /**
   * @return Returns the desired eigenvalues.
   */
  Eigen::VectorXd& getEigenvalues();

  /**
   * @return Returns the residual vectors
   */
  std::vector<Eigen::MatrixXd >& getResiduals();

  /**
   * @brief Solve eigenvalue problem.
   * @note  This function calls iterate() until it sets _nConverged to _nEigen, e.g.
   *        by calling checkForConvergence()
   */
  void solve();

  /**
   * @return Returns true if eigenvectors are converged
   */
  bool& isConverged();

protected:

  /*
   *
   * Input variables
   *
   */

  /*
   * The dimension of the eigenvalue problem
   */
  unsigned int _nDimension;

  /*
   * The number of eigenvalues to be determined. Note that it is not const unsigned int since one may
   * want to perform a pre-optimization where more eigenvectors are optimized.
   */
  unsigned int _nEigen;

  /*
   * A Vector holding the diagonal elements of the matrix for which the eigenvalues shall be computed.
   */
  Eigen::VectorXd _diagonalElements;

  /*
   * The convergence threshold. Note that it is not const double since one may want to perform a
   * pre-optimization until a higher convergence threshold.
   */
  double _convergenceCriterion;

  /*
   * The maximum number of iterations. Note that it is not const unsigned int since one may want to
   * perform a pre-optimization.
   */
  unsigned int _nMaxIterations;

  /*
   * Paths and ID for restart from file
   */
  std::string _loadPath;
  std::string _filePath;
  std::string _id;

  /*
   * True if restart information shall be written to file. Is set to true if filePath is different from "".
   */
  bool _restart;

  /*
   *
   * Further variables
   *
   */

  /*
   * A vector holding test eigenvector sets. Most eigenvalue solvers use only one set of test vectors.
   * However, some solvers for non-hermitian eigenvalue problems use two sets of test vectors for both
   * left and right eigenvectors.
   */
  std::vector<Eigen::MatrixXd > _guessVectors;

  /*
   * A vector holding the sigma vectors, i.e. the product of the matrix for which the eigenvalues shall
   * be calculated with the test vectors.
   */
  std::vector<Eigen::MatrixXd > _sigmaVectors;

  /*
   * A vector holding the residual vectors
   */
  std::vector<Eigen::MatrixXd > _residualVectors;

  /*
   * A vector holding the eigenvectors corresponding to the set of test vectors.
   */
  std::vector<Eigen::MatrixXd > _eigenvectors;

  /*
   * A vector holding the eigenvalues
   */
  Eigen::VectorXd _eigenvalues;

  /*
   * A vector holding the reduced matrices, i.e. the projection of the matrix for which the eigenvalues
   * shall be computed into a subspace spanned by the sets of guess vectors. Note that also for non-
   * hermitian matrices only one subspace matrix is used (which could be constructed from two sets of
   * guess vectors). I do not know any algorithm that would need two subspace matrices, at
   * least not for the most common eigenvalue problems in quantum chemistry.
   */
  Eigen::MatrixXd _subspaceMatrix;

  /*
   * A matrix holding the eigenvectors of the reduced eigenvalue problem, i.e. the expansion coefficients
   * for the sets of test vectors.
   */
  Eigen::MatrixXd _subspaceEigenvectors;

  /*
   * A vector holding the eigenvalues of the reduced eigenvalue problem, i.e. an approximation to the
   * true eigenvalues.
   */
  Eigen::VectorXd _subspaceEigenvalues;

  /*
   * The number of converged eigenvalues. Must be set to _nEigen if converged
   */
  unsigned int _nConverged = 0;

  /*
   * True if converged. Note that this variable is not private since one should have the
   * possibility to override the solve function.
   */
  bool _converged = false;

  /*
   * True if eigenpairs have been calculated. Note that this variable is not declared private since
   * one should have the possibility to override the solve function.
   */
  bool _hasBeenRun = false;

  /*
   * The number of iterations. Note that this variable is not declared private since one should have
   * the possibility to override the solve function.
   */
  unsigned int _nIter = 0;

  /*
   *
   * Common functions used by most eigenvalue solvers. Most of these functions are written in such
   * a way that they take a matrix as argument instead of a std::vector<Eigen::MatrixXd > because it
   * is not a priori clear how to handle both (or more) sets of guess vectors.
   *
   */
  /*
   * Perform one iteration. This function must be overridden. If convergence is achieved, it has to
   * set the member variable _nConverged to _nEigen.
   */
  virtual void iterate() = 0;

  /*
   * Post processing tasks, e.g. transformation of converged eigenvectors, etc.  This function must
   * be overridden in an actual implementation.
   */
  virtual void postProcessing() = 0;

  /*
   * Calculates unit vectors with ones in the position corresponding to the smallest diagonal elements.
   * Note, that the test vectors and diagonal elements are given as arguments to have the possibility
   * to work with edited diagonal elements and guess vectors.
   */
  void computeInitialGuessVectors(
      Eigen::MatrixXd& guessVectors,
      Eigen::VectorXd& diagonalElements);

  /*
   * Calculate the current approximation to the eigenvectors. The expansion coefficients are commonly given
   * by the subspace eigenvectors. Note that the trial vectors, expansion coefficients and guess vectors
   * are given as arguments to be able to edit them before calculating the current approximation.
   */
  void computeTrialVectors(
    Eigen::MatrixXd& trialVectors,
    Eigen::MatrixXd& expansionCoefficients,
    Eigen::MatrixXd& guessVectors);

  /*
   * Compute residual vectors. The residuals are calculated as
   *
   * R_i = \sum_j c_j * (\sigma_ij - \omega_i B_ij),
   *
   * where \sigma is the sigma vector and B the test vector matrix.
   * This function can (in principle) be used for both, hermitian and non-hermitian eigenvalue
   * problems.
   * Note that for non-hermitian eigenvalue problems the the guess vectors used to calculate
   * the sigma vectors are different from those needed here. See e.g. equation 27 and 28 in
   * J. Chem. Phys. 109, 8218 (1998), which is common for non-hermitian iterative eigenvalue
   * solvers.
   */
  void computeResidualVectors(
      Eigen::MatrixXd& residualVectors,
      Eigen::MatrixXd& sigmaVectors,
      Eigen::MatrixXd& guessVectors,
      Eigen::MatrixXd& expansionCoefficients,
      Eigen::VectorXd& subspaceEigenvalues);


  /*
   * Checks for convergence
   */
  void checkForConvergence(std::vector<Eigen::MatrixXd >& residualVectors);



  /*
   * Print some formatted data during the iterations. This function takes a std::vector of residual vectors
   * for both, left and right eigenvectors (order and length of the vector do not matter, it will only print
   * the maximum norm of all residuals.). The eigenvalues are specified as argument explicitly to have the
   * possibility of editing them befor printing.
   */
  void printIterationData(
      std::vector<Eigen::MatrixXd >& residualVectors,
      Eigen::VectorXd& eigenvalues);

  /*
   * Prints the caption of the eigenvalue solver.
   */
  void printHeader(std::string caption);

  /*
   * Writes guess vectors and sigma vectors to file for restart
   */
  void toHDF5();

  /*
   * Returns true if initialization from file was successful.
   */
  bool initializeFromHDF5();
};

} /* namespace Serenity */

#endif /* BASICS_MATH_LINEARALGEBRA_ITERATIVEEIGENVALUESOLVER_H_ */
