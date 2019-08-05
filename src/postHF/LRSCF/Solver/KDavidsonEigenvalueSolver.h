/**
 * @file   KDavidsonEigenvalueSolver.h
 *
 * @date   Mar 4, 2017
 * @author M. Boeckers
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

#ifndef ELECTRONICSTRUCTURECALCULATIONS_POSTHF_LRSCF_KDAVIDSONEIGENVALUESOLVER_H_
#define ELECTRONICSTRUCTURECALCULATIONS_POSTHF_LRSCF_KDAVIDSONEIGENVALUESOLVER_H_

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Solver/IterativeEigenvalueSolver.h"
#include "math/Matrix.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <functional>



namespace Serenity {
using namespace Eigen;
using namespace std;
/**
 * @brief Implementation of the K-Davidson iterative eigenvalue solver to obtain the k lowest
 *        eigenvalues of a linear response SCF problem. This is an adapted Davidson
 *        eigenvalue solver which uses the fact that the matrix \f$ (A+B)(A-B) = MK \f$ is symmetric
 *        (self-adjoint) with respect to the K-inner product
 *        \f[
 *          \left< v, MK \right>_K = v^T KMK v = (MKv)^T K v = \left< MK v, v \right>_K
 *        \f]
 *        where \f$ \left< X, Y \right>_K = X^T K Y \f$. Using the K-inner product and the fact that
 *        there is a relation between the X- and Y- component of the linear response eigenvalue problem,
 *        one only has to optimize one vector (either X or Y) and can produce the other vector in a
 *        postprocessing procedure.
 *        In contrast to the original Davidson(-Liu) algorithm, this algorithm only requires 2k
 *        multiplications of matrices with vectors per iteration instead of 4k multiplications.
 *        Thus, it should be two times faster than the original linear response Davidson solver.
 *        Furthermore, the subspace is constructed with only one set of the preconditioned residual
 *        vectors instead of two sets of vectors.
 *        However, here we cannot calculate the MKS sigma vectors in one step, which is different from
 *        other eigenvalue solvers where we can calculate the M and K sigma vectors simultaneously.
 *        Thus, when exact exchange is used, this eigenvalue solver is more expensive than the OJJ
 *        or SSF solver since the M matrix is no longer diagonal.
 *        This implementation follows E. Vecharynski et. al, JCC (submitted 2015). Some names of the
 *        variables are chosen according to the pseudocode of algorithm 5 in that report.
 *        This iterative eigenvalue solver has numerical problems when the convergence threshold
 *        is too low. A convergence threshold of 1.0e-6 is readily obtained, which should be sufficient
 *        for most applications. Convergence below that is not guaranteed. The reason for this seems
 *        to be the fact that the eigenvalues are calculated for the non-symmetric matrix (A+B)(A-B).
 *        Eigenvalue solvers for the symmetric LRSCF problem seem to be more stable.
 *        (see. Dissertation of F. Furche: Dichtefunktionalmethoden fuer angeregte Molekuele - Theorie,
 *        Implementierung und Anwendung, Cuvillier Verlag Goettingen, 2002)
 *
 *
 * @note  List of variables to understand this class:
 *
 *          * _guessVectors[0]: A test set S of guess vectors.
 *          * _sigmaVectors[0]: Sigma vectors (A-B)S
 *          * _sigmaVectors[1]: Sigma vectors (A+B)(A-B)S
 *          * _eigenvectors[0]: X part of the eigenvector
 *          * _eigenvectors[1]: Y part of the eigenvector
 */
class KDavidsonEigenvalueSolver: public IterativeEigenvalueSolver {
public:
  KDavidsonEigenvalueSolver(
      unsigned int nDimension,
      unsigned int nEigen,
      Eigen::VectorXd& diagonalElements,
      double convergenceCriterion,
      unsigned int nMaxIterations,
      std::function<Eigen::MatrixXd(Eigen::MatrixXd& guessVectors)> MSigmaCalculator,
      std::function<Eigen::MatrixXd(Eigen::MatrixXd& guessVectors)> KSigmaCalculator,
      unsigned int maxSubspaceDimension,
      std::string loadPath = "",
      std::string filePath = "",
      std::string id = "");

  virtual ~KDavidsonEigenvalueSolver() = default;

  /**
   *
   * @return Returns X
   */
  Eigen::MatrixXd& getX();
  /**
   *
   * @return Returns Y
   */
  Eigen::MatrixXd& getY();

private:
  /*
   *
   * Input variables not covered by IterativeEigenvalueSolver class
   */

  /*
   * A function to calculate the sigma vectors. Note that in and output arguments are
   * std::vector's to be consistent with the LRSCFSigmaVectorCalculator and
   * the IterativeEigenvalueSolver classes.
   */
  std::function<
      std::vector<Eigen::MatrixXd > (
          std::vector<Eigen::MatrixXd >& guessVector,
          std::vector<Options::SIGMAVECTOR_METHODS> method)> _sigmaVectorCalculator;


  /*
   * Functions to calculate the M (A+B) and K (A-B) sigma vectors
   */
  std::function<Eigen::MatrixXd(Eigen::MatrixXd& guessVector)> _MSigmaCalculator;
  std::function<Eigen::MatrixXd(Eigen::MatrixXd& guessVector)> _KSigmaCalculator;

  /*
   * Maximum subspace dimension
   */
  unsigned int _maxSubspaceDimension;

  /*
   *
   * Other member variables not covered by IterativeEigenvalueSolver class
   *
   */


  /*
   * Current eigenvectors of matrix MK. Needed for subspace collapse.
   */
  Eigen::MatrixXd _YM;

  /*
   * The method for which the sigma vector will  be calculated. Note that this is a std::vector
   * to be consistent with the LRSCFSigmaVectorCalculator, though only one element is needed.
   */
  std::vector<Options::SIGMAVECTOR_METHODS> _method;


  /*
   * A set of current correction vectors to be k-orthogonalized and appended to the set of guess vectors S
   */
  Eigen::MatrixXd _W;

  /*
   * The current dimension of the subspace eigenvalue problem
   */
  unsigned int _subspaceDimension;

  /*
   *
   * Member functions
   *
   */

  /*
   * Perform one iteration
   */
  void iterate() override final ;

  /*
   * Postprocessing: Obtain second component of the eigenvector
   */
  void postProcessing() override final;


  //Calculates correction vectors given the residual vectors and diagonal elements.
  //Currently, only the diagonal preconditioner is implemented.
  void computeCorrectionVectors(
      Eigen::MatrixXd& residualVectors,
      Eigen::MatrixXd& correctionVectors,
      Eigen::VectorXd& diagonalElements);


  /*
   * K-orthogonalizes S correction vectors iteratively using the SVQB algorithm of
   * A. Stathopoulos and K. Wu (SIAM J. Sci. Comput. 23, 2165 (2002).
   * In every iteration, one has to ensure that the correction vectors are k-orthonormal,
   * i.e.< S, S >_K = I .
   */
  void kOrthogonalize(
      Eigen::MatrixXd& S,
      Eigen::MatrixXd& KS,
      Eigen::MatrixXd& MKS);

  /*
   * Actions to be performed before the actual iterations begin
   */
  void initializeKDavidson();

  /*
   *  Collapses the subspace.
   */
  void subspaceCollapse();

};

} /* namespace Serenity */

#endif /* ELECTRONICSTRUCTURECALCULATIONS_POSTHF_LRSCF_KDAVIDSON_H_ */

