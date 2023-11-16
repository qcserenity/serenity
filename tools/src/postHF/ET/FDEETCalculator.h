/**
 * @file   FDEETCalculator.h
 *
 * @date   Apr. 20, 2020
 * @author Patrick Eschenbach
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
#ifndef FDEETCALCULATOR_H
#define FDEETCALCULATOR_H

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward Declarations */
class SystemController;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
/**
 * @class  FDEETCalculator FDEETCalculator.h
 * @brief  Stores and calculates all matrices needed in the FDEETTask. Is basically the "worker class".
 */
class FDEETCalculator {
 public:
  /**
   * @brief Constructor
   */
  FDEETCalculator();
  /**
   * @brief Default destructor.
   */
  virtual ~FDEETCalculator() = default;
  /**
   * @brief Calculates the transition MO overlap matrix for a given set of subsystems.
   * @param nStates Number of diabatic states
   * @param nSysPerState Number of subsystems each diabatic state is constructed from
   * @param jointIndex Indices of joined subsystems when disjoint approx. is used
   * @param allSystems Vector of all subsystems in the different states. First all systems of state 1 are listed
   *                    subsequently state 2 etc.
   * @param nOccState Information about total number of occupied orbitals in each state, separated by alpha and beta
   */
  void calculateMoOverlap(const unsigned nStates, const unsigned nSysPerState, const std::vector<unsigned> jointIndex,
                          const std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>>& stateCoeffs,
                          const std::vector<std::vector<std::shared_ptr<SystemController>>>& stateVector,
                          const std::shared_ptr<std::vector<unsigned>>& indexVector, const Eigen::MatrixXi nOccState);
  /**
   * @brief Calculates the Moore--Penrose pseudo inverse of blocks in a given matrix. Additionally calculates the
   * determinant of the transition overlap matrix blocks.
   * @param matrix The matrix whose blocks shall be inverted
   * @param nOccState Information about total number of occupied orbitals in each state, separated by alpha and beta
   * @param nStatesCouple Number of diabatic states that are coupled
   * @param invThreshold Threshold below which singular values are not inverted but set to 0
   */
  void calcBlockedMoorePenrose(const std::vector<Eigen::MatrixXd>& matrix, Eigen::MatrixXi nOccState,
                               unsigned nStatesCouple, double invThreshold);
  /**
   * @brief Calculates the transition electron density matrices of a given block matrix of pseudo inverses.
   * @param pseudoInverse A block matrix that contains pseudo inverses of the transition overlap matrix as blocks
   * @param superSystem The supersystem constructed from all subsystems. Only geometry and settings are necessary.
   *                     This supersystem shall not have any electronic structure
   * @param coeffs Stores the subsystems coefficients
   * @param nOccState Information about total number of occupied orbitals in each state, separated by alpha and beta
   * @param nStatesCouple Number of diabatic states that are coupled
   */
  void calcTransDensMats(const std::vector<Eigen::MatrixXd>& pseudoInverse, std::shared_ptr<SystemController> superSystem,
                         const std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>>& coeffs,
                         Eigen::MatrixXi nOccState, unsigned nStatesCouple);
  /**
   * @brief Calculates the transition electron density matrices of a given block matrix of pseudo inverses and writes
   *        them a HDF5 file.
   * @param pseudoInverse A block matrix that contains pseudo inverses of the transition overlap matrix as blocks
   * @param superSystem The supersystem constructed from all subsystems. Only geometry and settings are necessary.
   *                     This supersystem shall not have any electronic structure
   * @param coeffs Stores the subsystems coefficients
   * @param nOccState Information about total number of occupied orbitals in each state, separated by alpha and beta
   * @param nStatesCouple Number of diabatic states that are coupled
   */
  void calcTransDensMatsDisk(const std::vector<Eigen::MatrixXd>& pseudoInverse, std::shared_ptr<SystemController> superSystem,
                             const std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>>& coeffs,
                             Eigen::MatrixXi nOccState, unsigned nStatesCouple);
  /**
   * @brief Calculates the Hamilton matrix in the basis of given diabatic states. It therefore uses the transition
   *        electron density matrices
   * @param transDensMats Contains the transition electron density matrices
   * @param superSystem The supersystem constructed from all subsystems. Only geometry and settings are necessary.
   *                     This supersystem shall not have any electronic structure
   * @param nStatesCouple Number of diabatic states that are coupled
   */
  void calcHamiltonian(const std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>>& transDensMats,
                       std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple);
  /**
   * @brief Calculates the Hamilton matrix in the basis of given diabatic states in disk mode. It therefore uses the
   * transition electron density matrices
   * @param transDensMats Contains the transition electron density matrices
   * @param superSystem The supersystem constructed from all subsystems. Only geometry and settings are necessary.
   *                     This supersystem shall not have any electronic structure
   * @param nStatesCouple Number of diabatic states that are coupled
   */
  void calcHamiltonianDisk(std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple);
  /**
   * @brief Calculates the Hamilton matrix element in the basis of given diabatic states. It therefore uses the
   * transition electron density matrices. This routine is called by calcHamiltonian()
   * @param dMatC Contains the transition electron density matrices
   * @param superSystem The supersystem constructed from all subsystems. Only geometry and settings are necessary.
   *                     This supersystem shall not have any electronic structure
   */
  double calcDFTEnergy(std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatC,
                       std::shared_ptr<SystemController> superSystem);
  /**
   * @brief Solves a generalized eigenvalue problem in the basis of the diabatic states. Therefore uses the Hamilton
   * matrix and overlap of diabatic states (determinants).
   * @param H Hamilton matrix
   * @param determinants Overlap of diabatic states
   * @param nStatesCouple Number of diabatic states that are coupled
   */
  void solveEigenValueProblem(Eigen::MatrixXd H, Eigen::MatrixXd determinants, unsigned nStatesCouple);
  /**
   * @brief Prints the calculated transition electron density matrices into a cube file.
   * @param superSystem The supersystem constructed from all subsystems. Only geometry and settings are necessary.
   *                     This supersystem shall not have any electronic structure
   * @param nStatesCouple Number of diabatic states that are coupled
   */
  void printTransDensContributions(std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple);
  /**
   * @brief Prints the results of a FDE-ET/Diab calculation.
   * @param nStatesCouple Number of diabatic states that are coupled
   */
  void printResults(unsigned nStatesCouple);

  /**
   * Getter and Setter
   */

  /**
   * @return The transition MO overlap matrix
   */
  std::vector<Eigen::MatrixXd> getMoOverlap();
  /**
   * @brief Sets the transition MO overlap matrix
   * @param overlap A desired matrix
   */
  void setMoOverlap(std::vector<Eigen::MatrixXd> overlap);
  /**
   * @return The Moore--Penrose pseudo inverse matrix
   */
  std::vector<Eigen::MatrixXd> getPseudoInverse();
  /**
   * @brief Sets the Moore--Penrose pseudo inverse matrix
   * @param pseudo A desired matrix
   */
  void setPseudoInverse(std::vector<Eigen::MatrixXd> pseudo);
  /**
   * @return The overlap matrix of diabatic states
   */
  Eigen::MatrixXd getDeterminants();
  /**
   * @brief Sets the overlap matrix of diabatic states
   * @param dets A desired matrix
   */
  void setDeterminants(Eigen::MatrixXd dets);
  /**
   * @return The transition electron density matrix
   */
  std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> getTransDensMats();
  /**
   * @brief Sets the transition electron density matrix
   * @param mat A desired matrix
   */
  void setTransDensMats(std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> mat);
  /**
   * @return The Hamilton matrix
   */
  Eigen::MatrixXd getHamiltonian();
  /**
   * @brief Sets the Hamilton matrix
   * @param newH A desired matrix
   */
  void setHamiltonian(Eigen::MatrixXd newH);
  /**
   * @return The eigenvalues
   */
  Eigen::VectorXd getEigenValues();
  /**
   * @brief Sets the eigenvalues
   * @param newValues A desired matrix
   */
  void setEigenValues(Eigen::VectorXd newValues);
  /**
   * @return The linear combination coefficients
   */
  Eigen::MatrixXd getLinCoeffs();
  /**
   * @brief Sets the linear combination coefficients
   * @param newCoeffs A desired matrix
   */
  void setLinCoeffs(Eigen::MatrixXd newCoeffs);
  /**
   * @return The analytical electronic coupling
   */
  double getAnalyticalCoupling();
  /**
   * @brief Sets the analytical electronic coupling
   * @param newC A desired value
   */
  void setAnalyticalCoupling(double newC);
  /**
   * @return The excitation energy
   */
  double getLongRangeExEnergy();
  /**
   * @brief Sets the excitation energy
   * @param newE A desired value
   */
  void setLongRangeExEnergy(double newE);

  std::shared_ptr<std::vector<std::string>> getDensityMatrixFiles();

 private:
  ///@brief Transition MO overlap matrix
  std::vector<Eigen::MatrixXd> _overlapMO;
  ///@brief Moore-Penrose pseudo inverse
  std::vector<Eigen::MatrixXd> _pseudoInverse;
  ///@brief Overlap of diabatic states
  Eigen::MatrixXd _determinants;
  ///@brief Transition electron density matrices
  std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> _densMats;
  ///@brief Hamilton matrix
  Eigen::MatrixXd _hamiltonian;
  ///@brief Linear combination coefficients
  Eigen::MatrixXd _linCoeffs;
  ///@brief Eigenvalues of adiabatic states
  Eigen::VectorXd _eigenValues;
  ///@brief Analytical electronic coupling
  double _analyticalCoupling;
  ///@brief Excitation energy
  double _longRangeExEnergy;
  ///@brief densitymatrix files when disk mode is used
  std::shared_ptr<std::vector<std::string>> _dMatFiles;
};
} /* namespace Serenity */
#endif /* FDEETCALCULATOR_H */