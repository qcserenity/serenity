/**
 * @file   FDEDiabController.h
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
#ifndef FDEDIABCONTROLLER_H
#define FDEDIABCONTROLLER_H

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward Declarations */
class SystemController;

/**
 * @class  FDEDiabController FDEDiabController.h
 * @brief  Stores and calculates adiabatic wavefunctions and (spin)-densities after FDE-ET run.
 */
class FDEDiabController {
 public:
  /**
   * @brief Constructor
   * @param transDensMats Transition electron density matrices
   * @param linCoeffs Linear combination coefficients for eigenvalue problem in FDE-ET
   * @param determinants Overlap of diabatic states
   * @param superSystem The supersystem constructed from all subsystems. Only geometry and settings are necessary.
   *                     This supersystem shall not have any electronic structure
   * @param nStatesCouple Number of diabatic states that are coupled
   */
  FDEDiabController(const std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>>& transDensMats,
                    const Eigen::MatrixXd& linCoeffs, const Eigen::MatrixXd& determinants,
                    std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple, unsigned nStatesAdiab);

  FDEDiabController(std::shared_ptr<std::vector<std::string>> densityMatrixFiles, const Eigen::MatrixXd& linCoeffs,
                    const Eigen::MatrixXd& determinants, std::shared_ptr<SystemController> superSystem,
                    unsigned nStatesCouple, unsigned nStatesAdiab);
  /**
   * @brief Default destructor
   */
  virtual ~FDEDiabController() = default;
  /**
   * @return Vector of adiabatic electron density matrices
   */
  std::vector<std::shared_ptr<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> getAdiabDensities();
  /**
   * @brief Calculates the vector of adiabatic electron density matrices
   */
  void calcAdiabDensities();
  /**
   * @brief Prints all adiabatic electron density matrices to a cube file
   */
  void printAdiabDensities();
  /**
   * @return Vector of atom populations for each adiabatic state
   */
  std::vector<Eigen::VectorXd> getSortedAtomPopulations();
  /**
   * @brief Calculates vector of atomic spin populations for each adiabatic state
   * @param list A list of wanted adiabatic states
   */
  void calcPopulations(std::vector<unsigned int> list);
  /**
   * @brief Prints vector of atomic spin populations for each adiabatic state to output
   * @param list A list of wanted adiabatic states
   */
  void printPopulations(std::vector<unsigned int> list);
  /**
   * @brief Calculates the <S*S> value for a given adiabatic wave function according to
   *    [1] J. Wang, A. D. Becke, V.H. Smith. Evaluation of hS2i in restricted, unrestricted
   *        Hartree–Fock, and density functional based theories. J. Chem. Phys., 102(8) (1995)
   *        3477–3480.
   * @param iState Index of adiabatic state
   */
  double getS2(unsigned iState);
  /**
   * @brief Performs Becke Population analysis for adiabatic spin density and calculates <S*S> for each adiabatic state.
   *        Each of those quantities is printed to the output
   * @param list A list of wanted adiabatic states
   */
  void adiabaticWavefunctionAnalysis(std::vector<unsigned int> list);

 private:
  /// @brief Sets vector of adiabatic electron density matrices
  void setAdiabDensities(std::vector<std::shared_ptr<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> newD);
  /// @brief Sets vector of atom populations for each diabatic state
  void setSortedAtomPopulations(std::vector<Eigen::VectorXd> newPop);
  /// @brief The density matrices of the adiabatic wave functions
  std::vector<std::shared_ptr<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> _adiabDensities;
  /// @brief Transition electron density matrices
  const std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> _transDensMats;
  /// @brief Linear combination coefficients of states
  Eigen::MatrixXd _linCoeffs;
  /// @brief Overlap of diabatic states
  Eigen::MatrixXd _determinants;
  /// @brief The supersystem constructed from all subsystems
  std::shared_ptr<SystemController> _superSystem;
  /// @brief Number of diabatic states
  unsigned _nStates;
  /// @brief Number of adiabatic states
  unsigned _nAdiab;
  /// @brief Vector of atomic spin populations for each adiabatic state
  std::vector<Eigen::VectorXd> _sortedAtomPopulations;
  ///@brief densitymatrix files when disk mode is used
  std::shared_ptr<std::vector<std::string>> _densityMatrixFiles;
};
} /* namespace Serenity */
#endif /* FDEDIABCONTROLLER_H */