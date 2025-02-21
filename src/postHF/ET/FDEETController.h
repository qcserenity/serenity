/**
 * @file   FDEETController.h
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
#ifndef FDEETCONTROLLER_H
#define FDEETCONTROLLER_H

/* Include Serenity Internal Headers */
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward Declarations */
class SystemController;

/**
 * @class  FDEETController FDEETController.h
 * @brief  Stores and calculates all necessary dimensions and information needed in the ElectronTransferTask. Is
 * basically the "information class"
 */
class FDEETController {
 public:
  /**
   * @brief Constructor.
   * @param activeSystems  The active systems.
   * @param nStatesInput Total number of states specified in the input
   * @param nSysPerState Number of subsystems each diabatic state is constructed from
   */
  FDEETController(const std::vector<std::shared_ptr<SystemController>>& activeSystems, unsigned nStatesInput,
                  unsigned nSysPerState);
  /**
   * @brief Default destructor
   */
  virtual ~FDEETController() = default;
  /**
   * @brief Stores the information, which states in the input shall be used for coupling
   * @param couples List of diabatic states used in this coupling step
   */
  void setupIndexVector(std::vector<unsigned> couples);
  /**
   * @brief Provides the information, which states in the input shall be used for coupling
   * @return A vector containing the indices of the systemControllers to be used
   */
  std::shared_ptr<std::vector<unsigned>> getIndexVector();
  /**
   * @brief Set the information, which states in the input shall be used for coupling
   * @param states Set of indices for systemControllers to be used
   */
  void setIndexVector(std::shared_ptr<std::vector<unsigned>> states);
  /**
   * @brief Stores the information of occupied orbitals and basis functions for each diabatic state
   */
  void setupOrbitalInfo();
  /**
   * @return Information of occupied orbitals for each diabatic state
   */
  Eigen::MatrixXi getNOccState();
  /**
   * @brief Sets new nOccState vector
   * @param newOcc The new matrix
   */
  void setNOccState(Eigen::MatrixXi newOcc);
  /**
   * @return Information of basis functions for each diabatic state
   */
  Eigen::VectorXi getNBfState();
  /**
   * @brief Sets new nBfState vector
   * @param newBf The new matrix
   */
  void setNBfState(Eigen::VectorXi newBf);
  /**
   * @brief Provides the stateVector that stores all systemsControllers according to their diabatic state
   * @return The stateVector
   */
  std::vector<std::vector<std::shared_ptr<SystemController>>> getStateVector();
  /**
   * @brief Set the stateVector
   * @param newV A new stateVector
   */
  void setStateVector(std::vector<std::vector<std::shared_ptr<SystemController>>> newV);
  /**
   * @brief Provides all subsystems used in this coupling problem
   * @return all subsystems used in this coupling problem
   */
  std::vector<std::shared_ptr<SystemController>> getSubsystems();
  /**
   * @brief Set the subsystems
   * @param systems A new set of systemControllers
   */
  void setSubsystems(std::vector<std::shared_ptr<SystemController>> systems);
  /**
   * @brief Returns the subsystem coefficient matrices of each diabatic state
   * @return The blocked subsystem coefficient matrices of each diabatic state
   */
  std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>> getStateCoeffs();
  /**
   * @brief Returns the number of diabatic states to be coupled
   * @return The number of diabatic states to be coupled
   */
  unsigned getNStatesCouple();
  /**
   * @brief Sets the number of diabatic states to be coupled
   * @param nStates number of diabatic states to be coupled
   */
  void setNStatesCouple(unsigned nStates);
  /**
   * @brief Returns the total number of states specified in the input
   * @return Total number of states specified in the input
   */
  unsigned getNStatesInput();
  /**
   * @brief Returns the number of subsystems used to construct a diabatic state
   * @return The number of subsystems used to construct a diabatic state
   */
  unsigned getNSysPerState();
  /**
   * @brief Returns the systemControllers of systems used in this coupling problem
   * @return The systemControllers of systems used in this coupling problem
   */
  const std::vector<std::shared_ptr<SystemController>>& getActiveSystems();

 private:
  /// @brief Collects all systems used in this coupling problem and puts them into a vector
  void collectSubsystems();
  /// @brief Constructs vector from input that contains all systemControllers from the different states
  void setupStateVector();
  /// @brief Construct the blocked subsystem coefficient matrices for each diabatic state
  void setupStateCoeffs();
  /// @brief systemControllers of systems used in this coupling problem
  const std::vector<std::shared_ptr<SystemController>>& _activeSystems;
  /// @brief Contains indices which systemControllers will be used in this coupling problem
  std::shared_ptr<std::vector<unsigned>> _indexVector;
  /// @brief Contains all systemControllers from the different states
  std::vector<std::vector<std::shared_ptr<SystemController>>> _stateVector;
  /// @brief Systems used for coefficients construction
  std::vector<std::shared_ptr<SystemController>> _allSystems;
  /// @brief Blocked subsystem coefficient matrices for each diabatic state
  std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>> _stateCoeffs;
  /// @brief Total number of states specified in the input
  unsigned _nStatesInput;
  /// @brief Number of subsytems in each diabatic state
  unsigned _nSysPerState;
  /// @brief Number of states to couple
  unsigned _nStatesCouple;
  /// @brief Information about total number of occupied orbitals in each state, separated by alpha and beta
  Eigen::MatrixXi _nOccState;
  /// @brief Information about total number of basis functions in each state, separated by alpha and beta
  Eigen::VectorXi _nBfState;
};
} /* namespace Serenity */
#endif /* FDEETCONTROLLER_H */