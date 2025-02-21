/**
 * @file NTOCalculator.h
 *
 * @date Aug 26, 2019
 * @author Johannes Toelle
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

#ifndef NTOCALCULATOR_H_
#define NTOCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/matrices/MatrixInBasis.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

/* Forward declarations */
class SystemController;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class  NTOCalculator NTOCalculator.h
 * @brief  Calculates natural-transition-orbitals from linear response TDDFT/TDA/CIS/TDHF calculation
 *
 * The theory can be found in: J. Chem. Phys, 118, 4775 (2013) "Natural transition orbitals"
 * Hole and particle densities can be found in: J. Herbert, Phys. Chem. Chem. Phys., 26, 3755 (2024)
 */
template<Options::SCF_MODES SCFMode>
class NTOCalculator {
 public:
  /**
   * @brief Constructor.
   * @param activeSystems      A vector containing all active systems (in the form of SystemControllers).
   * @param environmentSystems A vector containing all environment systems
   *                           (in the form of SystemControllers).
   * @param plottingThreshold  Threshold which NTO contributes to a particular excitation
   */
  NTOCalculator(std::vector<std::shared_ptr<SystemController>> activeSystems,
                std::vector<std::shared_ptr<SystemController>> environmentSystems, const double plottingThreshold);
  /**
   * @brief Default destructor.
   */
  virtual ~NTOCalculator() = default;

  /// @returns the number of excited states
  unsigned int getNumberOfStates();
  /// @returns the folder in which the NTO of i gets written
  std::string getDir(unsigned int i, unsigned int iSubsystem);
  /// @returns the NTO matrix for the occupied orbitals of iState
  const SpinPolarizedData<SCFMode, Eigen::MatrixXd>& getOccNTOs(unsigned int iState, unsigned int iSubsystem) {
    makeAvailable(iState);
    return _occNTOs[iSubsystem];
  }
  /// @returns the NTO matrix for the virtual orbitals of iState
  const SpinPolarizedData<SCFMode, Eigen::MatrixXd>& getVirtNTOs(unsigned int iState, unsigned int iSubsystem) {
    makeAvailable(iState);
    return _virtNTOs[iSubsystem];
  }
  /// @returns the eigenvalues of the occupied NTO matrix of iState
  const SpinPolarizedData<SCFMode, Eigen::VectorXd>& getOccEigenvalues(unsigned int iState) {
    makeAvailable(iState);
    return _occEigenvalues;
  }
  /// @returns the eigenvalues of the virtual NTO matrix of iState
  const SpinPolarizedData<SCFMode, Eigen::VectorXd>& getVirtEigenvalues(unsigned int iState) {
    makeAvailable(iState);
    return _virtEigenvalues;
  }

  /// @returns the transition density for a specific state and subsystem in the AO basis of that subsystem
  const MatrixInBasis<SCFMode>& getTransitionDensity(unsigned iState, unsigned int iSubsystem) {
    makeAvailable(iState);
    return _transitionDensity[iSubsystem];
  }

  /// @returns the hole density for a specific state and subsystem in the AO basis of that subsystem
  const MatrixInBasis<SCFMode>& getHoleDensity(unsigned iState, unsigned int iSubsystem) {
    makeAvailable(iState);
    return _holeDensity[iSubsystem];
  }

  /// @returns the particle density for a specific state and subsystem in the AO basis of that subsystem
  const MatrixInBasis<SCFMode>& getParticleDensity(unsigned iState, unsigned int iSubsystem) {
    makeAvailable(iState);
    return _particleDensity[iSubsystem];
  }

  /// @returns the hole density correction for a specific state and subsystem in the AO basis of that subsystem. This is
  /// only relevant (non-zero) for FDEc.
  const MatrixInBasis<SCFMode>& getHoleDensityCorrection(unsigned iState, unsigned int iSubsystem) {
    makeAvailable(iState);
    return _holeDensityCorrection[iSubsystem];
  }

 private:
  /**
   * @brief prints information of a particular occupied and virtual NTO of state
   *
   * @param spin Additional output argument to distinguish between alpha/beta/RESTRICTED
   * @param iState defines the particular excited state
   * @param vEigenvalues eigenvalues of the virtual NTO matrix
   * @param oEigenvalues eigenvalues of the occupied NTO matrix
   * @param u eigenvectors of the occupied NTO matrix
   * @param v eigenvectors of the virtual NTO matrix
   */
  void printNTOInfo(std::string spin, const unsigned int iState, Eigen::VectorXd& oEigenvalues,
                    Eigen::VectorXd& vEigenvalues, Eigen::MatrixXd& u, Eigen::MatrixXd& v);
  // The active system
  std::vector<std::shared_ptr<SystemController>> _activeSystems;
  // The environment systems
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  // The NTO plotting threshold
  const double _ntoThreshold;
  // The lrscfController
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;
  // The combined excitation vector X+Y
  Eigen::MatrixXd _XPY;
  // The combined excitation vector X-Y
  Eigen::MatrixXd _XMY;
  // The excitation energies
  Eigen::VectorXd _eigenvalues;
  // The occupied NTO matrix
  std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _occNTOs;
  // The virtual NTO matrix
  std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _virtNTOs;
  // The occupied NTO matrix eigenvalues
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _occEigenvalues;
  // The virtual NTO matrix eigenvalues
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _virtEigenvalues;
  // The densities related to the excitation, vector length is the number of subsystems
  std::vector<MatrixInBasis<SCFMode>> _transitionDensity, _holeDensity, _particleDensity, _holeDensityCorrection;
  // Calculates the NTO matrices, eigenvalues, eigenvectors
  void calcNTOs(unsigned int iState);
  // Calculates the transition density
  void calcTransitionDensity(unsigned int iState);
  // Calculates the hole and particle density and the hole density correction
  void calcParticleAndHoleDensity(unsigned int iState);
  // Makes the NTOs and transition, hole and particle densities available for a specific state
  void makeAvailable(unsigned iState) {
    if (_state != iState) {
      calcNTOs(iState);
      calcTransitionDensity(iState);
      calcParticleAndHoleDensity(iState);
      _state = iState;
    }
  }
  // The state for which the calculation has been performed
  unsigned int _state = 999999;
};

} /* namespace Serenity */

#endif /* NTOCALCULATOR_H_ */