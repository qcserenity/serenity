/**
 * @file VirtualOrbitalSelectionAlgorithms.h
 *
 * @date May 28, 2020
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

#ifndef MISC_VIRTUALORBITALSELECTIONALGORITHMS_H_
#define MISC_VIRTUALORBITALSELECTIONALGORITHMS_H_

/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "settings/EmbeddingSettings.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declaration */
class SystemController;
/**
 * @class  VirtualOrbitalSelectionAlgorithms VirtualOrbitalSelectionAlgorithms.h
 * @brief  Performs different virtual orbital selection algorithms
 */
template<Options::SCF_MODES SCFMode>
class VirtualOrbitalSelectionAlgorithms {
 public:
  /**
   * @brief Default destructor.
   */
  virtual ~VirtualOrbitalSelectionAlgorithms() = default;

  /**
   * @brief Constructor.
   * @param activeSystem The system controller for the active system.
   * @param environmentSystems The systemcontroller of the environment systems
   */
  VirtualOrbitalSelectionAlgorithms(std::shared_ptr<SystemController> activeSystem,
                                    std::vector<std::shared_ptr<SystemController>> environmentSystems);

  /**
   * @brief Excludes artificially shifted occupied environmental orbitals in
   *        the virtual orbital space of the subsystem in case of projection-based embedding
   * @param coefs Coefficient matrix of the active subsystem
   * @param nOcc Number of occupied orbitals
   * @param indices The new indices of the partitioned orbitals
   */
  void excludeProjection(CoefficientMatrix<SCFMode>& coefs, SpinPolarizedData<SCFMode, unsigned int>& nOcc,
                         SpinPolarizedData<SCFMode, std::vector<unsigned int>>& indices);

  /**
   * @brief Partitiones the canonical virtual orbital space of a subsystem
   * @param coefs Coefficient matrix of the active subsystem
   * @param nOcc Number of occupied orbitals
   * @param nVirt Number of virtual orbitals
   * @param indices The new indices of the partitioned orbitals
   * @param localThresh Threshold for the selection of localized orbitals
   * @param envThres Threshold for the selection of environmental orbitals
   */
  void virtualCanonicalOrbitalSpaceSelection(CoefficientMatrix<SCFMode>& coefs, SpinPolarizedData<SCFMode, unsigned int>& nOcc,
                                             SpinPolarizedData<SCFMode, unsigned int>& nVirt,
                                             SpinPolarizedData<SCFMode, std::vector<unsigned int>>& indices,
                                             double& localThresh, double& envThres, bool onlyOne);
  /**
   * @brief Localizes the virtual orbital space of the active subsystem on itself or the environment
   * @param coefs Coefficient matrix of the active subsystem
   * @param eigenValues Eigenvalues of the active subsystem
   * @param nOcc Number of occupied orbitals
   * @param nVirt Number of virtual orbitals
   * @param embeddedFock The embedded Fock matrix of the active subsystem
   * @param local decides whether the the orbitals located on the active subsystem or on the environment should be
   * selected
   */
  void virtualOrbitalSpaceLocalization(CoefficientMatrix<SCFMode>& coefs,
                                       SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenValues,
                                       SpinPolarizedData<SCFMode, unsigned int>& nOcc,
                                       SpinPolarizedData<SCFMode, unsigned int>& nVirt,
                                       FockMatrix<SCFMode>& embeddedFock, bool local);

  /**
   * @brief Mixes the occupied states of env[0] with the virtuals of env[1]
   * @param relaxation Whether additional relaxation is performed
   * @param settings Embedding settings
   */
  void occupiedVirtualMixing(bool relaxation, EmbeddingSettings settings);

  /**
   * @brief calculates embedding Fock matrix
   * @param activeSystem system for which the Fock matrix is calculated
   * @param environmentSystem environment systems
   * @param settings Embedding settings
   * @return embedded Fock matrix
   */
  FockMatrix<SCFMode> calcEmbeddedFockMatrix(std::shared_ptr<SystemController>& activeSystem,
                                             std::vector<std::shared_ptr<SystemController>>& environmentSystem,
                                             EmbeddingSettings settings);

 private:
  /// @brief The active system
  std::shared_ptr<SystemController> _activeSystem;
  /// @brief The environment systems
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
};

} /* namespace Serenity */

#endif /* MISC_VIRTUALORBITALSELECTIONALGORITHMS_H_ */