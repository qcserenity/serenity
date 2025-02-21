/**
 * @file FDEPotentials.h
 *
 * @date Nov 24, 2016
 * @author Jan Unsleber
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

#ifndef POTENTIALS_BUNDLES_FDEPOTENTIALS_H_
#define POTENTIALS_BUNDLES_FDEPOTENTIALS_H_

/* Include Serenity Internal Headers */
#include "potentials/Potential.h"
#include "potentials/bundles/PotentialBundle.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class NAddFuncPotential;

/**
 * @class FDEPotentials FDEPotentials.h
 * @brief A class containing all the potentials relevant for an FDE-SCF.
 *
 * The potential bundle implies one active subsystem, embedded into
 * (possibly) multiple other subsystems.
 * This class is restricted to FDE calculations employing non-additive
 * functionals for both kinetic and exchange-correlation interaction energies.
 * For other frozen density type embedding schemes using projectors or optimized
 * effective potentials, see the other bundles.
 */
template<Options::SCF_MODES SCFMode>
class FDEPotentials : public PotentialBundle<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param activeSystemPot The active system (isolated) potential.
   *                         This is a general PotentialBundle to allow for any
   *                         method for the active system (HF, DFT, ...)
   * @param esiPot          The electrostatic interaction potential.
   *                         This is a general PotentialBundle in order to allow
   *                         for future implementations of i.e. potentials read
   *                         from disk.
   * @param naddXC          The non-additive exchange-correlation potential.
   *                         This potential is fixed to be a NAddFuncPotential
   *                         as this is part of what defines this particular
   *                         bundle.
   * @param naddKin         The non-additive kinetic energy potential.
   *                         This potentials is fixed to be a NAddFuncPotential
   *                         as this is part of what defines this particular
   *                         bundle.
   * @param ecp             The interaction between ECPs of different subsystems.
   * @param pcm             The interaction of the active subsystem with the implicit
   *                         solvent model
   */
  FDEPotentials(std::shared_ptr<PotentialBundle<SCFMode>> activeSystemPot, std::shared_ptr<PotentialBundle<SCFMode>> esiPot,
                std::vector<std::shared_ptr<NAddFuncPotential<SCFMode>>> naddXC,
                std::vector<std::shared_ptr<Potential<SCFMode>>> naddKin, std::shared_ptr<Potential<SCFMode>> ecp,
                std::shared_ptr<Potential<SCFMode>> pcm);

  /// @brief Default destructor.
  virtual ~FDEPotentials() = default;

  /**
   * @brief A function to get the entire Fock matrix.
   * @param P The density matrix.
   * @param energies The controller to add all the energy contributions to.
   * @return The current Fock matrix (rebuilt on every call).
   *         (Note that this does not imply that each Potential is rebuilt, they
   *          may very well be cached. But all potentials contributing are added
   *          together in every call)
   */
  FockMatrix<SCFMode> getFockMatrix(const DensityMatrix<SCFMode>& P,
                                    std::shared_ptr<EnergyComponentController> energies) override final;

  /**
   * @brief Returns added geometric gradients (cartesian) of all underlying potentials.
   *
   * The gradients generated here are in general not free of rotation and/or
   * translation.
   *
   * @return Returns the geometric gradients.
   */
  Eigen::MatrixXd getGradients() override final;

  std::shared_ptr<PotentialBundle<SCFMode>> getActiveSystemPotentials();
  std::shared_ptr<PotentialBundle<SCFMode>> getESIPotentials();
  std::vector<std::shared_ptr<NAddFuncPotential<SCFMode>>> getNaddXCPotentials();
  std::vector<std::shared_ptr<Potential<SCFMode>>> getNaddKinPotentials();
  std::shared_ptr<Potential<SCFMode>> getECPInteractionPotential();
  std::shared_ptr<Potential<SCFMode>> getPCMPotential();

  /**
   * @brief Getter for the gradients of external point charges.
   * @return The point charge gradients.
   */
  Eigen::MatrixXd getPointChargeGradients() override;

 private:
  ///@brief The active system potentials
  std::shared_ptr<PotentialBundle<SCFMode>> _activeSystemPot;
  ///@brief The clasical interaction terms
  std::shared_ptr<PotentialBundle<SCFMode>> _esiPot;
  ///@brief The non additive exchange correlation interaction.
  std::vector<std::shared_ptr<NAddFuncPotential<SCFMode>>> _naddXC;
  ///@brief The non additive kinetic interaction.
  std::vector<std::shared_ptr<Potential<SCFMode>>> _naddKin;
  ///@brief The contributions of the effective core potentials located in the environment.
  std::shared_ptr<Potential<SCFMode>> _ecp;
  ///@brief The contributions of the implicit solvation model.
  std::shared_ptr<Potential<SCFMode>> _pcm;
};

} /* namespace Serenity */

#endif /* POTENTIALS_BUNDLES_FDEPOTENTIALS_H_ */
