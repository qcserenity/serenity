/**
 * @file PBEPotentials.h
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

#ifndef POTENTIALS_BUNDLES_PBEPOTENTIALS_H_
#define POTENTIALS_BUNDLES_PBEPOTENTIALS_H_

/* Include Serenity Internal Headers */
#include "potentials/Potential.h"
#include "potentials/bundles/PotentialBundle.h"

namespace Serenity {

/**
 * @class PBEPotentials PBEPotentials.h
 * @brief A class containing all the potentials relevant for an SCF
 *        with projection based embedding potentials.
 */
template<Options::SCF_MODES SCFMode>
class PBEPotentials : public PotentialBundle<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param activeSystemPot Active fock matrix.
   * @param esiPot Electrostatic fock matrix.
   * @param naddXC Nadd-XC fock matrix.
   * @param projection Projection/Nadd-kinetic fock matrix.
   * @param ecp ECP interaction fock matrix.
   * @param pcm Implicit solvent model.
   */
  PBEPotentials(std::shared_ptr<PotentialBundle<SCFMode>> activeSystemPot, std::shared_ptr<PotentialBundle<SCFMode>> esiPot,
                std::shared_ptr<Potential<SCFMode>> naddXC, std::shared_ptr<Potential<SCFMode>> projection,
                std::shared_ptr<Potential<SCFMode>> ecp, std::shared_ptr<Potential<SCFMode>> pcm);
  /// @brief Default destructor.
  virtual ~PBEPotentials() = default;

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
   * @brief Returns gradients of all underlying potentials
   */
  Eigen::MatrixXd getGradients() override final;

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
  std::shared_ptr<Potential<SCFMode>> _naddXC;
  ///@brief The non additive kinetic interaction.
  std::shared_ptr<Potential<SCFMode>> _projection;
  ///@brief The ECP potential from the environment.
  std::shared_ptr<Potential<SCFMode>> _ecp;
  ///@brief The PCM potential from the implicit solvent model.
  std::shared_ptr<Potential<SCFMode>> _pcm;
};

} /* namespace Serenity */

#endif /* POTENTIALS_BUNDLES_PBEPOTENTIALS_H_ */
