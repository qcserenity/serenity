/**
 * @file ECPInteractionPotential.h
 *
 * @date Jun 11, 2018
 * @author Moritz Bensberg
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

#ifndef POTENTIALS_ECPINTERACTIONPOTENTIAL_H_
#define POTENTIALS_ECPINTERACTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "potentials/EffectiveCorePotential.h"
#include "potentials/Potential.h"

namespace Serenity {
/* Forward Declarations */
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
class Atom;

/**
 * @class ECPInteractionPotential ECPInteractionPotential.h
 * @brief Calculates the ECP interaction between subsystems.
 */
template<Options::SCF_MODES SCFMode>
class ECPInteractionPotential : public Potential<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param actSystem The system controller for the active system.
   * @param actAtoms The atoms of the active system.
   * @param envAtoms The atoms of the environment.
   * @param envDensities The densities of the environment.
   * @param basis The basis the Fock matrix should be expressed in.
   */
  ECPInteractionPotential(std::shared_ptr<SystemController> actSystem, std::vector<std::shared_ptr<Atom>> actAtoms,
                          std::vector<std::shared_ptr<Atom>> envAtoms,
                          std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensities,
                          std::shared_ptr<BasisController> basis);

  /**
   * @brief Default destructor.
   */
  virtual ~ECPInteractionPotential() = default;
  /**
   * @brief Getter for the actual potential.
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override final;
  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;
  /**
   * @brief Geometry gradient contribution from this Potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override final;

 private:
  std::weak_ptr<SystemController> _actSystem;
  ///@brief The active atoms.
  std::vector<std::shared_ptr<Atom>> _actAtoms;
  ///@brief The environment atoms.
  std::vector<std::shared_ptr<Atom>> _envAtoms;
  ///@brief The environment densities.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _envDensities;
  double _actEnvDensityECPEnergy;
  ///@brief The ECP of the environment atoms acting on the active density.
  std::shared_ptr<EffectiveCorePotential<SCFMode>> _envActDensECP;
  ///@brief The ECP of the active atoms acting on the environment density.
  std::vector<std::shared_ptr<EffectiveCorePotential<SCFMode>>> _actEnvDensECPs;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ECPINTERACTIONPOTENTIAL_H_ */
