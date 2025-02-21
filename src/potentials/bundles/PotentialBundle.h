/**
 * @file PotentialBundle.h
 *
 * @date Nov 22, 2016
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

#ifndef POTENTIALS_BUNDLES_POTENTIALBUNDLE_H_
#define POTENTIALS_BUNDLES_POTENTIALBUNDLE_H_

/* Include Serenity Internal Headers */

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "energies/EnergyComponentController.h"
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/**
 * @class PotentialBundle PotentialBundle.h
 *
 * @brief An interface for all the different ways to update a Fock matrix.
 *
 * Via bundles of potentials.
 * The difference between a potential and a bundle of potentials is obviously
 * that a bundle of potentials will consist of several potentials and bundle
 * them for convenience (hence the name).
 * However a bundle also knows in which field specific energies have to be
 * registered inside of the Serenity::EnergyComponentController .
 * Thus the function to get the Fock matrix requires a density and an
 * Serenity::EnergyComponentController and automatically adds all energy
 * contributions associated with the potentials bundled.
 */
template<Options::SCF_MODES SCFMode>
class PotentialBundle {
 public:
  ///@brief Default constructor.
  PotentialBundle() = default;
  /// @brief Default destructor.
  virtual ~PotentialBundle() = default;

  /**
   * @brief A function to get the entire Fock matrix.
   * @param P The density matrix.
   * @param energies The controller to add all the energy contributions to.
   * @return The current Fock matrix (rebuilt on every call).
   *         (Note that this does not imply that each Potential is rebuilt, they
   *          may very well be cached. But all potentials contributing are added
   *          together in every call)
   */
  virtual FockMatrix<SCFMode> getFockMatrix(const DensityMatrix<SCFMode>& P,
                                            std::shared_ptr<EnergyComponentController> energies) = 0;

  /**
   * @brief Returns gradients of all underlying potentials
   */
  virtual Eigen::MatrixXd getGradients() = 0;

  /**
   * @brief Getter for the gradients of external point charges.
   * @return The point charge gradients.
   */
  virtual Eigen::MatrixXd getPointChargeGradients() {
    throw SerenityError("No geometry gradients available.");
    Eigen::MatrixXd gradientContr(1, 3);
    return gradientContr;
  }
};

} /* namespace Serenity */

#endif /* POTENTIALS_BUNDLES_POTENTIALBUNDLE_H_ */
