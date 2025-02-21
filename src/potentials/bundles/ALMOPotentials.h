/**
 * @file ALMOPotentials.h
 *
 * @author Lukas Lampe
 * @date Nov 6, 2024
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

#ifndef POTENTIALS_ALMOPOTENTIALS_H_
#define POTENTIALS_ALMOPOTENTIALS_H_
/* Include Serenity Internal Headers */
#include "data/matrices/FockMatrix.h"           //Fock matrix definition.
#include "potentials/Potential.h"               //Potential definition.
#include "potentials/bundles/PotentialBundle.h" //Potential bundle definition.
/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {
/* Forward declarations */
class SystemController;
class EnergyComponentController;

/**
 * @class ALMOPotentials ALMOPotentials.h
 * @brief A class for the calculation of the potential for absolutely localized molecular orbitals.\n\n
 *
 * Absolutely localized molecular orbitals:\n
 * E.g., according to:\n
 *  J. Chem. Theory Comput. 5, 2702-2716 (2009) and \n
 *  J. Chem. Phys. 138, 134119 (2013).\n
 */
template<Options::SCF_MODES SCFMode>
class ALMOPotentials : public PotentialBundle<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param activeSystem         The active system.
   * @param environmentSystems   The environment systems.
   */
  ALMOPotentials(std::shared_ptr<SystemController> activeSystem,
                 std::vector<std::shared_ptr<SystemController>> environmentSystems);
  /**
   * @brief Defualt destrutor.
   */
  virtual ~ALMOPotentials() = default;

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

 private:
  /// @brief The active system controller.
  std::shared_ptr<SystemController> _activeSystem;
  /// @brief The environment system controllers.
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  /// @brief The supersystem controller.
  std::shared_ptr<SystemController> _supersystem;
  /// @brief The overlap matrix of the active system basis set with all environment basis sets.
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _s_ABs;
  /// @brief The potential in matrix representation.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  /// @brief The potential bundle of the supersystem.
  std::shared_ptr<PotentialBundle<SCFMode>> _superPotentialBundle;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ALMOPOTENTIALS_H_ */
