/**
 * @file EdmistonRuedenbergLocalization.h
 *
 * @date Nov 3, 2016
 * @author David Schnieders
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

#ifndef POSTSCF_LOCALIZATION_EDMISTONRUEDENBERGLOCALIZATION_H_
#define POSTSCF_LOCALIZATION_EDMISTONRUEDENBERGLOCALIZATION_H_

/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/Localization.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

class SystemController;
template<Options::SCF_MODES SCFMode>
class OrbitalController;

template<Options::SCF_MODES SCFMode>
class EdmistonRuedenbergLocalization : public Localization<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param The system to of which the orbitals are to be localized.
   */
  EdmistonRuedenbergLocalization(std::shared_ptr<SystemController> systemController);
  /**
   * @brief Destructor
   */
  virtual ~EdmistonRuedenbergLocalization() = default;

  /**
   * @brief Performs Jacobi rotations to localize the molecular orbitals as described in
   *        Edmiston, C.; Ruedenberg, K. Rev. Mod. Phys. 35, 457 (1963)
   * @param orbitals The orbitals to be localized (They are changed in place!)
   * @param maxSweeps The maximum number of iterations, where one iteration performs a
   *        Jacobi rotation on every orbital.
   * @param orbitalRange  The range of orbitals to be rotated in.
   */
  virtual void localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) override final;

 private:
  std::shared_ptr<SystemController> _systemController;
};

} /* namespace Serenity */

#endif /* POSTSCF_LOCALIZATION_EDMISTONRUEDENBERGLOCALIZATION_H_ */
