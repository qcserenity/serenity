/**
 * @file   Localization.h
 *
 * @date   Apr 22, 2014
 * @author Thomas Dresselhaus
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
#ifndef LOCALIZATION_H_
#define LOCALIZATION_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <vector>

namespace Serenity {
/* Forward declarations */
namespace Options {
enum class SCF_MODES;
}
template<Options::SCF_MODES SCFMode>
class OrbitalController;

/**
 * @class  Localization Localization.h
 * @brief  Interface for orbital localization routines
 */
template<Options::SCF_MODES SCFMode>
class Localization {
 public:
  Localization() = default;
  virtual ~Localization() = default;
  /**
   * @brief Localizes a set of orbitals.
   *
   * This method must be overridden in actual implementations. Localization of the orbitals should
   * not affect properties of the system, since they should be invariant under rotations among
   * the occupied orbitals only (or among the virtual orbitals only). This at least holds in
   * Hartree-Fock theory.
   *
   * @param orbitals The orbitals to localize. Will be mutated of course.
   * @param maxSweeps Maximum number of localization iterations.
   * @param orbitalRange Indices of the orbitals to be localized.
   */
  virtual void localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) = 0;
};

} /* namespace Serenity */

#endif /* LOCALIZATION_H_ */
