/**
 * @file FosterBoysLocalization.h
 *
 * @date Nov 5, 2015
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

#ifndef LOCALIZATION_FOSTERBOYSLOCALIZATION_H_
#define LOCALIZATION_FOSTERBOYSLOCALIZATION_H_

/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/Localization.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class OrbitalController;
class SystemController;

/**
 * @class FosterBoysLocalization FosterBoysLocalization.h
 *
 * @brief Localizes the Orbitals using the Foster-Boys algorithm.
 *
 * For more information about the two algorithms see Ref.:
 * J. Chem. Phys. 61, 3905 (1974)
 *
 */
template<Options::SCF_MODES SCFMode>
class FosterBoysLocalization : public Localization<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param systemController The system of which the orbitals are to be localized.
   */
  FosterBoysLocalization(std::shared_ptr<SystemController> systemController);
  /**
   * @brief Destructor
   */
  virtual ~FosterBoysLocalization() = default;

  /**
   * @brief @see Localization::localizeOrbitals()
   * @param orbitals The orbitals to be localized (They are changed in place!)
   * @param maxSweeps The maximum number of iterations, where one iteration performs a
   *        Jacobi rotation on every orbital.
   * @param orbitalRange  The range of orbitals to be rotated in.
   */
  void localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                        SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) override final;

 private:
  const std::shared_ptr<SystemController> _system;
};

} /* namespace Serenity */

#endif /* LOCALIZATION_FOSTERBOYSLOCALIZATION_H_ */
