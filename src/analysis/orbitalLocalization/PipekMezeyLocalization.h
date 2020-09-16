/**
 * @file   PipekMezeyLocalization.h
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
#ifndef PIPEKMEZEYLOCALIZATION_H_
#define PIPEKMEZEYLOCALIZATION_H_
/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/Localization.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <mutex>

namespace Serenity {
/* Forward declarations */
class SystemController;
/**
 * @class  Serenity::PipekMezeyLocalization PipekMezeyLocalization.h
 * @brief  Localizes orbitals using a Mulliken like measurement for the localization.
 *
 * Should scale with O(n^3)
 */
template<Options::SCF_MODES SCFMode>
class PipekMezeyLocalization : public Localization<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param systemController from which the information about geometry is taken
   */
  PipekMezeyLocalization(std::shared_ptr<SystemController> systemController);
  /**
   * @brief Destructor
   */
  virtual ~PipekMezeyLocalization() = default;

  /**
   * @brief Localizes a set of occupied orbitals by maximizing a Pipek-Mezey localization
   *        criterion using Jacobi rotations.
   * @param orbitals The OrbitalController containing the orbitals to be localized.
   * @param maxSweeps The maximum number of iterations, where one iteration performs a
   *        Jacobi rotation on every orbital.
   * @param orbitalRange  The range of orbitals to be rotated in.
   *
   * Reference: J. W. Boughton, P. Pulay; J. Comp. Chem. 14 (1993), 736 - 740.
   */
  void localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                        SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) override final;

 private:
  const std::shared_ptr<SystemController> _systemController;

  const double _convThreshold;
};

} /* namespace Serenity */

#endif /* PIPEKMEZEYLOCALIZATION_H_ */
