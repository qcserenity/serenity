/**
 * @file IBOLocalization.h
 *
 * @date Jun 15, 2016
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

#ifndef IBOLOCALIZATION_H_
#define IBOLOCALIZATION_H_

/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/Localization.h"
#include "data/matrices/CoefficientMatrix.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class OrbitalController;
class SystemController;

/**
 * @class IBOLocalization IBOLocalization.h
 * @brief Performs IBO and IAO Localizations
 *
 * Ref.:
 * Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical Concepts
 * Gerald Knizia, J. Chem. Theory Comput., 2013, 9 (11), pp 4834–4843
 */
template<Options::SCF_MODES SCFMode>
class IBOLocalization : public Localization<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param systemController The system to of which the orbitals are to be localized.
   * @param IAOsOnly Switch to stop after the generation of IAOs.
   * @param replaceVirtuals If true, the virtual valence orbitals are reconstructed to match the IAO orbital space.
   * @param enforceRestrictedOrbitals If true, restricted orbitals are enforced for alpha and beta. This only makes
   * sense if the original orbital set was already quasi restricted.
   */
  IBOLocalization(std::shared_ptr<SystemController> systemController, bool IAOsOnly = false,
                  bool replaceVirtuals = false, bool enforceRestrictedOrbitals = false);
  /**
   * @brief Destructor
   */
  virtual ~IBOLocalization() = default;

  /**
   * @brief @see Localization::localizeOrbitals()
   * @param orbitals The orbitals to be localized (They are changed in place!)
   * @param maxSweeps maximum number of iterations.
   * @param orbitalRange  The range of orbitals to be rotated in.
   */
  virtual void localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) override final;

 private:
  std::shared_ptr<SystemController> _system;
  const bool _IAOsOnly;
  const bool _replaceVirtuals;
  const bool _enforceRestrictedOrbitals;
  bool _virtualOrbitalsReplaced = false;

  void restrictOrbitals(CoefficientMatrix<SCFMode>& coefficientMatrix);
  SpinPolarizedData<SCFMode, unsigned int> restrictedOccupations();
};

} /* namespace Serenity */

#endif /* IBOLOCALIZATION_H_ */
