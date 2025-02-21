/**
 * @file NonOrthogonalLocalization.h
 *
 * @date Dec 7, 2015
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

#ifndef POSTSCF_LOCALIZATION_NONORTHOGONALLOCALIZATION_H_
#define POSTSCF_LOCALIZATION_NONORTHOGONALLOCALIZATION_H_

/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/Localization.h" //Localization
#include "settings/Options.h"                          // Options::SCF_MODES
/* Include Std and External Headers */
#include <Eigen/Dense> // Eigen::MatrixXd
#include <memory>      // shared_ptr

namespace Serenity {

/* Forward declarations */
class SystemController;

template<Options::SCF_MODES SCFMode>
class NonOrthogonalLocalization : public Localization<SCFMode> {
 public:
  NonOrthogonalLocalization(std::shared_ptr<SystemController> systemController);
  virtual ~NonOrthogonalLocalization() = default;

  /**
   * @brief @see Localization::localizeOrbitals()
   * @param orbitals The orbitals to be localized (They are changed in place!)
   * @param maxSweeps The maximum number of iterations, where one iteration performs a
   *        Jacobi rotation on every orbital.
   * @param orbitalRange  The range of orbitals to be rotated in.
   */
  virtual void localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) override final;

 private:
  /**
   * @brief Localizes the molecular orbitals as described in
   *        H. Feng, J. Bian, L. Li, W. Yang: An efficient method for constructing
   *        nonorthogonal localized molecular orbitals; J. Chem. Phys. 120, 9458 (2004).
   *        The resulting orbitals are not orthogonal anymore!
   * @param coefficients
   * @param nOccOrbs
   * @param maxSweeps
   */
  void doLocalization(Eigen::MatrixXd& coefficients, unsigned int nOccOrbs, unsigned int maxSweeps);

  const std::shared_ptr<SystemController> _system;
};

} /* namespace Serenity */

#endif /* POSTSCF_LOCALIZATION_NONORTHOGONALLOCALIZATION_H_ */
