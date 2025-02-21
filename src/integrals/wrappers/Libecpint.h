/**
 * @file   Libecpint.h
 *
 * @date   Dec 1, 2017
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
#ifndef LIBECPINT_H_
#define LIBECPINT_H_

/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h" // RESTRICTED/UNRESTRICTED
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace libecpint {
struct GaussianShell;
}

namespace Serenity {

// Forward declarations
class Atom;
class BasisController;
class Shell;
template<Options::SCF_MODES SCFMode>
class MatrixInBasis;
template<Options::SCF_MODES SCFMode>
class SPMatrix;
template<Options::SCF_MODES SCFMode>
class MatrixinBasis;
template<Options::SCF_MODES SCFMode>
using DensityMatrix = MatrixInBasis<SCFMode>;
class AtomCenteredBasisController;

/**
 * @class Libecpint Libecpint.h
 *
 * @brief Wrapper for the library with the same name for computing effective core potentials.
 */
class Libecpint {
 private:
  // Purely static, never instantiated
  Libecpint();

 public:
  virtual ~Libecpint();
  /**
   * @brief Method for the calculation of effective core potential integrals.
   *
   * Calculates sum_v <i|v|j>, where i and j are basis functions of basisController and v
   * are the gaussian functions representing the core potential. These functions are taken
   * from the atoms.
   *
   * @param   basisController The controller for the basis for which the integrals shall be
   *                          computed
   * @param   atoms The atoms from which the ECP functions are taken. Not every atom needs to
   *                actually have ECP functions.
   * @returns A matrix containing the effective core potential in basis form so that it can be
   *          added to, e.g., the core Hamiltonian.
   */
  static MatrixInBasis<Options::SCF_MODES::RESTRICTED>
  computeECPIntegrals(std::shared_ptr<BasisController> basisController, const std::vector<std::shared_ptr<Atom>>& atoms);
  /**
   * @brief Method for the calculation of effective core potential integrals.
   *
   * Calculates sum_v <i|v|j>, where i and j are basis functions of basisController A/B and v
   * are the gaussian functions representing the core potential. These functions are taken
   * from the atoms.
   *
   * @param   basisControllerA The controller for the basis on the LHS of the integrals that shall be
   *                           computed.
   * @param   basisControllerB The controller for the basis on the RHS of the integrals that shall be
   *                           computed.
   * @param   atoms The atoms from which the ECP functions are taken. Not every atom needs to
   *                actually have ECP functions.
   * @returns A matrix containing the effective core potential in basis form so that it can be
   *          added to, e.g., the core Hamiltonian.
   */
  static SPMatrix<Options::SCF_MODES::RESTRICTED> computeECPIntegrals(std::shared_ptr<BasisController> basisControllerA,
                                                                      std::shared_ptr<BasisController> basisControllerB,
                                                                      const std::vector<std::shared_ptr<Atom>>& atoms);
  /**
   * @brief Calculates the nuclear gradient contribution resulting from the ECPs.
   * @param basisController  The basis of to be used
   * @param atoms            All atoms present in the molecule.
   * @param density          The current electron density, basis must match the basisController.
   * @return Eigen::MatrixXd The gradient contributions.
   */
  static Eigen::MatrixXd computeECPGradientContribution(std::shared_ptr<AtomCenteredBasisController> basisController,
                                                        const std::vector<std::shared_ptr<Atom>>& atoms,
                                                        const DensityMatrix<RESTRICTED>& density);

 private:
  /**
   * @brief Constructs a libecpint style gaussian shell.
   * @param atoms All atoms in the molecule.
   * @param shell The Serenity::Shell
   * @return libecpint::GaussianShell The gaussian shell.
   */
  static libecpint::GaussianShell makeShell(const std::vector<std::shared_ptr<Atom>>& atoms,
                                            std::shared_ptr<const Serenity::Shell> shell);
};

} /* namespace Serenity */

#endif /* LIBECPINT_H_ */
