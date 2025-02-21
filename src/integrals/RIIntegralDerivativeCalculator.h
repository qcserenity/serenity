/**
 * @file   RIIntegralDerivativeCalculator.h
 *
 * @date   Jan 10, 2025
 * @author Anton Rikus
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef RIINTEGRALDERIVATIVECALCULATOR_H
#define RIINTEGRALDERIVATIVECALCULATOR_H
/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declarations */
class AtomCenteredBasisController;
class Geometry;
template<Options::SCF_MODES SCFMode>
class MatrixInBasis;

/**
 * @class RIIntegralDerivativeCalculator RIIntegralDerivativeCalculator.h
 *
 * This class calculates derivatives of Coulomb integrals using the resolution of the identity approximation and
 * contracts them with arbitrary density matrices.
 *
 * References:\n [1] D. Rappoport, F. Furche, J. Chem. Phys. 122, 064105 (2005)
 *
 */
class RIIntegralDerivativeCalculator {
 public:
  /**
   * @brief Constructor
   * @param basisController    The reference basis.
   * @param auxBasisController The auxiliary basis.
   */
  RIIntegralDerivativeCalculator(std::shared_ptr<AtomCenteredBasisController> basisController,
                                 std::shared_ptr<AtomCenteredBasisController> auxBasisController);

  /**
   * @brief This function contracts given densities density1, density2a and density2b with the derivatives of Coulomb
   * integrals using the resolution of the identity. It builds the 2-particle density matrix as 2 * density1 * density1
   * + density2a * density2b + density2b * density2a. The derivative of the metric is included as well (see equation
   * [1].19a and [1].19b). It only takes restricted MatrixInBasis as arguments, because the two-electron integrals in AO
   * basis are spin-independent. With an unrestricted density matrix, simply add both spin components using .total()
   * when calling the function.
   * @param density1 The density-like matrix to be contracted.
   * @param density2a The density-like matrix to be contracted.
   * @param density2b The density-like matrix to be contracted.
   * @returns The gradient contribution as a matrix of shape (number of atoms, 3).
   */
  const Eigen::MatrixXd getCoulombIntegralDerivatives(const MatrixInBasis<RESTRICTED>& density1,
                                                      const MatrixInBasis<RESTRICTED>& density2a,
                                                      const MatrixInBasis<RESTRICTED>& density2b);

 private:
  std::shared_ptr<AtomCenteredBasisController> _basisController;
  std::shared_ptr<AtomCenteredBasisController> _auxBasisController;
};

} /* namespace Serenity */
#endif /* RIINTEGRALDERIVATIVECALCULATOR_H */