/**
 * @file   OneElectronIntegralDerivativeCalculator.h
 *
 * @date   Jun 17, 2024
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
#ifndef ONEELECTRONINTEGRALDERIVATIVECALCULATOR_H
#define ONEELECTRONINTEGRALDERIVATIVECALCULATOR_H
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
class Point;
class ExternalChargeController;

/**
 * @class OneElectronIntegralDerivativeCalculator OneElectronIntegralDerivativeCalculator.h
 *
 * This class calculates derivatives of one-electron integrals and contracts them with arbitrary density matrices.
 *
 */
class OneElectronIntegralDerivativeCalculator {
 public:
  /**
   * @brief Constructor
   * @param basisController   The reference basis.
   * @param geometry          The reference geometry, needed for the nuclei-electron attraction integrals.
   * @param pointCharges      A list of point charges (including the effective atomic charges) with their charge and
   * position (optional).
   */
  OneElectronIntegralDerivativeCalculator(std::shared_ptr<AtomCenteredBasisController> basisController,
                                          std::shared_ptr<const Geometry> geometry,
                                          const std::vector<std::pair<double, Point>> pointCharges = {});

  /**
   * @brief Constructor
   * @param basisController            The reference basis.
   * @param geometry                   The reference geometry, needed for the nuclei-electron attraction integrals.
   * @param externalChargeController   The controller for the external charges.
   */
  OneElectronIntegralDerivativeCalculator(std::shared_ptr<AtomCenteredBasisController> basisController,
                                          std::shared_ptr<const Geometry> geometry,
                                          std::shared_ptr<ExternalChargeController> externalChargeController);

  virtual ~OneElectronIntegralDerivativeCalculator() = default;

  /**
   * @brief This function contracts a given density with the derivatives of the kinetic and nuclear integrals.
   *        It only takes a restricted MatrixInBasis as argument, because the one-electron integrals in AO basis are
   * spin-independent. With an unrestricted density matrix, simply add both spin components using .total() when calling
   * the function.
   * In case pointCharges were provided in the constructor, the derivatives have the shape (total number of charges, 3).
   * @param density The density-like matrix to be contracted.
   * @returns The gradient contribution as a matrix of shape (number of atoms, 3).
   */
  const Eigen::MatrixXd getNucKinDerivative(const MatrixInBasis<RESTRICTED>& density);

  /**
   * @brief This function contracts a given density with the derivatives of the overlap integrals.
   *        It only takes a restricted MatrixInBasis as argument, because the overlap integrals in AO basis are
   * spin-independent. With an unrestricted density matrix, simply add both spin components using .total() when calling
   * the function.
   * @param density The density-like matrix to be contracted.
   * @returns The gradient contribution as a matrix of shape (number of atoms, 3).
   */
  const Eigen::MatrixXd getOverlapDerivative(const MatrixInBasis<RESTRICTED>& density);

 private:
  static std::vector<std::pair<double, Point>> getAllCharges(std::shared_ptr<const Geometry> geometry,
                                                             std::shared_ptr<ExternalChargeController> externalChargeController);
  std::shared_ptr<AtomCenteredBasisController> _basisController;
  std::shared_ptr<const Geometry> _geometry;
  const std::vector<std::pair<double, Point>> _pointCharges;
};

} /* namespace Serenity */
#endif /* ONEELECTRONINTEGRADERIVATIVECALCULATOR_H */
