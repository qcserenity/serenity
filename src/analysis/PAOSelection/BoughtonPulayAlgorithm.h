/**
 * @file BoughtonPulayAlgorithm.h
 *
 * @date Dec 10, 2018
 * @author Moritz Bensberg
 *
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

#ifndef ANALYSIS_PAOSELECTION_BOUGHTONPULAYALGORITHM_H_
#define ANALYSIS_PAOSELECTION_BOUGHTONPULAYALGORITHM_H_

/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/PAOSelector.h" //Base class.
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {

/* Forward Declarations */
class AtomCenteredBasisController;
class Geometry;
class OneElectronIntegralController;

/**
 * @class BoughtonPulayAlgorithm BoughtonPulayAlgorithm.h
 * @brief Implements the Boughton--Pulay algorithm.
 *
 * For more details see:\n
 *   J. Comput. Chem. 143, 024105 (1993)\n
 */
class BoughtonPulayAlgorithm : public PAOSelector {
 public:
  /**
   * @brief Constructor.
   * @param oneElectronIntegralController The one electron integral controller.
   * @param atomCenteredBasisController The atom centered basis controller.
   * @param mullikenGrossCharges The Mulliken charges.
   * @param occupiedCoefficients The coefficients of the occupied orbitals.
   * @param completeness The completeness threshold.
   */
  BoughtonPulayAlgorithm(std::shared_ptr<OneElectronIntegralController> oneElectronIntegralController,
                         std::shared_ptr<AtomCenteredBasisController> atomCenteredBasisController,
                         std::shared_ptr<Eigen::MatrixXd> mullikenGrossCharges,
                         std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients, double completeness = 0.02);
  /**
   * @brief Destructor.
   */
  virtual ~BoughtonPulayAlgorithm() = default;
  /**
   * @brief Assigns atoms to an orbital.
   * @param mullikenGrossCharges The atom-wise mulliken gross charges of the orbital.
   * @param c_i The coefficients of the orbital.
   * @return The atom indices.
   */
  std::vector<unsigned int> assignAtoms(Eigen::VectorXd mullikenGrossCharges, const Eigen::VectorXd& c_i);
  /**
   * @brief Selects PAOs for each orbital.
   * @return The orbital selection.
   */
  std::shared_ptr<Eigen::SparseMatrix<int>> selectPAOs() override final;

 private:
  // The one electron integral controller for the overlap integrals.
  std::shared_ptr<OneElectronIntegralController> _oneElectronIntegralController;
  // The atom centered basis controller.
  std::shared_ptr<AtomCenteredBasisController> _atomCenteredBasisController;
  // The orbital-wise Mulliken charges.
  std::shared_ptr<Eigen::MatrixXd> _mullikenGrossCharges;
  // The coefficients of the occupied orbitals.
  std::shared_ptr<Eigen::MatrixXd> _occupiedCoefficients;
  // The completeness threshold.
  double _completeness;
};

} /* namespace Serenity */

#endif /* ANALYSIS_PAOSELECTION_BOUGHTONPULAYALGORITHM_H_ */
