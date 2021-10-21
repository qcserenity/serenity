/**
 * @file DOIBasedSelector.h
 *
 * @date Apr 3, 2019
 * @author Moritz Bensberg
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
#ifndef ANALYSIS_PAOSELECTION_DOIBASEDSELECTOR_H_
#define ANALYSIS_PAOSELECTION_DOIBASEDSELECTOR_H_

/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/PAOSelector.h" //Base class.
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dense matrices.
#include <Eigen/SparseCore> //Sparse matrices.

namespace Serenity {

/* Forward Declarations */
class PAOController;
class BasisFunctionOnGridController;
class AtomCenteredBasisController;

/**
 * @class DOIBasedSelector DOIBasedSelector.h
 * @brief PAO "selecter" based on the differential overlap between PAOs and occupied orbitals.
 *        As proposed by Neese and coworkers:\n
 *        J. Chem. Phys. 143, 034108 (2015).
 */
class DOIBasedSelector : public PAOSelector {
 public:
  /**
   * @brief Constructor. Saves attributes only.
   * @param occupiedCoefficients The coefficients of the occupied orbitals.
   * @param paoController The PAOController.
   * @param basOnGridController The basis function on grid controller.
   * @param atomCenteredBasisController The underlying atom centered basis controller.
   * @param doiThreshold The DOI selection threshold.
   * @param mnpPreThreshold Mulliken net population threshold for map-based prescreening.
   */
  DOIBasedSelector(std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients, std::shared_ptr<PAOController> paoController,
                   std::shared_ptr<BasisFunctionOnGridController> basOnGridController,
                   std::shared_ptr<AtomCenteredBasisController> atomCenteredBasisController,
                   Eigen::VectorXd doiThresholds, double mnpPreThreshold)
    : _occupiedCoefficients(occupiedCoefficients),
      _paoController(paoController),
      _basOnGridController(basOnGridController),
      _atomCenteredBasisController(atomCenteredBasisController),
      _doiThresholds(doiThresholds),
      _mnpPreThreshold(mnpPreThreshold){};
  /**
   * @brief Default destructor.
   */
  virtual ~DOIBasedSelector() = default;

  /**
   * @brief Selects PAOs for each orbital.
   * @return The orbital selection.
   */
  std::shared_ptr<Eigen::SparseMatrix<int>> selectPAOs() override final;

 private:
  ///@brief The coefficients of the occupied orbitals.
  std::shared_ptr<Eigen::MatrixXd> _occupiedCoefficients;
  ///@brief The PAOController.
  std::shared_ptr<PAOController> _paoController;
  ///@brief The basis function on grid controller.
  std::shared_ptr<BasisFunctionOnGridController> _basOnGridController;
  ///@brief The underlying atom centered basis controller.
  std::shared_ptr<AtomCenteredBasisController> _atomCenteredBasisController;
  ///@brief The DOI selection threshold.
  Eigen::VectorXd _doiThresholds;
  ///@brief Mulliken net population threshold for map-based prescreening.
  double _mnpPreThreshold;
};

} /* namespace Serenity */

#endif /* ANALYSIS_PAOSELECTION_DOIBASEDSELECTOR_H_ */
