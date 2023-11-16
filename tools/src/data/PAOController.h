/**
 * @file PAOController.h
 *
 * @date Jan 27, 2019
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

#ifndef DATA_PAOCONTROLLER_H_
#define DATA_PAOCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h" //Coefficient matrix.
#include "data/matrices/DensityMatrix.h"     //Density matrix.
#include "misc/HelperFunctions.h"            //constructProjectionMatrixFromSparse
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dense matrices.
#include <Eigen/SparseCore> //Sparse matrices.
#include <memory>           //smrt_ptr

namespace Serenity {
/**
 * @class PAOController PAOController.h
 * @brief A controller that holds and constructs the projected atomic orbital (PAO)
 *        coefficients
 *
 *        PAO definition:\n
 *          \f$ \chi^\mathrm{PAO}_\mu = (1-\hat{P}^\mathrm{occ})\chi^\mathrm{AO}_\mu~~~\mathrm{,} \f$ \n
 *          where \f$ \hat{P}^\mathrm{occ}\f$, \f$ \chi^\mathrm{AO}\f$ are the projection operator on the
 *          space of the occupied orbitals and an AO basis function, respectively.\n
 *          In matrix representation the coefficient matrix is thus given by\n
 *          \f$ C^\mathrm{PAO} = 1-DS~~~\mathrm{,}\f$\n
 *          where \f$ D \f$ and \f$ S \f$ are the density and overlap matrix, respectively.\n\n
 *
 *          If occupied environment orbitals are present, their density matrix is projected
 *          into the basis of the AOs by computing a least square fit:\n
 *          \f$ S_{AA}^{-1}S_{AB} D_{BB} * S_{BA}S_{AA}^{-1} ~~~\mathrm{,}\f$ \n
 *          where \f$ S_{AA} \f$, \f$ S_{AB} \f$ and \f$ D_{BB} \f$ are the overlap matrix of
 *          the active system (A), the overlap matrix between active and environment (B) and the
 *          density matrix of the environment, respectively.
 *
 */
class PAOController {
 public:
  /**
   * @brief Constructor.
   * @param activeDensity The density matrix of the active system.
   * @param overlapMatrix The overlap matrix of the active system.
   * @param paoNormalizationThreshold The threshold for the PAO normalization.
   * @param environmentDensities (optional) The density matrices of possible environment
   *        systems.
   */
  PAOController(std::shared_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>> activeDensity,
                std::shared_ptr<MatrixInBasis<Options::SCF_MODES::RESTRICTED>> overlapMatrix, double paoNormalizationThreshold,
                std::vector<std::shared_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>>> environmentDensities = {})
    : _activeDensity(activeDensity),
      _s(overlapMatrix),
      _paoNormalizationThreshold(paoNormalizationThreshold),
      _environmentDensities(environmentDensities) {
    // Even though everything is written for lazy evaluation,
    // I will construct the PAOs here in order to prevent
    // simultaneous access of multiple threads to the same controller.
    // This would lead to a memory leak, because the PAOs are
    // constructed twice (and the copy is not deleted) and probably
    // crash the program.
    constructPAOs();
    getS_PAO();
  };
  virtual ~PAOController() = default;

  /**
   * @brief Selects a set of PAOs, which corresponds to the given domain.
   * @param domain The domain.
   * @return The set of PAO coefficients.
   */
  inline Eigen::MatrixXd getPAOsFromDomain(const Eigen::SparseVector<int>& domain) {
    if (!_paoCoefficients)
      constructPAOs();
    const Eigen::MatrixXd& paoCoefficients = *_paoCoefficients;
    const Eigen::SparseMatrix<double> projection = constructProjectionMatrixFromSparse(domain);
    return (Eigen::MatrixXd)(paoCoefficients * projection).eval();
  }
  /**
   * @brief Returns all PAOs
   * @return All PAO coefficients.
   */
  const CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& getAllPAOs() {
    if (!_paoCoefficients)
      constructPAOs();
    return *_paoCoefficients;
  }
  /**
   * @brief Returns the number of PAOs.
   * @return The number of PAOs controlled by this controller.
   */
  unsigned int getNPAOs() {
    if (!_paoCoefficients)
      constructPAOs();
    return _paoCoefficients->cols();
  }
  /**
   * @brief Getter for the PAO--PAO overlap matrix.
   * @return The PAO--PAO overlap matrix.
   */
  const MatrixInBasis<RESTRICTED>& getS_PAO();

 private:
  /// Construct all PAOs.
  void constructPAOs();
  /// The PAO coefficients.
  std::unique_ptr<CoefficientMatrix<Options::SCF_MODES::RESTRICTED>> _paoCoefficients;
  /// The density matrix of the active system.
  std::shared_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>> _activeDensity;
  /// The overlap matrix of the acitve system.
  std::shared_ptr<MatrixInBasis<Options::SCF_MODES::RESTRICTED>> _s;
  /// The PAO--PAO overlap matrix.
  std::unique_ptr<MatrixInBasis<Options::SCF_MODES::RESTRICTED>> _s_pao;
  /// The PAO normalization threshold.
  double _paoNormalizationThreshold;
  /// The density matrices of the environment systems.
  std::vector<std::shared_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>>> _environmentDensities;
};

} /* namespace Serenity */

#endif /* DATA_PAOCONTROLLER_H_ */
