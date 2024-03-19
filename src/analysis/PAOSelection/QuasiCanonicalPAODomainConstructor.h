/**
 * @file QuasiCanonicalPAODomainConstructor.h
 *
 * @date May 13, 2019
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

#ifndef ANALYSIS_PAOSELECTION_QUASICANONICALPAODOMAINCONSTRUCTOR_H_
#define ANALYSIS_PAOSELECTION_QUASICANONICALPAODOMAINCONSTRUCTOR_H_

/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/FockMatrix.h"
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices.
#include <memory>      //smart ptr.
#include <vector>      //std::vector.

namespace Serenity {
/* Forward Declarations */
class OrbitalPair;
class PAOController;
class SystemController;

/**
 * @brief A class that constructs a set of quasi-canonical, non-redundant PAO for a given pair and the
 *        corresponding SC-MP2 amplitudes.
 *        The PAO domain within the pair needs to be initialized.\n
 *        Furthermore, this class calculates the overlap between the PAO domains defined in the pair
 *        or single.\n\n
 *
 *        Quasi-canonical, non-redundant PAOs denote a set of linear independent PAOs that diagonalize
 *        the given Fock matrix.
 */
class QuasiCanonicalPAODomainConstructor {
 public:
  /**
   * @brief Constructor.
   * @param coefficients The coefficient matrix.
   * @param overlapMatrix The AO overlap matrix.
   * @param f The Fock matrix.
   * @param paoController The PAO controller.
   * @param paoOrthogonalizationThreshold The norm threshold for the canonical orthogonalization of the PAOs.
   * @param environmentSystems Possible environment systems which are shifted in the fock matrix.
   * @param levelShiftParameter The shift of the occupied environment orbitals in the fock matrix.
   * @param ssScaling Same-spin scaling factor.
   * @param osScaling Opposite-spin scaling factor.
   * @param clear Clear integrals after pair energy calculation.
   */
  QuasiCanonicalPAODomainConstructor(const CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& coefficients,
                                     std::shared_ptr<const FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                                     std::shared_ptr<PAOController> paoController, double paoOrthogonalizationThreshold,
                                     std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                     double levelShiftParameter, double ssScaling = 1.0, double osScaling = 1.0,
                                     bool clear = false);

  /**
   * @brief Transforms the PAO basis of the given pair.
   * @param pair The orbital pair.
   */
  virtual void transformExternalBasis(std::shared_ptr<OrbitalPair> pair) {
    transformToQuasiCanonicalPAOBasis(pair);
  }
  /**
   * @brief Calculates the semi-canonical MP2 pair energy and initializes amplitudes.
   * @param pair The orbital pair.
   */
  virtual void postProcessing(std::shared_ptr<OrbitalPair> pair) {
    initializeAmplitudes(pair);
    if (_clear)
      clearIntegrals(pair);
  }
  /**
   * @brief Default destructor.
   */
  virtual ~QuasiCanonicalPAODomainConstructor() = default;

 protected:
  /**
   * @brief Transforms the PAO basis to its quasi-canonical, non-redundant representation.
   * @param pair The orbital pair.
   */
  void transformToQuasiCanonicalPAOBasis(std::shared_ptr<OrbitalPair> pair);
  /**
   * @brief Initializes amplitudes and calculate SC-MP2 pair energy.
   * @param pair The orbital pair.
   */
  void initializeAmplitudes(std::shared_ptr<OrbitalPair> pair);
  /**
   * @brief Clears the K_ij integrals of the given pair.
   * @param pair The orbital pair.
   */
  void clearIntegrals(std::shared_ptr<OrbitalPair> pair);
  ///@brief The internal block of the Fock matrix.
  Eigen::MatrixXd _f_MO;
  ///@brief Fock matrix with the right function transformed to the internal basis.
  Eigen::MatrixXd _f_PAO_MO;
  ///@brief The Fock matrix in PAO basis.
  Eigen::MatrixXd _f_pao;
  ///@brief The PAO controller.
  std::shared_ptr<PAOController> _paoController;

 private:
  /**
   * @brief The norm-threshold for canonical orthogonalization of the PAO domain.
   */
  double _paoOrthogonalizationThreshold;
  /**
   * @brief The coefficient matrix of the system (occupied and virtual).
   */
  Eigen::MatrixXd _coeff;
  // Same-spin scaling parameter.
  double _ssScaling;
  // Opposite-spin scaling parameter.
  double _osScaling;
  // Clear QC-PAO K_ij integrals.
  bool _clear;
};

} /* namespace Serenity */

#endif /* ANALYSIS_PAOSELECTION_QUASICANONICALPAODOMAINCONSTRUCTOR_H_ */
