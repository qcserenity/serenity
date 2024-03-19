/**
 * @file PNOConstructor.h
 *
 * @date Apr 11, 2019
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

#ifndef ANALYSIS_PAOSELECTION_PNOCONSTRUCTOR_H_
#define ANALYSIS_PAOSELECTION_PNOCONSTRUCTOR_H_

/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/QuasiCanonicalPAODomainConstructor.h" //QCPAO constructor, inherits from this class.
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/FockMatrix.h"
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices.
#include <memory>      //smart ptr.
#include <vector>      //std::vector

namespace Serenity {

/* Forward Declarations */
class OrbitalPair;
class PAOController;
class SystemController;

/**
 * @class PNOConstructor PNOConstructor.h
 * @brief Calculates the PNOs from initial amplitudes for the given pair and truncates them according to their
 * occupation. The PNOs are chosen to diagonalize their respective Fock matrix block. Amplitudes, exchange integrals,
 * transformation to PAO/PNO domain and uncoupled residual elements are updated for the pair in place. If no PNOs
 * survive the truncation procedure, the type of the pair is set to "VERY_DISTANT".\n See the following references if
 * you want to learn more:\n J. Chem. Phys. 143, 034108 (2015)\n J. Chem. Phys. 144, 024109 (2016)\n For the
 * contravariant amplitude (\f$ \tilde{T}^{ij} \f$) definition:\n J Chem. Phys. 139, 134101 (2013)\n
 *
 */
class PNOConstructor : public QuasiCanonicalPAODomainConstructor {
 public:
  /**
   * @brief Default destructor.
   */
  virtual ~PNOConstructor() = default;
  /**
   * @brief Constructor.
   * @param coefficients The coefficient matrix.
   * @param f The Fock matrix.
   * @param paoController The PAO-Controller.
   * @param paoOrthogonalizationThreshold Threshold for the canonical PAO orthogonalisation.
   * @param environmentSystems Optional environment systems.
   * @param levelShiftParameter Optional level-shift value for the occupied environment orbitals.
   * @param setFaiZero Force any F_ai block to be zero.
   * @param ssScaling Same spin scaling factor.
   * @param osScaling Opposite spin scaling factor.
   */
  PNOConstructor(const CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& coefficients,
                 std::shared_ptr<const FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                 std::shared_ptr<PAOController> paoController, double paoOrthogonalizationThreshold,
                 std::vector<std::shared_ptr<SystemController>> environmentSystems = {}, double levelShiftParameter = 1e+6,
                 bool setFaiZero = false, double ssScaling = 1.0, double osScaling = 1.0);

  /**
   * @brief Transform the external basis of the given pair.
   * @param orbitalPair The pair.
   */
  virtual void transformExternalBasis(std::shared_ptr<OrbitalPair> orbitalPair) override final;
  /**
   * @brief Constructs PNOs and transforms amplitudes and exchange integrals to the PNO basis.
   * @param pair The pair.
   */
  virtual void postProcessing(std::shared_ptr<OrbitalPair> pair) override;
  /**
   * @brief Transforms the PAO basis of the given pairs to their significant (quasi-canonical) PNO basis.
   * @param orbitalPairs  The orbital pairs.
   */
  void transformToPNOBasis(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs);

  /**
   * @brief Transforms the PAO basis of the given pair to its significant (quasi-canonical) PNO basis.
   * @param orbitalPairs  The orbital pair.
   */
  void transformToPNOBasis(std::shared_ptr<OrbitalPair> pair);

 private:
  ///@brief Force the F_ai block of the fock matrix to be zero.
  bool _setFaiZero;
  ///@brief Same-spin scaling parameter.
  double _ssScaling;
  ///@brief opposite-spin scaling parameter.
  double _osScaling;

  /**
   * @brief Calculate quasi canonical PNOs for the given pair.
   * @param d_ij_red The MP2 one particle density matrix.
   * @param pair The pair.
   * @return The eigenvalues and set of eigenvectors.
   */
  std::pair<Eigen::VectorXd, Eigen::MatrixXd> orthogonalizeFockMatrix(const Eigen::MatrixXd& d_ij_red,
                                                                      std::shared_ptr<OrbitalPair> pair);
};

} /* namespace Serenity */

#endif /* ANALYSIS_PAOSELECTION_PNOCONSTRUCTOR_H_ */
