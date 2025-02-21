/**
 * @file IAOPopulationCalculator.h
 *
 * @date Oct 15, 2018
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

#ifndef ANALYSIS_POPULATIONANALYSIS_IAOPOPULATIONCALCULATOR_H_
#define ANALYSIS_POPULATIONANALYSIS_IAOPOPULATIONCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h" //CoefficientMatrix.
#include "data/matrices/SPMatrix.h"          //SPMatrix.
/* Include Std and External Headers */
#include <memory> //smrt_ptr

namespace Serenity {

/* Forward Declarations */
class BasisController;
class AtomCenteredBasisController;
class SystemController;
class Geometry;
/**
 * @class IAOPopulationCalculator IAOPopulationCalculator.h
 * @brief A class to calculate the IAO coefficients and corresponding populations for any system.
 *
 * According to:\n
 *   Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical Concepts
 *   Gerald Knizia, J. Chem. Theory Comput., 2013, 9 (11), pp 4834–4843\n
 *   (better look it up on Knizia's homepage. The published article has some errors in the appendix.)\n\n
 *
 *   The IAO populations for an atom  \f$ A  \f$ are defined as:\n
 *      \f$ q_A = \sum_{\rho \in A} \langle\rho|\gamma|\rho\rangle~,\f$\n
 *   where  \f$ \gamma=k\sum_i|i\rangle\langle i| \f$ is the density matrix of the system with the occupied orbitals \f$
 * \{|i\rangle\} \f$ and the number of electrons per orbital  \f$ k  \f$. The set \f$ \{ |\rho\rangle \} \f$ is an
 * orthonormal basis which spans the space of the occupied orbitals while referencing a specific atomic orbital for each
 * \f$ |\rho\rangle \f$ . The final formula for \f$ q_A \f$ can be obtained by expanding the  \f$ |i\rangle \f$ in terms
 * of \f$ \{ |\rho\rangle \} \f$ . The expansion coefficients for this are called CIAO coefficients in this
 * implementation.
 */
template<Options::SCF_MODES SCFMode>
class IAOPopulationCalculator {
 private:
  /* purely static. Never instantiated. */
  IAOPopulationCalculator() = default;
  virtual ~IAOPopulationCalculator() = default;

 public:
  /**
   * @brief Calculates the IAO populations for a given system.
   * @param system The system.
   * @return The IAO populations.
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> static calculateIAOPopulations(std::shared_ptr<SystemController> system);
  /**
   * @brief Calculates the IAO populations on each atom orbital wise.
   * @param C The coefficient matrix.
   * @param S1 The orbital basis overlap integrals.
   * @param nOccOrbs The number of occupied orbitals.
   * @param B1 The basis controller of the system.
   * @param B2 The basis controller used for the IAO analysis.
   * @param geom The geometry.
   * @param withVirtuals  If true, the populations for virtual orbitals are calculated up to the maximum size in the
   *                      minimal IAO basis set.
   * @return The orbital populations on the system. No weighting for spin-polarization included.
   */
  SPMatrix<SCFMode> static calculateAtomwiseOrbitalPopulations(
      const CoefficientMatrix<SCFMode>& C, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& S1,
      const SpinPolarizedData<SCFMode, unsigned int> nOccOrbs, std::shared_ptr<BasisController> B1,
      std::shared_ptr<AtomCenteredBasisController> B2, const std::shared_ptr<Geometry> geom, bool withVirtuals = false);
  /**
   * @brief Calculates the IAO populations on each atom orbital wise.
   * @param system        The system controller.
   * @param withVirtuals  If true, the populations for virtual orbitals are calculated up to the maximum size in the
   *                      minimal IAO basis set.
   * @return The orbital populations on the system. No weighting for spin-polarization included.
   */
  SPMatrix<SCFMode> static calculateAtomwiseOrbitalPopulations(std::shared_ptr<SystemController> system,
                                                               bool withVirtuals = false);
  /**
   * @brief Calculate the shellwise orbital populations.
   * @param C The coefficient matrix.
   * @param S1 The orbital basis overlap integrals.
   * @param nOccOrbs The number of occupied orbitals.
   * @param B1 The basis controller of the system.
   * @param B2 The basis controller used for the IAO analysis.
   * @param withVirtuals  If true, the populations for virtual orbitals are calculated up to the maximum size in the
   *                      minimal IAO basis set.
   * @return The orbital populations.
   */
  SPMatrix<SCFMode> static calculateShellwiseOrbitalPopulations(const CoefficientMatrix<SCFMode>& C,
                                                                const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& S1,
                                                                const SpinPolarizedData<SCFMode, unsigned int> nOccOrbs,
                                                                std::shared_ptr<BasisController> B1,
                                                                std::shared_ptr<AtomCenteredBasisController> B2,
                                                                bool withVirtuals = false);
  /**
   * @brief Calculate the shellwise orbital populations.
   * @param system        The system controller.
   * @param withVirtuals  If true, the populations for virtual orbitals are calculated up to the maximum size in the
   *                      minimal IAO basis set.
   * @return The orbital populations.
   */
  SPMatrix<SCFMode> static calculateShellwiseOrbitalPopulations(std::shared_ptr<SystemController> system,
                                                                bool withVirtuals = false);
  /**
   * @brief Calculates the coefficients for the expansion of the orbitals in the IAO basis.
   * @param system The system controller.
   * @return The expansion coefficients in the IAO basis.
   */
  std::pair<SPMatrix<SCFMode>, SPMatrix<SCFMode>> static getCIAOCoefficients(std::shared_ptr<SystemController> system,
                                                                             bool withVirtuals);
  /**
   * @brief Calculates the coefficients for the expansion of the orbitals in the IAO basis.
   * @param C The coefficient matrix.
   * @param S1 The orbital basis overlap integrals.
   * @param nOccOrbs The number of occupied orbitals.
   * @param B1 The basis controller of the system.
   * @param B2 The basis controller used for the IAO analysis.
   * @param withVirtuals If true, the coefficients for the virtual orbitals are calculated as well.
   * @return The expansion coefficients in the IAO basis.
   */
  std::pair<SPMatrix<SCFMode>, SPMatrix<SCFMode>> static getCIAOCoefficients(
      const CoefficientMatrix<SCFMode>& C, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& S1,
      const SpinPolarizedData<SCFMode, unsigned int> nOccOrbs, std::shared_ptr<BasisController> B1,
      std::shared_ptr<BasisController> B2, bool withVirtuals = false);
  /**
   * @brief Calculates the orbital-wise 1S populations.
   * @param C The coefficient matrix.
   * @param S1 The orbital basis overlap integrals.
   * @param nOccOrbs The number of occupied orbitals.
   * @param B1 The basis controller of the system.
   * @param B2 The basis controller used for the IAO analysis.
   * @param geom The geometry.
   * @return The orbital populations.
   */
  SPMatrix<SCFMode> static calculate1SOrbitalPopulations(const CoefficientMatrix<SCFMode>& C,
                                                         const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& S1,
                                                         const SpinPolarizedData<SCFMode, unsigned int> nOccOrbs,
                                                         std::shared_ptr<BasisController> B1,
                                                         std::shared_ptr<AtomCenteredBasisController> B2,
                                                         const std::shared_ptr<Geometry> geom);
  /**
   * @brief Calculates the orbital-wise 1S populations.
   * @param system The system controller.
   * @return The orbital populations.
   */
  SPMatrix<SCFMode> static calculate1SOrbitalPopulations(std::shared_ptr<SystemController> system);
  /**
   * @brief Reconstruct the virtual valence orbitals to exactly span the virtual space of the IAOs.
   *        All remaining virtual orbitals are reconstructed by projection.
   * @param C The coefficient matrix.
   * @param S1 The orbital basis overlap integrals.
   * @param nOccOrbs The number of occupied orbitals.
   * @param B1 The basis controller of the system.
   * @param B2 The basis controller used for the IAO analysis.
   * @param removeSOMO If true, it is assumed that the orbitals are (quasi) restricted. The SOMO will then be separated
   * from the other occupied orbitals.
   */
  void static reconstructVirtualValenceOrbitalsInplace(CoefficientMatrix<SCFMode>& C,
                                                       const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& S1,
                                                       const SpinPolarizedData<SCFMode, unsigned int> nOccOrbs,
                                                       std::shared_ptr<BasisController> B1,
                                                       std::shared_ptr<BasisController> B2);
  /**
   * @brief Check whether the IAOs span at least the space of the virtual valence orbitals and all occupied orbitals.
   * @param C The coefficient matrix.
   * @param S1 The orbital basis overlap integrals.
   * @param nOccOrbs The number of occupied orbitals.
   * @param B1 The basis controller of the system.
   * @param B2 The basis controller used for the IAO analysis.
   * @return True, if the IAOs span the space of the orbitals. Otherwise, false.
   */
  static bool iaosSpanOrbitals(const CoefficientMatrix<SCFMode>& C, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& S1,
                               const SpinPolarizedData<SCFMode, unsigned int> nOccOrbs,
                               std::shared_ptr<BasisController> B1, std::shared_ptr<BasisController> B2);
};

} /* namespace Serenity */

#endif /* ANALYSIS_POPULATIONANALYSIS_IAOPOPULATIONCALCULATOR_H_ */
