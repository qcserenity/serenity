/**
 * @file LocalMP2.h
 *
 * @date Dec 7, 2018
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

#ifndef POSTHF_MPN_LOCALMP2_H_
#define POSTHF_MPN_LOCALMP2_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"

namespace Serenity {

/* Forward Declarations */
class DomainOverlapMatrixController;
class OrbitalPair;
class LocalCorrelationController;

/**
 * @class localMP2Settings LocalMP2.h
 * @brief The settings for LocalMP2.
 *
 *    maxResidual                --- Absolute maximum residual in the iterative optimization
 *    maxCycles                  --- Maximum number of amplitude optimization cycles.
 *    ssScaling                  --- Spin component scaling parameters.
 *    osScaling                  --- Spin component scaling parameters.
 *    useFourCenterIntegrals     --- Use full four center integrals (only sensible for testing).
 */
struct LocalMP2Settings {
  double maxResidual = 1e-5;           // Absolute maximum residual in the iterative optimization
  unsigned int maxCycles = 100;        // Maximum number of amplitude optimization cycles.
                                       // Usually 5-10 cycles are enough without the DIIS procedure.
  double ssScaling = 1.0;              // Spin component scaling parameters.
  double osScaling = 1.0;              // Spin component scaling parameters.
  bool useFourCenterIntegrals = false; // Use full four center integrals (only sensible for testing).
};

/**
 * @class LocalMP2 LocalMP2.h
 * @brief Implementation of an orbital invariant MP2 method.\n\n
 *
 * There is currently only a RESTRICTED variant available.\n
 *
 * The implementation of the Local-MP2 method is largely based on\n
 *    J. Chem. Phys. 143, 034108 (2015).\n\n
 *
 * In the local-MP2 method, the amplitudes for the double excitations are optimized iteratively.
 * For each orbital pair ij, with a virtual pair domain [ij] (\f$ \mu,\nu \in [ij] \f$) a residual
 * matrix is calculated:\n
 *        \f$ R^{ij}_{\mu\nu} = K^{ij}_{\mu\nu}+(\epsilon_\mu + \epsilon_\nu - F_{ii} - F{jj})T^{ij}_{\mu\nu}\\
 *                            - \sum_{k \neq i, \kappa\tau} F_{ik}S_{\mu\kappa}T_{\kappa\tau}^{kj}S_{\tau\nu}
 *                            - \sum_{k \neq j, \kappa\tau} F_{kj}S_{\mu\kappa}T_{\kappa\tau}^{ik}S_{\tau\nu}\f$ .\n
 * \f$ K^{ij},~ \epsilon_\mu,~F \f$ denote the matrix of the exchange integrals, expanded over [ij], the eigenvalue \f$
 * \mu \f$ of the Fock matrix in PAO basis and the Fock matrix in the basis of the occupied orbitals, respectively.\n\n
 *
 * The PAOs are selected for each orbital pair, linear dependencies removed via canonical orthogonalization and rotated
 * that they diagonalize the linear independent PAO-Fock matrix.\n\n
 *
 * For DLPNO-MP2 the initial amplitudes are used to construct pair density matrices and then select the most important
 * PNOs which diagonalize the pair density matrix. The PNOs are then rotated again that they diagonalize the fock
 * matrix. This reduces the number of PAOs/PNOs significantly and increases the computational speed.
 *
 * Orbital pairs are precreened based on their differential overlap and an upper bound to the dipole approximation of
 * the pair energies.
 *
 * If you are interesed in the theory, have a look into:\n
 *   J. Chem. Phys. 143, 034108 (2015) (DLPNO-MP2) \n
 *   J. Chem. Phys. 111 5691 (1999)    (PAO based approach by Werner et al.)\n
 *   J. Chem. Phys. 104, 6286 (1996)\n
 *
 */
class LocalMP2 {
 public:
  /**
   * @brief Constructor.
   * @param localCorrelationController The local correlation controller that
   *                                   defines the system.
   */
  LocalMP2(std::shared_ptr<LocalCorrelationController> localCorrelationController)
    : _localCorrelationController(localCorrelationController) {
  }

  /**
   * @brief Destructor.
   */
  virtual ~LocalMP2() = default;

  /**
   * @brief Calculates the energy corrections. Pairs are constructed on the fly.
   * @return The energy corrections order as:\n
   *         LMP2 pair energies,\n
   *         dipole approximation,\n
   *         PNO truncation
   */
  Eigen::VectorXd calculateEnergyCorrection();

  /**
   * @brief Calculates the energy corrections for a given set of pairs..
   * @param pairs The pairs.
   * @return The energy corrections order as:\n
   *         LMP2 pair energies,\n
   *         dipole approximation,\n
   *         PNO truncation
   */
  Eigen::VectorXd calculateEnergyCorrection(std::vector<std::shared_ptr<OrbitalPair>> pairs);

  /**
   * @brief Calculates the unrelaxed density correction.
   * @return Returns the local MP2 correction of the AO density matrix.
   */
  DensityMatrix<RESTRICTED> calculateDensityCorrection();
  /**
   * @brief The settings. More details are above.
   */
  LocalMP2Settings settings;

 private:
  // The local correlation controller.
  std::shared_ptr<LocalCorrelationController> _localCorrelationController;
  // Calculate (ia|jb) integrals for all significant pairs.
  void generateExchangeIntegrals(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs,
                                 std::vector<std::shared_ptr<OrbitalPair>> veryDistantPairs);
  // Set domain overlap matrix controller
  void setDomainOverlapMatrixController(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs);
  // Optimize the amplitudes.
  void optimizeAmplitudes(std::vector<std::shared_ptr<OrbitalPair>> closePairs,
                          std::vector<std::shared_ptr<OrbitalPair>> veryDistantPairs);
  // Calculate the energy.
  Eigen::VectorXd calculateEnergy(std::vector<std::shared_ptr<OrbitalPair>> closePairs,
                                  std::vector<std::shared_ptr<OrbitalPair>> veryDistantPairs);
};

} /* namespace Serenity */

#endif /* POSTHF_MPN_LOCALMP2_H_ */
