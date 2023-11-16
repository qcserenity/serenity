/**
 * @file DLPNO_CCSD.h
 *
 * @date May 16, 2019
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

#ifndef POSTHF_CC_DLPNO_CCSD_H_
#define POSTHF_CC_DLPNO_CCSD_H_

/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h" //Sigma vector in AO basis.
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices
#include <ctime>       //Printing the time for each update cycle.
#include <memory>      //smart ptr
#include <string>      //Printing the time for each update cycle.

namespace H5 {
class H5File;
} // namespace H5
namespace Serenity {

/* Forward Declarations */
class LocalCorrelationController;
class OrbitalPair;
class SingleSubstitution;
class CouplingOrbitalSet;
namespace HDF5 {
using H5File = H5::H5File;
} // namespace HDF5

/**
 * @class DLPNO_CCSD DLPNO_CCSD.h
 * @brief Performs a domain-localized pair-natural orbital (DLPNO) calculation
 *        on the Coupled-Cluster (CC) singles and doubles level.\n\n
 *
 *  DLPNO-CC is a, in principle, linear scaling CC-type method. This implementation is
 *  based on the following publications:\n
 *    [1]J. Chem. Phys. 138, 034106 (2013),\n
 *    [2]J. Chem. Phys. 131, 064103 (2009),\n
 *    [3]J. Chem. Phys. 104, 6286   (1996) and
 *    [4]J. Chem. Phys. 135, 144116 (2011).\n
 *  The equations are taken from Ref. 4 and were adopted for DLPNO-CCSD.\n
 *  The task fulfilled by this class is to optimize the CC-amplitudes and calculate the
 *  corresponding energies for a given set of orbital pairs, which are managed by the
 *  local-correlation controller.\n\n
 *
 *  Note that the equations given in Ref. 1 and 2 are full of typos. The equations which
 *  are implemented are tested carefully vs. canonical CCSD. There should be a pdf/tex
 *  document in the manual-directory with the implemented equations and some background information.
 *
 *  Note that the construction of the singles sigma vector is not fully linearized in this
 *  implementation!
 */
class DLPNO_CCSD {
 public:
  /**
   * @brief Constructor.
   * @param localCorrelationController The local correlation controller, which manages the given orbital
   *                                   pairs and singles.
   * @param maxResidual The maximum entry left in the residual matrices/vectors in order the assume
   *                    the amplitude optimization to be converged.
   * @param maxCycles The maximum number of iterations before the the amplitude optimization is cancelled.
   */
  DLPNO_CCSD(std::shared_ptr<LocalCorrelationController> localCorrelationController, double maxResidual, unsigned int maxCycles);
  /**
   * @brief Default destructor.
   */
  ~DLPNO_CCSD() = default;

  /**
   * @brief Main function. Runs DLPNO-CCSD calculation.
   * @return The DLPNO-CCSD correlation energies.
   */
  Eigen::VectorXd calculateElectronicEnergyCorrections();
  /**
   * @brief Flag for H2 test calculation. The integrals for the one orbital pair and single are
   *        substituted by the full AO four center integrals in order to check the error introduced
   *        by the RI-approximation.
   */
  bool testRun = false;

  struct DLPNO_CCSDSettings {
    /**
     * @brief Flag for skipping the crude SC-MP2 prescreening step in case
     *        they were set previously.
     */
    bool skipCrudePrescreening = false;
    /**
     * @brief Flag for keeping the MO3Center integrals.
     */
    bool keepMO3CenterIntegrals = false;
    /**
     * @brief Flag for keeping the integrals stored for each pair.
     */
    bool keepPairIntegrals = false;
  };
  /**
   * @brief DLPNO-CCSD specific settings.
   */
  DLPNO_CCSDSettings dlpnoCCSDSettings;

 private:
  /**
   * @brief The local-correlation controller. (The main hub for all local methods)
   */
  std::shared_ptr<LocalCorrelationController> _localCorrelationController;
  /**
   * @brief Maximum residual tolerated by amplitude optimization.
   */
  const double _maxResidual;
  /**
   * @brief Maximum number of cycles before the amplitude optimization cancels.
   */
  const unsigned int _maxCycles;
  /**
   * @brief Sigma vector in AO basis.
   */
  std::shared_ptr<MatrixInBasis<RESTRICTED>> _g_ao_ao;
  /**
   * @brief Sigma vector in mixed basis.
   */
  Eigen::MatrixXd _g_ao_occ;
  /**
   * @brief Sigma vector in the basis of the occupied orbitals.
   */
  Eigen::MatrixXd _g_occ_occ;
  /**
   * @brief If true, the integrals for the sigma vector construction are directly calculated in their PNO basis and
   * cached.
   */
  const bool _linearScalingSigmaVector;

  /**
   * @brief Time of the start of the last iteration.
   */
  timespec _time;
  /*
   * tilde X --> "dressed" X
   * tau     --> tau^ij_ab = t_ab^ij+(t_a^i*t_b^j)
   * t_ij    --> double substitution amplitudes
   * t_i     --> single substitution amplitudes
   *
   * All equations are given in the manual. Note that the indices of the equations refere
   * to the last digit: e.g. Eq. (3) may be Eq. 6.2.1.3 in the manual.
   * The equations are adopted (with slightly different notation) from the appendix of [4].
   */
  ///@brief Implements Eq. (3).
  inline void calculate_tau_ij(std::shared_ptr<OrbitalPair> pair);
  ///@brief Intermediate for pair--pair interaction. Implements Eq. (12).
  inline void calculate_Y(std::shared_ptr<OrbitalPair> pair);
  ///@brief Intermediate for pair--pair interaction. Implements Eq. (13).
  inline void calculate_Z(std::shared_ptr<OrbitalPair> pair);
  ///@brief Implements Eq. (11).
  inline void calculate_tilde_ikjl(std::shared_ptr<OrbitalPair> pair);
  ///@brief Implements Eq. (10).
  inline Eigen::MatrixXd calculate_K_tau(std::shared_ptr<OrbitalPair> pair);
  ///@brief Implements Eq. (9) for pq = ia.
  inline Eigen::VectorXd calculate_G_t1_ia(std::shared_ptr<SingleSubstitution> single);
  ///@brief Implements Eq. (9) for pq = ij.
  inline std::pair<double, double> calculate_G_t1_ij(std::shared_ptr<OrbitalPair> pair);
  ///@brief Implements Eq. (9) for pq = ab.
  inline Eigen::MatrixXd calculate_G_t1_ab(std::shared_ptr<OrbitalPair> pair);
  ///@brief Implements Eq. (4).
  inline void calculate_tildeF_ij(std::shared_ptr<OrbitalPair> pair);
  ///@brief Implements Eq. (5).
  inline void calculate_tildeF_ab(std::shared_ptr<OrbitalPair> pair);
  ///@brief Implements Eq. (6).
  inline void calculate_tildeF_ai(std::shared_ptr<SingleSubstitution> single);
  ///@brief Implements Eq. (7).
  inline void calculate_tildetildeF_ij(std::shared_ptr<OrbitalPair> pair);
  ///@brief Implements Eq. (8).
  inline Eigen::MatrixXd calculate_tildetildeF_ab(std::shared_ptr<OrbitalPair> pair);
  ///@brief Calculates the singles residuals for the given single. Implements Eq. (1).
  inline Eigen::VectorXd calculateSinglesResidual(std::shared_ptr<SingleSubstitution> single);
  ///@brief Calculates the double residuals for the given pair. Implements Eq. (2).
  inline Eigen::MatrixXd calculateDoublesResidual(std::shared_ptr<OrbitalPair> pair);
  ///@brief Optimizes the amplitudes for the pairs and singles.
  inline void optimizeAmplitudes();
  ///@brief Calculates the single fock matrix vectors tildeF_ij for all singles.
  inline void dressSingles();
  ///@brief Calculates the dressed quantities for the orbital pairs.
  inline void dressPairs();
  ///@brief Prints the singles norm, T1-diagnostic and the maximum singles coefficient to the output.
  void runDiagnostics();
  /**
   * @brief Prepare the orbital pairs for the following calculations.
   *
   *      Sets up: SC-MP2 amplitudes.
   *      PNO-Construction.
   *      Selection of distant and close orbital pairs.
   *      Pair-Coupling map.
   *      kl-Pairs assignment.
   *      PNO-Domain overlap integrals.
   *      Initialization of singles.
   */
  inline void prepareOrbitalPairs();
  /**
   * @brief Calculates the integrals needed for DLPNO-CCSD directly into
   *        the PNO-basis for each pair.
   */
  inline void calculateCCSDIntegrals();
  /**
   * @brief Calculates the DLPNO-CCSD energies for the pairs/singles with the
   *        current set of integrals and amplitudes.
   * @return The total DLPNO-CCSD energy corrections.
   *         Order: Pair energies, diagonal singles term, SC-MP2 correction,
   *                dipole correction, PNO truncation correction.
   */
  inline Eigen::VectorXd calculateEnergyCorrection();
  /**@brief Only used for tests. Switches the RI approximated integrals for the exact
   *        full four center integrals and forces the pairs to use the canonical virtual
   *        orbitals as their domain. Thus, the code will perform a canonical CCSD calculation.
   */
  void switchIntegrals();
  /**
   * @brief Deletes the pair integral files if they were written and are not supposed to be kept.
   */
  void deleteIntegralFiles();
  /**
   * @brief Removes the integrals stored in the orbital pairs.
   */
  void deleteIntegrals();
  /**
   * @brief Updates the sigma vector.
   */
  void updateSigmaVector();
  /**
   * @brief Getter for the time needed for the CCSD iteration.
   * @return The time as a string.
   */
  std::string getTimeString();
  /**
   * @brief Try to load as many integrals as possible. Some memory was probably freed previeosly.
   * @return A pointer to the open HDF5 file, if not all integrals can be loaded. Nullptr otherwise.
   */
  std::shared_ptr<HDF5::H5File> tryLoadingIntegrals();
};

} /* namespace Serenity */

#endif /* POSTHF_CC_DLPNO_CCSD_H_ */
