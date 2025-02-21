/**
 * @file CC2HelperFunctions.h
 * @author Niklas Niemeyer
 *
 * @date Jul 25, 2023
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

#ifndef LRSCF_CC2HELPERFUNCTIONS
#define LRSCF_CC2HELPERFUNCTIONS

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Tools/SigmaCalculator.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class DipoleIntegrals;

template<Options::SCF_MODES SCFMode>
/**
 * @class CC2HelperFunctions
 *
 * This is a class to clean up the CC2 part of the LRSCF Task so it does not look as messy.
 */
class CC2HelperFunctions {
 public:
  /**
   * @brief Calculates CIS eigenvectors as an initial guess.
   *
   * @param settings The LRSCF task settings.
   * @param eigenvectors The eigenvector pointer to pass the eigenvector solution to.
   * @param eigenvalues The vector to pass the eigenvalue solution to.
   * @param nDimension The dimension of the LRSCF problem (nv*no).
   * @param diagonal The orbital-energy differences.
   * @param initialSubspace The initial subspace size.
   * @param sigmaCalculator The CIS sigma vector calculation.
   */
  static void calculateCISEigenvectors(LRSCFTaskSettings& settings, std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                       Eigen::VectorXd& eigenvalues, unsigned nDimension, const Eigen::VectorXd& diagonal,
                                       unsigned initialSubspace, SigmaCalculator sigmaCalculator);

  /**
   * @brief Calculates right eigenpairs of the CC2 Jacoby matrix.
   *
   * @param settings The LRSCF task settings.
   * @param eigenvectors The eigenvector pointer to pass the eigenvector solution to.
   * @param eigenvalues The vector to pass the eigenvalue solution to.
   * @param diagonal The orbital-energy differences.
   * @param sigmaCalculator The right CC2 sigma vector calculation.
   * @param type Type of the LRSCF problem (iso, FDEu, FDEc).
   * @param writeToDisk Disk writer function to write temporary solution to disk.
   */
  static void calculateRightEigenvectors(LRSCFTaskSettings& settings, std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                         Eigen::VectorXd& eigenvalues, const Eigen::VectorXd& diagonal,
                                         NonlinearSigmaCalculator sigmaCalculator, Options::LRSCF_TYPE type,
                                         std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk);

  /**
   * @brief Calculates left eigenpairs of the CC2 Jacoby matrix.
   *
   * @param settings The LRSCF task settings.
   * @param eigenvectors The eigenvector pointer to pass the eigenvector solution to.
   * @param eigenvalues The vector to pass the eigenvalue solution to.
   * @param diagonal The orbital-energy differences.
   * @param sigmaCalculator The left CC2 sigma vector calculation.
   * @param type Type of the LRSCF problem (iso, FDEu, FDEc).
   */
  static void calculateLeftEigenvectors(LRSCFTaskSettings& settings, std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                        Eigen::VectorXd& eigenvalues, const Eigen::VectorXd& diagonal,
                                        NonlinearSigmaCalculator sigmaCalculator, Options::LRSCF_TYPE type);

  /**
   * @brief Normalizes the right eigenvectors (approximate normalization because the left eigenvectors are not
   * available).
   *
   * @param lrscf The LRSCF controllers.
   * @param settings The LRSCF task settings.
   * @param eigenvectors The eigenvectors (only right).
   * @param eigenvalues The eigenvalues/excitation energies.
   */
  static void normalizeRightEigenvectors(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                         const LRSCFTaskSettings& settings,
                                         std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                         Eigen::VectorXd& eigenvalues);

  /**
   * @brief Normalizes the right and left eigenvectors.
   *
   * @param lrscf The LRSCF controllers.
   * @param settings The LRSCF task settings.
   * @param eigenvectors The eigenvectors.
   * @param eigenvalues The eigenvalues/excitation energies.
   */
  static void normalizeBothEigenvectors(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                        const LRSCFTaskSettings& settings,
                                        std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                        Eigen::VectorXd& eigenvalues);

  /**
   * @brief Gets the dimensions of the objects to store the solutions in right.
   *
   * @param lrscf The LRSCF controllers.
   * @param settings The LRSCF task settings.
   * @param eigenvectors The eigenvectors (stores also the multipliers afterwards).
   * @param transitiondensities The transition and state densities.
   */
  static void prepareVectors(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                             const LRSCFTaskSettings& settings, std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                             std::shared_ptr<std::vector<Eigen::MatrixXd>>& transitiondensities);

  /**
   * @brief Calculates state multipliers (ground state and all excited states).
   *
   * @param lrscf The LRSCF controllers.
   * @param settings The LRSCF task settings.
   * @param eigenvectors The eigenvectors.
   * @param eigenvalues The eigenvalues/excitation energies.
   * @param sigmaCalculator The left CC2 sigma vector calculator.
   * @param diagonal The orbital-energy differences.
   */
  static void calculateStateMultipliers(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                        const LRSCFTaskSettings& settings,
                                        std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors, Eigen::VectorXd& eigenvalues,
                                        NonlinearSigmaCalculator sigmaCalculator, const Eigen::VectorXd& diagonal);

  /**
   * @brief Calculates state densities (ground state and all excited states).
   *
   * @param lrscf The LRSCF controllers.
   * @param transitiondensities The object to store the densities in.
   * @param eigenvectors The eigenvectors.
   * @param eigenvalues The eigenvalues/excitation energies.
   */
  static void calculateStateDensities(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                      std::shared_ptr<std::vector<Eigen::MatrixXd>>& transitiondensities,
                                      std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                      Eigen::VectorXd& eigenvalues);

  /**
   * @brief  Calculates transition multipliers (all excited states).
   *
   * @param lrscf The LRSCF controllers.
   * @param settings The LRSCF task settings.
   * @param eigenvectors The eigenvectors.
   * @param eigenvalues The eigenvalues/excitation energies.
   * @param sigmaCalculator The left CC2 sigma vector calculator.
   * @param diagonal The orbital-energy differences.
   */
  static void calculateTransitionMultipliers(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                             const LRSCFTaskSettings& settings,
                                             std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                             Eigen::VectorXd& eigenvalues, NonlinearSigmaCalculator sigmaCalculator,
                                             const Eigen::VectorXd& diagonal);

  /**
   * @brief Calculates transition densities (all excited states).
   *
   * @param lrscf The LRSCF controllers.
   * @param transitiondensities The object to store the densities in.
   * @param eigenvectors The eigenvectors.
   * @param eigenvalues The eigenvalues/excitation energies.
   */
  static void calculateTransitionDensities(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                           std::shared_ptr<std::vector<Eigen::MatrixXd>>& transitiondensities,
                                           std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                           Eigen::VectorXd& eigenvalues);

  /**
   * @brief Calculates perturbed amplitudes as solution to the CC2 response equation.
   *
   * @param lrscf The LRSCF controllers.
   * @param settings The LRSCF task settings.
   * @param solutionvectors The object to store the amplitudes in.
   * @param frequencies The frequencies to compute the perturbed amplitudes for.
   * @param sigmaCalculator The left CC2 sigma vector calculator.
   * @param diagonal The orbital-energy differences.
   * @param dipoles The dipole integrals: electric (three left columns, either length or velocity) and magnetic (three
   * right columns).
   */
  static void calculatePerturbedAmplitudes(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                           const LRSCFTaskSettings& settings,
                                           std::shared_ptr<std::vector<Eigen::MatrixXd>>& solutionvectors,
                                           const std::vector<double>& frequencies, NonlinearSigmaCalculator sigmaCalculator,
                                           const Eigen::VectorXd& diagonal, Eigen::MatrixXd& dipoles);

  /**
   * @brief Calculated perturbed densities and F-matrix contractions.
   *
   * @param lrscf The LRSCF controllers.
   * @param settings The LRSCF task settings.
   * @param solutionvectors The object to get the perturbed amplitudes from.
   * @param perturbeddensities The object to store the perturbed densities in.
   * @param frequencies The frequencies to compute the perturbed densities for.
   * @param dipoles The dipole integrals: electric (three left columns, either length or velocity) and magnetic (three
   * right columns).
   * @param Fdipdip F matrix contractions (both operators electric dipole operators).
   * @param Fdipmag F matrix contractions (one operator magnetic dipole, one operator electric dipole operators).
   */
  static void calculatePerturbedDensities(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                          const LRSCFTaskSettings& settings,
                                          std::shared_ptr<std::vector<Eigen::MatrixXd>>& solutionvectors,
                                          std::shared_ptr<std::vector<Eigen::MatrixXd>>& perturbeddensities,
                                          const std::vector<double>& frequencies, Eigen::MatrixXd& dipoles,
                                          std::vector<Eigen::Matrix3d>& Fdipdip, std::vector<Eigen::Matrix3d>& Fdipmag);
};

} /* namespace Serenity */
#endif /* LRSCF_CC2HELPERFUNCTIONS */
