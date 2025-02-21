/**
 * @file ExcitationSpectrum.h
 *
 * @date Dec. 17, 2018
 * @author Michael Boeckers
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

#ifndef LRSCF_EXCITATIONSPECTRUM
#define LRSCF_EXCITATIONSPECTRUM

/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/LRSCFOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class DipoleIntegrals;

/**
 * @class ExcitationSpectrum
 *
 * Prints the oscillator and rotatory strengths and corresponding transition dipole moments
 * for the excitation energies obtained earlier.
 */
template<Options::SCF_MODES SCFMode>
class ExcitationSpectrum {
 public:
  /**
   * @brief Calculates oscillator and rotatory strengths and transition dipole moments according to
   *        J. Am. Chem. Soc. 122, 1717 (2000)
   * @param method The linear-response method in use.
   * @param dipoles All property integrals needed (electric length/velocity and magnetic).
   * @param eigenvectors Eigenvectors of the response problem.
   * @param eigenvalues Eigenvalues of the response problem.
   * @param results A simple matrix to store eigenvalues, osc. and rot. strengths in.
   * @param fileName A string for the name of the file in which transition moments are to be stored in.
   */
  static void printTransitionMoments(Options::LR_METHOD method, const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                     const std::vector<Eigen::MatrixXd>& eigenvectors,
                                     const Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& results, std::string fileName);

  /**
   * @brief Calculates excited-state densities and dipole moments (currently only implemented for CC2).
   * @param method The linear-response method in use.
   * @param dipoles All property integrals needed (electric length/velocity and magnetic).
   * @param eigenvectors Eigenvectors of the response problem.
   * @param eigenvalues Eigenvalues of the response problem.
   * @param results A simple matrix to store eigenvalues, osc. and rot. strengths in.
   * @param nuclearPart The nulcear part of the ground-state dipole moment.
   */
  static void printStateMoments(Options::LR_METHOD method, const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                const std::vector<Eigen::MatrixXd>& eigenvectors, const Eigen::VectorXd& eigenvalues,
                                Eigen::MatrixXd& results, Eigen::Vector3d nuclearPart);

  static void printStateMoments(Options::LR_METHOD method, const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                const std::vector<MatrixInBasis<SCFMode>>& unrelaxedDensities,
                                const std::vector<MatrixInBasis<SCFMode>>& relaxedDensities, const Eigen::VectorXd& eigenvalues,
                                Eigen::Vector3d nuclearPart, SpinPolarizedData<SCFMode, unsigned int> nOcc);
};

} /* namespace Serenity */
#endif /* LRSCF_EXCITATIONSPECTRUM */
