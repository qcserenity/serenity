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
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Options.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
/**
 * @class ExcitationSpectrum
 *
 * Prints the oscillator and rotatory strengths and corresponding transition dipole moments
 * for the excitation energies obtained earlier.
 */
class ExcitationSpectrum {
 public:
  /**
   * @brief Calculates oscillator and rotatory strengths and transition dipole moments according to
   *        J. Am. Chem. Soc. 122, 1717 (2000)
   * @param dipoles All property integrals needed (electric length/velocity and magnetic).
   * @param eigenvecors Eigenvectors of the response problem.
   * @param eigenvalues Eigenvalues of the response problem.
   */
  static void printSpectrum(const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                            const std::vector<Eigen::MatrixXd>& eigenvectors, const Eigen::VectorXd& eigenvalues);
};

} /* namespace Serenity */
#endif /* LRSCF_EXCITATIONSPECTRUM */
