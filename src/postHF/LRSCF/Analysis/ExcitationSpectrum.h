/**
 * @file ExcitationSpectrum.h
 *
 * @date Nov 02, 2017
 * @author Michael Boeckers
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef POSTHF_LRSCF_ANALYSIS_EXCITATIONSPECTRUM_H_
#define POSTHF_LRSCF_ANALYSIS_EXCITATIONSPECTRUM_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES T> class ExcitationSpectrum {
public:
  /**
   * @brief Calculates oscillator strengths and transition dipole moments according to
   *        J. Am. Chem. Soc. 122, 1717 (2000)
   * @param systemController
   * @param eigenvectors
   * @param eigenvalues
   */
  ExcitationSpectrum(
      std::shared_ptr<SystemController> systemController,
      std::vector<Eigen::MatrixXd >& eigenvectors,
      Eigen::VectorXd& eigenvalues);

  virtual ~ExcitationSpectrum() = default;

  /**
   * @brief Prints excitation spectrum
   */
  void printSpectrum();
private:
  //Calculates MO dipole integrals
  Eigen::MatrixXd dipoleIntegrals();

  //Calculates MO momentum integrals
  Eigen::MatrixXd momentumIntegrals();

  //Transforms AO matrix to MO vector
  Eigen::MatrixXd ao2mo(std::vector<Eigen::MatrixXd>& ao_xyz);

  //The system controller
  std::shared_ptr<SystemController> _systemController;

  //A vector holding the CI coefficients
  std::vector<Eigen::MatrixXd> _eigenvectors;

  //A vector holding the excitation energies
  Eigen::VectorXd _eigenvalues;




};

} /* namespace Serenity */

#endif /* POSTHF_LRSCF_ANALYSIS_EXCITATIONSPECTRUM_H_ */
