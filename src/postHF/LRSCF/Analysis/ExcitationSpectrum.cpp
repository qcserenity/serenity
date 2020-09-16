/**
 * @file ExcitationSpectrum.cpp
 * @author: Michael Boeckers, Niklas Niemeyer
 *
 * @date Dec. 17, 2018
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
/* Include Class Header*/
#include "postHF/LRSCF/Analysis/ExcitationSpectrum.h"

/* Include Serenity Internal Headers */

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void ExcitationSpectrum<SCFMode>::printSpectrum(const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                                const std::vector<Eigen::MatrixXd>& eigenvectors,
                                                const Eigen::VectorXd& eigenvalues) {
  // Get dipole integrals
  auto dip_l = dipoles->getLengths();
  auto dip_v = dipoles->getVelocities();
  auto dip_m = dipoles->getMagnetics();

  // Calculate (X+Y) and (X-Y) excitation vector
  Eigen::MatrixXd xpy = eigenvectors[0];
  Eigen::MatrixXd xmy = eigenvectors[0];
  if (eigenvectors.size() == 2) {
    xpy += eigenvectors[1];
    xmy -= eigenvectors[1];
  }

  const double factor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? std::sqrt(2) : 1.0;

  // Transition dipole moments (<0|O|n>) where O in mu, p , m
  Eigen::MatrixXd tdm_l = factor * xpy.transpose() * (*dip_l);
  Eigen::MatrixXd tdm_v = factor * xmy.transpose() * (*dip_v);
  Eigen::MatrixXd tdm_m = factor / SPEEDOFLIGHT_AU * xmy.transpose() * (*dip_m);
  for (unsigned int iState = 0; iState < tdm_v.rows(); iState++)
    tdm_v.row(iState) /= eigenvalues(iState);

  // Oscillator and Rotatory Strenghts
  Eigen::VectorXd os_l(eigenvalues.rows());
  Eigen::VectorXd os_v(eigenvalues.rows());
  Eigen::VectorXd rs_l(eigenvalues.rows());
  Eigen::VectorXd rs_v(eigenvalues.rows());
  for (unsigned int iState = 0; iState < tdm_l.rows(); ++iState) {
    os_l(iState) = 2.0 / 3.0 * eigenvalues(iState) * tdm_l.row(iState).array().square().sum();
    os_v(iState) = 2.0 / 3.0 * eigenvalues(iState) * tdm_v.row(iState).array().square().sum();
    // R = Im<0|mu|n><n|m|0> = Im<0|mu|n>(<n|m|0>)^* = Im<0|mu|n>(-<0|m|n>)
    rs_l(iState) = tdm_l.row(iState).dot(-tdm_m.row(iState));
    rs_v(iState) = tdm_v.row(iState).dot(-tdm_m.row(iState));
  }

  // Print spectra
  printf("---------------------------------------------------------------------------------------\n");
  printf("                          Absorption Spectrum (dipole-length)                          \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state        energy        wavelength    fosc      mu**2     mu_x     mu_y     mu_z   \n");
  printf("          (eV)      (cm-1)     (nm)       (au)     (au**2)    (au)     (au)     (au)   \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned int iState = 0; iState < eigenvalues.rows(); ++iState) {
    printf(" %3i %10.4f %10.1f %8.1f %12.6f %9.4f %8.4f %8.4f %8.4f\n", iState + 1, eigenvalues(iState) * HARTREE_TO_EV,
           eigenvalues(iState) * HARTREE_TO_OOCM, HARTREE_TO_NM / eigenvalues(iState), os_l(iState),
           tdm_l.row(iState).array().square().sum(), tdm_l(iState, 0), tdm_l(iState, 1), tdm_l(iState, 2));
  }
  printf("---------------------------------------------------------------------------------------\n");
  printf("                         Absorption Spectrum (dipole-velocity)                         \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state        energy        wavelength    fosc      mu**2     mu_x     mu_y     mu_z   \n");
  printf("          (eV)      (cm-1)     (nm)       (au)     (au**2)    (au)     (au)     (au)   \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned int iState = 0; iState < eigenvalues.rows(); ++iState) {
    printf(" %3i %10.4f %10.1f %8.1f %12.6f %9.4f %8.4f %8.4f %8.4f\n", iState + 1, eigenvalues(iState) * HARTREE_TO_EV,
           eigenvalues(iState) * HARTREE_TO_OOCM, HARTREE_TO_NM / eigenvalues(iState), os_v(iState),
           tdm_v.row(iState).array().square().sum(), tdm_v(iState, 0), tdm_v(iState, 1), tdm_v(iState, 2));
  }
  printf("---------------------------------------------------------------------------------------\n");
  printf("                              CD Spectrum (dipole-length)                              \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state        energy        wavelength     R        mu**2     mu_x     mu_y     mu_z   \n");
  printf("          (eV)      (cm-1)     (nm)    (1e-40cgs)  (au**2)    (au)     (au)     (au)   \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned int iState = 0; iState < eigenvalues.rows(); ++iState) {
    printf(" %3i %10.4f %10.1f %8.1f %12.4f %9.4f %8.4f %8.4f %8.4f\n", iState + 1, eigenvalues(iState) * HARTREE_TO_EV,
           eigenvalues(iState) * HARTREE_TO_OOCM, HARTREE_TO_NM / eigenvalues(iState), AU_TO_CGS * rs_l(iState),
           tdm_m.row(iState).array().square().sum(), tdm_m(iState, 0), tdm_m(iState, 1), tdm_m(iState, 2));
  }
  printf("---------------------------------------------------------------------------------------\n");
  printf("                             CD Spectrum (dipole-velocity)                             \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state        energy        wavelength     R        mu**2     mu_x     mu_y     mu_z   \n");
  printf("          (eV)      (cm-1)     (nm)    (1e-40cgs)  (au**2)    (au)     (au)     (au)   \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned int iState = 0; iState < eigenvalues.rows(); ++iState) {
    printf(" %3i %10.4f %10.1f %8.1f %12.4f %9.4f %8.4f %8.4f %8.4f\n", iState + 1, eigenvalues(iState) * HARTREE_TO_EV,
           eigenvalues(iState) * HARTREE_TO_OOCM, HARTREE_TO_NM / eigenvalues(iState), AU_TO_CGS * rs_v(iState),
           tdm_m.row(iState).array().square().sum(), tdm_m(iState, 0), tdm_m(iState, 1), tdm_m(iState, 2));
  }
  printf("---------------------------------------------------------------------------------------\n\n\n");
}

template class ExcitationSpectrum<Options::SCF_MODES::RESTRICTED>;
template class ExcitationSpectrum<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
