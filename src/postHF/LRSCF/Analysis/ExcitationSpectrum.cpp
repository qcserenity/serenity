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
#include "parameters/Constants.h"
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"

/* Include External Headers */
#include <fstream>
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void ExcitationSpectrum<SCFMode>::printSpectrum(Options::LR_METHOD method,
                                                const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                                const std::vector<Eigen::MatrixXd>& densityMatrices,
                                                const Eigen::VectorXd& eigenvalues, Eigen::Ref<Eigen::MatrixXd> results,
                                                std::string fileName) {
  // Get dipole integrals.
  Eigen::MatrixXcd dip_l = std::complex<double>(1, 0) * (*dipoles->getLengths());
  Eigen::MatrixXcd dip_v = std::complex<double>(0, 1) * (*dipoles->getVelocities());
  Eigen::MatrixXcd dip_m = std::complex<double>(0, 1) * (*dipoles->getMagnetics());

  unsigned nEigen = eigenvalues.size();

  /**
   * Transition strength matrices for each excitation
   * Operators:
   *   l = electric dipole (length)
   *   p = electric dipole (velocity)
   *   m = magnetic dipole
   */
  std::vector<Eigen::Matrix3d> S_ll(nEigen);
  std::vector<Eigen::Matrix3d> S_lv(nEigen);
  std::vector<Eigen::Matrix3d> S_vv(nEigen);
  std::vector<Eigen::Matrix3d> S_lm(nEigen);
  std::vector<Eigen::Matrix3d> S_vm(nEigen);
  std::vector<Eigen::Matrix3d> S_lm_mod(nEigen);

  // Matrices to store transition moments in.
  Eigen::MatrixXcd len_R, vel_R, mag_R;
  Eigen::MatrixXcd len_L, vel_L, mag_L;

  double factor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? std::sqrt(2) : 1.0;

  // Calculate left and right transition moments.
  if (method == Options::LR_METHOD::CC2 || method == Options::LR_METHOD::CISDINF) {
    // <0|O|n>, O in (mu, p, m).
    len_R = factor * dip_l.transpose() * densityMatrices[0];
    vel_R = factor * dip_v.transpose() * densityMatrices[0] * eigenvalues.cwiseInverse().asDiagonal();
    mag_R = factor / SPEEDOFLIGHT_AU * dip_m.transpose() * densityMatrices[0];

    // <n|O|0>, O in (mu, p, m).
    len_L = factor * densityMatrices[1].transpose() * dip_l;
    vel_L = factor * eigenvalues.cwiseInverse().asDiagonal() * densityMatrices[1].transpose() * dip_v;
    mag_L = factor / SPEEDOFLIGHT_AU * densityMatrices[1].transpose() * dip_m;
  }
  else {
    // <0|O|n>, O in (mu, p, m).
    len_R = factor * dip_l.adjoint() * densityMatrices[0];
    vel_R = factor * dip_v.adjoint() * densityMatrices[1];
    mag_R = factor / SPEEDOFLIGHT_AU * dip_m.adjoint() * densityMatrices[1];

    // Write transition strengths to disk.
    if (fileName != "") {
      Eigen::MatrixXd exspectrum = Eigen::MatrixXd::Zero(10, nEigen);
      exspectrum.row(0) = eigenvalues;
      exspectrum.middleRows(1, 3) = len_R.real();
      exspectrum.middleRows(4, 3) = vel_R.imag();
      exspectrum.middleRows(7, 3) = mag_R.imag();

      std::ofstream file(fileName + ".exspectrum.txt");
      file << std::scientific << std::setprecision(15) << exspectrum;
      file.close();
    }

    // Scale velocity transition moments.
    vel_R *= eigenvalues.cwiseInverse().asDiagonal();

    // <n|O|0>, O in (mu, p, m).
    len_L = len_R.adjoint();
    vel_L = vel_R.adjoint();
    mag_L = mag_R.adjoint();
  }

  // Calculate transition strengths.
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    S_ll[iEigen] = 0.5 * ((len_R.col(iEigen) * len_L.row(iEigen)) + (len_R.col(iEigen) * len_L.row(iEigen)).adjoint()).real();
    S_lv[iEigen] = 0.5 * ((len_R.col(iEigen) * vel_L.row(iEigen)) + (vel_R.col(iEigen) * len_L.row(iEigen)).adjoint()).imag();
    S_vv[iEigen] = 0.5 * ((vel_R.col(iEigen) * vel_L.row(iEigen)) + (vel_R.col(iEigen) * vel_L.row(iEigen)).adjoint()).real();
    S_lm[iEigen] = 0.5 * ((len_R.col(iEigen) * mag_L.row(iEigen)) + (mag_R.col(iEigen) * len_L.row(iEigen)).adjoint()).imag();
    S_vm[iEigen] = 0.5 * ((vel_R.col(iEigen) * mag_L.row(iEigen)) + (mag_R.col(iEigen) * vel_L.row(iEigen)).adjoint()).real();

    // Calculate modified electric dipole--magnetic dipole transition strength
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(S_lv[iEigen], Eigen::ComputeThinU | Eigen::ComputeThinV);
    S_lm_mod[iEigen] = svd.matrixU().transpose() * S_lm[iEigen] * svd.matrixV();
  }

  // Oscillator and rotator strengths.
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    results(iEigen, 0) = eigenvalues(iEigen);
    // Oscillator Strength (length)
    results(iEigen, 1) = 2.0 / 3.0 * eigenvalues(iEigen) * S_ll[iEigen].trace();
    // Oscillator Strength (velocity)
    results(iEigen, 2) = 2.0 / 3.0 * eigenvalues(iEigen) * S_vv[iEigen].trace();
    // Rotator Strength (length)
    results(iEigen, 3) = S_lm[iEigen].trace();
    // Rotator Strength (velocity)
    results(iEigen, 4) = S_vm[iEigen].trace();
    // Rotator Strength (modified length)
    results(iEigen, 5) = S_lm_mod[iEigen].trace();
  }

  // Print spectra.
  printf("                          Absorption Spectrum (dipole-length)                          \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state       energy      wavelength        fosc          Sxx        Syy        Szz     \n");
  printf("              (eV)          (nm)           (au)                    (au)                \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    printf(" %3i %15.5f %12.1f %15.6f %12.5f %10.5f %10.5f\n", iEigen + 1, eigenvalues(iEigen) * HARTREE_TO_EV,
           HARTREE_TO_NM / eigenvalues(iEigen), results(iEigen, 1), S_ll[iEigen](0, 0), S_ll[iEigen](1, 1),
           S_ll[iEigen](2, 2));
  }
  printf("---------------------------------------------------------------------------------------\n");
  printf("                         Absorption Spectrum (dipole-velocity)                         \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state       energy      wavelength        fosc          Sxx        Syy        Szz     \n");
  printf("              (eV)          (nm)           (au)                    (au)                \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    printf(" %3i %15.5f %12.1f %15.6f %12.5f %10.5f %10.5f\n", iEigen + 1, eigenvalues(iEigen) * HARTREE_TO_EV,
           HARTREE_TO_NM / eigenvalues(iEigen), results(iEigen, 2), S_vv[iEigen](0, 0), S_vv[iEigen](1, 1),
           S_vv[iEigen](2, 2));
  }
  printf("---------------------------------------------------------------------------------------\n");
  printf("                              CD Spectrum (dipole-length)                              \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state       energy      wavelength         R            Sxx        Syy        Szz     \n");
  printf("              (eV)          (nm)        (1e-40cgs)               (1e-40cgs)            \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    printf(" %3i %15.5f %12.1f %15.4f %12.5f %10.5f %10.5f\n", iEigen + 1, eigenvalues(iEigen) * HARTREE_TO_EV,
           HARTREE_TO_NM / eigenvalues(iEigen), AU_TO_CGS * results(iEigen, 3), AU_TO_CGS * S_lm[iEigen](0, 0),
           AU_TO_CGS * S_lm[iEigen](1, 1), AU_TO_CGS * S_lm[iEigen](2, 2));
  }
  printf("---------------------------------------------------------------------------------------\n");
  printf("                             CD Spectrum (dipole-velocity)                             \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state       energy      wavelength         R            Sxx        Syy        Szz     \n");
  printf("              (eV)          (nm)        (1e-40cgs)               (1e-40cgs)            \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    printf(" %3i %15.5f %12.1f %15.4f %12.5f %10.5f %10.5f\n", iEigen + 1, eigenvalues(iEigen) * HARTREE_TO_EV,
           HARTREE_TO_NM / eigenvalues(iEigen), AU_TO_CGS * results(iEigen, 4), AU_TO_CGS * S_vm[iEigen](0, 0),
           AU_TO_CGS * S_vm[iEigen](1, 1), AU_TO_CGS * S_vm[iEigen](2, 2));
  }
  printf("---------------------------------------------------------------------------------------\n");
  printf("                            CD Spectrum (mod. dipole-length)                           \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state       energy      wavelength         R            Sxx        Syy        Szz     \n");
  printf("              (eV)          (nm)        (1e-40cgs)               (1e-40cgs)            \n");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    printf(" %3i %15.5f %12.1f %15.4f %12.5f %10.5f %10.5f\n", iEigen + 1, eigenvalues(iEigen) * HARTREE_TO_EV,
           HARTREE_TO_NM / eigenvalues(iEigen), AU_TO_CGS * results(iEigen, 5), AU_TO_CGS * S_lm_mod[iEigen](0, 0),
           AU_TO_CGS * S_lm_mod[iEigen](1, 1), AU_TO_CGS * S_lm_mod[iEigen](2, 2));
  }
  printf("---------------------------------------------------------------------------------------\n");
}

template class ExcitationSpectrum<Options::SCF_MODES::RESTRICTED>;
template class ExcitationSpectrum<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
