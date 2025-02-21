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
/* Include Std and External Headers */
#include <fstream>
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void ExcitationSpectrum<SCFMode>::printTransitionMoments(Options::LR_METHOD method,
                                                         const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                                         const std::vector<Eigen::MatrixXd>& densityMatrices,
                                                         const Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& results,
                                                         std::string fileName) {
  // Get dipole integrals.
  Eigen::MatrixXcd dip_l = std::complex<double>(1, 0) * (*dipoles->getLengths());
  Eigen::MatrixXcd dip_v = std::complex<double>(0, 1) * (*dipoles->getVelocities());
  Eigen::MatrixXcd dip_m = std::complex<double>(0, 1) * (*dipoles->getMagnetics());

  unsigned nEigen = densityMatrices[0].cols();

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
    // AR: changed this from a Matrix3d to a MatrixXd since Eigen debug says that thin U and V are only available for
    // dynamic-sized matrices
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S_lv[iEigen], Eigen::ComputeThinU | Eigen::ComputeThinV);
    S_lm_mod[iEigen] = svd.matrixU().transpose() * S_lm[iEigen] * svd.matrixV();
  }

  // Oscillator and rotator strengths.
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
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

  if (nEigen > 0) {
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
}

template<Options::SCF_MODES SCFMode>
void ExcitationSpectrum<SCFMode>::printStateMoments(Options::LR_METHOD, const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                                    const std::vector<Eigen::MatrixXd>& densityMatrices,
                                                    const Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& results,
                                                    Eigen::Vector3d nuclearPart) {
  // Get dipole integrals (in MO basis, dip_l has nocc * nvirt rows and 3 columns)
  Eigen::MatrixXcd dip_l = std::complex<double>(1, 0) * (*dipoles->getLengths());
  // Note: The first column in densityMatrices corresponds to the ground state density.
  unsigned nEigenExcitedStateDensities = densityMatrices[2].cols() - 1;
  Eigen::MatrixXd dipoleMoments = (dip_l.transpose() * densityMatrices[2]).real();
  dipoleMoments.colwise() += nuclearPart;
  // dipoleMoments also contains a column for the ground state, don't want to put that in the results matrix.
  results.conservativeResize(Eigen::NoChange, 7);
  results.col(6) = dipoleMoments.colwise().norm().tail(nEigenExcitedStateDensities);
  printf("                          State Dipole Moments (dipole-length)                         \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state       energy      wavelength        |mu|           x          y          z      \n");
  printf("              (eV)          (nm)           (au)                    (au)                \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" Nuclear part      :  %27.6f %12.6f %10.6f %10.6f\n", nuclearPart.norm(), nuclearPart(0), nuclearPart(1),
         nuclearPart(2));
  printf(" Ground-state part :  %27.6f %12.6f %10.6f %10.6f\n", dipoleMoments.col(0).norm(), dipoleMoments(0, 0),
         dipoleMoments(1, 0), dipoleMoments(2, 0));
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned iEigen = 0; iEigen < nEigenExcitedStateDensities; ++iEigen) {
    printf(" %3i %15.5f %12.1f %15.6f %12.6f %10.6f %10.6f\n", iEigen + 1, eigenvalues(iEigen) * HARTREE_TO_EV,
           HARTREE_TO_NM / eigenvalues(iEigen), dipoleMoments.col(iEigen + 1).norm(), dipoleMoments(0, iEigen + 1),
           dipoleMoments(1, iEigen + 1), dipoleMoments(2, iEigen + 1));
  }
  printf("---------------------------------------------------------------------------------------\n");
}

template<Options::SCF_MODES SCFMode>
void ExcitationSpectrum<SCFMode>::printStateMoments(Options::LR_METHOD, const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                                    const std::vector<MatrixInBasis<SCFMode>>& unrelaxedDensities,
                                                    const std::vector<MatrixInBasis<SCFMode>>& relaxedDensities,
                                                    const Eigen::VectorXd& eigenvalues, Eigen::Vector3d nuclearPart,
                                                    SpinPolarizedData<SCFMode, unsigned int> nOcc) {
  // Get dipole integrals (in MO basis, dip_l has nbasis * nbasis rows and 3 columns).
  dipoles->setFullSpace(true);
  Eigen::MatrixXcd dip_l = std::complex<double>(1, 0) * (*dipoles->getLengths());
  // reset dipole integrals size
  dipoles->setFullSpace(false);
  Eigen::MatrixXd dipoleMoments, relaxedDipoleMoments;
  unsigned nBasis = unrelaxedDensities[0].rows();
  Eigen::VectorXd groundStateVector((1 + (int)SCFMode) * nBasis * nBasis);
  unsigned nStart = 0;
  for_spin(nOcc) {
    Eigen::MatrixXd id(Eigen::MatrixXd::Zero(nBasis, nBasis));
    id.topLeftCorner(nOcc_spin, nOcc_spin) =
        ((SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0) * Eigen::MatrixXd::Identity(nOcc_spin, nOcc_spin);
    groundStateVector.segment(nStart, nBasis * nBasis) = Eigen::Map<const Eigen::VectorXd>(id.data(), nBasis * nBasis);
    nStart += nBasis * nBasis;
  };
  Eigen::VectorXd groundStateDipoleMoment = (dip_l.transpose() * groundStateVector).real();
  groundStateDipoleMoment += nuclearPart;
  dipoleMoments.resize(3, unrelaxedDensities.size());
  for (unsigned iEigen = 0; iEigen < unrelaxedDensities.size(); ++iEigen) {
    Eigen::VectorXd unrolledDens(groundStateVector.size());
    const MatrixInBasis<SCFMode>& temp = unrelaxedDensities[iEigen];
    nStart = 0;
    for_spin(temp) {
      unrolledDens.segment(nStart, nBasis * nBasis) = Eigen::Map<const Eigen::VectorXd>(temp_spin.data(), nBasis * nBasis);
      nStart += nBasis * nBasis;
    };
    dipoleMoments.col(iEigen) = (dip_l.transpose() * unrolledDens).real();
  }
  dipoleMoments.colwise() += groundStateDipoleMoment;
  relaxedDipoleMoments.resize(3, relaxedDensities.size());
  for (unsigned iEigen = 0; iEigen < relaxedDensities.size(); ++iEigen) {
    Eigen::VectorXd unrolledDens(groundStateVector.size());
    const MatrixInBasis<SCFMode>& temp = relaxedDensities[iEigen];
    nStart = 0;
    for_spin(temp) {
      unrolledDens.segment(nStart, nBasis * nBasis) = Eigen::Map<const Eigen::VectorXd>(temp_spin.data(), nBasis * nBasis);
      nStart += nBasis * nBasis;
    };
    relaxedDipoleMoments.col(iEigen) = (dip_l.transpose() * unrolledDens).real();
  }
  relaxedDipoleMoments.colwise() += groundStateDipoleMoment;
  printf("                         TDDFT State Dipole Moments (dipole-length)                    \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" state           energy     wavelength     |mu|           x          y          z      \n");
  printf("                  (eV)        (nm)         (au)                    (au)                \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" Nuclear part      :  %27.6f %12.6f %10.6f %10.6f\n", nuclearPart.norm(), nuclearPart(0), nuclearPart(1),
         nuclearPart(2));
  printf(" Ground-state part :  %27.6f %12.6f %10.6f %10.6f\n", groundStateDipoleMoment.col(0).norm(),
         groundStateDipoleMoment(0, 0), groundStateDipoleMoment(1, 0), groundStateDipoleMoment(2, 0));
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned iEigen = 0; iEigen < unrelaxedDensities.size(); ++iEigen) {
    printf(" %3i unrlx: %11.5f %10.1f %15.6f %12.6f %10.6f %10.6f\n", iEigen + 1, eigenvalues(iEigen) * HARTREE_TO_EV,
           HARTREE_TO_NM / eigenvalues(iEigen), dipoleMoments.col(iEigen).norm(), dipoleMoments(0, iEigen),
           dipoleMoments(1, iEigen), dipoleMoments(2, iEigen));
    if (iEigen < relaxedDensities.size())
      printf(" %3i rlx  : %11.5f %10.1f %15.6f %12.6f %10.6f %10.6f\n", iEigen + 1, eigenvalues(iEigen) * HARTREE_TO_EV,
             HARTREE_TO_NM / eigenvalues(iEigen), relaxedDipoleMoments.col(iEigen).norm(),
             relaxedDipoleMoments(0, iEigen), relaxedDipoleMoments(1, iEigen), relaxedDipoleMoments(2, iEigen));
  }
  printf("---------------------------------------------------------------------------------------\n");
}

template class ExcitationSpectrum<Options::SCF_MODES::RESTRICTED>;
template class ExcitationSpectrum<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
