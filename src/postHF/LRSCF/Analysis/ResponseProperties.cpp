/**
 * @file ResponseProperties.cpp
 * @author Niklas Niemeyer
 *
 * @date May 12, 2018
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
/* Include Class Header*/
#include "postHF/LRSCF/Analysis/ResponseProperties.h"
/* Include Serenity Internal Headers */
#include "parameters/Constants.h"
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"
/* Include Std and External Headers */
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void ResponseProperties<SCFMode>::printProperties(
    bool isNotCC2, const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
    const std::vector<Eigen::MatrixXd>& perturbeddensities, const std::vector<double>& frequencies, const double rotFactor,
    const double damping, Options::GAUGE gauge, std::vector<Eigen::Matrix3d> Fdipdip, std::vector<Eigen::Matrix3d> Fdipmag,
    std::vector<std::tuple<double, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d>>& results) {
  // Get dipole integrals.
  Eigen::MatrixXcd len = std::complex<double>(1, 0) * (*dipoles->getLengths());
  Eigen::MatrixXcd vel = std::complex<double>(0, 1) * (*dipoles->getVelocities());
  Eigen::MatrixXcd mag = std::complex<double>(0, 1) * (*dipoles->getMagnetics());

  // Set number of frequencies.
  unsigned nFreqs = frequencies.size();

  // Set dimension of perturbed densities (ia for TDDFT, pq for CC2).
  unsigned nDimension = len.rows();

  // Set damped or not.
  bool damped = (perturbeddensities.size() == 4 && isNotCC2);

  // Set perturbed density matrices.
  Eigen::MatrixXcd xpy = Eigen::MatrixXcd::Zero(nDimension, nFreqs * 3);
  Eigen::MatrixXcd xmy = Eigen::MatrixXcd::Zero(nDimension, nFreqs * 3);

  xpy.real() = perturbeddensities[0];
  xmy.real() = perturbeddensities[1];
  if (damped || !isNotCC2) {
    xpy.imag() = perturbeddensities[2];
    xmy.imag() = perturbeddensities[3];
  }

  Eigen::Matrix3cd static_or;

  // do gauge specific things and print some information about it.
  if (gauge == Options::GAUGE::LENGTH) {
    printf("\n  Note: These response properties were calculated in LENGTH gauge.\n\n");
    if (isNotCC2) {
      printf("        The origin-independent optical rotation is computed from a\n");
      printf("        singular-value decomposition of the mixed-gauge polarizability.\n");
      printf("        See J. Chem. Phys. 153 (2020) 151101 for reference.\n\n");
    }
    printf("        Tensors: alpha = -<<mu;mu>> on the left.\n");
    printf("                     G = -<<mu;m>> on the right.\n\n");
  }
  else if (gauge == Options::GAUGE::VELOCITY) {
    // The static limit is stored in the very last three columns of the
    // solution vectors as zero frequency was push backed in the beginning
    printf("\n  Note: These response properties were calculated in VELOCITY gauge.\n\n");
    printf("        The origin-independent optical rotation is computed by\n");
    printf("        subtracting the static limit from the optical rotation tensor.\n");
    printf("        See Chem. Phys. Lett. 393 (2004) 319-326 for reference.\n\n");
    printf("        Tensors: alpha = -1/w<<p;p>> on the left.\n");
    printf("                     G = -1/w(<<p;m>>-<<p;m>>_0) on the right.\n\n");

    // Don't loop over the static limit itself (corresponds to the last frequency)
    --nFreqs;
  }
  printf("        The Specific Rotation is given in deg/[dm(g/mL)].\n\n");
  if (damped) {
    printf("        Running in damped mode. Damping factor (eV): %6.3f\n", damping * HARTREE_TO_EV);
  }

  printf("\n --------------------------------------------------------------------------------- \n");
  for (unsigned iFreq = 0; iFreq < nFreqs; ++iFreq) {
    printf("    Frequency: %16.6f au %13.3f eV %14.2f nm\n", frequencies[iFreq], frequencies[iFreq] * HARTREE_TO_EV,
           HARTREE_TO_NM / frequencies[iFreq]);

    // Define prefactors.
    double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
    double gaugeFactor = (gauge == Options::GAUGE::LENGTH) ? 1.0 : 1.0 / frequencies[iFreq];
    // CC2 response function follows symmetrized formulation, the TDDFT one does not.
    double methodFactor = isNotCC2 ? 1.0 : -0.5;
    double factor = scfFactor * gaugeFactor * methodFactor;

    auto& dip = (gauge == Options::GAUGE::LENGTH) ? len : vel;

    // push back one entry for results
    auto tuple = std::make_tuple(frequencies[iFreq], Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero(),
                                 Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero());
    results.push_back(tuple);

    // Initialize tensors.
    Eigen::Matrix3cd poly = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd ordy = Eigen::Matrix3cd::Zero();
    if (isNotCC2) {
      if (gauge == Options::GAUGE::LENGTH) {
        poly = factor * dip.adjoint() * xpy.middleCols(iFreq * 3, 3);

        ordy = factor * mag.adjoint() * xmy.middleCols(iFreq * 3, 3);
        ordy.transposeInPlace();
      }
      else {
        poly = factor * dip.adjoint() * (xmy.middleCols(iFreq * 3, 3) - xmy.rightCols(3));
        poly *= std::complex<double>(0, 1) / frequencies[iFreq];

        ordy = factor * mag.adjoint() * xmy.middleCols(iFreq * 3, 3);
        ordy.transposeInPlace();
        static_or = factor * mag.adjoint() * xmy.rightCols(3);
        static_or.transposeInPlace();
      }
    }
    else {
      if (gauge == Options::GAUGE::LENGTH) {
        poly += xpy.real().middleCols(iFreq * 3, 3).transpose() * dip;
        poly += xmy.real().middleCols(iFreq * 3, 3).transpose() * dip;
        poly += dip.transpose() * xpy.real().middleCols(iFreq * 3, 3);
        poly += dip.transpose() * xmy.real().middleCols(iFreq * 3, 3);
        poly += Fdipdip[iFreq];
        poly *= factor;

        ordy -= xpy.real().middleCols(iFreq * 3, 3).transpose() * mag;
        ordy += xmy.real().middleCols(iFreq * 3, 3).transpose() * mag;
        ordy += std::complex<double>(0, 1) * dip.transpose() * xpy.imag().middleCols(iFreq * 3, 3);
        ordy -= std::complex<double>(0, 1) * dip.transpose() * xmy.imag().middleCols(iFreq * 3, 3);
        ordy += std::complex<double>(0, 1) * Fdipmag[iFreq];
        ordy *= -1 * factor;
      }
      else {
        poly += (xpy.real().middleCols(iFreq * 3, 3) - xpy.real().rightCols(3)).transpose() * dip;
        poly += (xmy.real().middleCols(iFreq * 3, 3) - xmy.real().rightCols(3)).transpose() * dip;
        poly += dip.transpose() * (xpy.real().middleCols(iFreq * 3, 3) - xpy.real().rightCols(3));
        poly += dip.transpose() * (xmy.real().middleCols(iFreq * 3, 3) - xmy.real().rightCols(3));
        poly *= factor * std::complex<double>(0, 1);
        poly -= factor * (Fdipdip[iFreq] - Fdipdip[nFreqs]);
        poly /= frequencies[iFreq];

        ordy += factor * xpy.real().middleCols(iFreq * 3, 3).transpose() * mag;
        ordy += factor * xmy.real().middleCols(iFreq * 3, 3).transpose() * mag;
        ordy += factor * dip.transpose() * xpy.imag().middleCols(iFreq * 3, 3);
        ordy += factor * dip.transpose() * xmy.imag().middleCols(iFreq * 3, 3);
        ordy += std::complex<double>(0, 1) * factor * Fdipmag[iFreq];

        static_or = factor * xpy.real().rightCols(3).transpose() * mag;
        static_or += factor * xmy.real().rightCols(3).transpose() * mag;
        static_or += factor * dip.transpose() * xpy.imag().rightCols(3);
        static_or += factor * dip.transpose() * xmy.imag().rightCols(3);
        static_or += std::complex<double>(0, 1) * factor * Fdipmag[nFreqs];
      }
    }

    double alpha = +1.0 / 3.0 * poly.real().trace();
    double beta0 = -1.0 / 3.0 / frequencies[iFreq] * ordy.imag().trace();
    double abs = 4.0 * PI * frequencies[iFreq] / 3.0 / SPEEDOFLIGHT_AU * poly.imag().trace();
    double cd0 = 6.533 * frequencies[iFreq] * ordy.real().trace();

    // Obtain origin-independent optical rotation.
    if (gauge == Options::GAUGE::LENGTH && isNotCC2) {
      Eigen::Matrix3cd mixedpoly = vel.adjoint() * xmy.middleCols(iFreq * 3, 3);
      // <<mu;p>> is needed for this decomposition.
      mixedpoly.adjointInPlace();

      // Imaginary part of G (optical rotation).
      Eigen::JacobiSVD<Eigen::Matrix3d> svd_imag(mixedpoly.imag(), Eigen::ComputeFullU | Eigen::ComputeFullV);
      ordy.imag() = (svd_imag.matrixU().transpose() * ordy.imag() * svd_imag.matrixV()).eval();

      // Real part of G (circular dichroism).
      Eigen::JacobiSVD<Eigen::Matrix3d> svd_real(mixedpoly.real(), Eigen::ComputeFullU | Eigen::ComputeFullV);
      ordy.real() = (svd_real.matrixU().transpose() * ordy.real() * svd_real.matrixV()).eval();
    }
    else if (gauge == Options::GAUGE::VELOCITY) {
      ordy -= static_or;
    }

    // Update optical rotation and circular dichroism.
    double beta1 = -1.0 / 3.0 / frequencies[iFreq] * ordy.imag().trace();
    double cd1 = 6.533 * frequencies[iFreq] * ordy.real().trace();

    // Store in results.
    std::get<1>(results[iFreq]) = poly.real();
    std::get<2>(results[iFreq]) = ordy.imag();
    std::get<3>(results[iFreq]) = poly.imag();
    std::get<4>(results[iFreq]) = ordy.real();

    printf("\n\t\t\t       ### Real Response ###\n");
    printf("\n          x           y           z      |          x           y          z    \n\n");
    printf("  x %11.6f %11.6f %11.6f  |  x %11.6f %11.6f %11.6f\n", poly(0, 0).real(), poly(0, 1).real(),
           poly(0, 2).real(), ordy(0, 0).imag(), ordy(0, 1).imag(), ordy(0, 2).imag());
    printf("  y %11.6f %11.6f %11.6f  |  y %11.6f %11.6f %11.6f\n", poly(1, 0).real(), poly(1, 1).real(),
           poly(1, 2).real(), ordy(1, 0).imag(), ordy(1, 1).imag(), ordy(1, 2).imag());
    printf("  z %11.6f %11.6f %11.6f  |  z %11.6f %11.6f %11.6f\n\n", poly(2, 0).real(), poly(2, 1).real(),
           poly(2, 2).real(), ordy(2, 0).imag(), ordy(2, 1).imag(), ordy(2, 2).imag());

    printf("  Isotropic Polarizability: %53.10f\n", alpha);
    printf("  Isotropic Optical Rotation: %51.10f\n", beta0);
    printf("  Specific Rotation: %60.7f\n", rotFactor * frequencies[iFreq] * frequencies[iFreq] * beta0);
    if (isNotCC2 || gauge == Options::GAUGE::VELOCITY) {
      printf("  Isotropic Optical Rotation (origin-independent): %30.10f\n", beta1);
      printf("  Specific Rotation (origin-independent): %39.10f\n", rotFactor * frequencies[iFreq] * frequencies[iFreq] * beta1);
    }

    if (damped) {
      printf("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - - - \n\n");
      printf("\t\t\t     *** Imaginary Response ***\n");
      printf("\n          x           y           z      |          x           y          z    \n\n");
      printf("  x %11.6f %11.6f %11.6f  |  x %11.6f %11.6f %11.6f\n", poly(0, 0).imag(), poly(0, 1).imag(),
             poly(0, 2).imag(), ordy(0, 0).real(), ordy(0, 1).real(), ordy(0, 2).real());
      printf("  y %11.6f %11.6f %11.6f  |  y %11.6f %11.6f %11.6f\n", poly(1, 0).imag(), poly(1, 1).imag(),
             poly(1, 2).imag(), ordy(1, 0).real(), ordy(1, 1).real(), ordy(1, 2).real());
      printf("  z %11.6f %11.6f %11.6f  |  z %11.6f %11.6f %11.6f\n\n", poly(2, 0).imag(), poly(2, 1).imag(),
             poly(2, 2).imag(), ordy(2, 0).real(), ordy(2, 1).real(), ordy(2, 2).real());

      // See Jiemchooroj J. Chem. Phys. 127, 165104 (2007).
      printf("  Isotropic Absorption: %57.10f\n", abs);
      printf("  Isotropic Circular Dichroism: %49.10f\n", cd0);
      printf("  Isotropic Circular Dichroism (origin-independent): %28.10f\n", cd1);
    }
    printf(" ---------------------------------------------------------------------------------\n\n");
  }
}

template class ResponseProperties<Options::SCF_MODES::RESTRICTED>;
template class ResponseProperties<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
