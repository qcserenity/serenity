/**
 * @file ResponseProperties.cpp
 * @author: Niklas Niemeyer
 *
 * @date May 12, 2019
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

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void ResponseProperties<SCFMode>::printProperties(
    const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles, const std::vector<Eigen::MatrixXd>& solutionvectors,
    const std::vector<double>& frequencies, const double rotFactor, const double damping, Options::GAUGE gauge,
    std::vector<std::tuple<double, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d>>& results) {
  // Get dipole integrals.
  Eigen::MatrixXcd dip_l = std::complex<double>(1, 0) * (*dipoles->getLengths());
  Eigen::MatrixXcd dip_v = std::complex<double>(0, 1) * (*dipoles->getVelocities());
  Eigen::MatrixXcd dip_m = std::complex<double>(0, 1) * (*dipoles->getMagnetics());

  // Set number of frequencies.
  unsigned nFreqs = frequencies.size();

  // Set dimension of response problem.
  unsigned nDimension = dip_l.rows();

  // Set damped or not.
  bool damped = (solutionvectors.size() == 4);

  // Set perturbed density matrices.
  Eigen::MatrixXcd xpy = Eigen::MatrixXcd::Zero(nDimension, nFreqs * 3);
  Eigen::MatrixXcd xmy = Eigen::MatrixXcd::Zero(nDimension, nFreqs * 3);
  xpy.real() = solutionvectors[0];
  xmy.real() = solutionvectors[1];
  if (damped) {
    xpy.imag() = solutionvectors[2];
    xmy.imag() = solutionvectors[3];
  }

  // Initialize tensors.
  Eigen::Matrix3cd poly;
  Eigen::Matrix3cd ordy;

  // do gauge specific things and print some information about it
  if (gauge == Options::GAUGE::LENGTH) {
    printf("\n  Note: These response properties were calculated in LENGTH gauge.\n\n");
    printf("        The origin-independent optical rotation is computed from a\n");
    printf("        singular-value decomposition of the mixed-gauge polarizability.\n");
    printf("        See J. Chem. Phys. 153 (2020) 151101 for reference.\n\n");
    printf("        Tensors: alpha = -<<mu;mu>> on the left, G = -<<mu;m>> on the right.\n\n");
  }
  else if (gauge == Options::GAUGE::VELOCITY) {
    // The static limit is stored in the very last three columns of the
    // solution vectors as zero frequency was push backed in the beginning
    printf("\n  Note: These response properties were calculated in VELOCITY gauge.\n\n");
    printf("        The origin-independent optical rotation is computed by\n");
    printf("        subtracting the static limit from the optical rotation tensor.\n");
    printf("        See Chem. Phys. Lett. 393 (2004) 319–326 for reference.\n\n");
    printf("        Tensors: alpha = -<<mu;p>> on the left, G = -<<p;m>> on the right.\n\n");

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
    double factor = scfFactor * gaugeFactor;

    // push back one entry for results
    std::tuple<double, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d> tuple =
        std::make_tuple(frequencies[iFreq], Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero(),
                        Eigen::Matrix3d::Zero());
    results.push_back(tuple);

    /*
     * Obtain polarizability and G tensors as minus the respective response-function
     *
     *                         alpha = -<<mu;mu>>,
     *                             G = -<<mu;m>>.
     *
     * the outer minus (from below) comes from here. The solution vectors here are solutions to the TDDFT
     * response equations |X,Y>. We note that the response function in TDDFT can be written as
     *
     *                         <<A;V>> = -<A,A|[M - wS]^(-1)|-(Q,R)>
     *                                 = -<A,A|X,Y>
     *                                 = -<A|X+Y> if A is real
     *                              or = -<A|X-Y> if A is imaginary.
     *
     * the inner minus comes from here. Concerning the notation, |X,Y> implies a supervector containing X and Y,
     * and <X,Y| its Hermitian transpose.
     *
     */
    poly = factor * dip_l.adjoint() * xpy.middleCols(iFreq * 3, 3);
    ordy = factor * dip_m.adjoint() * xmy.middleCols(iFreq * 3, 3);

    // We have calculated the <<m;mu>> response function, however, <<mu;m>> is requested.
    ordy.transposeInPlace();

    double alpha = +1.0 / 3.0 * poly.real().trace();
    double beta0 = -1.0 / 3.0 / frequencies[iFreq] * ordy.imag().trace();
    double abs = 4.0 * PI * frequencies[iFreq] / 3.0 / SPEEDOFLIGHT_AU * poly.imag().trace();
    // TODO check sign here.
    double cd0 = 6.533 * frequencies[iFreq] * ordy.real().trace();

    // Obtain origin-independent optical rotation.
    if (gauge == Options::GAUGE::LENGTH) {
      Eigen::Matrix3cd mixedpoly = factor / frequencies[iFreq] * dip_v.adjoint() * xmy.middleCols(iFreq * 3, 3);
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
      ordy -= factor * dip_m.adjoint() * xmy.rightCols(3);
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
    printf("  x %11.3e %11.3e %11.3e  |  x %11.3e %11.3e %11.3e\n", poly(0, 0).real(), poly(0, 1).real(),
           poly(0, 2).real(), ordy(0, 0).imag(), ordy(0, 1).imag(), ordy(0, 2).imag());
    printf("  y %11.3e %11.3e %11.3e  |  y %11.3e %11.3e %11.3e\n", poly(1, 0).real(), poly(1, 1).real(),
           poly(1, 2).real(), ordy(1, 0).imag(), ordy(1, 1).imag(), ordy(1, 2).imag());
    printf("  z %11.3e %11.3e %11.3e  |  z %11.3e %11.3e %11.3e\n\n", poly(2, 0).real(), poly(2, 1).real(),
           poly(2, 2).real(), ordy(2, 0).imag(), ordy(2, 1).imag(), ordy(2, 2).imag());

    printf("  Isotropic Polarizability: %53.7f\n", alpha);
    printf("  Isotropic Optical Rotation: %51.7f\n", beta0);
    printf("  Specific Rotation: %60.7f\n", rotFactor * frequencies[iFreq] * frequencies[iFreq] * beta0);
    printf("  Isotropic Optical Rotation (origin-independent): %30.7f\n", beta1);
    printf("  Specific Rotation (origin-independent): %39.7f\n", rotFactor * frequencies[iFreq] * frequencies[iFreq] * beta1);

    if (damped) {
      printf("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - - - \n\n");
      printf("\t\t\t     *** Imaginary Response ***\n");
      printf("\n          x           y           z      |          x           y          z    \n\n");
      printf("  x %11.3e %11.3e %11.3e  |  x %11.3e %11.3e %11.3e\n", poly(0, 0).imag(), poly(0, 1).imag(),
             poly(0, 2).imag(), ordy(0, 0).real(), ordy(0, 1).real(), ordy(0, 2).real());
      printf("  y %11.3e %11.3e %11.3e  |  y %11.3e %11.3e %11.3e\n", poly(1, 0).imag(), poly(1, 1).imag(),
             poly(1, 2).imag(), ordy(1, 0).real(), ordy(1, 1).real(), ordy(1, 2).real());
      printf("  z %11.3e %11.3e %11.3e  |  z %11.3e %11.3e %11.3e\n\n", poly(2, 0).imag(), poly(2, 1).imag(),
             poly(2, 2).imag(), ordy(2, 0).real(), ordy(2, 1).real(), ordy(2, 2).real());

      // See Jiemchooroj J. Chem. Phys. 127, 165104 (2007).
      printf("  Isotropic Absorption: %57.7f\n", abs);
      printf("  Isotropic Circular Dichroism: %49.7f\n", cd0);
      printf("  Isotropic Circular Dichroism (origin-independent): %28.7f\n", cd1);
    }
    printf(" ---------------------------------------------------------------------------------\n\n");
  }
}

template class ResponseProperties<Options::SCF_MODES::RESTRICTED>;
template class ResponseProperties<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
