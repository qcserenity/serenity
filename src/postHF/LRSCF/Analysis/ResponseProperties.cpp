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

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void ResponseProperties<SCFMode>::printProperties(const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                                                  const std::vector<Eigen::MatrixXd>& solutionvectors,
                                                  const std::vector<double>& frequencies, const double rotFactor,
                                                  const double damping, Options::GAUGE gauge) {
  // Get dipole integrals
  const auto& lengths = (*dipoles->getLengths());
  const auto& magnetics = (*dipoles->getMagnetics());

  // Set number of frequencies
  unsigned nFreqs = frequencies.size();

  // Set damped or not
  bool damped = (solutionvectors.size() == 4 ? true : false);

  // Initialize tensors
  Eigen::Matrix3d poly;
  Eigen::Matrix3d ordy;

  // Do gauge specific things and print some information about it
  if (gauge == Options::GAUGE::LENGTH) {
    printf("\n  Note: These response properties were calculated in LENGTH gauge and the\n");
    printf("        optical rotation tensor is therefore subject to origin dependence.\n\n");
    printf("        Tensors: alpha = -<<mu;mu>> on the left, G = -<<mu;m>> on the right.\n\n");
  }
  else if (gauge == Options::GAUGE::VELOCITY) {
    // The static limit is stored in the very last three columns of the
    // solution vectors as zero frequency was push backed in the beginning

    printf("\n  Note: These response properties were calculated in VELOCITY\n");
    printf("        gauge and are therefore gauge-origin independent.\n");
    printf("        The polarizability is given in mixed gauge.\n\n");
    printf("        Tensors: alpha = -<<mu;p>> on the left, G = -<<p;m>> on the right.\n\n");

    // Don't loop over the static limit itself (corresponds to the last frequency)
    --nFreqs;
  }
  else {
    assert(false);
  }
  printf("        The Specific Rotation is given in deg/[dm(g/mL)].\n\n");
  if (damped) {
    printf("        Running in damped mode. Damping factor (eV): %6.3f\n", damping * HARTREE_TO_EV);
  }

  printf("\n --------------------------------------------------------------------------------- \n\n");
  for (unsigned iFreq = 0; iFreq < nFreqs; ++iFreq) {
    printf("\n        Frequency: %16.6f au %13.3f eV %14.2f nm\n", frequencies[iFreq],
           frequencies[iFreq] * HARTREE_TO_EV, HARTREE_TO_NM / frequencies[iFreq]);

    // Define prefactors
    double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
    double gaugeFactor = (gauge == Options::GAUGE::LENGTH) ? 1.0 : 1.0 / frequencies[iFreq];
    double factor = scfFactor * gaugeFactor;

    /*
     * Obtain polarizability and G tensors as minus the respective response-function
     *
     *                         alpha = -<<mu;mu>>,
     *                         G = -<<mu;m>>.
     *
     * The outer minus (from below) comes from here. The solution vectors here are solutions to the TDDFT
     * response equations |X,Y>. We note that the response function in TDDFT can be written as
     *
     *                         <<A;V>> = -<A,A|[M - wS]^(-1)|-(Q,R)>
     *                                 = -<A,A|X,Y>
     *                                 = -<A|X+Y> if A is real
     *                              or = -<A|X-Y> if A is imaginary.
     *
     * The inner minus comes from here. Concerning the notation, |X,Y> implies a supervector containing X and Y,
     * and <X,Y| its Hermitian transpose.
     *
     */
    poly = -factor * (-lengths.transpose() * solutionvectors[0].middleCols(iFreq * 3, 3));
    ordy = -factor * (-magnetics.transpose() * solutionvectors[1].middleCols(iFreq * 3, 3));

    // We calculated the response function <<m,mu>>, but the actual definition of G can be found above
    // The result for the G tensor is thus transposed
    ordy.transposeInPlace();

    // beta = -1/(3w) TrG
    double beta0 = -1.0 / 3.0 / frequencies[iFreq] * ordy.trace();

    if (gauge == Options::GAUGE::VELOCITY) {
      // Replace tensor with velocity modified version
      ordy = factor * magnetics.transpose() * (solutionvectors[1].middleCols(iFreq * 3, 3) - solutionvectors[1].rightCols(3));
    }

    // Update optical rotation
    double beta1 = -1.0 / 3.0 / frequencies[iFreq] * ordy.trace();

    // Print out polarizability and optical rotation tensors and their isotropic components
    printf("\n\n\t\t\t         ### Real Part ###\n");
    printf("\n          x           y           z      |          x           y          z    \n\n");
    printf("  x %11.3e %11.3e %11.3e  |  x %11.3e %11.3e %11.3e\n", poly(0, 0), poly(0, 1), poly(0, 2), ordy(0, 0),
           ordy(0, 1), ordy(0, 2));
    printf("  y %11.3e %11.3e %11.3e  |  y %11.3e %11.3e %11.3e\n", poly(1, 0), poly(1, 1), poly(1, 2), ordy(1, 0),
           ordy(1, 1), ordy(1, 2));
    printf("  z %11.3e %11.3e %11.3e  |  z %11.3e %11.3e %11.3e\n\n", poly(2, 0), poly(2, 1), poly(2, 2), ordy(2, 0),
           ordy(2, 1), ordy(2, 2));

    printf("  Isotropic Polarizability: %53.6f\n\n", 1. / 3 * poly.trace());
    printf("  Isotropic Optical Rotation: %51.6f\n\n", beta1);
    if (gauge == Options::GAUGE::VELOCITY) {
      printf("  Isotropic Optical Rotation (unmod. velocity): %33.6f\n\n", beta0);
    }
    printf("  Specific Rotation: %60.6f\n\n", rotFactor * frequencies[iFreq] * frequencies[iFreq] * beta1);
    if (gauge == Options::GAUGE::VELOCITY) {
      printf("  Specific Rotation (unmod. velocity): %42.6f\n\n", rotFactor * frequencies[iFreq] * frequencies[iFreq] * beta0);
    }

    // Do the same for the imaginary part of the solution
    if (damped) {
      poly = -factor * (-lengths.transpose() * solutionvectors[2].middleCols(iFreq * 3, 3));
      ordy = -factor * (-magnetics.transpose() * solutionvectors[3].middleCols(iFreq * 3, 3));
      ordy.transposeInPlace();

      printf("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - - - \n\n");
      // Print out polarizability and optical rotation tensors and their isotropic components
      printf("\t\t\t       *** Imaginary Part ***\n");
      printf("\n          x           y           z      |          x           y          z    \n\n");
      printf("  x %11.3e %11.3e %11.3e  |  x %11.3e %11.3e %11.3e\n", poly(0, 0), poly(0, 1), poly(0, 2), ordy(0, 0),
             ordy(0, 1), ordy(0, 2));
      printf("  y %11.3e %11.3e %11.3e  |  y %11.3e %11.3e %11.3e\n", poly(1, 0), poly(1, 1), poly(1, 2), ordy(1, 0),
             ordy(1, 1), ordy(1, 2));
      printf("  z %11.3e %11.3e %11.3e  |  z %11.3e %11.3e %11.3e\n\n", poly(2, 0), poly(2, 1), poly(2, 2), ordy(2, 0),
             ordy(2, 1), ordy(2, 2));

      // See Jiemchooroj J. Chem. Phys. 127, 165104 (2007)
      printf("  Isotropic Absorption: %57.6f\n\n", 4. * PI * frequencies[iFreq] / 3. / 137.036 * poly.trace());
      printf("  Isotropic Circular Dichroism: %49.6f\n\n", -6.533 * frequencies[iFreq] * ordy.trace());
    }
    printf(" ---------------------------------------------------------------------------------\n\n");
  }
}

template class ResponseProperties<Options::SCF_MODES::RESTRICTED>;
template class ResponseProperties<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
