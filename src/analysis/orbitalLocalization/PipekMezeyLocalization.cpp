/**
 * @file   PipekMezeyLocalization.cpp
 *
 * @date   Apr 22, 2014
 * @author Thomas Dresselhaus
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
#include "analysis/orbitalLocalization/PipekMezeyLocalization.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/Matrix.h"
#include "math/linearAlgebra/JacobiRotation.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <cassert>
#include <cmath>
#include <limits>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
PipekMezeyLocalization<SCFMode>::PipekMezeyLocalization(std::shared_ptr<SystemController> systemController)
  : _systemController(systemController),
    // TODO this should be switchable via the input
    _convThreshold(1e-7) {
  assert(_systemController);
}

/*
 * Implementation according to
 * [1] J. W. Boughton, P. Pulay; J. Comp. Chem. 14 (1993), 736 - 740.
 */

template<Options::SCF_MODES SCFMode>
void PipekMezeyLocalization<SCFMode>::localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                                       SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) {
  CoefficientMatrix<SCFMode> coefficients = orbitals.getCoefficients();
  const auto& nOccOrbs = _systemController->getNOccupiedOrbitals<SCFMode>();

  /*
   * Do alpha and beta localization consecutively
   */
  for_spin(coefficients, nOccOrbs, orbitalRange) {
    /*
     * Get in / prepare some data for more convenient handling.
     */
    const unsigned int nAtoms = _systemController->getNAtoms();
    double oldLocalizationMeasure = std::numeric_limits<double>::lowest();
    // TODO this should still be specified by the incoming matrix.
    const auto basisController = _systemController->getAtomCenteredBasisController();
    const unsigned int nBasisFunctions = basisController->getNBasisFunctions();
    const std::vector<std::pair<unsigned int, unsigned int>>& basisIndices = basisController->getBasisIndices();
    const auto oneIntController = _systemController->getOneElectronIntegralController();
    assert(oneIntController->getBasisController() == basisController);
    const auto& overlaps = oneIntController->getOverlapIntegrals();

    unsigned int counter = 0;
    while (true) {
      /*
       * The gross Mulliken population of orbital i on atom a is Q(a,i).
       */
      Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(nAtoms, nOccOrbs_spin);

      /*
       * Calculate all Q(a,i) -> eq. [1].(5)
       */
      for (unsigned int iOrb = 0; iOrb < orbitalRange_spin.size(); ++iOrb) {
        unsigned int i = orbitalRange_spin[iOrb];
        for (unsigned int a = 0; a < nAtoms; ++a) {
          double& Q_ai = Q(a, i);
          // Basis functions on atom a
          for (unsigned int mu = basisIndices[a].first; mu < basisIndices[a].second; ++mu) {
            for (unsigned int nu = 0; nu < nBasisFunctions; ++nu) {
              Q_ai += coefficients_spin(mu, i) * coefficients_spin(nu, i) * overlaps(mu, nu);
            }
          }
        }
      }
      /*
       * Calculate all rotation angles gamma(s,t) and keep the largest rotation
       */
      Eigen::MatrixXd gammaST = Eigen::MatrixXd::Zero(nOccOrbs_spin, nOccOrbs_spin);

#pragma omp parallel
      {
#pragma omp for schedule(dynamic)
        for (unsigned int sOrb = 0; sOrb < orbitalRange_spin.size(); ++sOrb) {
          unsigned int s = orbitalRange_spin[sOrb];
          // Loop until t<s because a rotation of an orbital with itself won't make sense.
          for (unsigned int tOrb = 0; tOrb < sOrb; ++tOrb) {
            unsigned int t = orbitalRange_spin[tOrb];
            // Calculate the Q(a,s,t).
            Eigen::VectorXd Q_a_st = Eigen::VectorXd::Zero(nAtoms);
            for (unsigned int a = 0; a < nAtoms; ++a) {
              // Basis functions on atom a
              for (unsigned int mu = basisIndices[a].first; mu < basisIndices[a].second; ++mu) {
                double coefficientsMuS = coefficients_spin(mu, s);
                double coefficientsMuT = coefficients_spin(mu, t);
                for (unsigned int nu = 0; nu < nBasisFunctions; ++nu) {
                  Q_a_st[a] += (coefficientsMuS * coefficients_spin(nu, t) + coefficients_spin(nu, s) * coefficientsMuT) *
                               overlaps(nu, mu);
                }
              }
              Q_a_st[a] *= 0.5;
            }
            // Calculate B(s,t)
            double B_st = 0.0;
            for (unsigned int a = 0; a < nAtoms; ++a) {
              B_st += Q_a_st[a] * (Q(a, s) - Q(a, t));
            }
            // Calculate A(s,t)
            double A_st = 0;
            for (unsigned int a = 0; a < nAtoms; ++a) {
              A_st += Q_a_st[a] * Q_a_st[a] - 1.0 / 4.0 * (Q(a, s) - Q(a, t)) * (Q(a, s) - Q(a, t));
            }
            // Calculate gamma(s,t)
            double gamma_st = B_st / fabs(B_st) * 1.0 / 4.0 * acos(-A_st / sqrt(A_st * A_st + B_st * B_st));
            gammaST(s, t) = gamma_st;
            gammaST(t, s) = gamma_st;
          }
        }

      } /* pragma omp parallel */

      // Identify the largest (absolute) rotation angle
      unsigned int sOfMaxGamma = 0;
      unsigned int tOfMaxGamma = 0;
      gammaST.array().abs().maxCoeff(&sOfMaxGamma, &tOfMaxGamma);
      double maxGamma = gammaST(sOfMaxGamma, tOfMaxGamma);

      // If largest rotation is small enough: converged
      if (fabs(maxGamma) < _convThreshold) {
        break;
      }
      // Rotate orbitals with largest angle
      JacobiRotation::rotate(coefficients_spin.col(sOfMaxGamma), coefficients_spin.col(tOfMaxGamma), maxGamma);

      // Rotate other pairs of orbitals which were untouched so far
      std::vector<bool> wasRotated(nOccOrbs_spin, false);
      wasRotated[sOfMaxGamma] = true;
      wasRotated[tOfMaxGamma] = true;
      for (unsigned int sOrb = 0; sOrb < orbitalRange_spin.size(); ++sOrb) {
        unsigned int s = orbitalRange_spin[sOrb];
        if (wasRotated[s])
          continue;
        // Identify the largest (absolute) rotation angle
        unsigned int tOfMaxGammaThisS = 0;
        gammaST.row(s).array().abs().maxCoeff(&tOfMaxGammaThisS);
        double maxGammaThisS = gammaST(s, tOfMaxGammaThisS);
        // Do rotation
        JacobiRotation::rotate(coefficients_spin.col(s), coefficients_spin.col(tOfMaxGammaThisS), maxGammaThisS);
        // Mark as rotated
        wasRotated[s] = true;
        wasRotated[tOfMaxGammaThisS] = true;
      }
      // Check the optimization measure
      double P = Q.cwiseProduct(Q).sum();
      if (fabs(P - oldLocalizationMeasure) < _convThreshold or counter == maxSweeps)
        break;
      oldLocalizationMeasure = P;
      ++counter;
    }
  };

  orbitals.updateOrbitals(coefficients, orbitals.getEigenvalues());
}

template class PipekMezeyLocalization<Options::SCF_MODES::RESTRICTED>;
template class PipekMezeyLocalization<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
