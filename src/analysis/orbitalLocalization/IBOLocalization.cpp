/**
 * @file IBOLocalization.cpp
 *
 * @date Jun 15, 2016
 * @author Jan Unsleber
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
#include "analysis/orbitalLocalization/IBOLocalization.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/IAOPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/linearAlgebra/JacobiRotation.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
IBOLocalization<SCFMode>::IBOLocalization(std::shared_ptr<SystemController> systemController, bool IAOsOnly)
  : _system(systemController), _IAOsOnly(IAOsOnly){};

template<Options::SCF_MODES SCFMode>
void IBOLocalization<SCFMode>::localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                                SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) {
  /*
   * Most variables are named after their names in the
   *  paper describing the initial implementation.
   * Please refere to:
   * Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical Concepts
   * Gerald Knizia, J. Chem. Theory Comput., 2013, 9 (11), pp 4834–4843
   *
   * this link contains the reprint with corrected formulas in the appendix:
   * http://www.theochem.uni-stuttgart.de/~knizia/bin/iao_preprint.pdf
   */
  auto nOccOrbs = _system->getNOccupiedOrbitals<SCFMode>();
  // Create new basis
  std::shared_ptr<AtomCenteredBasisController> minaoBasis = AtomCenteredBasisControllerFactory::produce(
      _system->getGeometry(), _system->getSettings().basis.basisLibPath,
      _system->getSettings().basis.makeSphericalBasis, false, _system->getSettings().basis.firstECP, "MINAO");

  // store the basis
  _system->setBasisController(minaoBasis, Options::BASIS_PURPOSES::IAO_LOCALIZATION);

  // bases
  auto B1 = _system->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  auto B2 = _system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION);

  // coefficent matrix
  CoefficientMatrix<SCFMode> C = orbitals.getCoefficients();
  // Overlap integrals.
  const auto& S1 = _system->getOneElectronIntegralController()->getOverlapIntegrals();
  // Coefficients in IAO basis
  auto CIAOandOthoA = IAOPopulationCalculator<SCFMode>::getCIAOCoefficients(C, S1, nOccOrbs, B1, B2);
  auto& spinCIAO = CIAOandOthoA.first;
  const auto& spinOthoA = CIAOandOthoA.second;

  for_spin(spinOthoA, spinCIAO, C, nOccOrbs, orbitalRange) {
    auto& CIAO = spinCIAO_spin;
    const auto& othoA = spinOthoA_spin;
    if (_IAOsOnly) {
      /*
       * Stop here for IAOs
       */
      C_spin.leftCols(nOccOrbs_spin) = othoA * CIAO;
    }
    else {
      /*
       * Rotation from IAOs to IBOs
       */
      bool converged = false;
      unsigned int cycle = 0;
      while (!converged) {
        ++cycle;

        // In the reference implementation of Knizia it is explained that
        //  Bij is effectively the gradient, it is therefore used as convergence crit.
        double gradient = 0.0;

        /*
         * 2x2 rotation loop
         */
        for (unsigned int iOrb = 0; iOrb < orbitalRange_spin.size(); ++iOrb) {
          unsigned int i = orbitalRange_spin[iOrb];
          for (unsigned int jOrb = 0; jOrb < iOrb; ++jOrb) {
            unsigned int j = orbitalRange_spin[jOrb];
            double Aij = 0;
            double Bij = 0;
            for (unsigned int atomIndex = 0; atomIndex < _system->getGeometry()->getNAtoms(); ++atomIndex) {
              auto indices = minaoBasis->getBasisIndices();
              double Qii = 0;
              double Qij = 0;
              double Qjj = 0;
              for (unsigned int mu = indices[atomIndex].first; mu < indices[atomIndex].second; ++mu) {
                Qii += CIAO(mu, i) * CIAO(mu, i);
                Qij += CIAO(mu, j) * CIAO(mu, i);
                Qjj += CIAO(mu, j) * CIAO(mu, j);
              }
              // p=2 simmilar to PM
              // Aij += 4.0*Qij*Qij - (Qii - Qjj)*(Qii - Qjj);
              // Bij += 4.0*Qij*(Qii - Qjj);
              // Alternative p=4
              const double Qii3 = Qii * Qii * Qii;
              const double Qjj3 = Qjj * Qjj * Qjj;
              Aij += -Qii3 * Qii - Qjj3 * Qjj;
              Aij += 6.0 * (Qii * Qii + Qjj * Qjj) * Qij * Qij;
              Aij += Qii * Qjj3 + Qii3 * Qjj;
              Bij += 4.0 * Qij * (Qii3 - Qjj3);
            }
            gradient += Bij * Bij;
            /*
             * Calculate the angle:
             *   The atan2(y,x) function calculated the arcus tangent of
             *   y/x. Thus, it is better to check for numerical zeros, since
             *   atan(inf) = pi/2 and atan(-inf) = -pi/2.
             *   Otherwise the results may depend on the zero being +0 or -0.
             */
            double angle = 0.0;
            if (std::abs(Bij) > 1e-14) {
              if (std::abs(Aij) < 1e-14)
                Aij = std::abs(Aij);
              angle = 0.25 * atan2(Bij, -Aij);
            }
            // rotate orbital pair
            JacobiRotation::rotate(CIAO.col(i), CIAO.col(j), angle);
          }
        }
        // Convergence Check
        if (fabs(gradient) <= 1e-8) {
          std::cout << "    Converged after " << cycle << " orbital rotation cycles." << std::endl << std::endl;
          C_spin.leftCols(nOccOrbs_spin) = othoA * CIAO;
          converged = true;
        }
        else if (cycle == maxSweeps) {
          // The paper claims that 5-10 sweeps should converge the orbitals
          //  therefore the procedure will stop after 200 cycles.
          std::cout << "    ERROR: IBO procedure did not converged after " << cycle << " orbital rotation cycles."
                    << std::endl;
          std::cout << "           The orbitals will still be stored for error analysis." << std::endl;
          C_spin.leftCols(nOccOrbs_spin) = othoA * CIAO;
          converged = true;
        }
      }
    } /* rotations while loop */
  };  /* spin loop */
  orbitals.updateOrbitals(C, orbitals.getEigenvalues());
}

template class IBOLocalization<Options::SCF_MODES::RESTRICTED>;
template class IBOLocalization<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
