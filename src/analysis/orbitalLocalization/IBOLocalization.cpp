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
IBOLocalization<SCFMode>::IBOLocalization(std::shared_ptr<SystemController> systemController, bool IAOsOnly,
                                          bool replaceVirtuals, bool enforceRestrictedOrbitals)
  : _system(systemController),
    _IAOsOnly(IAOsOnly),
    _replaceVirtuals(replaceVirtuals),
    _enforceRestrictedOrbitals(enforceRestrictedOrbitals){};

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
  auto nOcc = _system->getNOccupiedOrbitals<SCFMode>();
  // Create new basis
  std::shared_ptr<AtomCenteredBasisController> minaoBasis = AtomCenteredBasisControllerFactory::produce(
      _system->getGeometry(), _system->getSettings().basis.basisLibPath,
      _system->getSettings().basis.makeSphericalBasis, false, _system->getSettings().basis.firstECP, "MINAO");

  // store the basis
  _system->setBasisController(minaoBasis, Options::BASIS_PURPOSES::IAO_LOCALIZATION);

  // bases
  auto B1 = _system->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  auto B2 = _system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION);

  bool localizeVirtuals = false;
  for_spin(orbitalRange, nOcc) {
    if (orbitalRange_spin.size() > 1) {
      unsigned int maxIndex = *std::max_element(orbitalRange_spin.begin(), orbitalRange_spin.end());
      if (maxIndex >= nOcc_spin)
        localizeVirtuals = true;
      if (maxIndex > B2->getNBasisFunctions())
        throw SerenityError("The IBO localization basis is too small to localize the given orbital selection.");
    }
  };

  // coefficent matrix
  CoefficientMatrix<SCFMode> C = orbitals.getCoefficients();
  // Overlap integrals.
  const auto& S1 = _system->getOneElectronIntegralController()->getOverlapIntegrals();
  /*
   * The following few lines are needed for the localization of virtual orbitals only. In the
   * case that the IAO basis set does not span the space of the virtual valence orbitals we must
   * first reconstruct the orbitals to do so. Otherwise, the orbital localization may lead to
   * non-orthogonal and non normalized virtual orbitals. In the case that "_replaceVirtuals" is
   * false, we explicitly check if the virtual valence orbital space is fine and replace the
   * virtual orbitals if required.
   */
  bool replaceVirtBeforeLoc =
      localizeVirtuals &&
      (_replaceVirtuals || !IAOPopulationCalculator<SCFMode>::iaosSpanOrbitals(
                               C, _system->getOneElectronIntegralController()->getOverlapIntegrals(),
                               _system->template getNOccupiedOrbitals<SCFMode>(), _system->getBasisController(),
                               _system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION)));
  if (replaceVirtBeforeLoc && !_virtualOrbitalsReplaced) {
    _virtualOrbitalsReplaced = true;
    /*
     * If we have a SOMO that should not mix with the remaining virtual orbital space, we must not include it in the
     * reconstruction of the virtual orbital space.
     */
    auto nOccValenceProjection =
        (_enforceRestrictedOrbitals) ? restrictedOccupations() : _system->template getNOccupiedOrbitals<SCFMode>();
    IAOPopulationCalculator<SCFMode>::reconstructVirtualValenceOrbitalsInplace(
        C, _system->getOneElectronIntegralController()->getOverlapIntegrals(), nOccValenceProjection,
        _system->getBasisController(), _system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION));
    orbitals.updateOrbitals(C, orbitals.getEigenvalues());
  }
  // Coefficients in IAO basis
  auto CIAOandOthoA = IAOPopulationCalculator<SCFMode>::getCIAOCoefficients(C, S1, nOcc, B1, B2, localizeVirtuals);

  auto& spinCIAO = CIAOandOthoA.first;
  const auto& spinOthoA = CIAOandOthoA.second;

  bool skipLoop = false;
  for_spin(spinOthoA, spinCIAO, C, orbitalRange) {
    if (skipLoop) {
      return;
    }

    auto& CIAO = spinCIAO_spin;
    const unsigned int nCIAOOrbs = CIAO.cols();
    const auto& othoA = spinOthoA_spin;
    if (_IAOsOnly) {
      /*
       * Stop here for IAOs
       */
      C_spin.leftCols(nCIAOOrbs) = othoA * CIAO;
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
             *
             *   We will do the following:
             *   If Aij = 0       --> angle = atan(0) = 0
             *   if |Bij| < 1e-14 --> Bij = |Bij| in order to fix the sing.
             *     and thus |angle| = pi/8
             */
            double angle = 0.0;
            if (std::abs(Aij) > 1e-14) {
              if (std::abs(Bij) < 1e-14)
                Bij = std::abs(Bij);
              angle = 0.25 * atan2(Bij, -Aij);
            }
            // rotate orbital pair
            JacobiRotation::rotate(CIAO.col(i), CIAO.col(j), angle);
          }
        }
        // Convergence Check
        if (fabs(gradient) <= 1e-8) {
          std::cout << "    Converged after " << cycle << " orbital rotation cycles." << std::endl << std::endl;
          C_spin.leftCols(nCIAOOrbs) = othoA * CIAO;
          converged = true;
        }
        else if (cycle == maxSweeps) {
          // The paper claims that 5-10 sweeps should converge the orbitals
          //  therefore the procedure will stop after 200 cycles.
          std::cout << "    ERROR: IBO procedure did not converged after " << cycle << " orbital rotation cycles."
                    << std::endl;
          std::cout << "           The orbitals will still be stored for error analysis." << std::endl;
          C_spin.leftCols(nCIAOOrbs) = othoA * CIAO;

          converged = true;
        }
        const double sanityCheck = (C_spin.leftCols(nCIAOOrbs).transpose() * S1 * C_spin.leftCols(nCIAOOrbs) -
                                    Eigen::MatrixXd::Identity(nCIAOOrbs, nCIAOOrbs))
                                       .array()
                                       .abs()
                                       .sum();
        if (sanityCheck > 1e-6)
          throw SerenityError("The orbitals are no longer orthogonal after IBO localization! Maximum overlap: " +
                              std::to_string(sanityCheck));
      }
    } /* rotations while loop */
    if (_enforceRestrictedOrbitals && SCFMode != RESTRICTED) {
      skipLoop = true;
      this->restrictOrbitals(C);
    }
  }; /* spin loop */
  orbitals.updateOrbitals(C, orbitals.getEigenvalues());
}
template<>
void IBOLocalization<UNRESTRICTED>::restrictOrbitals(CoefficientMatrix<UNRESTRICTED>& coefficientMatrix) {
  coefficientMatrix.beta = coefficientMatrix.alpha;
}
template<>
void IBOLocalization<RESTRICTED>::restrictOrbitals(CoefficientMatrix<RESTRICTED>&) {
  return;
}
template<>
SpinPolarizedData<RESTRICTED, unsigned int> IBOLocalization<RESTRICTED>::restrictedOccupations() {
  return _system->getNOccupiedOrbitals<RESTRICTED>();
}
template<>
SpinPolarizedData<UNRESTRICTED, unsigned int> IBOLocalization<UNRESTRICTED>::restrictedOccupations() {
  auto nOcc = _system->getNOccupiedOrbitals<UNRESTRICTED>();
  auto maxOcc = std::max(nOcc.alpha, nOcc.beta);
  nOcc.alpha = maxOcc;
  nOcc.beta = maxOcc;
  return nOcc;
}

template class IBOLocalization<Options::SCF_MODES::RESTRICTED>;
template class IBOLocalization<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
