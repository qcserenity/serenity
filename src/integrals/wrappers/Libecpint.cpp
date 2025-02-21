/**
 * @file   Libecpint.cpp
 *
 * @date   Dec 1, 2017
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
#include "integrals/wrappers/Libecpint.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/Basis.h" //Loop shells.
#include "basis/BasisController.h"
#include "basis/CartesianToSphericalTransformer.h"
#include "basis/Shell.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/SPMatrix.h"
#include "geometry/Atom.h"
#include "integrals/Normalization.h"
#include "misc/WarningTracker.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <array>
#include <iostream>
#include <libecpint/ecp.hpp>
#include <libecpint/ecpint.hpp>
#include <libecpint/gshell.hpp>
#include <libecpint/multiarr.hpp>

namespace Serenity {

MatrixInBasis<Options::SCF_MODES::RESTRICTED> Libecpint::computeECPIntegrals(std::shared_ptr<BasisController> basisController,
                                                                             const std::vector<std::shared_ptr<Atom>>& atoms) {
  MatrixInBasis<Options::SCF_MODES::RESTRICTED> result(basisController);
  result.block(0, 0, result.rows(), result.cols()) =
      computeECPIntegrals(basisController, basisController, atoms).block(0, 0, result.rows(), result.cols());
  return result;
}

SPMatrix<Options::SCF_MODES::RESTRICTED> Libecpint::computeECPIntegrals(std::shared_ptr<BasisController> basisControllerA,
                                                                        std::shared_ptr<BasisController> basisControllerB,
                                                                        const std::vector<std::shared_ptr<Atom>>& atoms) {
  bool singleBasisMode = basisControllerA == basisControllerB;
  unsigned int nBasisFunctionsA = basisControllerA->getNBasisFunctions();
  unsigned int nBasisFunctionsB = basisControllerB->getNBasisFunctions();

  SPMatrix<Options::SCF_MODES::RESTRICTED> result(nBasisFunctionsA, nBasisFunctionsB);
  const auto& basisA = basisControllerA->getBasis();
  const auto& basisB = basisControllerB->getBasis();
  const unsigned int nShellsA = basisControllerA->getReducedNBasisFunctions();
  const unsigned int nShellsB = basisControllerB->getReducedNBasisFunctions();
  // Initialize library
  unsigned int maxAngMomECP = 0;
  for (auto atom : atoms) {
    if (atom->usesECP() && (unsigned int)atom->getCorePotential()->getL() > maxAngMomECP) {
      maxAngMomECP = atom->getCorePotential()->getL();
      auto bar = atom->getCorePotential();
    }
  }
  unsigned int maxAngMomBasis = (basisControllerA->getMaxAngularMomentum() < basisControllerB->getMaxAngularMomentum())
                                    ? basisControllerB->getMaxAngularMomentum()
                                    : basisControllerA->getMaxAngularMomentum();
  // Parallel region
#pragma omp parallel
  {
    libecpint::ECPIntegral ecpint(maxAngMomBasis, maxAngMomECP);
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < nShellsA; ++i) {
      // Create shell in format for the library
      auto shellI = basisA[i];
      auto gshellI = makeShell(atoms, shellI);
      const unsigned int firstI = basisControllerA->extendedIndex(i);
      const unsigned int nI = N_SHELL_CART[shellI->getAngularMomentum()];
      unsigned int loopEnd = (singleBasisMode) ? i : nShellsB - 1;
      for (unsigned int j = 0; j <= loopEnd; ++j) {
        // Create shell in format for the library
        auto shellJ = basisB[j];
        auto gshellJ = makeShell(atoms, shellJ);
        const unsigned int firstJ = basisControllerB->extendedIndex(j);
        const unsigned int nJ = N_SHELL_CART[shellJ->getAngularMomentum()];
        // Loop over atoms
        for (const auto& atom : atoms) {
          if (atom->usesECP()) {
            libecpint::TwoIndex<double> shellResult(nI, nJ, 0.0);
            // calculate
            ecpint.compute_shell_pair(*atom->getCorePotential(), gshellI, gshellJ, shellResult);
            // Parse to Eigen::Matrix for easier handling.
            Eigen::MatrixXd eigenShellResult(nI, nJ);
            for (unsigned int ii = 0; ii < nI; ++ii) {
              for (unsigned int jj = 0; jj < nJ; ++jj) {
                double normFactor = 1.0;
                if (shellI->isCartesian())
                  normFactor *= shellI->getNormFactors()[ii];
                if (shellJ->isCartesian())
                  normFactor *= shellJ->getNormFactors()[jj];
                eigenShellResult(ii, jj) = normFactor * shellResult(ii, jj);
              }
            }
            // Transform to spherical harmonics if necessary.
            if (shellI->isSpherical()) {
              eigenShellResult =
                  CartesianToSphericalTransformer::getTransformationMatrix(shellI->getAngularMomentum()).transpose() *
                  eigenShellResult;
            }
            if (shellJ->isSpherical()) {
              eigenShellResult = eigenShellResult *
                                 CartesianToSphericalTransformer::getTransformationMatrix(shellJ->getAngularMomentum());
            }
            // unpack
            result.block(firstI, firstJ, eigenShellResult.rows(), eigenShellResult.cols()) += eigenShellResult;
            if (i != j && singleBasisMode)
              result.block(firstJ, firstI, eigenShellResult.cols(), eigenShellResult.rows()) += eigenShellResult.transpose();
          }
        } // atoms loop
      }   // j loop
    }     // i loop
  }       /* Parallel region */
  return result;
}

Eigen::MatrixXd Libecpint::computeECPGradientContribution(std::shared_ptr<AtomCenteredBasisController> basisController,
                                                          const std::vector<std::shared_ptr<Atom>>& atoms,
                                                          const DensityMatrix<RESTRICTED>& density) {
  Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(atoms.size(), 3);
  // Quick check if all of this can be skipped.
  bool skip = true;
  for (const auto& atom : atoms) {
    if (atom->usesECP()) {
      skip = false;
    }
  }
  if (skip) {
    return gradients;
  }
  auto mapping = basisController->getAtomIndicesOfBasisShells();
  const auto& basis = basisController->getBasis();
  const unsigned int nShells = basisController->getReducedNBasisFunctions();
  // Initialize library
  unsigned int maxAngMomECP = 0;
  for (auto atom : atoms) {
    if (atom->usesECP() && (unsigned int)atom->getCorePotential()->getL() > maxAngMomECP) {
      maxAngMomECP = atom->getCorePotential()->getL();
    }
  }
  unsigned int maxAngMomBasis = basisController->getMaxAngularMomentum();
  // Parallel region
#pragma omp parallel
  {
    Eigen::MatrixXd gradientContrPriv = Eigen::MatrixXd::Zero(atoms.size(), 3);
    libecpint::ECPIntegral ecpint(maxAngMomBasis, maxAngMomECP, 1);
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < nShells; ++i) {
      // Create shell in format for the library
      auto shellI = basis[i];
      auto gshellI = makeShell(atoms, shellI);
      const unsigned int firstI = basisController->extendedIndex(i);
      const unsigned int nI = N_SHELL_CART[shellI->getAngularMomentum()];
      for (unsigned int j = 0; j <= i; ++j) {
        // Create shell in format for the library
        auto shellJ = basis[j];
        auto gshellJ = makeShell(atoms, shellJ);
        const unsigned int firstJ = basisController->extendedIndex(j);
        const unsigned int nJ = N_SHELL_CART[shellJ->getAngularMomentum()];
        // Loop over atoms
        for (unsigned int iAtom = 0; iAtom < atoms.size(); iAtom++) {
          auto atom = atoms[iAtom];
          if (atom->usesECP()) {
            std::array<libecpint::TwoIndex<double>, 9> shellResult;
            ecpint.compute_shell_pair_derivative(*atom->getCorePotential(), gshellI, gshellJ, shellResult);

            for (unsigned int deriv = 0; deriv < 9; deriv++) {
              // Parse to Eigen::Matrix for easier handling.
              Eigen::MatrixXd eigenShellResult(nI, nJ);
              for (unsigned int ii = 0; ii < nI; ++ii) {
                for (unsigned int jj = 0; jj < nJ; ++jj) {
                  double normFactor = 1.0;
                  if (shellI->isCartesian())
                    normFactor *= shellI->getNormFactors()[ii];
                  if (shellJ->isCartesian())
                    normFactor *= shellJ->getNormFactors()[jj];
                  eigenShellResult(ii, jj) = normFactor * shellResult[deriv](ii, jj);
                }
              }
              // Transform to spherical harmonics if necessary.
              if (shellI->isSpherical()) {
                eigenShellResult =
                    CartesianToSphericalTransformer::getTransformationMatrix(shellI->getAngularMomentum()).transpose() *
                    eigenShellResult;
              }
              if (shellJ->isSpherical()) {
                eigenShellResult = eigenShellResult *
                                   CartesianToSphericalTransformer::getTransformationMatrix(shellJ->getAngularMomentum());
              }
              // unpack
              const unsigned int dir = deriv % 3;
              const double pref = (i != j) ? 2.0 : 1.0;
              if (deriv < 3) {
                const double contr =
                    pref * density.block(firstI, firstJ, eigenShellResult.rows(), eigenShellResult.cols())
                               .cwiseProduct(eigenShellResult)
                               .sum();
                gradientContrPriv(mapping[i], dir) += contr;
              }
              else if (deriv < 6) {
                const double contr =
                    pref * density.block(firstI, firstJ, eigenShellResult.rows(), eigenShellResult.cols())
                               .cwiseProduct(eigenShellResult)
                               .sum();
                gradientContrPriv(mapping[j], dir) += contr;
              }
              else {
                const double contr =
                    pref * density.block(firstI, firstJ, eigenShellResult.rows(), eigenShellResult.cols())
                               .cwiseProduct(eigenShellResult)
                               .sum();
                gradientContrPriv(iAtom, dir) += contr;
              }
            }
          }
        } // atoms loop
      }   // j loop
    }     // i loop
#pragma omp critical
    { gradients += gradientContrPriv; }
  } /* Parallel region */
  return gradients;
}

libecpint::GaussianShell Libecpint::makeShell(const std::vector<std::shared_ptr<Atom>>& atoms,
                                              std::shared_ptr<const Serenity::Shell> shell) {
  auto pos = shell->O;
  for (const auto& atom : atoms) {
    if (atom->usesECP()) {
      Eigen::Vector3d shellCoords;
      shellCoords << pos[0], pos[1], pos[2];
      Eigen::Vector3d atomCoords = atom->coords();
      const double dist = (shellCoords - atomCoords).norm();
      if (dist < 1e-3 && dist > 1e-10) {
        WarningTracker::printWarning("Warning: A center of an ECP and a basis function shell are close but not "
                                     "identical. I will assume identical centers.",
                                     iOOptions.printSCFCycleInfo);
        pos[0] = atomCoords(0);
        pos[1] = atomCoords(1);
        pos[2] = atomCoords(2);
      }
    }
  } // for atoms
  libecpint::GaussianShell gshell({pos[0], pos[1], pos[2]}, shell->getAngularMomentum());
  const auto& contractions = shell->contr[0].coeff;
  const auto& exponents = shell->alpha;
  for (unsigned int iPrim = 0; iPrim < shell->getNPrimitives(); ++iPrim) {
    gshell.addPrim(exponents[iPrim], contractions[iPrim]);
  }
  return gshell;
}

} /* namespace Serenity */
