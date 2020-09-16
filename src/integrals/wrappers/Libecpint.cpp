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
#include "integrals/wrappers/Libecpint.h"

// Serenity includes
#include "basis/BasisController.h"
#include "basis/CartesianToSphericalTransformer.h"
#include "basis/Shell.h"
#include "data/matrices/MatrixInBasis.h"
#include "data/matrices/SPMatrix.h"
#include "geometry/Atom.h"
#include "integrals/Normalization.h"
#include "misc/WarningTracker.h"
#include "parameters/Constants.h"
// External includes
#include <iostream>
#include <libecpint/ecp.hpp>
#include <libecpint/ecpint.hpp>
#include <libecpint/gshell.hpp>
#include <libecpint/multiarr.hpp>

namespace Serenity {
using namespace std;
using namespace Options;

MatrixInBasis<SCF_MODES::RESTRICTED> Libecpint::computeECPIntegrals(std::shared_ptr<BasisController> basisController,
                                                                    std::vector<std::shared_ptr<Atom>> atoms) {
  MatrixInBasis<SCF_MODES::RESTRICTED> result(basisController);
  result.block(0, 0, result.rows(), result.cols()) =
      computeECPIntegrals(basisController, basisController, atoms).block(0, 0, result.rows(), result.cols());
  return result;
}

SPMatrix<Options::SCF_MODES::RESTRICTED> Libecpint::computeECPIntegrals(std::shared_ptr<BasisController> basisControllerA,
                                                                        std::shared_ptr<BasisController> basisControllerB,
                                                                        std::vector<std::shared_ptr<Atom>> atoms) {
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
      auto posI = shellI->O;
      for (const auto& atom : atoms) {
        if (atom->usesECP()) {
          Eigen::Vector3d shellCoords;
          shellCoords << posI[0], posI[1], posI[2];
          Eigen::Vector3d atomCoords = atom->coords();
          Eigen::Vector3d diffVec = shellCoords - atomCoords;
          double distants = diffVec.norm();
          if (distants < 1e-3 && distants > 1e-10) {
            WarningTracker::printWarning("Warning: A center of an ECP and a basis function shell are close but not "
                                         "identical. I will assume identical centeres.",
                                         iOOptions.printSCFCycleInfo);
            posI[0] = atomCoords(0);
            posI[1] = atomCoords(1);
            posI[2] = atomCoords(2);
          }
        }
      } // for atoms
      libecpint::GaussianShell gshellI(posI.data(), shellI->getAngularMomentum());
      const auto& contractions = shellI->contr[0].coeff;
      const auto& exponents = shellI->alpha;
      for (unsigned int iPrim = 0; iPrim < shellI->getNPrimitives(); ++iPrim) {
        gshellI.addPrim(exponents[iPrim], contractions[iPrim]);
      }
      //      libecpint::GaussianShell gshellI = buildLibecpintGShell(shellI,atoms);
      const unsigned int firstI = basisControllerA->extendedIndex(i);
      const unsigned int nI = N_SHELL_CART[shellI->getAngularMomentum()];
      unsigned int loopEnd = (singleBasisMode) ? i : nShellsB - 1;
      for (unsigned int j = 0; j <= loopEnd; ++j) {
        // Create shell in format for the library
        auto shellJ = basisB[j];
        auto posJ = shellJ->O;
        for (const auto& atom : atoms) {
          if (atom->usesECP()) {
            Eigen::Vector3d shellCoords;
            shellCoords << posJ[0], posJ[1], posJ[2];
            Eigen::Vector3d atomCoords = atom->coords();
            Eigen::Vector3d diffVec = shellCoords - atomCoords;
            double distants = diffVec.norm();
            if (distants < 1e-3 && distants > 1e-10) {
              posJ[0] = atomCoords(0);
              posJ[1] = atomCoords(1);
              posJ[2] = atomCoords(2);
            }
          }
        } // for atoms
        libecpint::GaussianShell gshellJ(posJ.data(), shellJ->getAngularMomentum());
        const auto& contractionsJ = shellJ->contr[0].coeff;
        const auto& exponentsJ = shellJ->alpha;
        for (unsigned int jPrim = 0; jPrim < shellJ->getNPrimitives(); ++jPrim) {
          gshellJ.addPrim(exponentsJ[jPrim], contractionsJ[jPrim]);
        }
        //        libecpint::GaussianShell gshellJ = buildLibecpintGShell(shellJ,atoms);
        const unsigned int firstJ = basisControllerB->extendedIndex(j);
        const unsigned int nJ = N_SHELL_CART[shellJ->getAngularMomentum()];
        // Loop over atoms
        for (const auto atom : atoms) {
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

} /* namespace Serenity */
