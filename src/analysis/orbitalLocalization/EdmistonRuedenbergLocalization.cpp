/**
 * @file EdmistonRuedenbergLocalization.cpp
 *
 * @date Nov 3, 2016
 * @author David Schnieders
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
#include "analysis/orbitalLocalization/EdmistonRuedenbergLocalization.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/transformer/Ao2MoTransformer.h"
#include "math/RegularRankFourTensor.h"
#include "math/linearAlgebra/JacobiRotation.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <cmath>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
EdmistonRuedenbergLocalization<SCFMode>::EdmistonRuedenbergLocalization(std::shared_ptr<SystemController> systemController)
  : _systemController(systemController) {
}

template<Options::SCF_MODES SCFMode>
void EdmistonRuedenbergLocalization<SCFMode>::localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                                               SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) {
  const auto& nOccOrbs = _systemController->getNOccupiedOrbitals<SCFMode>();

  CoefficientMatrix<SCFMode> coefficients(orbitals.getCoefficients());
  auto basisController = _systemController->getBasisController();
  auto nBasisFunc = basisController->getNBasisFunctions();
  RegularRankFourTensor<double> eris(nBasisFunc, 0.0);
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, basisController, basisController->getPrescreeningThreshold());
  Ao2MoTransformer aoToMo(basisController);

  // Calculate and store 4 center integrals - Memory demanding!
  auto const storeERIS = [&eris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                 const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
    (void)threadId; // no warning, please
    eris(b, a, i, j) = integral(0);
    eris(b, a, j, i) = integral(0);
    eris(a, b, j, i) = integral(0);
    eris(a, b, i, j) = integral(0);
    eris(i, j, b, a) = integral(0);
    eris(i, j, a, b) = integral(0);
    eris(j, i, b, a) = integral(0);
    eris(j, i, a, b) = integral(0);
  };
  looper.loop(storeERIS);

  for_spin(coefficients, nOccOrbs, orbitalRange) {
    unsigned int cycle = 0;

    while (true) {
      // ao -> mo transformation of 4 center integrals
      RegularRankFourTensor<double> moERIS(nOccOrbs_spin, 0.0);
      aoToMo.transformTwoElectronIntegrals(eris, moERIS, coefficients_spin, nOccOrbs_spin);

      // A matrix to store a potential gain in SOS for pairwise rotations
      Eigen::MatrixXd selfRepIncrease = Eigen::MatrixXd::Zero(nOccOrbs_spin, nOccOrbs_spin);

      // A matrix to store the rotation Angles
      Eigen::MatrixXd rotAngles = Eigen::MatrixXd::Zero(nOccOrbs_spin, nOccOrbs_spin);

      // Calculate potential increase of self repulsion and rotation angles from electron repulsion integrals
      for (unsigned int iOrb = 0; iOrb < orbitalRange_spin.size(); ++iOrb) {
        unsigned int i = orbitalRange_spin[iOrb];
        for (unsigned int jOrb = 0; jOrb < iOrb; ++jOrb) {
          unsigned int j = orbitalRange_spin[jOrb];
          double factorA = 0.0;
          double factorB = 0.0;
          factorA = moERIS(i, j, i, j) - (moERIS(i, i, i, i) + moERIS(j, j, j, j) - 2 * moERIS(i, i, j, j)) / 4;
          factorB = moERIS(i, j, i, i) - moERIS(i, j, j, j);
          selfRepIncrease(i, j) = factorA + sqrt(factorA * factorA + factorB * factorB);
          rotAngles(i, j) = atan2(factorB, -factorA) / 4.0;
        }
      }

      // Find maximal increase in selfRepIncrease and its angle
      unsigned int rowOfMax, colOfMax;
      const double maximum = selfRepIncrease.maxCoeff(&rowOfMax, &colOfMax);
      const double maxAngle = rotAngles(rowOfMax, colOfMax);

      if (maximum < 1e-9)
        break;

      // Perform rotation on coefficients
      JacobiRotation::rotate(coefficients_spin.col(rowOfMax), coefficients_spin.col(colOfMax), maxAngle);

      // Rotate other pairs of orbitals which were untouched so far
      std::vector<bool> wasRotated(nOccOrbs_spin, false);
      wasRotated[rowOfMax] = true;
      wasRotated[colOfMax] = true;
      for (unsigned int i = 0; i < nOccOrbs_spin; ++i) {
        if (wasRotated[i])
          continue;
        // Identify maximum increase for this row
        unsigned int maxJ = 0;
        unsigned int iterTry = 0;
        while (iterTry < nOccOrbs_spin) {
          selfRepIncrease.row(i).maxCoeff(&maxJ);
          if (std::find(orbitalRange_spin.begin(), orbitalRange_spin.end(), maxJ) != orbitalRange_spin.end()) {
            break;
          }
          else {
            ++iterTry;
            selfRepIncrease(i, maxJ) = -99999;
          }
        }
        const double maxAngle = rotAngles(i, maxJ);
        // Do rotation
        JacobiRotation::rotate(coefficients_spin.col(i), coefficients_spin.col(maxJ), maxAngle);
        // Mark as rotated
        wasRotated[i] = true;
        wasRotated[maxJ] = true;
      }

      cycle++;

      if (cycle > maxSweeps)
        break;
    }
  };
  orbitals.updateOrbitals(coefficients, orbitals.getEigenvalues());
}

template class EdmistonRuedenbergLocalization<Options::SCF_MODES::RESTRICTED>;
template class EdmistonRuedenbergLocalization<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
