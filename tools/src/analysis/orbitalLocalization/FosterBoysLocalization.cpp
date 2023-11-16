/**
 * @file FosterBoysLocalization.cpp
 *
 * @date Nov 5, 2015
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
#include "analysis/orbitalLocalization/FosterBoysLocalization.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/Matrix.h"
#include "math/linearAlgebra/JacobiRotation.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <algorithm> // std::find

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FosterBoysLocalization<SCFMode>::FosterBoysLocalization(std::shared_ptr<SystemController> systemController)
  : _system(systemController) {
}

template<Options::SCF_MODES SCFMode>
void FosterBoysLocalization<SCFMode>::localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                                       SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) {
  CoefficientMatrix<SCFMode> coefficients = orbitals.getCoefficients();
  const auto& nOccOrbs = _system->getNOccupiedOrbitals<SCFMode>();

  auto basisController = _system->getBasisController();
  const unsigned int nBFs = basisController->getNBasisFunctions();
  unsigned int cycle = 0;

  auto dipoles = _system->getOneElectronIntegralController()->getDipoleLengths();

  for_spin(coefficients, nOccOrbs, orbitalRange) {
    while (true) {
      // AO dipole Integrals -> MO dipole Integrals
      std::vector<Eigen::MatrixXd> moDipoleInt(3, Eigen::MatrixXd(nBFs, nBFs));
      moDipoleInt[0] = -coefficients_spin.transpose() * dipoles[0] * coefficients_spin;
      moDipoleInt[1] = -coefficients_spin.transpose() * dipoles[1] * coefficients_spin;
      moDipoleInt[2] = -coefficients_spin.transpose() * dipoles[2] * coefficients_spin;

      // A matrix to store a potential gain in SOS for pairwise rotations
      Eigen::MatrixXd sosIncrease = Eigen::MatrixXd::Zero(nOccOrbs_spin, nOccOrbs_spin);
      // A matrix to store the rotation angles
      Eigen::MatrixXd rotAngles = Eigen::MatrixXd::Zero(nOccOrbs_spin, nOccOrbs_spin);

      // Calculate potential increase of SOS and rotation angles from dipole integrals
      for (unsigned int iOrb = 0; iOrb < orbitalRange_spin.size(); ++iOrb) {
        unsigned int i = orbitalRange_spin[iOrb];
        for (unsigned int jOrb = 0; jOrb < iOrb; ++jOrb) {
          unsigned int j = orbitalRange_spin[jOrb];
          double factorA = 0.0;
          double factorB = 0.0;
          for (unsigned int m = 0; m < 3; m++) {
            const double diff = moDipoleInt[m](i, i) - moDipoleInt[m](j, j);
            const double mij = moDipoleInt[m](i, j);
            factorA += mij * mij - 0.25 * diff * diff;
            factorB += diff * mij;
          }
          /*
           *  see: http://en.cppreference.com/w/cpp/numeric/math/atan2
           *  for information on how the numerics of atan2 work
           */
          if (fabs(factorA) < 1e-12)
            factorA = 0.0;
          if (fabs(factorB) < 1e-12)
            factorB = 0.0;
          sosIncrease(i, j) = factorA + sqrt(factorA * factorA + factorB * factorB);
          rotAngles(i, j) = 0.25 * atan2(factorB, -1.0 * factorA);
        }
      }

      // Find maximal increase in SOS and its angle
      unsigned int rowOfMax, colOfMax;
      const double maximum = sosIncrease.maxCoeff(&rowOfMax, &colOfMax);
      const double maxAngle = rotAngles(rowOfMax, colOfMax);

      if (maximum < 1e-9)
        break;

      // Perform rotation on coefficients
      JacobiRotation::rotate(coefficients_spin.col(rowOfMax), coefficients_spin.col(colOfMax), maxAngle);

      // Rotate other pairs of orbitals which were untouched so far
      std::vector<bool> wasRotated(nOccOrbs_spin, false);
      wasRotated[rowOfMax] = true;
      wasRotated[colOfMax] = true;
      for (unsigned int iOrb = 0; iOrb < orbitalRange_spin.size(); ++iOrb) {
        unsigned int i = orbitalRange_spin[iOrb];
        if (wasRotated[i])
          continue;
        // Identify maximum increase for this row
        unsigned int maxJ = 0;
        unsigned int iterTry = 0;
        while (iterTry < nOccOrbs_spin) {
          sosIncrease.row(i).maxCoeff(&maxJ);
          if (std::find(orbitalRange_spin.begin(), orbitalRange_spin.end(), maxJ) != orbitalRange_spin.end()) {
            break;
          }
          else {
            ++iterTry;
            sosIncrease(i, maxJ) = -99999;
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

  }; /* for_spin */

  // CMOs -> LMOs
  orbitals.updateOrbitals(coefficients, orbitals.getEigenvalues());
};

template class FosterBoysLocalization<Options::SCF_MODES::RESTRICTED>;
template class FosterBoysLocalization<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
