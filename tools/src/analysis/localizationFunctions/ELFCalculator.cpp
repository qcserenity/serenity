/**
 * @file ELFCalculator.cpp
 *
 * @date Apr 4, 2016, rework Jun 28, 2017
 * @author Melanie Börner, rework Jan Unsleber
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
#include "analysis/localizationFunctions/ELFCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/GridData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "grid/GridController.h"
#include "io/CubeFileWriter.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <cmath>
#include <vector>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ELFCalculator<SCFMode>::ELFCalculator(std::shared_ptr<SystemController> systemController)
  : _systemController(systemController), _basisController(systemController->getAtomCenteredBasisController()) {
  assert(_systemController);
  assert(_basisController);
}

template<Options::SCF_MODES SCFMode>
GridData<SCFMode> ELFCalculator<SCFMode>::calcKineticEnergyDensity(std::shared_ptr<GridController> gridController) {
  assert(gridController);

  // Definitions
  CoefficientMatrix<SCFMode> coefficientMatrix(
      _systemController->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->getCoefficients());
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(_systemController->getSettings(), _basisController, gridController);
  unsigned int nBlocks = basisFunctionOnGridController->getNBlocks();
  auto nOccOrbitals = _systemController->getNOccupiedOrbitals<SCFMode>();

  // Calculation of the kinetic energy density
  GridData<SCFMode> kineticEnergyDensity(gridController);
  for_spin(nOccOrbitals, coefficientMatrix, kineticEnergyDensity) {
    // Loop over the orbitals
    for (unsigned int i = 0; i < nOccOrbitals_spin; ++i) {
      // Loop over Grid. Data is accessed blockwise.
      for (unsigned int blockIndex = 0; blockIndex < nBlocks; ++blockIndex) {
        const auto& blockOnGridData = basisFunctionOnGridController->getBlockOnGridData(blockIndex);
        const auto& bfDerivatives = *blockOnGridData->derivativeValues;
        const unsigned int blockSize = blockOnGridData->functionValues.rows();

        // Loop over data entries in the current block
        for (unsigned int blockData = 0; blockData < blockSize; ++blockData) {
          const unsigned int pointIndex = basisFunctionOnGridController->getFirstIndexOfBlock(blockIndex) + blockData;
          Eigen::Vector3d gradPhi;
          gradPhi[0] = coefficientMatrix_spin.col(i).dot(bfDerivatives.x.row(blockData));
          gradPhi[1] = coefficientMatrix_spin.col(i).dot(bfDerivatives.y.row(blockData));
          gradPhi[2] = coefficientMatrix_spin.col(i).dot(bfDerivatives.z.row(blockData));
          kineticEnergyDensity_spin[pointIndex] += (SCFMode == RESTRICTED ? 1.0 : 0.5) * gradPhi.squaredNorm();
        } /* Loop over the entries in the current block */
      }   /* Loop over the grid blocks.*/
    }     /* Loop over the Orbitals */
  };
  return kineticEnergyDensity;
} /* calcKineticEnergyDensity */

template<Options::SCF_MODES SCFMode>
GridData<SCFMode> ELFCalculator<SCFMode>::calculateELFOnGrid(std::shared_ptr<GridController> gridController) {
  assert(gridController);

  // Definitions
  GridData<SCFMode> elf(gridController);
  auto electronicStructure = _systemController->getElectronicStructure<SCFMode>();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(_systemController->getSettings(), _basisController, gridController);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
      basisFunctionOnGridController, _systemController->getSettings().grid.blockAveThreshold);
  const unsigned int nGridPoints = gridController->getNGridPoints();
  auto densityMatrix(electronicStructure->getDensityMatrix());
  auto densityOnGrid = densOnGridCalculator->calcDensityOnGrid(electronicStructure->getDensityMatrix());
  auto densityGradientOnGrid = makeGradient<DensityOnGrid<SCFMode>>(gridController);

  // Calculation of the density and its derivatives on a grid
  densOnGridCalculator->calcDensityAndGradientOnGrid(densityMatrix, densityOnGrid, densityGradientOnGrid);

  GridData<SCFMode> nablaDensityOnGrid2(gridController);
  auto& gradientx = densityGradientOnGrid.x;
  auto& gradienty = densityGradientOnGrid.y;
  auto& gradientz = densityGradientOnGrid.z;
  for_spin(nablaDensityOnGrid2, gradientx, gradienty, gradientz) {
    nablaDensityOnGrid2_spin.array() = gradientx_spin.array() * gradientx_spin.array() +
                                       gradienty_spin.array() * gradienty_spin.array() +
                                       gradientz_spin.array() * gradientz_spin.array();
  };

  // Calculation of the kinetic energy density
  auto kineticEnergyDensity = this->calcKineticEnergyDensity(gridController);

  const double prefactor = (3.0 / 10.0) * pow(((SCFMode == RESTRICTED ? 3.0 : 6.0) * M_PI * M_PI), (2.0 / 3.0));

  for_spin(densityOnGrid, nablaDensityOnGrid2, kineticEnergyDensity, elf) {
    // Loop over all grid points
    for (unsigned int i = 0; i < nGridPoints; ++i) {
      // Calculation of d_0(r)
      const double d_0 = prefactor * pow(densityOnGrid_spin[i], (5.0 / 3.0));
      // Calculation of d
      const double d = (densityOnGrid_spin[i] > 1.0e-10)
                           ? kineticEnergyDensity_spin[i] - (nablaDensityOnGrid2_spin[i] / (8.0 * densityOnGrid_spin[i]))
                           : 0.0;
      // Calculation of elf
      elf_spin[i] = (d_0 > 1.0e-10) ? 1.0 / (1.0 + (d * d) / (d_0 * d_0)) : 0.0;
    }
  };
  return elf;
} /* calculateELFOnGrid */

template<>
GridData<RESTRICTED> ELFCalculator<RESTRICTED>::calculateTotalELFOnGrid(std::shared_ptr<GridController> gridController) {
  return this->calculateELFOnGrid(gridController);
}

template<>
GridData<RESTRICTED> ELFCalculator<UNRESTRICTED>::calculateTotalELFOnGrid(std::shared_ptr<GridController> gridController) {
  assert(gridController);

  // Definitions
  GridData<RESTRICTED> totalELF(gridController);
  auto electronicStructure = _systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(_systemController->getSettings(), _basisController, gridController);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>>(
      basisFunctionOnGridController, _systemController->getSettings().grid.blockAveThreshold);
  const unsigned int nGridPoints = gridController->getNGridPoints();
  auto densityMatrix(electronicStructure->getDensityMatrix());
  auto densityOnGrid = densOnGridCalculator->calcDensityOnGrid(electronicStructure->getDensityMatrix());
  auto densityGradientOnGrid = makeGradient<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>>(gridController);

  // Calculation of the density and its derivatives on a grid
  densOnGridCalculator->calcDensityAndGradientOnGrid(densityMatrix, densityOnGrid, densityGradientOnGrid);

  // Calculation of the kinetic energy density
  auto kineticEnergyDensity = this->calcKineticEnergyDensity(gridController);

  GridData<UNRESTRICTED> nablaDensityOnGrid2(gridController);
  auto& gradientx = densityGradientOnGrid.x;
  auto& gradienty = densityGradientOnGrid.y;
  auto& gradientz = densityGradientOnGrid.z;
  for_spin(nablaDensityOnGrid2, gradientx, gradienty, gradientz) {
    nablaDensityOnGrid2_spin.array() = gradientx_spin.array() * gradientx_spin.array() +
                                       gradienty_spin.array() * gradienty_spin.array() +
                                       gradientz_spin.array() * gradientz_spin.array();
  };

  const double prefactor = (3.0 / 10.0) * pow((6.0 * M_PI * M_PI), (2.0 / 3.0));

  // Loop over all grid points
  for (unsigned int i = 0; i < nGridPoints; ++i) {
    // Calculation of d_0(r)
    const double d_0 = prefactor * (pow(densityOnGrid.alpha[i], (5.0 / 3.0)) + pow(densityOnGrid.beta[i], (5.0 / 3.0)));
    // Calculation of d
    double d = 0.0;
    if (densityOnGrid.alpha[i] > 1.0e-10 && densityOnGrid.beta[i] > 1.0e-10) {
      d = kineticEnergyDensity.alpha[i] + kineticEnergyDensity.beta[i] -
          (nablaDensityOnGrid2.alpha[i] / (8.0 * densityOnGrid.alpha[i])) -
          (nablaDensityOnGrid2.beta[i] / (8.0 * densityOnGrid.beta[i]));
    }
    // Calculation of elf
    totalELF[i] = (d_0 > 1.0e-10) ? 1.0 / (1.0 + (d * d) / (d_0 * d_0)) : 0.0;
  }
  return totalELF;
} /* calculateTotalELFOnGrid */

template<Options::SCF_MODES SCFMode>
GridData<SCFMode> ELFCalculator<SCFMode>::calculateELFTSOnGrid(std::shared_ptr<GridController> gridController) {
  assert(gridController);

  // Definitions
  GridData<SCFMode> elf_TS(gridController);
  auto electronicStructure = _systemController->getElectronicStructure<SCFMode>();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(_systemController->getSettings(), _basisController, gridController);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
      basisFunctionOnGridController, _systemController->getSettings().grid.blockAveThreshold);
  const unsigned int nGridPoints = gridController->getNGridPoints();
  auto densityMatrix(electronicStructure->getDensityMatrix());
  auto densityOnGrid = densOnGridCalculator->calcDensityOnGrid(electronicStructure->getDensityMatrix());
  auto densityGradientOnGrid = makeGradient<DensityOnGrid<SCFMode>>(gridController);
  auto densityHessianOnGrid = makeHessian<DensityOnGrid<SCFMode>>(gridController);

  // Calculation of the density and its derivatives on grid
  densOnGridCalculator->calcDensityAndDerivativesOnGrid(densityMatrix, densityOnGrid, densityGradientOnGrid,
                                                        densityHessianOnGrid);

  // Calculation of the x, y and z parts of the gradient and the laplacian of the electron density on a grid
  auto& densityGradientOnGridx = densityGradientOnGrid.x;
  auto& densityGradientOnGridy = densityGradientOnGrid.y;
  auto& densityGradientOnGridz = densityGradientOnGrid.z;
  auto& densityHessianOnGridxx = densityHessianOnGrid.xx;
  auto& densityHessianOnGridyy = densityHessianOnGrid.yy;
  auto& densityHessianOnGridzz = densityHessianOnGrid.zz;

  // Calculation of |grad(rho)|^2
  GridData<SCFMode> nablaDensityOnGrid2(gridController);
  for_spin(nablaDensityOnGrid2, densityGradientOnGridx, densityGradientOnGridy, densityGradientOnGridz) {
    nablaDensityOnGrid2_spin.array() = (densityGradientOnGridx_spin.array() * densityGradientOnGridx_spin.array()) +
                                       (densityGradientOnGridy_spin.array() * densityGradientOnGridy_spin.array()) +
                                       (densityGradientOnGridz_spin.array() * densityGradientOnGridz_spin.array());
  };

  // Calculation of laplace(roh)
  GridData<SCFMode> laplaceDensityOnGrid(gridController);
  for_spin(laplaceDensityOnGrid, densityHessianOnGridxx, densityHessianOnGridyy, densityHessianOnGridzz) {
    laplaceDensityOnGrid_spin.array() =
        densityHessianOnGridxx_spin.array() + densityHessianOnGridyy_spin.array() + densityHessianOnGridzz_spin.array();
  };

  // Loop over all grid points
  for_spin(laplaceDensityOnGrid, densityOnGrid, nablaDensityOnGrid2, elf_TS) {
    for (unsigned int i = 0; i < nGridPoints; ++i) {
      // Calculation of d_0(r)
      const double d_0 = 0.3 * pow(((SCFMode == RESTRICTED ? 3.0 : 6.0) * M_PI * M_PI), (2.0 / 3.0)) *
                         pow(densityOnGrid_spin[i], (5.0 / 3.0));

      double kineticEnergyDensity_TS = 0.0;
      double d_TS = 0.0;
      if (densityOnGrid_spin[i] > 1.0e-10) {
        // Calculation of the approximated kinetic energy density (index _TS, Tsirelson and Stash)
        kineticEnergyDensity_TS = (d_0 + (1.0 / 72.0) * (nablaDensityOnGrid2_spin[i] / densityOnGrid_spin[i]) +
                                   (1.0 / 6.0) * laplaceDensityOnGrid_spin[i]);
        // Calculation of d_TS
        d_TS = (kineticEnergyDensity_TS) - ((1.0 / 8.0) * (nablaDensityOnGrid2_spin[i] / densityOnGrid_spin[i]));
      }
      // Calculation of elf_TS
      if (d_0 > 1.0e-10)
        elf_TS_spin[i] = 1.0 / (1.0 + (d_TS * d_TS) / (d_0 * d_0));
    }
  };
  return elf_TS;
} /* calculateELFTSOnGrid */

template<>
GridData<RESTRICTED> ELFCalculator<RESTRICTED>::calculateTotalELFTSOnGrid(std::shared_ptr<GridController> gridController) {
  return this->calculateELFTSOnGrid(gridController);
}

template<>
GridData<RESTRICTED> ELFCalculator<UNRESTRICTED>::calculateTotalELFTSOnGrid(std::shared_ptr<GridController> gridController) {
  assert(gridController);

  // Definitions
  GridData<RESTRICTED> totalELF_TS(gridController);
  auto electronicStructure = _systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(_systemController->getSettings(), _basisController, gridController);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>>(
      basisFunctionOnGridController, _systemController->getSettings().grid.blockAveThreshold);
  const unsigned int nGridPoints = gridController->getNGridPoints();
  auto densityMatrix(electronicStructure->getDensityMatrix());
  auto densityOnGrid = densOnGridCalculator->calcDensityOnGrid(electronicStructure->getDensityMatrix());
  auto densityGradientOnGrid = makeGradient<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>>(gridController);
  auto densityHessianOnGrid = makeHessian<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>>(gridController);

  // Calculation of the density and its derivatives on a grid
  densOnGridCalculator->calcDensityAndDerivativesOnGrid(densityMatrix, densityOnGrid, densityGradientOnGrid,
                                                        densityHessianOnGrid);

  // Calculation of the x, y and z parts of the gradient and the laplacian of the electron density on a grid
  auto& densityGradientOnGridx = densityGradientOnGrid.x;
  auto& densityGradientOnGridy = densityGradientOnGrid.y;
  auto& densityGradientOnGridz = densityGradientOnGrid.z;
  auto& densityHessianOnGridxx = densityHessianOnGrid.xx;
  auto& densityHessianOnGridyy = densityHessianOnGrid.yy;
  auto& densityHessianOnGridzz = densityHessianOnGrid.zz;

  // Calculation of |grad(rho)|^2
  GridData<UNRESTRICTED> nablaDensityOnGrid2(gridController);
  for_spin(nablaDensityOnGrid2, densityGradientOnGridx, densityGradientOnGridy, densityGradientOnGridz) {
    nablaDensityOnGrid2_spin.array() = (densityGradientOnGridx_spin.array() * densityGradientOnGridx_spin.array()) +
                                       (densityGradientOnGridy_spin.array() * densityGradientOnGridy_spin.array()) +
                                       (densityGradientOnGridz_spin.array() * densityGradientOnGridz_spin.array());
  };

  // Calculation of laplace(roh)
  GridData<UNRESTRICTED> laplaceDensityOnGrid(gridController);
  for_spin(laplaceDensityOnGrid, densityHessianOnGridxx, densityHessianOnGridyy, densityHessianOnGridzz) {
    laplaceDensityOnGrid_spin.array() =
        densityHessianOnGridxx_spin.array() + densityHessianOnGridyy_spin.array() + densityHessianOnGridzz_spin.array();
  };

  // Loop over all grid points
  for (unsigned int i = 0; i < nGridPoints; ++i) {
    // Calculation of d_0(r)
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> d_0;
    for_spin(d_0, densityOnGrid) {
      d_0_spin = (3.0 / 10.0) * pow((6.0 * M_PI * M_PI), (2.0 / 3.0)) * pow(densityOnGrid_spin[i], (5.0 / 3.0));
    };

    // Calculation of the approximated kinetic energy density: (index _TS, Tsirelson and Stash)
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> kineticEnergyDensity_TS;
    for_spin(d_0, kineticEnergyDensity_TS, densityOnGrid, nablaDensityOnGrid2, laplaceDensityOnGrid) {
      if (densityOnGrid_spin[i] > 1.0e-10) {
        kineticEnergyDensity_TS_spin = (d_0_spin + (1.0 / 72.0) * (nablaDensityOnGrid2_spin[i] / densityOnGrid_spin[i]) +
                                        (1.0 / 6.0) * laplaceDensityOnGrid_spin[i]);
      }
    };
    // Calculation of elf
    if (densityOnGrid.alpha[i] > 1.0e-10 && densityOnGrid.beta[i] > 1.0e-10) {
      totalELF_TS[i] =
          1.0 / (1.0 + pow((kineticEnergyDensity_TS.alpha + kineticEnergyDensity_TS.beta -
                            (1.0 / 8.0) * (nablaDensityOnGrid2.alpha[i] / densityOnGrid.alpha[i]) -
                            (1.0 / 8.0) * (nablaDensityOnGrid2.beta[i] / densityOnGrid.beta[i])) /
                               ((3.0 / 10.0) * pow((6.0 * M_PI * M_PI), (2.0 / 3.0)) *
                                (pow(densityOnGrid.alpha[i], (5.0 / 3.0)) + pow(densityOnGrid.beta[i], (5.0 / 3.0)))),
                           2.0));
    }
  }
  return totalELF_TS;
} /* calculateTotalELFTSOnGrid */

template class ELFCalculator<Options::SCF_MODES::RESTRICTED>;
template class ELFCalculator<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
