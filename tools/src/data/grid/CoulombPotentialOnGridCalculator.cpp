/**
 * @file CoulombPotentialOnGridCalculator.cpp
 *
 * @date Mar 31, 2016
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
#include "data/grid/CoulombPotentialOnGridCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/SpinPolarizedData.h"
#include "data/grid/DensityOnGridController.h"
#include "data/grid/GridData.h"
#include "data/matrices/DensityMatrix.h"
#include "geometry/Atom.h"
#include "grid/Grid.h"
#include "grid/GridController.h"
#include "integrals/wrappers/Libint.h"
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <Eigen/Core>
#include <vector>

namespace Serenity {

template<Options::SCF_MODES SCF_MODE>
void CoulombPotentialOnGridCalculator::calculateElectronElectron(GridPotential<RESTRICTED>& result,
                                                                 const DensityMatrix<SCF_MODE>& densMat) {
  auto gridController = result.getGridController();
  auto nGridPoints = gridController->getNGridPoints();
  unsigned int nBlocks = omp_get_max_threads();
  auto basisController = densMat.getBasisController();
  const Eigen::MatrixXd totalDensity = densMat.total();
  auto basis = basisController->getBasis();
  const auto& gridPoints = gridController->getGridPoints();

  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::nuclear, 0, 2);
  Eigen::setNbThreads(1);
  // go through grid blockwise
#pragma omp for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    unsigned int n = (unsigned int)(nGridPoints / nBlocks);
    const unsigned int start = block * n;
    if (block == nBlocks - 1)
      n += nGridPoints % nBlocks;
    unsigned int blockEnd = start + n;
    for (unsigned int gridpoint = start; gridpoint < blockEnd; gridpoint++) {
      /*
       * integrals should be calculated to a point charge of -1
       * at the position of the current grid point
       */
      std::vector<std::pair<double, std::array<double, 3>>> point = {
          {-1.0, {{gridPoints(0, gridpoint), gridPoints(1, gridpoint), gridPoints(2, gridpoint)}}}};

      Eigen::MatrixXd ints = libint.compute1eInts(LIBINT_OPERATOR::nuclear, basisController, point);

      double pot = totalDensity.cwiseProduct(ints).sum();
      //... and store in the GridPotential
      result[gridpoint] += pot;
    } // gridPoint
  }   // block
  Eigen::setNbThreads(0);
}

void CoulombPotentialOnGridCalculator::calculateElectronNuclei(GridPotential<RESTRICTED>& result,
                                                               const std::vector<std::shared_ptr<Atom>>& atoms) {
  auto gridController = result.getGridController();
  auto nGridPoints = gridController->getNGridPoints();
  unsigned int nBlocks = omp_get_max_threads();
  const auto& gridPoints = gridController->getGridPoints();

  Eigen::setNbThreads(1);
  // go through grid blockwise
#pragma omp for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    unsigned int n = (unsigned int)(nGridPoints / nBlocks);
    const unsigned int start = block * n;
    if (block == nBlocks - 1)
      n += nGridPoints % nBlocks;
    unsigned int blockEnd = start + n;
    for (unsigned int gridpoint = start; gridpoint < blockEnd; gridpoint++) {
      // Eigen vector class for easy vector calculation
      Eigen::Vector3d gridCoord = gridPoints.col(gridpoint);
      // calculate sum_A [-charge_A/(|r-R_A|)] on every gridpoint
      for (auto atom : atoms) {
        // Eigen vector class for easy vector calculation
        Eigen::Vector3d atomCoord(atom->getX(), atom->getY(), atom->getZ());
        /*
         * precalculate the potential before storing it in the
         * spin polarized GridPotential (to prevent double calculation)
         */
        double tmpPot = -atom->getEffectiveCharge() / ((gridCoord - atomCoord).norm());
        // store
        result[gridpoint] += tmpPot;
      }
    } // gridPoint
  }   // block
  Eigen::setNbThreads(0);
}

template void CoulombPotentialOnGridCalculator::calculateElectronElectron<Options::SCF_MODES::RESTRICTED>(
    GridPotential<Options::SCF_MODES::RESTRICTED>& result, const DensityMatrix<Options::SCF_MODES::RESTRICTED>& densMat);
template void CoulombPotentialOnGridCalculator::calculateElectronElectron<Options::SCF_MODES::UNRESTRICTED>(
    GridPotential<Options::SCF_MODES::RESTRICTED>& result, const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& densMat);

} /* namespace Serenity */
