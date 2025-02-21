/**
 * @file   GridFactory.cpp
 *
 * @date   Mar 14, 2014
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
#include "grid/construction/GridFactory.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "geometry/AtomType.h"
#include "geometry/Geometry.h"
#include "geometry/Point.h"
#include "grid/construction/AtomGrid.h"
#include "grid/construction/AtomGridFactory.h"
#include "io/FormattedOutput.h"
#include "io/IOOptions.h"
#include "math/Matrix.h"
#include "misc/Timing.h"
#include "settings/GridOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <vector>

namespace Serenity {

GridFactory::GridFactory(Options::GRID_TYPES flavour, unsigned int smoothing, const Options::RADIAL_GRID_TYPES radialType,
                         const Options::SPHERICAL_GRID_TYPES sphericalType, unsigned int accuracyLevel, double weightThreshold)
  : _flavour(flavour),
    _smoothing(smoothing),
    _radialType(radialType),
    _sphericalType(sphericalType),
    _accuracyLevel(accuracyLevel),
    _weightThreshold(weightThreshold) {
}
/*
 * This implementation is according to: A.D. Becke, J.Chem.Phys. 88 (1988), 2547.
 */
std::unique_ptr<AtomCenteredGrid> GridFactory::produce(std::shared_ptr<const Geometry> geometry) {
  Timings::takeTime("Tech. -     Grid Constructions");
  assert(geometry);

  if (geometry->hasIdenticalAtoms()) {
    throw SerenityError("Grid construction with duplicated atoms was detected!\n\
                         You have to delete all identical atoms from the geometry before\n\
                         entering the grid construction!");
  }

  takeTime("Setup Grid");
  if (iOOptions.printGridInfo)
    printSmallCaption("Constructing an integration grid (accuracy: " + std::to_string(_accuracyLevel) + ")");

  /* ===================================
   *   Tabulate Inter Atomic Distances
   * ===================================*/
  const auto& atoms = geometry->getAtoms();
  auto coords = geometry->getCoordinates();
  const unsigned int nAtoms = atoms.size();
  /* Precompute atomic distances */
  Matrix<double> atomicDistances(nAtoms, nAtoms);
  auto atomDistPtr = atomicDistances.data();
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      atomDistPtr[i + nAtoms * j] = distance(*atoms[i], *atoms[j]);
      atomDistPtr[j + nAtoms * i] = atomDistPtr[i + nAtoms * j];
    }
  }

  /* =============================================
   *   Calculate Size Adjustment For Fuzzy Cells
   * =============================================*/
  /* This parameter is explained in the appendix of the reference. It modifies the atom sizes. */
  Eigen::MatrixXd aOfAtomPair = Eigen::MatrixXd::Zero(nAtoms, nAtoms);
  if (_flavour == Options::GRID_TYPES::BECKE) {
    /*
     * Calculate the atomic size adjustments parameter for each atom pair. (
     * See the appendix of the reference.)
     */
#pragma omp parallel for ordered schedule(static)
    for (unsigned int i = 0; i < nAtoms; ++i) {
      const auto& otherAtom = atoms[i];
      for (unsigned int j = 0; j < nAtoms; ++j) {
        const auto& currentAtom = atoms[j];
        double braggSlaterQuot =
            sqrt(currentAtom->getAtomType()->getBraggSlaterRadius() / otherAtom->getAtomType()->getBraggSlaterRadius());
        double u = (braggSlaterQuot - 1) / (braggSlaterQuot + 1);
        aOfAtomPair(j, i) = u / (u * u - 1);
        if (aOfAtomPair(j, i) < -0.5) {
          aOfAtomPair(j, i) = -0.5;
        }
        else if (aOfAtomPair(j, i) > 0.5) {
          aOfAtomPair(j, i) = 0.5;
        }
      }
    }
  }
  auto a = aOfAtomPair.data();

  /* ===================================
   *   Construct Initial Grid Per Atom
   * ===================================*/
  // First atom explicitly, to have the Singletons initialized
  std::vector<std::shared_ptr<const AtomGrid>> atomGrids(nAtoms, nullptr);
  const auto& currentGridAtm = atoms[0];
  /* get the reference grid of the current atom type */
  atomGrids[0] = AtomGridFactory::produce(_radialType, _sphericalType, currentGridAtm->getAtomType(), _accuracyLevel);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int k = 1; k < nAtoms; ++k) {
    const auto& currentGridAtom = atoms[k];
    /* get the reference grid of the current atom type */
    atomGrids[k] = AtomGridFactory::produce(_radialType, _sphericalType, currentGridAtom->getAtomType(), _accuracyLevel);
  }

  /* =====================================
   *   Construct/Cut Final Grid Per Atom
   * =====================================*/
  std::vector<std::shared_ptr<std::vector<double>>> atomGridPoints(nAtoms, nullptr);
  std::vector<std::shared_ptr<std::vector<double>>> atomWeights(nAtoms, nullptr);
  // parallelize over atoms
#pragma omp parallel for schedule(dynamic)
  for (unsigned int k = 0; k < nAtoms; ++k) {
    atomGridPoints[k] = std::make_shared<std::vector<double>>();
    atomWeights[k] = std::make_shared<std::vector<double>>();
    auto pts = atomGridPoints[k];
    auto w = atomWeights[k];

    std::vector<unsigned int> significantAtoms;
    double minAtomDist = 999999999.9;
    for (unsigned int l = 0; l < nAtoms; ++l) {
      if (atomDistPtr[l + nAtoms * k] < 40.0)
        significantAtoms.push_back(l);
      if (l != k)
        minAtomDist = std::min(minAtomDist, atomDistPtr[l + nAtoms * k]);
    }
    /* get the reference grid of the current atom type */
    const auto& currentAtomGrid = atomGrids[k];
    /* reassign reference gridPoints and reference weights of the current atom to usable variables */
    const auto& refGridPoints = currentAtomGrid->getGridPoints();
    const auto& refWeights = currentAtomGrid->getWeights();
    // Reserve space to efficiently push_back data
    atomWeights[k]->reserve(currentAtomGrid->nPoints());
    atomGridPoints[k]->reserve(currentAtomGrid->nPoints() * 3);

    /*
     * Loop over all grid points in currentAtomGrid
     */
    for (unsigned int i = 0; i < currentAtomGrid->nPoints(); i++) {
      /* get a grid point and shift it according to the atom coordinates */
      Eigen::Vector3d gridPoint(refGridPoints.col(i));
      gridPoint += coords.row(k);
      /* copy the reference weight in order to manipulate it */
      double weight(refWeights[i]);
      /* reset cellFunctionSum */
      double cellFunctionSum = 0.0;

      /* Precalculate distances between atoms and this gridPoint */
      Eigen::VectorXd distancesAtomP((coords.rowwise() - gridPoint.transpose()).rowwise().norm());
      auto rdist = distancesAtomP.data();

      /*
       * Calculate the weight of this grid point in the final grid as:
       *   w(r) = P_k(r) / \sum_l P_l(r),
       * where P_K is the cell function on atom k
       *
       * all P_l are calculated in one double loop for all nu_lj
       */
      if (_flavour == Options::GRID_TYPES::BECKE) {
        for (auto& l : significantAtoms) {
          /* reset cellFunction */
          double cellFunction = 1.0;
          for (auto& j : significantAtoms) {
            if (l == j)
              continue;
            /* Eq. (11) and Eq. (A2)*/
            const double mu = (rdist[l] - rdist[j]) / atomicDistances(l, j);
            const double nu = mu + a[j + nAtoms * l] * (1.0 - (mu * mu));
            /* Smooth the grid to make a difference to the strict Voronoi cut */
            cellFunction *= beckeSmoothFunction(nu, _smoothing);
          }
          if (l == k)
            weight *= cellFunction;
          cellFunctionSum += cellFunction;
        }
        weight /= cellFunctionSum;
      }
      else if (_flavour == Options::GRID_TYPES::VORONOI) {
        for (auto& l : significantAtoms) {
          const double mu = (rdist[l] - rdist[k]) / atomDistPtr[l + nAtoms * k];
          const double nu = mu + a[k + nAtoms * l] * (1.0 - (mu * mu));
          if (nu < 0.0) {
            weight = 0.0;
            break;
          }
        }
      }
      else if (_flavour == Options::GRID_TYPES::SSF) {
        /*
         * Reference:
         * Stratmann, Scuseria and Frisch, Chem. Phys. Lett 257, (1996) 213-223
         */

        // Screen according to Eq. (15)
        if (rdist[k] >= 0.5 * (1.0 - 0.64) * minAtomDist) {
          // Screen if nu_kj>0.64 => w_k(r)=0 => p_k(r)=0
          bool done = false;
          for (unsigned int j = 0; j < nAtoms; j++) {
            double nuij = (rdist[k] - rdist[j]) / atomDistPtr[k + nAtoms * j];
            if (nuij >= 0.64) {
              weight = 0.0;
              done = true;
              break;
            }
          }
          if (done)
            continue;

          for (unsigned int i = 0; i < nAtoms; i++) {
            /* reset cellFunction */
            double cellFunction = 1.0;
            for (unsigned int j = 0; j < nAtoms; j++) {
              if (i == j)
                continue;
              // Eq. (4)
              double nuij = (rdist[i] - rdist[j]) / atomDistPtr[i + nAtoms * j];
              if (nuij <= -0.64) {
                // Screen if nu_j<-0.64 => s(nuij)=1
                continue;
              }
              else if (nuij >= 0.64) {
                // Screen if nu_kj>0.64 => s(nuij)=0 => w_k(r)=0
                cellFunction = 0.0;
                break;
              }
              else {
                // Solve the default polynomial
                nuij /= 0.64;
                const double n3 = nuij * nuij * nuij;
                const double poly = (-5.0 * n3 * n3 * nuij + 21.0 * n3 * nuij * nuij - 35.0 * n3 + 35.0 * nuij) / 16.0;
                cellFunction *= 0.5 * (1 - poly);
              }
            }
            if (i == k)
              weight *= cellFunction;
            cellFunctionSum += cellFunction;
          }
          weight /= cellFunctionSum;
        }
      }
      //      assert(weight >= 0.0);
      /* Check for significance of this grid point. Throw away if insignificant. */
      if (weight > _weightThreshold) {
        /* This grid point is significant. Store. */
        w->push_back(weight);
        pts->push_back(gridPoint[0]);
        pts->push_back(gridPoint[1]);
        pts->push_back(gridPoint[2]);
      }
    } /*loop points*/
  }   /*loop atoms*/

  /* ============================
   *   Set indices needed later
   * ============================*/
  std::vector<std::pair<unsigned int, unsigned int>> gridIndicesOfAtoms(nAtoms);
  int sum = 0;
  unsigned int nGridPointsTot = atomGrids[nAtoms - 1]->nPoints();
  for (unsigned int i = 0; i < nAtoms - 1; i++) {
    nGridPointsTot += atomGrids[i]->nPoints();
    sum += atomWeights[i]->size();
    gridIndicesOfAtoms[i].second = sum;
    gridIndicesOfAtoms[i + 1].first = sum;
  }
  unsigned int nGridPointsRelevant = gridIndicesOfAtoms[nAtoms - 1].first + atomWeights[nAtoms - 1]->size();
  gridIndicesOfAtoms[nAtoms - 1].second = nGridPointsRelevant;
  gridIndicesOfAtoms[0].first = 0;

  /* ============================
   *   Merge atomwise grid data
   * ============================*/
  auto gridPoints = std::unique_ptr<Eigen::Matrix3Xd>(new Eigen::Matrix3Xd(3, nGridPointsRelevant));
  auto weights = std::unique_ptr<Eigen::VectorXd>(new Eigen::VectorXd(nGridPointsRelevant));

#pragma omp parallel for schedule(dynamic)
  for (unsigned int k = 0; k < nAtoms; ++k) {
    unsigned int n = gridIndicesOfAtoms[k].second - gridIndicesOfAtoms[k].first;
    gridPoints->block(0, gridIndicesOfAtoms[k].first, 3, n) = Eigen::MatrixXd::Map(atomGridPoints[k]->data(), 3, n);
    weights->segment(gridIndicesOfAtoms[k].first, n) = Eigen::VectorXd::Map(atomWeights[k]->data(), n);
    atomGridPoints[k] = nullptr;
    atomWeights[k] = nullptr;
  }

  /* ===================================
   *   Sanity Check, Output and Return
   * ===================================*/
  assert(nGridPointsRelevant == gridPoints->cols());
  assert(nGridPointsRelevant == weights->size());

  if (iOOptions.printGridInfo) {
    print((std::string) "Total number of grid points (before weight cutoff):    " + nGridPointsTot);
    print((std::string) "Total number of grid points (after weight cutoff):     " + nGridPointsRelevant);
    print("");
  }
  timeTaken(2, "Setup Grid");

  auto grid =
      std::unique_ptr<AtomCenteredGrid>(new AtomCenteredGrid(std::move(gridPoints), std::move(weights), gridIndicesOfAtoms));
  Timings::timeTaken("Tech. -     Grid Constructions");
  return grid;
}
/* Eq. (19) and (20) */
constexpr double GridFactory::pRecursive(double nu, unsigned int k) {
  /* Slow recursive implementation
  if (k > 1) {
    double tmp;
    tmp = pRecursive(nu, k - 1);
    return (3.0 * tmp - tmp * tmp * tmp) / 2.0;
  }
  return (3.0 * nu - nu * nu * nu) / 2.0;
   */
  /* Alternative, WAY more efficient than the recursive implementation!
  for (unsigned int i=k; i>1; i--) {
    nu = (3.0 * nu - nu * nu * nu) / 2.0;
  }
  return (3.0 * nu - nu * nu * nu) / 2.0;
   */
  /* Functional programming style to make the function constexpr -> pure. Even more efficient. */
  return (k > 1) ? pRecursive((3.0 - nu * nu) * nu / 2.0, k - 1) : (3.0 - nu * nu) * nu / 2.0;
}
/* Eq. (21) */
constexpr double GridFactory::beckeSmoothFunction(double nu, unsigned int k) {
  return 0.5 * (1 - pRecursive(nu, k));
}

} /* namespace Serenity */
