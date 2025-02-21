/**
 * @file   AtomGridFactory.cpp
 *
 * @date   Mar 14, 2014
 * @author Dennis Barton
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
#include "grid/construction/AtomGridFactory.h"
/* Include Serenity Internal Headers */
#include "geometry/AtomType.h"
#include "geometry/Point.h"
#include "grid/construction/AtomGrid.h"
#include "grid/construction/sphere_lebedev_rule.h"
#include "math/FloatMaths.h"
#include "misc/SerenityError.h"
#include "settings/GridOptions.h"
/* Include Std and External Headers */
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace Serenity {
/*
 * Initialize static members
 */
const std::array<double, N_ELEMENTS_IN_PERIODIC_TABLE + 1> AtomGridFactory::_ahlrichsAlphaValues = {
    {1.1, // Dummy entry to fit the indices to the nuclear charges starting at 1
     0.8, 0.9, 1.8, 1.4, 1.3, 1.1, 0.9, 0.9, 0.9, 0.9, 1.4, 1.3, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0, 1.5, 1.4, 1.3, 1.2, 1.2,
     1.2, 1.2, 1.2, 1.2, 1.1, 1.1, 1.1, 1.1, 1.0, 0.9, 0.9, 0.9, 0.9,
     /*
      * The original paper [J. Chem. Phys. (1995) 102, 346] only has values until Krypton
      * (i.e. until here). However, the authors
      * also state that there are "essential changes for the hydrogen atom only" upon optimization
      * of the alpha value. Thus, we use the numbers given for the fourth period also for the
      * following ones. For the lanthanoids and actinoids we use the unoptimized value of 1.0.
      */
     1.5, 1.4, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.1, 1.1, 1.1, 1.1, 1.0, 0.9, 0.9, 0.9, 0.9, 1.5, 1.4, 1.3, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.1, 1.1, 1.1, 1.1, 1.0,
     0.9, 0.9, 0.9, 0.9, 1.5, 1.4, 1.3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.2, 1.2,
     1.2, 1.2, 1.2, 1.2, 1.1, 1.1, 1.1, 1.1, 1.0, 0.9, 0.9, 0.9, 0.9}};
const std::array<double, N_ELEMENTS_IN_PERIODIC_TABLE + 1> AtomGridFactory::_clementisRadii = {
    {1.00, // Dummy atoms will create a H like grid
     1.00, 0.59, 3.16, 2.12, 1.64, 1.27, 1.06, 0.91, 0.79, 0.72, 3.59, 2.74, 2.23, 2.10, 1.85, 1.66, 1.49, 1.34, 4.59,
     3.67, 3.48, 3.33, 3.23, 3.14, 3.04, 2.95, 2.87, 2.82, 2.74, 2.68, 2.57, 2.36, 2.15, 1.95, 1.78, 1.66, 5.01, 4.14,
     4.01, 3.89, 3.74, 3.59, 3.46, 3.36, 3.27, 3.19, 3.12, 3.04, 2.95, 2.74, 2.51, 2.32, 2.17, 2.04, 5.63, 4.78, 3.68,
     3.50, 4.67, 3.89, 3.87, 4.50, 4.37, 4.40, 4.25, 4.31, 4.27, 4.27, 4.20, 4.20, 4.10, 3.93, 3.78, 3.65, 3.55, 3.50,
     3.40, 3.34, 3.29, 3.23, 2.95, 2.91, 2.70, 2.55, 2.40, 2.27, 6.00, 5.70, 3.68, // Guessed from
                                                                                   // here on
     3.40, 3.40, 3.31, 3.31, 3.31, 3.31, 3.31, 3.31, 3.31, 3.31, 3.31, 3.31, 3.31, 3.31, 4.10, 4.10, 4.10, 4.10, 4.10,
     4.10, 4.10, 4.10, 4.10, 3.00, 2.90, 2.80, 2.70, 2.60, 2.50}};

std::unique_ptr<AtomGridFactory> AtomGridFactory::_instance;

std::shared_ptr<const AtomGrid> AtomGridFactory::produce(Options::RADIAL_GRID_TYPES radType,
                                                         Options::SPHERICAL_GRID_TYPES sphType,
                                                         const std::shared_ptr<const AtomType> atomType, unsigned int acc) {
  if (!_instance)
    _instance.reset(new AtomGridFactory());
  return _instance->getOrProduce(radType, sphType, atomType, acc);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wswitch"

std::unique_ptr<const AtomGrid> AtomGridFactory::produceNew(const Options::RADIAL_GRID_TYPES radialType,
                                                            const Options::SPHERICAL_GRID_TYPES sphericalType,
                                                            const std::shared_ptr<const AtomType> atomType,
                                                            const unsigned int acc) {
  /*===============
   *  Radial Part
   *===============*/
  if (!(0 < acc and acc < 8))
    throw SerenityError("Invalid accuracy level chosen, must be between [1,7].");
  constexpr std::array<int, 7> radAcc = {{13, 13, 13, 14, 15, 16, 17}};
  const unsigned int nRadial = (unsigned int)(5.0 * (radAcc[acc - 1] + atomType->getRow() - 8));

  double alpha = 0.0;
  std::vector<double> points;
  std::vector<double> weights;
  switch (radialType) {
    case Options::RADIAL_GRID_TYPES::BECKE:
      // Becke uses Bragg-Slater radii as
      //   parameter to optimize radial integration
      if (atomType->getElementSymbol() == "H") {
        alpha = atomType->getBraggSlaterRadius();
      }
      else {
        alpha = 0.5 * atomType->getBraggSlaterRadius();
      }
      break;
    case Options::RADIAL_GRID_TYPES::HANDY:
      /*
       *  We will use the scheme Handy describes.
       *  Because of reasons stated at the declaration of
       *  the Ahlrichs alpha values above, we will use
       *  the Bragg-Slater radii here.
       */

      if (atomType->getElementSymbol() == "H") {
        alpha = atomType->getBraggSlaterRadius();
      }
      else {
        alpha = 0.5 * atomType->getBraggSlaterRadius();
      }
      break;
    case Options::RADIAL_GRID_TYPES::AHLRICHS:
      alpha = _ahlrichsAlphaValues[atomType->getPSEPosition()];
      assert(alpha > 0.0);
      break;
    case Options::RADIAL_GRID_TYPES::KNOWLES:
      if (atomType->getElementSymbol() == "H" or atomType->getElementSymbol() == "He" or
          atomType->getElementSymbol() == "Li" or atomType->getElementSymbol() == "Be" or
          atomType->getElementSymbol() == "Na" or atomType->getElementSymbol() == "Mg" or
          atomType->getElementSymbol() == "K" or atomType->getElementSymbol() == "Ca") {
        alpha = 5.0;
      }
      else {
        alpha = 7.0;
      }
      break;
    case Options::RADIAL_GRID_TYPES::EQUI:
      // No alpha needed
      alpha = 0.0;
      break;
      // No default --> compiler should throw an error if a case is not covered!
  }
  // Call radial grid function
  // type will be determined internally
  std::vector<double> radPoints(nRadial);
  std::vector<double> radWeights(nRadial);
  _radialGrid(alpha, nRadial, radPoints, radWeights, radialType);

  /*================
   *  Angular Part
   *================*/
  switch (sphericalType) {
    // For now the only one we have
    // Lebedev
    case Options::SPHERICAL_GRID_TYPES::LEBEDEV:

      // Lebedev rules according to accuracy,
      //  following the data in the ORCA manual (4.0) p291 ff.
      constexpr std::array<std::array<unsigned int, 5>, 7> ldval = {{{{4, 4, 4, 4, 4}},
                                                                     {{4, 4, 4, 7, 4}},
                                                                     {{4, 4, 7, 10, 7}},
                                                                     {{4, 7, 10, 13, 10}},
                                                                     {{7, 10, 13, 15, 13}},
                                                                     {{10, 13, 15, 16, 15}},
                                                                     {{13, 15, 16, 17, 16}}}};

      // definition of pruning zones according to SG1 grid
      //  used with Clementi radii following the description
      //  in the ORCA manual (4.0) p291 ff.
      constexpr std::array<std::array<double, 4>, 3> ranges = {
          {{{0.25, 0.5, 1.0, 4.5}}, {{0.1667, 0.5, 0.9, 3.5}}, {{0.1, 0.4, 0.8, 2.5}}}};

      unsigned int sphericalAcc = 0;

      const unsigned int row = (atomType->getRow() > 3) ? 2 : atomType->getRow() - 1;
      const unsigned int redp1 = (atomType->getRow() == 1 and acc > 1) ? 2 : 1;
      std::array<double, 4> r = ranges[row];
      for (auto& i : r) {
        i *= _clementisRadii[atomType->getPSEPosition()];
      }
      for (unsigned int i = 0; i < nRadial; ++i) {
        unsigned int sphN;
        std::vector<double> sphX;
        std::vector<double> sphY;
        std::vector<double> sphZ;
        std::vector<double> sphW;
        if (sphericalAcc < 4)
          if (radPoints[i] > r[sphericalAcc])
            ++sphericalAcc;
        _lebedevSphericalGrid(ldval[acc - redp1][sphericalAcc], sphN, sphX, sphY, sphZ, sphW);

        for (unsigned int k = 0; k < sphN; k++) {
          weights.push_back(sphW[k] * radWeights[i] * 4.0 * M_PI);
          points.push_back(sphX[k] * radPoints[i]);
          points.push_back(sphY[k] * radPoints[i]);
          points.push_back(sphZ[k] * radPoints[i]);
        }
      }
      break;
  }
  Eigen::VectorXd weights_e = Eigen::VectorXd::Map(weights.data(), weights.size());
  Eigen::Matrix3Xd points_e = Eigen::MatrixXd::Map(points.data(), 3, int(points.size() / 3));
  return std::unique_ptr<AtomGrid>(new AtomGrid(points_e, weights_e));
}

#pragma GCC diagnostic pop // -Wswitch

/*
 * Function for the creation of a radial grid
 */
void AtomGridFactory::_radialGrid(double alpha, unsigned int nRadial, std::vector<double>& radPoints,
                                  std::vector<double>& radWeights, const Options::RADIAL_GRID_TYPES& radType) {
  assert(radPoints.size() == nRadial);
  assert(radWeights.size() == nRadial);

  switch (radType) {
    case Options::RADIAL_GRID_TYPES::BECKE:
      /*
       * BECKE grid with Chebyshev quadrature of the second kind according to
       *     ref.: J. Comput. Chem. (2003) 24 p732-740
       *
       */
      for (unsigned int i = 1; i <= nRadial; i++) {
        double xi = cos((double)i * M_PI / (nRadial + 1));
        radWeights[nRadial - i] =
            sqrt(fipow((1.0 + xi), 5) / fipow((1.0 - xi), 7)) * (2.0 * M_PI) * fipow(alpha, 3) / (nRadial + 1);
        radPoints[nRadial - i] = alpha * (1.0 + xi) / (1.0 - xi);
      }
      break;

    case Options::RADIAL_GRID_TYPES::HANDY:
      /*
       *  HANDY grid according to
       *     ref.: J. Comput. Chem. (2003) 24 p732-740
       */
      for (unsigned int i = 1; i <= nRadial; i++) {
        double xi = i / (nRadial + 1.0);
        radPoints[i - 1] = alpha * xi * xi / ((1 - xi) * (1 - xi));
        radWeights[i - 1] = 2 * pow(alpha, 3) * pow(xi, 5) / (nRadial + 1) / pow((1 - xi), 7);
      }
      break;

    case Options::RADIAL_GRID_TYPES::AHLRICHS:
      /*
       * AHLRICHS grid with Chebyshev quadrature of the second kind according to
       *     ref.: J. Chem. Phys. (1998) 108, 3226
       *
       *     Original paper:
       *     ref.:  J. Chem. Phys. (1995) 102, 346
       */
      for (unsigned int i = 1; i <= nRadial; i++) {
        // M3 with m (called alpha in Ahlrichs paper) = 0.6
        const double log2 = log(2.0);
        const double tmp = alpha / log2;
        const double xi = cos(i * M_PI / (nRadial + 1.0));
        const double ri = tmp * pow((xi + 1.0), 0.6) * log(2.0 / (1.0 - xi));
        const double ln = log((1.0 - xi) / 2.0);
        const double sq = sqrt((1.0 + xi) / (1.0 - xi));
        const double wi = (M_PI / (nRadial + 1.0)) * pow((1.0 + xi), (3 * 0.6)) * tmp * tmp * tmp *
                          (sq * ln * ln - 0.6 * ln * ln * ln / sq);
        radWeights[nRadial - i] = wi;
        radPoints[nRadial - i] = ri;
      }
      break;

    case Options::RADIAL_GRID_TYPES::KNOWLES:
      /*
       * KNOWLES grid according to
       *     ref.: J. Comput. Chem. (2003) 24 p732-740
       */
      for (unsigned int i = 1; i <= nRadial; i++) {
        double xi = i / (nRadial + 1.0);
        double tmp = 1 - xi * xi * xi;
        double ln = log(tmp);
        radPoints[nRadial - i] = -alpha * ln;
        radWeights[nRadial - i] = 3 * xi * xi * ln * ln * alpha * alpha * alpha / ((nRadial + 1) * tmp);
      }
      break;

    case Options::RADIAL_GRID_TYPES::EQUI:
      /*
       * EQUIDISTANT grid
       *     for test purposes.
       */
      for (unsigned int i = 1; i <= nRadial; i++) {
        radPoints[i - 1] = i * 5.0 / nRadial;
        radWeights[i - 1] = 1.0 / nRadial;
      }
      break;
  }
  return;
}

// function for spherical grid
void AtomGridFactory::_lebedevSphericalGrid(unsigned int i, unsigned int& sphN, std::vector<double>& sphX,
                                            std::vector<double>& sphY, std::vector<double>& sphZ, std::vector<double>& sphW) {
  switch (i) {
    case 0:
      sphN = 6;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0006(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 1:
      sphN = 14;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0014(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 2:
      sphN = 26;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0026(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 3:
      sphN = 38;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0038(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 4:
      sphN = 50;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0050(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 5:
      sphN = 74;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0074(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 6:
      sphN = 86;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0086(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 7:
      sphN = 110;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0110(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 8:
      sphN = 146;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0146(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 9:
      sphN = 170;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0170(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 10:
      sphN = 194;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0194(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 11:
      sphN = 230;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0230(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 12:
      sphN = 266;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0266(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 13:
      sphN = 302;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0302(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 14:
      sphN = 350;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0350(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 15:
      sphN = 434;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0434(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 16:
      sphN = 590;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0590(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 17:
      sphN = 770;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0770(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 18:
      sphN = 974;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld0974(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 19:
      sphN = 1202;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld1202(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    case 20:
      sphN = 1454;
      sphX.resize(sphN);
      sphY.resize(sphN);
      sphZ.resize(sphN);
      sphW.resize(sphN);
      ld1454(sphX.data(), sphY.data(), sphZ.data(), sphW.data());
      break;
    default:
      throw SerenityError("Spherical grid (Lebedev) with unsupported accuracy was requested.");
  }
  return;
}
} /* namespace Serenity */
