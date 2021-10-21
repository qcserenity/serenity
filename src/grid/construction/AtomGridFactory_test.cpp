/**
 * @file   AtomGridFactory_test.cpp
 *
 * @date   Mar 15, 2014
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
/* Include Serenity Internal Headers */
#include "grid/construction/AtomGridFactory.h"
#include "geometry/AtomType.h"
#include "grid/construction/AtomGrid.h"
#include "settings/GridOptions.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>

namespace Serenity {

/**
 * @class AtomGridFactoryTest
 * @brief Sets everything up for the tests of Matrix.h .
 */
class AtomGridFactoryTest : public ::testing::Test {};

/**
 * @test
 * @brief Tests radial integration with a simple set of Gaussian functions.
 *
 * Does a one-dimensional integration. The weights are expected to already
 * contain a factor r^2, which is needed for an integration in spherical
 * coordinates. This factor is here calculated out again.
 * Tested are different alphas and grid sized (nPoints)
 */
TEST_F(AtomGridFactoryTest, CheckRadialPointsByIntegration) {
  std::vector<Options::RADIAL_GRID_TYPES> radTypes = {
      Options::RADIAL_GRID_TYPES::AHLRICHS //,
                                           //    RADIAL_GRID_TYPES::BECKE//,
                                           //    RADIAL_GRID_TYPES::HANDY,
                                           //    RADIAL_GRID_TYPES::KNOWLES
  };
  for (unsigned int radI = 0; radI < radTypes.size(); ++radI) {
    for (double c = 5; c > 0.1; c *= 0.2) {
      for (unsigned int nRad = 3000; nRad <= 100000; nRad *= 6.37) {
        std::vector<double> radPoints(nRad);
        std::vector<double> radWeights(nRad);
        for (double alpha = 0.8; alpha <= 2.6; alpha += 0.8) {
          AtomGridFactory::_radialGrid(alpha, nRad, radPoints, radWeights, radTypes[radI]);
          double integral = 0;
          for (unsigned int i = 1; i < nRad; i++) {
            if (radPoints[i] != 0.0) {
              const double x = radPoints[i];
              integral += exp(-(x * x) / (2 * c * c)) * radWeights[i] / (x * x);
            }
          }
          EXPECT_NEAR(integral, 0.5 * sqrt(2.0 * M_PI) * c, 1e-6);
        }
      }
    }
  }
}
/**
 * @test
 * @brief Makes sure the spherical grid integrates to one on each accuracy level
 */
TEST_F(AtomGridFactoryTest, CheckSphericalPointsByIntegration) {
  for (double i = 0; i < 20; i++) {
    unsigned int nPoints;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> weights;
    AtomGridFactory::_lebedevSphericalGrid(i, nPoints, x, y, z, weights);
    ASSERT_EQ(nPoints, x.size());
    ASSERT_EQ(nPoints, y.size());
    ASSERT_EQ(nPoints, z.size());
    ASSERT_EQ(nPoints, weights.size());
    double integral = 0;
    for (unsigned int j = 0; j < nPoints; j++) {
      integral += weights[j];
      EXPECT_NEAR(x[i] * x[i] + y[i] * y[i] + z[i] * z[i], 1.0, 1e-8);
    }
    EXPECT_NEAR(integral, 1.0, 1e-8);
  }
}

/*
double sphereFunction(const Point& point) {
  if (point.distanceToOrigin() <= 1.0) {
    return 1.0;
  } else {
    return 0.0;
  }
}

TEST_F(AtomGridFactoryTest, IntegrateSphereAroundDummyAtom) {
  AtomType dummyType("C", 1, 1.0, 1.0);

  const vector<unsigned int> allNRad = {
      1000
  };
  const vector<vector<unsigned int> > allIncrSphAcc = {
    {100,100,100,100,100,100,100,100}
  };

  for (unsigned int testIndex=0; testIndex < allNRad.size(); testIndex++) {
    const AtomGrid dummyAtomGrid = tAtomGridFactory.produce(
        &dummyType, allNRad[testIndex], allIncrSphAcc[testIndex]);

    const vector<Point> points = dummyAtomGrid.getGridPoints();
    const vector<double> weights = dummyAtomGrid.getWeights();
    double volume = 0.0;
    double volumeShifted = 0.0;
    Point shift(1.0, 0.0, 0.0);
    for (unsigned int i=0; i<points.size(); i++) {
      volume += sphereFunction(points[i]) * weights[i];
      volumeShifted += sphereFunction(points[i]+shift) * weights[i];
    }
    EXPECT_NEAR(volume, 4.0/3.0*M_PI, 1e-3);
    EXPECT_NEAR(volumeShifted, 4.0/3.0*M_PI, 0.4);
  }
}
*/
} // namespace Serenity
