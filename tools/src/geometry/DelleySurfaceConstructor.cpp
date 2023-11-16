/**
 * @file DelleySurfaceConstructor.cpp
 *
 * @date   May 27, 2020
 * @author Moritz Bensberg
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
#include "geometry/DelleySurfaceConstructor.h"
/* Include Serenity Internal Headers */
#include "geometry/Ellipse.h"                  //Ellipse definition.
#include "geometry/Line.h"                     //Line definition.
#include "geometry/MolecularSurface.h"         //The final result: The molecular surface.
#include "geometry/Plane.h"                    //Plane definition.
#include "geometry/Sphere.h"                   //Sphere definition.
#include "grid/GridController.h"               //Surface grid controller.
#include "grid/construction/AtomGridFactory.h" //Lebedev spherical grid.
#include "math/FloatMaths.h"                   //isEqual();
#include "misc/SerenityError.h"                //throw errors
#include "misc/Timing.h"                       //Timings.
/* Include Std and External Headers */
#include <cmath>

namespace Serenity {

DelleySurfaceConstructor::DelleySurfaceConstructor(std::vector<Sphere> spheres, double r_s, double alpha, double projectionCutOff,
                                                   double minimalDistance, bool oneCavity, double connectivityFactor)
  : _spheres(spheres),
    _r_s(r_s),
    _alpha(alpha),
    _projectionCutOff(projectionCutOff),
    _minimalDistance(minimalDistance),
    _oneCavity(oneCavity),
    _connectivityFactor(connectivityFactor),
    _cutOff(4.0 / _alpha) {
  double eM4 = exp(-4);
  _a = 6.0 / 4.0 * _alpha * eM4;
  _b = -5.0 / 16.0 * eM4 * _alpha * _alpha;
  for (const auto& sphere : _spheres)
    getUnitSphere(sphere.getAngularMomentum());
  calculateCylinderRadiiAndParameters();
}

DelleySurfaceConstructor::~DelleySurfaceConstructor() = default;

Plane DelleySurfaceConstructor::getSphereSphereTouchingPlane(const Sphere& sphere_i, const Sphere& sphere_j) {
  const double totalWeight = sphere_i.getRadius() + sphere_j.getRadius();
  const Eigen::Vector3d r_i = sphere_i.getCenterCoords();
  const Eigen::Vector3d r_j = sphere_j.getCenterCoords();
  const Eigen::Vector3d weightedCenter = sphere_i.getRadius() / totalWeight * r_i + sphere_j.getRadius() / totalWeight * r_j;
  const Eigen::Vector3d normalVector = r_i - r_j;
  return Plane(weightedCenter, normalVector);
}

std::unique_ptr<MolecularSurface> DelleySurfaceConstructor::getMolecularSurface() {
  /*
   * The surface construction is done sphere-wise.
   * For every sphere the points are projected from the unit-sphere
   * to the molecular surface. The weights are set to the weights
   * that correspond to a sphere with the radius being the distance
   * to the original center.
   * Weights are adjusted if they are on the bond part between spheres
   * and are constrained to remain on the side of the original sphere.
   * If a point falls into the area close to the plane orthogonal to
   * the bond, the point is interpreted as an ellipse and it is checked
   * if the plane cuts the ellipse. If so, the weight of the point is
   * adjusted by ratio that is cut and the point is moved to the new
   * center of gravity of the ellipse fragment.
   */
  std::vector<std::vector<Eigen::Vector3d>> allFinalCoordinates;
  std::vector<std::vector<Eigen::Vector3d>> allFinalNorms;
  std::vector<std::vector<double>> allFinalWeights;
  std::vector<std::vector<unsigned int>> allSphereIndices;
  Eigen::setNbThreads(1);
  unsigned int nThreads = 1;
#ifdef _OPENMP
  nThreads = omp_get_max_threads();
#endif
  for (unsigned int iThread = 0; iThread < nThreads; ++iThread) {
    std::vector<Eigen::Vector3d> newCoordVector = {};
    std::vector<double> newWeightVector = {};
    std::vector<unsigned int> newIndexVector = {};
    allFinalCoordinates.push_back(newCoordVector);
    allFinalNorms.push_back(newCoordVector);
    allFinalWeights.push_back(newWeightVector);
    allSphereIndices.push_back(newIndexVector);
  }
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iSphere = 0; iSphere < _spheres.size(); ++iSphere) {
    unsigned int threadId = 0;
#ifdef _OPENMP
    threadId = omp_get_thread_num();
#endif
    auto& finalCoordinates = allFinalCoordinates[threadId];
    auto& finalNorms = allFinalNorms[threadId];
    auto& finalWeights = allFinalWeights[threadId];
    auto& finalSphereIndices = allSphereIndices[threadId];

    const Sphere& sphere = _spheres[iSphere];
    const UnitSphere& unitSphere = getUnitSphere(sphere.getAngularMomentum());
    const Eigen::Vector3d r_i = sphere.getCenterCoords();
    const double R_i = sphere.getRadius();
    const Eigen::Vector3d center = sphere.getCenterCoords();
    Eigen::Matrix3Xd sCoords = unitSphere.coordinates;
    const Eigen::VectorXd& sWeights = unitSphere.weights;
    sCoords.colwise() += center;
    const unsigned int nPoints = sCoords.cols();
    for (unsigned int i = 0; i < nPoints; ++i) {
      Eigen::Vector3d r, grad;
      if (nThreads == 1)
        Timings::takeTime("Tech. -  mol. surf. projection");
      const bool skip = projectCenterOntoSurface(center, sCoords.col(i), R_i, r, grad);
      if (nThreads == 1)
        Timings::timeTaken("Tech. -  mol. surf. projection");
      if (skip)
        continue;
      const double radius = (r - center).norm();
      double weight = sWeights(i) * 4.0 * M_PI;
      if (weight < 0.0)
        continue;
      grad *= 1.0 / grad.norm();
      if (std::abs(radius - 1.0) > 1e-6) {
        const Eigen::Vector3d& centerOnUnit = sCoords.col(i);
        Eigen::Vector3d normalOnUnit = (centerOnUnit - center);
        normalOnUnit /= normalOnUnit.norm();
        std::shared_ptr<Ellipse> ellipse = projectCircleOntoSurface(sCoords.col(i), normalOnUnit, r, grad, weight);
        weight *= radius * radius;
        // Set the area to the actual area of the spherical part.
        // This will scale the radii accordingly.
        ellipse->scaleByArea(weight);

        // Check if the point is on the bond-part of the surface.
        // If so, weight reduction and shuffling by cutting with sphere--sphere planes is necessary.
        if (nThreads == 1)
          Timings::takeTime("Tech. -  mol. surf. scale area");
        for (unsigned int jSphere = 0; jSphere < _spheres.size(); ++jSphere) {
          if (iSphere == jSphere)
            continue; // No bond with itself.
          const Sphere& sphere_j = _spheres[jSphere];
          const Eigen::Vector3d r_j = sphere_j.getCenterCoords();
          Eigen::Vector3d rB = getBondProjection(r_i, r_j, r);
          const bool purelyOnReferenceSide = (rB - r_i).norm() < 1e-6;
          const bool purelyDistant = (rB - r_j).norm() < 1e-6;
          const double dij = R_i + sphere_j.getRadius() + 2.0 * _r_s;
          const bool bondIsFarAway = (rB - r).norm() > dij || (r - r_j).norm() > dij;
          if (purelyDistant) {
            // The bond is projected past the second sphere.
            // It can not contribute.
            weight = 0.0;
            break;
          }
          // The point is not on the bond or the bond is very far away, Thus, the j-sphere is not relevant.
          if (purelyOnReferenceSide || bondIsFarAway)
            continue;

          const Plane touchingPlane = getSphereSphereTouchingPlane(sphere, sphere_j);
          // Scale by skew projection on the touching plane.
          /*
           *        |
           *        A             Perforation point of skew projection with touching plane.
           *       /|
           *   a  / | b           scaling = b/a = sin[<(rA,rB)].
           *     /  |
           *    r---B             Orthogonal projection.
           *   /    |
           *  /     |
           * I -------------- J   Sphere centers I and J.
           *        |
           */
          const double angle = acos(normalOnUnit.dot(touchingPlane.getNormalVector()));
          const double scaling = sin(angle);
          ellipse->scaleByFactor(scaling);
          weight = ellipse->getArea();
          if (weight < 1e-9)
            break;

          // Check for cuts by the touching plane.
          // If cut, scale the area and move the center.
          CutEllipse farSideFragment;
          CutEllipse atomSideFragment;
          // Move the point in case the plane cuts the ellipse.
          bool movePoint = ellipse->cutWithPlane(touchingPlane, r_i, atomSideFragment, farSideFragment);
          if (atomSideFragment.area < 1e-9) {
            // The ellipse is completely hidden by the plane.
            // Stop here.
            weight = 0.0;
            break;
          }
          if (movePoint) {
            // If the center of gravity changed, we construct a new ellipse at the
            // shifted center with the given area/scaled radii.
            ellipse->scaleByArea(atomSideFragment.area);
            double guessDistance = (center - atomSideFragment.centerOfGravity).norm();
            // Shift the point back to the surface.
            projectCenterOntoSurface(center, atomSideFragment.centerOfGravity, guessDistance, r, grad);
            ellipse = std::make_shared<Ellipse>(r, ellipse->getR1(), ellipse->getR2());
          } // of movePoint
          weight = ellipse->getArea();
          if (weight < 1e-9)
            break;
        } // for jSphere
      }
      else {
        // Unlikely to happen, but not impossible. The surface is on the unit sphere.
        weight = sWeights(i);
      }
      if (nThreads == 1)
        Timings::timeTaken("Tech. -  mol. surf. scale area");
      // Omit points with small weights or that are very close to already existing points.
      // Save all significant points.
      if (weight > 1e-9) {
        finalCoordinates.push_back(r);
        finalNorms.push_back(grad);
        finalWeights.push_back(weight);
        finalSphereIndices.push_back(iSphere);
      } // if weight > 1e-9
    }   // for i
  }     // for sphere
  // Gather the points from each thread and check if they are to close to another.
  std::vector<Eigen::Vector3d> coords;
  std::vector<Eigen::Vector3d> norms;
  std::vector<double> weights;
  std::vector<unsigned int> sphereIndices;
  mergeGridPoints(allFinalCoordinates, allFinalNorms, allFinalWeights, allSphereIndices, coords, norms, weights,
                  sphereIndices, _spheres);

  if (_oneCavity)
    groupPoints(coords, norms, weights, sphereIndices);
  // Build grid controller and norm vectors.
  const unsigned int nCavityPoints = coords.size();
  std::unique_ptr<Eigen::Matrix3Xd> gridCoods = std::make_unique<Eigen::Matrix3Xd>(3, nCavityPoints);
  auto normVectors = std::make_unique<Eigen::Matrix3Xd>(3, nCavityPoints);
  std::unique_ptr<Eigen::VectorXd> gridWeights = std::make_unique<Eigen::VectorXd>(nCavityPoints);
  auto sphereIndexPairs = collectSphereIndices(sphereIndices, _spheres.size());
  for (unsigned int i = 0; i < nCavityPoints; ++i) {
    gridCoods->col(i) = coords[i];
    (*gridWeights)[i] = weights[i];
    normVectors->col(i) = norms[i];
  } // for i
  Eigen::setNbThreads(0);
  return std::make_unique<MolecularSurface>(std::move(gridCoods), std::move(gridWeights), std::move(normVectors),
                                            "DELLEY", sphereIndexPairs, sphereIndices, _r_s, _spheres);
}

void DelleySurfaceConstructor::mergeGridPoints(const std::vector<std::vector<Eigen::Vector3d>>& allCoordinates,
                                               const std::vector<std::vector<Eigen::Vector3d>>& allNorms,
                                               const std::vector<std::vector<double>>& allWeights,
                                               const std::vector<std::vector<unsigned int>>& allSphereIndices,
                                               std::vector<Eigen::Vector3d>& coords, std::vector<Eigen::Vector3d>& norms,
                                               std::vector<double>& weights, std::vector<unsigned int>& sphereIndices,
                                               const std::vector<Sphere>& spheres) {
  for (unsigned int threadId = 0; threadId < allCoordinates.size(); ++threadId) {
    const unsigned int nPointsThread = allCoordinates[threadId].size();
    for (unsigned int i = 0; i < nPointsThread; ++i) {
      bool skip = false;
      const Eigen::Vector3d& r_i = allCoordinates[threadId][i];
      if (!skip && _minimalDistance > 0.0) {
        for (unsigned int j = 0; j < coords.size(); ++j) {
          const double distance = (r_i - coords[j]).norm();
          if (distance < _minimalDistance) {
            skip = true;
            const Eigen::Vector3d newCoord = 0.5 * (coords[j] + r_i);
            const double newWeight = allWeights[threadId][i] + weights[j];
            Eigen::Vector3d newNorm = 0.5 * (norms[j] + allNorms[threadId][i]);
            newNorm /= newNorm.norm();
            unsigned int newIndex = sphereIndices[j];
            // Chose the index as the closest sphere center.
            // Note that this makes it necessary to shuffle the points a little bit
            // to insure that all points of a sphere are in consecutive order.
            const Eigen::Vector3d& center_i = spheres[allSphereIndices[threadId][i]].getCenterCoords();
            const Eigen::Vector3d& center_j = spheres[sphereIndices[j]].getCenterCoords();
            const double distToRi = (newCoord - center_i).norm();
            const double distToRj = (newCoord - center_j).norm();
            bool equalDistance = isEqual(distToRi, distToRj, NORMAL_D);
            if (equalDistance) {
              newIndex = std::min(newIndex, allSphereIndices[threadId][i]);
            }
            else if (distToRi < distToRj) {
              newIndex = allSphereIndices[threadId][i];
            }
            if (newIndex == sphereIndices[j]) {
              // Set the new values.
              coords[j] = newCoord;
              weights[j] = newWeight;
              norms[j] = newNorm;
            }
            else {
              // Append the new point and remove point i.
              coords.push_back(newCoord);
              norms.push_back(newNorm);
              weights.push_back(newWeight);
              sphereIndices.push_back(newIndex);
              // Remove the former point i.
              coords.erase(coords.begin() + j);
              norms.erase(norms.begin() + j);
              weights.erase(weights.begin() + j);
              sphereIndices.erase(sphereIndices.begin() + j);
            }
            break;
          } // if distance < _minimalDistance
        }   // for j
      }     // if !skip

      if (not skip) {
        coords.push_back(r_i);
        norms.push_back(allNorms[threadId][i]);
        weights.push_back(allWeights[threadId][i]);
        sphereIndices.push_back(allSphereIndices[threadId][i]);
      }
    } // for i
  }   // for threadId
}

std::shared_ptr<Ellipse> DelleySurfaceConstructor::projectCircleOntoSurface(const Eigen::Vector3d& centerOnUnit,
                                                                            const Eigen::Vector3d& normalOnUnit,
                                                                            const Eigen::Vector3d& centerOnSurface,
                                                                            const Eigen::Vector3d& normalOnSurface,
                                                                            double weightOnUnit) {
  Plane unitPlane = Plane(centerOnUnit, normalOnUnit);
  Plane surfacePlane = Plane(centerOnSurface, normalOnSurface);
  std::shared_ptr<Line> planeIntersection = unitPlane.calculateIntersection(surfacePlane);
  const double circleRadius = sqrt(weightOnUnit / M_PI);
  Eigen::Vector3d r1;
  Eigen::Vector3d r2;
  if (planeIntersection) {
    // Not parallel.
    // r1 will be parallel to the intersection direction
    // and r2 will be orthogonal to the intersection
    // and the normal vector.
    r1 = planeIntersection->getDirection(); // Direction is normalized by Line.
    r2 = r1.cross(normalOnSurface);
    // Length of r1 is given by the projection from the circle on
    // the unit sphere to the surface.
    Eigen::Vector3d r1ProjectUnit = centerOnUnit + circleRadius * r1;
    const Eigen::Vector3d projectionDirection = centerOnSurface - centerOnUnit;
    Line rayToSurface(r1ProjectUnit, projectionDirection);
    std::shared_ptr<Eigen::Vector3d> r1ProjectSurface = surfacePlane.calculateIntersection(rayToSurface);
    const double r1Length = (centerOnSurface - *r1ProjectSurface).norm();
    r1 *= r1Length;

    Eigen::Vector3d r2ProjectUnit = r1.cross(normalOnUnit);
    r2ProjectUnit = 1.0 / r2ProjectUnit.norm() * circleRadius * r2ProjectUnit + centerOnUnit;
    Line rayToSurface2(r2ProjectUnit, projectionDirection);
    std::shared_ptr<Eigen::Vector3d> r2ProjectSurface = surfacePlane.calculateIntersection(rayToSurface2);
    const double r2Length = (centerOnSurface - *r2ProjectSurface).norm();
    r2 = 1.0 / r2.norm() * r2Length * r2;
  }
  else {
    // It will be circle! Both radii have the same norm and are orthogonal.
    // Construct them by guessing. They only have to be orthogonal to the
    // normal vector.
    Eigen::Vector3d xAxis = {1, 0, 0};
    r1 = xAxis.cross(normalOnUnit);
    if (r1.norm() < 1e-6) {
      // xAxis and normalOnUnit are parallel.
      Eigen::Vector3d yAxis = {0, 1, 0};
      r1 = yAxis.cross(normalOnUnit);
    }
    r1 *= 1.0 / r1.norm() * circleRadius;
    r2 = r1.cross(normalOnUnit);
    r2 *= 1.0 / r2.norm() * circleRadius;
  }
  return std::make_shared<Ellipse>(centerOnSurface, r1, r2);
}

bool DelleySurfaceConstructor::projectCenterOntoSurface(const Eigen::Vector3d& center, const Eigen::Vector3d& unitPoint,
                                                        double initialGuess, Eigen::Vector3d& r, Eigen::Vector3d& grad) {
  Eigen::Vector3d direction = unitPoint - center;
  direction *= 1.0 / direction.norm();
  // Get starting point for the right and left limits.
  // The isosurface of a single sphere--> left limit.
  Eigen::Vector3d r_l = center + initialGuess * direction;
  Eigen::Vector3d r_r;
  Eigen::Vector3d grad_l;
  Eigen::Vector3d grad_r;
  double value_l = 1;
  double value_r = 1;
  bool found_l = false;
  bool found_r = false;
  // Check the most likely case: The isosurface is
  // dominated by the sphere associated to the point.
  grad_l = calculateGradientAndValue(r_l, value_l);
  double value = 1;
  if (std::abs(value_l) < 1e-6) {
    r = r_l;
    grad = grad_l;
    return false;
  }
  if (value_l <= 1e-7) {
    found_l = true;

    // Do a linear search along the direction until we hit
    // a non negative function value.
    // As a step length I assume an increase of the surface function
    // as exp(-alpha delta R/R_s). This is roughly the case for a single
    // atom with radius 1 in this model. A minimal step-length of 0.2 is ensured.
    double prefactor = 0;
    r_r = r_l;
    do {
      value_r = value_l;
      double deltaR = _r_s * log(-value_r) / (_alpha);
      deltaR = (deltaR < 0.2) ? 0.2 : deltaR;
      r_r += deltaR * direction;
      value_r = calculateFunctionValue(r_r);
      if (value_r >= 0.0 || std::abs(value_r) < 1e-7) {
        found_r = true;
        break;
      } // if
      prefactor += deltaR;
      r_l = r_r;
      value_l = value_r;
    } while (prefactor < _projectionCutOff);
  }
  else {
    found_r = true;
    r_r = r_l;
    value_r = value_l;
    r_l = r_r - 1.0 * direction;
    value_l = calculateFunctionValue(r_l);
    if (value_l <= 1e-7) {
      found_l = true;
    }
    else {
      r_l = r_r - 2.0 * direction;
      value_l = calculateFunctionValue(r_l);
      if (value_l <= 1e-7)
        found_l = true;
    }
  }
  // Stop here if the point is covered within the molecule.
  if (not found_l || not found_r)
    return true;
  // If not: Perform binary search for the isosurface.
  // The distance between upper and lower border is halfed in
  // each iteration. Thus, we know in advance how long it will
  // take to converge the coordinate the accuracy given below.
  // This converged the PCM energy to around 1e-6.
  const double coordinateAccuracy = 1e-4;
  const double initialDistance = (r_l - r_r).norm();
  const unsigned int nSteps = log2(initialDistance / coordinateAccuracy) + 1;
  unsigned int counter = 0;
  do {
    r = 0.5 * (r_l + r_r);
    value = calculateFunctionValue(r);
    if (std::abs(value) < 1e-6)
      break;

    if (value < 0) {
      r_l = r;
      value_l = value;
    }
    else if (value > 0) {
      r_r = r;
      value_r = value;
    }
    else {
      if (value != value)
        throw SerenityError("Nan detected!");
      // is zero.
      break;
    }
    ++counter;
  } while (counter < nSteps);
  grad = calculateGradientAndValue(r, value);
  return false;
}

double DelleySurfaceConstructor::calculateFunctionValue(const Eigen::Vector3d& r) {
  double value = 1.0;
  const unsigned int nSpheres = _spheres.size();
  for (unsigned int i = 0; i < nSpheres; ++i) {
    const Sphere& sphere_i = _spheres[i];
    const Eigen::Vector3d& r_i = sphere_i.getCenterCoords();
    const double R_i = sphere_i.getRadius();
    const double R_i2 = R_i * R_i;
    const double r_sRi2 = 2.0 * _r_s * R_i;
    const Eigen::Vector3d r_r_i = r - r_i;
    const double arg = (r_r_i.squaredNorm() - R_i2) / r_sRi2;
    if (arg < _cutOff)
      value += f(arg);
    for (Eigen::SparseMatrix<double>::InnerIterator itJ(_cylinderRadii, i); itJ; ++itJ) {
      const unsigned int j = itJ.row();
      const Sphere& sphere_j = _spheres[j];
      const Eigen::Vector3d& r_j = sphere_j.getCenterCoords();
      const Eigen::Vector3d r_z = getBondProjection(r_i, r_j, r);
      const double& R_z = itJ.value();
      const double& r_sRz2 = _cylinderParameters(i, j);
      const Eigen::Vector3d r_r_z = r - r_z;
      const double arg2 = (r_r_z.squaredNorm() - R_z * R_z) / r_sRz2;
      if (arg2 < _cutOff)
        value += f(arg2);
    } // for j
  }   // for i
  return value;
}

Eigen::Vector3d DelleySurfaceConstructor::calculateGradientAndValue(const Eigen::Vector3d& r, double& value) {
  Eigen::Vector3d grad = Eigen::Vector3d::Zero(3);
  value = 1.0;
  const unsigned int nSpheres = _spheres.size();
  for (unsigned int i = 0; i < nSpheres; ++i) {
    const Sphere& sphere_i = _spheres[i];
    const Eigen::Vector3d& r_i = sphere_i.getCenterCoords();
    const double R_i = sphere_i.getRadius();
    const double R_i2 = R_i * R_i;
    const double r_sRi2 = 2.0 * _r_s * R_i;
    const Eigen::Vector3d r_r_i = r - r_i;
    const double arg = (r_r_i.squaredNorm() - R_i2) / r_sRi2;
    if (arg < _cutOff) {
      const Eigen::Vector3d d_arg = 2.0 / r_sRi2 * r_r_i;
      grad += df(arg) * d_arg;
      value += f(arg);
    }
    for (Eigen::SparseMatrix<double>::InnerIterator itJ(_cylinderRadii, i); itJ; ++itJ) {
      const unsigned int j = itJ.row();
      const Sphere& sphere_j = _spheres[j];
      const Eigen::Vector3d& r_j = sphere_j.getCenterCoords();
      const Eigen::Vector3d r_z = getBondProjection(r_i, r_j, r);
      const double& R_z = itJ.value();
      const double& r_sRz2 = _cylinderParameters(i, j);
      const Eigen::Vector3d r_r_z = r - r_z;
      const double arg2 = (r_r_z.squaredNorm() - R_z * R_z) / r_sRz2;
      if (arg2 < _cutOff) {
        const Eigen::Vector3d d_arg2 = -2.0 / r_sRz2 * r_r_z;
        grad += df(arg2) * d_arg2;
        value += f(arg2);
      }
    } // for j
  }   // for i
  return grad;
}

std::shared_ptr<UnitSphere> DelleySurfaceConstructor::buildUnitSphere(unsigned int l) {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> w;
  unsigned int n = 0;
  auto unitSphere = std::make_shared<UnitSphere>();
  AtomGridFactory::_lebedevSphericalGrid(l, n, x, y, z, w);
  unitSphere->coordinates = Eigen::Matrix3Xd::Zero(3, n);
  unitSphere->weights = Eigen::VectorXd::Zero(n);
  for (unsigned int i = 0; i < n; ++i) {
    unitSphere->coordinates(0, i) = x[i];
    unitSphere->coordinates(1, i) = y[i];
    unitSphere->coordinates(2, i) = z[i];
    unitSphere->weights(i) = w[i];
  } // for i
  return unitSphere;
}

Eigen::Vector3d DelleySurfaceConstructor::getBondProjection(const Eigen::Vector3d& r_i, const Eigen::Vector3d& r_j,
                                                            const Eigen::Vector3d& r) {
  // Bond line between i and j starting.
  const Eigen::Vector3d r_ij = r_j - r_i;
  const double r_ijNorm = r_ij.norm();
  // Check if r_i and r_j are identical.
  if (r_ijNorm < 1e-6)
    return r_i;
  // Unit vector in r_ij direction.
  const Eigen::Vector3d e_rij = 1.0 / r_ijNorm * r_ij;
  // Relative coordinate of r with respect to i.
  const Eigen::Vector3d tmp = r - r_i;
  // Projection of r along the ij bond line.
  const double projection = tmp.dot(e_rij);
  // If r falls in between i and j, projection will be larger than 0 and smaller than |r_ij|.
  // Then: use the point of the projection on the bond line of i and j.
  // Else: use the bond end.
  if (projection > 0.0 && projection < r_ijNorm) {
    return r_i + projection * e_rij;
  }
  if (projection <= 0.0)
    return r_i;
  //(projection >= r_ijNorm) has to be true.
  return r_j;
}

const UnitSphere& DelleySurfaceConstructor::getUnitSphere(unsigned int l) {
  if (!_unitSpheres[l])
    _unitSpheres[l] = buildUnitSphere(l);
  return *_unitSpheres[l];
}

void DelleySurfaceConstructor::calculateCylinderRadiiAndParameters() {
  std::vector<Eigen::Triplet<double>> radii;
  const unsigned int nSpheres = _spheres.size();
  _cylinderParameters = Eigen::MatrixXd::Constant(nSpheres, nSpheres, -1);
  for (unsigned int iSphere = 0; iSphere < nSpheres; ++iSphere) {
    const Sphere& sphere_i = _spheres[iSphere];
    const Eigen::Vector3d& r_i = sphere_i.getCenterCoords();
    const double R_i = sphere_i.getRadius();
    const double S_i = R_i + _r_s;
    for (unsigned int jSphere = 0; jSphere < iSphere; ++jSphere) {
      const Sphere& sphere_j = _spheres[jSphere];
      const Eigen::Vector3d& r_j = sphere_j.getCenterCoords();
      const double R_j = sphere_j.getRadius();
      const double dij = (r_i - r_j).norm();
      const double S_j = R_j + _r_s;
      if (dij <= S_i + S_j) {
        const double R_z = R_z_ij(R_i, R_j, dij);
        radii.push_back(Eigen::Triplet<double>(jSphere, iSphere, R_z));
        radii.push_back(Eigen::Triplet<double>(iSphere, jSphere, R_z));
        const double r_sRz2 = 2.0 * _r_s * std::max(R_z, 1.0);
        _cylinderParameters(jSphere, iSphere) = r_sRz2;
        _cylinderParameters(iSphere, jSphere) = r_sRz2;
      } // if dij <= S_i + S_j
    }   // for jSphere
  }     // for iSphere
  _cylinderRadii.resize(0, 0);
  _cylinderRadii.resize(nSpheres, nSpheres);
  _cylinderRadii.setFromTriplets(radii.begin(), radii.end());
}

void DelleySurfaceConstructor::groupPoints(std::vector<Eigen::Vector3d>& coords, std::vector<Eigen::Vector3d>& norms,
                                           std::vector<double>& weights, std::vector<unsigned int>& sphereIndices) {
  unsigned int nGridPoints = coords.size();
  unsigned int extremeIndex = nGridPoints + 1;
  unsigned int extremeCoordinate = 0;
  for (unsigned int i = 0; i < nGridPoints; ++i) {
    if (std::abs(coords[i](0)) > extremeCoordinate) {
      extremeCoordinate = std::abs(coords[i](0));
      extremeIndex = i;
    }
  }
  // Build connection matrix. I will use the probe diameter as
  // a cut off. This eliminates points completely within the
  // cavity.
  std::vector<Eigen::Triplet<int>> triplets;
  const double cutOff = _connectivityFactor * _r_s;
  for (unsigned int i = 0; i < nGridPoints; ++i) {
    const Eigen::Vector3d& r_i = coords[i];
    triplets.push_back(Eigen::Triplet<int>(i, i, 1));
    for (unsigned int j = 0; j < i; ++j) {
      const Eigen::Vector3d& r_j = coords[j];
      const double distanceSquared = (r_i - r_j).squaredNorm();
      if (distanceSquared < cutOff) {
        triplets.push_back(Eigen::Triplet<int>(i, j, 1));
        triplets.push_back(Eigen::Triplet<int>(j, i, 1));
      }
    } // for j
  }   // for i

  Eigen::SparseMatrix<int> spareConnections(nGridPoints, nGridPoints);
  spareConnections.setFromTriplets(triplets.begin(), triplets.end());
  std::vector<unsigned int> pointCloud = buildPointCloud(spareConnections, extremeIndex);

  unsigned int nSpheres = _spheres.size();
  std::vector<std::vector<Eigen::Vector3d>> newCoords(nSpheres);
  std::vector<std::vector<Eigen::Vector3d>> newNorms(nSpheres);
  std::vector<std::vector<double>> newWeights(nSpheres);
  for (auto point : pointCloud) {
    const unsigned int index = sphereIndices[point];
    newCoords[index].push_back(coords[point]);
    newNorms[index].push_back(norms[point]);
    newWeights[index].push_back(weights[point]);
  }
  // Resort the points.
  coords.clear();
  norms.clear();
  weights.clear();
  sphereIndices.clear();
  for (unsigned int iSphere = 0; iSphere < nSpheres; ++iSphere) {
    coords.insert(coords.end(), newCoords[iSphere].begin(), newCoords[iSphere].end());
    norms.insert(norms.end(), newNorms[iSphere].begin(), newNorms[iSphere].end());
    weights.insert(weights.end(), newWeights[iSphere].begin(), newWeights[iSphere].end());
    std::vector<unsigned int> dummyIndices(newCoords[iSphere].size(), iSphere);
    sphereIndices.insert(sphereIndices.end(), dummyIndices.begin(), dummyIndices.end());
  }
}

std::vector<unsigned int> DelleySurfaceConstructor::buildPointCloud(Eigen::SparseMatrix<int>& connections, unsigned int seed) {
  std::vector<unsigned int> pointSet = {seed};
  unsigned int nPoints = 1;
  unsigned int i = 0;
  /*
   * The algorithms starts at the seed and selects all vertices connected to it
   * and adds them to the <pointSet> if not already contained in it. It then selects
   * the next entry of the <pointSet> (if any) and does the same. This process is
   * repeated until all vertices in the <pointSet> where checked for their connecting
   * vertices.
   */
  while (i < nPoints) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(connections, pointSet[i]); it; ++it) {
      unsigned int newIndex = it.row();
      if (std::find(pointSet.begin(), pointSet.end(), newIndex) == pointSet.end()) {
        pointSet.push_back(newIndex);
      }
    } // for it
    ++i;
    nPoints = pointSet.size();
  } // while i
  return pointSet;
}

std::vector<std::pair<unsigned int, unsigned int>>
DelleySurfaceConstructor::collectSphereIndices(std::vector<unsigned int> sphereIndices, unsigned int nSpheres) {
  unsigned int nPoints = sphereIndices.size();
  std::vector<std::pair<unsigned int, unsigned int>> indexPairs(nSpheres, {nPoints + 1, nPoints});
  for (unsigned int iPoint = 0; iPoint < nPoints;) {
    const unsigned int sphereIndex = sphereIndices[iPoint];
    if (indexPairs[sphereIndex].first <= indexPairs[sphereIndex].second) {
      throw SerenityError("ERROR: This sphere has already points assigned to it."
                          " The initial point list was not ordered as expected.");
    }
    const unsigned int firstPoint = iPoint;
    ++iPoint;
    // Loop until we get to the next sphere index.
    while (iPoint < nPoints && sphereIndex == sphereIndices[iPoint])
      ++iPoint;
    // The second index will always be larger than the first!
    indexPairs[sphereIndex].first = firstPoint;
    indexPairs[sphereIndex].second = iPoint;
  } // for iPoint
  return indexPairs;
}

std::map<unsigned int, std::shared_ptr<UnitSphere>> DelleySurfaceConstructor::_unitSpheres = {
    {0, nullptr},  {1, nullptr},  {2, nullptr},  {3, nullptr},  {4, nullptr},  {5, nullptr},  {6, nullptr},
    {7, nullptr},  {8, nullptr},  {9, nullptr},  {10, nullptr}, {11, nullptr}, {12, nullptr}, {13, nullptr},
    {14, nullptr}, {15, nullptr}, {16, nullptr}, {17, nullptr}, {18, nullptr}, {19, nullptr}, {20, nullptr}};

} /* namespace Serenity */
