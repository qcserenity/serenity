/**
 * @file GEPOLSurfaceConstructor.cpp
 *
 * @date Sep 27, 2016
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
#include "geometry/GEPOLSurfaceConstructor.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "geometry/MolecularSurface.h"
#include "geometry/Point.h"
#include "geometry/Sphere.h"
#include "geometry/Triangle.h"
#include "grid/GridController.h"
#include "io/FormattedOutputStream.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

GEPOLSurfaceConstructor::GEPOLSurfaceConstructor(std::vector<Sphere> spheres, bool isSAS, unsigned int patchLevel,
                                                 double minDistance, double minRadius, double probeRadius, double overlapFactor)
  : _initialized(false),
    _isSAS(isSAS),
    _calveLevel(patchLevel),
    _minDistance(minDistance),
    _minRad(minRadius),
    _solvRad(probeRadius),
    _overlapFactor(overlapFactor),
    _spheres(spheres) {
}

GEPOLSurfaceConstructor::~GEPOLSurfaceConstructor() = default;

std::unique_ptr<MolecularSurface> GEPOLSurfaceConstructor::getMolecularSurface() {
  if (!_initialized) {
    initializeSurface();
  }
  if (!_surface)
    throw SerenityError(
        "ERROR: Logic error in GEPOL Surface construction. The molecular surface is not available any more!");
  return std::move(_surface);
}

void GEPOLSurfaceConstructor::checkTriangleCenter(bool& withinSurface, bool& intersected, const Triangle& triangle,
                                                  const double& effSolvRad, const std::vector<Sphere>& spheres) {
  const Point& center = triangle.getCenter();
  withinSurface = false;
  intersected = false;
  for (auto sphere : spheres) {
    double dist = (center - sphere.getCenter()).distanceToOrigin();
    double cutOff = sphere.getRadius() + effSolvRad;
    if ((dist - cutOff) < -1e-9) {
      withinSurface = true;
      intersected = false;
      break;
    }
    auto trianglePoints = triangle.getPoints();
    for (auto point : trianglePoints) {
      double pointDist = (point - sphere.getCenter()).distanceToOrigin();
      if ((pointDist - cutOff) < -1e-9) {
        intersected = true;
        break;
      }
    }
  }
}

std::vector<Triangle> GEPOLSurfaceConstructor::patchIntersection(const Triangle& triangle, double effSolvRad,
                                                                 int patchLevel, const std::vector<Sphere>& spheres) {
  std::vector<Triangle> triangles = {};
  bool withinSurface, intersected;
  checkTriangleCenter(withinSurface, intersected, triangle, effSolvRad, spheres);
  if (withinSurface)
    return triangles;
  if (intersected and patchLevel > 0) {
    auto calvedTriangles = triangle.calve(true);
    for (auto tri : calvedTriangles) {
      auto patchedTriangles = patchIntersection(tri, effSolvRad, patchLevel - 1, spheres);
      for (auto patch : patchedTriangles)
        triangles.push_back(patch);
    }
  } // if intersected
  else {
    triangles = {triangle};
  }
  return triangles;
}

void GEPOLSurfaceConstructor::buildSurface(const std::vector<Sphere>& spheres) {
  std::vector<Triangle> allTriangles;
  std::vector<double> allRadii;
  std::vector<Point> allSphereCenters;
  std::vector<double> allAreas;
  double effSolvRad = (_isSAS) ? _solvRad : 0.0;

  for (auto sphere : spheres) {
    std::vector<Triangle> triangles = sphere.getSolventAccessibleSurface(effSolvRad);
    allTriangles.insert(allTriangles.end(), triangles.begin(), triangles.end());
  }
  std::vector<Triangle> screenedTriangles = {};
  for (auto triangle : allTriangles) {
    auto patchedTriangles = patchIntersection(triangle, effSolvRad, _calveLevel, spheres);
    for (auto patch : patchedTriangles)
      screenedTriangles.push_back(patch);
  }
  for (unsigned int i = 0; i < screenedTriangles.size(); ++i) {
    const Triangle& tri = screenedTriangles[i];
    const Point& i_center = tri.getCenter();
    bool nonClose = true;
    for (unsigned int j = 0; j < i; ++j) {
      const Point& j_center = screenedTriangles[j].getCenter();
      double dist = (i_center - j_center).distanceToOrigin();
      if (dist < _minDistance) {
        nonClose = false;
        break;
      }
    }
    if (nonClose)
      _molecularSurface.push_back(tri);
  }
  auto centerCoordinates = std::make_unique<Eigen::Matrix3Xd>(3, _molecularSurface.size());
  auto triangleAreas = std::make_unique<Eigen::VectorXd>(_molecularSurface.size());
  auto normalVectors = std::make_unique<Eigen::Matrix3Xd>(3, _molecularSurface.size());
  auto& coordinates = *centerCoordinates;
  auto& areas = *triangleAreas;
  auto& normVectors = *normalVectors;
  for (unsigned int i = 0; i < _molecularSurface.size(); ++i) {
    Triangle& triangle = _molecularSurface[i];
    const Point& center = triangle.getCenter();
    coordinates(0, i) = center.getX();
    coordinates(1, i) = center.getY();
    coordinates(2, i) = center.getZ();
    normVectors.col(i) = triangle.getNormVector();
    areas(i) = triangle.getArea(4);
  }
  std::string label = (_isSAS) ? "GEPOL-SAS" : "GEPOL-SES";
  // TODO  Keep track of the triangle centers ...
  std::vector<std::pair<unsigned int, unsigned int>> sphereIndices = {std::make_pair(0, centerCoordinates->cols())};
  std::vector<unsigned int> pointWiseIndicesOfSpheres = {0};
  _surface = std::make_unique<MolecularSurface>(std::move(centerCoordinates), std::move(triangleAreas),
                                                std::move(normalVectors), label, sphereIndices,
                                                pointWiseIndicesOfSpheres, _solvRad, spheres);
}

void GEPOLSurfaceConstructor::initializeSurface() {
  auto spheres = _spheres;
  addSpheres(true, spheres);
  addSpheres(false, spheres);
  buildSurface(spheres);
  _initialized = true;
}

void GEPOLSurfaceConstructor::addSpheres(bool coarse, std::vector<Sphere>& spheres) {
  while (true) {
    unsigned int nSpheresOld = 0;
    unsigned int nSpheres = spheres.size();
    for (unsigned int A = nSpheresOld; A < nSpheres; A++) {
      /*
       * Spheres with zero area do not contribute to
       * surface generation.
       */
      if (spheres[A].getSphereType() == SphereType::ZeroArea)
        continue;
      for (unsigned int B = 0; B < A; B++) {
        /*
         * Spheres with zero area do not contribute to
         * to surface generation.
         */
        if (spheres[B].getSphereType() == SphereType::ZeroArea)
          continue;
        /*
         * Get the data of the two spheres. Distinguish
         * between smaller and bigger sphere.
         */
        double dist = spheres[A].distanceTo(spheres[B]);
        double smallRad = 0;
        double bigRad = 0;
        Point smallCenter(0, 0, 0);
        Point bigCenter(0, 0, 0);
        if (spheres[A].getRadius() <= spheres[B].getRadius()) {
          smallRad = spheres[A].getRadius();
          bigRad = spheres[B].getRadius();
          smallCenter = spheres[A].getCenter();
          bigCenter = spheres[B].getCenter();
        }
        else {
          smallRad = spheres[B].getRadius();
          bigRad = spheres[A].getRadius();
          smallCenter = spheres[B].getCenter();
          bigCenter = spheres[A].getCenter();
        }
        /*
         * Check if solvent can pass through pair of spheres.
         */
        if (dist >= smallRad + bigRad + 2 * _solvRad)
          continue;
        /*
         * If not: check if spheres are overlapped enough
         */
        if (dist < bigRad + smallRad * (1 - 2 * _overlapFactor))
          continue;
        /*
         * If both tests were passed: Create new sphere with
         * radius bigRad centered right between the surfaces of
         * the two spheres. This step is the main difference
         * between the Coarse and the Fine generation of spheres
         * which is why a switch is used at this step.
         */
        double lambda = (dist + bigRad - smallRad) / (dist + smallRad - bigRad);
        Point newCenter = (smallCenter * lambda + bigCenter) / (1 + lambda);
        double newRad = 0.0;
        Sphere newSphere(newCenter, newRad, SphereType::Ghost, 0);
        /*
         * Boolean to check whether the newSphere should be kept
         * or not
         */
        bool save = true;
        if (coarse) {
          /*
           * No different spheretypes needed here
           */
          newSphere = Sphere(newCenter, bigRad, SphereType::Ghost, 0);
        }
        else {
          newRad = bigRad + _solvRad;
          newRad *= newRad;
          double omega = (dist + bigRad - smallRad) / 2;
          newRad += omega * (omega - (newRad + dist * dist - (smallRad + _solvRad) * (smallRad + _solvRad)) / dist);
          if (dist < bigRad + smallRad) {
            /*
             * Create sphere of type A
             */
            newRad = sqrt(newRad);
            newRad -= _solvRad;
            newSphere = Sphere(newCenter, newRad, SphereType::Ghost, 0);
          }
          else {
            if (newRad > (_solvRad + _minRad) * (_solvRad + _minRad)) {
              /*
               * Create sphere of type B
               */
              newRad = sqrt(newRad);
              newRad -= _solvRad;
              newSphere = Sphere(newCenter, newRad, SphereType::Ghost, 0);
            }
            else {
              /*
               * Create sphere of type C
               */
              newRad = bigRad + _solvRad;
              newRad *= newRad;
              newRad += bigRad * (bigRad - (newRad + dist * dist - (smallRad + _solvRad) * (smallRad + _solvRad)) / dist);
              newRad = sqrt(newRad);
              newRad -= _solvRad;
              lambda = bigRad / (dist - bigRad);
              newCenter = (bigCenter + smallCenter * lambda) / (1 + lambda);
              newSphere = Sphere(newCenter, newRad, SphereType::Ghost, 0);
            }
          }
        }
        if (newSphere.getRadius() < _minRad)
          continue;
        /*
         * Check if the newSphere is valid by several criteria
         */
        for (unsigned C = 0; C < nSpheres; C++) {
          Sphere otherSphere = spheres[C];
          dist = newSphere.distanceTo(otherSphere);
          /*
           * Boolean to check whether a sphere from the original
           * set needs to be deleted (which is the case if it
           * is engulfed by the newSphere).
           */
          bool markedToDeath = false;
          if (newSphere.getRadius() <= otherSphere.getRadius()) {
            smallRad = newSphere.getRadius();
            bigRad = otherSphere.getRadius();
          }
          else {
            markedToDeath = true;
            smallRad = otherSphere.getRadius();
            bigRad = newSphere.getRadius();
          }
          /*
           * Check if one sphere engulfs the other. If a sphere
           * of the original set is engulfed (i.e. if it is smaller
           * than newSphere) it will be deleted later. If newSphere
           * is engulfed, it will be discarded directly.
           */
          if ((dist * dist) < ((smallRad - bigRad) * (smallRad - bigRad))) {
            if (markedToDeath)
              spheres[C].setSphereType(SphereType::Engulfed);
            else {
              save = false;
              break;
            }
          }
          /*
           * Check if the new sphere overlaps strongly with any
           * of the existing spheres. If so: discard newSphere.
           */
          if (dist <= bigRad + smallRad * (1 - 2 * _overlapFactor)) {
            save = false;
            break;
          }
        }
        if (save) {
          /*
           * Discretize surface and check whether points lie
           * outside of the accessible surface. If so: discard
           * newSphere
           */
          auto surfaceTriangles = newSphere.getSolventAccessibleSurface(_solvRad);
          std::vector<Point> surfacePoints;
          for (auto triangle : surfaceTriangles) {
            surfacePoints.push_back(triangle.getCenter());
          }
          bool hasZeroArea = true;
          bool accessibleSurfaceContr = false;
          for (auto point : surfacePoints) {
            bool isOutsideSurface = true;
            for (unsigned C = 0; C < nSpheres; C++) {
              Sphere otherSphere = spheres[C];
              /*
               * Check if point is within allowed distance (=sphereRadius+solvRadius).
               * If true, this point is within the surface.
               */
              if ((point - otherSphere.getCenter()).distanceToOrigin() < otherSphere.getRadius() + _solvRad) {
                isOutsideSurface = false;
              }
              /*
               * Check if newSphere generates any new volume. If not, the newSphere
               * will participate in the generation of new spheres (but still it
               * will be kept).
               */
              if ((point - otherSphere.getCenter()).distanceToOrigin() > otherSphere.getRadius()) {
                hasZeroArea = false;
              }
            }
            if (isOutsideSurface) {
              accessibleSurfaceContr = true;
              break;
            }
          }
          if (hasZeroArea) {
            newSphere.setSphereType(SphereType::ZeroArea);
          }
          /*
           * If newSphere survived until this point it
           * deserves to live...
           *
           * I am a merciful god.
           */
          if (!accessibleSurfaceContr and save)
            spheres.push_back(newSphere);
        }
      }
    }
    /*
     * Set nSpheresOld and nSpheres so that the outermost for-loop
     * only covers the new spheres.
     */

    /*
     * Delete engulfed spheres.
     */
    unsigned int i = 0;
    for (std::vector<Sphere>::iterator it = spheres.begin(); it < spheres.end(); ++it) {
      if (spheres[i].getSphereType() == SphereType::Engulfed) {
        spheres.erase(it);
      }
      i++;
    }
    nSpheresOld = nSpheres;
    nSpheres = spheres.size();

    /*
     * Check if new spheres were added. If not: Stop
     */
    if (nSpheresOld == nSpheres) {
      break;
    }
  }
}

} /* namespace Serenity */
