/**
 * @file   EFieldPlates.cpp
 *
 * @date   Mai 25, 2020
 * @author Eric Niehoff
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
#include "geometry/efields/EFieldPlates.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "parameters/Constants.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace Serenity {

EFieldPlates::EFieldPlates(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double distance, unsigned nRings, double radius,
                           double fieldStrength, std::string nameOutput) {
  printSmallCaption("electric field plate generation");
  printf("  field strength (au): %4.2e\n\n", fieldStrength);
  printf("  reference points (Angstrom):\n");
  printf("  %7.3f %7.3f %7.3f\n", pos1(0) * BOHR_TO_ANGSTROM, pos1(1) * BOHR_TO_ANGSTROM, pos1(2) * BOHR_TO_ANGSTROM);
  printf("  %7.3f %7.3f %7.3f\n\n", pos2(0) * BOHR_TO_ANGSTROM, pos2(1) * BOHR_TO_ANGSTROM, pos2(2) * BOHR_TO_ANGSTROM);

  // convert input (Angstrom) to Bohr
  distance *= ANGSTROM_TO_BOHR;
  radius *= ANGSTROM_TO_BOHR;

  // vector connecting both reference points
  Eigen::Vector3d conVec = pos2 - pos1;

  // point in the middle of that vector
  Eigen::Vector3d midPoint = pos1 + 0.5 * conVec;

  conVec.normalize();

  // get vector in plane
  Eigen::Vector3d inPlane = this->getPerpendicularVector(conVec);

  // vector containing positive/negative point charges
  std::vector<Eigen::Vector3d> posPoints;
  std::vector<Eigen::Vector3d> negPoints;

  // center of plate 1
  Eigen::Vector3d centerPlate = midPoint - conVec * distance;
  posPoints.push_back(centerPlate);

  // create plate 1
  for (unsigned iRing = 1; iRing < nRings + 1; ++iRing) {
    // number of points on each ring
    unsigned nPoints = std::round(1.0 + PI / std::asin(1.0 / (2.0 * iRing)));
    // angle increment
    double angle = 2.0 * PI / nPoints;

    for (unsigned iPoint = 0; iPoint < nPoints; ++iPoint) {
      posPoints.push_back(centerPlate + Eigen::AngleAxisd(angle * iPoint, conVec) * inPlane * radius * iRing);
    }
  }
  // create plate2
  for (auto point : posPoints) {
    negPoints.push_back(point + 2.0 * conVec * distance);
  }

  // get initial start efield with charge = +- 1.0
  Eigen::Vector3d initEfield(0.0, 0.0, 0.0);
  for (auto point : posPoints) {
    Eigen::Vector3d vecPToP = midPoint - point;
    double dist = vecPToP.norm();
    initEfield += vecPToP.normalized() / (dist * dist);
  }
  for (auto point : negPoints) {
    Eigen::Vector3d vecPToP = midPoint - point;
    double dist = vecPToP.norm();
    initEfield -= vecPToP.normalized() / (dist * dist);
  }

  // scale point charges to satisfy efield strength
  _charge = fieldStrength / conVec.dot(initEfield);

  // write both plates to disk
  if (nameOutput != "") {
    std::ofstream file;
    file.open(nameOutput + ".xyz", std::ofstream::out | std::ofstream::trunc);
    file << posPoints.size() + negPoints.size() << std::endl;
    file << "ID: " << nameOutput << std::endl;
    file << std::fixed << std::setprecision(6);
    for (auto point : posPoints) {
      point *= BOHR_TO_ANGSTROM;
      file << "Xe"
           << "    " << point[0] << "    " << point[1] << "    " << point[2] << std::endl;
    }
    for (auto point : negPoints) {
      point *= BOHR_TO_ANGSTROM;
      file << "B"
           << "    " << point[0] << "    " << point[1] << "    " << point[2] << std::endl;
    }
    file.close();
  }

  // setup pair list
  for (auto point : posPoints) {
    _pairList.push_back(std::make_pair(+_charge, std::array<double, 3>{point[0], point[1], point[2]}));
  }
  for (auto point : negPoints) {
    _pairList.push_back(std::make_pair(-_charge, std::array<double, 3>{point[0], point[1], point[2]}));
  }
}

Eigen::Vector3d EFieldPlates::getPerpendicularVector(Eigen::Vector3d inVec) {
  Eigen::Vector3d x(1.0, 0.0, 0.0);
  Eigen::Vector3d y(0.0, 1.0, 0.0);

  if ((inVec - x).norm() > 0.05 && (inVec - x).norm() < 1.95) {
    Eigen::Vector3d outVector = inVec.cross(x);
    return outVector.normalized();
  }
  else {
    Eigen::Vector3d outVector = inVec.cross(y);
    return outVector.normalized();
  }
}

} /* namespace Serenity */
