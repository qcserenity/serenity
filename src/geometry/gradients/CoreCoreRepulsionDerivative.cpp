/**
 * @file   CoreCoreRepulsionDerivative.cpp
 *
 * @date   Nov 21, 2014
 * @author k_klah01
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
#include "geometry/gradients/CoreCoreRepulsionDerivative.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"

namespace Serenity {

std::vector<std::pair<double, Point>>
CoreCoreRepulsionDerivative::convertAtomsToCharges(const std::vector<std::shared_ptr<Atom>>& atoms) {
  std::vector<std::pair<double, Point>> convertedPointCharges;
  for (const auto& atom : atoms) {
    convertedPointCharges.emplace_back(atom->getEffectiveCharge(), *atom);
  }
  return convertedPointCharges;
}

Matrix<double> CoreCoreRepulsionDerivative::calculateDerivative(const std::vector<std::shared_ptr<Atom>>& atoms) {
  return calculateDerivative(convertAtomsToCharges(atoms));
}

Matrix<double> CoreCoreRepulsionDerivative::calculateDerivative(const std::vector<std::shared_ptr<Atom>>& atomsAct,
                                                                const std::vector<std::shared_ptr<Atom>>& atomsEnv) {
  return calculateDerivative(convertAtomsToCharges(atomsAct), convertAtomsToCharges(atomsEnv));
}

Matrix<double> CoreCoreRepulsionDerivative::calculateDerivative(const std::vector<std::pair<double, Point>>& charges) {
  const unsigned int nCharges = charges.size();
  Matrix<double> gradient = Eigen::MatrixXd::Zero(nCharges, 3);
  for (unsigned int i = 0; i < nCharges; ++i) {
    for (unsigned int j = 0; j < nCharges; ++j) {
      if (i == j) {
        // Do nothing for repulsion with self
        continue;
      }
      const double dist = distance(charges[i].second, charges[j].second);
      // Multiplication of Charge(i) at the end!
      const double jChargeOverDistCube = charges[j].first / (dist * dist * dist);
      // X
      gradient(i, 0) += jChargeOverDistCube * (charges[j].second.getX() - charges[i].second.getX());
      // Y
      gradient(i, 1) += jChargeOverDistCube * (charges[j].second.getY() - charges[i].second.getY());
      // Z
      gradient(i, 2) += jChargeOverDistCube * (charges[j].second.getZ() - charges[i].second.getZ());
    }
    gradient(i, 0) *= charges[i].first;
    gradient(i, 1) *= charges[i].first;
    gradient(i, 2) *= charges[i].first;
  }
  return gradient;
}
Matrix<double> CoreCoreRepulsionDerivative::calculateDerivative(const std::vector<std::pair<double, Point>>& chargesAct,
                                                                const std::vector<std::pair<double, Point>>& chargesEnv) {
  const unsigned int nChargesAct = chargesAct.size();
  const unsigned int nChargesEnv = chargesEnv.size();
  Matrix<double> gradient(nChargesAct, 3);
  gradient.setZero();

  for (unsigned int i = 0; i < nChargesAct; ++i) {
    for (unsigned int j = 0; j < nChargesEnv; ++j) {
      const double dist = distance(chargesAct[i].second, chargesEnv[j].second);
      const double jChargeOverDistCube = chargesEnv[j].first / (dist * dist * dist);
      // X
      gradient(i, 0) += jChargeOverDistCube * (chargesEnv[j].second.getX() - chargesAct[i].second.getX());
      // Y
      gradient(i, 1) += jChargeOverDistCube * (chargesEnv[j].second.getY() - chargesAct[i].second.getY());
      // Z
      gradient(i, 2) += jChargeOverDistCube * (chargesEnv[j].second.getZ() - chargesAct[i].second.getZ());
    }
    gradient(i, 0) *= chargesAct[i].first;
    gradient(i, 1) *= chargesAct[i].first;
    gradient(i, 2) *= chargesAct[i].first;
  }
  return gradient;
}
Matrix<double> CoreCoreRepulsionDerivative::calculateDerivative(const std::vector<std::pair<double, Point>>& chargesAct,
                                                                const std::vector<std::shared_ptr<Atom>>& atomsEnv) {
  return calculateDerivative(chargesAct, convertAtomsToCharges(atomsEnv));
}
Matrix<double> CoreCoreRepulsionDerivative::calculateDerivative(const std::vector<std::shared_ptr<Atom>>& atomsAct,
                                                                const std::vector<std::pair<double, Point>>& chargesEnv) {
  return calculateDerivative(convertAtomsToCharges(atomsAct), chargesEnv);
}

} /* namespace Serenity */
