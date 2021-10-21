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

Matrix<double> CoreCoreRepulsionDerivative::calculateDerivative(const std::vector<std::shared_ptr<Atom>>& atoms) {
  const unsigned int nAtoms = atoms.size();
  Matrix<double> gradient = Eigen::MatrixXd::Zero(nAtoms, 3);

  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = 0; j < nAtoms; ++j) {
      if (i == j) {
        // Do nothing for repulsion with self
        continue;
      }
      else {
        const double dist = distance(*atoms[i], *atoms[j]);
        // Multiplication of Charge(i) at the end!
        const double jChargeOverDistCube = atoms[j]->getEffectiveCharge() / (dist * dist * dist);
        // X
        gradient(i, 0) += jChargeOverDistCube * (atoms[j]->getX() - atoms[i]->getX());
        // Y
        gradient(i, 1) += jChargeOverDistCube * (atoms[j]->getY() - atoms[i]->getY());
        // Z
        gradient(i, 2) += jChargeOverDistCube * (atoms[j]->getZ() - atoms[i]->getZ());
      }
    }
    gradient(i, 0) *= atoms[i]->getEffectiveCharge();
    gradient(i, 1) *= atoms[i]->getEffectiveCharge();
    gradient(i, 2) *= atoms[i]->getEffectiveCharge();
  }
  return gradient;
}

Matrix<double> CoreCoreRepulsionDerivative::calculateDerivative(const std::vector<std::shared_ptr<Atom>>& atomsAct,
                                                                const std::vector<std::shared_ptr<Atom>>& atomsEnv) {
  const unsigned int nAtomsAct = atomsAct.size();
  const unsigned int nAtomsEnv = atomsEnv.size();
  Matrix<double> gradient(nAtomsAct, 3);
  gradient.setZero();

  for (unsigned int i = 0; i < nAtomsAct; ++i) {
    for (unsigned int j = 0; j < nAtomsEnv; ++j) {
      const double dist = distance(*atomsAct[i], *atomsEnv[j]);
      const double jChargeOverDistCube = atomsEnv[j]->getEffectiveCharge() / (dist * dist * dist);
      // X
      gradient(i, 0) += jChargeOverDistCube * (atomsEnv[j]->getX() - atomsAct[i]->getX());
      // Y
      gradient(i, 1) += jChargeOverDistCube * (atomsEnv[j]->getY() - atomsAct[i]->getY());
      // Z
      gradient(i, 2) += jChargeOverDistCube * (atomsEnv[j]->getZ() - atomsAct[i]->getZ());
    }
    gradient(i, 0) *= atomsAct[i]->getEffectiveCharge();
    gradient(i, 1) *= atomsAct[i]->getEffectiveCharge();
    gradient(i, 2) *= atomsAct[i]->getEffectiveCharge();
  }
  return gradient;
}

} /* namespace Serenity */
