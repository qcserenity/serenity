/**
 * @file CM5PopulationCalculator.cpp
 *
 * @date October 1, 2024
 * @author: Thorben Wiegmann
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
#include "analysis/populationAnalysis/CM5PopulationCalculator.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "parameters/AtomicParameters.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <cmath>
#include <iostream>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CM5PopulationCalculator<SCFMode>::CM5PopulationCalculator(std::shared_ptr<SystemController> system,
                                                          std::shared_ptr<HirshfeldPopulationCalculator<SCFMode>> hirshFeld)
  : _system(system), _hirshFeld(hirshFeld), _atomPopulations(nullptr) {
  _atoms = _system->getGeometry()->getAtoms();
  _nAtoms = _atoms.size();
  _atSymbols = _system->getGeometry()->getAtomSymbols();
}

template<Options::SCF_MODES SCFMode>
void CM5PopulationCalculator<SCFMode>::calculateCM5Populations() {
  Eigen::MatrixXd paulingBondOrder = this->getPaulingBondOrderMatrix();
  Eigen::MatrixXd parameterMatrix = this->getParameterMatrix();
  // calculate correction to Hirshfeld charges (eq. 1 in [1])
  Eigen::MatrixXd correction = paulingBondOrder.array() * parameterMatrix.array();
  auto q = _hirshFeld->getAtomPopulations();
  _atomPopulations = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(Eigen::VectorXd::Zero(_nAtoms));
  SpinPolarizedData<SCFMode, Eigen::VectorXd>& atomPopulations = *_atomPopulations;
  double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 1.0 : 0.5;
  for_spin(atomPopulations, q) {
    atomPopulations_spin = q_spin - scfFactor * correction.rowwise().sum();
  };
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd CM5PopulationCalculator<SCFMode>::getPaulingBondOrderMatrix() {
  Eigen::MatrixXd paulingBondOrder = Eigen::MatrixXd(_nAtoms, _nAtoms).setZero();
  for (unsigned int i = 0; i < _nAtoms; i++) {
    for (unsigned int j = 0; j < _nAtoms; j++) {
      if (i == j) {
        continue;
      }
      paulingBondOrder(i, j) = std::exp(-_alpha * (distance(*_atoms[i], *_atoms[j]) -
                                                   COVALENT_RADII[_atoms[i]->getNuclearCharge()] * ANGSTROM_TO_BOHR -
                                                   COVALENT_RADII[_atoms[j]->getNuclearCharge()] * ANGSTROM_TO_BOHR));
    }
  }
  return paulingBondOrder;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd CM5PopulationCalculator<SCFMode>::getParameterMatrix() {
  Eigen::MatrixXd parameterMatrix = Eigen::MatrixXd(_nAtoms, _nAtoms).setZero();
  for (unsigned int i = 0; i < _nAtoms; i++) {
    for (unsigned int j = 0; j < _nAtoms; j++) {
      std::map<std::string, double>::iterator iter = _bondParameters.find(_atSymbols[i] + _atSymbols[j]);
      if (iter != _bondParameters.end()) {
        parameterMatrix(i, j) = iter->second;
      }
      else {
        parameterMatrix(i, j) = _atomicParameters[_atSymbols[i]] - _atomicParameters[_atSymbols[j]];
      }
    }
  }
  return parameterMatrix;
}

template class CM5PopulationCalculator<Options::SCF_MODES::RESTRICTED>;
template class CM5PopulationCalculator<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity