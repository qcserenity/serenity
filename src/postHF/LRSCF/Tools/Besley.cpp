/**
 * @file Besley.cpp
 *
 * @date Dec 06, 2018
 * @author Johannes Toelle, Niklas Niemeyer
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
#include "postHF/LRSCF/Tools/Besley.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "integrals/OneElectronIntegralController.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
Besley<SCFMode>::Besley(std::shared_ptr<SystemController> system, unsigned int nAtoms, std::vector<double> besleyCutoff)
  : _system(system), _nAtoms(nAtoms), _besleyCutoff(besleyCutoff) {
  _nAtoms = _system->getNAtoms() < _nAtoms ? _system->getNAtoms() : _nAtoms;
  if (_nAtoms == _system->getNAtoms())
    printf("ATTENTION: Using all atoms for Besley restriction!");
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, std::vector<unsigned int>> Besley<SCFMode>::getWhiteList() {
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> whiteList(0);

  auto mulliken = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
      _system->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->getCoefficients(),
      _system->getOneElectronIntegralController()->getOverlapIntegrals(),
      _system->getAtomCenteredBasisController()->getBasisIndices());
  auto nOccupied = _system->getNOccupiedOrbitals<SCFMode>();
  auto coefficients = _system->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->getCoefficients();
  auto basisIndices = _system->getAtomCenteredBasisController()->getBasisIndices();
  unsigned int nActAOs = basisIndices[_nAtoms - 1].second;

  // Determine which occupied orbitals will be further used (eq. 9 of reference given in header file)
  for_spin(nOccupied, mulliken, whiteList) {
    for (unsigned int iMo = 0; iMo < nOccupied_spin; ++iMo) {
      double value = 0.0;
      for (unsigned int iAtom = 0; iAtom < _nAtoms; ++iAtom) {
        value += mulliken_spin(iAtom, iMo);
      }
      if (value > _besleyCutoff[0])
        whiteList_spin.push_back(iMo);
    }
  }; /* Occupied */

  // Determine which virtual orbitals will be further used (eq. 10 of reference given in header file)
  for_spin(nOccupied, coefficients, whiteList) {
    for (unsigned int iMo = nOccupied_spin; iMo < coefficients_spin.cols(); ++iMo) {
      double value = coefficients_spin.col(iMo).segment(0, nActAOs).array().square().sum();
      double norm = coefficients_spin.col(iMo).array().square().sum();
      if (value / norm > _besleyCutoff[1])
        whiteList_spin.push_back(iMo);
    }
  }; /* Virtual */

  return whiteList;
}

template class Besley<Options::SCF_MODES::RESTRICTED>;
template class Besley<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
