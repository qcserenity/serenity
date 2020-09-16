/**
 * @file SparseMapsController.cpp
 *
 * @date May 9, 2019
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
#include "data/SparseMapsController.h"
/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/PAOSelecter.h"                        //PAO selection
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h" //Occ to Atom
#include "basis/AtomCenteredBasisController.h"                        //Mulliken populations.
#include "data/ElectronicStructure.h"                                 //Coefficients.
#include "data/OrbitalController.h"                                   //Coefficients.
#include "data/OrbitalPair.h"                                         //OrbitalPair definition
#include "data/PAOController.h"                                       //PAOs
#include "geometry/Geometry.h"                                        //Defintion of getNAtoms().
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"               //Definition of a k-set.

namespace Serenity {

SparseMapsController::SparseMapsController(std::shared_ptr<SystemController> system, std::shared_ptr<PAOController> paoController,
                                           std::shared_ptr<SparseMap> occupiedToPAOOrbitalMap,
                                           std::vector<std::shared_ptr<OrbitalPair>> closeOrbitalPairs,
                                           std::vector<std::shared_ptr<OrbitalPair>> distantOrbitalPairs,
                                           double mullikenThreshold, double orbitalToShellThreshold,
                                           double strongTripletMullikenThreshold, double weakTripletMullikenThreshold)
  : _system(system),
    _paoController(paoController),
    _occupiedToPAOOrbitalMap(occupiedToPAOOrbitalMap),
    _closeOrbitalPairs(closeOrbitalPairs),
    _distantOrbitalPairs(distantOrbitalPairs),
    _mullikenThreshold(mullikenThreshold),
    _orbitalToShellThreshold(orbitalToShellThreshold),
    _strongTripletMullikenThreshold(strongTripletMullikenThreshold),
    _weakTripletMullikenThreshold(weakTripletMullikenThreshold) {
  _orbitalPairs = _closeOrbitalPairs;
  _orbitalPairs.insert(_orbitalPairs.end(), _distantOrbitalPairs.begin(), _distantOrbitalPairs.end());
  assert(_system);
  assert(_system->getLastSCFMode() == Options::SCF_MODES::RESTRICTED &&
         "Sparse maps are only implemented for RESTRICTED.");
  assert(_occupiedToPAOOrbitalMap);
  assert(_paoController);
}

const Eigen::MatrixXd& SparseMapsController::getOrbitalWiseMullikenPopulations() {
  if (!_orbitalWiseMullikenPopulations) {
    auto coefficients = _system->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
    auto occupiedMullikenPops = MullikenPopulationCalculator<RESTRICTED>::calculateAtomwiseOrbitalPopulations(
        coefficients, _system->getOneElectronIntegralController()->getOverlapIntegrals(),
        _system->getAtomCenteredBasisController()->getBasisIndices());
    _orbitalWiseMullikenPopulations = std::make_shared<Eigen::MatrixXd>(occupiedMullikenPops);
  }
  return *_orbitalWiseMullikenPopulations;
}

SparseMap SparseMapsController::buildOccToAtomMap(double threshold) {
  unsigned int nAtoms = _system->getGeometry()->getNAtoms();
  auto nOcc = _system->getNOccupiedOrbitals<RESTRICTED>();
  auto coefficients = _system->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  const auto& occupiedMullikenPops = getOrbitalWiseMullikenPopulations();
  SparseMap toReturn(nAtoms, nOcc);
  std::vector<Eigen::Triplet<int>> tripletList;
  for (unsigned int iOcc = 0; iOcc < nOcc; ++iOcc) {
    unsigned int nAtomsAssigned = 0;
    // Use the Mulliken charges to assign atom indices to the orbitals for integral prescreening.
    for (unsigned int atomIndex = 0; atomIndex < occupiedMullikenPops.rows(); ++atomIndex) {
      if (std::fabs(occupiedMullikenPops(atomIndex, iOcc)) > threshold) {
        tripletList.push_back(Eigen::Triplet<int>(atomIndex, iOcc, 1));
        ++nAtomsAssigned;
      } // if mulliken
    }   // for atomIndex
    if (nAtomsAssigned == 0)
      throw SerenityError("Non existing fitting domain detected. Please adjust the fitting thresholds!");
  } // for iOcc
  toReturn.setFromTriplets(tripletList.begin(), tripletList.end());
  return toReturn;
}

const SparseMap& SparseMapsController::getOccToAtomMap() {
  if (!_occupiedOrbitalToAtomMap) {
    _occupiedOrbitalToAtomMap = std::make_shared<SparseMap>(buildOccToAtomMap(_mullikenThreshold));
  }
  return *_occupiedOrbitalToAtomMap;
}

const SparseMap& SparseMapsController::getAtomToPAOMap() {
  if (!_atomToPAOMap) {
    unsigned int nPAOs = _paoController->getNPAOs();
    unsigned int nAtoms = _system->getGeometry()->getNAtoms();
    auto restMullikenGrossPops = MullikenPopulationCalculator<RESTRICTED>::calculateAtomwiseOrbitalPopulations(
        _paoController->getAllPAOs(), _system->getOneElectronIntegralController()->getOverlapIntegrals(),
        _system->getAtomCenteredBasisController()->getBasisIndices());
    _atomToPAOMap = std::make_shared<SparseMap>(nPAOs, nAtoms);
    std::vector<Eigen::Triplet<int>> tripletList;
    for (unsigned int iPAO = 0; iPAO < nPAOs; ++iPAO) {
      // Use the mulliken charges to assign atom indices to the orbitals for integral prescreening.
      for (unsigned int atomIndex = 0; atomIndex < restMullikenGrossPops.rows(); ++atomIndex) {
        if (std::fabs(restMullikenGrossPops(atomIndex, iPAO)) > _mullikenThreshold) {
          tripletList.push_back(Eigen::Triplet<int>(iPAO, atomIndex, 1));
        }
      }
    }
    _atomToPAOMap->setFromTriplets(tripletList.begin(), tripletList.end());
  }
  return *_atomToPAOMap;
}

const SparseMap& SparseMapsController::getOccToPAOMap() {
  return *_occupiedToPAOOrbitalMap;
}

const SparseMap& SparseMapsController::getShellToOccMap() {
  if (!_shellToOccMap) {
    auto nOcc = _system->getNOccupiedOrbitals<RESTRICTED>();
    Eigen::MatrixXd coefficients = _system->getActiveOrbitalController<RESTRICTED>()->getCoefficients().leftCols(nOcc).eval();
    _shellToOccMap = constructShellToOrbitalMap(coefficients);
  }
  return *_shellToOccMap;
}

const SparseMap& SparseMapsController::getShellToPAOMap() {
  if (!_shellToPAOMap) {
    _shellToPAOMap = constructShellToOrbitalMap(_paoController->getAllPAOs());
  }
  return *_shellToPAOMap;
}

const SparseMap& SparseMapsController::getAtomToAuxShellMap() {
  if (!_atomToAuxShellMap) {
    auto auxBasis = _system->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
    unsigned int nAtoms = _system->getGeometry()->getNAtoms();
    _atomToAuxShellMap = std::make_shared<SparseMap>(auxBasis->getBasis().size(), nAtoms);
    std::vector<Eigen::Triplet<int>> tripletList;
    // Loop over atoms and their index pairs.
    auto idx = auxBasis->getBasisIndicesRed();
    unsigned int atomIndex = 0;
    for (const auto& index : idx) {
      // Loop over basis indices of this atom and add to the triplet list
      for (unsigned int iShell = index.first; iShell < index.second; ++iShell) {
        tripletList.push_back(Eigen::Triplet<int>(iShell, atomIndex, 1));
      }
      ++atomIndex;
    }
    // Build the map from triplets
    _atomToAuxShellMap->setFromTriplets(tripletList.begin(), tripletList.end());
  }
  return *_atomToAuxShellMap;
}

const SparseMap& SparseMapsController::getOccToAuxShellMap() {
  if (!_occToKMap) {
    _occToKMap = std::make_shared<SparseMap>(this->getAtomToAuxShellMap().rows(), this->getOccToAtomMap().cols());
    *_occToKMap = (this->getAtomToAuxShellMap() * this->getOccToAtomMap()).pruned().eval();
  }
  return *_occToKMap;
}

const SparseMap& SparseMapsController::getTripletOccToAuxShellMap(bool weak) {
  if (weak) {
    if (!_weakTripletOccToKMap) {
      const SparseMap occToAtom = buildOccToAtomMap(_weakTripletMullikenThreshold);
      _weakTripletOccToKMap = std::make_shared<SparseMap>(this->getAtomToAuxShellMap().rows(), occToAtom.cols());
      *_weakTripletOccToKMap = (this->getAtomToAuxShellMap() * occToAtom).pruned().eval();
    } // if !_weakTripletOccToKMap
    return *_weakTripletOccToKMap;
  }
  else {
    if (!_strongTripletOccToKMap) {
      const SparseMap occToAtom = buildOccToAtomMap(_strongTripletMullikenThreshold);
      _strongTripletOccToKMap = std::make_shared<SparseMap>(this->getAtomToAuxShellMap().rows(), occToAtom.cols());
      *_strongTripletOccToKMap = (this->getAtomToAuxShellMap() * occToAtom).pruned().eval();
    } // if !_strongTripletOccToKMap
    return *_strongTripletOccToKMap;
  } // else weak
}

const SparseMap& SparseMapsController::getExtendedOccToAuxShellMap() {
  if (!_extendedOccToK) {
    _extendedOccToK = buildExtendedMap(this->getOccToAuxShellMap());
  }
  return *_extendedOccToK;
}

const SparseMap& SparseMapsController::getCloseExtendedOccToAuxShellMap() {
  if (!_closeExtendedOccToK) {
    _closeExtendedOccToK = buildCloseExtendedMap(this->getOccToAuxShellMap());
  }
  return *_closeExtendedOccToK;
}

const SparseMap& SparseMapsController::getExtendedOccToPAOMap() {
  if (!_extendedOccToPAO) {
    _extendedOccToPAO = buildExtendedMap(this->getOccToPAOMap());
  }
  return *_extendedOccToPAO;
}

const SparseMap& SparseMapsController::getExtendedKtoRhoMap() {
  if (!_extendedAuxShellToRho) {
    const SparseMap& shellToOcc = this->getShellToOccMap();
    const SparseMap& extOccToAux = this->getExtendedOccToAuxShellMap();
    _extendedAuxShellToRho = std::make_shared<SparseMap>(shellToOcc.cols(), extOccToAux.rows());
    *_extendedAuxShellToRho = SparseMap((extOccToAux * shellToOcc).transpose()).pruned().eval();
  }
  return *_extendedAuxShellToRho;
}

const SparseMap& SparseMapsController::getExtendedKtoPAOMap() {
  if (!_extendedAuxShellToPAO) {
    const SparseMap& extOccToAux = this->getExtendedOccToAuxShellMap();
    const SparseMap& occToPAOMap = this->getExtendedOccToPAOMap();
    _extendedAuxShellToPAO = std::make_shared<SparseMap>(occToPAOMap.rows(), extOccToAux.rows());
    *_extendedAuxShellToPAO = SparseMap((extOccToAux * occToPAOMap.transpose()).transpose()).pruned().eval();
  }
  return *_extendedAuxShellToPAO;
}

const SparseMap& SparseMapsController::getExtendedKtoSigmaMap() {
  if (!_extendedAuxShellToSigma) {
    const SparseMap& extendedKtoPAOMap = getExtendedKtoPAOMap();
    const SparseMap& shellToPAO = this->getShellToPAOMap();
    _extendedAuxShellToSigma = std::make_shared<SparseMap>(shellToPAO.cols(), extendedKtoPAOMap.rows());
    *_extendedAuxShellToSigma = SparseMap((extendedKtoPAOMap.transpose() * shellToPAO).transpose()).pruned().eval();
  }
  return *_extendedAuxShellToSigma;
}

std::shared_ptr<SparseMap> SparseMapsController::constructShellToOrbitalMap(const Eigen::MatrixXd& coefficients) {
  auto basisController = _system->getAtomCenteredBasisController();
  Eigen::MatrixXd squaredC = coefficients.cwiseProduct(coefficients);
  squaredC.array() = squaredC.array().sqrt();
  Eigen::SparseMatrix<int> shellToOrbitalMap(coefficients.cols(), basisController->getReducedNBasisFunctions());
  std::vector<Eigen::Triplet<int>> tripletList;
  auto idx = basisController->getBasisIndices();
  auto idxRed = basisController->getBasisIndicesRed();
  for (unsigned int i = 0; i < squaredC.cols(); ++i) {
    for (unsigned int iAtom = 0; iAtom < idx.size(); ++iAtom) {
      unsigned int nFunc = idx[iAtom].second - idx[iAtom].first;
      int nImportant = (squaredC.block(idx[iAtom].first, i, nFunc, 1).array() > _orbitalToShellThreshold).count();
      if (nImportant > 0) {
        for (unsigned int iShell = idxRed[iAtom].first; iShell < idxRed[iAtom].second; ++iShell) {
          tripletList.push_back(Eigen::Triplet<int>(i, iShell, 1));
        } // for iShell
      }   // if nImportant > 0
    }     // for iAtom
  }       // for i
  // Build sparse map
  shellToOrbitalMap.setFromTriplets(tripletList.begin(), tripletList.end());
  return std::make_shared<SparseMap>(shellToOrbitalMap);
}

std::shared_ptr<SparseMap> SparseMapsController::buildCloseExtendedMap(const Eigen::SparseMatrix<int>& initialMap) {
  auto newMap = std::make_shared<SparseMap>(initialMap.rows(), initialMap.cols());
  assert(initialMap.cols() == _system->getNOccupiedOrbitals<RESTRICTED>());
  for (const auto& pair : _closeOrbitalPairs) {
    Eigen::SparseVector<int> sumVector = (initialMap.col(pair->i) + initialMap.col(pair->j)).eval();
    newMap->col(pair->i) += sumVector.eval();
    newMap->col(pair->j) += sumVector.eval();
  }
  newMap->pruned();
  return newMap;
}

std::shared_ptr<SparseMap> SparseMapsController::buildExtendedMap(const SparseMap& initialMap) {
  auto newMap = std::make_shared<SparseMap>(initialMap.rows(), initialMap.cols());
  assert(initialMap.cols() == _system->getNOccupiedOrbitals<RESTRICTED>());
  for (const auto& pair : _orbitalPairs) {
    Eigen::SparseVector<int> sumVector = (initialMap.col(pair->i) + initialMap.col(pair->j)).eval();
    newMap->col(pair->i) += sumVector.eval();
    newMap->col(pair->j) += sumVector.eval();
  }
  newMap->pruned();
  return newMap;
}

} /* namespace Serenity */
