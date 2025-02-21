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
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h" //Occ to Atom
#include "basis/AtomCenteredBasisController.h"                        //Mulliken populations.
#include "data/OrbitalController.h"                                   //Coefficients.
#include "data/OrbitalPair.h"                                         //OrbitalPair definition
#include "data/OrbitalTriple.h"                                       //OrbitalTriple definition.
#include "data/PAOController.h"                                       //PAOs
#include "geometry/Geometry.h"                                        //Defintion of getNAtoms().
#include "integrals/OneElectronIntegralController.h"                  //Overlap interals.
#include "io/FormattedOutputStream.h"                                 //Filtered output streams.
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"               //Definition of a k-set.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"                     //Map extension by kl sets.
#include "system/SystemController.h"                                  //System controller definition.

namespace Serenity {

SparseMapsController::SparseMapsController(
    std::shared_ptr<SystemController> system, std::shared_ptr<PAOController> paoController,
    std::shared_ptr<SparseMap> occupiedToPAOOrbitalMap, std::vector<std::shared_ptr<OrbitalPair>> closeOrbitalPairs,
    std::vector<std::shared_ptr<OrbitalPair>> distantOrbitalPairs, Eigen::VectorXd orbitalWiseMullikenThresholds,
    Eigen::VectorXd orbitalToShellThresholds, double strongTripletMullikenScaling, double weakTripletMullikenScaling,
    bool klListExtension, std::vector<std::shared_ptr<OrbitalTriple>> triples)
  : _system(system),
    _virtualCoefficients(std::make_shared<Eigen::MatrixXd>(paoController->getAllPAOs())),
    _occupiedToPAOOrbitalMap(occupiedToPAOOrbitalMap),
    _closeOrbitalPairs(closeOrbitalPairs),
    _distantOrbitalPairs(distantOrbitalPairs),
    _orbitalWiseMullikenThresholds(orbitalWiseMullikenThresholds),
    _orbitalToShellThresholds(orbitalToShellThresholds),
    _strongTripletMullikenScaling(strongTripletMullikenScaling),
    _weakTripletMullikenScaling(weakTripletMullikenScaling),
    _klListExtension(klListExtension),
    _triples(triples) {
  _orbitalPairs = _closeOrbitalPairs;
  _orbitalPairs.insert(_orbitalPairs.end(), _distantOrbitalPairs.begin(), _distantOrbitalPairs.end());
  if (system->getLastSCFMode() != RESTRICTED) {
    throw SerenityError("This sparse map controller constructor is implemented for RESTRICTED orbitals."
                        " If you are using UNRESTRICTED orbitals, please use a different constructor");
  }
  unsigned int nOcc = _system->getNOccupiedOrbitals<RESTRICTED>();
  const Eigen::MatrixXd coefficients = _system->getActiveOrbitalController<RESTRICTED>()->getCoefficients().leftCols(nOcc);
  _occupiedCoefficients = std::make_shared<Eigen::MatrixXd>(coefficients);
}

SparseMapsController::SparseMapsController(
    std::shared_ptr<SystemController> system, std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients,
    std::shared_ptr<Eigen::MatrixXd> virtualCoefficients, std::shared_ptr<SparseMap> occupiedToPAOOrbitalMap,
    std::vector<std::shared_ptr<OrbitalPair>> closeOrbitalPairs,
    std::vector<std::shared_ptr<OrbitalPair>> distantOrbitalPairs, Eigen::VectorXd orbitalWiseMullikenThresholds,
    Eigen::VectorXd orbitalToShellThresholds, double strongTripletMullikenScaling, double weakTripletMullikenScaling,
    bool klListExtension, std::vector<std::shared_ptr<OrbitalTriple>> triples)
  : _system(system),
    _occupiedCoefficients(occupiedCoefficients),
    _virtualCoefficients(virtualCoefficients),
    _occupiedToPAOOrbitalMap(occupiedToPAOOrbitalMap),
    _closeOrbitalPairs(closeOrbitalPairs),
    _distantOrbitalPairs(distantOrbitalPairs),
    _orbitalWiseMullikenThresholds(orbitalWiseMullikenThresholds),
    _orbitalToShellThresholds(orbitalToShellThresholds),
    _strongTripletMullikenScaling(strongTripletMullikenScaling),
    _weakTripletMullikenScaling(weakTripletMullikenScaling),
    _klListExtension(klListExtension),
    _triples(triples) {
  _orbitalPairs = _closeOrbitalPairs;
  _orbitalPairs.insert(_orbitalPairs.end(), _distantOrbitalPairs.begin(), _distantOrbitalPairs.end());
}

const Eigen::MatrixXd& SparseMapsController::getOrbitalWiseMullikenPopulations() {
  if (!_orbitalWiseMullikenPopulations) {
    const Eigen::MatrixXd& coefficients = *_occupiedCoefficients;
    auto occupiedMullikenPops = MullikenPopulationCalculator<RESTRICTED>::calculateAtomwiseOrbitalPopulations(
        coefficients, _system->getOneElectronIntegralController()->getOverlapIntegrals(),
        _system->getAtomCenteredBasisController()->getBasisIndices());
    _orbitalWiseMullikenPopulations = std::make_shared<Eigen::MatrixXd>(occupiedMullikenPops);
  }
  return *_orbitalWiseMullikenPopulations;
}

SparseMap SparseMapsController::buildOccToAtomMap(Eigen::VectorXd thresholds) {
  unsigned int nAtoms = _system->getGeometry()->getNAtoms();
  auto nOcc = _occupiedCoefficients->cols();
  const auto& occupiedMullikenPops = getOrbitalWiseMullikenPopulations();
  SparseMap toReturn(nAtoms, nOcc);
  std::vector<Eigen::Triplet<int>> tripletList;
  for (unsigned int iOcc = 0; iOcc < nOcc; ++iOcc) {
    unsigned int nAtomsAssigned = 0;
    // Use the Mulliken charges to assign atom indices to the orbitals for integral prescreening.
    for (unsigned int atomIndex = 0; atomIndex < occupiedMullikenPops.rows(); ++atomIndex) {
      if (std::fabs(occupiedMullikenPops(atomIndex, iOcc)) > thresholds(iOcc)) {
        tripletList.push_back(Eigen::Triplet<int>(atomIndex, iOcc, 1));
        ++nAtomsAssigned;
      } // if mulliken
    }   // for atomIndex
    if (nAtomsAssigned == 0) {
      // Assign atom with the largest contribution.
      unsigned int maxAtom;
      occupiedMullikenPops.col(iOcc).array().abs().maxCoeff(&maxAtom);
      tripletList.push_back(Eigen::Triplet<int>(maxAtom, iOcc, 1));
      OutputControl::dOut << "Fitting domain selection lead to tiny fitting domain for orbital " << iOcc << "." << std::endl;
    }
  } // for iOcc
  toReturn.setFromTriplets(tripletList.begin(), tripletList.end());
  return toReturn;
}

const SparseMap& SparseMapsController::getOccToAtomMap() {
  if (!_occupiedOrbitalToAtomMap) {
    _occupiedOrbitalToAtomMap = std::make_shared<SparseMap>(buildOccToAtomMap(_orbitalWiseMullikenThresholds));
  }
  return *_occupiedOrbitalToAtomMap;
}

const SparseMap& SparseMapsController::getAtomToPAOMap() {
  if (!_atomToPAOMap) {
    unsigned int nPAOs = _virtualCoefficients->cols();
    unsigned int nAtoms = _system->getGeometry()->getNAtoms();
    Eigen::VectorXd atomWiseMullikenThresholds = convertToAtomWiseThresholds(_orbitalWiseMullikenThresholds);
    auto restMullikenGrossPops = MullikenPopulationCalculator<RESTRICTED>::calculateAtomwiseOrbitalPopulations(
        *_virtualCoefficients, _system->getOneElectronIntegralController()->getOverlapIntegrals(),
        _system->getAtomCenteredBasisController()->getBasisIndices());
    _atomToPAOMap = std::make_shared<SparseMap>(nPAOs, nAtoms);
    std::vector<Eigen::Triplet<int>> tripletList;
    for (unsigned int iPAO = 0; iPAO < nPAOs; ++iPAO) {
      // Use the mulliken charges to assign atom indices to the orbitals for integral prescreening.
      for (unsigned int atomIndex = 0; atomIndex < restMullikenGrossPops.rows(); ++atomIndex) {
        if (std::fabs(restMullikenGrossPops(atomIndex, iPAO)) > atomWiseMullikenThresholds(atomIndex)) {
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
    const Eigen::MatrixXd& coefficients = *_occupiedCoefficients;
    _shellToOccMap = constructShellToOrbitalMap(coefficients, _orbitalToShellThresholds);
  }
  return *_shellToOccMap;
}

const SparseMap& SparseMapsController::getShellToPAOMap() {
  if (!_shellToPAOMap) {
    const SparseMap paoToAtomMap = this->getAtomToPAOMap().transpose().eval();
    const unsigned int nPAOs = paoToAtomMap.cols();
    Eigen::VectorXd paoWiseThresholds = Eigen::VectorXd::Constant(nPAOs, 1e-3);
    const Eigen::MatrixXd atomWiseShellThresholds = convertToAtomWiseThresholds(_orbitalToShellThresholds);
    for (unsigned int iPAO = 0; iPAO < nPAOs; ++iPAO) {
      double minThreshold = 1e-3;
      for (Eigen::SparseMatrix<int>::InnerIterator itAtom(paoToAtomMap, iPAO); itAtom; ++itAtom) {
        minThreshold = std::min(minThreshold, atomWiseShellThresholds(itAtom.row()));
      }
      paoWiseThresholds(iPAO) = minThreshold;
    }
    _shellToPAOMap = constructShellToOrbitalMap(*_virtualCoefficients, paoWiseThresholds);
  }
  return *_shellToPAOMap;
}

const SparseMap& SparseMapsController::getAtomToAuxShellMap() {
  if (!_atomToAuxShellMap) {
    auto auxBasis = _system->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
    unsigned int nAtoms = _system->getGeometry()->getNAtoms();
    _atomToAuxShellMap = std::make_shared<SparseMap>(auxBasis->getReducedNBasisFunctions(), nAtoms);
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
      const SparseMap occToAtom = buildOccToAtomMap(_weakTripletMullikenScaling * _orbitalWiseMullikenThresholds);
      _weakTripletOccToKMap = std::make_shared<SparseMap>(this->getAtomToAuxShellMap().rows(), occToAtom.cols());
      *_weakTripletOccToKMap = (this->getAtomToAuxShellMap() * occToAtom).pruned().eval();
    } // if !_weakTripletOccToKMap
    return *_weakTripletOccToKMap;
  }
  else {
    if (!_strongTripletOccToKMap) {
      const SparseMap occToAtom = buildOccToAtomMap(_strongTripletMullikenScaling * _orbitalWiseMullikenThresholds);
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

const SparseMap& SparseMapsController::getExtendedOccToAuxShellMap_triples() {
  if (!_extendedOccToK_triples) {
    // Ensure existence of the weak occ->K map.
    this->getTripletOccToAuxShellMap(true);
    // Parse weak and strong map to the extended map construction.
    _extendedOccToK_triples =
        buildTriplesExtendedMap(this->getTripletOccToAuxShellMap(false), this->getTripletOccToAuxShellMap(true));
  }
  return *_extendedOccToK_triples;
}

const SparseMap& SparseMapsController::getExtendedOccToPAOMap() {
  if (!_extendedOccToPAO) {
    _extendedOccToPAO = buildExtendedMap(this->getOccToPAOMap());
  }
  return *_extendedOccToPAO;
}

const SparseMap& SparseMapsController::getExtendedOccToPAOMap_triples() {
  if (!_extendedOccToPAO_triples) {
    _extendedOccToPAO_triples = buildTriplesExtendedMap(this->getOccToPAOMap());
  }
  return *_extendedOccToPAO_triples;
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

const SparseMap& SparseMapsController::getExtendedKtoRhoMap_triples() {
  if (!_extendedAuxShellToRho_triples) {
    const SparseMap& shellToOcc = this->getShellToOccMap();
    const SparseMap& extOccToAux = this->getExtendedOccToAuxShellMap_triples();
    _extendedAuxShellToRho_triples = std::make_shared<SparseMap>(shellToOcc.cols(), extOccToAux.rows());
    *_extendedAuxShellToRho_triples = SparseMap((extOccToAux * shellToOcc).transpose()).pruned().eval();
  }
  return *_extendedAuxShellToRho_triples;
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

const SparseMap& SparseMapsController::getExtendedKtoPAOMap_triples() {
  if (!_extendedAuxShellToPAO_triples) {
    const SparseMap& extOccToAux = this->getExtendedOccToAuxShellMap_triples();
    const SparseMap& occToPAOMap = this->getExtendedOccToPAOMap_triples();
    _extendedAuxShellToPAO_triples = std::make_shared<SparseMap>(occToPAOMap.rows(), extOccToAux.rows());
    *_extendedAuxShellToPAO_triples = SparseMap((extOccToAux * occToPAOMap.transpose()).transpose()).pruned().eval();
  }
  return *_extendedAuxShellToPAO_triples;
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
const SparseMap& SparseMapsController::getExtendedKtoSigmaMap_triples() {
  if (!_extendedAuxShellToSigma_triples) {
    const SparseMap& extendedKtoPAOMap = getExtendedKtoPAOMap_triples();
    const SparseMap& shellToPAO = this->getShellToPAOMap();
    _extendedAuxShellToSigma_triples = std::make_shared<SparseMap>(shellToPAO.cols(), extendedKtoPAOMap.rows());
    *_extendedAuxShellToSigma_triples = SparseMap((extendedKtoPAOMap.transpose() * shellToPAO).transpose()).pruned().eval();
  }
  return *_extendedAuxShellToSigma_triples;
}

std::shared_ptr<SparseMap> SparseMapsController::constructShellToOrbitalMap(const Eigen::MatrixXd& coefficients,
                                                                            const Eigen::VectorXd& thresholds) {
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
      int nImportant = (squaredC.block(idx[iAtom].first, i, nFunc, 1).array() > thresholds(i)).count();
      if (nImportant > 0) {
        for (unsigned int iShell = idxRed[iAtom].first; iShell < idxRed[iAtom].second; ++iShell) {
          tripletList.emplace_back(i, iShell, 1);
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
  *newMap += initialMap;
  for (const auto& pair : _closeOrbitalPairs) {
    Eigen::SparseVector<int> sumVector = (initialMap.col(pair->i) + initialMap.col(pair->j)).eval();
    for (auto kSet : pair->coupledPairs)
      sumVector += initialMap.col(kSet->getK()).eval();
    newMap->col(pair->i) += sumVector.eval();
    newMap->col(pair->j) += sumVector.eval();
  }
  newMap->pruned();
  return newMap;
}

std::shared_ptr<SparseMap> SparseMapsController::buildExtendedMap(const SparseMap& initialMap) {
  auto newMap = std::make_shared<SparseMap>(initialMap.rows(), initialMap.cols());
  *newMap += initialMap;
  for (const auto& pair : _orbitalPairs) {
    Eigen::SparseVector<int> sumVector = (initialMap.col(pair->i) + initialMap.col(pair->j)).eval();
    if (_klListExtension) {
      for (const auto& klSet : pair->klPairSets) {
        newMap->col(klSet->getKLPair()->i) += sumVector.eval();
        newMap->col(klSet->getKLPair()->j) += sumVector.eval();
      }
    }
    for (auto kSet : pair->coupledPairs)
      sumVector += initialMap.col(kSet->getK()).eval();
    newMap->col(pair->i) += sumVector.eval();
    newMap->col(pair->j) += sumVector.eval();
  }
  newMap->pruned();
  return newMap;
}

std::shared_ptr<SparseMap> SparseMapsController::buildTriplesExtendedMap(const Eigen::SparseMatrix<int>& initialMap) {
  if (_triples.size() < 1)
    throw SerenityError(
        (std::string) "ERROR: This SparseMapsController is unable to construct extended triples maps!\n" +
        "       Triples were never parsed!");
  auto newMap = std::make_shared<SparseMap>(initialMap.rows(), initialMap.cols());
  *newMap += initialMap;
  for (const auto& triple : _triples) {
    Eigen::SparseVector<int> sumVector =
        (initialMap.col(triple->getI()) + initialMap.col(triple->getJ()) + initialMap.col(triple->getK())).eval();
    newMap->col(triple->getI()) += sumVector.eval();
    newMap->col(triple->getJ()) += sumVector.eval();
    newMap->col(triple->getK()) += sumVector.eval();
  }
  newMap->pruned();
  return newMap;
}

std::shared_ptr<SparseMap> SparseMapsController::buildTriplesExtendedMap(const Eigen::SparseMatrix<int>& strongMap,
                                                                         const Eigen::SparseMatrix<int>& weakMap) {
  if (_triples.size() < 1)
    throw SerenityError(
        (std::string) "ERROR: This SparseMapsController is unable to construct extended triples maps!\n" +
        "       Triples were never parsed!");
  auto newMap = std::make_shared<SparseMap>(strongMap.rows(), strongMap.cols());
  *newMap += strongMap;
  *newMap += weakMap;
  for (const auto& triple : _triples) {
    bool notWeak = not triple->isWeak();
    Eigen::SparseVector<int> sumVector =
        (notWeak) ? strongMap.col(triple->getI()) + strongMap.col(triple->getJ()) + strongMap.col(triple->getK())
                  : weakMap.col(triple->getI()) + weakMap.col(triple->getJ()) + weakMap.col(triple->getK());
    newMap->col(triple->getI()) += sumVector.eval();
    newMap->col(triple->getJ()) += sumVector.eval();
    newMap->col(triple->getK()) += sumVector.eval();
  }
  newMap->pruned();
  return newMap;
}

Eigen::VectorXd SparseMapsController::convertToAtomWiseThresholds(Eigen::VectorXd orbitalWiseThresholds) {
  const SparseMap& occToAtomMap = this->getOccToAtomMap(); // nAtoms x nOcc
  const SparseMap atomToOccMap = occToAtomMap.transpose(); // nOcc x nAtoms
  Eigen::VectorXd atomWiseThresholds = Eigen::VectorXd::Zero(atomToOccMap.cols());
  // Loop Atoms. Select tightest threshold for each atom!
  for (unsigned int iAtom = 0; iAtom < atomToOccMap.cols(); ++iAtom) {
    double threshold = 1;
    for (Eigen::SparseMatrix<int>::InnerIterator itIOcc(atomToOccMap, iAtom); itIOcc; ++itIOcc) {
      const double testThreshold = orbitalWiseThresholds(itIOcc.row());
      if (threshold > testThreshold)
        threshold = testThreshold;
    } // for itIOcc
    atomWiseThresholds(iAtom) = threshold;
  } // for iAtom
  return atomWiseThresholds;
}

} /* namespace Serenity */
