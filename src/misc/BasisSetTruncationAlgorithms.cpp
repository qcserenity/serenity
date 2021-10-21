/**
 * @file BasisSetTruncationAlgorithms.cpp
 *
 * @date Oct 30, 2017
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
#include "misc/BasisSetTruncationAlgorithms.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "basis/Shell.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrixController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutputStream.h" //Filtered output streams.
#include "misc/WarningTracker.h"
#include "settings/MiscOptions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <assert.h>
#include <algorithm>
#include <vector>

namespace Serenity {

BasisSetTruncationAlgorithms::BasisSetTruncationAlgorithms(std::shared_ptr<SystemController> activeSystem)
  : _activeSystem(activeSystem) {
}

Eigen::VectorXi BasisSetTruncationAlgorithms::netPopulationTruncation(const DensityMatrix<Options::SCF_MODES::RESTRICTED> dMatActiveSystem,
                                                                      double netPopThreshold) {
  /*
   * Keeping all AOs centered on an active atom and all environmental AOs "a" with a Mulliken net
   * population q_a >= netPopTheshold calculated from the density matrix of the active system.
   *
   * q_a = P_aa^Active * S_aa
   *
   * where P_aa is the density matrix entry aa of the active density matrix and S_aa the corresponding overlap matrix
   * entry.
   */
  if (netPopThreshold <= 0.0) {
    WarningTracker::printWarning("WARNING: You selected the 'net population truncation' with a threshold <= 0.",
                                 iOOptions.printSCFCycleInfo);
  }

  unsigned int nBFstot = dMatActiveSystem.cols();
  unsigned int nShellsTot = _activeSystem->getBasisController()->getReducedNBasisFunctions();
  Eigen::VectorXd importance = Eigen::VectorXd::Zero(nBFstot);
  auto sMatSuper(_activeSystem->getOneElectronIntegralController()->getOverlapIntegrals());
  // generate measures for importance of the environment AOs (for more information see above.)
  // Basis functions are organized like: ...active...environmental... in the array
  Eigen::VectorXi importantShells = Eigen::VectorXi::Zero(nShellsTot);

  auto idxRed = _activeSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
  auto idxNotRed = _activeSystem->getAtomCenteredBasisController()->getBasisIndices();
  auto actAtoms = _activeSystem->getGeometry()->getAtoms();
  unsigned int atomIndex = 0;
  for (const auto& atom : actAtoms) {
    // All non-dummy shells are important.
    if (not atom->isDummy()) {
      for (unsigned int j = idxRed[atomIndex].first; j < idxRed[atomIndex].second; j++) {
        importantShells(j) = true;
      } // for j
    }   // if not atom->isDummy()
    else {
      // Use net-population criterion for dummy atom shells.
      for (unsigned int mu = idxNotRed[atomIndex].first; mu < idxNotRed[atomIndex].second; ++mu) {
        double netPop = dMatActiveSystem(mu, mu) * sMatSuper(mu, mu);
        if (netPop >= netPopThreshold) {
          importantShells[_activeSystem->getBasisController()->reducedIndex(mu)] = 1;
        }
      }
    }
    ++atomIndex;
  } // for atom
  return importantShells;
}

Eigen::VectorXi BasisSetTruncationAlgorithms::primitiveNetPopTruncation(const DensityMatrix<Options::SCF_MODES::RESTRICTED> dMatActiveSystem,
                                                                        double truncationFactor) {
  assert(truncationFactor >= 0.0 && "A negative ratio of kept basis functions does not make any sense.");
  unsigned int nBFstot = dMatActiveSystem.cols();
  unsigned int nShellsTot = _activeSystem->getBasisController()->getReducedNBasisFunctions();
  Eigen::VectorXd importance = Eigen::VectorXd::Zero(nBFstot);
  auto sMatSuper(_activeSystem->getOneElectronIntegralController()->getOverlapIntegrals());

  // generate measures for importance of the environment AOs
  for (unsigned int i = 0; i < nBFstot; i++) {
    importance[i] = dMatActiveSystem(i, i) * sMatSuper(i, i);
  }
  // Select all non-dummy shells as important.
  Eigen::VectorXi importantShells = Eigen::VectorXi::Zero(nShellsTot);
  auto idxRed = _activeSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
  unsigned int nNonDummyBasisFunctions = 0;
  auto actAtoms = _activeSystem->getGeometry()->getAtoms();
  unsigned int atomIndex = 0;
  for (const auto& atom : actAtoms) {
    // All non-dummy shells are important.
    if (not atom->isDummy()) {
      auto nBasisFunctions = 0;
      const auto& shells = atom->getBasisFunctions();
      for (const auto& shell : shells)
        nBasisFunctions += shell->getNContracted();
      nNonDummyBasisFunctions += nBasisFunctions;
      for (unsigned int j = idxRed[atomIndex].first; j < idxRed[atomIndex].second; j++) {
        importantShells(j) = true;
      } // for j
    }   // if not atom->isDummy()
    ++atomIndex;
  } // for atom

  // group importance per Shell
  // keeping only (environment) basis functions associated to a shell
  // with at least one basis function being in the
  // top 1-truncation most "important" environment basis functions.
  unsigned int i = 0;
  while (i < (unsigned int)((nBFstot - nNonDummyBasisFunctions) * (1.0 - truncationFactor))) {
    int imp;
    importance.maxCoeff(&imp);
    importance[imp] = 0.0;
    auto reducedIndex = _activeSystem->getBasisController()->reducedIndex(imp);
    if (not importantShells[reducedIndex]) {
      importantShells[reducedIndex] = 1;
      ++i;
    }
  } // while i<(unsigned int)((nBFstot-_nCoreBasisFunctionsActive)*(1.0-truncationFactor))
  return importantShells;
}

void BasisSetTruncationAlgorithms::truncateBasis(Options::BASIS_SET_TRUNCATION_ALGORITHMS truncationAlgorithm,
                                                 double truncationFactor,
                                                 std::shared_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>> dMatActiveSystem,
                                                 double netPopThreshold) {
  Eigen::VectorXi importantShells;
  switch (truncationAlgorithm) {
    case Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE:
      assert(false && "Do not use a truncation algorithm if you do not want to truncate anything!");
      break;

    case Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION:
      assert(dMatActiveSystem &&
             "The 'net population' truncation was called without the density matrix of the active system.");
      importantShells = this->netPopulationTruncation(*dMatActiveSystem, netPopThreshold);
      break;

    case Options::BASIS_SET_TRUNCATION_ALGORITHMS::PRIMITIVE_NET_POPULATION:
      assert(dMatActiveSystem &&
             "The 'primitive net population' truncation was called without the density matrix of the active system.");
      importantShells = this->primitiveNetPopTruncation(*dMatActiveSystem, truncationFactor);
      break;

    default:
      assert(false && "Unknown truncation algorithm.");
  }
  // build basis
  this->buildTruncatedBasis(importantShells);
}

/* ==============================================================
 *                       Helper Functions
 * ============================================================== */
inline void BasisSetTruncationAlgorithms::buildTruncatedBasis(Eigen::VectorXi importantShells) {
  // get active system and supersystem atoms
  auto activeAtoms = _activeSystem->getGeometry()->getAtoms();
  int nBFstot = _activeSystem->getBasisController()->getNBasisFunctions();
  auto idx = _activeSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
  // loop over all environmental atoms
  std::vector<Eigen::VectorXi> activeShells;
  std::vector<bool> deleteAtom = std::vector<bool>(activeAtoms.size(), true);
  unsigned int atomIndex = 0;
  for (const auto& atom : activeAtoms) {
    Eigen::VectorXi tmp;
    if (atom->isDummy()) {
      tmp = Eigen::VectorXi::Zero(idx[atomIndex].second - idx[atomIndex].first);
      // generate list of shells for this atom
      for (unsigned int j = idx[atomIndex].first; j < idx[atomIndex].second; j++) {
        if (importantShells[j]) {
          tmp[j - idx[atomIndex].first] = true;
        }
      }
    } /* atom->isDummy() */
    else {
      tmp = Eigen::VectorXi::Ones(idx[atomIndex].second - idx[atomIndex].first);
    }
    if (tmp.sum() > 0) {
      activeShells.push_back(tmp);
      deleteAtom[atomIndex] = false;
    }
    else {
      OutputControl::vOut << "NOTE: Deleting dummy atom! No shells left. Atom index " << atomIndex + 1 << std::endl;
    }
    ++atomIndex;
  } // for i

  // Delete all obsolete atoms of the old geometry.
  unsigned int nDeleted = 0;
  for (unsigned int iAtom = 0; iAtom < activeAtoms.size(); ++iAtom) {
    if (deleteAtom[iAtom]) {
      _activeSystem->getGeometry()->deleteAtom(iAtom - nDeleted);
      ++nDeleted;
    }
  }
  _activeSystem->getGeometry()->printToFile(_activeSystem->getHDF5BaseName(), _activeSystem->getSystemIdentifier());
  // set new basis
  std::string label = _activeSystem->getAtomCenteredBasisController()->getBasisLabel();
  auto newBas =
      std::make_shared<AtomCenteredBasisController>(_activeSystem->getGeometry(), _activeSystem->getSettings().basis.basisLibPath,
                                                    _activeSystem->getSettings().basis.makeSphericalBasis, true, label,
                                                    999999, // No new ECP assignment
                                                    activeShells);
  newBas->toHDF5(_activeSystem->getHDF5BaseName(), _activeSystem->getSystemIdentifier());
  auto newAuxBas = std::make_shared<AtomCenteredBasisController>(
      _activeSystem->getGeometry(), _activeSystem->getSettings().basis.basisLibPath,
      _activeSystem->getSettings().basis.makeSphericalBasis, false, _activeSystem->getSettings().basis.auxJLabel);
  _activeSystem->setBasisController(newBas);
  _activeSystem->setBasisController(newAuxBas, Options::BASIS_PURPOSES::AUX_COULOMB);
  int nBFsAct = _activeSystem->getBasisController()->getNBasisFunctions();

  OutputControl::nOut << "Supersystem basis: " << nBFstot << " basis functions" << std::endl;
  OutputControl::nOut << "Active basis: " << nBFsAct << " basis functions" << std::endl;
  OutputControl::nOut << "Truncated supersystem basis by " << nBFstot - nBFsAct << " basis functions." << std::endl
                      << std::endl;
}

} /* namespace Serenity */
