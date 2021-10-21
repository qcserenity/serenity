/**
 * @file BasisExtension.cpp
 *
 * @date Jan 11, 2018
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
#include "basis/BasisExtension.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/matrices/MatrixInBasis.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <assert.h>

namespace Serenity {

void BasisExtension::extendAllSystems(std::vector<std::shared_ptr<SystemController>> systems, double overlapThreshold) {
  std::cout << "Extending the basis for all active systems" << std::endl << std::endl;
  /*
   * 1. Evaluate the how important the other shells are for all systems.
   * 2. Build the new atom vectors --> core-atoms + ghost-atoms
   * 3. Update the basis
   */
  /* ====================== */
  /*   1. Evaluate Shells   */
  /* ====================== */
  // activeSystem->envSystem->VectorWithImportantShellsPerAtom
  std::vector<std::vector<std::vector<Eigen::VectorXi>>> importantShellPerSystem;
  for (unsigned int i = 0; i < systems.size(); ++i) {
    auto act = systems[i];
    std::vector<std::vector<Eigen::VectorXi>> vec;
    importantShellPerSystem.push_back(vec);
    // Loop over environment systems
    for (unsigned int j = 0; j < systems.size(); ++j) {
      auto env = systems[j];
      std::vector<Eigen::VectorXi> importantShells;
      if (i == j) {
        auto actAtoms = act->getGeometry()->getAtoms();
        for (unsigned int iAtom = 0; iAtom < actAtoms.size(); ++iAtom) {
          Eigen::VectorXi impShellsOnAtom = Eigen::VectorXi::Ones(actAtoms[iAtom]->getBasisFunctions().size());
          importantShells.push_back(impShellsOnAtom);
        }
      }
      else {
        importantShells = getListOfImportantShells(env, act, overlapThreshold);
      }
      importantShellPerSystem[i].push_back(importantShells);
    } // j
  }   // i
  /* ================================= */
  /*   2. Build the new atom vectors   */
  /* ================================= */
  for (unsigned int i = 0; i < systems.size(); ++i) {
    auto activeSystem = systems[i];
    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "Extending basis of system: " << activeSystem->getSystemName() << std::endl;
    auto activeAtoms = buildActiveSystemAtoms(activeSystem);
    // copy!
    std::vector<Eigen::VectorXi> importantActiveShells = importantShellPerSystem[i][i];
    for (unsigned int j = 0; j < systems.size(); ++j) {
      if (i != j) {
        auto env = systems[j];
        auto importantShells = importantShellPerSystem[i][j];
        addGhostAtoms(importantShells, env, importantActiveShells, activeAtoms);
      }
    } // j
    _extendedAtomSets.push_back(activeAtoms);
    _importantShells.push_back(importantActiveShells);
  } // i
  /* ================================ */
  /*   3. Update Basis and Geometry   */
  /* ================================ */
  for (unsigned int i = 0; i < systems.size(); ++i) {
    auto activeSystem = systems[i];
    unsigned int nBFActOld = activeSystem->getBasisController()->getNBasisFunctions();
    const auto& activeAtoms = _extendedAtomSets[i];
    const auto& importantActiveShells = _importantShells[i];
    buildBasisFromGeometry(activeSystem, activeAtoms, importantActiveShells);
    unsigned int nBFAct = activeSystem->getBasisController()->getNBasisFunctions();
    std::cout << "Extended the active basis by " << nBFAct - nBFActOld << " basis functions" << std::endl << std::endl;
  }
}

/* ==================== */
/*   Helper Functions   */
/* ==================== */

inline void BasisExtension::addGhostAtoms(std::vector<Eigen::VectorXi> newImportantShells,
                                          std::shared_ptr<SystemController> environmentSystem,
                                          std::vector<Eigen::VectorXi>& activeImportantShells,
                                          std::vector<std::shared_ptr<Atom>>& activeAtoms) {
  // add dummy atoms which include the environment basis
  std::cout << "Shells from environment system " << environmentSystem->getSystemName() << std::endl;
  auto envAtoms = environmentSystem->getGeometry()->getAtoms();
  for (unsigned int iAtom = 0; iAtom < envAtoms.size(); ++iAtom) {
    const Eigen::VectorXi& impShellsOnAtom = newImportantShells[iAtom];
    // if there were important shells, generate a ghost atom and add the shells to it
    if (impShellsOnAtom.sum() > 0) {
      auto newGhost = std::make_shared<Atom>(envAtoms[iAtom]->getAtomType()->getElementSymbol() + ":",
                                             envAtoms[iAtom]->getX(), envAtoms[iAtom]->getY(), envAtoms[iAtom]->getZ());
      std::cout << "------------" << std::endl;
      std::cout << "Atom: " << envAtoms[iAtom]->getAtomType()->getElementSymbol() << std::endl;
      auto shells = envAtoms[iAtom]->getBasisFunctions();
      unsigned int nShells = shells.size();
      for (unsigned int iShell = 0; iShell < nShells; ++iShell) {
        if (impShellsOnAtom(iShell)) {
          auto shell = shells[iShell];
          std::string spherical = (!shell->isCartesian()) ? "YES" : "NO";
          std::cout << "l: " << shell->getAngularMomentum() << "      Spherical? " << spherical
                    << "    N contracted: " << shell->getNContracted() << std::endl;
        }
      }
      // add the new ghost atom to the list of  atoms for the active system.
      activeAtoms.push_back(newGhost);
      activeImportantShells.push_back(impShellsOnAtom);
    }
  } // iAtom
}

inline std::vector<Eigen::VectorXi>
BasisExtension::getListOfImportantShells(std::shared_ptr<SystemController> environmentSystem,
                                         std::shared_ptr<SystemController> activeSystem, double overlapThreshold) {
  // build the combined basis of the active and environment system
  // TODO: I expect the label of the basis to be identically for all systems.
  assert(activeSystem->getAtomCenteredBasisController()->getBasisLabel() ==
         environmentSystem->getAtomCenteredBasisController()->getBasisLabel());
  auto envBasisController = environmentSystem->getBasisController();
  auto& libint = Libint::getInstance();
  Eigen::MatrixXd s_AB =
      libint.compute1eInts(LIBINT_OPERATOR::overlap, envBasisController, activeSystem->getBasisController());
  // vector of important environment shells
  std::vector<Eigen::VectorXi> importantShells; // =
                                                // Eigen::VectorXi::Zero(envBasisController->getReducedNBasisFunctions())
  // loop over environment atoms
  const auto& envAtoms = environmentSystem->getGeometry()->getAtoms();
  auto idx = environmentSystem->getAtomCenteredBasisController()->getBasisIndices();
  auto idxRed = environmentSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
  for (unsigned int iAtom = 0; iAtom < envAtoms.size(); ++iAtom) {
    Eigen::VectorXi importantShellsOnAtom = Eigen::VectorXi::Zero(envAtoms[iAtom]->getBasisFunctions().size());
    // loop over basis functions centered on the atom
    unsigned int firstShell = idxRed[iAtom].first;
    for (unsigned int j = idx[iAtom].first; j < idx[iAtom].second; ++j) {
      double totalOverlap = s_AB.col(j).norm();
      if (totalOverlap > overlapThreshold) {
        unsigned int reducedIndex = envBasisController->reducedIndex(j) - firstShell;
        importantShellsOnAtom(reducedIndex) = 1;
      } // if overlap
    }   // for j
    importantShells.push_back(importantShellsOnAtom);
  } // for iAtom
  return importantShells;
}

inline std::vector<std::shared_ptr<Atom>> BasisExtension::buildActiveSystemAtoms(std::shared_ptr<SystemController> activeSystem) {
  // get active system atoms
  auto atoms(activeSystem->getGeometry()->getAtoms());
  // rename/add old basis as "EMB_TRUNK" basis
  for (auto& atom : atoms) {
    auto bfs = std::make_pair(std::string("EMB_EXTENDED"), atom->getBasisFunctions());
    atom->addBasis(bfs, true);
  }
  return atoms;
}

inline void BasisExtension::buildBasisFromGeometry(std::shared_ptr<SystemController> activeSystem,
                                                   std::vector<std::shared_ptr<Atom>> activeAtoms,
                                                   std::vector<Eigen::VectorXi> importantActiveShells) {
  // build the new active basis via the "new" geometry (active atoms + ghost atoms)
  // generate geometry
  auto geom = std::make_shared<Geometry>(activeAtoms);
  // set new basis
  auto newBas = std::make_shared<AtomCenteredBasisController>(
      geom, activeSystem->getSettings().basis.basisLibPath, activeSystem->getSettings().basis.makeSphericalBasis, true,
      activeSystem->getAtomCenteredBasisController()->getBasisLabel(), 37, importantActiveShells);
  auto newAuxBas = std::make_shared<AtomCenteredBasisController>(geom, activeSystem->getSettings().basis.basisLibPath,
                                                                 activeSystem->getSettings().basis.makeSphericalBasis,
                                                                 false, activeSystem->getSettings().basis.auxJLabel);
  *activeSystem->getGeometry() += Geometry(activeAtoms);
  activeSystem->getGeometry()->deleteIdenticalAtoms();
  activeSystem->getGeometry()->printToFile(activeSystem->getHDF5BaseName(), activeSystem->getSystemIdentifier());
  activeSystem->setBasisController(newBas);
  // Write the basis to disk. Enables restarts with the same basis.
  newBas->toHDF5(activeSystem->getHDF5BaseName(), activeSystem->getSystemIdentifier());
  activeSystem->setBasisController(newAuxBas, Options::BASIS_PURPOSES::AUX_COULOMB);
}

} /* namespace Serenity */
