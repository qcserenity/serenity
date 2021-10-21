/**
 * @file AtomicDensityGuessCalculator.cpp
 *
 * @date Jul 12, 2014
 * @author Thomas Dresselhaus, David Schnieders
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
#include "scf/initialGuess/AtomicDensityGuessCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "basis/Basis.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "io/Filesystem.h"
#include "io/HDF5.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"
#include "system/System.h"
#include "system/SystemController.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <fstream>
#include <vector>

namespace Serenity {

Eigen::MatrixXd AtomicDensityGuessCalculator::performAtomInitialGuess(Settings settings, std::shared_ptr<Atom> atom) {
  Eigen::MatrixXd atomDensMat;
  bool usesECPs = atom->usesECP();
  // Copy system settings
  settings.name = atom->getAtomType()->getElementSymbol() + "_FREE";
  // Uncharged atoms
  settings.charge = 0;
  if (atom->getNuclearCharge() % 2 != 0) {
    settings.spin = 1;
  }
  else {
    settings.spin = 0;
  }
  settings.scf.initialguess = (!usesECPs) ? Options::INITIAL_GUESSES::ATOM_SCF : Options::INITIAL_GUESSES::H_CORE;
  settings.scf.energyThreshold = 1e-6;
  settings.scf.rmsdThreshold = 1e-6;
  settings.scf.diisThreshold = 1e-6;
  settings.scf.damping = Options::DAMPING_ALGORITHMS::STATIC;
  settings.scf.staticDampingFactor = 0.7;
  settings.scf.useOffDiagLevelshift = false;
  settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::NONE;
  std::vector<std::shared_ptr<Atom>> atomVec = {atom};
  auto atomGeom = std::make_shared<Geometry>(atomVec);
  auto atomSys = std::make_shared<SystemController>(atomGeom, settings);
  if (settings.spin == 1) {
    ScfTask<UNRESTRICTED> scf(atomSys);
    scf.settings.fractionalDegeneracy = true;
    scf.run();
    atomDensMat = atomSys->getElectronicStructure<UNRESTRICTED>()->getDensityMatrix().total();
  }
  else {
    ScfTask<RESTRICTED> scf(atomSys);
    scf.settings.fractionalDegeneracy = true;
    scf.run();
    atomDensMat = atomSys->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  }
  return atomDensMat;
}

std::unique_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>>
AtomicDensityGuessCalculator::calculateInitialDensity(std::shared_ptr<SystemController> systemController,
                                                      bool keepMinimalBasis, bool scale) {
  assert(systemController);

  std::string pathToGuessDensities;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToGuessDensities = env_p;
  }
  else {
    std::cout << "ERROR Environment variable SERENITY_RESOURCES not set." << std::endl;
    assert(false);
  }
  if (systemController->getSettings().basis.makeSphericalBasis) {
    pathToGuessDensities += "initialGuess/spherical/";
  }
  else {
    pathToGuessDensities += "initialGuess/cartesian/";
  }
  if (_scf == GUESSMODES::SCF) {
    if (systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      pathToGuessDensities += "dft_scf/";
    }
    else {
      pathToGuessDensities += "hf_scf/";
    }
  }
  else {
    pathToGuessDensities += "ao_occupation/";
  }

  /*
   * Get in some data locally
   */
  const auto& atoms = systemController->getAtoms();
  std::shared_ptr<BasisController> minimalBasisController;
  if (_scf != GUESSMODES::SCF_INPLACE) {
    minimalBasisController = systemController->getAtomCenteredBasisController(
        (_scf == GUESSMODES::SCF) ? Options::BASIS_PURPOSES::SCF_DENS_GUESS : Options::BASIS_PURPOSES::MINBAS);
  }

  auto basisController = systemController->getBasisController();

  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::overlap, 0, 2);
  // these two variables are needed in the outer scope
  Eigen::MatrixXd targetBasisOverlap = systemController->getOneElectronIntegralController()->getOverlapIntegrals();

  auto newGuessDensity = std::unique_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>>(
      new DensityMatrix<Options::SCF_MODES::RESTRICTED>(keepMinimalBasis ? minimalBasisController : basisController));
  auto& newDens = *newGuessDensity;
  unsigned int blockstart = 0;
  for (auto atom : atoms) {
    if (atom->isDummy()) {
      for (auto& shell : atom->getBasisFunctions()) {
        blockstart += shell->getNContracted();
      }
      continue;
    }
    bool usesECPs = atom->usesECP();
    /*
     * If not already available: Get from file
     */
    if (_atomDensities.find(atom->getAtomType()->getElementSymbol()) == _atomDensities.end()) {
      Eigen::MatrixXd atomDensMat;
      if (_scf == GUESSMODES::SCF_INPLACE || usesECPs) {
        std::string freeAtomName = atom->getAtomType()->getElementSymbol() + "_FREE";
        std::string atomFilesPath = systemController->getSettings().path + freeAtomName + "/";
        std::string pathRes = atomFilesPath + freeAtomName + ".dmat.res.h5";
        std::string pathUnres = atomFilesPath + freeAtomName + ".dmat.unres.h5";
        std::string settingsPath = atomFilesPath + freeAtomName + ".settings";

        // Check if inital guess files already exist
        struct stat buffer;
        const bool resExists = stat(pathRes.c_str(), &buffer) == 0;
        const bool unresExists = (stat(pathUnres.c_str(), &buffer) == 0);
        const bool settingsExist = stat(settingsPath.c_str(), &buffer) == 0;
        if ((resExists || unresExists) && settingsExist) {
          // Check the settings of the old files.
          // Constuct the settings from file.
          std::ifstream settingsFile(settingsPath);
          Settings oldSettings(settingsFile);
          settingsFile.close();
          const auto& newSettings = systemController->getSettings();
          // Check if every thing as reqiured.
          const bool labelIsFine = oldSettings.basis.label == newSettings.basis.label;
          const bool angularIsFine = oldSettings.basis.makeSphericalBasis == newSettings.basis.makeSphericalBasis;
          const bool basisSetLibIsFine = oldSettings.basis.basisLibPath == newSettings.basis.basisLibPath;
          if (not labelIsFine || not angularIsFine || not basisSetLibIsFine) {
            // Something did not match. Reconstruct!
            WarningTracker::printWarning(
                "WARNING: Old initial guess files were found on diks, but the basis label/library\n"
                "         or angular momentum mode did not match the required settings!\n"
                "         The old initial guess files will be deleted and replaced!",
                true);
            removeSystemFiles(atomFilesPath, freeAtomName);
            atomDensMat = performAtomInitialGuess(systemController->getSettings(), atom);
          }
          else {
            if (resExists) {
              // Load old files.
              HDF5::Filepath name(pathRes);
              HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
              HDF5::dataset_exists(file, "densityMatrix");
              HDF5::load(file, "densityMatrix", atomDensMat);
              file.close();
            }
            else {
              Eigen::MatrixXd dummy;
              HDF5::Filepath name(pathUnres);
              HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
              HDF5::dataset_exists(file, "densityMatrix_alpha");
              HDF5::load(file, "densityMatrix_alpha", dummy);
              atomDensMat = dummy;
              HDF5::dataset_exists(file, "densityMatrix_beta");
              HDF5::load(file, "densityMatrix_beta", dummy);
              atomDensMat += dummy;
              file.close();
            }
          }
        }
        // If initial guess does not already exists, run atom calculations
        else {
          atomDensMat = performAtomInitialGuess(systemController->getSettings(), atom);
        }
      }
      else {
        HDF5::Filepath name(pathToGuessDensities + atom->getAtomType()->getElementSymbol() + ".dmat.res.h5");
        HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        HDF5::dataset_exists(file, "densityMatrix");
        HDF5::load(file, "densityMatrix", atomDensMat);
        file.close();
      }
      /*
       * If basis change is wanted: Project to new basis
       */
      if (!keepMinimalBasis and _scf != GUESSMODES::SCF_INPLACE and !usesECPs) {
        std::vector<std::shared_ptr<Atom>> dummyVec = {atom};
        auto dummyGeom = std::make_shared<Geometry>(dummyVec);
        /*
         * Build Projection operation
         */
        auto atomMinBas = AtomCenteredBasisControllerFactory::produce(
            dummyGeom, systemController->getSettings().basis.basisLibPath,
            systemController->getSettings().basis.makeSphericalBasis, false,
            systemController->getSettings().basis.firstECP, (_scf == GUESSMODES::SCF) ? "DEF2-QZVP" : "STO-3G");
        auto atomTargetBas = AtomCenteredBasisControllerFactory::produce(
            dummyGeom, systemController->getSettings().basis.basisLibPath,
            systemController->getSettings().basis.makeSphericalBasis, false,
            systemController->getSettings().basis.firstECP, systemController->getSettings().basis.label);
        unsigned int nBasFunc = atomTargetBas->getNBasisFunctions();
        Eigen::MatrixXd overlapB = targetBasisOverlap.block(blockstart, blockstart, nBasFunc, nBasFunc);
        auto overlapAB = libint.compute1eInts(LIBINT_OPERATOR::overlap, atomMinBas, atomTargetBas);
        // Calculate inverse of AO overlap integrals in basis B. Use SVD,
        // i.e. S_B = U * D * V^T, since overlapB could be ill-conditioned
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(overlapB, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.setThreshold(1e-6);
        // Calculate projection operator P_{BA}
        Eigen::MatrixXd projectionOperator;
        projectionOperator = svd.solve(overlapAB);
        auto atomdens = (projectionOperator * atomDensMat * projectionOperator.transpose()).eval();
        const double factor = scale ? atom->getEffectiveCharge() / overlapB.cwiseProduct(atomdens).sum() : 1.0;
        atomDensMat = atomdens * factor;
      }
      /*
       * Save matrix for this atom
       */
      _atomDensities[atom->getAtomType()->getElementSymbol()] = atomDensMat;
    }
    /*
     * Get atom density matrix and fill systems density matrix
     */
    auto atomDensMat = _atomDensities[atom->getAtomType()->getElementSymbol()];
    unsigned int nBasFunc = atomDensMat.rows();
    newDens.block(blockstart, blockstart, nBasFunc, nBasFunc) = atomDensMat;
    blockstart += nBasFunc;
  } // for atom

  libint.freeEngines(LIBINT_OPERATOR::overlap, 0, 2);

  return newGuessDensity;
}
} /* namespace Serenity */
