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
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "io/Filesystem.h"
#include "io/HDF5.h"
#include "misc/SerenityError.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"
#include "system/System.h"
#include "system/SystemController.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <vector>

namespace Serenity {

Eigen::MatrixXd AtomicDensityGuessCalculator::performAtomInitialGuess(Settings settings, std::shared_ptr<Atom> atom) {
  Eigen::MatrixXd atomDensMat;
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
  settings.scf.initialguess = Options::INITIAL_GUESSES::H_CORE;
  settings.scf.energyThreshold = 1e-6;
  settings.scf.rmsdThreshold = 1e-6;
  settings.scf.diisThreshold = 1e-6;
  settings.scf.damping = Options::DAMPING_ALGORITHMS::STATIC;
  settings.scf.staticDampingFactor = 0.7;
  settings.scf.useOffDiagLevelshift = false;
  settings.scf.allowNotConverged = true;
  settings.scf.degeneracyThreshold = 0.1;
  settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::NONE;
  settings.extCharges.externalChargesFile = "";

  std::vector<std::shared_ptr<Atom>> atomVec = {atom};
  auto atomGeom = std::make_shared<Geometry>(atomVec);
  auto atomSys = std::make_shared<SystemController>(atomGeom, settings);
  if (settings.spin == 1) {
    ScfTask<UNRESTRICTED> scf(atomSys);
    scf.run();
    atomDensMat = atomSys->getElectronicStructure<UNRESTRICTED>()->getDensityMatrix().total();
  }
  else {
    ScfTask<RESTRICTED> scf(atomSys);
    scf.run();
    atomDensMat = atomSys->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
  }
  return atomDensMat;
}

std::unique_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>>
AtomicDensityGuessCalculator::calculateInitialDensity(std::shared_ptr<SystemController> systemController) {
  if (!systemController)
    throw SerenityError("AtomicDensityGuessCalculator does not have a proper systemController!");

  std::string pathToGuessDensities;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToGuessDensities = env_p;
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  if (systemController->getSettings().basis.makeSphericalBasis) {
    pathToGuessDensities += "initialGuess/spherical/";
  }
  else {
    pathToGuessDensities += "initialGuess/cartesian/";
  }

  const auto& atoms = systemController->getAtoms();
  auto basisController = systemController->getBasisController();

  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::overlap, 0, 2);
  // these two variables are needed in the outer scope
  Eigen::MatrixXd targetBasisOverlap = systemController->getOneElectronIntegralController()->getOverlapIntegrals();

  auto minimalBasisController = systemController->getAtomCenteredBasisController(Options::BASIS_PURPOSES::SCF_DENS_GUESS);
  auto newGuessDensity = std::make_unique<DensityMatrix<Options::SCF_MODES::RESTRICTED>>(basisController);
  auto& newDens = (*newGuessDensity);
  unsigned int blockstart = 0;
  for (auto atom : atoms) {
    if (atom->isDummy()) {
      for (auto& shell : atom->getBasisFunctions()) {
        blockstart += shell->getNContracted();
      }
      continue;
    }
    bool usesECPs = atom->usesECP();

    if (_atomDensities.find(atom->getAtomType()->getElementSymbol()) == _atomDensities.end()) {
      Eigen::MatrixXd atomDensMat;
      if (_scf == GUESSMODES::SCF_INPLACE || usesECPs) {
        atomDensMat = performAtomInitialGuess(systemController->getSettings(), atom);
      }
      else {
        try {
          HDF5::Filepath name(pathToGuessDensities + atom->getAtomType()->getElementSymbol() + ".dmat.res.h5");
          HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
          HDF5::dataset_exists(file, "densityMatrix");
          HDF5::load(file, "densityMatrix", atomDensMat);
          file.close();
        }
        catch (...) {
          WarningTracker::printWarning((std::string) "  Warning: Was not able to read " + pathToGuessDensities +
                                           atom->getAtomType()->getElementSymbol() +
                                           ".dmat.res.h5\n"
                                           "  Performing an atom SCF now.",
                                       true);
          atomDensMat = performAtomInitialGuess(systemController->getSettings(), atom);
        }

        std::vector<std::shared_ptr<Atom>> dummyVec = {atom};
        auto dummyGeom = std::make_shared<Geometry>(dummyVec);

        // Get projection matrix minbasis -> actual basis.
        auto atomMinBas =
            AtomCenteredBasisControllerFactory::produce(dummyGeom, systemController->getSettings().basis.basisLibPath,
                                                        systemController->getSettings().basis.makeSphericalBasis, false,
                                                        systemController->getSettings().basis.firstECP, "MINAO");
        auto atomTargetBas = AtomCenteredBasisControllerFactory::produce(
            dummyGeom, systemController->getSettings().basis.basisLibPath,
            systemController->getSettings().basis.makeSphericalBasis, false,
            systemController->getSettings().basis.firstECP, systemController->getSettings().basis.label);
        unsigned int nBasFunc = atomTargetBas->getNBasisFunctions();
        Eigen::MatrixXd overlapB = targetBasisOverlap.block(blockstart, blockstart, nBasFunc, nBasFunc);

        auto overlapAB = libint.compute1eInts(LIBINT_OPERATOR::overlap, atomMinBas, atomTargetBas);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(overlapB, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.setThreshold(1e-6);

        Eigen::MatrixXd projectionOperator = svd.solve(overlapAB);
        atomDensMat = projectionOperator * atomDensMat * projectionOperator.transpose();
      }

      // Keep matrix for this atom type.
      _atomDensities[atom->getAtomType()->getElementSymbol()] = atomDensMat;
    }

    // Get atom density matrix and fill systems density matrix
    auto atomDensMat = _atomDensities[atom->getAtomType()->getElementSymbol()];
    unsigned int nBasFunc = atomDensMat.rows();
    newDens.block(blockstart, blockstart, nBasFunc, nBasFunc) = atomDensMat;
    blockstart += nBasFunc;
  }

  libint.freeEngines(LIBINT_OPERATOR::overlap, 0, 2);

  return newGuessDensity;
}
} /* namespace Serenity */
