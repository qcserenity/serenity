/**
 * @file   VirtualOrbitalSpaceSelectionTask.cpp
 *
 * @date   Aug 7, 2019
 * @author J. Toelle
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "tasks/VirtualOrbitalSpaceSelectionTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "io/FormattedOutputStream.h"
#include "io/HDF5.h"
#include "misc/VirtualOrbitalSelectionAlgorithms.h"
#include "parameters/Constants.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/ElectronicStructureCopyTask.h"
/* Include Std and External Headers */
#include <iomanip>
#include <numeric>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
VirtualOrbitalSpaceSelectionTask<SCFMode>::VirtualOrbitalSpaceSelectionTask(
    const std::vector<std::shared_ptr<SystemController>>& activeSystems,
    const std::vector<std::shared_ptr<SystemController>>& passiveSystems)
  : _act(activeSystems), _env(passiveSystems) {
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSpaceSelectionTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _act, _env);
  printSectionTitle("VirtualOrbitalSpaceSelectionTask");
  checkInput();
  if (settings.excludeProjection) {
    printf("  ------- Exclude Projection -------\n");
    // Some initial initialization
    CoefficientMatrix<SCFMode> coefs = _act[0]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenValues =
        _act[0]->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
    auto nOcc = _act[0]->template getNOccupiedOrbitals<SCFMode>();
    // The Virtual orbital space selector
    auto selectorAlgorithm = std::make_shared<VirtualOrbitalSelectionAlgorithms<SCFMode>>(_act[0], _env);
    // Creates Reference List of orbital indices
    SpinPolarizedData<SCFMode, std::vector<unsigned int>> indices;
    for_spin(indices) {
      indices_spin.resize(_act[0]->getBasisController()->getNBasisFunctions());
      std::iota(std::begin(indices_spin), std::end(indices_spin), 0);
    };
    selectorAlgorithm->excludeProjection(coefs, nOcc, indices);
    updateOrbitals(coefs, eigenValues, indices);
    // Update Orbitals
    _act[0]->template getActiveOrbitalController<SCFMode>()->updateOrbitals(coefs, eigenValues);
    selectorAlgorithm = nullptr;
    writeOrbitalsToHDF5(coefs, eigenValues, _act[0]);
  } /* Exclude Projection */

  if (settings.localCanonicalVirtuals != 0.0 || settings.envCanonicalVirtuals != 0.0) {
    printf("\n  ------- Choose Local Canonical Orbitals -------\n");
    _act[1]->getGeometry()->addDummy(*_act[0]->getGeometry());
    _act[1]->getGeometry()->deleteIdenticalAtoms();
    _act[1]->getGeometry()->print();
    _act[1]->getGeometry()->printToFile(_act[1]->getHDF5BaseName(), _act[1]->getSystemIdentifier());
    ElectronicStructureCopyTask<SCFMode> copytask(_act[0], {_act[1]});
    copytask.run();
    auto nOcc = _act[1]->template getNOccupiedOrbitals<SCFMode>();
    auto nVirt = _act[1]->template getNVirtualOrbitalsTruncated<SCFMode>();
    CoefficientMatrix<SCFMode> coeffs = _act[1]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenvals =
        _act[1]->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
    auto selectorAlgorithm = std::make_shared<VirtualOrbitalSelectionAlgorithms<SCFMode>>(_act[1], _env);
    auto es = _act[0]->template getElectronicStructure<SCFMode>();
    std::shared_ptr<FockMatrix<SCFMode>> fock = nullptr;
    if (es->checkFock() && settings.recalculateFockMatrix == false) {
      std::cout << "   Fock matrix from a previous embedding/supermolecular calculation used!" << std::endl;
      fock = std::make_shared<FockMatrix<SCFMode>>(es->getFockMatrix());
      _act[1]->template getElectronicStructure<SCFMode>()->setFockMatrix(*fock);
      _act[1]->template getElectronicStructure<SCFMode>()->toHDF5(_act[1]->getHDF5BaseName(), _act[1]->getSettings().identifier);
    }
    SpinPolarizedData<SCFMode, std::vector<unsigned int>> indices;
    for_spin(indices, nOcc) {
      indices_spin.resize(nOcc_spin);
      std::iota(std::begin(indices_spin), std::end(indices_spin), 0);
    };
    selectorAlgorithm->virtualCanonicalOrbitalSpaceSelection(coeffs, nOcc, nVirt, indices, settings.localCanonicalVirtuals,
                                                             settings.envCanonicalVirtuals, settings.onlyOne);
    updateOrbitals(coeffs, eigenvals, indices);
    // Update Orbitals
    _act[1]->template getActiveOrbitalController<SCFMode>()->updateOrbitals(coeffs, eigenvals);
    printNewOrbitals(nOcc, eigenvals);
    writeOrbitalsToHDF5(coeffs, eigenvals, _act[1]);
  } /* localCanonicalVirutalOrbitals */

  if (settings.localizedVirtualorbitals == true || settings.localizedEnvVirtualorbitals == true) {
    printf("\n  ------- Choose Localized Virtual Orbitals -------\n");
    _act[1]->getGeometry()->addDummy(*_act[0]->getGeometry());
    _act[1]->getGeometry()->deleteIdenticalAtoms();
    _act[1]->getGeometry()->print();
    _act[1]->getGeometry()->printToFile(_act[1]->getHDF5BaseName(), _act[1]->getSystemIdentifier());
    ElectronicStructureCopyTask<SCFMode> copytask(_act[0], {_act[1]});
    copytask.run();
    auto nOcc = _act[1]->template getNOccupiedOrbitals<SCFMode>();
    auto nVirt = _act[1]->template getNVirtualOrbitalsTruncated<SCFMode>();
    CoefficientMatrix<SCFMode> coeffs = _act[1]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenvals =
        _act[1]->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
    // Initialize VirtualOrbitalSelection
    auto selectorAlgorithm = std::make_shared<VirtualOrbitalSelectionAlgorithms<SCFMode>>(_act[1], _env);
    // Check if a Fock matrix exists and set Fockmatrix
    auto es = _act[0]->template getElectronicStructure<SCFMode>();
    std::shared_ptr<FockMatrix<SCFMode>> fock = nullptr;
    if (es->checkFock() && settings.recalculateFockMatrix == false) {
      std::cout << "   Fock matrix from a previous embedding/supermolecular calculation used!" << std::endl;
      fock = std::make_shared<FockMatrix<SCFMode>>(es->getFockMatrix());
    }
    else {
      fock = std::make_shared<FockMatrix<SCFMode>>(selectorAlgorithm->calcEmbeddedFockMatrix(_act[1], _env, settings.embedding));
    }
    _act[1]->template getElectronicStructure<SCFMode>()->setFockMatrix(*fock);
    _act[1]->template getElectronicStructure<SCFMode>()->toHDF5(_act[1]->getHDF5BaseName(), _act[1]->getSettings().identifier);
    // Localization
    bool localized = false;
    if (settings.localizedVirtualorbitals)
      localized = true;
    selectorAlgorithm->virtualOrbitalSpaceLocalization(coeffs, eigenvals, nOcc, nVirt, *fock, localized);
    // New orbitalenergies
    for_spin(eigenvals, coeffs) {
      auto temp = eigenvals_spin;
      eigenvals_spin.conservativeResize(coeffs_spin.cols());
      eigenvals_spin.segment(0, temp.size()) = temp;
      eigenvals_spin.segment(temp.size(), eigenvals_spin.size() - temp.size()).array() =
          std::numeric_limits<double>::infinity();
    };
    _act[1]->template getActiveOrbitalController<SCFMode>()->updateOrbitals(coeffs, eigenvals);
    printNewOrbitals(nOcc, eigenvals);
    writeOrbitalsToHDF5(coeffs, eigenvals, _act[1]);
  } /* localizedVirtualOrbitals */

  if (settings.mixingOccAndVirtualOrbitals) {
    printf("\n  ------- Mixing virtual and occupied orbitals of two subsystems -------\n");
    auto selectorAlgorithm = std::make_shared<VirtualOrbitalSelectionAlgorithms<SCFMode>>(_act[0], _env);
    selectorAlgorithm->occupiedVirtualMixing(settings.relaxation, settings.embedding);
    CoefficientMatrix<SCFMode> coefs = _act[0]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenValues =
        _act[0]->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
    auto nOcc = _act[0]->template getNOccupiedOrbitals<SCFMode>();
    printNewOrbitals(nOcc, eigenValues);
    writeOrbitalsToHDF5(coefs, eigenValues, _act[0]);
  } /*OccVirtMixing*/
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSpaceSelectionTask<SCFMode>::writeOrbitalsToHDF5(CoefficientMatrix<SCFMode>& coefficients,
                                                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues,
                                                                    std::shared_ptr<SystemController>& system) {
  std::string mode = (SCFMode == RESTRICTED) ? "res" : "unres";
  std::string fName = system->getSystemPath() + system->getSettings().name + settings.identifier + ".orbs." + mode + ".h5";
  HDF5::H5File file(fName.c_str(), H5F_ACC_TRUNC);
  std::string id = system->getSettings().identifier;
  unsigned int iCount = 0;
  for_spin(eigenvalues, coefficients) {
    std::string spin = (SCFMode == RESTRICTED) ? "" : "_alpha";
    if (iCount > 0)
      spin = "_beta";
    HDF5::save(file, "eigenvalues" + spin, eigenvalues_spin);
    HDF5::save(file, "coefficients" + spin, coefficients_spin);
    iCount += 1;
  };
  HDF5::save_scalar_attribute(file, "ID", id);
  printf("\n  Writes the new coefficients and orbital energies to file : %20s \n", fName.c_str());
  file.close();
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSpaceSelectionTask<SCFMode>::updateOrbitals(CoefficientMatrix<SCFMode>& coefs,
                                                               SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues,
                                                               SpinPolarizedData<SCFMode, std::vector<unsigned int>> indices) {
  for_spin(coefs, eigenvalues, indices) {
    // Update orbitals
    auto oldCoefficients = coefs_spin;
    auto oldOrbitalEnergies = eigenvalues_spin;
    // Resize Eigenvalues
    for (unsigned int iMO = 0; iMO < indices_spin.size(); ++iMO) {
      coefs_spin.col(iMO) = oldCoefficients.col(indices_spin[iMO]);
      eigenvalues_spin(iMO) = oldOrbitalEnergies(indices_spin[iMO]);
    }
    eigenvalues_spin.segment(indices_spin.size(), eigenvalues_spin.size() - indices_spin.size()).array() =
        std::numeric_limits<double>::infinity();
    for (unsigned int iMO = indices_spin.size(); iMO < oldCoefficients.cols(); iMO++) {
      coefs_spin.col(iMO).setZero();
    }
  };
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSpaceSelectionTask<SCFMode>::printNewOrbitals(SpinPolarizedData<SCFMode, unsigned int>& nOcc,
                                                                 SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues) {
  if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
    printf("Alpha:\n");
  }
  unsigned int spin_counter = 0;
  for_spin(nOcc, eigenvalues) {
    if (spin_counter == 1)
      printf("Beta:\n");
    int firstToPrint = nOcc_spin;
    unsigned int nToPrint = 0;
    switch (GLOBAL_PRINT_LEVEL) {
      case Options::GLOBAL_PRINT_LEVELS::MINIMUM:
        nToPrint = 10;
        break;
      case Options::GLOBAL_PRINT_LEVELS::NORMAL:
        nToPrint = 20;
        break;
      case Options::GLOBAL_PRINT_LEVELS::VERBOSE:
        nToPrint = 30;
        break;
      case Options::GLOBAL_PRINT_LEVELS::DEBUGGING:
        nToPrint = eigenvalues_spin.size() - nOcc_spin;
        break;
    }
    if (nToPrint > (eigenvalues_spin.size() - nOcc_spin))
      nToPrint = eigenvalues_spin.size() - nOcc_spin;
    if (firstToPrint < 0)
      firstToPrint = nOcc_spin;

    printf("%4s %6s %12s %19s\n", "", " # ", " Hartree ", "   eV   ");
    printf("%4s %6s %12s %19s\n", "", "---", "---------", "--------");
    for (unsigned int i = firstToPrint; i < nToPrint + firstToPrint; i++) {
      if (eigenvalues_spin[i] == std::numeric_limits<double>::infinity())
        continue;
      printf("%4s %5d   %+15.10f %+19.10f\n", "", (i + 1), eigenvalues_spin[i], eigenvalues_spin[i] * HARTREE_TO_EV);
    }
    spin_counter += 1;
  };
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSpaceSelectionTask<SCFMode>::checkInput() {
  if (_act.size() == 0 || _env.size() == 0) {
    throw SerenityError("You need to specify at least one active and one environment system!");
  }
  if (settings.localCanonicalVirtuals != 0.0 || settings.envCanonicalVirtuals != 0.0) {
    if (settings.localCanonicalVirtuals != 0.0 && settings.envCanonicalVirtuals != 0.0)
      throw SerenityError("You can not select local end environment virtual orbitals at the same time!");
    if (settings.localCanonicalVirtuals < 0.0)
      throw SerenityError("Negative localCanonicalVirtuals threshold not defined!");
    if (_env.size() > 1 && settings.envCanonicalVirtuals != 0.0)
      throw SerenityError(
          "You can select environment canonical orbitals only associated with one environment subsystem!");
  }
  if (settings.localizedVirtualorbitals == true || settings.localizedEnvVirtualorbitals == true) {
    if (settings.localizedVirtualorbitals == true && settings.localizedEnvVirtualorbitals == true)
      throw SerenityError("You can not select local end environment virtual orbitals at the same time!");
    if (_act.size() != 2)
      throw SerenityError("You need to specify two active subsystems, First: The original orbital space; Second: The "
                          "system to whom the new orbital space belongs");
  }
  if (settings.mixingOccAndVirtualOrbitals) {
    if (_env.size() != 2)
      throw SerenityError("Only two subsystems allowed for orbital mixing!");
    if (_act.size() != 1)
      throw SerenityError("Only one active system allowed for orbital mixing!");
    auto atoms_act = _act[0]->getGeometry()->getAtoms();
    auto atoms_env1 = _env[0]->getGeometry()->getAtoms();
    auto atoms_env2 = _env[1]->getGeometry()->getAtoms();
    bool sameGeo = true;
    if (atoms_act.size() == (atoms_env1.size() + atoms_env2.size())) {
      for (unsigned int i = 0; i < atoms_act.size(); i++) {
        if (i < atoms_env1.size()) {
          if (*atoms_act[i] == *atoms_env1[i]) {
            sameGeo = true;
          }
          else {
            sameGeo = false;
            break;
          }
        }
        else if (!atoms_act[i]->isDummy()) {
          sameGeo = false;
          break;
        }
      }
    }
    else {
      sameGeo = false;
    }
    if (sameGeo == false) {
      throw SerenityError("The active systems does not have the geometry of the environment system 1 and the ghost "
                          "atoms of the environment system 2!");
    }
  }
}

template class VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED>;
template class VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
