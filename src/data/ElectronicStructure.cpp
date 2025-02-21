/**
 * @file   ElectronicStructure.cpp
 *
 * @date   Aug 2, 2016
 * @author Jan Unsleber
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
#include "data/ElectronicStructure.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "integrals/OneIntControllerFactory.h"
#include "io/FormattedOutputStream.h"
#include "io/HDF5.h"
#include "parameters/Constants.h"
#include "potentials/bundles/PotentialBundle.h"

namespace Serenity {

/* ===============================
 *   Constructors and Destructor
 * =============================== */

template<Options::SCF_MODES SCFMode>
ElectronicStructure<SCFMode>::ElectronicStructure(std::shared_ptr<OneElectronIntegralController> oneEIntController,
                                                  const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals,
                                                  const SpinPolarizedData<SCFMode, unsigned int> nCoreElectrons)
  : state(ES_STATE::INITIAL),
    _diskmode(false),
    _oneEIntController(oneEIntController),
    _nOccupiedOrbitals(nOccupiedOrbitals),
    _molecularOrbitals(new OrbitalController<SCFMode>(oneEIntController->getBasisController(), nCoreElectrons)),
    _densityMatrixController(nullptr),
    _energyComponentController(new EnergyComponentController),
    _naddKinPotential(nullptr) {
  this->_densityMatrixController.reset(new DensityMatrixController<SCFMode>(_molecularOrbitals, nOccupiedOrbitals));
};

template<Options::SCF_MODES SCFMode>
ElectronicStructure<SCFMode>::ElectronicStructure(std::shared_ptr<BasisController> basisController,
                                                  std::shared_ptr<const Geometry> geometry,
                                                  std::shared_ptr<ExternalChargeController> externalCharges,
                                                  const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals,
                                                  const SpinPolarizedData<SCFMode, unsigned int> nCoreElectrons)
  : state(ES_STATE::INITIAL),
    _diskmode(false),
    _oneEIntController(OneIntControllerFactory::getInstance().produce(basisController, geometry, externalCharges)),
    _nOccupiedOrbitals(nOccupiedOrbitals),
    _molecularOrbitals(new OrbitalController<SCFMode>(basisController, nCoreElectrons)),
    _densityMatrixController(nullptr),
    _energyComponentController(new EnergyComponentController),
    _naddKinPotential(nullptr) {
  this->_densityMatrixController.reset(new DensityMatrixController<SCFMode>(_molecularOrbitals, nOccupiedOrbitals));
};

template<Options::SCF_MODES SCFMode>
ElectronicStructure<SCFMode>::ElectronicStructure(std::shared_ptr<OrbitalController<SCFMode>> molecularOrbitals,
                                                  std::shared_ptr<OneElectronIntegralController> oneEIntController,
                                                  const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals)
  : state(ES_STATE::GUESS),
    _diskmode(false),
    _oneEIntController(oneEIntController),
    _nOccupiedOrbitals(nOccupiedOrbitals),
    _molecularOrbitals(molecularOrbitals),
    _densityMatrixController(new DensityMatrixController<SCFMode>(molecularOrbitals, nOccupiedOrbitals)),
    _energyComponentController(new EnergyComponentController),
    _naddKinPotential(nullptr) {
  assert(oneEIntController->getBasisController() == molecularOrbitals->getBasisController());
};

template<Options::SCF_MODES SCFMode>
ElectronicStructure<SCFMode>::ElectronicStructure(std::string fBaseName, std::shared_ptr<BasisController> basis,
                                                  std::shared_ptr<const Geometry> geometry,
                                                  std::shared_ptr<ExternalChargeController> externalCharges, std::string id)
  : state(ES_STATE::GUESS),
    _diskmode(false),
    _oneEIntController(OneIntControllerFactory::getInstance().produce(basis, geometry, externalCharges)),
    _molecularOrbitals(new OrbitalController<SCFMode>(fBaseName, basis, id)),
    _densityMatrixController(new DensityMatrixController<SCFMode>(fBaseName, basis, id)),
    _energyComponentController(new EnergyComponentController),
    _naddKinPotential(nullptr) {
  _densityMatrixController->attachOrbitals(_molecularOrbitals, _densityMatrixController->getOccupations(), false);
  _molecularOrbitals->fromHDF5(fBaseName, id);
  _densityMatrixController->fromHDF5(fBaseName, id);
  try {
    this->fockFromHDF5(fBaseName, id);
  }
  catch (...) {
    _fockMatrix = nullptr;
  }
  if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
    fBaseName = fBaseName + ".energies.unres";
  }
  else {
    fBaseName = fBaseName + ".energies.res";
  }
  SpinPolarizedData<SCFMode, Eigen::VectorXd> occupations = _densityMatrixController->getOccupations();
  for_spin(occupations, _nOccupiedOrbitals) {
    _nOccupiedOrbitals_spin = (occupations_spin.array() > 0.0).count();
  };
  _energyComponentController->fromFile(fBaseName, id);
  auto& factory = OneIntControllerFactory::getInstance();
  _oneEIntController = factory.produce(_molecularOrbitals->getBasisController(), geometry, externalCharges);
};

/* ==============================
 *   Potentials and Fock Matrix
 * ============================== */

template<Options::SCF_MODES SCFMode>
void ElectronicStructure<SCFMode>::setFockMatrix(FockMatrix<SCFMode>& fock) {
  _fockMatrix = std::make_shared<FockMatrix<SCFMode>>(fock);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> ElectronicStructure<SCFMode>::getFockMatrix() {
  if (_fockMatrix) {
    return *_fockMatrix;
  }
  else {
    assert(this->potentialsAvailable() && "Tried to generate Fock matrix without any potentials present.");
    return this->getPotentialBundle()->getFockMatrix(this->getDensityMatrix(), _energyComponentController);
  }
}

/* ==============================
 *      I/O and Storage Mode
 * ============================== */

template<Options::SCF_MODES SCFMode>
void ElectronicStructure<SCFMode>::setDiskMode(bool diskmode, std::string fBaseName, std::string id) {
  _fBaseName = fBaseName;
  _id = id;
  if (diskmode) {
    if (_fBaseName.empty()) {
      throw SerenityError("Need to set file path before settings orbital controller to disk mode.");
    }
    if (_id.empty()) {
      throw SerenityError("Need to set file ID before settings orbital controller to disk mode.");
    }
  }
  _potentials = nullptr;
  _molecularOrbitals->setDiskMode(diskmode, fBaseName, id);
  _densityMatrixController->setDiskMode(diskmode, fBaseName, id);
  _diskmode = diskmode;
  if (diskmode) {
    _oneEIntController->clearOneInts();
  }
}

template<Options::SCF_MODES SCFMode>
void ElectronicStructure<SCFMode>::toHDF5(std::string fBaseName, std::string id) {
  _molecularOrbitals->toHDF5(fBaseName, id);
  _densityMatrixController->toHDF5(fBaseName, id);
  if (_fockMatrix) {
    std::string name = fBaseName;
    if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
      name += ".FockMatrix.unres.h5";
    }
    else {
      name += ".FockMatrix.res.h5";
    }
    HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
    const auto& fock = *_fockMatrix;
    unsigned int spinCounter = 0;
    for_spin(fock) {
      if (SCFMode == Options::SCF_MODES::RESTRICTED) {
        HDF5::save(file, "FockMatrix", fock_spin);
      }
      else if (SCFMode == Options::SCF_MODES::UNRESTRICTED && spinCounter == 0) {
        HDF5::save(file, "FockMatrix_alpha", fock_spin);
      }
      else if (SCFMode == Options::SCF_MODES::UNRESTRICTED && spinCounter == 1) {
        HDF5::save(file, "FockMatrix_beta", fock_spin);
      }
      spinCounter++;
    };
    HDF5::save_scalar_attribute(file, "ID", id);
    file.close();
  }
  if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
    fBaseName = fBaseName + ".energies.unres";
  }
  else {
    fBaseName = fBaseName + ".energies.res";
  }
  _energyComponentController->toFile(fBaseName, id);
}

template<>
void ElectronicStructure<Options::SCF_MODES::RESTRICTED>::fockFromHDF5(std::string fBaseName, std::string id) {
  HDF5::Filepath name(fBaseName + ".FockMatrix.res.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "FockMatrix");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  _fockMatrix = std::make_shared<FockMatrix<Options::SCF_MODES::RESTRICTED>>(
      FockMatrix<Options::SCF_MODES::RESTRICTED>(_oneEIntController->getBasisController()));
  HDF5::load(file, "FockMatrix", *_fockMatrix);
  file.close();
}
template<>
void ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>::fockFromHDF5(std::string fBaseName, std::string id) {
  HDF5::Filepath name(fBaseName + ".FockMatrix.unres.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "FockMatrix_alpha");
  HDF5::dataset_exists(file, "FockMatrix_beta");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  _fockMatrix = std::make_shared<FockMatrix<Options::SCF_MODES::UNRESTRICTED>>(
      FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_oneEIntController->getBasisController()));
  HDF5::load(file, "FockMatrix_alpha", _fockMatrix->alpha);
  HDF5::load(file, "FockMatrix_beta", _fockMatrix->beta);
  file.close();
}

template<>
void ElectronicStructure<Options::SCF_MODES::RESTRICTED>::printMOEnergies() const {
  const auto& eigenvalues = _molecularOrbitals->getEigenvalues();
  auto occ = _densityMatrixController->getOccupations();

  /* Get the orbital eigenvalue indices that should be printed.
   * By restricted the interval.--> First and number of eigenvalues required.
   */
  int firstToPrint = 0;
  unsigned int nToPrint = 0;
  unsigned int nOcc = (occ.array() > 0.0).count();
  switch (GLOBAL_PRINT_LEVEL) {
    case Options::GLOBAL_PRINT_LEVELS::MINIMUM:
      firstToPrint = nOcc - 1;
      nToPrint = 2;
      break;
    case Options::GLOBAL_PRINT_LEVELS::NORMAL:
      firstToPrint = nOcc - 10;
      nToPrint = 20;
      break;
    case Options::GLOBAL_PRINT_LEVELS::VERBOSE:
      firstToPrint = 0;
      nToPrint = nOcc + 10;
      break;
    case Options::GLOBAL_PRINT_LEVELS::DEBUGGING:
      firstToPrint = 0;
      nToPrint = eigenvalues.size();
      break;
  }
  if (nToPrint > eigenvalues.size())
    nToPrint = eigenvalues.size();
  if (firstToPrint < 0)
    firstToPrint = 0;

  printf("%4s %5s  %6s %12s %19s\n", "", " # ", " Occ. ", " Hartree ", "   eV   ");
  printf("%4s %5s  %6s %12s %19s\n", "", "---", "------", "---------", "--------");
  for (unsigned int i = firstToPrint; i < nToPrint + firstToPrint; i++) {
    if (eigenvalues[i] == std::numeric_limits<double>::infinity())
      continue;
    printf("%4s %5d   %4.2f %+15.10f %+19.10f\n", "", (i + 1), occ[i], eigenvalues[i], eigenvalues[i] * HARTREE_TO_EV);
  }
}

template<>
void ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>::printMOEnergies() const {
  const auto& eigenvalues = _molecularOrbitals->getEigenvalues();
  auto occ = _densityMatrixController->getOccupations();

  int firstToPrintAlpha = 0;
  unsigned int nToPrintAlpha = 0;
  unsigned int nOccAlpha = (occ.alpha.array() > 0.0).count();
  int firstToPrintBeta = 0;
  unsigned int nToPrintBeta = 0;
  unsigned int nOccBeta = (occ.beta.array() > 0.0).count();
  switch (GLOBAL_PRINT_LEVEL) {
    case Options::GLOBAL_PRINT_LEVELS::MINIMUM:
      firstToPrintAlpha = nOccAlpha - 1;
      nToPrintAlpha = 2;
      firstToPrintBeta = nOccBeta - 1;
      nToPrintBeta = 2;
      break;
    case Options::GLOBAL_PRINT_LEVELS::NORMAL:
      firstToPrintAlpha = nOccAlpha - 10;
      nToPrintAlpha = 20;
      firstToPrintBeta = nOccBeta - 10;
      nToPrintBeta = 20;
      break;
    case Options::GLOBAL_PRINT_LEVELS::VERBOSE:
      firstToPrintAlpha = 0;
      nToPrintAlpha = nOccAlpha + 10;
      firstToPrintBeta = 0;
      nToPrintBeta = nOccBeta + 10;
      break;
    case Options::GLOBAL_PRINT_LEVELS::DEBUGGING:
      firstToPrintAlpha = 0;
      nToPrintAlpha = eigenvalues.alpha.size();
      firstToPrintBeta = 0;
      nToPrintBeta = eigenvalues.beta.size();
      break;
  }
  if (nToPrintAlpha > eigenvalues.alpha.size())
    nToPrintAlpha = eigenvalues.alpha.size();
  if (firstToPrintAlpha < 0)
    firstToPrintAlpha = 0;
  if (nToPrintBeta > eigenvalues.beta.size())
    nToPrintBeta = eigenvalues.beta.size();
  if (firstToPrintBeta < 0)
    firstToPrintBeta = 0;

  printf("Alpha:\n");
  printf("%4s %5s  %6s %12s %19s\n", "", " # ", " Occ. ", " Hartree ", "   eV   ");
  printf("%4s %5s  %6s %12s %19s\n", "", "---", "------", "---------", "--------");
  for (unsigned int i = firstToPrintAlpha; i < nToPrintAlpha + firstToPrintAlpha; i++) {
    printf("%4s %5d   %4.2f %+15.10f %+19.10f\n", "", (i + 1), occ.alpha[i], eigenvalues.alpha[i],
           eigenvalues.alpha[i] * HARTREE_TO_EV);
  }
  printf("Beta:\n");
  printf("%4s %5s  %6s %12s %19s\n", "", " # ", " Occ. ", " Hartree ", "   eV   ");
  printf("%4s %5s  %6s %12s %19s\n", "", "---", "------", "---------", "--------");
  for (unsigned int i = firstToPrintBeta; i < nToPrintBeta + firstToPrintBeta; i++) {
    printf("%4s %5d   %4.2f %+15.10f %+19.10f\n", "", (i + 1), occ.beta[i], eigenvalues.beta[i],
           eigenvalues.beta[i] * HARTREE_TO_EV);
  }
}

template<>
ElectronicStructure<Options::SCF_MODES::RESTRICTED>::ElectronicStructure(std::shared_ptr<ElectronicStructure<RESTRICTED>> other)
  : _diskmode(other->getDiskMode()),
    _oneEIntController(other->getOneElectronIntegralController()),
    _nOccupiedOrbitals(other->getNOccupiedOrbitals()),
    _molecularOrbitals(other->getMolecularOrbitals()),
    _densityMatrixController(other->getDensityMatrixController()),
    _energyComponentController(other->getEnergyComponentController()) {
}

template<>
ElectronicStructure<Options::SCF_MODES::RESTRICTED>::ElectronicStructure(std::shared_ptr<ElectronicStructure<UNRESTRICTED>> other) {
  (void)other;
  throw SerenityError("ERROR: Forbidden conversion from UNRESTRICTED to RESTRICTED\n"
                      "       ElectronicStructure!");
}

template<>
ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>::ElectronicStructure(std::shared_ptr<ElectronicStructure<UNRESTRICTED>> other)
  : _diskmode(other->getDiskMode()),
    _oneEIntController(other->getOneElectronIntegralController()),
    _nOccupiedOrbitals(other->getNOccupiedOrbitals()),
    _molecularOrbitals(other->getMolecularOrbitals()),
    _densityMatrixController(other->getDensityMatrixController()),
    _energyComponentController(other->getEnergyComponentController()) {
}

template<>
ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>::ElectronicStructure(std::shared_ptr<ElectronicStructure<RESTRICTED>> other)
  : _diskmode(other->getDiskMode()),
    _oneEIntController(other->getOneElectronIntegralController()),
    _nOccupiedOrbitals(other->getNOccupiedOrbitals()),
    _energyComponentController(other->getEnergyComponentController()) {
  auto basisController = other->getMolecularOrbitals()->getBasisController();
  auto coefficients = std::make_unique<CoefficientMatrix<UNRESTRICTED>>(other->getMolecularOrbitals()->getCoefficients());
  auto eigenvalues =
      std::make_unique<SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>>(other->getMolecularOrbitals()->getEigenvalues());
  _molecularOrbitals = std::make_shared<OrbitalController<UNRESTRICTED>>(
      std::move(coefficients), basisController, *eigenvalues, other->getMolecularOrbitals()->getNCoreOrbitals());
  // Else, use the coefficients.
  _densityMatrixController = std::make_shared<DensityMatrixController<UNRESTRICTED>>(_molecularOrbitals, _nOccupiedOrbitals);
  if (other->checkFock())
    _fockMatrix = std::make_shared<FockMatrix<UNRESTRICTED>>(other->getFockMatrix());
}

template class ElectronicStructure<Options::SCF_MODES::RESTRICTED>;
template class ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
