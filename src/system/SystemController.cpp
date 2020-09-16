/**
 * @file   SystemController.cpp
 * @author Thomas Dresselhaus, Jan Unsleber
 *
 * @date   20. Juli 2015, 16:51
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
#include "system/SystemController.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/ElectrostaticPotentialOnGridController.h" //Potential on grid construction.
#include "data/matrices/CoefficientMatrix.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "geometry/Geometry.h"
#include "geometry/Geometry.h"                   //Geometry definition.
#include "geometry/MolecularSurfaceController.h" //Cavity generation.
#include "geometry/MolecularSurfaceFactory.h"
#include "geometry/XyzFileToGeometryConverter.h"
#include "grid/AtomCenteredGridController.h"
#include "grid/GridControllerFactory.h"
#include "integrals/OneIntControllerFactory.h" //One electron integral controller construction.
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "math/IntegerMaths.h"
#include "potentials/ERIPotential.h"
#include "potentials/FuncPotential.h"
#include "potentials/HCorePotential.h"
#include "potentials/PCMPotential.h"
#include "potentials/bundles/DFTPotentials.h"
#include "potentials/bundles/HFPotentials.h"
#include "scf/initialGuess/InitialGuessFactory.h"
#include "system/System.h"            //System definition.
#include "tasks/ScfTask.h"            //ScfTask construction.
#include "tasks/SystemAdditionTask.h" //needed in operator+

namespace Serenity {
using namespace std;

SystemController::~SystemController() = default;

SystemController::SystemController(Settings settings) : _system(nullptr) {
  // If system should be loaded: Load new settings in load/name/...
  bool read = !settings.load.empty();
  std::string loadPath;
  if (read) {
    // Set load path
    if (settings.load.substr(settings.load.length() - 1) != "/")
      settings.load = settings.load + "/";
    // Prevent overwriting of loaded systems
    if (settings.path == settings.load)
      throw SerenityError(
          "ERROR: Write path and load path are the same. Stopping to prevent overwriting of loaded system");
    loadPath = settings.load + settings.name + "/";
    std::ifstream newSettingsFile;
    // open file
    newSettingsFile.open(loadPath + settings.name + ".settings", std::ifstream::in);
    if (!newSettingsFile.good()) {
      throw SerenityError("File " + loadPath + ".settings not existent or corrupted");
    }
    Settings newSettings(newSettingsFile);
    // set the geometry to the one that should be loaded
    newSettings.geometry = loadPath + settings.name + ".xyz";
    // set write path to value from old settings
    newSettings.path = settings.path;
    // set load path
    newSettings.load = loadPath;
    // set settings
    settings = newSettings;
    newSettingsFile.close();
  }
  else if (settings.geometry.empty()) {
    throw SerenityError("ERROR: neither .xyz path ('geometry') nor 'load' path given for system: " + settings.name);
  }
  // Set basis path if not specified in input
  if (settings.basis.basisLibPath.empty()) {
    if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
      settings.basis.basisLibPath = (std::string)env_p + "basis/";
    }
    else {
      throw SerenityError("ERROR: Neither BasisLibPath nor environment variable SERENITY_RESOURCES set.");
    }
  }
  // Set path to write files to: path/name/ and create folder
  int start = settings.path.length() - settings.name.length() - 1;
  if (start >= 0) {
    if (settings.path.substr(start).compare(settings.name + "/")) {
      settings.path = settings.path + settings.name + "/";
    }
  }
  else {
    settings.path = settings.path + settings.name + "/";
  }
  if (not makePath(settings.path))
    throw SerenityError((std::string) "Failed to create directory: " + settings.path + "\n" +
                        "The directory does not already exist. A file with the same name may be\n" +
                        "preventing the directory creation.");
  // print settings
  settings.printSettings();
  XyzFileToGeometryConverter reader(settings.geometry);
  std::shared_ptr<Geometry> geom;
  auto geo = reader.readGeometry();
  geom.reset(geo.release());
  _system.reset(new System(geom, settings));
  assert(_system);
  if (read)
    fromHDF5(loadPath + settings.name);
  setCharge(_system->_settings.charge);
  geom->printToFile(getHDF5BaseName(), settings.identifier);
  print();
}

SystemController::SystemController(std::shared_ptr<Geometry> geometry, Settings settings) : _system(nullptr) {
  if (settings.basis.basisLibPath.empty()) {
    if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
      settings.basis.basisLibPath = (std::string)env_p + "basis/";
    }
    else {
      throw SerenityError("ERROR: Neither BasisLibPath nor environment variable SERENITY_RESOURCES set.");
    }
  }
  int start = settings.path.length() - settings.name.length() - 1;
  if (start >= 0) {
    if (settings.path.substr(start).compare(settings.name + "/")) {
      settings.path = settings.path + settings.name + "/";
    }
  }
  else {
    settings.path = settings.path + settings.name + "/";
  }
  if (not makePath(settings.path))
    throw SerenityError((std::string) "Failed to create directory: " + settings.path + "\n" +
                        "The directory does not already exist. A file with the same name may be\n" +
                        "preventing the directory creation.");
  _system.reset(new System(geometry, settings));
  assert(_system);
  setCharge(_system->_settings.charge);
  // print settings
  settings.printSettings();
  geometry->printToFile(getHDF5BaseName(), settings.identifier);
  print();
}

/********************/
/* Getter functions */
/********************/
const Settings& SystemController::getSettings() const {
  return _system->_settings;
}

std::shared_ptr<BasisController> SystemController::getBasisController(Options::BASIS_PURPOSES basisPurpose) const {
  if (!_system->_basisControllers[basisPurpose])
    produceBasisController(basisPurpose);
  return _system->_basisControllers[basisPurpose];
}

template<>
void SystemController::produceScfTask<RESTRICTED>() {
  _restrictedScfTask.reset(new ScfTask<RESTRICTED>(shared_from_this()));
}
template<>
void SystemController::produceScfTask<UNRESTRICTED>() {
  _unrestrictedScfTask.reset(new ScfTask<UNRESTRICTED>(shared_from_this()));
}
template<>
std::shared_ptr<ElectronicStructure<RESTRICTED>> SystemController::getElectronicStructure<RESTRICTED>() {
  if (!_system->_restrictedElectronicStructure) {
    if (!_restrictedScfTask)
      produceScfTask<RESTRICTED>();
    _restrictedScfTask->run();
  }
  return _system->_restrictedElectronicStructure;
}
template<>
std::shared_ptr<ElectronicStructure<UNRESTRICTED>> SystemController::getElectronicStructure<UNRESTRICTED>() {
  if (!_system->_unrestrictedElectronicStructure) {
    if (!_unrestrictedScfTask)
      produceScfTask<UNRESTRICTED>();
    _unrestrictedScfTask->run();
  }
  return _system->_unrestrictedElectronicStructure;
}

template<>
bool SystemController::hasElectronicStructure<RESTRICTED>() {
  return (nullptr != _system->_restrictedElectronicStructure);
}
template<>
bool SystemController::hasElectronicStructure<UNRESTRICTED>() {
  return (nullptr != _system->_unrestrictedElectronicStructure);
}

template<>
void SystemController::setElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(
    std::shared_ptr<ElectrostaticPotentialOnGridController<RESTRICTED>> potential, MOLECULAR_SURFACE_TYPES surfaceType) {
  _system->_restrictedElectrostaticPotentialsOnGridControllers[surfaceType] = potential;
}

template<>
void SystemController::setElectrostaticPotentialOnMolecularSurfaceController<UNRESTRICTED>(
    std::shared_ptr<ElectrostaticPotentialOnGridController<UNRESTRICTED>> potential, MOLECULAR_SURFACE_TYPES surfaceType) {
  _system->_unrestrictedElectrostaticPotentialsOnGridControllers[surfaceType] = potential;
}

template<>
void SystemController::produceElectrostaticPotentialOnMolecularSurfaceController<Options::SCF_MODES::RESTRICTED>(
    MOLECULAR_SURFACE_TYPES surfaceType) {
  auto gridController = this->getMolecularSurface(surfaceType)->getGridController();
  auto densityMatrixController = this->getElectronicStructure<RESTRICTED>()->getDensityMatrixController();
  this->setElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(
      std::make_shared<ElectrostaticPotentialOnGridController<Options::SCF_MODES::RESTRICTED>>(
          gridController, densityMatrixController, this->getGeometry(), this->getHDF5BaseName(), _system->_settings.pcm.cacheSize),
      surfaceType);
}

template<>
void SystemController::produceElectrostaticPotentialOnMolecularSurfaceController<UNRESTRICTED>(MOLECULAR_SURFACE_TYPES surfaceType) {
  auto gridController = this->getMolecularSurface(surfaceType)->getGridController();
  auto densityMatrixController = this->getElectronicStructure<UNRESTRICTED>()->getDensityMatrixController();
  this->setElectrostaticPotentialOnMolecularSurfaceController<UNRESTRICTED>(
      std::make_shared<ElectrostaticPotentialOnGridController<Options::SCF_MODES::UNRESTRICTED>>(
          gridController, densityMatrixController, this->getGeometry(), this->getHDF5BaseName(), _system->_settings.pcm.cacheSize),
      surfaceType);
}

template<>
std::shared_ptr<ElectrostaticPotentialOnGridController<RESTRICTED>>
SystemController::getElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(MOLECULAR_SURFACE_TYPES surfaceType) {
  if (!_system->_restrictedElectrostaticPotentialsOnGridControllers[surfaceType])
    produceElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(surfaceType);
  return _system->_restrictedElectrostaticPotentialsOnGridControllers[surfaceType];
}

template<>
std::shared_ptr<ElectrostaticPotentialOnGridController<UNRESTRICTED>>
SystemController::getElectrostaticPotentialOnMolecularSurfaceController<UNRESTRICTED>(MOLECULAR_SURFACE_TYPES surfaceType) {
  if (!_system->_unrestrictedElectrostaticPotentialsOnGridControllers[surfaceType])
    produceElectrostaticPotentialOnMolecularSurfaceController<UNRESTRICTED>(surfaceType);
  return _system->_unrestrictedElectrostaticPotentialsOnGridControllers[surfaceType];
}

template<Options::SCF_MODES T>
std::shared_ptr<OrbitalController<T>> SystemController::getActiveOrbitalController() {
  return this->getElectronicStructure<T>()->getMolecularOrbitals();
};
template std::shared_ptr<OrbitalController<RESTRICTED>> SystemController::getActiveOrbitalController<RESTRICTED>();
template std::shared_ptr<OrbitalController<UNRESTRICTED>> SystemController::getActiveOrbitalController<UNRESTRICTED>();

template<>
SpinPolarizedData<RESTRICTED, unsigned int> SystemController::getNElectrons<RESTRICTED>() const {
  return SpinPolarizedData<RESTRICTED, unsigned int>(_system->_nElectrons);
}
template<>
SpinPolarizedData<UNRESTRICTED, unsigned int> SystemController::getNElectrons<UNRESTRICTED>() const {
  assert(isEven(_system->_nElectrons + _system->_settings.spin));
  assert(isEven(_system->_nElectrons - _system->_settings.spin));
  assert((int)_system->_nElectrons + _system->_settings.spin >= 0);
  assert((int)_system->_nElectrons - _system->_settings.spin >= 0);
  return makeUnrestrictedFromPieces<unsigned int>((unsigned int)(_system->_nElectrons + _system->_settings.spin) / 2,
                                                  (unsigned int)(_system->_nElectrons - _system->_settings.spin) / 2);
}

template<>
SpinPolarizedData<RESTRICTED, unsigned int> SystemController::getNOccupiedOrbitals<RESTRICTED>() const {
  return SpinPolarizedData<RESTRICTED, unsigned int>(_system->_nElectrons / 2);
}

template<>
SpinPolarizedData<UNRESTRICTED, unsigned int> SystemController::getNOccupiedOrbitals<UNRESTRICTED>() const {
  assert(isEven(_system->_nElectrons + _system->_settings.spin));
  assert(isEven(_system->_nElectrons - _system->_settings.spin));
  assert((int)_system->_nElectrons + _system->_settings.spin >= 0);
  assert((int)_system->_nElectrons - _system->_settings.spin >= 0);
  return makeUnrestrictedFromPieces<unsigned int>((unsigned int)(_system->_nElectrons + _system->_settings.spin) / 2,
                                                  (unsigned int)(_system->_nElectrons - _system->_settings.spin) / 2);
}

template<>
SpinPolarizedData<RESTRICTED, unsigned int> SystemController::getNVirtualOrbitals<RESTRICTED>() {
  unsigned int nMolecularOrbitals = getActiveOrbitalController<RESTRICTED>()->getNOrbitals();
  return SpinPolarizedData<RESTRICTED, unsigned int>(nMolecularOrbitals - _system->_nElectrons / 2);
}

template<>
SpinPolarizedData<UNRESTRICTED, unsigned int> SystemController::getNVirtualOrbitals<UNRESTRICTED>() {
  assert(isEven(_system->_nElectrons + _system->_settings.spin));
  assert(isEven(_system->_nElectrons - _system->_settings.spin));
  assert((int)_system->_nElectrons + _system->_settings.spin >= 0);
  assert((int)_system->_nElectrons - _system->_settings.spin >= 0);
  unsigned int nMolecularOrbitals = getActiveOrbitalController<UNRESTRICTED>()->getNOrbitals();
  return makeUnrestrictedFromPieces<unsigned int>(
      (unsigned int)(nMolecularOrbitals - (_system->_nElectrons + _system->_settings.spin) / 2),
      (unsigned int)(nMolecularOrbitals - (_system->_nElectrons - _system->_settings.spin) / 2));
}

/**
 * @brief Forworded getter for _settings.pcm.use.
 */
bool SystemController::getSystemContinuumModelMode() {
  return _system->_settings.pcm.use;
}
/********************/
/* Setter functions */
/********************/

void SystemController::setBasisController(std::shared_ptr<AtomCenteredBasisController> basisController,
                                          Options::BASIS_PURPOSES basisPurpose) {
  /*
   * It should not be possible to override an existing basis controller
   * assign it before it is ever needed, or don't assign it!
   */
  // assert(!_system->_basisControllers[basisPurpose]);
  _system->_basisControllers[basisPurpose] = basisController;
  if (basisPurpose == Options::BASIS_PURPOSES::DEFAULT) {
    _system->_restrictedElectronicStructure = nullptr;
    _system->_unrestrictedElectronicStructure = nullptr;
  }
}

template<>
void SystemController::setElectronicStructure(std::shared_ptr<ElectronicStructure<RESTRICTED>> electronicStructure) {
  _system->_restrictedElectronicStructure = electronicStructure;
  _system->_lastSCFMode = RESTRICTED;
}
template<>
void SystemController::setElectronicStructure(std::shared_ptr<ElectronicStructure<UNRESTRICTED>> electronicStructure) {
  _system->_unrestrictedElectronicStructure = electronicStructure;
  _system->_lastSCFMode = UNRESTRICTED;
}

void SystemController::setCharge(const int charge) {
  _system->_settings.charge = charge;
  /*
   * Make sure the default basis is built as the effective charge of the nuclei and
   * thus also the number of electrons are dependent on whether effective core potentials
   * are used. However, whether effective core potentials are used is determined by the
   * basis set.
   */
  this->getBasisController();
  /*
   * Adjust nElectrons
   */
  _system->_nElectrons = 0;
  for (const auto& atom : _system->_geometry->getAtoms()) {
    _system->_nElectrons += atom->getEffectiveCharge();
  }
  assert((int)_system->_nElectrons >= _system->_settings.charge);
  _system->_nElectrons -= _system->_settings.charge;
  if (!isEven(_system->_nElectrons + _system->_settings.spin)) {
    throw SerenityError("ERROR: Number of electron + spin yields an odd number."
                        " Please check your charge/spin!");
  }
  _system->_settings.printSettings();
}

void SystemController::setSpin(const int spin) {
  _system->_settings.spin = spin;
  _system->_settings.printSettings();
}

/*************/
/* Disk Mode */
/*************/
void SystemController::setDiskMode(bool diskmode) {
  if (diskmode) {
    // Delete RI integral data
    auto& fac = RI_J_IntegralControllerFactory::getInstance();
    auto cont = fac.produce(this->getAtomCenteredBasisController(),
                            this->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));
    cont->notify();
    //    // Delete grids
    //    auto ptr1 =
    //    dynamic_pointer_cast<AtomCenteredGridController>(this->getGridController(Options::GRID_PURPOSES::DEFAULT));
    //    auto ptr2 =
    //    dynamic_pointer_cast<AtomCenteredGridController>(this->getGridController(Options::GRID_PURPOSES::SMALL));
    //    ptr1->notify();
    //    ptr2->notify();
    // Delete 1e Ints
    this->getOneElectronIntegralController()->notify();
  }
  // Switch modes for electronic structures
  if (_system->_restrictedElectronicStructure) {
    _system->_restrictedElectronicStructure->setDiskMode(diskmode, this->getHDF5BaseName(), this->getSettings().identifier);
  }
  if (_system->_unrestrictedElectronicStructure) {
    _system->_unrestrictedElectronicStructure->setDiskMode(diskmode, this->getHDF5BaseName(), this->getSettings().identifier);
  }
}

void SystemController::setSystemContinuumModelMode(bool newMode) {
  this->_system->_settings.pcm.use = newMode;
}

/*****************************/
/* Default Potential getters */
/*****************************/

template<>
std::shared_ptr<PotentialBundle<RESTRICTED>>
SystemController::getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>(Options::GRID_PURPOSES grid) {
  (void)grid;
  if (!_system->_restrictedElectronicStructure) {
    this->setElectronicStructure(std::shared_ptr<ElectronicStructure<RESTRICTED>>(
        InitialGuessFactory::produce<RESTRICTED>(_system->_settings.scf.initialguess)->calculateInitialGuess(this->getSharedPtr())));
  }
  auto hcore = std::make_shared<HCorePotential<RESTRICTED>>(this->getSharedPtr());
  auto hf = std::make_shared<ERIPotential<RESTRICTED>>(
      this->getSharedPtr(), _system->_restrictedElectronicStructure->getDensityMatrixController(), 1.0,
      this->getSettings().basis.integralThreshold, this->getSettings().basis.integralIncrementThresholdStart,
      this->getSettings().basis.integralIncrementThresholdEnd, this->getSettings().basis.incrementalSteps);
  bool usesPCM = _system->_settings.pcm.use;
  auto pcm = std::make_shared<PCMPotential<RESTRICTED>>(
      _system->_settings.pcm, this->getBasisController(), this->getGeometry(),
      (usesPCM) ? this->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE) : nullptr,
      (usesPCM) ? this->getElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(MOLECULAR_SURFACE_TYPES::ACTIVE)
                : nullptr);
  return std::make_shared<HFPotentials<RESTRICTED>>(hcore, hf, pcm, this->getGeometry());
}
template<>
std::shared_ptr<PotentialBundle<UNRESTRICTED>>
SystemController::getPotentials<UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>(Options::GRID_PURPOSES grid) {
  (void)grid;
  if (!_system->_unrestrictedElectronicStructure) {
    this->setElectronicStructure(std::shared_ptr<ElectronicStructure<UNRESTRICTED>>(
        InitialGuessFactory::produce<UNRESTRICTED>(_system->_settings.scf.initialguess)->calculateInitialGuess(this->getSharedPtr())));
  }
  auto hcore = std::make_shared<HCorePotential<UNRESTRICTED>>(this->getSharedPtr());
  auto hf = std::make_shared<ERIPotential<UNRESTRICTED>>(
      this->getSharedPtr(), _system->_unrestrictedElectronicStructure->getDensityMatrixController(), 1.0,
      this->getSettings().basis.integralThreshold, this->getSettings().basis.integralIncrementThresholdStart,
      this->getSettings().basis.integralIncrementThresholdEnd, this->getSettings().basis.incrementalSteps);
  bool usesPCM = _system->_settings.pcm.use;
  auto pcm = std::make_shared<PCMPotential<UNRESTRICTED>>(
      _system->_settings.pcm, this->getBasisController(), this->getGeometry(),
      (usesPCM) ? this->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE) : nullptr,
      (usesPCM) ? this->getElectrostaticPotentialOnMolecularSurfaceController<UNRESTRICTED>(MOLECULAR_SURFACE_TYPES::ACTIVE)
                : nullptr);
  return std::make_shared<HFPotentials<UNRESTRICTED>>(hcore, hf, pcm, this->getGeometry());
}
template<>
std::shared_ptr<PotentialBundle<RESTRICTED>>
SystemController::getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(Options::GRID_PURPOSES grid) {
  if (!_system->_restrictedElectronicStructure) {
    this->setElectronicStructure(std::shared_ptr<ElectronicStructure<RESTRICTED>>(
        InitialGuessFactory::produce<RESTRICTED>(_system->_settings.scf.initialguess)->calculateInitialGuess(this->getSharedPtr())));
  }
  // Hcore
  auto hcore = std::make_shared<HCorePotential<RESTRICTED>>(this->getSharedPtr());

  // XC Func
  auto functional = resolveFunctional(this->getSettings().dft.functional);
  auto Vxc = std::make_shared<FuncPotential<RESTRICTED>>(
      this->getSharedPtr(), _system->_restrictedElectronicStructure->getDensityMatrixController(),
      this->getGridController(grid), functional);
  // J
  std::shared_ptr<Potential<RESTRICTED>> J;
  double thresh = this->getSettings().basis.integralThreshold;
  J = std::make_shared<ERIPotential<RESTRICTED>>(
      this->getSharedPtr(), _system->_restrictedElectronicStructure->getDensityMatrixController(),
      functional.getHfExchangeRatio(), thresh, this->getSettings().basis.integralIncrementThresholdStart,
      this->getSettings().basis.integralIncrementThresholdEnd, this->getSettings().basis.incrementalSteps, true,
      functional.getLRExchangeRatio(), functional.getRangeSeparationParameter());
  bool usesPCM = _system->_settings.pcm.use;
  auto pcm = std::make_shared<PCMPotential<RESTRICTED>>(
      _system->_settings.pcm, this->getBasisController(), this->getGeometry(),
      (usesPCM) ? this->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE) : nullptr,
      (usesPCM) ? this->getElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(MOLECULAR_SURFACE_TYPES::ACTIVE)
                : nullptr);
  // Bundle
  return std::make_shared<DFTPotentials<RESTRICTED>>(hcore, J, Vxc, pcm, this->getGeometry(),
                                                     _system->_restrictedElectronicStructure->getDensityMatrixController(),
                                                     this->getSettings().basis.integralThreshold);
}
template<>
std::shared_ptr<PotentialBundle<UNRESTRICTED>>
SystemController::getPotentials<UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(Options::GRID_PURPOSES grid) {
  if (!_system->_unrestrictedElectronicStructure) {
    this->setElectronicStructure(std::shared_ptr<ElectronicStructure<UNRESTRICTED>>(
        InitialGuessFactory::produce<UNRESTRICTED>(_system->_settings.scf.initialguess)->calculateInitialGuess(this->getSharedPtr())));
  }
  // Hcore
  auto hcore = std::make_shared<HCorePotential<UNRESTRICTED>>(this->getSharedPtr());
  // XC Func
  auto functional = resolveFunctional(this->getSettings().dft.functional);
  auto Vxc = std::make_shared<FuncPotential<UNRESTRICTED>>(
      this->getSharedPtr(), _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
      this->getGridController(grid), functional);
  // J
  std::shared_ptr<Potential<UNRESTRICTED>> J;
  double thresh = this->getSettings().basis.integralThreshold;
  J = std::make_shared<ERIPotential<UNRESTRICTED>>(
      this->getSharedPtr(), _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
      functional.getHfExchangeRatio(), thresh, this->getSettings().basis.integralIncrementThresholdStart,
      this->getSettings().basis.integralIncrementThresholdEnd, this->getSettings().basis.incrementalSteps, true,
      functional.getLRExchangeRatio(), functional.getRangeSeparationParameter());
  bool usesPCM = _system->_settings.pcm.use;
  auto pcm = std::make_shared<PCMPotential<UNRESTRICTED>>(
      _system->_settings.pcm, this->getBasisController(), this->getGeometry(),
      (usesPCM) ? this->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE) : nullptr,
      (usesPCM) ? this->getElectrostaticPotentialOnMolecularSurfaceController<UNRESTRICTED>(MOLECULAR_SURFACE_TYPES::ACTIVE)
                : nullptr);
  // Bundle
  return std::make_shared<DFTPotentials<UNRESTRICTED>>(
      hcore, J, Vxc, pcm, this->getGeometry(), _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
      this->getSettings().basis.integralThreshold);
}

void SystemController::fromHDF5(std::string loadPath) {
  try {
    this->getAtomCenteredBasisController()->fromHDF5(loadPath, this->getSettings().identifier);
    this->getAtomCenteredBasisController()->getBasis();
  }
  catch (...) {
  }
  try {
    setElectronicStructure(std::make_shared<ElectronicStructure<UNRESTRICTED>>(
        loadPath, this->getBasisController(), this->getGeometry(), _system->_settings.identifier));
  }
  catch (...) {
  }
  try {
    setElectronicStructure(std::make_shared<ElectronicStructure<RESTRICTED>>(
        loadPath, this->getBasisController(), this->getGeometry(), _system->_settings.identifier));
  }
  catch (...) {
  }
}

std::shared_ptr<SystemController> SystemController::operator+(SystemController& rhs) {
  // Set the geometry of the joint system manually in order
  // to have the constructor of SystemController show the final
  // geometry when printing it to the output.
  Settings settings = this->getSettings();
  settings.charge = this->getCharge() + rhs.getCharge();
  settings.spin = this->getSpin() + rhs.getSpin();
  settings.name = this->getSystemName() + "+" + rhs.getSystemName();
  settings.path = settings.path.substr(0, settings.path.size() - (this->getSystemName().size() + 1));
  auto geom = std::make_shared<Geometry>();
  *geom += *(this->getGeometry());
  *geom += *(rhs.getGeometry());
  geom->deleteIdenticalAtoms();
  auto newSys = std::make_shared<SystemController>(geom, settings);
  if (this->getSettings().scfMode == UNRESTRICTED or rhs.getSettings().scfMode == UNRESTRICTED) {
    SystemAdditionTask<UNRESTRICTED> addTask(newSys, {this->getSharedPtr(), rhs.getSharedPtr()});
    addTask.settings.checkSuperGeom = true;
    addTask.settings.checkSuperCharge = true;
    addTask.settings.addOccupiedOrbitals = true;
    addTask.run();
  }
  else {
    SystemAdditionTask<RESTRICTED> addTask(newSys, {this->getSharedPtr(), rhs.getSharedPtr()});
    addTask.settings.checkSuperGeom = true;
    addTask.settings.checkSuperCharge = true;
    addTask.settings.addOccupiedOrbitals = true;
    addTask.run();
  }
  return newSys;
}

void SystemController::print() {
  std::string name = "System " + getSettings().name;
  printSubSectionTitle(name);
  auto m = getSettings().method;
  std::string method;
  Options::resolve<Options::ELECTRONIC_STRUCTURE_THEORIES>(method, m);
  printf("%4s Method:                %15s\n", "", method.c_str());
  if (getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    std::string functional;
    auto func = getSettings().dft.functional;
    Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(functional, func);
    printf("%4s Functional:            %15s\n", "", functional.c_str());
    std::string fitting;
  }
  printf("%4s Basis Set:             %15s\n", "", getSettings().basis.label.c_str());
  if (getGeometry()->hasAtomsWithECPs()) {
    printf("%4s ECP Start:             %15d\n", "", getSettings().basis.firstECP);
  }
  getGeometry()->print();
}

void SystemController::setMolecularSurface(std::shared_ptr<MolecularSurfaceController> surface,
                                           MOLECULAR_SURFACE_TYPES surfaceType) {
  _system->_molecularSurfaces[surfaceType] = surface;
}

std::shared_ptr<MolecularSurfaceController> SystemController::getMolecularSurface(MOLECULAR_SURFACE_TYPES surfaceType) {
  if (!_system->_molecularSurfaces[surfaceType]) {
    if (surfaceType == MOLECULAR_SURFACE_TYPES::ACTIVE) {
      this->produceMolecularSurface();
    }
    else {
      throw SerenityError("Logic error in molecular surface creation. Molecular surfaces containing more than one "
                          "molecule can not be created from the SystemController.");
    }
  }
  return _system->_molecularSurfaces[surfaceType];
}

/*********************/
/* Private functions */
/*********************/
void SystemController::produceMolecularSurface() {
  _system->_molecularSurfaces[MOLECULAR_SURFACE_TYPES::ACTIVE] =
      MolecularSurfaceFactory::produce(this->getGeometry(), _system->_settings.pcm);
}

void SystemController::produceBasisController(const Options::BASIS_PURPOSES basisPurpose) const {
  std::string label;
  if (basisPurpose == Options::BASIS_PURPOSES::DEFAULT) {
    label = _system->_settings.basis.label;
  }
  else if (basisPurpose == Options::BASIS_PURPOSES::AUX_COULOMB) {
    label = _system->_settings.basis.auxJLabel;
  }
  else if (basisPurpose == Options::BASIS_PURPOSES::SCF_DENS_GUESS) {
    label = "DEF2-QZVP";
  }
  else if (basisPurpose == Options::BASIS_PURPOSES::MINBAS) {
    label = "STO-3G";
  }
  else if (basisPurpose == Options::BASIS_PURPOSES::HUECKEL) {
    label = "STO-6G";
  }
  else if (basisPurpose == Options::BASIS_PURPOSES::IAO_LOCALIZATION) {
    label = "MINAO";
  }
  else if (basisPurpose == Options::BASIS_PURPOSES::AUX_CORREL) {
    if (_system->_settings.basis.auxCLabel == "") {
      label = _system->_settings.basis.label;
      std::transform(label.begin(), label.end(), label.begin(), ::toupper);
      label += "-RI-C";
      auto auxCPath = _system->_settings.basis.basisLibPath + label;
      ifstream f(auxCPath.c_str());
      if (!f.good())
        throw SerenityError(
            (std::string) "No default auxiliary (correlation optimized) basis set for chosen basis. Set AuxCLabel: " + label);
    }
    else {
      label = _system->_settings.basis.auxCLabel;
    }
  }
  _system->_basisControllers[basisPurpose] = AtomCenteredBasisControllerFactory::produce(
      _system->_geometry, _system->_settings.basis.basisLibPath, _system->_settings.basis.makeSphericalBasis,
      (basisPurpose == Options::BASIS_PURPOSES::DEFAULT) ? true : false, _system->_settings.basis.firstECP, label);
}

void SystemController::produceGridController(const Options::GRID_PURPOSES gridPurpose) const {
  _system->_geometry->deleteIdenticalAtoms();
  _system->_gridControllers[gridPurpose] = GridControllerFactory::produce(_system->_geometry, _system->_settings, gridPurpose);
}

void SystemController::setSCFMode(Options::SCF_MODES mode) {
  _system->_settings.scfMode = mode;
}

std::string SystemController::getSystemName() {
  return _system->_settings.name;
}

std::string SystemController::getSystemIdentifier() {
  return _system->_settings.identifier;
}

std::string SystemController::getSystemPath() {
  return _system->_settings.path;
}

std::string SystemController::getHDF5BaseName() {
  return _system->_settings.path + _system->_settings.name;
}

int SystemController::getCharge() const {
  return _system->_settings.charge;
}

enum Options::SCF_MODES SystemController::getLastSCFMode() const {
  return _system->_lastSCFMode;
};

int SystemController::getSpin() const {
  return _system->_settings.spin;
}

SpinPolarizedData<UNRESTRICTED, unsigned int> SystemController::getNAlphaAndBetaElectrons() const {
  assert(_system->_settings.spin % 2 == (int)_system->_nElectrons % 2);
  assert(_system->_nElectrons + _system->_settings.spin > 0 && _system->_nElectrons - _system->_settings.spin > 0);
  return makeUnrestrictedFromPieces<unsigned int>((unsigned int)(_system->_nElectrons + _system->_settings.spin) / 2,
                                                  (unsigned int)(_system->_nElectrons - _system->_settings.spin) / 2);
}

bool SystemController::isOpenShell() const {
  return (_system->_settings.spin != 0);
}

std::shared_ptr<AtomCenteredBasisController>
SystemController::getAtomCenteredBasisController(Options::BASIS_PURPOSES basisPurpose) const {
  if (!_system->_basisControllers[basisPurpose]) {
    produceBasisController(basisPurpose);
  }
  return _system->_basisControllers[basisPurpose];
}

std::shared_ptr<OneElectronIntegralController>
SystemController::getOneElectronIntegralController(Options::BASIS_PURPOSES basisPurpose) const {
  auto& factory = OneIntControllerFactory::getInstance();
  return factory.produce(this->getBasisController(basisPurpose), this->getGeometry());
}

std::shared_ptr<Geometry> SystemController::getGeometry() const {
  return _system->_geometry;
}

std::shared_ptr<GridController> SystemController::getGridController(Options::GRID_PURPOSES gridPurpose) const {
  if (!_system->_gridControllers[gridPurpose])
    produceGridController(gridPurpose);
  return _system->_gridControllers[gridPurpose];
}

void SystemController::setGridController(std::shared_ptr<GridController> gridController,
                                         Options::GRID_PURPOSES gridPurpose) const {
  _system->_gridControllers[gridPurpose] = gridController;
}

std::shared_ptr<AtomCenteredGridController>
SystemController::getAtomCenteredGridController(Options::GRID_PURPOSES gridPurpose) const {
  if (!_system->_gridControllers[gridPurpose])
    produceGridController(gridPurpose);
  // TODO hack, an ugly dynamic pointer cast. Should be cleaned away.
  std::shared_ptr<AtomCenteredGridController> result =
      std::dynamic_pointer_cast<AtomCenteredGridController>(_system->_gridControllers[gridPurpose]);
  assert(result);
  return result;
}

const std::vector<std::shared_ptr<Atom>>& SystemController::getAtoms() const {
  return _system->_geometry->getAtoms();
}

unsigned int SystemController::getNAtoms() const {
  return _system->_geometry->getAtoms().size();
}

} /* namespace Serenity */
