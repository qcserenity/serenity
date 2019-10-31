/**
 * @file   SystemController.cpp
 * @author Thomas Dresselhaus, Jan Unsleber
 *
 * @date   20. Juli 2015, 16:51
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
#include "system/SystemController.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "basis/BasisController.h"
#include "basis/BasisFunctionProvider.h"
#include "data/matrices/CoefficientMatrix.h"
#include "potentials/CoulombPotential.h"
#include "potentials/bundles/DFTPotentials.h"
#include "potentials/EffectiveCorePotential.h"
#include "data/ElectronicStructure.h"
#include "potentials/FuncPotential.h"
#include "input/FunctionalClassResolver.h"
#include "grid/GridControllerFactory.h"
#include "potentials/HCorePotential.h"
#include "io/HDF5.h"
#include "io/Filesystem.h"
#include "potentials/HFPotential.h"
#include "potentials/bundles/HFPotentials.h"
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "scf/initialGuess/InitialGuessFactory.h"
#include "math/IntegerMaths.h"
#include "math/Matrix.h"
#include "data/OrbitalController.h"
#include "basis/Shell.h"
#include "geometry/XyzFileToGeometryConverter.h"


namespace Serenity {
using namespace std;

SystemController::SystemController(Settings settings) :
      _system(nullptr){
  //If system should be loaded: Load new settings in load/name/...
  bool read=!settings.load.empty();
  std::string loadPath;
  if (read){
    //Set load path
    if(settings.load.substr(settings.load.length()-1)!="/") settings.load=settings.load+"/";
    //Prevent overwriting of loaded systems
    if(settings.path==settings.load)
      throw SerenityError("ERROR: Write path and load path are the same. Stopping to prevent overwriting of loaded system");
    loadPath=settings.load+settings.name+"/";
    std::ifstream newSettingsFile;
    //open file
    newSettingsFile.open (loadPath+settings.name+".settings", std::ifstream::in);
    if(!newSettingsFile.good()){
      throw SerenityError("File "+loadPath+".settings not existent or corrupted");
    }
    Settings newSettings(newSettingsFile);
    //set the geometry to the one that should be loaded
    newSettings.geometry=loadPath+settings.name+".xyz";
    //set write path to value from old settings
    newSettings.path=settings.path;
    //set load path
    newSettings.load=loadPath;
    //set settings
    settings=newSettings;
    newSettingsFile.close();
  } else if (settings.geometry.empty()){
    throw SerenityError("ERROR: neither .xyz path ('geometry') nor 'load' path given for system: "+settings.name);
  }
  //Set basis path if not specified in input
  if(settings.basis.basisLibPath.empty()){
    if(const char* env_p = std::getenv("SERENITY_RESOURCES")){
      settings.basis.basisLibPath = (std::string)env_p+"basis/";
    }else{
      throw SerenityError("ERROR: Neither BasisLibPath nor environment variable SERENITY_RESOURCES set.");
    }
  }
  //Set path to write files to: path/name/ and create folder
  int start = settings.path.length()-settings.name.length()-1;
  if (start>=0){
   if (settings.path.substr(start).compare(settings.name+"/")){
     settings.path=settings.path+settings.name+"/";
   }
  } else {
    settings.path=settings.path+settings.name+"/";
  }
  makePath(settings.path);
  //print settings
  settings.printSettings();
  XyzFileToGeometryConverter reader(settings.geometry);
  std::shared_ptr<Geometry> geom;
  auto geo = reader.readGeometry();
  geom.reset(geo.release());
  _system.reset(new System(geom, settings));
  assert(_system);
  if(read) fromHDF5(loadPath+settings.name);
  setCharge(_system->_settings.charge);
  geom->printToFile(getHDF5BaseName(),settings.identifier);
  print();
}

SystemController::SystemController(std::shared_ptr<Geometry> geometry,Settings settings) :
      _system(nullptr){
  if(settings.basis.basisLibPath.empty()){
    if(const char* env_p = std::getenv("SERENITY_RESOURCES")){
      settings.basis.basisLibPath = (std::string)env_p + "basis/";
    }else{
      throw SerenityError("ERROR: Neither BasisLibPath nor environment variable SERENITY_RESOURCES set.");
    }
  }
  int start = settings.path.length()-settings.name.length()-1;
  if (start>=0){
    if (settings.path.substr(start).compare(settings.name+"/")){
      settings.path=settings.path+settings.name+"/";
    }
  } else {
    settings.path=settings.path+settings.name+"/";
  }
  makePath(settings.path);
  _system.reset(new System(geometry,settings));
  assert(_system);
  setCharge(_system->_settings.charge);
  //print settings
  settings.printSettings();
  geometry->printToFile(getHDF5BaseName(),settings.identifier);
  print();
}


/********************/
/* Getter functions */
/********************/
const Settings& SystemController::getSettings() const {
  return _system->_settings;
}

std::shared_ptr<BasisController> SystemController::getBasisController(
      Options::BASIS_PURPOSES basisPurpose) const {
  if (!_system->_basisControllers[basisPurpose]) produceBasisController(basisPurpose);
  return _system->_basisControllers[basisPurpose];
}


template<>void SystemController::produceScfTask<Options::SCF_MODES::RESTRICTED>() {
  _restrictedScfTask.reset(
          new ScfTask<Options::SCF_MODES::RESTRICTED>(shared_from_this()));
}
template<>void SystemController::produceScfTask<Options::SCF_MODES::UNRESTRICTED>() {
  _unrestrictedScfTask.reset(
          new ScfTask<Options::SCF_MODES::UNRESTRICTED>(shared_from_this()));
}
template<>std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED> >
      SystemController::getElectronicStructure<Options::SCF_MODES::RESTRICTED>() {
  if (!_system->_restrictedElectronicStructure){
    if (!_restrictedScfTask) produceScfTask<Options::SCF_MODES::RESTRICTED>();
    _restrictedScfTask->run();
  }
  return _system->_restrictedElectronicStructure;
}
template<>std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED> >
      SystemController::getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>() {
  if (!_system->_unrestrictedElectronicStructure){
    if (!_unrestrictedScfTask) produceScfTask<Options::SCF_MODES::UNRESTRICTED>();
    _unrestrictedScfTask->run();
  }
  return _system->_unrestrictedElectronicStructure;
}

template<> bool
SystemController::hasElectronicStructure<Options::SCF_MODES::RESTRICTED>() {
  return (nullptr!=_system->_restrictedElectronicStructure);
}
template<> bool
SystemController::hasElectronicStructure<Options::SCF_MODES::UNRESTRICTED>() {
  return (nullptr!=_system->_unrestrictedElectronicStructure);
}

template <Options::SCF_MODES T>
std::shared_ptr<OrbitalController<T> > SystemController::getActiveOrbitalController() {
  return this->getElectronicStructure<T>()->getMolecularOrbitals();
};
template std::shared_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED> >
SystemController::getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>();
template std::shared_ptr<OrbitalController<Options::SCF_MODES::UNRESTRICTED> >
SystemController::getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>();

template <>SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>
    SystemController::getNElectrons<Options::SCF_MODES::RESTRICTED>() const {
  return SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>(_system->_nElectrons);
}
template <>SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int>
    SystemController::getNElectrons<Options::SCF_MODES::UNRESTRICTED>() const {
  assert(isEven(_system->_nElectrons+_system->_settings.spin));
  assert(isEven(_system->_nElectrons-_system->_settings.spin));
  assert((int)_system->_nElectrons+_system->_settings.spin >= 0);
  assert((int)_system->_nElectrons-_system->_settings.spin >= 0);
    return makeUnrestrictedFromPieces<unsigned int>(
        (unsigned int)(_system->_nElectrons+_system->_settings.spin)/2,
        (unsigned int)(_system->_nElectrons-_system->_settings.spin)/2);
}

template <>SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>
    SystemController::getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>() const {
  return SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>(_system->_nElectrons /2);
}

template <>SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int>
    SystemController::getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>() const {
  assert(isEven(_system->_nElectrons+_system->_settings.spin));
  assert(isEven(_system->_nElectrons-_system->_settings.spin));
  assert((int)_system->_nElectrons+_system->_settings.spin >= 0);
  assert((int)_system->_nElectrons-_system->_settings.spin >= 0);
    return makeUnrestrictedFromPieces<unsigned int>(
        (unsigned int)(_system->_nElectrons+_system->_settings.spin)/2,
        (unsigned int)(_system->_nElectrons-_system->_settings.spin)/2);
}

template <>SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>
    SystemController::getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>()  {
  unsigned int nMolecularOrbitals = getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getNOrbitals();
  return SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>(nMolecularOrbitals -_system->_nElectrons /2);
}

template <>SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int>
    SystemController::getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>() {
  assert(isEven(_system->_nElectrons+_system->_settings.spin));
  assert(isEven(_system->_nElectrons-_system->_settings.spin));
  assert((int)_system->_nElectrons+_system->_settings.spin >= 0);
  assert((int)_system->_nElectrons-_system->_settings.spin >= 0);
  unsigned int nMolecularOrbitals = getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getNOrbitals();
    return makeUnrestrictedFromPieces<unsigned int>(
        (unsigned int)(nMolecularOrbitals -(_system->_nElectrons+_system->_settings.spin)/2),
        (unsigned int)(nMolecularOrbitals -(_system->_nElectrons-_system->_settings.spin)/2));
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
  //assert(!_system->_basisControllers[basisPurpose]);
    _system->_basisControllers[basisPurpose] = basisController;
    if (basisPurpose==Options::BASIS_PURPOSES::DEFAULT){
      _system->_restrictedElectronicStructure = nullptr;
      _system->_unrestrictedElectronicStructure = nullptr;
    }
}

template<>void SystemController::setElectronicStructure(
    std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED> > electronicStructure) {
  _system->_restrictedElectronicStructure = electronicStructure;
  _system->_lastSCFMode = Options::SCF_MODES::RESTRICTED;
}
template<>void SystemController::setElectronicStructure(
    std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED> > electronicStructure) {
  _system->_unrestrictedElectronicStructure = electronicStructure;
  _system->_lastSCFMode = Options::SCF_MODES::UNRESTRICTED;
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
  if(!isEven(_system->_nElectrons+_system->_settings.spin)){
    throw SerenityError("ERROR: Number of electron + spin yields an odd number."
        " Please check your charge/spin!");
  }
}

void SystemController::setSpin(const int spin) {
  _system->_settings.spin = spin;
}

/*************/
/* Disk Mode */
/*************/
void SystemController::setDiskMode(bool diskmode){
  if (diskmode) {
    // Delete RI integral data
    auto& fac = RI_J_IntegralControllerFactory::getInstance();
    auto cont = fac.produce(this->getAtomCenteredBasisController(),
        this->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));
    cont->notify();
//    // Delete grids
//    auto ptr1 = dynamic_pointer_cast<AtomCenteredGridController>(this->getGridController(Options::GRID_PURPOSES::DEFAULT));
//    auto ptr2 = dynamic_pointer_cast<AtomCenteredGridController>(this->getGridController(Options::GRID_PURPOSES::SMALL));
//    ptr1->notify();
//    ptr2->notify();
    // Delete 1e Ints
    this->getOneElectronIntegralController()->notify();
  }
  // Switch modes for electronic structures
  if (_system->_restrictedElectronicStructure){
    _system->_restrictedElectronicStructure->setDiskMode(diskmode,
        this->getHDF5BaseName(),
        this->getSettings().identifier);
  }
  if (_system->_unrestrictedElectronicStructure){
    _system->_unrestrictedElectronicStructure->setDiskMode(diskmode,
        this->getHDF5BaseName(),
        this->getSettings().identifier);
  }
}

/*****************************/
/* Default Potential getters */
/*****************************/

template<> std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> >
SystemController::getPotentials<Options::SCF_MODES::RESTRICTED,
                                Options::ELECTRONIC_STRUCTURE_THEORIES::HF>(
                                    Options::GRID_PURPOSES grid ){
  (void)grid;
  if (!_system->_restrictedElectronicStructure){
    this->setElectronicStructure(
        std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED> >(
            InitialGuessFactory::produce<Options::SCF_MODES::RESTRICTED>(
                _system->_settings.scf.initialguess)
                          ->calculateInitialGuess(this->getSharedPtr())
                          )
    );
  }
  std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED> > hcore
     (new HCorePotential<Options::SCF_MODES::RESTRICTED>(this->getSharedPtr()));
  std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED> > hf
     (new HFPotential<Options::SCF_MODES::RESTRICTED>(this->getSharedPtr(),
         _system->_restrictedElectronicStructure->getDensityMatrixController(),
         1.0,
         this->getSettings().basis.integralThreshold));
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > pot(
      new HFPotentials<Options::SCF_MODES::RESTRICTED>(hcore,hf,this->getGeometry()));
  return pot;
}
template<> std::shared_ptr<PotentialBundle<Options::SCF_MODES::UNRESTRICTED> >
SystemController::getPotentials<Options::SCF_MODES::UNRESTRICTED,
                                Options::ELECTRONIC_STRUCTURE_THEORIES::HF>(
                                    Options::GRID_PURPOSES grid ){
  (void)grid;
  if (!_system->_unrestrictedElectronicStructure){
    this->setElectronicStructure(
        std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED> >(
            InitialGuessFactory::produce<Options::SCF_MODES::UNRESTRICTED>(
                _system->_settings.scf.initialguess)
                          ->calculateInitialGuess(this->getSharedPtr())
                          )
    );
  }
  std::shared_ptr<Potential<Options::SCF_MODES::UNRESTRICTED> > hcore
     (new HCorePotential<Options::SCF_MODES::UNRESTRICTED>(this->getSharedPtr()));
  std::shared_ptr<Potential<Options::SCF_MODES::UNRESTRICTED> > hf
     (new HFPotential<Options::SCF_MODES::UNRESTRICTED>(this->getSharedPtr(),
         _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
         1.0,
         this->getSettings().basis.integralThreshold));
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::UNRESTRICTED> > pot(
      new HFPotentials<Options::SCF_MODES::UNRESTRICTED>(hcore,hf,this->getGeometry()));
  return pot;
}
template<> std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> >
SystemController::getPotentials<Options::SCF_MODES::RESTRICTED,
                                Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(
                                    Options::GRID_PURPOSES grid ){
  if (!_system->_restrictedElectronicStructure){
    this->setElectronicStructure(
        std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED> >(
            InitialGuessFactory::produce<Options::SCF_MODES::RESTRICTED>(
                _system->_settings.scf.initialguess)
                          ->calculateInitialGuess(this->getSharedPtr())
                          )
    );
  }
  // Hcore
  std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED> > hcore
     (new HCorePotential<Options::SCF_MODES::RESTRICTED>(this->getSharedPtr()));

  // XC Func
   auto functional = FunctionalClassResolver::resolveFunctional(this->getSettings().dft.functional);
   std::shared_ptr<FuncPotential<Options::SCF_MODES::RESTRICTED> > Vxc
      (new FuncPotential<Options::SCF_MODES::RESTRICTED>(
          this->getSharedPtr(),
          _system->_restrictedElectronicStructure->getDensityMatrixController(),
          this->getGridController(grid),
          functional));
   // J
   std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED> > J;
   if (!functional.isHybrid()){
     double thresh = this->getSettings().basis.integralThreshold;
     if (this->getSettings().dft.densityFitting == Options::DENS_FITS::RI){
       J = std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED> >(
           new CoulombPotential<Options::SCF_MODES::RESTRICTED>(this->getSharedPtr(),
               _system->_restrictedElectronicStructure->getDensityMatrixController(),
               RI_J_IntegralControllerFactory::getInstance().produce(
                   this->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                   this->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)),
                   thresh));
     } else {
       J = std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED> >(
           new HFPotential<Options::SCF_MODES::RESTRICTED>(this->getSharedPtr(),
               _system->_restrictedElectronicStructure->getDensityMatrixController(),
               0.0,
               thresh));
     }
   } else {
     double thresh = this->getSettings().basis.integralThreshold;
     J = std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED> >(
         new HFPotential<Options::SCF_MODES::RESTRICTED>(this->getSharedPtr(),
             _system->_restrictedElectronicStructure->getDensityMatrixController(),
             functional.getHfExchangeRatio(),
             thresh));
   }
  // Bundle
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > pot(
      new DFTPotentials<Options::SCF_MODES::RESTRICTED>(
          hcore,
          J,
          Vxc,
          this->getGeometry(),
          _system->_restrictedElectronicStructure->getDensityMatrixController(),
          this->getSettings().basis.integralThreshold));
  return pot;
}
template<> std::shared_ptr<PotentialBundle<Options::SCF_MODES::UNRESTRICTED> >
SystemController::getPotentials<Options::SCF_MODES::UNRESTRICTED,
                                Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(
                                    Options::GRID_PURPOSES grid ){
  if (!_system->_unrestrictedElectronicStructure){
    this->setElectronicStructure(
        std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED> >(
            InitialGuessFactory::produce<Options::SCF_MODES::UNRESTRICTED>(
                _system->_settings.scf.initialguess)
                          ->calculateInitialGuess(this->getSharedPtr())
                          )
    );
  }
  // Hcore
  std::shared_ptr<Potential<Options::SCF_MODES::UNRESTRICTED> > hcore
     (new HCorePotential<Options::SCF_MODES::UNRESTRICTED>(this->getSharedPtr()));
  // XC Func
  auto functional = FunctionalClassResolver::resolveFunctional(this->getSettings().dft.functional);
  std::shared_ptr<FuncPotential<Options::SCF_MODES::UNRESTRICTED> > Vxc
     (new FuncPotential<Options::SCF_MODES::UNRESTRICTED>(
         this->getSharedPtr(),
         _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
         this->getGridController(grid),
         functional));
  // J
  std::shared_ptr<Potential<Options::SCF_MODES::UNRESTRICTED> > J;
  if (!functional.isHybrid()){
    double thresh = this->getSettings().basis.integralThreshold;
    if (this->getSettings().dft.densityFitting == Options::DENS_FITS::RI){
      J = std::shared_ptr<Potential<Options::SCF_MODES::UNRESTRICTED> >(
          new CoulombPotential<Options::SCF_MODES::UNRESTRICTED>(this->getSharedPtr(),
              _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
              RI_J_IntegralControllerFactory::getInstance().produce(
                  this->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
                  this->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)),
                  thresh));
    } else {
      J = std::shared_ptr<Potential<Options::SCF_MODES::UNRESTRICTED> >(
          new HFPotential<Options::SCF_MODES::UNRESTRICTED>(this->getSharedPtr(),
              _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
              0.0,
              thresh));
    }
  } else {
    double thresh = this->getSettings().basis.integralThreshold;
    J = std::shared_ptr<Potential<Options::SCF_MODES::UNRESTRICTED> >(
        new HFPotential<Options::SCF_MODES::UNRESTRICTED>(this->getSharedPtr(),
            _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
            functional.getHfExchangeRatio(),
            thresh));
  }
  // Bundle
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::UNRESTRICTED> > pot(
      new DFTPotentials<Options::SCF_MODES::UNRESTRICTED>(
          hcore,
          J,
          Vxc,
          this->getGeometry(),
          _system->_unrestrictedElectronicStructure->getDensityMatrixController(),
          this->getSettings().basis.integralThreshold));
  return pot;
}



void SystemController::fromHDF5(std::string loadPath){
  try{
    this->getAtomCenteredBasisController()->fromHDF5(loadPath,this->getSettings().identifier);
    this->getAtomCenteredBasisController()->getBasis();
  }catch(...){
  }
  try{
  setElectronicStructure(
      std::make_shared<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>>(
          loadPath,
          this->getBasisController(),
          this->getGeometry(),
          _system->_settings.identifier));
  }catch(...){
  }
  try{
  setElectronicStructure(
      std::make_shared<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>(
          loadPath,
          this->getBasisController(),
          this->getGeometry(),
          _system->_settings.identifier));
  }catch(...){
  }
}

std::shared_ptr<SystemController> SystemController::operator+(
    SystemController& rhs) {
  Settings settings=this->getSettings();
  //New charge, spin
  settings.charge=this->getCharge()+rhs.getCharge();
  settings.spin=this->getSpin()+rhs.getSpin();
  settings.name=this->getSystemName()+"+"+rhs.getSystemName();
  settings.path = settings.path.substr(0, settings.path.size()-(this->getSystemName().size()+1));
  auto geom=std::make_shared<Geometry>();
  *geom+=*(this->getGeometry());
  *geom+=*(rhs.getGeometry());
  geom->deleteIdenticalAtoms();
  auto newSys=std::make_shared<SystemController>(geom,settings);
  auto basisController = newSys->getBasisController();
  auto lhsBasisFunc = this->getBasisController()->getNBasisFunctions();
  auto rhsBasisFunc = rhs.getBasisController()->getNBasisFunctions();
  unsigned int rhsBlockStart = lhsBasisFunc;
  if(rhsBlockStart==basisController->getNBasisFunctions()) rhsBlockStart=0;
  //Add electronic structures: Orbitals and energies.
  if (this->getSettings().scfMode == UNRESTRICTED
       or rhs.getSettings().scfMode == UNRESTRICTED) {
    auto nOccOrbsLeft=this->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
    auto nOccOrbsRight=rhs.getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
    CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> lhsCoeffs = this->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
            ->getMolecularOrbitals()->getCoefficients();
    CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> rhsCoeffs = rhs.getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
            ->getMolecularOrbitals()->getCoefficients();
    //OrbitalCoefficients
    std::unique_ptr<CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>> coeffs
    (new CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>(basisController));
    auto& newCoeffs=*coeffs;
    for_spin(newCoeffs,nOccOrbsLeft,nOccOrbsRight,lhsCoeffs,rhsCoeffs){
    newCoeffs_spin.block(0,0,lhsBasisFunc,nOccOrbsLeft_spin)
              =lhsCoeffs_spin.block(0,0,lhsBasisFunc,nOccOrbsLeft_spin);
    newCoeffs_spin.block(rhsBlockStart,nOccOrbsLeft_spin,rhsBasisFunc,
        nOccOrbsRight_spin)
              =rhsCoeffs_spin.block(0,0,rhsBasisFunc,nOccOrbsRight_spin);
    //Now the virtual orbitals
    if(rhsBlockStart!=0){
      newCoeffs_spin.block(0,nOccOrbsLeft_spin+nOccOrbsRight_spin,
          lhsBasisFunc,lhsBasisFunc-nOccOrbsLeft_spin)
                  =lhsCoeffs_spin.block(0,nOccOrbsLeft_spin,lhsBasisFunc,lhsBasisFunc-nOccOrbsLeft_spin);
      newCoeffs_spin.block(lhsBasisFunc,lhsBasisFunc+nOccOrbsRight_spin,
          rhsBasisFunc,rhsBasisFunc-nOccOrbsRight_spin)
                  =rhsCoeffs_spin.block(0,nOccOrbsRight_spin,rhsBasisFunc,rhsBasisFunc-nOccOrbsRight_spin);
    }
    };
    //Orbitalenergies
    std::unique_ptr<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>> newEsPtr
    (new SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>
    (basisController->getNBasisFunctions()));
    auto& newEs=*newEsPtr;
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd > lhsEs = this->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
            ->getMolecularOrbitals()->getEigenvalues();
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd > rhsEs = rhs.getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
                ->getMolecularOrbitals()->getEigenvalues();
    for_spin(newEs,nOccOrbsLeft,nOccOrbsRight,lhsEs,rhsEs){
    newEs_spin.segment(0,nOccOrbsLeft_spin)=lhsEs_spin.segment(0,nOccOrbsLeft_spin);
    newEs_spin.segment(nOccOrbsLeft_spin,nOccOrbsRight_spin)
                      =rhsEs_spin.segment(0,nOccOrbsRight_spin);
    if(rhsBlockStart!=0){
      newEs_spin.segment(nOccOrbsLeft_spin+nOccOrbsRight_spin,lhsBasisFunc-nOccOrbsLeft_spin)
                          =lhsEs_spin.segment(nOccOrbsLeft_spin,lhsBasisFunc-nOccOrbsLeft_spin);
      newEs_spin.segment(lhsBasisFunc+nOccOrbsRight_spin,rhsBasisFunc-nOccOrbsRight_spin)
                          =rhsEs_spin.segment(nOccOrbsRight_spin,rhsBasisFunc-nOccOrbsRight_spin);
    }
    };
    //New Orbitals
    auto newOrbs=std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>
        (std::move(coeffs),basisController,std::move(newEsPtr));
    //New electronic structure
    auto oneEIntController=OneIntControllerFactory::getInstance().produce(basisController,geom);
    auto newElecStruc=std::make_shared<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>>
      (newOrbs,oneEIntController,newSys->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>());
    //Set electronic structure
    newSys->setElectronicStructure<Options::SCF_MODES::UNRESTRICTED>(newElecStruc);
  }else{
    auto nOccOrbsLeft=this->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
    auto nOccOrbsRight=rhs.getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
    CoefficientMatrix<Options::SCF_MODES::RESTRICTED> lhsCoeffs = this->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
        ->getMolecularOrbitals()->getCoefficients();
    CoefficientMatrix<Options::SCF_MODES::RESTRICTED> rhsCoeffs = rhs.getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
            ->getMolecularOrbitals()->getCoefficients();
    std::unique_ptr<CoefficientMatrix<Options::SCF_MODES::RESTRICTED>> coeffs
    (new CoefficientMatrix<Options::SCF_MODES::RESTRICTED>(basisController));
    auto& newCoeffs=*coeffs;
    //Occupied Coeffs
    newCoeffs.block(0,0,lhsBasisFunc,nOccOrbsLeft)
        =lhsCoeffs.block(0,0,lhsBasisFunc,nOccOrbsLeft);
    newCoeffs.block(rhsBlockStart,nOccOrbsLeft,rhsBasisFunc,
        nOccOrbsRight)
        =rhsCoeffs.block(0,0,rhsBasisFunc,nOccOrbsRight);
    //Virtual Coeffs
    if(rhsBlockStart!=0){
      newCoeffs.block(0,nOccOrbsLeft+nOccOrbsRight,
          lhsBasisFunc,lhsBasisFunc-nOccOrbsLeft)
            =lhsCoeffs.block(0,nOccOrbsLeft,lhsBasisFunc,lhsBasisFunc-nOccOrbsLeft);
      newCoeffs.block(lhsBasisFunc,lhsBasisFunc+nOccOrbsRight,
          rhsBasisFunc,rhsBasisFunc-nOccOrbsRight)
            =rhsCoeffs.block(0,nOccOrbsRight,rhsBasisFunc,rhsBasisFunc-nOccOrbsRight);
    }
    std::unique_ptr<SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>> newEsPtr
        (new SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>
        (basisController->getNBasisFunctions()));
    auto& newEs=*newEsPtr;
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd > lhsEs = this->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
            ->getMolecularOrbitals()->getEigenvalues();
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd > rhsEs = rhs.getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
            ->getMolecularOrbitals()->getEigenvalues();
    //Occupied energies
    newEs.segment(0,nOccOrbsLeft)=lhsEs.segment(0,nOccOrbsLeft);
    newEs.segment(nOccOrbsLeft,nOccOrbsRight)
        =rhsEs.segment(0,nOccOrbsRight);
    //Virtual energies
    if(rhsBlockStart!=0){
      newEs.segment(nOccOrbsLeft+nOccOrbsRight,lhsBasisFunc-nOccOrbsLeft)
            =lhsEs.segment(nOccOrbsLeft,lhsBasisFunc-nOccOrbsLeft);
      newEs.segment(lhsBasisFunc+nOccOrbsRight,rhsBasisFunc-nOccOrbsRight)
            =rhsEs.segment(nOccOrbsRight,rhsBasisFunc-nOccOrbsRight);
    }
    //Build new Orbitals
    auto newOrbs=std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>
        (std::move(coeffs),basisController,std::move(newEsPtr));
    //New electronic structure
    auto oneEIntController=OneIntControllerFactory::getInstance().produce(basisController,geom);
    auto newElecStruc=std::make_shared<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>
        (newOrbs,oneEIntController,newSys->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>());
    //Set electronic structure
    newSys->setElectronicStructure<Options::SCF_MODES::RESTRICTED>(newElecStruc);
  }
  return newSys;
}

void SystemController::print(){
  std::string name = "System "+getSettings().name;
  printSubSectionTitle(name);
      auto m = getSettings().method;
      std::string method;
      Options::resolve<Options::ELECTRONIC_STRUCTURE_THEORIES>(method,m);
      printf("%4s Method:                %15s\n","",method.c_str());\
      if (getSettings().method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
        std::string functional;
        auto func = getSettings().dft.functional;
        Options::resolve<Options::XCFUNCTIONALS>(functional,func);
        printf("%4s Functional:            %15s\n","",functional.c_str());
        std::string fitting;
      }
      printf("%4s Basis Set:             %15s\n","",getSettings().basis.label.c_str());
      if(getGeometry()->hasAtomsWithECPs()) {
        printf("%4s ECP Start:             %15d\n","",getSettings().basis.firstECP);
      }
      getGeometry()->print();
}

/*********************/
/* Private functions */
/*********************/
void SystemController::produceBasisController(const Options::BASIS_PURPOSES basisPurpose) const {
  std::string label;
  if (basisPurpose==Options::BASIS_PURPOSES::DEFAULT){
    label = _system->_settings.basis.label;
  }else if(basisPurpose==Options::BASIS_PURPOSES::AUX_COULOMB){
    label = _system->_settings.basis.auxJLabel;
  }else if(basisPurpose==Options::BASIS_PURPOSES::SCF_DENS_GUESS){
    label = "DEF2-QZVP";
  }else if(basisPurpose==Options::BASIS_PURPOSES::MINBAS){
    label = "STO-3G";
  }else if(basisPurpose==Options::BASIS_PURPOSES::HUECKEL){
    label = "STO-6G";
  }else if(basisPurpose==Options::BASIS_PURPOSES::IAO_LOCALIZATION){
    label = "MINAO";
  }else if(basisPurpose==Options::BASIS_PURPOSES::AUX_CORREL){
    if (_system->_settings.basis.auxCLabel == ""){
      label = _system->_settings.basis.label;
      std::transform(label.begin(), label.end(),label.begin(), ::toupper);
      label += "-RI-C";
      auto auxCPath = _system->_settings.basis.basisLibPath + label;
      ifstream f(auxCPath.c_str());
      if (!f.good()) throw SerenityError("No default auxiliary (correlation optimized) basisset for chosen basis. Set AuxCLabel.");
    }else{
      label = _system->_settings.basis.auxCLabel;
    }
  }
  _system->_basisControllers[basisPurpose] = AtomCenteredBasisControllerFactory::produce(
          _system->_geometry,
          _system->_settings.basis.basisLibPath,
          _system->_settings.basis.makeSphericalBasis,
          (basisPurpose==Options::BASIS_PURPOSES::DEFAULT) ? true : false,
          _system->_settings.basis.firstECP,
          label);
}

void SystemController::produceGridController(const Options::GRID_PURPOSES gridPurpose) const {
  _system->_geometry->deleteIdenticalAtoms();
  _system->_gridControllers[gridPurpose] = GridControllerFactory::produce(_system->_geometry, _system->_settings, gridPurpose);
}

} /* namespace Serenity */
