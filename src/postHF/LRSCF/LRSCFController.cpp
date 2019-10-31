/**
 * @file LRSCFController.cpp
 *
 * @date Dec 06, 2018
 * @author Michael Boeckers, Niklas Niemeyer, Johannes Toelle
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
#include "postHF/LRSCF/LRSCFController.h"

/* Include Serenity Internal Headers */
#include "tasks/LRSCFTask.h"
#include "postHF/LRSCF/Tools/Besley.h"
#include "data/OrbitalController.h"
#include "io/HDF5.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "basis/AtomCenteredBasisController.h"

namespace Serenity {

  template<Options::SCF_MODES SCFMode>
  LRSCFController<SCFMode>::LRSCFController(
      std::shared_ptr<SystemController> system,
      LRSCFTaskSettings& settings) :
      _system(system),
      _settings(settings),
      _coefficients(_system->getActiveOrbitalController<SCFMode>()->getCoefficients()),
      _orbitalEnergies(_system->getActiveOrbitalController<SCFMode>()->getEigenvalues()),
      _excitationVectors(nullptr){

    //Setup initial occupation vector corresponding to orbitals stored in system
    auto nOcc = _system->getNOccupiedOrbitals<SCFMode>();
    auto nVirt = _system->getNVirtualOrbitals<SCFMode>();
    for_spin(_occupation,nOcc,nVirt) {
      _occupation_spin = Eigen::VectorXi::Zero(nOcc_spin+nVirt_spin);
      _occupation_spin.segment(0,nOcc_spin) =  Eigen::VectorXi::Ones(nOcc_spin);
      _occupation_spin *= (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2 : 1;
    };

    //Creates Reference List of orbital indices
    for_spin(_indexWhiteList,nOcc,nVirt){
      _indexWhiteList_spin.resize(nOcc_spin+nVirt_spin);
      std::iota(std::begin(_indexWhiteList_spin), std::end(_indexWhiteList_spin), 0);
    };

    // Restricts orbital space according to the user's input and adjusts occupation vector
    this->editReference();
  }

  template<Options::SCF_MODES SCFMode>
  std::shared_ptr<std::vector<Eigen::MatrixXd> > LRSCFController<SCFMode>::getExcitationVectors(
      Options::LRSCF_TYPE type) {
    if (!_excitationVectors || _type != type) loadFromH5(type);
    return _excitationVectors;
  }

  template<Options::SCF_MODES SCFMode>
  std::shared_ptr<Eigen::VectorXd > LRSCFController<SCFMode>::getExcitationEnergies(
      Options::LRSCF_TYPE type) {
    if (!_excitationEnergies || _type != type) loadFromH5(type);
    return _excitationEnergies;
  }

  template<Options::SCF_MODES SCFMode>
  void LRSCFController<SCFMode>::loadFromH5(Options::LRSCF_TYPE type) {
    try {
      //Try to initialize from h5
      std::string mode = (SCFMode == RESTRICTED) ? "res":"unres";
      std::string type_str = (type==Options::LRSCF_TYPE::ISOLATED) ? ".iso" : (type==Options::LRSCF_TYPE::UNCOUPLED) ? ".fdeu" : ".fdec";
      std::string fName = _system->getSettings().path+_system->getSettings().name+type_str+"_lrscf." + mode + ".h5";
      HDF5::Filepath name(fName);
      HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      HDF5::dataset_exists(file,"X");
      HDF5::dataset_exists(file,"Y");
      std::vector<Eigen::MatrixXd> XY(2);
      HDF5::load(file,"X",XY[0]);
      HDF5::load(file,"Y",XY[1]);
      HDF5::dataset_exists(file,"EIGENVALUES");
      Eigen::VectorXd e;
      HDF5::load(file,"EIGENVALUES",e);
      file.close();
      //Found solution vectors of type, reset _excitationEnergies and _excitationVectors
      if(e.rows() != XY[0].cols()) throw SerenityError("The number of eigenvectors (X) and eigenvalues do not match!");
      if(e.rows() != XY[1].cols()) throw SerenityError("The number of eigenvectors (Y) and eigenvalues do not match!");
      _excitationVectors.reset();
      _excitationEnergies.reset();
      _excitationEnergies = std::make_shared<Eigen::VectorXd>(e);
      _excitationVectors = std::make_shared<std::vector<Eigen::MatrixXd> >(XY);
      _type = type;
      printf("\nFound %4i excitation vectors in %20s with eigenvalues:\n",(unsigned int)e.rows(),
          fName.c_str());
    } catch (...) {
      throw SerenityError("Cannot find the excitation vectors and excitation energies from HDF5.");
    }
  }

  template<Options::SCF_MODES SCFMode>
  void LRSCFController<SCFMode>::setSolution(
      std::shared_ptr<std::vector<Eigen::MatrixXd> > eigenvectors,
      std::shared_ptr<Eigen::VectorXd> eigenvalues,
      Options::LRSCF_TYPE type) {
  _excitationVectors = eigenvectors;
  _excitationEnergies = eigenvalues;
  _type = type;
  //write to h5
  std::string mode = (SCFMode==RESTRICTED) ? "res":"unres";
  std::string type_str = (type==Options::LRSCF_TYPE::ISOLATED) ? ".iso" : (type==Options::LRSCF_TYPE::UNCOUPLED) ? ".fdeu" : ".fdec";
  std::string fName = _system->getSettings().path+_system->getSettings().name+type_str+"_lrscf." + mode + ".h5";
  HDF5::H5File file(fName.c_str(), H5F_ACC_TRUNC);
  HDF5::save_scalar_attribute(file,"ID",_system->getSettings().identifier);
  HDF5::save(file,"X",(*eigenvectors)[0]);
  if((*eigenvectors).size() == 2 && !_settings.tda) {
    HDF5::save(file,"Y",(*eigenvectors)[1]);
  } else if ((*eigenvectors).size() == 1 && _settings.tda) {
    Eigen::MatrixXd zero = Eigen::MatrixXd::Zero((*eigenvectors)[0].rows(),(*eigenvectors)[0].cols());
    HDF5::save(file,"Y",zero);
  } else {
    assert(false);
  }
  HDF5::save(file,"EIGENVALUES",(*eigenvalues));
  file.close();
}


  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode,unsigned int> LRSCFController<SCFMode>::getNOccupied(){
    SpinPolarizedData<SCFMode,unsigned int> nOcc;
    for_spin(nOcc,_occupation) {
      nOcc_spin = (SCFMode == Options::SCF_MODES::RESTRICTED) ? _occupation_spin.sum() / 2 : _occupation_spin.sum();
    };
    return nOcc;
  }

  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode,unsigned int> LRSCFController<SCFMode>::getNVirtual(){
    SpinPolarizedData<SCFMode,unsigned int> nVirt;
    auto nOcc = this->getNOccupied();
    for_spin(nVirt,nOcc,_occupation) {
      nVirt_spin = _occupation_spin.size() - nOcc_spin;
    };
    return nVirt;
  }

  template<Options::SCF_MODES SCFMode>
  CoefficientMatrix<SCFMode>& LRSCFController<SCFMode>::getCoefficients() {
    return _coefficients;
  }

  template<Options::SCF_MODES SCFMode>
  std::shared_ptr<BasisController> LRSCFController<SCFMode>::getBasisController() {
    return _coefficients.getBasisController();
  }

  template<Options::SCF_MODES SCFMode>
    std::shared_ptr<MatrixInBasis<SCFMode> > LRSCFController<SCFMode>::getFockNonCanon() {
      return _fockNonCanon;
  }

  template<Options::SCF_MODES SCFMode>
    void LRSCFController<SCFMode>::setFockNonCanon(std::shared_ptr<MatrixInBasis<SCFMode> > fockNonCanon) {
      _fockNonCanon = fockNonCanon;
  }

  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode,Eigen::VectorXd> LRSCFController<SCFMode>::getEigenvalues() {
    return _orbitalEnergies;
  }

  template<Options::SCF_MODES SCFMode>
  std::shared_ptr<GridController> LRSCFController<SCFMode>::getGridController() {
    return _system->getGridController();
  }

  template<Options::SCF_MODES SCFMode>
  const Settings& LRSCFController<SCFMode>::getSysSettings() {
    return _system->getSettings();
  }

  template<Options::SCF_MODES SCFMode>
  std::shared_ptr<SystemController> LRSCFController<SCFMode>::getSys() {
    return _system;
  }

  template<Options::SCF_MODES SCFMode>
  const LRSCFTaskSettings& LRSCFController<SCFMode>::getLRSCFSettings(){
    return _settings;
  }

  template<>
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED,std::vector<unsigned int > >
  LRSCFController<Options::SCF_MODES::RESTRICTED>::getSet() {
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED,std::vector<unsigned int> > set;
    set = _settings.setAlpha;
    return set;
  }

  template<>
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,std::vector<unsigned int > >
  LRSCFController<Options::SCF_MODES::UNRESTRICTED>::getSet() {
    auto set = makeUnrestrictedFromPieces<std::vector<unsigned int> >(_settings.setAlpha, _settings.setBeta);
    return set;
  }

  template<>
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED,std::vector<unsigned int > >
  LRSCFController<Options::SCF_MODES::RESTRICTED>::getExclude() {
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED,std::vector<unsigned int> > exclude;
    exclude = _settings.excludeAlpha;
    return exclude;
  }

  template<>
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,std::vector<unsigned int > >
  LRSCFController<Options::SCF_MODES::UNRESTRICTED>::getExclude() {
    auto exclude = makeUnrestrictedFromPieces<std::vector<unsigned int> >(_settings.excludeAlpha, _settings.excludeBeta);
    return exclude;
  }

  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode,std::vector<unsigned int> > LRSCFController<SCFMode>::getReferenceOrbitals() {
    return _indexWhiteList;
  }

  template<Options::SCF_MODES SCFMode>
  void LRSCFController<SCFMode>::editReference(){
    auto set = this->getSet();
    auto exclude = this->getExclude();
    bool editReference = false;
    for_spin(set,exclude) {
      if (!set_spin.empty()) editReference = true;
      if (!exclude_spin.empty()) editReference = true;
    };
    if (_settings.besleyAtoms > 0) editReference = true;
    if (!_settings.energyInclusion.empty()) editReference = true;
    if (!_settings.energyExclusion.empty()) editReference = true;
    if (_settings.localVirtualOrbitals != 0.0) editReference = true;
    if (_settings.envVirtualOrbitals != 0.0) editReference = true;

    //Orbitals will not be edited. Work with reference set from system.
    if (!editReference) return;

    //Some variables needed
    auto nOcc = _system->getNOccupiedOrbitals<SCFMode>();
 
    //Restrict Orbitals according to Besley's criterion for occupied and virtual orbitals
    if(_settings.besleyAtoms > 0) {
      if(_settings.besleyCutoff.size() != 2) throw SerenityError("Keyword besleyCutoff needs two arguments!");
      Besley<SCFMode> besley(_system, _settings.besleyAtoms, _settings.besleyCutoff);
      _indexWhiteList = besley.getWhiteList();
    }

    //Virtual Orbital Partitioning (located on the subsysten/ environment subsystem)
    //or virtual orbitals located on the environment subsystem
    if(_settings.localVirtualOrbitals != 0 || _settings.envVirtualOrbitals != 0){
     //Get atom indices for non-dummy atoms
      auto atoms = _system ->getAtoms();
     
      //Calculate Overlap between Virtual MOs in basis of atoms of active subsystem and entire basis
      auto nBasisFunc = _system->getBasisController()->getNBasisFunctions();
      auto basisIndices = _system->getAtomCenteredBasisController()->getBasisIndices(); 
      auto coeff = this->getCoefficients(); 
      const auto& oneIntController = _system->getOneElectronIntegralController();
      const auto& overlaps = oneIntController->getOverlapIntegrals();
      auto nOccupied = this->getNOccupied();
      auto nVirtual = this->getNVirtual();
      //find the atom indices, which are not dummy
      unsigned int counter = 0;
      std::vector<unsigned int > index;
      for (auto atom : atoms){
        if (atom -> isDummy() == false){
          index.push_back(counter);
        }
        counter += 1;
      }
   
      //Entire Dimensions
      for_spin(nOccupied,nVirtual,coeff,_indexWhiteList){
        _indexWhiteList_spin.clear();
        for(unsigned int iOcc = 0;  iOcc < nOccupied_spin; iOcc++){
          _indexWhiteList_spin.push_back(iOcc);
        }
        //ToDo: Assumes that the atoms of the subsystem are behind each other
        unsigned int start = basisIndices[index[0]].first;
        unsigned int end = basisIndices[index[index.size()-1]].second;
         
        Eigen::MatrixXd overlapMOModified = coeff_spin.block(start,nOccupied_spin,end - start,nVirtual_spin).transpose() *
                    overlaps.block(start, 0, end - start, nBasisFunc) *
                    coeff_spin.block(0,nOccupied_spin,nBasisFunc,nVirtual_spin);

        for (unsigned int col = 0; col < overlapMOModified.cols(); col ++){
          if(_settings.localVirtualOrbitals != 0)
          if (overlapMOModified(col,col) > _settings.localVirtualOrbitals){
            _indexWhiteList_spin.push_back(col+nOccupied_spin);
          }
          if(_settings.envVirtualOrbitals != 0)
          if ((1-overlapMOModified(col,col)) >= _settings.envVirtualOrbitals){
            _indexWhiteList_spin.push_back(col+nOccupied_spin);
          }
        }
      };
    }/* Virtual Orbital Partitioning */

    unsigned int iSpin = 0;
    for_spin(_coefficients,_orbitalEnergies,_occupation,_indexWhiteList,set,exclude,nOcc){
      //Set cut off
      if (!set_spin.empty()) {
        _indexWhiteList_spin.clear();
        if(!exclude_spin.empty()) throw SerenityError("Keywords 'set' and 'exclude' cannot be used at the same time!");
        for(auto index : set_spin) {
          _indexWhiteList_spin.push_back(index);
        }
      }

      // Exclude cut off
      if (!exclude_spin.empty()) {
        if(!set_spin.empty()) throw SerenityError("Keywords 'set' and 'exclude' cannot be used at the same time!");
        for (int index = _indexWhiteList_spin.size() - 1; index >= 0; --index) {
          if (std::find(exclude_spin.begin(), exclude_spin.end(), index) != exclude_spin.end()) {
            _indexWhiteList_spin.erase(_indexWhiteList_spin.begin() + index);
          }
        }
      }

      // Energy cut off
      if(!_settings.energyExclusion.empty()) {
        if(_settings.energyExclusion.size() % 2 != 0) throw SerenityError("Keyword energyExclusion needs an even number of arguments!");
        if(!(_settings.energyInclusion.empty())) throw SerenityError("Keywords 'energyExclusion' and 'energyInclusion' cannot be used at the same time!");
        for(int index = _indexWhiteList_spin.size() - 1; index >= 0; --index){
          for(unsigned int pair = 0; pair < _settings.energyExclusion.size(); pair += 2) {
            if(_settings.energyExclusion[pair + 0] >= _settings.energyExclusion[pair + 1]){
              throw SerenityError("The first energyExclusion parameter must be smaller than the second!");
            }
            if(_orbitalEnergies_spin(_indexWhiteList_spin[index])* HARTREE_TO_EV > _settings.energyExclusion[pair + 0] &&
               _orbitalEnergies_spin(_indexWhiteList_spin[index])* HARTREE_TO_EV < _settings.energyExclusion[pair + 1]){
              _indexWhiteList_spin.erase(_indexWhiteList_spin.begin() + index);
            }
          }
        }
      }

      if(!_settings.energyInclusion.empty()) {
        if(_settings.energyInclusion.size() % 2 != 0) throw SerenityError("Keyword energyInclusion needs an even number of arguments!");
        if(!(_settings.energyExclusion.empty())) throw SerenityError("Keywords 'energyExclusion' and 'energyInclusion' cannot be used at the same time!");
        for(int index = _indexWhiteList_spin.size() - 1; index >= 0; --index){
          bool removeIndex = true;
          for(unsigned int pair = 0; pair < _settings.energyInclusion.size(); pair += 2) {
            if(_settings.energyInclusion[pair + 0] >= _settings.energyInclusion[pair + 1]){
              throw SerenityError("The first energyExclusion parameter must be smaller than the second!");
            }
            if (_orbitalEnergies_spin(_indexWhiteList_spin[index])* HARTREE_TO_EV > _settings.energyInclusion[pair + 0] &&
                _orbitalEnergies_spin(_indexWhiteList_spin[index])* HARTREE_TO_EV < _settings.energyInclusion[pair + 1]) {
              removeIndex = false;
              break;
            }
          }
          // If the energy of the orbital with a given index is not found in any of the intervals, it gets removed from the list
          if (removeIndex) _indexWhiteList_spin.erase(_indexWhiteList_spin.begin() + index);
        }
      }/* Energy cut off */

      // Update orbitals
      auto oldCoefficients = _coefficients_spin;
      auto oldOrbitalEnergies = _orbitalEnergies_spin;
      _coefficients_spin.setZero();

      //Adapt occupation vector and orbital energies
      _occupation_spin.resize(_indexWhiteList_spin.size());
      _orbitalEnergies_spin.resize(_indexWhiteList_spin.size());
      _occupation_spin.setZero();

      for (unsigned int iMO = 0; iMO < _indexWhiteList_spin.size(); ++iMO) {
        _coefficients_spin.col(iMO) = oldCoefficients.col(_indexWhiteList_spin[iMO]);
        _orbitalEnergies_spin(iMO) = oldOrbitalEnergies(_indexWhiteList_spin[iMO]);
        _occupation_spin(iMO) = (_indexWhiteList_spin[iMO] < nOcc_spin) ?
            ((SCFMode == Options::SCF_MODES::RESTRICTED) ? 2 : 1) : 0;
      }/* Update orbitals */

      //Print info
      if (SCFMode == Options::SCF_MODES::RESTRICTED) {
        printf("\n Reference orbitals : \n");
      } else {
        printf("\n %s Reference orbitals : \n",(iSpin == 0) ? "Alpha" : "Beta");
      }
      for (unsigned int iMO = 0; iMO < _indexWhiteList_spin.size(); ++iMO) {
        printf("%4i",_indexWhiteList_spin[iMO]+1);
        if ( (iMO+1) % 10 == 0) printf("\n");
      }
      printf("\n");
      iSpin += 1;
    };

    //Using isolated as default since it will not matter in this case
    this->indexToH5(Options::LRSCF_TYPE::ISOLATED);
  }

  //Edit reference for exclude projection where the individual calculation is carried out in LRSCFTask
  template<Options::SCF_MODES SCFMode>
  void LRSCFController<SCFMode>::editReference(SpinPolarizedData<SCFMode,std::vector<unsigned int> > indexWhiteList, 
                                                unsigned int system, Options::LRSCF_TYPE type){
    //Set new indices
    _indexWhiteList = indexWhiteList;
    //Some information
    auto nOcc = _system->getNOccupiedOrbitals<SCFMode>();
    //The new index list is referring to not earlier changed dimensions in the coefficient matrix
    //Therefore is a reinitialization of the coefficients / orbital energies neccessary
    _coefficients = _system->getActiveOrbitalController<SCFMode>()->getCoefficients();
    _orbitalEnergies = _system->getActiveOrbitalController<SCFMode>()->getEigenvalues();

    unsigned int iSpin = 0;
    for_spin(_coefficients,_orbitalEnergies,_occupation,_indexWhiteList,nOcc){
      // Update orbitals
      auto oldCoefficients = _coefficients_spin;
      auto oldOrbitalEnergies = _orbitalEnergies_spin;
      _coefficients_spin.setZero();

      //Adapt occupation vector and orbital energies
      _occupation_spin.resize(_indexWhiteList_spin.size());
      _orbitalEnergies_spin.resize(_indexWhiteList_spin.size());
      _occupation_spin.setZero();

      for (unsigned int iMO = 0; iMO < _indexWhiteList_spin.size(); ++iMO) {
        _coefficients_spin.col(iMO) = oldCoefficients.col(_indexWhiteList_spin[iMO]);
        _orbitalEnergies_spin(iMO) = oldOrbitalEnergies(_indexWhiteList_spin[iMO]);
        _occupation_spin(iMO) = (_indexWhiteList_spin[iMO] < nOcc_spin) ?
            ((SCFMode == Options::SCF_MODES::RESTRICTED) ? 2 : 1) : 0;
      }/* Update orbitals */

      //Print info
      if (SCFMode == Options::SCF_MODES::RESTRICTED) {
        std::cout<<" System: "<<system+1<<" \n";
        printf(" NEW Reference orbitals : \n");
      } else {
        std::cout<<" System: "<<system+1<<" \n";
        printf("%s NEW Reference orbitals : \n",(iSpin == 0) ? "Alpha" : "Beta");
      }
      for (unsigned int iMO = 0; iMO < _indexWhiteList_spin.size(); ++iMO) {
        printf("%4i",_indexWhiteList_spin[iMO]+1);
        if ( (iMO+1) % 10 == 0) printf("\n");
      }
      printf("\n");
      iSpin += 1;
    };

    this->indexToH5(type);
  }

  template<Options::SCF_MODES SCFMode>
  void LRSCFController<SCFMode>::indexToH5(Options::LRSCF_TYPE type){

    std::string mode = (SCFMode==RESTRICTED) ? "res":"unres";
    std::string type_str = (type==Options::LRSCF_TYPE::ISOLATED || type==Options::LRSCF_TYPE::UNCOUPLED) ? "" : ".fdec";
    std::string name = _system->getSettings().path + _system->getSettings().name + type_str + ".lrscfSpace." + mode + ".h5";
    HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
    HDF5::save_scalar_attribute(file,"ID",_system->getSettings().identifier);
    unsigned int iCount = 0;
    for_spin(_indexWhiteList) {
      std::string spin = "alpha";
      if (iCount > 0) spin = "beta";
      iCount += 1;

      //Cannot use unsigned vector for Eigen::Map
      std::vector<int> tmp(_indexWhiteList_spin.begin(),_indexWhiteList_spin.end());
      EigenHDF5::save(file,spin,Eigen::VectorXi::Map(tmp.data(), tmp.size()));
    };
    file.close();
  }

  template class LRSCFController<Options::SCF_MODES::RESTRICTED>;
  template class LRSCFController<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
