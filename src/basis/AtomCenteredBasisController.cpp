/**
 * @file AtomCenteredBasisController.cpp
 * @author Thomas Dresselhaus
 *
 * @date Jul 30, 2015
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
#include "basis/AtomCenteredBasisController.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisFunctionProvider.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "io/HDF5.h"
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {

AtomCenteredBasisController::AtomCenteredBasisController(std::shared_ptr<Geometry> geometry,
                                                         const std::string basisLibrary, bool makeSphericalBasis,
                                                         bool makePrimary, const std::string basisLabel, int firstECP,
                                                         std::vector<Eigen::VectorXi> importantShells)
  : BasisController(basisLabel),
    _geometry(geometry),
    _basisLibrary(basisLibrary),
    _basisLabel(basisLabel),
    _makeSphericalBasis(makeSphericalBasis),
    _makePrimary(makePrimary),
    _firstECP(firstECP) {
  assert(_geometry);
  _geometry->addSensitiveObject(this->ObjectSensitiveClass<Geometry>::_self);
  resetBasisLoadingData(importantShells);
  produceBasis();
}

std::unique_ptr<Basis> AtomCenteredBasisController::produceBasisFunctionVector() {
  auto basis = std::unique_ptr<Basis>(new Basis());
  if (_basisLoadingData.size() == 0)
    resetBasisLoadingData();
  // Add new basis functions to the atoms
  unsigned int atomCount = 0;
  for (const auto& atom : _geometry->getAtoms()) {
    if (!atom->basisFunctionsExist(_basisLabel)) {
      BasisFunctionProvider::provideAtomWithBasisFunction(*atom, _basisLibrary, _basisLoadingData[atomCount].label,
                                                          _makeSphericalBasis, _makePrimary, _firstECP);
    }
    // Gather them again to form a Basis.
    auto& thisAtomBasis = atom->getBasisFunctions(_basisLoadingData[atomCount].label);
    if (_basisLoadingData[atomCount].shells.size() != 0) {
      for (unsigned int shell = 0; shell < thisAtomBasis.size(); shell++) {
        if (_basisLoadingData[atomCount].shells[shell])
          basis->push_back(thisAtomBasis[shell]);
      }
    }
    else {
      basis->insert(basis->end(), thisAtomBasis.begin(), thisAtomBasis.end());
      _basisLoadingData[atomCount].shells = Eigen::VectorXi::Ones(thisAtomBasis.size());
    }
    atomCount++;
  }
  return basis;
}

void AtomCenteredBasisController::postConstruction() {
  /*
   * To keep track of the connection between atoms and their basis functions the limiting indices
   * of the basis functions of each atom is calculated and stored here.
   */
  const unsigned int nAtoms = _geometry->getNAtoms();
  _basisIndicesOfAtom.resize(0);
  _basisIndicesRedOfAtom.resize(0);
  _basisIndicesOfAtom.reserve(nAtoms);
  _basisIndicesRedOfAtom.reserve(nAtoms);
  unsigned int basisIndex = 0;
  for (unsigned int atomIndex = 0; atomIndex < _geometry->getNAtoms(); atomIndex++) {
    auto shells = _basisLoadingData[atomIndex].shells;
    unsigned int usedShells = shells.sum();
    unsigned int firstIndexRed = basisIndex;
    // If the current index is larger than the number of shells:
    // 1.) There cannot be any shells left!
    // 2.) The shell index will be identical to the index of the last shell.
    if (firstIndexRed >= _basis->size()) {
      assert(usedShells == 0);
      firstIndexRed = basisIndex - 1;
    }
    unsigned int firstIndex = this->extendedIndex(firstIndexRed);
    unsigned int endIndexRed = 0;
    unsigned int endIndex = 0;
    /*
     * If no shells are left on the atom, set the end index to the first index.
     */
    if (usedShells == 0) {
      endIndex = firstIndex;
      endIndexRed = firstIndexRed;
    }
    else {
      basisIndex += usedShells;
      endIndexRed = basisIndex;
      endIndex = this->extendedIndex(basisIndex - 1) + (*this->_basis)[basisIndex - 1]->getNContracted();
    }
    assert(firstIndex <= endIndex);
    assert(firstIndexRed <= endIndexRed);
    assert(endIndexRed <= _basis->size());
    assert(endIndex <= this->getNBasisFunctions());
    _basisIndicesOfAtom.push_back(std::pair<unsigned int, unsigned int>(firstIndex, endIndex));
    _basisIndicesRedOfAtom.push_back(std::pair<unsigned int, unsigned int>(firstIndexRed, endIndexRed));
  }
  _atomIndicesOfBasisFunctions = std::vector<unsigned int>(this->getNBasisFunctions(), _geometry->getNAtoms());
  _atomBasisProjection = Eigen::MatrixXd::Zero(_geometry->getNAtoms(), this->getNBasisFunctions());
  _atomIndicesOfShells = std::vector<unsigned int>(this->getReducedNBasisFunctions(), _geometry->getNAtoms());
  for (unsigned int atomIndex = 0; atomIndex < _geometry->getNAtoms(); atomIndex++) {
    for (unsigned int iBas = _basisIndicesOfAtom[atomIndex].first; iBas < _basisIndicesOfAtom[atomIndex].second; ++iBas) {
      _atomIndicesOfBasisFunctions[iBas] = atomIndex;
      _atomBasisProjection(atomIndex, iBas) = 1.0;
    } // for iBas
    for (unsigned int iShell = _basisIndicesRedOfAtom[atomIndex].first;
         iShell < _basisIndicesRedOfAtom[atomIndex].second; ++iShell) {
      _atomIndicesOfShells[iShell] = atomIndex;
    } // for iShell
  }   // for atomIndex
  for (auto ibas : _atomIndicesOfBasisFunctions) {
    if (ibas >= _geometry->getNAtoms())
      throw SerenityError("ERROR: While constructing basis function to atom map.");
  }
  for (auto ibas : _atomIndicesOfShells) {
    if (ibas >= _geometry->getNAtoms())
      throw SerenityError("ERROR: While constructing basis function shell to atom map.");
  }
  if (_geometry->hasAtomsWithECPs())
    _geometry->updateCoreCoreRepulsion();
}

void AtomCenteredBasisController::toHDF5(std::string fBaseName, std::string id) {
  std::string name = fBaseName + ".basis.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  for (unsigned int i = 0; i < _geometry->getNAtoms(); i++) {
    std::string count = std::to_string(i);
    HDF5::save_attribute(file, "label_" + count, _basisLoadingData[i].label);
    HDF5::save(file, "shells_" + count, _basisLoadingData[i].shells);
    HDF5::save_attribute(file, "ecp_" + count, _basisLoadingData[i].ecp);
  }
  HDF5::save_scalar_attribute(file, "ID", id);
  file.close();
}

void AtomCenteredBasisController::fromHDF5(std::string fBaseName, std::string id) {
  HDF5::Filepath name(fBaseName + ".basis.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  for (unsigned int i = 0; i < _geometry->getNAtoms(); i++) {
    std::string count = std::to_string(i);
    HDF5::attribute_exists(file, "label_" + count);
    HDF5::dataset_exists(file, "shells_" + count);
    HDF5::attribute_exists(file, "ecp_" + count);
    HDF5::load_scalar_attribute(file, "label_" + count, _basisLoadingData[i].label);
    HDF5::load(file, "shells_" + count, _basisLoadingData[i].shells);
    HDF5::load_scalar_attribute(file, "ecp_" + count, _basisLoadingData[i].ecp);
  }
  // Delete old basis vector in order to work with the basis read from disk.
  this->_basis = nullptr;
  this->BasisController::notify();
  file.close();
}

void AtomCenteredBasisController::notify() {
  this->_basis = nullptr;
  this->BasisController::notify();
  resetBasisLoadingData();
}

void AtomCenteredBasisController::resetBasisLoadingData(const std::vector<Eigen::VectorXi>& importantShells) {
  _basisLoadingData.clear();
  unsigned int nAtom = 0;
  for (auto atom : _geometry->getAtoms()) {
    BasisLoadingData tmp;
    tmp.label = _basisLabel.empty() ? atom->getPrimaryBasisLabel() : _basisLabel;
    tmp.ecp = false;
    tmp.shells = importantShells.size() == 0 ? Eigen::VectorXi(0) : importantShells[nAtom];
    _basisLoadingData.push_back(tmp);
    nAtom++;
  }
}

} /* namespace Serenity */
