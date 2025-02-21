/**
 * @file   Geometry.cpp
 *
 * @date   Mar 19, 2013
 * @author Thomas Dresselhaus
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
#include "geometry/Geometry.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "geometry/Point.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutput.h"
#include "io/HDF5.h"
#include "math/Matrix.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace Serenity {

Geometry::Geometry(std::vector<std::shared_ptr<Atom>> atoms)
  : _atoms(atoms),
    _minX(std::numeric_limits<double>::max()),
    _minY(std::numeric_limits<double>::max()),
    _minZ(std::numeric_limits<double>::max()),
    _maxX(std::numeric_limits<double>::lowest()),
    _maxY(std::numeric_limits<double>::lowest()),
    _maxZ(std::numeric_limits<double>::lowest()),
    _coreCoreRepulsion(0.0),
    _centerOfMass(Point(0.0, 0.0, 0.0)),
    _coreCoreRepulsionUpToDate(false),
    _centerOfMassUpToDate(false) {
  calcCoreCoreRepulsion();
  for (unsigned int i = 0; i < _atoms.size(); ++i) {
    _atoms[i]->addSensitiveObject(this->_self);
    if (_atoms[i]->getX() < _minX || i == 0)
      _minX = _atoms[i]->getX();
    if (_atoms[i]->getY() < _minY || i == 0)
      _minY = _atoms[i]->getY();
    if (_atoms[i]->getZ() < _minZ || i == 0)
      _minZ = _atoms[i]->getZ();
    if (_atoms[i]->getX() > _maxX || i == 0)
      _maxX = _atoms[i]->getX();
    if (_atoms[i]->getY() > _maxY || i == 0)
      _maxY = _atoms[i]->getY();
    if (_atoms[i]->getZ() > _maxZ || i == 0)
      _maxZ = _atoms[i]->getZ();
  }
}

Geometry::Geometry()
  : _atoms(std::vector<std::shared_ptr<Atom>>(0, nullptr)),
    _minX(std::numeric_limits<double>::max()),
    _minY(std::numeric_limits<double>::max()),
    _minZ(std::numeric_limits<double>::max()),
    _maxX(std::numeric_limits<double>::lowest()),
    _maxY(std::numeric_limits<double>::lowest()),
    _maxZ(std::numeric_limits<double>::lowest()),
    _coreCoreRepulsion(0.0),
    _centerOfMass(Point(0.0, 0.0, 0.0)),
    _coreCoreRepulsionUpToDate(false),
    _centerOfMassUpToDate(false) {
}

Geometry::Geometry(std::vector<std::string> atomSymbols, Matrix<double> atomPositions) : Geometry() {
  assert(atomPositions.cols() == 3);
  assert(atomPositions.rows() == (int)atomSymbols.size());
  for (unsigned int i = 0; i < atomSymbols.size(); i++) {
    _atoms.push_back(std::make_shared<Atom>(atomSymbols[i], atomPositions(i, 0), atomPositions(i, 1), atomPositions(i, 2)));
  }
  calcCoreCoreRepulsion();
  for (unsigned int i = 0; i < _atoms.size(); ++i) {
    _atoms[i]->addSensitiveObject(this->_self);
    if (_atoms[i]->getX() < _minX || i == 0)
      _minX = _atoms[i]->getX();
    if (_atoms[i]->getY() < _minY || i == 0)
      _minY = _atoms[i]->getY();
    if (_atoms[i]->getZ() < _minZ || i == 0)
      _minZ = _atoms[i]->getZ();
    if (_atoms[i]->getX() > _maxX || i == 0)
      _maxX = _atoms[i]->getX();
    if (_atoms[i]->getY() > _maxY || i == 0)
      _maxY = _atoms[i]->getY();
    if (_atoms[i]->getZ() > _maxZ || i == 0)
      _maxZ = _atoms[i]->getZ();
  }
}

double Geometry::getCoreCoreRepulsion() const {
  if (!_coreCoreRepulsionUpToDate) {
    calcCoreCoreRepulsion();
  }
  return _coreCoreRepulsion;
}

Point Geometry::getCenterOfMass() const {
  if (!_centerOfMassUpToDate) {
    calcCenterOfMass();
  }
  return _centerOfMass;
}

std::vector<std::string> Geometry::getAtomSymbols() const {
  std::vector<std::string> symbols;
  for (auto& atom : this->_atoms) {
    symbols.push_back(atom->getAtomType()->getElementSymbol());
  }
  return symbols;
}

void Geometry::deleteAtom(unsigned int i) {
  std::vector<std::shared_ptr<Atom>>::iterator it = _atoms.begin();
  std::advance(it, i);
  _atoms.erase(it);
  notify();
};

bool Geometry::hasAtomsWithECPs() const {
  bool atomsWithECPs = false;
  for (auto& atom : _atoms) {
    if (atom->usesECP())
      atomsWithECPs = true;
  }
  return atomsWithECPs;
}

void Geometry::operator+=(const Geometry& rhs) {
  for (auto atom : rhs.getAtoms()) {
    this->_atoms.push_back(atom);
    atom->addSensitiveObject(this->_self);
    if (atom->getX() < this->_minX) {
      this->_minX = atom->getX();
    }
    if (atom->getY() < this->_minY) {
      this->_minY = atom->getY();
    }
    if (atom->getZ() < this->_minZ) {
      this->_minZ = atom->getZ();
    }
    if (atom->getX() > this->_maxX) {
      this->_maxX = atom->getX();
    }
    if (atom->getY() > this->_maxY) {
      this->_maxY = atom->getY();
    }
    if (atom->getZ() > this->_maxZ) {
      this->_maxZ = atom->getZ();
    }
  }
}

bool Geometry::operator==(const Geometry& rhs) {
  bool same;
  if (_atoms.size() == rhs.getAtoms().size()) {
    for (unsigned iAtom = 0; iAtom < _atoms.size(); ++iAtom) {
      same = (*(_atoms[iAtom]) == *(rhs.getAtoms()[iAtom])) ? true : false;
      if (!same)
        break;
    }
  }
  else {
    same = false;
  }
  return same;
}

void Geometry::calcCoreCoreRepulsion() const {
  _coreCoreRepulsion = 0.0;
  for (unsigned int i = 0; i < _atoms.size(); ++i) {
    if (_atoms[i]->isDummy())
      continue;
    for (unsigned int j = 0; j < i; ++j) {
      if (!_atoms[j]->isDummy()) {
        _coreCoreRepulsion +=
            _atoms[i]->getEffectiveCharge() * _atoms[j]->getEffectiveCharge() / distance(*_atoms[i], *_atoms[j]);
      }
    }
  }
  _coreCoreRepulsionUpToDate = true;
}

void Geometry::calcCenterOfMass() const {
  _centerOfMass = Point(0.0, 0.0, 0.0);
  double mass = 0.0;
  double massTimesX = 0.0;
  double massTimesY = 0.0;
  double massTimesZ = 0.0;
  for (auto atom : _atoms) {
    double atomMass = atom->getAtomType()->getMass();
    massTimesX += atomMass * atom->getX();
    massTimesY += atomMass * atom->getY();
    massTimesZ += atomMass * atom->getZ();
    mass += atomMass;
  }
  _centerOfMass.setX(massTimesX / mass);
  _centerOfMass.setY(massTimesY / mass);
  _centerOfMass.setZ(massTimesZ / mass);
  _coreCoreRepulsionUpToDate = true;
}

void Geometry::notify() {
  for (unsigned int i = 0; i < _atoms.size(); i++) {
    if (_atoms[i]->getX() < _minX || i == 0)
      _minX = _atoms[i]->getX();
    if (_atoms[i]->getY() < _minY || i == 0)
      _minY = _atoms[i]->getY();
    if (_atoms[i]->getZ() < _minZ || i == 0)
      _minZ = _atoms[i]->getZ();
    if (_atoms[i]->getX() > _maxX || i == 0)
      _maxX = _atoms[i]->getX();
    if (_atoms[i]->getY() > _maxY || i == 0)
      _maxY = _atoms[i]->getY();
    if (_atoms[i]->getZ() > _maxZ || i == 0)
      _maxZ = _atoms[i]->getZ();
  }
  _coreCoreRepulsionUpToDate = false;
  _centerOfMassUpToDate = false;
  this->notifyObjects();
}

void Geometry::print() const {
  int counter = 0;
  printSmallCaption("Current Geometry (Angstrom)");
  for (auto atom : _atoms) {
    ++counter;
    printf("%4s %4d %2s %+15.10f %+15.10f %+15.10f\n", "", counter, atom->getAtomType()->getName().c_str(),
           atom->getX() * BOHR_TO_ANGSTROM, atom->getY() * BOHR_TO_ANGSTROM, atom->getZ() * BOHR_TO_ANGSTROM);
  }
  std::cout << "\n" << std::endl;
}
void Geometry::printToFile(std::string baseName, std::string id) const {
  std::ofstream file;
  file.open(baseName + ".xyz", std::ofstream::out | std::ofstream::trunc);
  file << _atoms.size() << std::endl;
  file << "ID: " << id << std::endl;
  for (auto atom : _atoms) {
    //    std::string dummy=atom->isDummy()? ":" : ""; // The name already contains the information whether it is a dummy.
    file << atom->getAtomType()->getName() //<< dummy
         << "  " << std::fixed << std::setprecision(12) << atom->getX() * BOHR_TO_ANGSTROM << "  " << std::fixed
         << std::setprecision(12) << atom->getY() * BOHR_TO_ANGSTROM << "  " << std::fixed << std::setprecision(12)
         << atom->getZ() * BOHR_TO_ANGSTROM << std::endl;
  }
  file.close();
}

void Geometry::updateTrajFile(std::string baseName, double energy, double gradNorm) const {
  std::ofstream file;
  file.open(baseName + ".trj", std::ofstream::out | std::ofstream::app);
  file << _atoms.size() << std::endl;
  file << "Energy = " << std::fixed << energy << " Gradient = " << gradNorm << std::endl;
  for (auto atom : _atoms) {
    //    std::string dummy=atom->isDummy()? ":" : ""; // The name already contains the information whether it is a dummy.
    file << atom->getAtomType()->getName() //<< dummy
         << "  " << std::fixed << atom->getX() * BOHR_TO_ANGSTROM << "  " << std::fixed
         << atom->getY() * BOHR_TO_ANGSTROM << "  " << std::fixed << atom->getZ() * BOHR_TO_ANGSTROM << std::endl;
  }
  file.close();
}

Matrix<double> Geometry::getCoordinates() const {
  Matrix<double> coords(this->getNAtoms(), 3);
  for (unsigned int i = 0; i != this->getNAtoms(); ++i) {
    coords(i, 0) = this->_atoms[i]->getX();
    coords(i, 1) = this->_atoms[i]->getY();
    coords(i, 2) = this->_atoms[i]->getZ();
  }
  return coords;
}

Matrix<double> Geometry::getAlignedCoordinates() const {
  // Get coordinates as Eigen3 object
  Matrix<double> coords(this->getNAtoms(), 3);
  for (unsigned int i = 0; i != this->getNAtoms(); ++i) {
    coords(i, 0) = this->_atoms[i]->getX();
    coords(i, 1) = this->_atoms[i]->getY();
    coords(i, 2) = this->_atoms[i]->getZ();
  }

  /* ------------------------
   *   Get 3 Heaviest Atoms
   * ------------------------*/
  auto copy = _atoms;
  std::vector<std::shared_ptr<Atom>> maxAt(3, nullptr);

  for (unsigned int i = 0; i < 3; ++i) {
    unsigned int maxNuc = 0;
    std::vector<std::shared_ptr<Atom>> maxAtTmp;
    for (const auto& atom : copy) {
      unsigned int nuc = atom->getNuclearCharge();
      if (nuc > maxNuc) {
        maxAtTmp.clear();
        maxAtTmp.push_back(atom);
        maxNuc = nuc;
      }
      else if (nuc == maxNuc) {
        maxAtTmp.push_back(atom);
      }
    }
    assert(maxAtTmp.size() > 0);
    double minDist = std::numeric_limits<double>::infinity();
    std::shared_ptr<Atom> minAtom = nullptr;
    for (const auto& atom : maxAtTmp) {
      double dist = distance(*atom, this->getCenterOfMass());
      if (dist < minDist) {
        minAtom = atom;
        minDist = dist;
      }
    }
    assert(minAtom);
    // remove atom from copy
    copy.erase(std::remove(copy.begin(), copy.end(), minAtom), copy.end());
    // add to final list
    maxAt[i] = minAtom;
  }
  Eigen::Vector3d at1;
  Eigen::Vector3d at2;
  Eigen::Vector3d at3;
  at1 << maxAt[0]->getX(), maxAt[0]->getY(), maxAt[0]->getZ();
  at2 << maxAt[1]->getX(), maxAt[1]->getY(), maxAt[1]->getZ();
  at3 << maxAt[2]->getX(), maxAt[2]->getY(), maxAt[2]->getZ();

  /* ---------------------
   *   Move #1 to Origin
   * ---------------------*/
  coords.col(0).array() -= at1[0];
  coords.col(1).array() -= at1[1];
  coords.col(2).array() -= at1[2];
  at2 -= at1;
  at3 -= at1;

  /* --------------------------------
   *   Prep. Rotation of #2 into -z
   * --------------------------------*/
  Eigen::Vector3d mz(0.0, 0.0, -1.0);
  Eigen::Vector3d s = at2.cross(mz);
  double w1 = at2.norm() + mz.dot(at2);
  Eigen::Quaterniond q1(w1, s[0], s[1], s[2]);
  q1.normalize();
  at3 = q1.toRotationMatrix() * at3;

  /* --------------------------------
   *   Prep. Rotation #3 into -z,-y
   * --------------------------------*/
  Eigen::Vector3d my(0.0, -1.0, -0.0);
  at3[2] = 0.0;
  s = at3.cross(my);
  double w2 = at3.norm() + my.dot(at3);
  Eigen::Quaterniond q2(w2, s[0], s[1], s[2]);
  q2.normalize();

  /* -------------------
   *   Apply Rotations
   * ------------------- */
  Eigen::Matrix3d r = q2.toRotationMatrix() * q1.toRotationMatrix();
  for (unsigned int i = 0; i < coords.rows(); ++i) {
    Eigen::Vector3d tmp = coords.row(i);
    coords.row(i) = r * tmp;
  }

  return coords;
}

void Geometry::setCoordinates(const Eigen::MatrixXd& newCoords) const {
  assert(newCoords.rows() == this->getNAtoms());
  assert(newCoords.cols() == 3);
  for (unsigned int i = 0; i < this->getNAtoms(); ++i) {
    this->_atoms[i]->setX(newCoords(i, 0));
    this->_atoms[i]->setY(newCoords(i, 1));
    this->_atoms[i]->setZ(newCoords(i, 2));
  }
  auto& libint = Libint::getInstance();
  libint.clearAllEngines();
}

Matrix<double> Geometry::getGradients() const {
  Matrix<double> grads(this->getNAtoms(), 3);
  for (unsigned int i = 0; i != this->getNAtoms(); ++i) {
    grads.row(i) = this->_atoms[i]->_gradient;
  }
  return grads;
}

void Geometry::setGradients(const Eigen::MatrixXd& newGradients) const {
  assert(newGradients.rows() == this->getNAtoms());
  assert(newGradients.cols() == 3);
  for (unsigned int i = 0; i != this->getNAtoms(); ++i) {
    this->_atoms[i]->_gradient = newGradients.row(i);
    _atoms[i]->_gradientsUpToDate = true;
  }
}

void Geometry::printGradients() const {
  int counter = 0;
  printSmallCaption("Current Geometry Gradients (a.u.)");
  for (auto atom : _atoms) {
    ++counter;
    printf("%4s %4d %2s %+15.10f %+15.10f %+15.10f\n", "", counter, atom->getAtomType()->getElementSymbol().c_str(),
           atom->getGradient()[0], atom->getGradient()[1], atom->getGradient()[2]);
  }
}

bool Geometry::isLinear() {
  bool linear = false;
  if (this->getNAtoms() <= 2) {
    linear = true;
  }
  else if (this->getNAtoms() > 2) {
    linear = true;
    auto coord = this->getCoordinates();
    for (unsigned int i = 2; i < this->getNAtoms(); i++) {
      Eigen::Vector3d vecAB(coord.row(0) - coord.row(1));
      Eigen::Vector3d vecAC(coord.row(0) - coord.row(i));
      double angle = vecAB.dot(vecAC) / (vecAB.norm() * vecAC.norm());
      if (abs(angle) < 0.98) {
        linear = false;
        break;
      }
    }
  }
  return linear;
}

Eigen::MatrixXd Geometry::getTransModes() {
  //  generate translational modes
  auto atoms = this->getAtoms();
  unsigned int nAtoms = this->getNAtoms();
  Eigen::MatrixXd transModes(3 * nAtoms, 3);
  transModes.setZero();
  for (unsigned int i = 0; i != nAtoms; ++i) {
    for (unsigned int j = 0; j != 3; ++j) {
      transModes(i * 3 + j, j) = sqrt(atoms[i]->getAtomType()->getMass());
    }
  }
  // normalizing transModes
  for (unsigned int j = 0; j != 3; ++j) {
    Eigen::Map<Eigen::VectorXd> colVec(transModes.col(j).data(), transModes.rows());
    double normFactor = 1.0 / (sqrt(colVec.dot(colVec)));
    transModes.col(j) *= normFactor;
  }
  return transModes;
}

Eigen::MatrixXd Geometry::getRotModes() {
  //  generate rotational modes
  auto atoms = this->getAtoms();
  unsigned int nAtoms = this->getNAtoms();

  Eigen::MatrixXd IDMat(3, 3);
  IDMat = Eigen::MatrixXd::Identity(3, 3);

  // construct inertia matrix

  Eigen::MatrixXd inertiaMat(3, 3);
  inertiaMat.setZero();
  for (unsigned int i = 0; i != 3; ++i) {
    for (unsigned int j = 0; j != 3; ++j) {
      for (unsigned int k = 0; k != nAtoms; ++k) {
        inertiaMat(i, j) += atoms[k]->getAtomType()->getMass() *
                            (IDMat(i, j) * (atoms[k]->getX() * atoms[k]->getX() + atoms[k]->getY() * atoms[k]->getY() +
                                            atoms[k]->getZ() * atoms[k]->getZ()) -
                             ((i == 0 ? atoms[k]->getX() : (i == 1 ? atoms[k]->getY() : atoms[k]->getZ())) *
                              (j == 0 ? atoms[k]->getX() : (j == 1 ? atoms[k]->getY() : atoms[k]->getZ()))));
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> iDiagSolver(3);
  iDiagSolver.compute(inertiaMat);
  Eigen::VectorXd inertiaEigenvalues = iDiagSolver.eigenvalues();
  auto eVecs = iDiagSolver.eigenvectors();

  // construct rotational normal modes

  Eigen::MatrixXd rotModes(3 * nAtoms, 3);
  rotModes.setZero();
  Eigen::Vector3d x_a;
  x_a.setZero();

  for (unsigned int i = 0; i != 3; ++i) {
    if (inertiaEigenvalues[i] > 0.75) {
      for (unsigned int k = 0; k != nAtoms; ++k) {
        x_a[0] = atoms[k]->getX() -
                 ((atoms[k]->getX() * eVecs(0, i) + atoms[k]->getY() * eVecs(1, i) + atoms[k]->getZ() * eVecs(2, i)) /
                  (eVecs(0, i) * eVecs(0, i) + eVecs(1, i) * eVecs(1, i) + eVecs(2, i) * eVecs(2, i)) * eVecs(0, i));

        x_a[1] = atoms[k]->getY() -
                 ((atoms[k]->getX() * eVecs(0, i) + atoms[k]->getY() * eVecs(1, i) + atoms[k]->getZ() * eVecs(2, i)) /
                  (eVecs(0, i) * eVecs(0, i) + eVecs(1, i) * eVecs(1, i) + eVecs(2, i) * eVecs(2, i)) * eVecs(1, i));

        x_a[2] = atoms[k]->getZ() -
                 ((atoms[k]->getX() * eVecs(0, i) + atoms[k]->getY() * eVecs(1, i) + atoms[k]->getZ() * eVecs(2, i)) /
                  (eVecs(0, i) * eVecs(0, i) + eVecs(1, i) * eVecs(1, i) + eVecs(2, i) * eVecs(2, i)) * eVecs(2, i));

        Eigen::Vector3d colVec;
        colVec = eVecs.col(i);

        auto crossProd = x_a.cross(colVec);
        rotModes((k + 1) * 3 - 3, i) = sqrt(atoms[k]->getAtomType()->getMass()) * crossProd[0];
        rotModes((k + 1) * 3 - 2, i) = sqrt(atoms[k]->getAtomType()->getMass()) * crossProd[1];
        rotModes((k + 1) * 3 - 1, i) = sqrt(atoms[k]->getAtomType()->getMass()) * crossProd[2];
      }
    }
    else {
      for (unsigned int j = 0; j != 3 * nAtoms; ++j) {
        rotModes(j, i) = 0.0;
      }
    }
  }

  for (unsigned int j = 0; j != 3; ++j) {
    Eigen::Map<Eigen::VectorXd> colVec(rotModes.col(j).data(), rotModes.rows());
    double normFactor = (sqrt(colVec.dot(colVec)));
    if (normFactor > 0.00001)
      normFactor = 1.0 / normFactor;
    rotModes.col(j) *= normFactor;
  }

  return rotModes;
}

bool Geometry::hasIdenticalAtoms() const {
  for (unsigned int atom1 = 0; atom1 < this->_atoms.size(); atom1++) {
    for (unsigned int atom2 = 0; atom2 < atom1; atom2++) {
      double dist = distance(*(this->_atoms[atom1]), *(this->_atoms[atom2]));
      if (dist < 1e-5)
        return true;
    }
  }
  return false;
}

void Geometry::deleteIdenticalAtoms() {
  bool deleted = false;
  Eigen::VectorXi deleteAtom = Eigen::VectorXi::Zero(this->_atoms.size());
  for (unsigned int atom1 = 0; atom1 < this->_atoms.size(); atom1++) {
    for (unsigned int atom2 = 0; atom2 < atom1; atom2++) {
      double dist = distance(*(this->_atoms[atom1]), *(this->_atoms[atom2]));
      if (dist < 1e-5) {
        deleted = true;
        if (_atoms[atom1]->isDummy()) {
          deleteAtom[atom1] = 1;
        }
        else if (_atoms[atom2]->isDummy()) {
          deleteAtom[atom2] = 1;
        }
        else {
          deleteAtom[atom1] = 1;
        }
      }
    }
  }
  std::vector<std::shared_ptr<Atom>> newAtomSet;
  for (unsigned int i = 0; i < this->_atoms.size(); i++) {
    if (!deleteAtom[i])
      newAtomSet.push_back(this->_atoms[i]);
  }
  this->_atoms = newAtomSet;
  if (deleted)
    this->notify();
}

void Geometry::deleteGhostAtoms() {
  bool deleted = false;
  Eigen::VectorXi deleteAtom = Eigen::VectorXi::Zero(this->_atoms.size());
  for (unsigned int iAtom = 0; iAtom < this->_atoms.size(); iAtom++) {
    if (this->_atoms[iAtom]->isDummy()) {
      deleteAtom[iAtom] = 1;
      deleted = true;
    }
  }
  std::vector<std::shared_ptr<Atom>> newAtomSet;
  for (unsigned int i = 0; i < this->_atoms.size(); i++) {
    if (!deleteAtom[i])
      newAtomSet.push_back(this->_atoms[i]);
  }
  this->_atoms = newAtomSet;
  if (deleted)
    this->notify();
}

void Geometry::makeGradientsTranslationallyInvariant() {
  auto totalGrad = this->getGradients();
  auto atoms = this->getAtoms();
  unsigned int nAtoms = this->getNAtoms();

  Eigen::MatrixXd transModes = this->getTransModes();

  for (unsigned int j = 0; j != 3; ++j) {
    Eigen::MatrixXd mode = Eigen::Map<const Eigen::MatrixXd>(transModes.col(j).data(), nAtoms, 3);
    totalGrad -= (mode.cwiseProduct(totalGrad).sum() / transModes.col(j).squaredNorm()) * mode;
  }
  this->setGradients(totalGrad);
}

void Geometry::addAsDummy(const Geometry& add, bool toFront) {
  auto addAtoms = add.getAtoms();

  std::vector<std::shared_ptr<Atom>> dummyAtoms;
  for (auto atom : addAtoms) {
    // Dummy atoms should not be added
    if (atom->isDummy())
      continue;
    auto oldAtomType = atom->getAtomType();
    // Create new dummy atomtype
    auto atomType = std::make_shared<AtomType>(
        oldAtomType->getName() + ":", oldAtomType->getNuclearCharge(), oldAtomType->getMass(),
        oldAtomType->getBraggSlaterRadius(), oldAtomType->getVanDerWaalsRadius(), oldAtomType->getUFFRadius(),
        oldAtomType->getNCoreElectrons(), oldAtomType->getOccupations(), oldAtomType->getChemicalHardness(), true);

    std::string basisLabel = atom->getPrimaryBasisLabel();
    auto basisFunctions = atom->getBasisFunctions();

    std::pair<std::string, std::vector<std::shared_ptr<Shell>>> atomBasisFunctions(basisLabel, basisFunctions);

    auto dummyAtom = std::make_shared<Atom>(atomType, atom->getX(), atom->getY(), atom->getZ(), atomBasisFunctions);

    dummyAtoms.push_back(dummyAtom);
  }

  if (toFront) {
    _atoms.insert(_atoms.begin(), dummyAtoms.begin(), dummyAtoms.end());
  }
  else {
    _atoms.insert(_atoms.end(), dummyAtoms.begin(), dummyAtoms.end());
  }
  this->notifyObjects();
}

void Geometry::addDummy(const Geometry& add, bool toFront) {
  auto addAtoms = add.getAtoms();

  std::vector<std::shared_ptr<Atom>> dummyAtoms;
  for (auto atom : addAtoms) {
    // Dummy atoms should not be added
    if (!atom->isDummy())
      continue;
    dummyAtoms.push_back(atom);
  }

  if (toFront) {
    _atoms.insert(_atoms.begin(), dummyAtoms.begin(), dummyAtoms.end());
  }
  else {
    _atoms.insert(_atoms.end(), dummyAtoms.begin(), dummyAtoms.end());
  }
  this->notifyObjects();
}

int Geometry::getTotalEffectiveCharge() {
  int effCharge = 0;
  for (const auto& atom : _atoms)
    effCharge += atom->getEffectiveCharge();
  return effCharge;
}

void Geometry::updateCoreCoreRepulsion() {
  _coreCoreRepulsionUpToDate = false;
}

unsigned int Geometry::getNumberOfCoreElectrons() {
  unsigned int nCoreElectrons = 0;
  for (const auto& atom : _atoms) {
    nCoreElectrons += atom->getNCoreElectrons();
  }
  return nCoreElectrons;
}

double Geometry::getMinimumDistance(const Geometry& rhs) {
  double minimumDistance = std::numeric_limits<double>::infinity();

  const auto& atomsA = this->getAtoms();
  const auto& atomsB = rhs.getAtoms();

  for (const auto& a : atomsA) {
    for (const auto& b : atomsB) {
      minimumDistance = std::min(minimumDistance, distance(*a, *b));
    }
  }

  return minimumDistance;
}

unsigned int Geometry::getNMinimalBasisFunctions(bool excludeDummyAtoms) const {
  unsigned int nMinimalBasisFunctions = 0;
  for (const auto& atom : _atoms) {
    if (excludeDummyAtoms && atom->isDummy())
      continue;
    nMinimalBasisFunctions += atom->getAtomType()->getMinimalBasisSize();
  }
  return nMinimalBasisFunctions;
}

} /* namespace Serenity */
