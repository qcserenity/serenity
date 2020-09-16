/**
 * @file OrbitalPair.cpp
 *
 * @date May 3, 2019
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
#include "data/OrbitalPair.h"
/* Include Serenity Internal Headers */
#include "io/HDF5.h"                                               //Writing to file.
#include "io/IOOptions.h"                                          //Warnings.
#include "misc/WarningTracker.h"                                   //Warnings
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"            //Definition of a k-set.
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h" //Overlap matrices.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"                  //Definition of a kl-set.

namespace Serenity {

OrbitalPair::~OrbitalPair() {
}

void OrbitalPair::deleteIntegrals() {
  if (_integralsOnDisk) {
    std::remove(("PairIntegrals_" + std::to_string(i) + "_" + std::to_string(j) + ".h5").c_str());
  }
  else {
    WarningTracker::printWarning("WARNING: You are trying to delete integrals which were never written to disk!\n",
                                 iOOptions.printSCFCycleInfo);
  }
}

void OrbitalPair::write_acbd(HDF5::H5File& file) {
  unsigned int nPNOs = (*ac_bd)(0, 0).rows();
  assert(nPNOs != 0);
  Eigen::MatrixXd tmpRep(nPNOs * nPNOs, nPNOs * (nPNOs + 1) / 2);
  unsigned int counter = 0;
  for (unsigned int a = 0; a < nPNOs; ++a) {
    for (unsigned int b = 0; b <= a; ++b) {
      tmpRep.col(counter) = Eigen::Map<Eigen::VectorXd>((*ac_bd)(a, b).data(), nPNOs * nPNOs);
      ++counter;
    } // for b
  }   // for a
  std::string groupName = _groupNames[FOUR_CENTER_INTEGRAL_TYPE::ac_bd];
  EigenHDF5::save(file, groupName, tmpRep);
}

void OrbitalPair::load_acbd(HDF5::H5File& file) {
  std::string groupName = _groupNames[FOUR_CENTER_INTEGRAL_TYPE::ac_bd];
  Eigen::MatrixXd tmpRep;
  EigenHDF5::load(file, groupName, tmpRep);
  HDF5::dataset_exists(file, groupName);
  unsigned int nPNOs = toPAODomain.cols();
  ac_bd.reset(new Matrix<Eigen::MatrixXd>(nPNOs, nPNOs, Eigen::MatrixXd(nPNOs, nPNOs)));
  unsigned int counter = 0;
  for (unsigned int a = 0; a < nPNOs; ++a) {
    for (unsigned int b = 0; b <= a; ++b) {
      (*ac_bd)(a, b) = Eigen::Map<Eigen::MatrixXd>(tmpRep.col(counter).data(), nPNOs, nPNOs);
      (*ac_bd)(b, a) = (*ac_bd)(a, b).transpose().eval();
      ++counter;
    } // for b
  }   // for a
  assert(ac_bd->cols() != 0);
}

void OrbitalPair::writeIntegralsToFile() {
  if (type != OrbitalPairTypes::CLOSE)
    return;
  std::string name = "PairIntegrals_" + std::to_string(i) + "_" + std::to_string(j) + ".h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  write_acbd(file);
  for (auto kSet : coupledPairs)
    kSet->writeIntegralsToFile(file, _groupNames);
  file.close();
  _integralsOnDisk = true;
}

void OrbitalPair::loadIntegralsFromFile() {
  if (type != OrbitalPairTypes::CLOSE)
    return;
  assert(_integralsOnDisk && "There were no integrals written to disk!");
  HDF5::Filepath name("PairIntegrals_" + std::to_string(i) + "_" + std::to_string(j) + ".h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  load_acbd(file);
  for (auto kSet : coupledPairs)
    kSet->loadIntegralsFromFile(file, _groupNames);
  file.close();
}

void OrbitalPair::cleanUp() {
  flushIntegrals();
  ic_ab = {};
  jc_ab = {};
  ij_ab.resize(0, 0);
  tau_ij.resize(0, 0);
  y_ij.resize(0, 0);
  y_ji.resize(0, 0);
  z_ij.resize(0, 0);
  z_ji.resize(0, 0);
  tf_ab.resize(0, 0);
  residual.resize(0, 0);
}

void OrbitalPair::flushIntegrals() {
  ac_bd.reset(nullptr);
  for (auto& kSet : coupledPairs)
    kSet->flushIntegrals();
}

double OrbitalPair::getMemoryRequirement(bool mp2Memory) {
  double memorySize = 0.0;
  const unsigned int nPNOs = toPAODomain.cols();
  //(ia|jb)
  memorySize += nPNOs * nPNOs;
  if (!mp2Memory) {
    const unsigned int nPNOs_j = singles_j->toPAODomain.cols();
    const unsigned int nPNOs_i = singles_i->toPAODomain.cols();
    //(ac|bd)
    memorySize += nPNOs * nPNOs * nPNOs * nPNOs;
    //(ij|ab),(ia|bc),(ja|bc),(jc|ab),(ic|ab)
    memorySize += nPNOs * nPNOs * (1 + 2 * nPNOs_i + 2 * nPNOs_j);
    memorySize += klPairSets.size() * (2 * nPNOs_j + 2 * nPNOs_i + 2);
  } // if !mp2Memory
  memorySize *= sizeof(double);
  if (!mp2Memory)
    for (const auto& kSet : coupledPairs)
      memorySize += kSet->getMemoryRequirement();
  return memorySize;
}
const Eigen::MatrixXd& OrbitalPair::getS_ij_i() {
  if (!_s_ij_i)
    _s_ij_i = _domainSController.lock()->getS(*this, singles_i);
  return *_s_ij_i;
}
const Eigen::MatrixXd& OrbitalPair::getS_ij_j() {
  if (!_s_ij_j)
    _s_ij_j = _domainSController.lock()->getS(*this, singles_j);
  return *_s_ij_j;
}
const Eigen::MatrixXd& OrbitalPair::getS_i_j() {
  if (!_s_i_j)
    _s_i_j = _domainSController.lock()->getS(singles_i, singles_j);
  return *_s_i_j;
}
void OrbitalPair::setOverlapMatrixController(std::shared_ptr<DomainOverlapMatrixController> domainSController) {
  _domainSController = domainSController;
}
} /* namespace Serenity */
