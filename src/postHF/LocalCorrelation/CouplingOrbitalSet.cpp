/**
 * @file CouplingOrbitalSet.cpp
 *
 * @date 1 Jun 2019
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
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"
/* Include Serenity Internal Headers */
#include "data/SingleSubstitution.h"
#include "io/HDF5.h"
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h"
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {

CouplingOrbitalSet::CouplingOrbitalSet(std::shared_ptr<OrbitalPair> parentPair, std::shared_ptr<OrbitalPair> ikPair,
                                       std::shared_ptr<OrbitalPair> kjPair, unsigned int k)
  : _k(k), _ijPair(parentPair), _ikPair(ikPair), _kjPair(kjPair), _ik_eq_ij(parentPair == ikPair), _kj_eq_ij(parentPair == kjPair) {
}

void CouplingOrbitalSet::flushIntegrals() {
  ia_kc.resize(0, 0);
  ja_kc.resize(0, 0);
  ik_ca.resize(0, 0);
  jk_ca.resize(0, 0);
  std::vector<Eigen::MatrixXd>().swap(ka_bc);
  std::vector<Eigen::MatrixXd>().swap(ab_kcX2_M_ak_bc);
}

void CouplingOrbitalSet::loadIntegralsFromFile(HDF5::H5File& file,
                                               std::map<FOUR_CENTER_INTEGRAL_TYPE, std::string>& groupNames, std::string id) {
  std::string kString = std::to_string(_k) + "_";
  HDF5::dataset_exists(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ia_kc] + id);
  HDF5::dataset_exists(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ja_kc] + id);
  HDF5::dataset_exists(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ik_ca] + id);
  HDF5::dataset_exists(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::jk_ca] + id);

  EigenHDF5::load(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ia_kc] + id, ia_kc);
  EigenHDF5::load(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ja_kc] + id, ja_kc);
  EigenHDF5::load(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ik_ca] + id, ik_ca);
  EigenHDF5::load(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::jk_ca] + id, jk_ca);
  unsigned int nPNOs_ij = _ijPair.lock()->t_ij.rows();
  load_vectorSet(file, ka_bc, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ka_bc] + id, nPNOs_ij, nPNOs_ij);
  if (_loadSinglesSigmaInts) {
    load_vectorSet(file, ab_kcX2_M_ak_bc, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ab_kcX2_M_ak_bc] + id,
                   nPNOs_ij, nPNOs_ij);
  }
}

void CouplingOrbitalSet::writeIntegralsToFile(HDF5::H5File& file,
                                              std::map<FOUR_CENTER_INTEGRAL_TYPE, std::string>& groupNames, std::string id) {
  std::string kString = std::to_string(_k) + "_";
  assert(ia_kc.rows() > 0);
  assert(ja_kc.rows() > 0);
  assert(ik_ca.rows() > 0);
  assert(jk_ca.rows() > 0);
  assert(ka_bc.size() > 0);

  EigenHDF5::save(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ia_kc] + id, ia_kc);
  EigenHDF5::save(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ja_kc] + id, ja_kc);
  EigenHDF5::save(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ik_ca] + id, ik_ca);
  EigenHDF5::save(file, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::jk_ca] + id, jk_ca);
  write_vectorSet(file, ka_bc, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ka_bc] + id);
  if (ab_kcX2_M_ak_bc.size() > 0) {
    _loadSinglesSigmaInts = true;
    write_vectorSet(file, ab_kcX2_M_ak_bc, kString + groupNames[FOUR_CENTER_INTEGRAL_TYPE::ab_kcX2_M_ak_bc] + id);
  }
}

void CouplingOrbitalSet::write_vectorSet(HDF5::H5File& file, std::vector<Eigen::MatrixXd>& ints, std::string name) {
  unsigned int nRows = ints[0].rows() * ints[0].cols();
  unsigned int nCols = ints.size();
  Eigen::MatrixXd tmpRep(nRows, nCols);
  for (unsigned int iCol = 0; iCol < nCols; ++iCol) {
    tmpRep.col(iCol) = Eigen::Map<Eigen::VectorXd>(ints[iCol].data(), nRows);
  } // for iCol
  EigenHDF5::save(file, name, tmpRep);
}

void CouplingOrbitalSet::load_vectorSet(HDF5::H5File& file, std::vector<Eigen::MatrixXd>& ints, std::string name,
                                        unsigned int nRows, unsigned int nCols) {
  assert(ints.size() == 0);
  HDF5::dataset_exists(file, name);
  Eigen::MatrixXd tmpRep;
  EigenHDF5::load(file, name, tmpRep);
  unsigned int nVectors = tmpRep.cols();
  for (unsigned int iCol = 0; iCol < nVectors; ++iCol) {
    Eigen::MatrixXd newMatrix = Eigen::Map<Eigen::MatrixXd>(tmpRep.col(iCol).data(), nRows, nCols);
    ints.push_back(newMatrix);
  }
}

double CouplingOrbitalSet::getMemoryRequirement(bool sigmaInts) {
  double memorySize = 0.0;
  const unsigned int nPNOs_ij = _ijPair.lock()->toPAODomain.cols();
  const unsigned int nPNOs_kj = _ikPair.lock()->toPAODomain.cols();
  const unsigned int nPNOs_ik = _kjPair.lock()->toPAODomain.cols();
  const unsigned int nPNOs_k = this->getKSingles()->toPAODomain.cols();
  //(ia|kc),(ja|kc),(ik|ca),(jk|ca)
  memorySize += 2 * nPNOs_ij * (nPNOs_kj + nPNOs_ik);
  //(ij|ak),(ja|ik),(ia|jk)
  memorySize += 3 * nPNOs_ij;
  //(ka|bc)
  memorySize += nPNOs_ij * nPNOs_ij * nPNOs_ij;
  if (sigmaInts) {
    // 2(ij|ak)-(ia|jk) and 2(ij|ak)-(ja|ik)
    memorySize += 2 * nPNOs_k;
    // 2(ab|kc)-(ak|bc)
    memorySize += nPNOs_ij * nPNOs_ij * nPNOs_k;
  }
  return memorySize * sizeof(double);
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_ij_kj() {
  if (!_s_ij_kj) {
    _s_ij_kj = _domainSController.lock()->getS(_ijPair.lock(), _kjPair.lock());
  }
  return *_s_ij_kj;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_ij_ik() {
  if (!_s_ij_ik) {
    _s_ij_ik = _domainSController.lock()->getS(_ijPair.lock(), _ikPair.lock());
  }
  return *_s_ij_ik;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_ij_k() {
  if (!_s_ij_k) {
    _s_ij_k = _domainSController.lock()->getS(_ijPair.lock(), getKSingles());
  }
  return *_s_ij_k;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_kj_ik() {
  if (!_s_kj_ik) {
    _s_kj_ik = _domainSController.lock()->getS(_kjPair.lock(), _ikPair.lock());
  }
  return *_s_kj_ik;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_ik_kj() {
  if (!_s_ik_kj) {
    _s_ik_kj = _domainSController.lock()->getS(_ikPair.lock(), _kjPair.lock());
  }
  return *_s_ik_kj;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_ik_j() {
  if (!_s_ik_j) {
    _s_ik_j = _domainSController.lock()->getS(_ikPair.lock(), _ijPair.lock()->singles_j);
  }
  return *_s_ik_j;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_ik_i() {
  if (!_s_ik_i) {
    _s_ik_i = _domainSController.lock()->getS(_ikPair.lock(), _ijPair.lock()->singles_i);
  }
  return *_s_ik_i;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_kj_i() {
  if (!_s_kj_i) {
    _s_kj_i = _domainSController.lock()->getS(_kjPair.lock(), _ijPair.lock()->singles_i);
  }
  return *_s_kj_i;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_kj_j() {
  if (!_s_kj_j) {
    _s_kj_j = _domainSController.lock()->getS(_kjPair.lock(), _ijPair.lock()->singles_j);
  }
  return *_s_kj_j;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_i_k() {
  if (!_s_i_k) {
    _s_i_k = _domainSController.lock()->getS(_ijPair.lock()->singles_i, getKSingles());
  }
  return *_s_i_k;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_j_k() {
  if (!_s_j_k) {
    _s_j_k = _domainSController.lock()->getS(_ijPair.lock()->singles_j, getKSingles());
  }
  return *_s_j_k;
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_kj_k() {
  auto kjPair = _kjPair.lock();
  return (kjPair->i == _k) ? kjPair->getS_ij_i() : kjPair->getS_ij_j();
}
const Eigen::MatrixXd& CouplingOrbitalSet::getS_ik_k() {
  auto ikPair = _ikPair.lock();
  return (ikPair->i == _k) ? ikPair->getS_ij_i() : ikPair->getS_ij_j();
}

void CouplingOrbitalSet::setOverlapMatrixController(std::shared_ptr<DomainOverlapMatrixController> domainSController) {
  _domainSController = domainSController;
}

} /* namespace Serenity */
