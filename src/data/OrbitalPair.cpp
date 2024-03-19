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
#include "data/SingleSubstitution.h"
#include "io/HDF5.h"                                               //Writing to file.
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"            //Definition of a k-set.
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h" //Overlap matrices.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"                  //Definition of a kl-set.

namespace Serenity {

OrbitalPair::OrbitalPair(unsigned int orbital_i, unsigned int orbital_j, double pnoThreshold,
                         double mainPairEnergyThreshold, double collinearScaling, double sparseMapsConstructionPairScaling)
  : i(orbital_i),
    j(orbital_j),
    _pnoThreshold(pnoThreshold),
    _mainPairEnergyThreshold(mainPairEnergyThreshold),
    _collinearDipoleApproxThreshold(mainPairEnergyThreshold * collinearScaling),
    _sparseMapsConstructionPairThreshold(sparseMapsConstructionPairScaling * mainPairEnergyThreshold) {
  if (sparseMapsConstructionPairScaling == std::numeric_limits<double>::infinity()) {
    _sparseMapsConstructionPairThreshold = mainPairEnergyThreshold;
  }
};

OrbitalPair::~OrbitalPair() {
}

void OrbitalPair::write_acbd(HDF5::H5File& file, std::string id) {
  unsigned int nPNOs = (*ac_bd)(0, 0).rows();
  assert(nPNOs != 0);
  Eigen::MatrixXd tmpRep(nPNOs * nPNOs, nPNOs * (nPNOs + 1) / 2);
  unsigned int counter = 0;
  for (unsigned int a = 0; a < nPNOs; ++a) {
    for (unsigned int b = 0; b <= a; ++b) {
      tmpRep.col(counter) = Eigen::Map<Eigen::VectorXd>((*ac_bd)(b, a).data(), nPNOs * nPNOs);
      ++counter;
    } // for b
  }   // for a
  std::string groupName = _groupNames[FOUR_CENTER_INTEGRAL_TYPE::ac_bd] + id;
  EigenHDF5::save(file, groupName, tmpRep);
}

void OrbitalPair::load_acbd(HDF5::H5File& file, std::string id) {
  std::string groupName = _groupNames[FOUR_CENTER_INTEGRAL_TYPE::ac_bd] + id;
  Eigen::MatrixXd tmpRep;
  EigenHDF5::load(file, groupName, tmpRep);
  HDF5::dataset_exists(file, groupName);
  unsigned int nPNOs = toPAODomain.cols();
  ac_bd = std::make_unique<Matrix<Eigen::MatrixXd>>(nPNOs, nPNOs, Eigen::MatrixXd(nPNOs, nPNOs));
  unsigned int counter = 0;
  for (unsigned int a = 0; a < nPNOs; ++a) {
    for (unsigned int b = 0; b <= a; ++b) {
      (*ac_bd)(b, a) = Eigen::Map<Eigen::MatrixXd>(tmpRep.col(counter).data(), nPNOs, nPNOs);
      ++counter;
    } // for b
  }   // for a
  assert(ac_bd->cols() != 0);
}

void OrbitalPair::writeIntegralsToFile(HDF5::H5File& file) {
  if (type != OrbitalPairTypes::CLOSE)
    return;
  _id = std::to_string(i) + "_" + std::to_string(j);
  write_acbd(file, _id);
  for (auto kSet : coupledPairs)
    kSet->writeIntegralsToFile(file, _groupNames, _id);
  _integralsOnDisk = true;
}

void OrbitalPair::loadIntegralsFromFile(HDF5::H5File& file) {
  if (type != OrbitalPairTypes::CLOSE)
    return;
  assert(_integralsOnDisk && "There were no integrals written to disk!");
  load_acbd(file, _id);
  for (auto kSet : coupledPairs)
    kSet->loadIntegralsFromFile(file, _groupNames, _id);
}

void OrbitalPair::cleanUp() {
  flushIntegrals();
  std::vector<Eigen::MatrixXd>().swap(ic_ab);
  std::vector<Eigen::MatrixXd>().swap(jc_ab);
  ij_ab.resize(0, 0);
  tau_ij.resize(0, 0);
  y_ij.resize(0, 0);
  y_ji.resize(0, 0);
  z_ij.resize(0, 0);
  z_ji.resize(0, 0);
  tf_ab.resize(0, 0);
  residual.resize(0, 0);
  for (auto& klSet : klPairSets)
    klSet->cleanUp();
}

void OrbitalPair::flushIntegrals() {
  if (ac_bd) {
    for (unsigned int a = 0; a < ac_bd->cols(); ++a) {
      for (unsigned int b = 0; b < ac_bd->rows(); ++b) {
        (*ac_bd)(b, a).resize(0, 0);
      }
    }
    ac_bd.reset(nullptr);
  }
  for (auto& kSet : coupledPairs)
    kSet->flushIntegrals();
}

double OrbitalPair::getMemoryRequirement(bool mp2Memory, bool sigmaInts) {
  double memorySize = 0.0;
  const unsigned int nPNOs = toPAODomain.cols();
  //(ia|jb)
  memorySize += nPNOs * nPNOs;
  if (!mp2Memory) {
    if (this->type != OrbitalPairTypes::CLOSE)
      return 0;
    const unsigned int nPNOs_j = singles_j->toPAODomain.cols();
    const unsigned int nPNOs_i = singles_i->toPAODomain.cols();
    //(ac|bd)
    memorySize += nPNOs * (nPNOs + 1) / 2 * nPNOs * nPNOs;
    //(ij|ab),(ia|bc),(ja|bc),(jc|ab),(ic|ab)
    memorySize += nPNOs * nPNOs * (1 + 2 * nPNOs_i + 2 * nPNOs_j);
    memorySize += klPairSets.size() * (2 * nPNOs_j + 2 * nPNOs_i + 2);
    // 2(ia|jb)-(ij|ab)
    if (sigmaInts)
      memorySize += nPNOs_i * nPNOs_j;
  } // if !mp2Memory
  memorySize *= sizeof(double);
  if (!mp2Memory) {
    for (const auto& kSet : coupledPairs) {
      memorySize += kSet->getMemoryRequirement(sigmaInts);
    }
  }
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

double OrbitalPair::getLMP2PairEnergy() {
  double pairEnergy = 0.0;
  if (this->lMP2PairEnergy != 0.0) {
    pairEnergy = this->lMP2PairEnergy;
    pairEnergy += this->deltaPNO;
  }
  else {
    if (this->scMP2PairEnergy != 0.0) {
      pairEnergy = this->scMP2PairEnergy;
    }
    else {
      pairEnergy = this->dipolePairEnergy;
    }
  }
  return pairEnergy;
}

double OrbitalPair::getCCSDPairEnergy() {
  double pairEnergy = 0.0;
  if (this->dlpnoCCSDPairEnergy != 0.0) {
    pairEnergy = this->dlpnoCCSDPairEnergy;
    pairEnergy += this->deltaPNO;
  }
  else {
    if (this->scMP2PairEnergy != 0.0) {
      pairEnergy = this->scMP2PairEnergy;
    }
    else {
      pairEnergy = this->dipolePairEnergy;
    }
  }
  return pairEnergy;
}

double OrbitalPair::getPairEnergy() {
  return (this->dlpnoCCSDPairEnergy != 0.0) ? this->getCCSDPairEnergy() : this->getLMP2PairEnergy();
}

double OrbitalPair::getPNOThreshold() {
  return _pnoThreshold;
}

double OrbitalPair::getPairEnergyThreshold() {
  return _mainPairEnergyThreshold;
}
double OrbitalPair::getCollinearDipolePairThreshold() {
  return _collinearDipoleApproxThreshold;
}

double OrbitalPair::getSparseMapConstructionPairThreshold() {
  return _sparseMapsConstructionPairThreshold;
}

const Eigen::SparseVector<int>& OrbitalPair::getFittingDomain() {
  return _extendedAuxDomain;
}

void OrbitalPair::setFittingDomain(const Eigen::SparseVector<int>& fittingDomain) {
  _extendedAuxDomain = fittingDomain;
}

} /* namespace Serenity */
