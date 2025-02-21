/**
 * @file MO3CenterIntegralController.cpp
 *
 * @date May 7, 2019
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
#include "integrals/MO3CenterIntegralController.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"               //Loop shells.
#include "basis/BasisController.h"     //Basis controller - AO integral calculation.
#include "data/PAOController.h"        //PAOController definition.
#include "integrals/wrappers/Libint.h" //Integral calculation.
#include "io/FormattedOutputStream.h"  //Filtered output.
#include "io/HDF5.h"                   //Dumping integrals to disk.
#include "io/IOOptions.h"              //Warnings.
#include "misc/HelperFunctions.h"      //Sparse projections from sparse maps.
#include "misc/WarningTracker.h"       //Warnings.
/* Include Std and External Headers */
#if __linux__ || __unix__ || __unix
#include <malloc.h> //Free unused memory.
#endif
#include <sys/stat.h> //Check for an already existing file on disk.

namespace Serenity {

MO3CenterIntegralController::MO3CenterIntegralController(std::shared_ptr<BasisController> auxiliaryBasisController,
                                                         std::shared_ptr<BasisController> basisController,
                                                         const std::shared_ptr<SparseMapsController> sparseMaps,
                                                         std::shared_ptr<PAOController> paoController,
                                                         std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients,
                                                         std::string fBaseName, std::string id, bool triplesMode)
  : _auxiliaryBasisController(auxiliaryBasisController),
    _basisController(basisController),
    _sparseMaps(sparseMaps),
    _virtualCoefficients(std::make_shared<Eigen::MatrixXd>(paoController->getAllPAOs())),
    _occupiedCoefficients(occupiedCoefficients),
    _fBaseName(fBaseName),
    _id(id),
    _triplesMode(triplesMode) {
}

MO3CenterIntegralController::MO3CenterIntegralController(std::shared_ptr<BasisController> auxiliaryBasisController,
                                                         std::shared_ptr<BasisController> basisController,
                                                         const std::shared_ptr<SparseMapsController> sparseMaps,
                                                         std::shared_ptr<Eigen::MatrixXd> virtualCoefficients,
                                                         std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients,
                                                         std::string fBaseName, std::string id, bool triplesMode)
  : _auxiliaryBasisController(auxiliaryBasisController),
    _basisController(basisController),
    _sparseMaps(sparseMaps),
    _virtualCoefficients(virtualCoefficients),
    _occupiedCoefficients(occupiedCoefficients),
    _fBaseName(fBaseName),
    _id(id),
    _triplesMode(triplesMode) {
}

MO3CenterIntegralController::~MO3CenterIntegralController() {
  OutputControl::dOut << "Deleting Cached MO-3-Center integrals!" << std::endl;
  auto types = {MO3CENTER_INTS::ia_K, MO3CENTER_INTS::kl_K, MO3CENTER_INTS::ab_K};
  for (const auto& type : types) {
    if (checkDisk(type)) {
      std::string name = _fBaseName + ".mo3c." + _fileNames[type] + ".h5";
      std::remove(name.c_str());
    } // if checkDisk(type)
  }   // for type
  // Free unused memory.
#if __linux__ || __unix__ || __unix
  malloc_trim(0);
#endif
}

void MO3CenterIntegralController::removeFromMemory(MO3CENTER_INTS mo3CenterType) {
  _diskMode = true;
  if (_integrals[mo3CenterType].first) {
    writeToDisk(mo3CenterType);
    flushIntegrals(mo3CenterType);
    _integrals[mo3CenterType].first = nullptr;
    _integrals[mo3CenterType].second = Eigen::SparseVector<int>(0);
  } // if _integrals[type]
  else {
    WarningTracker::printWarning("Warning: You are trying to remove a set of integrals that is not stored in memory.\n",
                                 iOOptions.printSCFCycleInfo);
  }
}

void MO3CenterIntegralController::flushIntegrals(MO3CENTER_INTS mo3CenterType) {
  *_integrals[mo3CenterType].first = {};
  for (auto& ints : *_integrals[mo3CenterType].first)
    ints.resize(0, 0);
  _integrals[mo3CenterType].first = nullptr;
  _integrals[mo3CenterType].second = Eigen::SparseVector<int>(0);
  // Free unused memory.
#if __linux__ || __unix__ || __unix
  malloc_trim(0);
#endif
}

void MO3CenterIntegralController::setDiskMode(bool diskMode) {
  _diskMode = diskMode;
  if (_diskMode) {
    auto types = {MO3CENTER_INTS::ia_K, MO3CENTER_INTS::kl_K, MO3CENTER_INTS::ab_K};
    for (const auto& type : types) {
      removeFromMemory(type);
    } // for type
  }   // if _diskMode
}

bool MO3CenterIntegralController::checkDisk(MO3CENTER_INTS mo3CenterType) {
  std::string name = _fBaseName + ".mo3c." + _fileNames[mo3CenterType] + ".h5";
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}

void MO3CenterIntegralController::writeToDisk(MO3CENTER_INTS mo3CenterType) {
  MO3CenterIntegrals& ints = *_integrals[mo3CenterType].first;
  std::string name = _fBaseName + ".mo3c." + _fileNames[mo3CenterType] + ".h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save_scalar_attribute(file, "ID", _id);
  for (Eigen::SparseVector<int>::InnerIterator itK(_integrals[mo3CenterType].second); itK; ++itK) {
    unsigned int K = itK.row();
    std::string groupName = std::to_string(K);
    EigenHDF5::save(file, groupName, ints[K]);
  } // for K
  file.close();
}

Eigen::SparseVector<int> MO3CenterIntegralController::getMissingDomain(MO3CENTER_INTS mo3CenterType,
                                                                       const Eigen::SparseVector<int>& kDomain) {
  if (_integrals[mo3CenterType].second.size() == 0 || not _integrals[mo3CenterType].first)
    return kDomain;
  Eigen::SparseMatrix<int> missing(kDomain.rows(), 1);
  Eigen::VectorXi tmp = (Eigen::VectorXi)_integrals[mo3CenterType].second;
  std::vector<Eigen::Triplet<int>> tripletList;
  for (Eigen::SparseVector<int>::InnerIterator itK(kDomain); itK; ++itK) {
    if (tmp(itK.row()) < 1)
      tripletList.push_back(Eigen::Triplet<int>(itK.row(), 0, 1));
  }
  missing.setFromTriplets(tripletList.begin(), tripletList.end());
  return missing;
}

Eigen::SparseVector<int> MO3CenterIntegralController::getUnusedDomain(const Eigen::SparseVector<int>& kDomain) {
  Eigen::VectorXi unsedDomain = Eigen::VectorXi::Constant(kDomain.rows(), 1);
  std::vector<Eigen::Triplet<int>> tripletList;
  for (Eigen::SparseVector<int>::InnerIterator itK(kDomain); itK; ++itK) {
    unsedDomain(itK.row()) = 0;
  }
  return unsedDomain.sparseView();
}

void MO3CenterIntegralController::loadIntegrals(MO3CENTER_INTS mo3CenterType, const Eigen::SparseVector<int>& kDomain) {
  unsigned int nAuxFunc = _auxiliaryBasisController->getNBasisFunctions();
  if (!_integrals[mo3CenterType].first) {
    _integrals[mo3CenterType].first.reset(new MO3CenterIntegrals(nAuxFunc, Eigen::SparseMatrix<double>(0, 0)));
    _integrals[mo3CenterType].second = kDomain;
  }
  else {
    _integrals[mo3CenterType].second += kDomain;
  }
  HDF5::Filepath name(_fBaseName + ".mo3c." + _fileNames[mo3CenterType] + ".h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", _id);
  for (Eigen::SparseVector<int>::InnerIterator itK(kDomain); itK; ++itK) {
    unsigned int K = itK.row();
    Eigen::MatrixXd& ints = (*_integrals[mo3CenterType].first)[K];
    const std::string groupName = std::to_string(K);
    HDF5::dataset_exists(file, groupName);
    EigenHDF5::load(file, groupName, ints);
  } // for K
  file.close();
}

void MO3CenterIntegralController::calculateIntegrals(MO3CENTER_INTS mo3CenterType, const Eigen::SparseVector<int>& kDomain,
                                                     const Eigen::SparseVector<int>& paoDomain) {
  takeTime("AO2MO Exchange Integral Transformation -- Sparse Maps");
  const SparseMap& kToRhoMap = *selectPrescreeningMaps(mo3CenterType).first;
  const SparseMap& kToSigmaMap = *selectPrescreeningMaps(mo3CenterType).second;
  const SparseMap& shellToFunction = _basisController->getFunctionToShellMap();
  const SparseMap kToRhoFuncMap = shellToFunction * kToRhoMap;
  const SparseMap kToSigmaFuncMap = shellToFunction * kToSigmaMap;

  SparseMap kToPAOMap =
      (this->_triplesMode) ? _sparseMaps->getExtendedKtoPAOMap_triples() : _sparseMaps->getExtendedKtoPAOMap();
  if (paoDomain.size() != 0) {
    for (unsigned int col = 0; col < kToPAOMap.cols(); ++col) {
      kToPAOMap.col(col) = (kToPAOMap.col(col).cwiseProduct(paoDomain)).eval();
    }
    _projection_virt = getProjectionMatrices(kToPAOMap);
  }
  const SparseMap kToOccMap = (this->_triplesMode) ? _sparseMaps->getExtendedOccToAuxShellMap_triples().transpose()
                                                   : _sparseMaps->getExtendedOccToAuxShellMap().transpose();
  if (!_integrals[mo3CenterType].first) {
    _integrals[mo3CenterType].first = std::make_shared<MO3CenterIntegrals>();
  }
  const Eigen::MatrixXd& paoCoefficients = *_virtualCoefficients;
  const Eigen::MatrixXd& aoCoefficients = *_occupiedCoefficients;
  this->printInfo(kToRhoMap, kToSigmaMap, kToOccMap, kToPAOMap, mo3CenterType, kDomain);

  auto& libint = Libint::getInstance();
  libint.initialize_plain(
      LIBINT_OPERATOR::coulomb, 3, std::numeric_limits<double>::epsilon(),
      std::max(paoCoefficients.lpNorm<Eigen::Infinity>(), aoCoefficients.lpNorm<Eigen::Infinity>()),
      std::max(_auxiliaryBasisController->getMaxNumberOfPrimitives(), _basisController->getMaxNumberOfPrimitives()));

  unsigned int m = _auxiliaryBasisController->getNBasisFunctions();
  // Resize the result matrices.
  auto& mo3CenterInts = *_integrals[mo3CenterType].first;
  if (mo3CenterInts.empty()) {
    mo3CenterInts = std::vector<Eigen::MatrixXd>(m, Eigen::MatrixXd(0, 0));
  }

  auto auxShells_K = _auxiliaryBasisController->getBasis();
  auto aoShells = _basisController->getBasis();

  bool normAux = !(_auxiliaryBasisController->isAtomicCholesky());

  // Construction of the projection matrices is not allowed within a parallel region.
  // Thus, we will ensure that they are available upon integral calculation.
  this->getProjection(ORBITAL_TYPE::OCCUPIED);
  this->getProjection(ORBITAL_TYPE::VIRTUAL);

#ifdef _OPENMP
  std::vector<Eigen::MatrixXd> integrals(omp_get_max_threads());
  unsigned int nThreads = omp_get_max_threads();
  Eigen::setNbThreads(1);
#else
  std::vector<Eigen::MatrixXd> integrals(1);
#endif

  std::vector<unsigned int> nonZeroIndices;
  for (Eigen::SparseVector<int>::InnerIterator itK(kDomain); itK; ++itK)
    nonZeroIndices.push_back(itK.row());
    // Loop over K
#pragma omp parallel for schedule(dynamic)
  for (unsigned int kCounter = 0; kCounter < nonZeroIndices.size(); ++kCounter) {
    const unsigned int indexK = nonZeroIndices[kCounter];
    const auto& shellK = *auxShells_K[indexK];
    const unsigned int nK = shellK.getNContracted();
    const Eigen::SparseMatrix<double>& paoProjection = *_projection_virt[indexK];
    const Eigen::SparseMatrix<double>& occProjection = *_projection_occ[indexK];
    if (occProjection.nonZeros() == 0 && (mo3CenterType == MO3CENTER_INTS::kl_K || mo3CenterType == MO3CENTER_INTS::ia_K))
      continue;
    if (paoProjection.nonZeros() == 0 && (mo3CenterType == MO3CENTER_INTS::ab_K || mo3CenterType == MO3CENTER_INTS::ia_K))
      continue;
    const Eigen::MatrixXd signPAOCoefficients = (paoCoefficients * paoProjection).eval(); // nBasFunc x nSigPAO
    const Eigen::MatrixXd signOccCoefficients = (aoCoefficients * occProjection).eval();  // nBasFunc x nSigOcc

#ifdef _OPENMP
    const unsigned int threadId = omp_get_thread_num();
#else
    const unsigned int threadId = 0;
#endif
    // Dimensions of storage matrices.
    unsigned int nSignRho = 0;
    unsigned int nSignSigma = 0;
    for (Eigen::SparseMatrix<int>::InnerIterator itR(kToRhoMap, indexK); itR; ++itR) {
      unsigned int shellIndexR = itR.row();
      nSignRho += aoShells[shellIndexR]->getNContracted();
    }
    for (Eigen::SparseMatrix<int>::InnerIterator itS(kToSigmaMap, indexK); itS; ++itS) {
      unsigned int shellIndexS = itS.row();
      nSignSigma += aoShells[shellIndexS]->getNContracted();
    }
    // Storage matrix I_K
    std::vector<Eigen::MatrixXd> i_Ks;
    for (unsigned int K = 0; K < nK; ++K) {
      i_Ks.push_back(Eigen::MatrixXd::Zero(nSignSigma, nSignRho));
    }
    unsigned int shiftedRhoIndex = 0;
    // Calculate AO integrals
    if (omp_get_max_threads() == 1)
      Timings::takeTime(" Local Cor. -    Int. (rs|K) AO");
    for (Eigen::SparseMatrix<int>::InnerIterator itR(kToRhoMap, indexK); itR; ++itR) {
      unsigned int shellIndexR = itR.row();
      const unsigned int nRho = aoShells[shellIndexR]->getNContracted();
      unsigned int shiftedSigmaIndex = 0;
      for (Eigen::SparseMatrix<int>::InnerIterator itS(kToSigmaMap, indexK); itS; ++itS) {
        unsigned int shellIndexS = itS.row();
        const unsigned int nSigma = aoShells[shellIndexS]->getNContracted();
        // Calculate integral, normalize and store in i_K
        if (libint.compute(LIBINT_OPERATOR::coulomb, 0, shellK, *aoShells[shellIndexR], *aoShells[shellIndexS],
                           integrals[threadId], normAux)) {
          for (unsigned int K = 0; K < nK; ++K) {
            unsigned int blockLength = nRho * nSigma;
            unsigned int startRow = K * blockLength;
            Eigen::MatrixXd kResultBlock =
                Eigen::Map<Eigen::MatrixXd>(integrals[threadId].block(startRow, 0, blockLength, 1).data(), nSigma, nRho);
            // Store the block
            i_Ks[K].block(shiftedSigmaIndex, shiftedRhoIndex, nSigma, nRho) = kResultBlock;
          } /* primitives of k -> K */
        }   // if compute
        shiftedSigmaIndex += nSigma;
      } // for itS
      shiftedRhoIndex += nRho;
    } // for itR
    if (omp_get_max_threads() == 1)
      Timings::timeTaken(" Local Cor. -    Int. (rs|K) AO");
    // Perform transformation of the AO coefficients and save results
    if (omp_get_max_threads() == 1)
      Timings::takeTime(" Local Cor. - Int. (rs|K) Tran.");
    // Extract significant AO indices for coefficient matrices.
    const Eigen::SparseMatrix<double> proRho = constructProjectionMatrixFromSparse(kToRhoFuncMap.col(indexK));
    const Eigen::SparseMatrix<double> proSig = constructProjectionMatrixFromSparse(kToSigmaFuncMap.col(indexK));
    const Eigen::MatrixXd c_K = (proRho.transpose() * signOccCoefficients).eval();   // nSigRho x nSigOcc
    const Eigen::MatrixXd pao_K = (proSig.transpose() * signPAOCoefficients).eval(); // nSigSigma x nSigPAO

    unsigned int counter = 0;
    for (const auto& i_K : i_Ks) {
      unsigned int extendedK = _auxiliaryBasisController->extendedIndex(indexK) + counter;
      transformCoefficients(mo3CenterType, mo3CenterInts[extendedK], c_K, pao_K, i_K);
      ++counter;
    } // for i_K
    if (omp_get_max_threads() == 1)
      Timings::timeTaken(" Local Cor. - Int. (rs|K) Tran.");
  } // for indexK
#ifdef _OPENMP
  Eigen::setNbThreads(nThreads);
#endif

  libint.finalize(LIBINT_OPERATOR::coulomb, 0, 3);
  if (_diskMode)
    writeToDisk(mo3CenterType);
  OutputControl::nOut << " done" << std::endl;
  timeTaken(2, "AO2MO Exchange Integral Transformation -- Sparse Maps");
}

inline void MO3CenterIntegralController::transformCoefficients(MO3CENTER_INTS mo3CenterType, Eigen::MatrixXd& result,
                                                               const Eigen::MatrixXd& c_K, const Eigen::MatrixXd& pao_K,
                                                               const Eigen::MatrixXd& i_K) {
  switch (mo3CenterType) {
    case MO3CENTER_INTS::ia_K:
      result.noalias() = (pao_K.transpose() * i_K * c_K).eval();
      break;
    case MO3CENTER_INTS::kl_K:
      result.noalias() = (c_K.transpose() * i_K * c_K).eval();
      break;
    case MO3CENTER_INTS::ab_K:
      result.noalias() = (pao_K.transpose() * i_K * pao_K).eval();
      break;
  } // switch mo3CenterType
}

std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>
MO3CenterIntegralController::getProjectionMatrices(const SparseMap& map) {
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> projectionMatrices;
  for (unsigned int col = 0; col < map.cols(); ++col) {
    projectionMatrices.push_back(
        std::make_shared<Eigen::SparseMatrix<double>>(constructProjectionMatrixFromSparse(map.col(col))));
  } // for col
  return projectionMatrices;
}

std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>
MO3CenterIntegralController::getReducedIndices(const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_redToFullMaps) {
  const unsigned int nAux = k_redToFullMaps.size();
  std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>> indices;
  for (unsigned int kIndex = 0; kIndex < nAux; ++kIndex) {
    const auto& sparseMap = *k_redToFullMaps[kIndex];
    indices.push_back(std::make_shared<std::map<unsigned int, unsigned int>>());
    for (unsigned int iSmall = 0; iSmall < sparseMap.cols(); ++iSmall) {
      for (Eigen::SparseMatrix<double>::InnerIterator itOrb(sparseMap, iSmall); itOrb; ++itOrb) {
        indices[kIndex]->insert(std::make_pair(itOrb.row(), iSmall));
      } // for itOrb
    }   // for iSmall
  }     // for kIndex
  assert(indices.size() == k_redToFullMaps.size());
  return indices;
}

std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> MO3CenterIntegralController::getProjection(ORBITAL_TYPE orbitalType) {
  switch (orbitalType) {
    case ORBITAL_TYPE::OCCUPIED: {
      if (_projection_occ.size() < 1) {
        const SparseMap& map = (this->_triplesMode) ? _sparseMaps->getExtendedOccToAuxShellMap_triples().transpose()
                                                    : _sparseMaps->getExtendedOccToAuxShellMap().transpose();
        _projection_occ = getProjectionMatrices(map);
      }
      return _projection_occ;
    }
    case ORBITAL_TYPE::VIRTUAL: {
      if (_projection_virt.size() < 1) {
        const SparseMap& map =
            (this->_triplesMode) ? _sparseMaps->getExtendedKtoPAOMap_triples() : _sparseMaps->getExtendedKtoPAOMap();
        _projection_virt = getProjectionMatrices(map);
      }
      return _projection_virt;
    }
  }
  return {};
}

std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>> MO3CenterIntegralController::getIndices(ORBITAL_TYPE orbitalType) {
  switch (orbitalType) {
    case ORBITAL_TYPE::OCCUPIED: {
      if (_indices_occ.size() < 1) {
        _indices_occ = this->getReducedIndices(this->getProjection(orbitalType));
      }
      return _indices_occ;
    }
    case ORBITAL_TYPE::VIRTUAL: {
      if (_indices_virt.size() < 1) {
        _indices_virt = this->getReducedIndices(this->getProjection(orbitalType));
      }
      return _indices_virt;
    }
  }
  return {};
}

std::pair<const std::shared_ptr<SparseMap>, const std::shared_ptr<SparseMap>>
MO3CenterIntegralController::selectPrescreeningMaps(MO3CENTER_INTS mo3CenterType) {
  if (!_kToRhoMap)
    _kToRhoMap = (this->_triplesMode) ? std::make_shared<SparseMap>(_sparseMaps->getExtendedKtoRhoMap_triples())
                                      : std::make_shared<SparseMap>(_sparseMaps->getExtendedKtoRhoMap());
  if (!_kToSigmaMap)
    _kToSigmaMap = (this->_triplesMode) ? std::make_shared<SparseMap>(_sparseMaps->getExtendedKtoSigmaMap_triples())
                                        : std::make_shared<SparseMap>(_sparseMaps->getExtendedKtoSigmaMap());
  switch (mo3CenterType) {
    case MO3CENTER_INTS::ia_K:
      return std::make_pair(_kToRhoMap, _kToSigmaMap);
    case MO3CENTER_INTS::kl_K:
      return std::make_pair(_kToRhoMap, _kToRhoMap);
    case MO3CENTER_INTS::ab_K:
      return std::make_pair(_kToSigmaMap, _kToSigmaMap);
  }
  // No default so that the compiler gives a warning.
  assert(false);
  return std::make_pair(nullptr, nullptr);
}

double MO3CenterIntegralController::getMemoryRequirement(MO3CENTER_INTS type, std::shared_ptr<SparseMapsController> sparseMaps,
                                                         std::shared_ptr<BasisController> auxBasisController,
                                                         const Eigen::SparseVector<int>& kDomain, bool triplesMode) {
  const SparseMap& kToPAOMap =
      (triplesMode) ? sparseMaps->getExtendedKtoPAOMap_triples() : sparseMaps->getExtendedKtoPAOMap();
  const SparseMap kToOccMap = (triplesMode) ? sparseMaps->getExtendedOccToAuxShellMap_triples().transpose()
                                            : sparseMaps->getExtendedOccToAuxShellMap().transpose();
  const SparseMap& ktoFirst = (type == MO3CENTER_INTS::kl_K || type == MO3CENTER_INTS::ia_K) ? kToOccMap : kToPAOMap;
  const SparseMap& ktoSecond = (type == MO3CENTER_INTS::ab_K || type == MO3CENTER_INTS::ia_K) ? kToPAOMap : kToOccMap;
  double memory = 0.0;
  auto auxShells_K = auxBasisController->getBasis();
  for (Eigen::SparseVector<int>::InnerIterator itK(kDomain); itK; ++itK) {
    const unsigned int k = itK.row();
    const auto& shellK = *auxShells_K[k];
    const unsigned int nK = shellK.getNContracted();
    memory += ktoFirst.col(k).nonZeros() * ktoSecond.col(k).nonZeros() * nK;
  }
  memory *= sizeof(double);
  return memory;
}

double MO3CenterIntegralController::getTotalMemoryRequirement(std::shared_ptr<SparseMapsController> sparseMaps,
                                                              std::shared_ptr<BasisController> auxBasisController,
                                                              const Eigen::SparseVector<int>& kDomain) {
  double memory = getMemoryRequirement(MO3CENTER_INTS::ia_K, sparseMaps, auxBasisController, kDomain);
  memory += getMemoryRequirement(MO3CENTER_INTS::ab_K, sparseMaps, auxBasisController, kDomain);
  memory += getMemoryRequirement(MO3CENTER_INTS::kl_K, sparseMaps, auxBasisController, kDomain);
  return memory;
}

void MO3CenterIntegralController::printInfo(const SparseMap& kToRhoMap, const SparseMap& kToSigmaMap,
                                            const SparseMap& kToOccMap, const SparseMap& kToPAOMap, MO3CENTER_INTS type,
                                            const Eigen::SparseVector<int>& kDomain) {
  // print general information about the prescreening.
  double nTotalCoeff_AO = kToRhoMap.rows() * kToRhoMap.cols();
  double densRho = ((double)kToRhoMap.nonZeros()) / nTotalCoeff_AO;
  double densSig = ((double)kToSigmaMap.nonZeros()) / nTotalCoeff_AO;
  const SparseMap& ktoFirst = (type == MO3CENTER_INTS::kl_K || type == MO3CENTER_INTS::ia_K) ? kToOccMap : kToPAOMap;
  const SparseMap& ktoSecond = (type == MO3CENTER_INTS::ab_K || type == MO3CENTER_INTS::ia_K) ? kToPAOMap : kToOccMap;
  double nTotalCoeffFirst = ktoFirst.rows() * ktoFirst.cols();
  double nTotalCoeffSecond = ktoSecond.rows() * ktoSecond.cols();
  double densMOFirst = ((double)ktoFirst.nonZeros()) / nTotalCoeffFirst;
  double densMOSecond = ((double)ktoSecond.nonZeros()) / nTotalCoeffSecond;
  std::string integralType = (type == MO3CENTER_INTS::kl_K)   ? "(kl|K)"
                             : (type == MO3CENTER_INTS::ab_K) ? "(ab|K)"
                                                              : "(ia|K)";
  double memory = this->getMemoryRequirement(type, _sparseMaps, _auxiliaryBasisController, kDomain, _triplesMode) * 1e-9;
  OutputControl::vOut << std::string(40, '-') << std::setprecision(4) << std::endl;
  OutputControl::vOut << "         Sparse map densities " << integralType << std::endl;
  OutputControl::vOut << "            " << std::setw(12) << "AO-Indices" << std::setw(12) << "MO-Indices" << std::endl;
  OutputControl::vOut << std::string(40, '-') << std::endl;
  OutputControl::vOut << "K --> first " << std::setw(12) << densRho << std::setw(12) << densMOFirst << std::endl;
  OutputControl::vOut << "K --> second" << std::setw(12) << densSig << std::setw(12) << densMOSecond << std::endl;
  OutputControl::vOut << std::string(40, '-') << std::endl;
  OutputControl::vOut << "Total       " << std::setw(12) << densRho * densSig << std::setw(12)
                      << densMOFirst * densMOSecond << std::endl;
  OutputControl::vOut << std::string(40, '-') << std::endl;
  OutputControl::nOut << "  Memory requirement for integrals " << integralType << ": " << memory << " GB" << std::endl;
  OutputControl::vOut << std::setprecision(6);
  OutputControl::nOut << "  Entering linear scaling integral transformation " << ((_triplesMode) ? "(triples) " : "")
                      << integralType << " ...";
  OutputControl::nOut.flush();
}

void MO3CenterIntegralController::removeIntegralsByDomain(const Eigen::SparseVector<int>& kDomain, MO3CENTER_INTS mo3CenterType) {
  auto auxShells_K = _auxiliaryBasisController->getBasis();
  if (_integrals[mo3CenterType].first) {
    for (Eigen::SparseVector<int>::InnerIterator itK(kDomain); itK; ++itK) {
      const unsigned int kShellIndex = itK.row();
      const auto& shellK = *auxShells_K[kShellIndex];
      const unsigned int nK = shellK.getNContracted();
      unsigned int extendedK = _auxiliaryBasisController->extendedIndex(kShellIndex);
      for (unsigned int k = 0; k < nK; ++k) {
        const unsigned int kFctIndex = extendedK + k;
        (*_integrals[mo3CenterType].first)[kFctIndex].resize(0, 0);
      }
    } // for K
  }
  // Free unused memory.
#if __linux__ || __unix__ || __unix
  malloc_trim(0);
#endif
}

} /* namespace Serenity */
