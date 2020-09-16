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
#include "basis/BasisController.h"     //Basis controller - AO integral calculation.
#include "data/PAOController.h"        //PAOController definition.
#include "integrals/wrappers/Libint.h" //Integral calculation.
#include "io/FormattedOutputStream.h"  //Filtered output.
#include "io/HDF5.h"                   //Dumping integrals to disk.
#include "io/IOOptions.h"              //Warnings.
#include "misc/HelperFunctions.h"      //Sparse prjections from sparse maps.
#include "misc/WarningTracker.h"       //Warnings.
/* Include Std and External Headers */
#include <sys/stat.h> //Check for an already existing file on disk.

namespace Serenity {

MO3CenterIntegralController::~MO3CenterIntegralController() {
  auto types = {MO3CENTER_INTS::ia_K, MO3CENTER_INTS::kl_K, MO3CENTER_INTS::ab_K};
  for (const auto& type : types) {
    if (checkDisk(type)) {
      std::string name = _fBaseName + ".mo3c." + _fileNames[type] + ".h5";
      std::remove(name.c_str());
    } // if checkDisk(type)
  }   // for type
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
  _integrals[mo3CenterType].first = nullptr;
  _integrals[mo3CenterType].second = Eigen::SparseVector<int>(0);
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
  if (not _integrals[mo3CenterType].first)
    calculateIntegrals(mo3CenterType);
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

void MO3CenterIntegralController::loadOrCalculateIntegrals(MO3CENTER_INTS mo3CenterType, const Eigen::SparseVector<int>& kDomain) {
  if (checkDisk(mo3CenterType)) {
    loadIntegrals(mo3CenterType, kDomain);
  }
  else {
    calculateIntegrals(mo3CenterType);
  }
}

Eigen::SparseVector<int> MO3CenterIntegralController::getMissingDomain(MO3CENTER_INTS mo3CenterType,
                                                                       const Eigen::SparseVector<int>& kDomain) {
  if (_integrals[mo3CenterType].second.size() == 0 || not _integrals[mo3CenterType].first)
    return kDomain;
  Eigen::SparseMatrix<int> missing(kDomain.rows(), 1);
  Eigen::SparseVector<int> tmp = _integrals[mo3CenterType].second - kDomain;
  std::vector<Eigen::Triplet<int>> tripletList;
  for (Eigen::SparseVector<int>::InnerIterator itK(tmp); itK; ++itK) {
    if (itK.value() < 0)
      tripletList.push_back(Eigen::Triplet<int>(itK.row(), 0, 1));
  }
  missing.setFromTriplets(tripletList.begin(), tripletList.end());
  return missing;
}

void MO3CenterIntegralController::loadIntegrals(MO3CENTER_INTS mo3CenterType, const Eigen::SparseVector<int>& kDomain) {
  unsigned int nAuxFunc = _auxilliaryBasisController->getNBasisFunctions();
  if (!_integrals[mo3CenterType].first) {
    _integrals[mo3CenterType].first.reset(new MO3CenterIntegrals(nAuxFunc, Eigen::SparseMatrix<double>(0, 0)));
    _integrals[mo3CenterType].second = kDomain;
  }
  else {
    _integrals[mo3CenterType].second += kDomain;
  }
  HDF5::Filepath name(_fBaseName + ".mo3c." + _fileNames[mo3CenterType] + ".h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
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

void MO3CenterIntegralController::calculateIntegrals(MO3CENTER_INTS mo3CenterType) {
  takeTime("AO2MO Exchange Integral Transformation -- Sparse Maps");
  const SparseMap& kToRhoMap = *selectPrescreeningMaps(mo3CenterType).first;
  const SparseMap& kToSigmaMap = *selectPrescreeningMaps(mo3CenterType).second;
  const SparseMap& kToPAOMap = _sparseMaps->getExtendedKtoPAOMap();
  const SparseMap kToOccMap = _sparseMaps->getExtendedOccToAuxShellMap().transpose();
  auto mo3CenterInts = std::make_shared<MO3CenterIntegrals>();
  _integrals[mo3CenterType].first = mo3CenterInts;
  _integrals[mo3CenterType].second = Eigen::VectorXi::Constant(_auxilliaryBasisController->getBasis().size(), 1).sparseView();
  const Eigen::MatrixXd& paoCoefficients = _paoController->getAllPAOs();
  Eigen::MatrixXd& aoCoefficients = *_occupiedCoefficients;

  // print general information about the prescreening.
  double nTotalCoeff = kToRhoMap.rows() * kToRhoMap.cols();
  OutputControl::vOut << "Sparse map densities" << std::endl;
  OutputControl::vOut << " Sparse map K --> Rho    (first index): " << ((double)kToRhoMap.nonZeros()) / nTotalCoeff
                      << std::endl;
  OutputControl::vOut << " Sparse map K --> Sigma (second index): " << ((double)kToSigmaMap.nonZeros()) / nTotalCoeff
                      << std::endl;

  OutputControl::nOut << "  Entering linear scaling integral transformation ...";
  OutputControl::nOut.flush();

  auto& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::coulomb, 0, 3);

  unsigned int m = _auxilliaryBasisController->getNBasisFunctions();
  // Resize the result matrices.
  *mo3CenterInts = std::vector<Eigen::MatrixXd>(m, Eigen::MatrixXd(paoCoefficients.cols(), aoCoefficients.cols()));

  auto auxShells_K = _auxilliaryBasisController->getBasis();
  auto aoShells = _basisController->getBasis();

#ifdef _OPENMP
  std::vector<Eigen::MatrixXd> integrals(omp_get_max_threads());
  unsigned int nThreads = omp_get_max_threads();
  Eigen::setNbThreads(1);
#else
  std::vector<Eigen::MatrixXd> integrals(1);
#endif

  // Loop over K
#pragma omp parallel for schedule(dynamic)
  for (unsigned int indexK = 0; indexK < auxShells_K.size(); ++indexK) {
    const auto& shellK = *auxShells_K[indexK];
    const unsigned int nK = shellK.getNContracted();
    const Eigen::SparseMatrix<double> paoProjection = constructProjectionMatrixFromSparse(kToPAOMap.col(indexK));
    const Eigen::SparseMatrix<double> occProjection = constructProjectionMatrixFromSparse(kToOccMap.col(indexK));
    const Eigen::MatrixXd signPAOCoefficients = (paoCoefficients * paoProjection).eval();
    const Eigen::MatrixXd signOccCoefficients = (aoCoefficients * occProjection).eval();

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
    Eigen::MatrixXd c_K(nSignRho, signOccCoefficients.cols());
    Eigen::MatrixXd pao_K(nSignSigma, signPAOCoefficients.cols());
    // Loop over basis functions associated with K
    // Storage matrix I_K
    std::vector<Eigen::MatrixXd> i_Ks;
    for (unsigned int K = 0; K < nK; ++K) {
      i_Ks.push_back(Eigen::MatrixXd::Zero(nSignSigma, nSignRho));
    }
    unsigned int shiftetRhoIndex = 0;
    for (Eigen::SparseMatrix<int>::InnerIterator itR(kToRhoMap, indexK); itR; ++itR) {
      unsigned int shellIndexR = itR.row();
      const unsigned int nRho = aoShells[shellIndexR]->getNContracted();
      unsigned int shiftetSigmaIndex = 0;
      for (Eigen::SparseMatrix<int>::InnerIterator itS(kToSigmaMap, indexK); itS; ++itS) {
        unsigned int shellIndexS = itS.row();
        const unsigned int nSigma = aoShells[shellIndexS]->getNContracted();
        // Calculate integral, normalize and store in i_K
        if (libint.compute(libint2::Operator::coulomb, 0, shellK, *aoShells[shellIndexR], *aoShells[shellIndexS],
                           integrals[threadId])) {
          for (unsigned int K = 0; K < nK; ++K) {
            // const unsigned int kk = auxBasis->extendedIndex(indexK)+K;
            unsigned int blockLength = nRho * nSigma;
            unsigned int startRow = K * blockLength;
            Eigen::MatrixXd kResultBlock = Eigen::Map<Eigen::MatrixXd>( //<double, Eigen::Dynamic, Eigen::Dynamic,
                                                                        // Eigen::RowMajor> > (
                integrals[threadId].block(startRow, 0, blockLength, 1).data(), nSigma, nRho);
            // Store the block
            i_Ks[K].block(shiftetSigmaIndex, shiftetRhoIndex, nSigma, nRho) = kResultBlock;
          } /* primitives of k -> K */
        }   // if compute
        // Extract the rows of the coefficient matrices belonging to the
        // significant basis functions. --> PAOs
        unsigned int pao_KRowStart = _basisController->extendedIndex(shellIndexS);
        pao_K.block(shiftetSigmaIndex, 0, nSigma, pao_K.cols()) =
            signPAOCoefficients.block(pao_KRowStart, 0, nSigma, pao_K.cols());
        shiftetSigmaIndex += nSigma;
      } // for itS
      // Extract the rows of the coefficient matrices belonging to the
      // significant basis functions. --> occupied orbitals
      unsigned int c_KRowStart = _basisController->extendedIndex(shellIndexR);
      c_K.block(shiftetRhoIndex, 0, nRho, c_K.cols()) = signOccCoefficients.block(c_KRowStart, 0, nRho, c_K.cols());
      shiftetRhoIndex += nRho;
    } // for itR
    // Perform transformation of the AO coefficients and save results
    unsigned int counter = 0;
    for (const auto& i_K : i_Ks) {
      unsigned int extendedK = _auxilliaryBasisController->extendedIndex(indexK) + counter;
      transformCoefficients(mo3CenterType, (*mo3CenterInts)[extendedK], c_K, pao_K, i_K);
      ++counter;
    } // for i_K
  }   // for indexK
#ifdef _OPENMP
  Eigen::setNbThreads(nThreads);
#endif
  libint.finalize(libint2::Operator::coulomb, 0, 3);
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
      result = (pao_K.transpose() * i_K * c_K).eval();
      break;
    case MO3CENTER_INTS::kl_K:
      result = (c_K.transpose() * i_K * c_K).eval();
      break;
    case MO3CENTER_INTS::ab_K:
      result = (pao_K.transpose() * i_K * pao_K).eval();
      break;
  } // switch mo3CenterType
}

std::vector<Eigen::SparseMatrix<double>> MO3CenterIntegralController::getProjectionMatrices(const SparseMap& map) {
  std::vector<Eigen::SparseMatrix<double>> projectionMatrices;
  for (unsigned int col = 0; col < map.cols(); ++col) {
    projectionMatrices.push_back(constructProjectionMatrixFromSparse(map.col(col)));
  } // for col
  return projectionMatrices;
}

std::vector<Eigen::VectorXi>
MO3CenterIntegralController::getReducedIndices(const std::vector<Eigen::SparseMatrix<double>>& k_redToFullMaps) {
  std::vector<Eigen::VectorXi> indicies;
  for (unsigned int kIndex = 0; kIndex < k_redToFullMaps.size(); ++kIndex) {
    const Eigen::MatrixXd map = Eigen::MatrixXd(k_redToFullMaps[kIndex].transpose());
    Eigen::VectorXi reducedIndices = Eigen::VectorXi::Constant(map.cols(), -1);
    if (map.rows() < 1) {
      indicies.push_back(reducedIndices);
      continue;
    }
    for (unsigned int col = 0; col < map.cols(); ++col) {
      int row;
      double maxCoeff = map.col(col).maxCoeff(&row);
      if (maxCoeff != 0)
        reducedIndices(col) = row;
    } // for col
    indicies.push_back(reducedIndices);
  } // for kIndex
  assert(indicies.size() == k_redToFullMaps.size());
  return indicies;
}

const std::pair<std::vector<Eigen::SparseMatrix<double>>, std::vector<Eigen::VectorXi>>&
MO3CenterIntegralController::getProjectionAndIndices(ORBITAL_TYPE orbitalType) {
  if (!_calculatedIntsProjectionsAndIndices[orbitalType]) {
    SparseMap map;
    switch (orbitalType) {
      case ORBITAL_TYPE::OCCUPIED:
        map = _sparseMaps->getExtendedOccToAuxShellMap().transpose();
        break;
      case ORBITAL_TYPE::VIRTUAL:
        map = _sparseMaps->getExtendedKtoPAOMap();
        break;
      default:
        throw SerenityError{"ORBITAL TYPE was not handled in switch."};
    }
    std::vector<Eigen::SparseMatrix<double>> projection = getProjectionMatrices(map);
    std::vector<Eigen::VectorXi> indices = getReducedIndices(projection);
    _calculatedIntsProjectionsAndIndices[orbitalType] =
        std::make_shared<std::pair<std::vector<Eigen::SparseMatrix<double>>, std::vector<Eigen::VectorXi>>>(projection, indices);
  }
  return *_calculatedIntsProjectionsAndIndices[orbitalType];
}

} /* namespace Serenity */
