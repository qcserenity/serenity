/**
 * @file DomainOverlapMatrixController.cpp
 *
 * @author Moritz Bensberg
 * @date Dec 11, 2019
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
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h"
/* Include Std and External Headers */
#include "data/OrbitalPair.h"                           //OrbitalPair definition.
#include "data/PAOController.h"                         //PAOController
#include "data/SingleSubstitution.h"                    //SingleSubstitution definition.
#include "misc/SerenityError.h"                         //Error messages.
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h" //K-Set definition.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"       //KL-Set definition.

namespace Serenity {

DomainOverlapMatrixController::DomainOverlapMatrixController(MatrixInBasis<RESTRICTED> s,
                                                             std::shared_ptr<PAOController> paoController,
                                                             std::vector<std::shared_ptr<OrbitalPair>> closeOrbitalPairs,
                                                             std::vector<std::shared_ptr<SingleSubstitution>> singles,
                                                             const Eigen::MatrixXi closeOrbitalPairIndices, unsigned int nOcc)
  : _s(s),
    _paoController(paoController),
    _closeOrbitalPairs(closeOrbitalPairs),
    _singles(singles),
    _closeOrbitalPairIndices(closeOrbitalPairIndices) {
  unsigned int nPairs = _closeOrbitalPairs.size();
  unsigned int nSingles = _singles.size();
  bool calculateSingles = nSingles > 0;
  _s_ij_kl = std::move(std::unique_ptr<Matrix<std::shared_ptr<Eigen::MatrixXd>>>(
      new Matrix<std::shared_ptr<Eigen::MatrixXd>>(nPairs, nPairs, nullptr)));
  _s_ij_k = std::move(std::unique_ptr<Matrix<std::shared_ptr<Eigen::MatrixXd>>>(
      new Matrix<std::shared_ptr<Eigen::MatrixXd>>(nPairs, nSingles, nullptr)));
  _s_k_l = std::move(std::unique_ptr<Matrix<std::shared_ptr<Eigen::MatrixXd>>>(
      new Matrix<std::shared_ptr<Eigen::MatrixXd>>(nSingles, nSingles, nullptr)));
  _singlesIndices = Eigen::VectorXi::Constant(nOcc, -1);
  for (unsigned int iSingle = 0; iSingle < _singles.size(); ++iSingle) {
    const auto& single = _singles[iSingle];
    _singlesIndices[single->i] = iSingle;
  }
  auto& s_ij_kl = *_s_ij_kl;
  auto& s_ij_k = *_s_ij_k;
  auto& s_k_l = *_s_k_l;
#ifdef _OPENMP
  Eigen::setNbThreads(1);
#endif
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iPair = 0; iPair < _closeOrbitalPairs.size(); ++iPair) {
    std::shared_ptr<OrbitalPair> pair = _closeOrbitalPairs[iPair];
    const Eigen::MatrixXd p_ij_T_S =
        ((_paoController->getPAOsFromDomain(pair->paoDomain) * pair->toPAODomain).transpose() * s);
    unsigned int ijIndex = _closeOrbitalPairIndices(pair->i, pair->j);
    for (auto& kSet : pair->coupledPairs) {
      const std::shared_ptr<OrbitalPair>& ikPair = kSet->getIKPair();
      const std::shared_ptr<OrbitalPair>& kjPair = kSet->getKJPair();
      unsigned int ikIndex = _closeOrbitalPairIndices(ikPair->i, ikPair->j);
      unsigned int kjIndex = _closeOrbitalPairIndices(kjPair->i, kjPair->j);
      const Eigen::MatrixXd p_ik = _paoController->getPAOsFromDomain(ikPair->paoDomain) * ikPair->toPAODomain;
      const Eigen::MatrixXd p_kj = _paoController->getPAOsFromDomain(kjPair->paoDomain) * kjPair->toPAODomain;
      s_ij_kl(ijIndex, ikIndex) = std::make_shared<Eigen::MatrixXd>((p_ij_T_S * p_ik).eval());
      s_ij_kl(ijIndex, kjIndex) = std::make_shared<Eigen::MatrixXd>((p_ij_T_S * p_kj).eval());
      if (calculateSingles) {
        const std::shared_ptr<SingleSubstitution>& kSingle = kSet->getKSingles();
        const Eigen::MatrixXd p_k =
            _paoController->getPAOsFromDomain(kSingle->getDiagonalPair()->paoDomain) * kSingle->toPAODomain;
        s_ij_k(ijIndex, _singlesIndices[kSet->getK()]) = std::make_shared<Eigen::MatrixXd>((p_ij_T_S * p_k).eval());
      }
    } // for kSet
    for (auto klSet : pair->klPairSets) {
      const std::shared_ptr<OrbitalPair>& klPair = klSet->getKLPair();
      int klIndex = _closeOrbitalPairIndices(klPair->i, klPair->j);
      assert(klIndex >= 0);
      if (!s_ij_kl(ijIndex, klIndex)) {
        const Eigen::MatrixXd p_kl = _paoController->getPAOsFromDomain(klPair->paoDomain) * klPair->toPAODomain;
        s_ij_kl(ijIndex, klIndex) = std::make_shared<Eigen::MatrixXd>((p_ij_T_S * p_kl).eval());
      }
      if (calculateSingles) {
        if (!s_ij_k(ijIndex, _singlesIndices[klPair->i])) {
          const Eigen::MatrixXd p_k = _paoController->getPAOsFromDomain(klPair->singles_i->getDiagonalPair()->paoDomain) *
                                      klPair->singles_i->toPAODomain;
          s_ij_k(ijIndex, _singlesIndices[klPair->i]) = std::make_shared<Eigen::MatrixXd>((p_ij_T_S * p_k).eval());
        }
        if (!s_ij_k(ijIndex, _singlesIndices[klPair->j])) {
          const Eigen::MatrixXd p_l = _paoController->getPAOsFromDomain(klPair->singles_j->getDiagonalPair()->paoDomain) *
                                      klPair->singles_j->toPAODomain;
          s_ij_k(ijIndex, _singlesIndices[klPair->j]) = std::make_shared<Eigen::MatrixXd>((p_ij_T_S * p_l).eval());
        }
      }
    } // for klSet
  }   // for pair
  if (calculateSingles) {
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iSingleIndex = 0; iSingleIndex < _singles.size(); ++iSingleIndex) {
      const std::shared_ptr<SingleSubstitution>& iSingle = _singles[iSingleIndex];
      const Eigen::MatrixXd p_i_T_S =
          ((_paoController->getPAOsFromDomain(iSingle->getDiagonalPair()->paoDomain) * iSingle->toPAODomain).transpose() * s);
      s_k_l(_singlesIndices[iSingle->i], _singlesIndices[iSingle->i]) =
          std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Identity(p_i_T_S.rows(), p_i_T_S.rows()));
      for (auto ijPair_ptr : iSingle->orbitalPairs) {
        auto ijPair = ijPair_ptr.lock();
        if (ijPair->i == ijPair->j)
          continue;
        const std::shared_ptr<SingleSubstitution>& jSingle = (ijPair->i == iSingle->i) ? ijPair->singles_j : ijPair->singles_i;
        const Eigen::MatrixXd p_j =
            _paoController->getPAOsFromDomain(jSingle->getDiagonalPair()->paoDomain) * jSingle->toPAODomain;
        const Eigen::MatrixXd s_i_j = p_i_T_S * p_j;
        s_k_l(_singlesIndices[iSingle->i], _singlesIndices[jSingle->i]) = std::make_shared<Eigen::MatrixXd>(s_i_j);
        s_k_l(_singlesIndices[jSingle->i], _singlesIndices[iSingle->i]) =
            std::make_shared<Eigen::MatrixXd>(s_i_j.transpose());
      }
    } // for iSingleIndex
  }
#ifdef _OPENMP
  Eigen::setNbThreads(0);
#endif
}

double DomainOverlapMatrixController::getOverlapMatrixSize() {
  unsigned int memorySize = 0;
  const auto& s_ij_kl = *_s_ij_kl;
  const auto& s_ij_k = *_s_ij_k;
  const auto& s_k_l = *_s_k_l;
  for (unsigned int iCol = 0; iCol < s_ij_kl.cols(); ++iCol) {
    for (unsigned int iRow = 0; iRow < s_ij_kl.rows(); ++iRow) {
      if (s_ij_kl(iRow, iCol))
        memorySize += s_ij_kl(iRow, iCol)->rows() * s_ij_kl(iRow, iCol)->cols();
    } // for iRow
  }   // for iCol
  for (unsigned int iCol = 0; iCol < s_ij_k.cols(); ++iCol) {
    for (unsigned int iRow = 0; iRow < s_ij_k.rows(); ++iRow) {
      if (s_ij_k(iRow, iCol))
        memorySize += s_ij_k(iRow, iCol)->rows() * s_ij_k(iRow, iCol)->cols();
    } // for iRow
  }   // for iCol
  for (unsigned int iCol = 0; iCol < s_k_l.cols(); ++iCol) {
    for (unsigned int iRow = 0; iRow < s_k_l.rows(); ++iRow) {
      if (s_k_l(iRow, iCol))
        memorySize += s_k_l(iRow, iCol)->rows() * s_k_l(iRow, iCol)->cols();
    } // for iRow
  }   // for iCol
  return memorySize * sizeof(double);
}

const std::shared_ptr<Eigen::MatrixXd> DomainOverlapMatrixController::getS(std::shared_ptr<OrbitalPair> ijPair,
                                                                           std::shared_ptr<OrbitalPair> klPair) {
  if (_debugIdentity)
    return _debugIdentity;
  auto& s_ij_kl = *_s_ij_kl;
  int ijIndex = _closeOrbitalPairIndices(ijPair->i, ijPair->j);
  int klIndex = _closeOrbitalPairIndices(klPair->i, klPair->j);
  if (ijIndex < 0 || klIndex < 0)
    throw SerenityError("non existing orbital pair");
  assert(s_ij_kl(ijIndex, klIndex) && "This should have been precomputed S_ij_kl!");
  return s_ij_kl(ijIndex, klIndex);
}
const std::shared_ptr<Eigen::MatrixXd> DomainOverlapMatrixController::getS(std::shared_ptr<OrbitalPair> ijPair,
                                                                           std::shared_ptr<SingleSubstitution> kSingle) {
  if (_debugIdentity)
    return _debugIdentity;
  auto& s_ij_k = *_s_ij_k;
  unsigned int ijIndex = _closeOrbitalPairIndices(ijPair->i, ijPair->j);
  assert(s_ij_k(ijIndex, _singlesIndices[kSingle->i]) && "This should have been precomputed S_ij_k!");
  return s_ij_k(ijIndex, _singlesIndices[kSingle->i]);
}
const std::shared_ptr<Eigen::MatrixXd> DomainOverlapMatrixController::getS(const OrbitalPair& ijPair,
                                                                           std::shared_ptr<SingleSubstitution> kSingle) {
  if (_debugIdentity)
    return _debugIdentity;
  auto& s_ij_k = *_s_ij_k;
  unsigned int ijIndex = _closeOrbitalPairIndices(ijPair.i, ijPair.j);
  assert(s_ij_k(ijIndex, _singlesIndices[kSingle->i]) && "This should have been precomputed S_ij_k!");
  return s_ij_k(ijIndex, _singlesIndices[kSingle->i]);
}
const std::shared_ptr<Eigen::MatrixXd> DomainOverlapMatrixController::getS(std::shared_ptr<SingleSubstitution> kSingle,
                                                                           std::shared_ptr<SingleSubstitution> lSingle) {
  if (_debugIdentity)
    return _debugIdentity;
  auto& s_k_l = *_s_k_l;
  assert(s_k_l(_singlesIndices[kSingle->i], _singlesIndices[lSingle->i]) && "This should have been precomputed S_k_l!");
  return s_k_l(_singlesIndices[kSingle->i], _singlesIndices[lSingle->i]);
}
void DomainOverlapMatrixController::setIdentity(std::shared_ptr<Eigen::MatrixXd> identity) {
  _debugIdentity = identity;
}

} /* namespace Serenity */
