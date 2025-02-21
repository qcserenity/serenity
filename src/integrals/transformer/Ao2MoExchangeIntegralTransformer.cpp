/**
 * @file   Ao2MoExchangeIntegralTransformer.cpp
 *
 * @date   Dec. 23, 2018
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
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h"
/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/QuasiCanonicalPAODomainConstructor.h" //Memory efficient, inplace PNO/PAO construction
#include "basis/AtomCenteredBasisController.h"
#include "data/OrbitalPair.h"
#include "data/OrbitalPairSet.h"
#include "data/PAOController.h"
#include "data/SingleSubstitution.h"
#include "integrals/RI_J_IntegralController.h" //RI_J_IntegralController for schwartz based prescreening transformation.
#include "integrals/transformer/Ao2MoHalfTransformer.h" //Transformation via full four center integrals
#include "integrals/wrappers/Libint.h"                  //Integrals
#include "io/FormattedOutputStream.h"                   //Filtered output.
#include "io/HDF5.h"                                    //Open/close HDF5 file.
#include "memory/MemoryManager.h"                       //Memory handling
#include "misc/SystemSplittingTools.h"                  //Block selection via atom indices
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h" //Definition of an k-set.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"       //KL-pair set definition.

namespace Serenity {

void Ao2MoExchangeIntegralTransformer::calculateTwoCenterIntegrals(MatrixInBasis<Options::SCF_MODES::RESTRICTED>& metric) {
  takeTime("Calc 2-center ints");
  /*
   * Build (Q|1/r|K) = S_aux(Q,K)
   */
  metric.setZero();
  auto auxBasisController = metric.getBasisController();
  bool normAux = !(auxBasisController->isAtomicCholesky());
  auto libint = Libint::getSharedPtr();
  libint->initialize(LIBINT_OPERATOR::coulomb, 0, 2);
  const auto& auxBasis = auxBasisController->getBasis();
  unsigned int nAuxFunctionsRed = auxBasisController->getReducedNBasisFunctions();
#pragma omp parallel for schedule(dynamic)
  for (unsigned int i = 0; i < nAuxFunctionsRed; ++i) {
    const unsigned int pStart = auxBasisController->extendedIndex(i);
    const unsigned int nI = auxBasis[i]->getNContracted();
    const auto& shellA = *auxBasis[i];
    for (unsigned int j = 0; j <= i; ++j) {
      const unsigned int qStart = auxBasisController->extendedIndex(j);
      const unsigned int nJ = auxBasis[j]->getNContracted();
      const auto& shellB = *auxBasis[j];

      Eigen::MatrixXd ints;
      if (libint->compute(LIBINT_OPERATOR::coulomb, 0, shellA, shellB, ints, normAux)) {
        Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(), nJ, nI);
        metric.block(qStart, pStart, nJ, nI) = tmp;
        metric.block(pStart, qStart, nI, nJ) = tmp.transpose();
      }
    }
  }
  libint->finalize(LIBINT_OPERATOR::coulomb, 0, 2);
  timeTaken(3, "Calc 2-center ints");
}

std::vector<Eigen::MatrixXd> Ao2MoExchangeIntegralTransformer::get_abK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
    const Eigen::SparseMatrix<double>& projectionMatrix_a, const Eigen::SparseMatrix<double>& projectionMatrix_b,
    const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps, const Eigen::MatrixXd& toPNO_a,
    const Eigen::MatrixXd& toPNO_b, const MO3CenterIntegrals& abK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nPNOs_a = toPNO_a.cols();
  const unsigned int nPNOs_b = toPNO_b.cols();
  const Eigen::MatrixXd toPNO_b_T = toPNO_b.transpose();
  std::vector<Eigen::MatrixXd> m_abks(nPNOs_a, Eigen::MatrixXd::Zero(nPNOs_b, nLocalAux));
  Eigen::MatrixXd m_abk_tmp(toPNO_b_T.rows(), toPNO_a.cols());
  unsigned int kCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itK(auxDomain); itK; ++itK) {
    unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
    unsigned int nContK = auxBasis[itK.row()]->getNContracted();
    if (k_PAOToFullPAOMaps[itK.row()]->cols() < 1 || abK[extendedK].rows() < 1) {
      kCounter += nContK;
      continue;
    }
    // Transform both indices.
    const Eigen::SparseMatrix<double> finalProjection_b = (projectionMatrix_b * *k_PAOToFullPAOMaps[itK.row()]).pruned();
    const Eigen::SparseMatrix<double> finalProjection_a = (projectionMatrix_a * *k_PAOToFullPAOMaps[itK.row()]).pruned();
    for (unsigned int kk = 0; kk < nContK; ++kk) {
      const Eigen::MatrixXd tmp = (finalProjection_b * abK[extendedK + kk] * finalProjection_a.transpose());
      m_abk_tmp = toPNO_b_T * tmp * toPNO_a;
      // Resort the integrals.
      for (unsigned int aPNO = 0; aPNO < nPNOs_a; ++aPNO) {
        m_abks[aPNO].col(kCounter) = m_abk_tmp.col(aPNO);
      } // for aPNO
      kCounter++;
    } // for kk
  }   // for itK
  return m_abks;
}

std::vector<Eigen::MatrixXd> Ao2MoExchangeIntegralTransformer::get_abK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
    const Eigen::SparseMatrix<double>& projectionMatrix_a, const Eigen::SparseMatrix<double>& projectionMatrix_b,
    const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps, const Eigen::MatrixXd& toPNO_a,
    const MO3CenterIntegrals& abK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nPNOs_a = toPNO_a.cols();
  const unsigned int nPAOs_b = projectionMatrix_b.rows();
  std::vector<Eigen::MatrixXd> m_abks(nPNOs_a, Eigen::MatrixXd::Zero(nPAOs_b, nLocalAux));
  Eigen::MatrixXd m_abk_tmp(nPAOs_b, toPNO_a.cols());
  unsigned int kCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itK(auxDomain); itK; ++itK) {
    unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
    unsigned int nContK = auxBasis[itK.row()]->getNContracted();
    if (k_PAOToFullPAOMaps[itK.row()]->cols() < 1 || abK[extendedK].rows() < 1) {
      kCounter += nContK;
      continue;
    }
    // Transform both indices.
    const Eigen::SparseMatrix<double> finalProjection_b = (projectionMatrix_b * *k_PAOToFullPAOMaps[itK.row()]).pruned();
    const Eigen::SparseMatrix<double> finalProjection_a = (projectionMatrix_a * *k_PAOToFullPAOMaps[itK.row()]).pruned();
    for (unsigned int kk = 0; kk < nContK; ++kk) {
      const Eigen::MatrixXd tmp = (finalProjection_b * abK[extendedK + kk] * finalProjection_a.transpose());
      m_abk_tmp = tmp * toPNO_a;
      // Resort the integrals.
      for (unsigned int aPNO = 0; aPNO < nPNOs_a; ++aPNO) {
        m_abks[aPNO].col(kCounter) = m_abk_tmp.col(aPNO);
      } // for aPNO
      kCounter++;
    } // for kk
  }   // for itK
  return m_abks;
}

std::vector<Eigen::MatrixXd> Ao2MoExchangeIntegralTransformer::get_abK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
    const Eigen::SparseMatrix<double>& projectionMatrix_a, const Eigen::SparseMatrix<double>& projectionMatrix_b,
    const Eigen::SparseMatrix<double>& projectionMatrix_c, const Eigen::SparseMatrix<double>& projectionMatrix_d,
    const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps, const Eigen::MatrixXd& toPNO_a,
    const Eigen::MatrixXd& toPNO_c, const Eigen::MatrixXd& toPNO_d, unsigned int nLocalAuxQ,
    const Eigen::SparseVector<int>& auxDomainQ, const MO3CenterIntegrals& abK, MO3CenterIntegrals& cdK) {
  const Eigen::VectorXi auxDomainQ_dens = (Eigen::VectorXi)auxDomainQ;
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nPNOs_a = toPNO_a.cols();
  const unsigned int nPAOs_b = projectionMatrix_b.rows();
  const unsigned int nPNOs_c = toPNO_c.cols();
  const unsigned int nPNOs_d = toPNO_d.cols();
  const Eigen::SparseMatrix<double> projection_bd = projectionMatrix_d * projectionMatrix_b.transpose();
  const Eigen::SparseMatrix<double> projection_ac = projectionMatrix_c * projectionMatrix_a.transpose();
  std::vector<Eigen::MatrixXd> m_abks(nPNOs_a, Eigen::MatrixXd::Zero(nPAOs_b, nLocalAux));
  cdK = std::vector<Eigen::MatrixXd>(nPNOs_c, Eigen::MatrixXd::Zero(nPNOs_d, nLocalAuxQ));
  Eigen::MatrixXd m_abk_tmp(nPAOs_b, nPNOs_a);
  Eigen::MatrixXd m_abq_tmp(nPNOs_d, nPNOs_c);
  unsigned int kCounter = 0;
  unsigned int qCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itK(auxDomain); itK; ++itK) {
    bool getBlockForQ = auxDomainQ_dens[itK.row()] != 0;
    unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
    unsigned int nContK = auxBasis[itK.row()]->getNContracted();
    if (k_PAOToFullPAOMaps[itK.row()]->cols() < 1 || abK[extendedK].rows() < 1) {
      kCounter += nContK;
      if (getBlockForQ)
        qCounter += nContK;
      continue;
    }
    // Transform both indices.
    const Eigen::SparseMatrix<double> finalProjection_b = (projectionMatrix_b * *k_PAOToFullPAOMaps[itK.row()]).pruned();
    const Eigen::SparseMatrix<double> finalProjection_a = (projectionMatrix_a * *k_PAOToFullPAOMaps[itK.row()]).pruned();
    for (unsigned int kk = 0; kk < nContK; ++kk) {
      const Eigen::MatrixXd tmp = (finalProjection_b * abK[extendedK + kk] * finalProjection_a.transpose());
      m_abk_tmp = tmp * toPNO_a;
      // Resort the integrals.
      for (unsigned int aPNO = 0; aPNO < nPNOs_a; ++aPNO) {
        m_abks[aPNO].col(kCounter) = m_abk_tmp.col(aPNO);
      } // for aPNO
      if (getBlockForQ) {
        const Eigen::MatrixXd tmp2 = projection_bd * tmp * projection_ac.transpose();
        m_abq_tmp = toPNO_d.transpose() * tmp2 * toPNO_c;
        for (unsigned int cPNO = 0; cPNO < nPNOs_c; ++cPNO) {
          cdK[cPNO].col(qCounter) = m_abq_tmp.col(cPNO);
        }
        ++qCounter;
      }
      kCounter++;
    } // for kk
  }   // for itK
  return m_abks;
}

Eigen::MatrixXd Ao2MoExchangeIntegralTransformer::get_iaK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux,
    const Eigen::SparseVector<int>& auxDomain, const Eigen::SparseMatrix<double>& projectionMatrix_a,
    const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
    const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
    const Eigen::MatrixXd& toPNO_a, const unsigned int i, const MO3CenterIntegrals& iaK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nLocalPAOs_a = toPNO_a.rows();
  Eigen::MatrixXd m_iaK = Eigen::MatrixXd::Zero(nLocalPAOs_a, nLocalAux);
  unsigned int qCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itQ(auxDomain); itQ; ++itQ) {
    unsigned int extendedQ = auxBasisController->extendedIndex(itQ.row());
    unsigned int nContQ = auxBasis[itQ.row()]->getNContracted();
    const auto& indices = reducedOccIndices[itQ.row()];
    const auto it_i_red = indices->find(i);
    if (it_i_red == indices->end() || k_PAOToFullPAOMaps[itQ.row()]->cols() < 1 || iaK[extendedQ].rows() < 1) {
      qCounter += nContQ;
      continue;
    }
    const Eigen::SparseMatrix<double> finalPAOProjection = (projectionMatrix_a * *k_PAOToFullPAOMaps[itQ.row()]).eval();
    for (unsigned int qq = 0; qq < nContQ; ++qq) {
      unsigned int totQIndex = extendedQ + qq;
      m_iaK.col(qCounter) = finalPAOProjection * iaK[totQIndex].col(it_i_red->second);
      qCounter++;
    } // for qq
  }   // for itQ
  m_iaK = toPNO_a.transpose() * m_iaK;
  return m_iaK;
}

std::vector<Eigen::MatrixXd> Ao2MoExchangeIntegralTransformer::get_iaK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux,
    const Eigen::SparseVector<int>& auxDomain, const Eigen::SparseMatrix<double>& projectionMatrix_a,
    const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
    const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
    const Eigen::MatrixXd& toPNO_a, const std::map<unsigned int, unsigned int>& iIndices, const MO3CenterIntegrals& iaK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nLocalPAOs_a = toPNO_a.rows();
  std::vector<Eigen::MatrixXd> m_iaKs(iIndices.size(), Eigen::MatrixXd::Zero(nLocalPAOs_a, nLocalAux));

  unsigned int qStart = 0;
  for (Eigen::SparseVector<int>::InnerIterator itQ(auxDomain); itQ; ++itQ) {
    unsigned int extendedQ = auxBasisController->extendedIndex(itQ.row());
    unsigned int nContQ = auxBasis[itQ.row()]->getNContracted();
    const auto& indices = reducedOccIndices[itQ.row()];
    const Eigen::SparseMatrix<double> finalPAOProjection = (projectionMatrix_a * *k_PAOToFullPAOMaps[itQ.row()]).eval();
    for (unsigned int qq = 0; qq < nContQ; ++qq) {
      for (const auto& i : iIndices) {
        // Skip if i is not included in integral set.
        const auto it_i_red = indices->find(i.first);
        if (it_i_red != indices->end() && k_PAOToFullPAOMaps[itQ.row()]->cols() > 0 && iaK[extendedQ].rows() > 0) {
          const unsigned int totQIndex = extendedQ + qq;
          m_iaKs[i.second].col(qStart + qq) = finalPAOProjection * iaK[totQIndex].col(it_i_red->second);
        }
      } // for i
    }   // for qq
    qStart += nContQ;
  } // for itQ
  for (auto& m_iak : m_iaKs) {
    m_iak = toPNO_a.transpose() * m_iak;
  }
  return m_iaKs;
}

Eigen::MatrixXd Ao2MoExchangeIntegralTransformer::get_ikK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
    const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
    const std::map<unsigned int, unsigned int>& kIndices, const unsigned int i, const MO3CenterIntegrals& klK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nK = kIndices.size();
  Eigen::MatrixXd m_ikK = Eigen::MatrixXd::Zero(nK, nLocalAux);
  auto m_ikKptr = m_ikK.data();
  unsigned int qCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itQ(auxDomain); itQ; ++itQ) {
    const auto& indices = reducedOccIndices[itQ.row()];
    const auto it_i_red = indices->find(i);
    unsigned int nContQ = auxBasis[itQ.row()]->getNContracted();
    if (it_i_red == indices->end()) {
      qCounter += nContQ;
      continue;
    }
    unsigned int extendedQ = auxBasisController->extendedIndex(itQ.row());
    for (unsigned int qq = 0; qq < nContQ; ++qq) {
      unsigned int totQIndex = extendedQ + qq;
      const Eigen::VectorXd tmp_i = Eigen::VectorXd(klK[totQIndex].col(it_i_red->second));
      for (auto& k : kIndices) {
        const auto it_k_red = indices->find(k.first);
        if (it_k_red == indices->end()) {
          continue;
        }
        m_ikKptr[qCounter * nK + k.second] = tmp_i(it_k_red->second);
      } // for k
      qCounter++;
    } // for qq
  }   // for itQ
  return m_ikK;
}

Eigen::MatrixXd Ao2MoExchangeIntegralTransformer::get_ikK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
    const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
    const std::vector<unsigned int> kIndices, const unsigned int i, const MO3CenterIntegrals& klK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const size_t nK = kIndices.size();
  Eigen::MatrixXd m_ikK = Eigen::MatrixXd::Zero(nK, nLocalAux);
  auto m_ikKptr = m_ikK.data();
  size_t qCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itQ(auxDomain); itQ; ++itQ) {
    const auto& indices = reducedOccIndices[itQ.row()];
    const auto it_i_red = indices->find(i);
    unsigned int nContQ = auxBasis[itQ.row()]->getNContracted();
    if (it_i_red == indices->end()) {
      qCounter += nContQ;
      continue;
    }
    unsigned int extendedQ = auxBasisController->extendedIndex(itQ.row());
    for (unsigned int qq = 0; qq < nContQ; ++qq) {
      unsigned int totQIndex = extendedQ + qq;
      const Eigen::VectorXd tmp_i = Eigen::VectorXd(klK[totQIndex].col(it_i_red->second));
      size_t kCounter = 0;
      for (auto& k : kIndices) {
        const auto it_k_red = indices->find(k);
        if (it_k_red == indices->end()) {
          ++kCounter;
          continue;
        }
        m_ikKptr[qCounter * nK + kCounter] = tmp_i(it_k_red->second);
        ++kCounter;
      } // for k
      qCounter++;
    } // for qq
  }   // for itQ
  return m_ikK;
}

Eigen::SparseMatrix<double>
Ao2MoExchangeIntegralTransformer::buildSparseAuxProjection(std::shared_ptr<BasisController> auxBasisController,
                                                           const Eigen::SparseVector<int>& auxDomain) {
  const auto& auxBasis = auxBasisController->getBasis();
  unsigned int row = 0;
  std::vector<Eigen::Triplet<double>> projectionTriplets;
  for (Eigen::SparseVector<int>::InnerIterator itQ(auxDomain); itQ; ++itQ) {
    unsigned int nContQ = auxBasis[itQ.row()]->getNContracted();
    unsigned int extendedQ = auxBasisController->extendedIndex(itQ.row());
    for (unsigned int qq = 0; qq < nContQ; ++qq) {
      unsigned int totQIndex = extendedQ + qq;
      projectionTriplets.push_back(Eigen::Triplet<double>(row, totQIndex, 1.0));
      ++row;
    } // for qq
  }   // for itQ
  Eigen::SparseMatrix<double> auxProjectionMatrix(projectionTriplets.size(), auxBasisController->getNBasisFunctions());
  auxProjectionMatrix.setFromTriplets(projectionTriplets.begin(), projectionTriplets.end());
  return auxProjectionMatrix;
}

void Ao2MoExchangeIntegralTransformer::calculateMixedIntegrals(
    std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
    const MatrixInBasis<RESTRICTED>& metric, const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
    const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices, const SparseMap& occToK,
    const MO3CenterIntegrals& iaK, const MO3CenterIntegrals& klK, const Eigen::SparseMatrix<double>& extendedPAODomainProjection,
    const Eigen::SparseVector<int>& extendedAuxDomain, const Eigen::SparseMatrix<double>& extendedAuxProjectionMatrix_T,
    const std::vector<Eigen::MatrixXd>& mEx_amuQ, bool calcSigmaVectorInts) {
  bool oneThread = omp_get_num_threads() == 1;
  const Eigen::MatrixXd& toPNO_ij = pair->toPAODomain;

  /*
   * Extract all integrals ikQ/jkQ, iaQ/jaQ which will be used at some point during
   * the k-Set loop. The integrals used for each k-Set can easily be extracted by
   * sparse projection. The projection should only lead to an insignificant overhead,
   * while the early integral extraction is slightly faster than extracting in the loop
   * itself.
   */
  const Eigen::MatrixXd extendedBlock =
      SystemSplittingTools<RESTRICTED>::getMatrixBlockShellWise(metric, extendedAuxDomain, extendedAuxDomain);

  std::map<unsigned int, unsigned int> kTokIndexMap; // Key: k, value: index in kIndices list.
  unsigned int kIndex = 0;
  for (const auto& kSet : pair->coupledPairs) {
    unsigned int k = kSet->getK();
    if (kTokIndexMap.find(k) == kTokIndexMap.end()) {
      kTokIndexMap.insert(std::make_pair(k, kIndex));
      ++kIndex;
    }
  }
  const unsigned int nExtendedLocalAux = extendedBlock.cols();
  std::vector<Eigen::MatrixXd> tmp_mEx_iaQ;
  std::vector<Eigen::MatrixXd> tmp_mEx_ikQ_T;
  tmp_mEx_iaQ.push_back(get_iaK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, pair->domainProjection,
                                k_PAOToFullPAOMaps, reducedOccIndices, toPNO_ij, pair->i, iaK));
  tmp_mEx_ikQ_T.push_back(
      get_ikK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, reducedOccIndices, kTokIndexMap, pair->i, klK)
          .transpose()
          .eval());
  if (pair->i != pair->j) {
    tmp_mEx_iaQ.push_back(get_iaK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, pair->domainProjection,
                                  k_PAOToFullPAOMaps, reducedOccIndices, toPNO_ij, pair->j, iaK));
    tmp_mEx_ikQ_T.push_back(
        get_ikK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, reducedOccIndices, kTokIndexMap, pair->j, klK)
            .transpose()
            .eval());
  }
  const Eigen::MatrixXd& mEx_iaQ = tmp_mEx_iaQ[0];
  const Eigen::MatrixXd& mEx_jaQ = (pair->i != pair->j) ? tmp_mEx_iaQ[1] : tmp_mEx_iaQ[0];
  const Eigen::MatrixXd& mEx_ikQ_T = tmp_mEx_ikQ_T[0];
  const Eigen::MatrixXd& mEx_jkQ_T = (pair->i != pair->j) ? tmp_mEx_ikQ_T[1] : tmp_mEx_ikQ_T[0];

  MO3CenterIntegrals kaK;
  if (calcSigmaVectorInts) {
    kaK = get_iaK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, pair->domainProjection, k_PAOToFullPAOMaps,
                  reducedOccIndices, toPNO_ij, kTokIndexMap, iaK);
  }

  for (unsigned int ik = 0; ik < pair->coupledPairs.size(); ++ik) {
    auto kSet = pair->coupledPairs[ik];
    const unsigned int k = kSet->getK();
    unsigned int kIndex = kTokIndexMap.find(k)->second;
    if (kTokIndexMap.find(k) == kTokIndexMap.end())
      throw SerenityError("ERROR: Wrong index map! 4");
    const Eigen::SparseVector<int> auxDomain = (occToK.col(pair->i) + occToK.col(pair->j) + occToK.col(k)).pruned().eval();

    const Eigen::SparseMatrix<double> auxProjectionMatrix = buildSparseAuxProjection(auxBasisController, auxDomain);
    const Eigen::SparseMatrix<double> auxIJKtoExtended =
        (auxProjectionMatrix * extendedAuxProjectionMatrix_T).pruned().eval();

    const Eigen::MatrixXd localBlock = (auxIJKtoExtended * extendedBlock * auxIJKtoExtended.transpose()).eval();
    const unsigned int nLocalAux = localBlock.cols();
    const Eigen::LLT<Eigen::MatrixXd> llt = localBlock.llt();
    Eigen::MatrixXd m_iaQ = mEx_iaQ * auxIJKtoExtended.transpose();
    Eigen::MatrixXd m_jaQ = mEx_jaQ * auxIJKtoExtended.transpose();
    Eigen::VectorXd m_ikQ = auxIJKtoExtended * mEx_ikQ_T.col(kIndex);
    Eigen::VectorXd m_jkQ = auxIJKtoExtended * mEx_jkQ_T.col(kIndex);
    m_ikQ = llt.solve(m_ikQ).eval();
    m_jkQ = llt.solve(m_jkQ).eval();
    m_iaQ = llt.solve(m_iaQ.transpose()).eval();
    m_jaQ = llt.solve(m_jaQ.transpose()).eval();

    const Eigen::MatrixXd& toPNO_ik = kSet->getIKPair()->toPAODomain;
    const Eigen::MatrixXd& toPNO_kj = kSet->getKJPair()->toPAODomain;
    const Eigen::SparseMatrix<double>& domainProjection_ik = kSet->getIKPair()->domainProjection;
    const Eigen::SparseMatrix<double>& domainProjection_kj = kSet->getKJPair()->domainProjection;
    const Eigen::SparseMatrix<double> domainProjection_ijTo_ik =
        domainProjection_ik * extendedPAODomainProjection.transpose();
    const Eigen::SparseMatrix<double> domainProjection_ijTo_kj =
        domainProjection_kj * extendedPAODomainProjection.transpose();

    // Coulomb type integrals.
    kSet->ik_ca = Eigen::MatrixXd::Zero(toPNO_kj.cols(), toPNO_ij.cols());
    kSet->jk_ca = Eigen::MatrixXd::Zero(toPNO_ik.cols(), toPNO_ij.cols());
    for (unsigned int aPNO = 0; aPNO < mEx_amuQ.size(); ++aPNO) {
      const Eigen::MatrixXd m_amuQ = mEx_amuQ[aPNO] * auxIJKtoExtended.transpose();
      kSet->ik_ca.col(aPNO) = toPNO_kj.transpose() * domainProjection_ijTo_kj * m_amuQ * m_ikQ;
      kSet->jk_ca.col(aPNO) = toPNO_ik.transpose() * domainProjection_ijTo_ik * m_amuQ * m_jkQ;
    } // for aPNO

    // Exchange type integrals
    const Eigen::MatrixXd m_kcQ_ik = get_iaK(auxBasisController, nLocalAux, auxDomain, domainProjection_ik,
                                             k_PAOToFullPAOMaps, reducedOccIndices, toPNO_ik, k, iaK);
    const Eigen::MatrixXd m_kcQ_kj = get_iaK(auxBasisController, nLocalAux, auxDomain, domainProjection_kj,
                                             k_PAOToFullPAOMaps, reducedOccIndices, toPNO_kj, k, iaK);
    kSet->ia_kc = (m_kcQ_kj * m_iaQ).transpose().eval();
    kSet->ja_kc = (m_kcQ_ik * m_jaQ).transpose().eval();

    if (calcSigmaVectorInts) {
      if (oneThread)
        Timings::takeTime("sigma ints. over k");
      const auto& kSingles = kSet->getKSingles();
      const Eigen::MatrixXd& toPNO_k = kSingles->toPAODomain;
      const unsigned int nPNO_k = toPNO_k.cols();
      const unsigned int nPNO_ij = toPNO_ij.cols();
      const Eigen::SparseMatrix<double>& domainProjection_k = kSingles->getDiagonalPair()->domainProjection;
      const Eigen::SparseMatrix<double> domainProjection_ExTo_k =
          domainProjection_k * extendedPAODomainProjection.transpose();
      const Eigen::SparseMatrix<double> domainProjection_ExTo_ij =
          pair->domainProjection * extendedPAODomainProjection.transpose();
      const Eigen::MatrixXd m_kcQ_k = get_iaK(auxBasisController, nLocalAux, auxDomain, domainProjection_k,
                                              k_PAOToFullPAOMaps, reducedOccIndices, toPNO_k, k, iaK);
      const Eigen::MatrixXd llt_m_kcK = llt.solve(m_kcQ_k.transpose().eval());
      kSet->ab_kcX2_M_ak_bc = std::vector<Eigen::MatrixXd>(nPNO_k, Eigen::MatrixXd::Zero(nPNO_ij, nPNO_ij));
      // 2(ab|kc)
      for (unsigned int aPNO = 0; aPNO < nPNO_ij; ++aPNO) {
        const Eigen::MatrixXd m_amuQ = mEx_amuQ[aPNO] * auxIJKtoExtended.transpose();
        const Eigen::MatrixXd ab_kcX2 = 2 * toPNO_ij.transpose() * domainProjection_ExTo_ij * m_amuQ * llt_m_kcK;
        for (unsigned int cPNO = 0; cPNO < nPNO_k; ++cPNO) {
          kSet->ab_kcX2_M_ak_bc[cPNO].col(aPNO) = ab_kcX2.col(cPNO);
        } // for cPNO
      }   // for aPNO
      // -(ak|bc)
      const auto m_bcK_k = flipAndTransform(mEx_amuQ, toPNO_k, domainProjection_ExTo_k, auxIJKtoExtended);
      const Eigen::MatrixXd m_kaK = kaK[kIndex] * auxIJKtoExtended.transpose();
      const Eigen::MatrixXd llt_m_kaK = llt.solve(m_kaK.transpose().eval());
      for (unsigned int cPNO = 0; cPNO < nPNO_k; ++cPNO) {
        const Eigen::MatrixXd ak_bc = (m_bcK_k[cPNO] * llt_m_kaK).transpose();
        kSet->ab_kcX2_M_ak_bc[cPNO] -= ak_bc;
      } // for cPNO
      // 2(ij|ak) - (ia|jk)  and 2(ij|ak) - (ja|ik) a in [k]
      std::map<unsigned int, unsigned int> ijMap;
      ijMap.insert(std::make_pair(pair->i, 0));
      if (pair->i != pair->j)
        ijMap.insert(std::make_pair(pair->j, 1));
      const auto m_iaqs = get_iaK(auxBasisController, nLocalAux, auxDomain, domainProjection_k, k_PAOToFullPAOMaps,
                                  reducedOccIndices, toPNO_k, ijMap, iaK);
      const Eigen::MatrixXd& m_iaq = m_iaqs[ijMap.find(pair->i)->second];
      const Eigen::MatrixXd& m_jaq = m_iaqs[ijMap.find(pair->j)->second];
      unsigned int jIndex = kTokIndexMap.find(pair->j)->second;
      const Eigen::VectorXd ijK = auxIJKtoExtended * mEx_ikQ_T.col(jIndex);
      const Eigen::VectorXd invV_ijK = llt.solve(ijK).eval();
      const Eigen::VectorXd ij_ak = m_kcQ_k * invV_ijK;
      kSet->ij_akX2_M_ia_jk = 2 * ij_ak - (m_iaq * m_jkQ).eval();
      kSet->ij_akX2_M_ja_ik = 2 * ij_ak - (m_jaq * m_ikQ).eval();
      if (oneThread)
        Timings::timeTaken("sigma ints. over k");
    }
  } // for kSet
}

void Ao2MoExchangeIntegralTransformer::calculate_3Occ_Integrals(
    std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
    const Eigen::LLT<Eigen::MatrixXd>& llt_metric, const Eigen::SparseVector<int>& pairDomainToK,
    const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
    const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
    const MO3CenterIntegrals& iaK, const MO3CenterIntegrals& klK) {
  const unsigned int nLocalAux = llt_metric.cols();
  Eigen::SparseMatrix<double>& domainProjection_i = pair->singles_i->getDiagonalPair()->domainProjection;
  Eigen::SparseMatrix<double>& domainProjection_j = pair->singles_j->getDiagonalPair()->domainProjection;
  const Eigen::MatrixXd& toPNO = pair->toPAODomain;
  const Eigen::MatrixXd& toPNO_i = pair->singles_i->toPAODomain;
  const Eigen::MatrixXd& toPNO_j = pair->singles_j->toPAODomain;

  // Build list of orbitals k, which we will need in the integral extraction and keep track of their ordering
  // and remove any identical k's. We do this using a #-map.
  std::map<unsigned int, unsigned int> kTokIndexMap; // Key: k, value: index in kIndices list.
  unsigned int kIndex = 0;
  for (const auto& klSet : pair->klPairSets) {
    unsigned int k = klSet->getKLPair()->i;
    unsigned int l = klSet->getKLPair()->j;
    if (kTokIndexMap.find(k) == kTokIndexMap.end()) {
      kTokIndexMap.insert(std::make_pair(k, kIndex));
      ++kIndex;
    }
    if (kTokIndexMap.find(l) == kTokIndexMap.end()) {
      kTokIndexMap.insert(std::make_pair(l, kIndex));
      ++kIndex;
    }
  }

  // These are rather large, so I will use the possible symmetry for i==j.
  std::vector<Eigen::MatrixXd> tmp_mikq_T;
  std::vector<std::vector<Eigen::MatrixXd>> tmp_mkaqs;
  tmp_mikq_T.push_back(
      get_ikK(auxBasisController, nLocalAux, pairDomainToK, reducedOccIndices, kTokIndexMap, pair->i, klK).transpose());
  tmp_mkaqs.push_back(get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_i, k_PAOToFullPAOMaps,
                              reducedOccIndices, toPNO_i, kTokIndexMap, iaK));
  if (pair->i != pair->j) {
    tmp_mikq_T.push_back(
        get_ikK(auxBasisController, nLocalAux, pairDomainToK, reducedOccIndices, kTokIndexMap, pair->j, klK).transpose());
    tmp_mkaqs.push_back(get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_j, k_PAOToFullPAOMaps,
                                reducedOccIndices, toPNO_j, kTokIndexMap, iaK));
  }
  const Eigen::MatrixXd& m_ikq_T = tmp_mikq_T[0];
  const Eigen::MatrixXd& m_jkq_T = (pair->i != pair->j) ? tmp_mikq_T[1] : tmp_mikq_T[0];
  const std::vector<Eigen::MatrixXd>& m_kaqs_i = tmp_mkaqs[0];
  const std::vector<Eigen::MatrixXd>& m_kaqs_j = (pair->i != pair->j) ? tmp_mkaqs[1] : tmp_mkaqs[0];
  /* (ik|jl) */
  for (auto& klSet : pair->klPairSets) {
    unsigned int k = klSet->getKLPair()->i;
    unsigned int l = klSet->getKLPair()->j;
    const auto itK = kTokIndexMap.find(k);
    const auto itL = kTokIndexMap.find(l);
    if (itK == kTokIndexMap.end() || itL == kTokIndexMap.end())
      throw SerenityError("ERROR: Wrong index map!");
    const Eigen::VectorXd m_ikQ_V_T = llt_metric.solve(m_ikq_T.col(itK->second));
    const Eigen::VectorXd m_ilQ_V_T = llt_metric.solve(m_ikq_T.col(itL->second));
    const Eigen::VectorXd m_jkQ_T = m_jkq_T.col(itK->second);
    const Eigen::VectorXd m_jlQ_T = m_jkq_T.col(itL->second);
    klSet->ik_jl = (m_ikQ_V_T.transpose() * m_jlQ_T).eval()(0, 0);
    klSet->il_jk = (m_ilQ_V_T.transpose() * m_jkQ_T).eval()(0, 0);

    //(ki|la), (kj|la), (li|ka) and (lj|ka)
    const Eigen::MatrixXd& m_kaq_i = m_kaqs_i[itK->second];
    const Eigen::MatrixXd& m_kaq_j = m_kaqs_j[itK->second];
    const Eigen::MatrixXd& m_laq_j = m_kaqs_j[itL->second];
    const Eigen::MatrixXd& m_laq_i = m_kaqs_i[itL->second];
    klSet->ki_la = m_laq_j * m_ikQ_V_T;
    klSet->kj_la = m_laq_i * llt_metric.solve(m_jkQ_T);
    klSet->li_ka = m_kaq_j * m_ilQ_V_T;
    klSet->lj_ka = m_kaq_i * llt_metric.solve(m_jlQ_T);
  }
  // Use the symmetry, because we can!
  std::map<unsigned int, unsigned int> ijMap;
  ijMap.insert(std::make_pair(pair->i, 0));
  if (pair->i != pair->j)
    ijMap.insert(std::make_pair(pair->j, 1));
  std::vector<Eigen::MatrixXd> m_iaqs = get_iaK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                                k_PAOToFullPAOMaps, reducedOccIndices, toPNO, ijMap, iaK);
  for (auto& m_iaq : m_iaqs) {
    m_iaq = llt_metric.solve(m_iaq.transpose()).transpose().eval();
  }
  const Eigen::MatrixXd& m_iaq_V = m_iaqs[ijMap.find(pair->i)->second];
  const Eigen::MatrixXd& m_jaq_V = m_iaqs[ijMap.find(pair->j)->second];
  if (ijMap.find(pair->i) == ijMap.end() || ijMap.find(pair->j) == ijMap.end())
    throw SerenityError("ERROR: Wrong index map! 5");
  for (auto& kSet : pair->coupledPairs) {
    const auto itK = kTokIndexMap.find(kSet->getK());
    if (itK == kTokIndexMap.end())
      throw SerenityError("ERROR: Wrong index map! 2");
    /* (ja|ik) and (ia|jk) */
    kSet->ja_ik = m_jaq_V * m_ikq_T.col(itK->second);
    kSet->ia_jk = m_iaq_V * m_jkq_T.col(itK->second);
  } // for kSet
}

Eigen::VectorXd Ao2MoExchangeIntegralTransformer::calculate_invV_ijK(
    std::shared_ptr<OrbitalPair> pair, const Eigen::SparseVector<int>& pairDomainToK,
    std::shared_ptr<BasisController> auxBasisController, const Eigen::LLT<Eigen::MatrixXd>& llt_metric,
    const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
    const MO3CenterIntegrals& klK) {
  /*
   * TODO: The sorting of the klK integrals is not ideal.
   *       Change it if possible.
   *       Currently this is not the time limiting step (by far).
   *       Thus, this is something to keep in mind but currently not
   *       important.
   */
  const auto& auxBasis = auxBasisController->getBasis();
  Eigen::MatrixXd m_ijk = Eigen::VectorXd::Zero(llt_metric.cols());
  unsigned int kCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itK(pairDomainToK); itK; ++itK) {
    const auto& indices = reducedOccIndices[itK.row()];
    const auto it_i_red = indices->find(pair->i);
    const auto it_j_red = indices->find(pair->j);
    unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
    unsigned int nContK = auxBasis[itK.row()]->getNContracted();
    for (unsigned int kk = 0; kk < nContK; ++kk) {
      m_ijk(kCounter, 0) = klK[extendedK + kk](it_i_red->second, it_j_red->second);
      kCounter++;
    } // for kk
  }   // for itK
  return llt_metric.solve(m_ijk);
}

MO3CenterIntegrals Ao2MoExchangeIntegralTransformer::flipAndTransform(const MO3CenterIntegrals& ints,
                                                                      const Eigen::MatrixXd& toPNO,
                                                                      const Eigen::SparseMatrix<double>& sparseProjection) {
  unsigned int nPNOs = ints.size();
  unsigned int nTargetPNO = toPNO.cols();
  unsigned int nAux = ints[0].cols();
  MO3CenterIntegrals returnInts(nTargetPNO, Eigen::MatrixXd::Zero(nAux, nPNOs));
  for (unsigned int bPNO = 0; bPNO < nPNOs; ++bPNO) {
    // Transform
    const Eigen::MatrixXd m_bcK = (toPNO.transpose() * sparseProjection * ints[bPNO]).transpose();
    // Sort
    for (unsigned int cPNO = 0; cPNO < nTargetPNO; ++cPNO) {
      returnInts[cPNO].col(bPNO) = m_bcK.col(cPNO);
    }
  }
  // Flip
  for (auto& m_bcK : returnInts)
    m_bcK = m_bcK.transpose().eval();
  return returnInts;
}

MO3CenterIntegrals Ao2MoExchangeIntegralTransformer::flipAndTransform(const MO3CenterIntegrals& ints,
                                                                      const Eigen::MatrixXd& toPNO,
                                                                      const Eigen::SparseMatrix<double>& sparseProjection,
                                                                      const Eigen::SparseMatrix<double>& auxProjection) {
  unsigned int nPNOs = ints.size();
  unsigned int nTargetPNO = toPNO.cols();
  unsigned int nAux = auxProjection.rows();
  MO3CenterIntegrals returnInts(nTargetPNO, Eigen::MatrixXd::Zero(nAux, nPNOs));
  for (unsigned int bPNO = 0; bPNO < nPNOs; ++bPNO) {
    // Transform
    const Eigen::MatrixXd m_bcK = (toPNO.transpose() * sparseProjection * ints[bPNO] * auxProjection.transpose()).transpose();
    // Sort
    for (unsigned int cPNO = 0; cPNO < nTargetPNO; ++cPNO) {
      returnInts[cPNO].col(bPNO) = m_bcK.col(cPNO);
    }
  }
  // Flip
  for (auto& m_bcK : returnInts)
    m_bcK = m_bcK.transpose().eval();
  return returnInts;
}

void Ao2MoExchangeIntegralTransformer::calculate_2Virt_integrals(
    std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
    const Eigen::LLT<Eigen::MatrixXd>& llt_metric, const Eigen::SparseVector<int>& pairDomainToK,
    const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
    const MO3CenterIntegrals& m_abKs, const Eigen::SparseMatrix<double>& extendedPAODomainProjection,
    const Eigen::SparseMatrix<double>& extendedAuxProjectionMatrix_T, const Eigen::VectorXd& invV_ijK,
    const MO3CenterIntegrals& iaK, const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
    bool calcSigmaVectorInts, const MO3CenterIntegrals& abK) {
  bool oneThread = omp_get_num_threads() == 1;
  const unsigned int nLocalAux = llt_metric.cols();
  const Eigen::MatrixXd& toPNO = pair->toPAODomain;
  const Eigen::MatrixXd& toPNO_i = pair->singles_i->toPAODomain;
  const Eigen::MatrixXd& toPNO_j = pair->singles_j->toPAODomain;
  const unsigned int nPNOs = toPNO.cols();
  // A lot of sparse projections.
  // pair aux <--> all aux
  const Eigen::SparseMatrix<double> auxProjectionMatrix = buildSparseAuxProjection(auxBasisController, pairDomainToK);
  // pair aux <--> pair ext. aux.
  const Eigen::SparseMatrix<double> auxIJtoExtended = (auxProjectionMatrix * extendedAuxProjectionMatrix_T).pruned().eval();
  // pair PAO <--> pair ext. PAO
  const Eigen::SparseMatrix<double> domainProjection_ExTo_ij =
      pair->domainProjection * extendedPAODomainProjection.transpose();
  // sing. i PAO <--> all aux.
  const Eigen::SparseMatrix<double>& domainProjection_i = pair->singles_i->getDiagonalPair()->domainProjection;
  // sing. j PAO <--> all aux.
  const Eigen::SparseMatrix<double>& domainProjection_j = pair->singles_j->getDiagonalPair()->domainProjection;
  // sing. i PAO <--> pair ext. PAO
  const Eigen::SparseMatrix<double> domainProjection_ExTo_i = domainProjection_i * extendedPAODomainProjection.transpose();
  // sing. j PAO <--> pair ext. PAO
  const Eigen::SparseMatrix<double> domainProjection_ExTo_j = domainProjection_j * extendedPAODomainProjection.transpose();

  // Extract aux domain.
  std::vector<Eigen::MatrixXd> m_acks_KEx;
  for (const auto& m_abk : m_abKs) {
    m_acks_KEx.push_back((m_abk * auxIJtoExtended.transpose()).eval());
  }
  // Extract from extended domain integral list and transform the second index.
  std::vector<Eigen::MatrixXd> m_acks;
  for (const auto& m_abk : m_acks_KEx) {
    const Eigen::MatrixXd k_bd = toPNO.transpose() * domainProjection_ExTo_ij * m_abk;
    m_acks.push_back(k_bd);
  }
  // (ij|ab) integals a,b in [ij]
  pair->ij_ab = Eigen::MatrixXd::Zero(toPNO.cols(), toPNO.cols());
  for (unsigned int bPNO = 0; bPNO < toPNO.cols(); ++bPNO) {
    pair->ij_ab.col(bPNO) = m_acks[bPNO] * invV_ijK;
  } // for bPNO

  std::map<unsigned int, unsigned int> kTokIndexMap;
  unsigned int kIndex = 0;
  for (auto& kSet : pair->coupledPairs) {
    unsigned int k = kSet->getK();
    if (kTokIndexMap.find(k) == kTokIndexMap.end()) {
      kTokIndexMap.insert(std::make_pair(k, kIndex));
      ++kIndex;
    }
  } // for kSet
  const auto m_kaKs = get_iaK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection, k_PAOToFullPAOMaps,
                              reducedOccIndices, toPNO, kTokIndexMap, iaK);
  for (unsigned int ik = 0; ik < pair->coupledPairs.size(); ++ik) {
    std::shared_ptr<CouplingOrbitalSet> coupledPair = pair->coupledPairs[ik];
    unsigned int k = coupledPair->getK();
    const auto itK = kTokIndexMap.find(k);
    if (itK == kTokIndexMap.end())
      throw SerenityError("ERROR: Wrong index map. 6");
    const Eigen::MatrixXd& m_kaK = m_kaKs[itK->second];
    const Eigen::MatrixXd llt_m_kaK = llt_metric.solve(m_kaK.transpose());
    //(ka|bc) a,b,c in [ij]
    for (auto& m_bck : m_acks)
      coupledPair->ka_bc.push_back((m_bck * llt_m_kaK).transpose().eval());
    //(ij|ak) a in [ij]
    coupledPair->ij_ak = m_kaK * invV_ijK;
  } // for ik
  const auto itI = kTokIndexMap.find(pair->i);
  const auto itJ = kTokIndexMap.find(pair->j);
  if (itI == kTokIndexMap.end() || itJ == kTokIndexMap.end())
    throw SerenityError("ERROR: Wrong index map! 3");
  //(ja|bc) a,b in [ij] and c in [i]
  {
    const Eigen::MatrixXd& m_jaK = m_kaKs[itJ->second];
    const Eigen::MatrixXd llt_m_jaK = llt_metric.solve(m_jaK.transpose());
    const auto m_bcK_i = flipAndTransform(m_acks_KEx, toPNO_i, domainProjection_ExTo_i);
    pair->ja_bc = {};
    for (auto& m_bck : m_bcK_i)
      pair->ja_bc.push_back((m_bck * llt_m_jaK).transpose());
  }
  //(ia|bc) a,b in [ij] and c in [j]
  {
    const Eigen::MatrixXd& m_iaK = m_kaKs[itI->second];
    const Eigen::MatrixXd llt_m_iaK = llt_metric.solve(m_iaK.transpose());
    const auto m_bcK_j = flipAndTransform(m_acks_KEx, toPNO_j, domainProjection_ExTo_j);
    pair->ia_bc = {};
    for (auto& m_bck : m_bcK_j)
      pair->ia_bc.push_back((m_bck * llt_m_iaK).transpose());
  }

  if (calcSigmaVectorInts) {
    if (oneThread)
      Timings::takeTime("sigma ints. over ij");
    std::vector<Eigen::MatrixXd> tmp_m_iaK;
    tmp_m_iaK.push_back(get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_i, k_PAOToFullPAOMaps,
                                reducedOccIndices, toPNO_i, pair->i, iaK));
    if (pair->i != pair->j) {
      tmp_m_iaK.push_back(get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_j, k_PAOToFullPAOMaps,
                                  reducedOccIndices, toPNO_j, pair->j, iaK));
    }
    const Eigen::MatrixXd& m_iaK = tmp_m_iaK[0];
    const Eigen::MatrixXd& m_jbK = (pair->i != pair->j) ? tmp_m_iaK[1] : tmp_m_iaK[0];
    pair->iaS_jbSX2_M_ij_aSbS = 2 * (m_iaK * llt_metric.solve(m_jbK.transpose()).eval());
    // -(ij|ab)
    // This here is expensive ... However, I do not really see ho else to do this.
    Eigen::MatrixXd ij_ab = Eigen::MatrixXd::Zero(toPNO_i.cols(), toPNO_j.cols());
    for (unsigned int bPNO = 0; bPNO < toPNO_j.cols(); ++bPNO) {
      ij_ab.col(bPNO) = abK[bPNO] * invV_ijK;
    } // for bPNO
    pair->iaS_jbSX2_M_ij_aSbS -= ij_ab;
    if (oneThread)
      Timings::timeTaken("sigma ints. over ij");
  }

  //(jc|ab) c in [i], a,b in [ij]
  {
    const Eigen::MatrixXd m_jck_i = get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_i,
                                            k_PAOToFullPAOMaps, reducedOccIndices, toPNO_i, pair->j, iaK);
    const Eigen::MatrixXd llt_m_jck_i = llt_metric.solve(m_jck_i.transpose());
    for (unsigned int cPNO = 0; cPNO < toPNO_i.cols(); ++cPNO) {
      Eigen::MatrixXd ab_jc = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
      for (unsigned int bPNO = 0; bPNO < nPNOs; ++bPNO) {
        ab_jc.col(bPNO) = m_acks[bPNO] * llt_m_jck_i.col(cPNO);
      } // for bPNO
      assert((ab_jc - ab_jc.transpose()).array().abs().sum() < 1e-7);
      pair->jc_ab.push_back(ab_jc);
    } // for cPNO
  }
  //(ic|ab) c in [j], a,b in [ij]
  {
    const Eigen::MatrixXd m_ick_j = get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_j,
                                            k_PAOToFullPAOMaps, reducedOccIndices, toPNO_j, pair->i, iaK);
    const Eigen::MatrixXd llt_m_ick_j = llt_metric.solve(m_ick_j.transpose());
    for (unsigned int cPNO = 0; cPNO < toPNO_j.cols(); ++cPNO) {
      Eigen::MatrixXd ab_ic = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
      for (unsigned int bPNO = 0; bPNO < nPNOs; ++bPNO) {
        ab_ic.col(bPNO) = m_acks[bPNO] * llt_m_ick_j.col(cPNO);
      } // for bPNO
      assert((ab_ic - ab_ic.transpose()).array().abs().sum() < 1e-7);
      pair->ic_ab.push_back(ab_ic);
    } // for cPNO
  }

  //(ac|bd) integrals
  pair->ac_bd = std::make_unique<Matrix<Eigen::MatrixXd>>(nPNOs, nPNOs, Eigen::MatrixXd(0, 0));
  for (unsigned int bPNO = 0; bPNO < nPNOs; ++bPNO) {
    const Eigen::MatrixXd invV_k_bd = llt_metric.solve(m_acks[bPNO].transpose().eval());
    for (unsigned int aPNO = 0; aPNO <= bPNO; ++aPNO) {
      const Eigen::MatrixXd ac_bd = m_acks[aPNO] * invV_k_bd;
      (*pair->ac_bd)(aPNO, bPNO) = ac_bd;
    } // for bPNO
  }   // for aPNO
}

void Ao2MoExchangeIntegralTransformer::transformExchangeIntegrals(
    std::shared_ptr<BasisController> basisController, const Eigen::MatrixXd& aoCoefficients,
    std::shared_ptr<PAOController> paoController, std::vector<std::shared_ptr<OrbitalPair>>& orbitalPairs,
    std::shared_ptr<QuasiCanonicalPAODomainConstructor> pnoConstructor) {
  Timings::takeTime("Local Cor. -   Int. Transform.");
  Ao2MoHalfTransformer halfTransformer(basisController, basisController);
#ifdef _OPENMP
  Eigen::setNbThreads(1);
#endif
  for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
    auto& pair = orbitalPairs[iPair];
    pnoConstructor->transformExternalBasis(pair);
    Eigen::MatrixXd pairDensityMatrix = aoCoefficients.col(pair->i) * aoCoefficients.col(pair->j).transpose();
    Eigen::MatrixXd result;
    halfTransformer.transformTwoElectronIntegrals(result, pairDensityMatrix);
    Eigen::MatrixXd p_ij = paoController->getPAOsFromDomain(pair->paoDomain) * pair->toPAODomain;
    pair->k_ij = p_ij.transpose() * result * p_ij;
    pnoConstructor->postProcessing(pair);
  } // for pair
#ifdef _OPENMP
  Eigen::setNbThreads(0);
#endif
  OutputControl::nOut << " done" << std::endl;
  Timings::timeTaken("Local Cor. -   Int. Transform.");
}

void Ao2MoExchangeIntegralTransformer::transformExchangeIntegrals(
    std::shared_ptr<BasisController> auxBasisController, std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
    std::vector<std::shared_ptr<OrbitalPair>>& orbitalPairs,
    std::shared_ptr<QuasiCanonicalPAODomainConstructor> pnoConstructor) {
  Timings::takeTime("Local Cor. -   Int. Transform.");
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
#ifdef _OPENMP
  unsigned int nThreads = omp_get_max_threads();
  Eigen::setNbThreads(1);
#endif
  MatrixInBasis<scfMode> M(auxBasisController);
  takeTime("2-center Coulomb metric");
  calculateTwoCenterIntegrals(M);
  timeTaken(2, "2-center Coulomb metric");
  takeTime("Prepare Sparse Maps");
  const auto sparseMaps = mo3CenterIntegralController->getSparseMapsController();
  const SparseMap& occToK = sparseMaps->getOccToAuxShellMap();
  const auto& k_redToFullMapsPAO = mo3CenterIntegralController->getProjection(ORBITAL_TYPE::VIRTUAL);
  const auto& reducedOccIndices = mo3CenterIntegralController->getIndices(ORBITAL_TYPE::OCCUPIED);

  timeTaken(2, "Prepare Sparse Maps");

  const auto& auxBasis = auxBasisController->getBasis();
  // Temporary precalculation of all integrals.
  // This could be changed to a block wise calculation and writing.
  Eigen::SparseVector<int> auxSuperDomain = Eigen::VectorXi::Constant(auxBasis.size(), 1).sparseView();
  const MO3CenterIntegrals& iaK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ia_K, auxSuperDomain);
  OutputControl::nOut << "  Performing (ia|jb) integral transformation             ...";
  OutputControl::nOut.flush();
  takeTime("3-center int. normal.");
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
    auto& pair = orbitalPairs[iPair];
    pnoConstructor->transformExternalBasis(pair);
    if (pair->type == OrbitalPairTypes::VERY_DISTANT || pair->type == OrbitalPairTypes::DISTANT)
      continue;
    assert(pair->paoDomain.nonZeros() > 0);
    // Get the fitting domain indices
    unsigned int i = pair->i;
    unsigned int j = pair->j;
    /*
     * Build local metric for each pair by selecting a block from S_aux.
     */
    const Eigen::SparseVector<int> pairDomainToK = (occToK.col(i) + occToK.col(j)).pruned();
    const Eigen::MatrixXd localBlock = SystemSplittingTools<scfMode>::getMatrixBlockShellWise(M, pairDomainToK, pairDomainToK);

    unsigned int nLocalAux = localBlock.cols();
    unsigned int nLocalPAOs = pair->paoDomain.nonZeros();
    Eigen::MatrixXd m_i = Eigen::MatrixXd::Zero(nLocalPAOs, nLocalAux);
    Eigen::MatrixXd m_j = Eigen::MatrixXd::Zero(nLocalPAOs, nLocalAux);
    unsigned int kCounter = 0;
    for (Eigen::SparseVector<int>::InnerIterator itK(pairDomainToK); itK; ++itK) {
      unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
      unsigned int nContK = auxBasis[itK.row()]->getNContracted();
      const auto indices = reducedOccIndices[itK.row()];
      const auto it_i_red = indices->find(i);
      const auto it_j_red = indices->find(j);
      if (it_i_red == indices->end() || it_j_red == indices->end() || k_redToFullMapsPAO[itK.row()]->cols() < 1 ||
          iaK[extendedK].rows() < 1 || iaK[extendedK].cols() < 1) {
        kCounter += nContK;
        continue;
      }
      const Eigen::SparseMatrix<double> finalPAOProjection = (pair->domainProjection * *k_redToFullMapsPAO[itK.row()]).eval();
      for (unsigned int kk = 0; kk < nContK; ++kk) {
        m_i.col(kCounter).noalias() = finalPAOProjection * iaK[extendedK + kk].col(it_i_red->second);
        m_j.col(kCounter).noalias() = finalPAOProjection * iaK[extendedK + kk].col(it_j_red->second);
        kCounter++;
      } // for kk
    }   // for itK
    // The actual Cholesky decomposition and solving of the set of linear equations.
    pair->nAuxFunctions = m_i.cols();
    pair->k_ij.noalias() = pair->toPAODomain.transpose() * m_i *
                           localBlock.llt().solve((pair->toPAODomain.transpose() * m_j).eval().transpose());
    pnoConstructor->postProcessing(pair);
  } // for pair
#ifdef _OPENMP
  Eigen::setNbThreads(nThreads);
#endif
  OutputControl::nOut << " done" << std::endl;
  timeTaken(1, "3-center int. normal.");
  Timings::timeTaken("Local Cor. -   Int. Transform.");
}

void Ao2MoExchangeIntegralTransformer::transformAllIntegrals(std::shared_ptr<BasisController> auxBasisController,
                                                             std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                                                             std::vector<std::shared_ptr<OrbitalPair>>& orbitalPairs,
                                                             bool calcSigmaVectorInts, bool lowMemory,
                                                             double memoryDemandForSet, bool ignoreMemoryHandling) {
  Timings::takeTime("Local Cor. -   Int. Transform.");
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  const auto sparseMaps = mo3CenterIntegralController->getSparseMapsController();
  const SparseMap& occToK = sparseMaps->getOccToAuxShellMap();
  const SparseMap& extendedOccToPAO = sparseMaps->getExtendedOccToPAOMap();
  // Temporary precalculation of all integrals.
  // This could be changed to a block wise calculation and writing.
  Eigen::SparseVector<int> auxSuperDomain(auxBasisController->getBasis().size());
  for (auto& pair : orbitalPairs) {
    auxSuperDomain += pair->getFittingDomain();
  }
  Eigen::SparseVector<int> extendedPAODomain(0);
  if (lowMemory) {
    auto firstPair = orbitalPairs[0];
    extendedPAODomain = (extendedOccToPAO.col(firstPair->i) + extendedOccToPAO.col(firstPair->j)).eval();
    for (unsigned int i = 1; i < orbitalPairs.size(); ++i) {
      auto currentPair = orbitalPairs[i];
      extendedPAODomain += (extendedOccToPAO.col(currentPair->i) + extendedOccToPAO.col(currentPair->j)).eval();
    }
  }

  const MO3CenterIntegrals& iaK =
      mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ia_K, auxSuperDomain, extendedPAODomain);
  const MO3CenterIntegrals& abK =
      mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ab_K, auxSuperDomain, extendedPAODomain);
  const MO3CenterIntegrals& klK =
      mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::kl_K, auxSuperDomain, extendedPAODomain);

  auto& libint = Libint::getInstance();
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);

  const auto& k_redToFullMapsPAO = mo3CenterIntegralController->getProjection(ORBITAL_TYPE::VIRTUAL);
  const auto& reducedOccIndices = mo3CenterIntegralController->getIndices(ORBITAL_TYPE::OCCUPIED);

  unsigned int nThreads = 1;
  unsigned int nMaxThreads = 1;
#ifdef _OPENMP
  /*
   * Reserve 1.5 GB memory for each thread and make sure that the number of threads is at least 1.
   */
  nMaxThreads = omp_get_max_threads();
  if (ignoreMemoryHandling) {
    OutputControl::nOut
        << "  Note: Assuming that sufficient memory is available and ignoring any constraints by the operating system."
        << std::endl;
    nThreads = nMaxThreads;
  }
  else {
    double availableMemory = MemoryManager::getInstance()->getAvailableSystemMemory() - memoryDemandForSet;
    availableMemory = (availableMemory < 0) ? 0 : availableMemory;
    double memoryPerThread = 1.5e+9;
    nThreads = std::min((unsigned int)(availableMemory / memoryPerThread), nMaxThreads);
    nThreads = (nThreads == 0) ? 1 : nThreads;
    omp_set_num_threads(nThreads);
  }
  Eigen::setNbThreads(1);
  OutputControl::vOut << "Number of threads used in the integral transformation: " << nThreads << std::endl;
#endif
  MatrixInBasis<scfMode> M(auxBasisController);
  calculateTwoCenterIntegrals(M);
  OutputControl::nOut << "  Entering final Integral Transformation                 ...";
  OutputControl::nOut.flush();
  takeTime("3-center int. normal.");
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
    auto& pair = orbitalPairs[iPair];
    if (pair->type == OrbitalPairTypes::VERY_DISTANT || pair->type == OrbitalPairTypes::DISTANT)
      assert(false);
    assert(pair->k_ij.cols() > 0 && "(ia|jb) integrals are not initialized!");
    // Get the fitting domain indices
    unsigned int i = pair->i;
    unsigned int j = pair->j;

    /*
     * Build local metric for each pair by selecting a block from S_aux.
     */
    const Eigen::SparseVector<int> pairDomainToK = (occToK.col(i) + occToK.col(j)).pruned();
    const Eigen::SparseVector<int>& extendedAuxDomain = pair->getFittingDomain();
    const Eigen::SparseMatrix<double> extendedAuxProjectionMatrix_T =
        buildSparseAuxProjection(auxBasisController, extendedAuxDomain).transpose();

    const Eigen::SparseVector<int> extendedPairPAODomain = (extendedOccToPAO.col(i) + extendedOccToPAO.col(j)).pruned();
    const auto extendedPAODomainProjection_T = constructProjectionMatrixFromSparse_T(extendedPairPAODomain);
    const Eigen::MatrixXd localBlock = SystemSplittingTools<scfMode>::getMatrixBlockShellWise(M, pairDomainToK, pairDomainToK);
    // Calculate the Cholesky decomposition of the local metric.
    const Eigen::LLT<Eigen::MatrixXd> llt = localBlock.llt();
    // Calculate V^(-1) * (K|ij) <-- needed in at least two transformations.
    const Eigen::VectorXd invV_ijK = calculate_invV_ijK(pair, pairDomainToK, auxBasisController, llt, reducedOccIndices, klK);

    // This extraction step is the most expensive part. Thus, we will start with it and reuse the result.
    if (nThreads == 1)
      Timings::takeTime("(ac|K) integrals");
    // Extract (ac|K) integrals, a in [ij], c in ext. PAO [ij]
    MO3CenterIntegrals mEx_amuQ, cdQ;
    const unsigned int nExtendedLocalAux = extendedAuxProjectionMatrix_T.cols();
    if (calcSigmaVectorInts) {
      const auto singles_i = pair->singles_i;
      const auto singles_j = pair->singles_j;
      const auto domainProjection_i = singles_i->getDiagonalPair()->domainProjection;
      const auto domainProjection_j = singles_j->getDiagonalPair()->domainProjection;
      mEx_amuQ = get_abK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, pair->domainProjection,
                         extendedPAODomainProjection_T, domainProjection_j, domainProjection_i, k_redToFullMapsPAO,
                         pair->toPAODomain, singles_j->toPAODomain, singles_i->toPAODomain, llt.cols(), pairDomainToK,
                         abK, cdQ);
    }
    else {
      mEx_amuQ = get_abK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, pair->domainProjection,
                         extendedPAODomainProjection_T, k_redToFullMapsPAO, pair->toPAODomain, abK);
    }
    if (nThreads == 1)
      Timings::timeTaken("(ac|K) integrals");

    if (nThreads == 1)
      Timings::takeTime("2-4 virt. integrals");
    // Calculate all integrals with at least two virtual indices over non-mixed PNO domains.
    calculate_2Virt_integrals(pair, auxBasisController, llt, pairDomainToK, k_redToFullMapsPAO, mEx_amuQ,
                              extendedPAODomainProjection_T, extendedAuxProjectionMatrix_T, invV_ijK, iaK,
                              reducedOccIndices, calcSigmaVectorInts, cdQ);
    if (nThreads == 1)
      Timings::timeTaken("2-4 virt. integrals");

    if (nThreads == 1)
      Timings::takeTime("3,4 occ. integrals");
    // Calculate all integrals with at least 3 occupied indices.
    calculate_3Occ_Integrals(pair, auxBasisController, llt, pairDomainToK, k_redToFullMapsPAO, reducedOccIndices, iaK, klK);
    if (nThreads == 1)
      Timings::timeTaken("3,4 occ. integrals");

    if (nThreads == 1)
      Timings::takeTime("Mixed PNO-basis integrals");
    // Calculate integrals over mixed PNO domains.
    calculateMixedIntegrals(pair, auxBasisController, M, k_redToFullMapsPAO, reducedOccIndices, occToK, iaK, klK,
                            extendedPAODomainProjection_T, extendedAuxDomain, extendedAuxProjectionMatrix_T, mEx_amuQ,
                            calcSigmaVectorInts);
    if (nThreads == 1)
      Timings::timeTaken("Mixed PNO-basis integrals");
  } // for pair
#ifdef _OPENMP
  Eigen::setNbThreads(nMaxThreads);
  omp_set_num_threads(nMaxThreads);
#endif
  OutputControl::nOut << " done" << std::endl;
  timeTaken(1, "3-center int. normal.");
  Timings::timeTaken("Local Cor. -   Int. Transform.");
}

void Ao2MoExchangeIntegralTransformer::transformAllIntegrals(std::shared_ptr<BasisController> auxBasisController,
                                                             std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                                                             std::vector<std::shared_ptr<OrbitalPairSet>> orbitalPairSets,
                                                             bool dumpIntegrals, std::string pairIntegralFileName,
                                                             bool calcSigmaVectorInts, bool lowMemory,
                                                             bool ignoreMemoryHandling) {
  std::shared_ptr<HDF5::H5File> file;
  if (dumpIntegrals) {
    file = std::make_shared<HDF5::H5File>(pairIntegralFileName.c_str(), H5F_ACC_TRUNC);
  }
  for (unsigned int iSet = 0; iSet < orbitalPairSets.size(); ++iSet) {
    auto& orbitalPairSet = orbitalPairSets[iSet];
    // Transform integrals
    transformAllIntegrals(auxBasisController, mo3CenterIntegralController, *orbitalPairSet, calcSigmaVectorInts,
                          lowMemory, orbitalPairSet->memoryDemand(calcSigmaVectorInts), ignoreMemoryHandling);
    orbitalPairSet->setInMemory(true);
    // Write to file
    if (dumpIntegrals) {
      takeTime("Writing Integral Sets to File");
      orbitalPairSet->toHDF5(*file);
      file->flush(H5F_SCOPE_GLOBAL);
      if (iSet != orbitalPairSets.size() - 1)
        orbitalPairSet->removeInteralsFromMemory();
      timeTaken(2, "Writing Integral Sets to File");
    }
  } // for orbitalPairSets
  if (file)
    file->close();
}

} /* namespace Serenity */
