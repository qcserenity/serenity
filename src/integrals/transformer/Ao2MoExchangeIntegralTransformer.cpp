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
#include "basis/AtomCenteredBasisController.h"                        //AtomCenteredBasisController
#include "data/OrbitalPair.h"                                         //OrbitalPair definition
#include "data/PAOController.h"                                       //PAO controller.
#include "data/SingleSubstitution.h"                                  //Definition of a SingleSubstitution
#include "integrals/RI_J_IntegralController.h" //RI_J_IntegralController for schwartz based prescreening transformation.
#include "integrals/transformer/Ao2MoHalfTransformer.h" //Transformation via full four center integrals
#include "integrals/wrappers/Libint.h"                  //Integrals
#include "io/FormattedOutputStream.h"                   //Filtered output.
#include "math/linearAlgebra/MatrixFunctions.h"         //Symmetric orthogonalization of the aux. basis.
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
  auto libint = Libint::getSharedPtr();
  libint->initialize(libint2::Operator::coulomb, 0, 2);
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
      if (libint->compute(libint2::Operator::coulomb, 0, shellA, shellB, ints)) {
        Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(), nJ, nI);
        metric.block(qStart, pStart, nJ, nI) = tmp;
        metric.block(pStart, qStart, nI, nJ) = tmp.transpose();
      }
    }
  }
  libint->finalize(libint2::Operator::coulomb, 0, 2);
  timeTaken(3, "Calc 2-center ints");
}

std::pair<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>> Ao2MoExchangeIntegralTransformer::get_abK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux,
    const Eigen::SparseVector<int>& auxDomain, const Eigen::SparseMatrix<double>& projectionMatrix_a,
    const Eigen::SparseMatrix<double>& projectionMatrix_b, const Eigen::SparseMatrix<double>& projectionMatrix_c,
    const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps, const Eigen::MatrixXd& toPNO_a,
    const Eigen::MatrixXd& toPNO_b, const Eigen::MatrixXd& toPNO_c, const MO3CenterIntegrals& abK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nPNOs_a = toPNO_a.cols();
  const unsigned int nPNOs_b = toPNO_b.cols();
  const unsigned int nPNOs_c = toPNO_c.cols();
  const Eigen::MatrixXd toPNO_a_T = toPNO_a.transpose();
  const Eigen::MatrixXd toPNO_b_T = toPNO_b.transpose();
  const Eigen::MatrixXd toPNO_c_T = toPNO_c.transpose();
  std::pair<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>> result =
      std::make_pair(std::vector<Eigen::MatrixXd>(nPNOs_a, Eigen::MatrixXd::Zero(nPNOs_b, nLocalAux)),
                     std::vector<Eigen::MatrixXd>(nPNOs_a, Eigen::MatrixXd::Zero(nPNOs_c, nLocalAux)));
  std::vector<Eigen::MatrixXd>& m_abks = result.first;
  std::vector<Eigen::MatrixXd>& m_acks = result.second;
  Eigen::MatrixXd m_abk_tmp(toPNO_b_T.rows(), toPNO_a.cols());
  Eigen::MatrixXd m_ack_tmp(toPNO_c_T.rows(), toPNO_a.cols());
  unsigned int kCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itK(auxDomain); itK; ++itK) {
    unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
    unsigned int nContK = auxBasis[itK.row()]->getNContracted();
    if (k_PAOToFullPAOMaps[itK.row()].cols() < 1 || abK[extendedK].rows() < 1) {
      kCounter += nContK;
      continue;
    }
    // Transform both indices.
    const Eigen::SparseMatrix<double> finalProjection_b = (projectionMatrix_b * k_PAOToFullPAOMaps[itK.row()]).pruned();
    const Eigen::SparseMatrix<double> finalProjection_c = (projectionMatrix_c * k_PAOToFullPAOMaps[itK.row()]).pruned();
    const Eigen::SparseMatrix<double> finalProjection_a = (projectionMatrix_a * k_PAOToFullPAOMaps[itK.row()]).pruned();
    for (unsigned int kk = 0; kk < nContK; ++kk) {
      const Eigen::MatrixXd tmp_b = (finalProjection_b * abK[extendedK + kk] * finalProjection_a.transpose());
      const Eigen::MatrixXd tmp_c = (finalProjection_c * abK[extendedK + kk] * finalProjection_a.transpose());
      m_abk_tmp = toPNO_b_T * tmp_b * toPNO_a;
      m_ack_tmp = toPNO_c_T * tmp_c * toPNO_a;
      // Resort the integrals.
      for (unsigned int aPNO = 0; aPNO < nPNOs_a; ++aPNO) {
        m_abks[aPNO].col(kCounter) = m_abk_tmp.col(aPNO);
        m_acks[aPNO].col(kCounter) = m_ack_tmp.col(aPNO);
      } // for aPNO
      kCounter++;
    } // for kk
  }   // for itK
  return result;
}

std::vector<Eigen::MatrixXd> Ao2MoExchangeIntegralTransformer::get_abK(
    std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
    const Eigen::SparseMatrix<double>& projectionMatrix_a, const Eigen::SparseMatrix<double>& projectionMatrix_b,
    const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps, const Eigen::MatrixXd& toPNO_a,
    const Eigen::MatrixXd& toPNO_b, const MO3CenterIntegrals& abK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nPNOs_a = toPNO_a.cols();
  const unsigned int nPNOs_b = toPNO_b.cols();
  const Eigen::MatrixXd toPNO_a_T = toPNO_a.transpose();
  const Eigen::MatrixXd toPNO_b_T = toPNO_b.transpose();
  std::vector<Eigen::MatrixXd> m_abks(nPNOs_a, Eigen::MatrixXd::Zero(nPNOs_b, nLocalAux));
  Eigen::MatrixXd m_abk_tmp(toPNO_b_T.rows(), toPNO_a.cols());
  unsigned int kCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itK(auxDomain); itK; ++itK) {
    unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
    unsigned int nContK = auxBasis[itK.row()]->getNContracted();
    if (k_PAOToFullPAOMaps[itK.row()].cols() < 1 || abK[extendedK].rows() < 1) {
      kCounter += nContK;
      continue;
    }
    // Transform both indices.
    const Eigen::SparseMatrix<double> finalProjection_b = (projectionMatrix_b * k_PAOToFullPAOMaps[itK.row()]).pruned();
    const Eigen::SparseMatrix<double> finalProjection_a = (projectionMatrix_a * k_PAOToFullPAOMaps[itK.row()]).pruned();
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

Eigen::MatrixXd Ao2MoExchangeIntegralTransformer::get_iaK(std::shared_ptr<BasisController> auxBasisController,
                                                          unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
                                                          const Eigen::SparseMatrix<double>& projectionMatrix_a,
                                                          const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
                                                          const std::vector<Eigen::VectorXi>& reducedOccIndices,
                                                          const Eigen::MatrixXd& toPNO_a, const unsigned int i,
                                                          const MO3CenterIntegrals& iaK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nLocalPAOs_a = toPNO_a.rows();
  Eigen::MatrixXd m_iaK = Eigen::MatrixXd::Zero(nLocalPAOs_a, nLocalAux);
  unsigned int qCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itQ(auxDomain); itQ; ++itQ) {
    unsigned int extendedQ = auxBasisController->extendedIndex(itQ.row());
    unsigned int nContQ = auxBasis[itQ.row()]->getNContracted();
    const int i_red = reducedOccIndices[itQ.row()](i);
    if (i_red < 0 || k_PAOToFullPAOMaps[itQ.row()].cols() < 1 || iaK[extendedQ].rows() < 1) {
      qCounter += nContQ;
      continue;
    }
    const Eigen::SparseMatrix<double> finalPAOProjection = (projectionMatrix_a * k_PAOToFullPAOMaps[itQ.row()]).eval();
    for (unsigned int qq = 0; qq < nContQ; ++qq) {
      unsigned int totQIndex = extendedQ + qq;
      m_iaK.col(qCounter) = finalPAOProjection * iaK[totQIndex].col(i_red);
      qCounter++;
    } // for qq
  }   // for itQ
  m_iaK = toPNO_a.transpose() * m_iaK;
  return m_iaK;
}

Eigen::MatrixXd Ao2MoExchangeIntegralTransformer::get_ikK(std::shared_ptr<BasisController> auxBasisController,
                                                          unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
                                                          const std::vector<Eigen::VectorXi>& reducedOccIndices,
                                                          std::vector<unsigned int> kIndices, const unsigned int i,
                                                          const MO3CenterIntegrals& klK) {
  const auto& auxBasis = auxBasisController->getBasis();
  const unsigned int nK = kIndices.size();
  Eigen::MatrixXd m_ikK = Eigen::MatrixXd::Zero(nK, nLocalAux);
  auto m_ikKptr = m_ikK.data();
  unsigned int qCounter = 0;
  for (Eigen::SparseVector<int>::InnerIterator itQ(auxDomain); itQ; ++itQ) {
    const int i_red = reducedOccIndices[itQ.row()](i);
    unsigned int nContQ = auxBasis[itQ.row()]->getNContracted();
    if (i_red < 0) {
      qCounter += nContQ;
      continue;
    }
    unsigned int extendedQ = auxBasisController->extendedIndex(itQ.row());
    for (unsigned int qq = 0; qq < nContQ; ++qq) {
      unsigned int totQIndex = extendedQ + qq;
      const Eigen::VectorXd tmp_i = Eigen::VectorXd(klK[totQIndex].col(i_red));
      unsigned int kCounter = 0;
      for (auto k : kIndices) {
        const int k_red = reducedOccIndices[itQ.row()](k);
        if (k_red < 0) {
          ++kCounter;
          continue;
        }
        m_ikKptr[qCounter * nK + kCounter] = tmp_i(k_red);
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
    const MatrixInBasis<RESTRICTED>& metric, const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
    const std::vector<Eigen::VectorXi>& reducedOccIndices, const SparseMap& occToK, const MO3CenterIntegrals& iaK,
    const MO3CenterIntegrals& abK, const MO3CenterIntegrals& klK) {
  const Eigen::MatrixXd& toPNO_ij = pair->toPAODomain;
  bool takeTime = false;
#ifdef _OPENMP
  takeTime = omp_get_max_threads() == 1;
#endif

  /*
   * Extract all integrals ikQ/jkQ, iaQ/jaQ which will be used at some point during
   * the k-Set loop. The integrals used for each k-Set can easily be extracted by
   * sparse projection. The projection should only lead to an insignificant overhead,
   * while the early integral extraction is slightly faster than extracting in the loop
   * itself.
   */
  Eigen::SparseVector<int> extendedAuxDomain = (occToK.col(pair->i) + occToK.col(pair->j)).pruned();
  for (auto& kSet : pair->coupledPairs)
    extendedAuxDomain += occToK.col(kSet->getK());
  extendedAuxDomain.pruned().eval();
  const Eigen::SparseMatrix<double> extendedAuxProjectionMatrix_T =
      buildSparseAuxProjection(auxBasisController, extendedAuxDomain).transpose();
  const Eigen::MatrixXd extendedBlock =
      SystemSplittingTools<RESTRICTED>::getMatrixBlockShellWise(metric, extendedAuxDomain, extendedAuxDomain);

  std::vector<unsigned int> kIndices;
  for (auto& kSet : pair->coupledPairs)
    kIndices.push_back(kSet->getK());

  const unsigned int nExtendedLocalAux = extendedBlock.cols();
  Eigen::MatrixXd mEx_iaQ = get_iaK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, pair->domainProjection,
                                    k_PAOToFullPAOMaps, reducedOccIndices, toPNO_ij, pair->i, iaK);
  Eigen::MatrixXd mEx_jaQ = get_iaK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, pair->domainProjection,
                                    k_PAOToFullPAOMaps, reducedOccIndices, toPNO_ij, pair->j, iaK);
  Eigen::MatrixXd mEx_ikQ_T =
      get_ikK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, reducedOccIndices, kIndices, pair->i, klK)
          .transpose()
          .eval();
  Eigen::MatrixXd mEx_jkQ_T =
      get_ikK(auxBasisController, nExtendedLocalAux, extendedAuxDomain, reducedOccIndices, kIndices, pair->j, klK)
          .transpose()
          .eval();

  for (unsigned int ik = 0; ik < pair->coupledPairs.size(); ++ik) {
    auto kSet = pair->coupledPairs[ik];
    const unsigned int k = kSet->getK();
    const Eigen::SparseVector<int> auxDomain = (occToK.col(pair->i) + occToK.col(pair->j) + occToK.col(k)).pruned().eval();

    const Eigen::SparseMatrix<double> auxProjectionMatrix = buildSparseAuxProjection(auxBasisController, auxDomain);
    const Eigen::SparseMatrix<double> auxIJKtoExtended =
        (auxProjectionMatrix * extendedAuxProjectionMatrix_T).pruned().eval();

    const Eigen::MatrixXd localBlock = (auxIJKtoExtended * extendedBlock * auxIJKtoExtended.transpose()).eval();
    const unsigned int nLocalAux = localBlock.cols();
    const Eigen::LLT<Eigen::MatrixXd> llt = localBlock.llt();
    Eigen::MatrixXd m_iaQ = mEx_iaQ * auxIJKtoExtended.transpose();
    Eigen::MatrixXd m_jaQ = mEx_jaQ * auxIJKtoExtended.transpose();
    Eigen::MatrixXd m_ikQ = auxIJKtoExtended * mEx_ikQ_T.col(ik);
    Eigen::MatrixXd m_jkQ = auxIJKtoExtended * mEx_jkQ_T.col(ik);
    m_ikQ = llt.solve(m_ikQ).eval();
    m_jkQ = llt.solve(m_jkQ).eval();
    m_iaQ = llt.solve(m_iaQ.transpose()).eval();
    m_jaQ = llt.solve(m_jaQ.transpose()).eval();

    const Eigen::MatrixXd& toPNO_ik = kSet->getIKPair()->toPAODomain;
    const Eigen::MatrixXd& toPNO_kj = kSet->getKJPair()->toPAODomain;
    const Eigen::SparseMatrix<double>& domainProjection_ik = kSet->getIKPair()->domainProjection;
    const Eigen::SparseMatrix<double>& domainProjection_kj = kSet->getKJPair()->domainProjection;

    if (takeTime)
      Timings::takeTime("Mixed. K-Set extraction");
    // Coulomb type integrals.
    // This extraction steps cost 80% of the computational time ...
    const auto mVirtVirtAuxPair = get_abK(auxBasisController, nLocalAux, auxDomain, pair->domainProjection, domainProjection_ik,
                                          domainProjection_kj, k_PAOToFullPAOMaps, toPNO_ij, toPNO_ik, toPNO_kj, abK);
    const std::vector<Eigen::MatrixXd>& m_acQ_ik = mVirtVirtAuxPair.first;
    const std::vector<Eigen::MatrixXd>& m_acQ_kj = mVirtVirtAuxPair.second;
    if (takeTime)
      Timings::timeTaken("Mixed. K-Set extraction");

    kSet->ik_ca = Eigen::MatrixXd::Zero(toPNO_kj.cols(), toPNO_ij.cols());
    kSet->jk_ca = Eigen::MatrixXd::Zero(toPNO_ik.cols(), toPNO_ij.cols());
    for (unsigned int aPNO = 0; aPNO < m_acQ_ik.size(); ++aPNO) {
      kSet->ik_ca.col(aPNO) = m_acQ_kj[aPNO] * m_ikQ;
      kSet->jk_ca.col(aPNO) = m_acQ_ik[aPNO] * m_jkQ;
    } // for aPNO

    // Exchange type integrals
    const Eigen::MatrixXd m_kcQ_ik = get_iaK(auxBasisController, nLocalAux, auxDomain, domainProjection_ik,
                                             k_PAOToFullPAOMaps, reducedOccIndices, toPNO_ik, k, iaK);
    const Eigen::MatrixXd m_kcQ_kj = get_iaK(auxBasisController, nLocalAux, auxDomain, domainProjection_kj,
                                             k_PAOToFullPAOMaps, reducedOccIndices, toPNO_kj, k, iaK);
    kSet->ia_kc = (m_kcQ_kj * m_iaQ).transpose().eval();
    kSet->ja_kc = (m_kcQ_ik * m_jaQ).transpose().eval();
  } // for kSet
}

void Ao2MoExchangeIntegralTransformer::calculate_ijka_integrals(
    std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
    const Eigen::SparseVector<int>& pairDomainToK, const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
    const std::vector<Eigen::VectorXi>& reducedOccIndices, const MO3CenterIntegrals& iaK, const Eigen::VectorXd& invV_ijK) {
  const unsigned int nLocalAux = invV_ijK.rows();
  const Eigen::MatrixXd& toPNO = pair->toPAODomain;

  for (auto coupledPair : pair->coupledPairs) {
    unsigned int k = coupledPair->getK();
    Eigen::MatrixXd m_kaQ = get_iaK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                    k_PAOToFullPAOMaps, reducedOccIndices, toPNO, k, iaK);
    coupledPair->ij_ak = m_kaQ * invV_ijK;
  } // for k
}

void Ao2MoExchangeIntegralTransformer::calculate_ikjaANDikjlANDkilaANDkjla_integrals(
    std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
    const Eigen::LLT<Eigen::MatrixXd>& llt_metric, const Eigen::SparseVector<int>& pairDomainToK,
    const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
    const std::vector<Eigen::VectorXi>& reducedOccIndices, const MO3CenterIntegrals& iaK, const MO3CenterIntegrals& klK) {
  const unsigned int nLocalAux = llt_metric.cols();
  Eigen::SparseMatrix<double>& domainProjection_i = pair->singles_i->getDiagonalPair()->domainProjection;
  Eigen::SparseMatrix<double>& domainProjection_j = pair->singles_j->getDiagonalPair()->domainProjection;
  const Eigen::MatrixXd& toPNO = pair->toPAODomain;
  const Eigen::MatrixXd& toPNO_i = pair->singles_i->toPAODomain;
  const Eigen::MatrixXd& toPNO_j = pair->singles_j->toPAODomain;

  std::vector<unsigned int> kIndices;
  for (const auto& kSet : pair->coupledPairs)
    kIndices.push_back(kSet->getK());

  const Eigen::MatrixXd m_ikq =
      get_ikK(auxBasisController, nLocalAux, pairDomainToK, reducedOccIndices, kIndices, pair->i, klK).transpose();
  const Eigen::MatrixXd m_jkq =
      get_ikK(auxBasisController, nLocalAux, pairDomainToK, reducedOccIndices, kIndices, pair->j, klK).transpose();
  /* (ik|jl) */
  for (auto& klSet : pair->klPairSets) {
    unsigned int k = klSet->getKLPair()->i;
    unsigned int l = klSet->getKLPair()->j;
    const Eigen::MatrixXd m_ikQ =
        get_ikK(auxBasisController, nLocalAux, pairDomainToK, reducedOccIndices, {k, l}, pair->i, klK);
    const Eigen::MatrixXd m_jlQ =
        get_ikK(auxBasisController, nLocalAux, pairDomainToK, reducedOccIndices, {k, l}, pair->j, klK);
    const Eigen::MatrixXd m_ikQ_V = llt_metric.solve(m_ikQ.transpose());
    const Eigen::MatrixXd m_jlQ_V = llt_metric.solve(m_jlQ.transpose());
    klSet->ik_jl = (m_ikQ_V.col(0).transpose() * m_jlQ.row(1).transpose()).eval()(0, 0);
    klSet->il_jk = (m_ikQ_V.col(1).transpose() * m_jlQ.row(0).transpose()).eval()(0, 0);

    //(ki|la), (kj|la), (li|ka) and (lj|ka)
    const Eigen::MatrixXd m_kaq_i = get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_i,
                                            k_PAOToFullPAOMaps, reducedOccIndices, toPNO_i, k, iaK);
    const Eigen::MatrixXd m_kaq_j = get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_j,
                                            k_PAOToFullPAOMaps, reducedOccIndices, toPNO_j, k, iaK);
    const Eigen::MatrixXd m_laq_j = get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_j,
                                            k_PAOToFullPAOMaps, reducedOccIndices, toPNO_j, l, iaK);
    const Eigen::MatrixXd m_laq_i = get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_i,
                                            k_PAOToFullPAOMaps, reducedOccIndices, toPNO_i, l, iaK);
    klSet->ki_la = m_laq_j * m_ikQ_V.col(0);
    klSet->kj_la = m_laq_i * m_jlQ_V.col(0);
    klSet->li_ka = m_kaq_j * m_ikQ_V.col(1);
    klSet->lj_ka = m_kaq_i * m_jlQ_V.col(1);
  }

  Eigen::MatrixXd m_iaq = get_iaK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                  k_PAOToFullPAOMaps, reducedOccIndices, toPNO, pair->i, iaK);
  Eigen::MatrixXd m_jaq = get_iaK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                  k_PAOToFullPAOMaps, reducedOccIndices, toPNO, pair->j, iaK);
  m_iaq = llt_metric.solve(m_iaq.transpose()).transpose().eval();
  m_jaq = llt_metric.solve(m_jaq.transpose()).transpose().eval();

  unsigned int iK = 0;
  for (auto& kSet : pair->coupledPairs) {
    /* (ja|ik) and (ia|jk) */
    kSet->ja_ik = m_jaq * m_ikq.col(iK);
    kSet->ia_jk = m_iaq * m_jkq.col(iK);
    ++iK;
  } // for kSet
}

Eigen::VectorXd Ao2MoExchangeIntegralTransformer::calculate_invV_ijK(std::shared_ptr<OrbitalPair> pair,
                                                                     const Eigen::SparseVector<int>& pairDomainToK,
                                                                     std::shared_ptr<BasisController> auxBasisController,
                                                                     const Eigen::LLT<Eigen::MatrixXd>& llt_metric,
                                                                     const std::vector<Eigen::VectorXi>& reducedOccIndices,
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
    const unsigned int i_red = reducedOccIndices[itK.row()](pair->i);
    const unsigned int j_red = reducedOccIndices[itK.row()](pair->j);
    unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
    unsigned int nContK = auxBasis[itK.row()]->getNContracted();
    for (unsigned int kk = 0; kk < nContK; ++kk) {
      m_ijk(kCounter, 0) = klK[extendedK + kk](i_red, j_red);
      kCounter++;
    } // for kk
  }   // for itK
  return llt_metric.solve(m_ijk);
}

void Ao2MoExchangeIntegralTransformer::calculate_ijab_integrals(
    std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
    const Eigen::SparseVector<int>& pairDomainToK, const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
    const MO3CenterIntegrals& abK, const Eigen::VectorXd& invV_ijK) {
  const unsigned int nLocalAux = invV_ijK.rows();
  const Eigen::MatrixXd& toPNO = pair->toPAODomain;
  const std::vector<Eigen::MatrixXd> m_abK_ij = get_abK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                                        pair->domainProjection, k_PAOToFullPAOMaps, toPNO, toPNO, abK);
  pair->ij_ab = Eigen::MatrixXd::Zero(toPNO.cols(), toPNO.cols());
  for (unsigned int bPNO = 0; bPNO < toPNO.cols(); ++bPNO) {
    pair->ij_ab.col(bPNO) = m_abK_ij[bPNO] * invV_ijK;
  } // for bPNO
}

void Ao2MoExchangeIntegralTransformer::calculate_kabcANDib_acANDjb_ac_integrals(
    std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
    const Eigen::LLT<Eigen::MatrixXd>& llt_metric, const Eigen::SparseVector<int>& pairDomainToK,
    const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
    const std::vector<Eigen::VectorXi>& reducedOccIndices, const MO3CenterIntegrals& abK, const MO3CenterIntegrals& iaK) {
  const auto& domainProjection_i = pair->singles_i->getDiagonalPair()->domainProjection;
  const auto& domainProjection_j = pair->singles_j->getDiagonalPair()->domainProjection;
  const unsigned int nLocalAux = llt_metric.cols();
  const Eigen::MatrixXd& toPNO = pair->toPAODomain;
  const Eigen::MatrixXd& toPNO_i = pair->singles_i->toPAODomain;
  const Eigen::MatrixXd& toPNO_j = pair->singles_j->toPAODomain;
  const Eigen::MatrixXd toPNO_T = toPNO.transpose().eval();
  const Eigen::MatrixXd toPNO_i_T = toPNO_i.transpose().eval();
  const Eigen::MatrixXd toPNO_j_T = toPNO_j.transpose().eval();
  const unsigned int nPNOs = toPNO.cols();

  // Calculate m_bck in PNO basis!
  // It would be faster to extract all three sets in one function,
  // but this step is far from being a bottleneck.
  std::vector<Eigen::MatrixXd> m_abK_ij = get_abK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                                  pair->domainProjection, k_PAOToFullPAOMaps, toPNO, toPNO, abK);
  std::vector<Eigen::MatrixXd> m_bcK_i = get_abK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_i,
                                                 pair->domainProjection, k_PAOToFullPAOMaps, toPNO_i, toPNO, abK);
  std::vector<Eigen::MatrixXd> m_bcK_j = get_abK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_j,
                                                 pair->domainProjection, k_PAOToFullPAOMaps, toPNO_j, toPNO, abK);
  // Calculate kak and construct ka_bc from kak and bck.
  for (unsigned int ik = 0; ik < pair->coupledPairs.size(); ++ik) {
    std::shared_ptr<CouplingOrbitalSet> coupledPair = pair->coupledPairs[ik];
    unsigned int k = coupledPair->getK();
    const Eigen::MatrixXd m_kaK = get_iaK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                          k_PAOToFullPAOMaps, reducedOccIndices, toPNO, k, iaK);
    const Eigen::MatrixXd llt_m_kaK = llt_metric.solve(m_kaK.transpose());
    for (auto& m_bck : m_abK_ij)
      coupledPair->ka_bc.push_back((m_bck * llt_m_kaK).transpose().eval());
  } // for ik
  // calculate ib_ac, ib_ac_i, jb_ac and jb_ac_j
  Eigen::MatrixXd m_iaK = get_iaK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                  k_PAOToFullPAOMaps, reducedOccIndices, toPNO, pair->i, iaK);
  Eigen::MatrixXd m_jaK = get_iaK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                  k_PAOToFullPAOMaps, reducedOccIndices, toPNO, pair->j, iaK);
  Eigen::MatrixXd llt_m_iaK = llt_metric.solve(m_iaK.transpose());
  Eigen::MatrixXd llt_m_jaK = llt_metric.solve(m_jaK.transpose());
  pair->ia_bc = {};
  pair->ja_bc = {};

  for (auto& m_bck : m_bcK_i)
    pair->ja_bc.push_back((m_bck * llt_m_jaK).transpose());
  for (auto& m_bck : m_bcK_j)
    pair->ia_bc.push_back((m_bck * llt_m_iaK).transpose());
  // Calculate (jc|ab) and (ic|ab) integrals.
  const Eigen::MatrixXd m_jck_i = get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_i,
                                          k_PAOToFullPAOMaps, reducedOccIndices, toPNO_i, pair->j, iaK);
  const Eigen::MatrixXd llt_m_jck_i = llt_metric.solve(m_jck_i.transpose());
  for (unsigned int cPNO = 0; cPNO < toPNO_i.cols(); ++cPNO) {
    Eigen::MatrixXd ab_jc = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
    for (unsigned int bPNO = 0; bPNO < nPNOs; ++bPNO) {
      ab_jc.col(bPNO) = m_abK_ij[bPNO] * llt_m_jck_i.col(cPNO);
    } // for bPNO
    assert((ab_jc - ab_jc.transpose()).array().abs().sum() < 1e-7);
    pair->jc_ab.push_back(ab_jc);
  } // for cPNO
  const Eigen::MatrixXd m_ick_j = get_iaK(auxBasisController, nLocalAux, pairDomainToK, domainProjection_j,
                                          k_PAOToFullPAOMaps, reducedOccIndices, toPNO_j, pair->i, iaK);
  const Eigen::MatrixXd llt_m_ick_j = llt_metric.solve(m_ick_j.transpose());
  for (unsigned int cPNO = 0; cPNO < toPNO_j.cols(); ++cPNO) {
    Eigen::MatrixXd ab_ic = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
    for (unsigned int bPNO = 0; bPNO < nPNOs; ++bPNO) {
      ab_ic.col(bPNO) = m_abK_ij[bPNO] * llt_m_ick_j.col(cPNO);
    } // for bPNO
    assert((ab_ic - ab_ic.transpose()).array().abs().sum() < 1e-7);
    pair->ic_ab.push_back(ab_ic);
  } // for cPNO
}

void Ao2MoExchangeIntegralTransformer::calculate_acbd_integrals(
    std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
    const Eigen::LLT<Eigen::MatrixXd>& llt_metric, const Eigen::SparseVector<int>& pairDomainToK,
    const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps, const MO3CenterIntegrals& abK) {
  const unsigned int nLocalAux = llt_metric.cols();
  const Eigen::MatrixXd& toPNO = pair->toPAODomain;
  const unsigned int nPNOs = toPNO.cols();
  pair->ac_bd = std::make_unique<Matrix<Eigen::MatrixXd>>(nPNOs, nPNOs, Eigen::MatrixXd::Zero(nPNOs, nPNOs));

  const std::vector<Eigen::MatrixXd> m_acks = get_abK(auxBasisController, nLocalAux, pairDomainToK, pair->domainProjection,
                                                      pair->domainProjection, k_PAOToFullPAOMaps, toPNO, toPNO, abK);
  for (unsigned int bPNO = 0; bPNO < nPNOs; ++bPNO) {
    const Eigen::MatrixXd invV_k_bd = llt_metric.solve(m_acks[bPNO].transpose().eval());
    for (unsigned int aPNO = 0; aPNO <= bPNO; ++aPNO) {
      const Eigen::MatrixXd ac_bd = m_acks[aPNO] * invV_k_bd;
      (*pair->ac_bd)(aPNO, bPNO) = ac_bd;
      if (bPNO != aPNO)
        (*pair->ac_bd)(bPNO, aPNO) = ac_bd.transpose();
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
  calculateTwoCenterIntegrals(M);
  const auto sparseMaps = mo3CenterIntegralController->getSparseMapsController();
  const SparseMap& occToK = sparseMaps->getOccToAuxShellMap();

  const auto paoMapAndIndices = mo3CenterIntegralController->getProjectionAndIndices(ORBITAL_TYPE::VIRTUAL);
  const auto occMapAndIndices = mo3CenterIntegralController->getProjectionAndIndices(ORBITAL_TYPE::OCCUPIED);

  const auto& k_redToFullMapsPAO = paoMapAndIndices.first;
  const auto& reducedOccIndices = occMapAndIndices.second;

  const auto& auxBasis = auxBasisController->getBasis();
  // Temporary precalculation of all integrals.
  // This could be change to a block wise calculation and writing.
  Eigen::SparseVector<int> auxSuperDomain = Eigen::VectorXi::Constant(auxBasis.size(), 1).sparseView();
  const MO3CenterIntegrals& iaK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ia_K, auxSuperDomain);
  OutputControl::nOut << "  Performing (ia|jb) integral transformation      ...";
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

    std::vector<Eigen::Triplet<double>> projectionTriplets;
    unsigned int row = 0;
    for (Eigen::SparseVector<int>::InnerIterator itPAO(pair->paoDomain); itPAO; ++itPAO) {
      projectionTriplets.push_back(Eigen::Triplet<double>(row, itPAO.row(), 1.0));
      row++;
    } // for itMu
    Eigen::SparseMatrix<double> projectionMatrix(pair->paoDomain.nonZeros(), pair->paoDomain.rows());
    projectionMatrix.setFromTriplets(projectionTriplets.begin(), projectionTriplets.end());

    pair->domainProjection = projectionMatrix;

    /*
     * Build local metric for each pair by selecting a block from S_aux.
     */
    Eigen::SparseVector<int> pairDomainToK = (occToK.col(i) + occToK.col(j)).pruned();
    Eigen::MatrixXd localBlock = SystemSplittingTools<scfMode>::getMatrixBlockShellWise(M, pairDomainToK, pairDomainToK);

    unsigned int nLocalAux = localBlock.cols();
    unsigned int nLocalPAOs = pair->paoDomain.nonZeros();
    Eigen::MatrixXd m_i = Eigen::MatrixXd::Zero(nLocalPAOs, nLocalAux);
    Eigen::MatrixXd m_j = Eigen::MatrixXd::Zero(nLocalPAOs, nLocalAux);
    unsigned int kCounter = 0;
    for (Eigen::SparseVector<int>::InnerIterator itK(pairDomainToK); itK; ++itK) {
      unsigned int extendedK = auxBasisController->extendedIndex(itK.row());
      unsigned int nContK = auxBasis[itK.row()]->getNContracted();
      const int i_red = reducedOccIndices[itK.row()](i);
      const int j_red = reducedOccIndices[itK.row()](j);
      if (i_red < 0 || j_red < 0 || k_redToFullMapsPAO[itK.row()].cols() < 1 || iaK[extendedK].rows() < 1 ||
          iaK[extendedK].cols() < 1) {
        kCounter += nContK;
        continue;
      }
      const Eigen::SparseMatrix<double> finalPAOProjection = (projectionMatrix * k_redToFullMapsPAO[itK.row()]).eval();
      for (unsigned int kk = 0; kk < nContK; ++kk) {
        m_i.col(kCounter) = finalPAOProjection * iaK[extendedK + kk].col(i_red);
        m_j.col(kCounter) = finalPAOProjection * iaK[extendedK + kk].col(j_red);
        kCounter++;
      } // for kk
    }   // for itK
    // The actual Cholesky decomposition and solving of the set of linear equations.
    pair->nAuxFunctions = m_i.cols();
    pair->k_ij = pair->toPAODomain.transpose() * m_i *
                 localBlock.llt().solve((pair->toPAODomain.transpose() * m_j).eval().transpose());
    ;
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
                                                             std::vector<std::shared_ptr<OrbitalPair>>& orbitalPairs) {
  Timings::takeTime("Local Cor. -   Int. Transform.");
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  // Temporary precalculation of all integrals.
  // This could be change to a block wise calculation and writing.
  Eigen::SparseVector<int> auxSuperDomain = Eigen::VectorXi::Constant(auxBasisController->getBasis().size(), 1).sparseView();
  const MO3CenterIntegrals& iaK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ia_K, auxSuperDomain);
  const MO3CenterIntegrals& abK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ab_K, auxSuperDomain);
  const MO3CenterIntegrals& klK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::kl_K, auxSuperDomain);

  const auto sparseMaps = mo3CenterIntegralController->getSparseMapsController();
  const SparseMap& occToK = sparseMaps->getOccToAuxShellMap();
  const SparseMap& extendedOccToK = sparseMaps->getCloseExtendedOccToAuxShellMap();

  const auto paoMapAndIndices = mo3CenterIntegralController->getProjectionAndIndices(ORBITAL_TYPE::VIRTUAL);
  const auto occMapAndIndices = mo3CenterIntegralController->getProjectionAndIndices(ORBITAL_TYPE::OCCUPIED);
  const auto& k_redToFullMapsPAO = paoMapAndIndices.first;
  const auto& reducedOccIndices = occMapAndIndices.second;

  unsigned int nThreads = 1;
#ifdef _OPENMP
  nThreads = omp_get_max_threads();
  Eigen::setNbThreads(1);
#endif
  MatrixInBasis<scfMode> M(auxBasisController);
  calculateTwoCenterIntegrals(M);
  OutputControl::nOut << "  Entering final Integral Transformation          ...";
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
    const Eigen::SparseVector<int> extendedPairDomainToK = (extendedOccToK.col(i) + extendedOccToK.col(j)).pruned();
    const Eigen::MatrixXd localBlock = SystemSplittingTools<scfMode>::getMatrixBlockShellWise(M, pairDomainToK, pairDomainToK);
    // Calculate the Cholesky decomposition of the local metric.
    const Eigen::LLT<Eigen::MatrixXd> llt = localBlock.llt();
    // Calculate V^(-1) * (K|ij) <-- needed in at least two transformations.
    const Eigen::VectorXd invV_ijK = calculate_invV_ijK(pair, pairDomainToK, auxBasisController, llt, reducedOccIndices, klK);
    /*
     * Calculate all integrals:
     * Note: There are a lot of integrals ...
     * Notation:
     * i,j            -> occ. orbitals of the given pair.
     * k,l            -> occ. orbitals of pairs coupling to the pair ij.
     * a,b,c,d        -> PNOs of the pair ij if not said otherwise.
     */
    if (nThreads == 1)
      Timings::takeTime("(ac|bd) integrals");
    calculate_acbd_integrals(pair, auxBasisController, llt, pairDomainToK, k_redToFullMapsPAO, abK);
    if (nThreads == 1)
      Timings::timeTaken("(ac|bd) integrals");
    if (nThreads == 1)
      Timings::takeTime("(ka|bc) 3 virt. integrals");
    calculate_kabcANDib_acANDjb_ac_integrals(pair, auxBasisController, llt, pairDomainToK, k_redToFullMapsPAO,
                                             reducedOccIndices, abK, iaK);
    if (nThreads == 1)
      Timings::timeTaken("(ka|bc) 3 virt. integrals");
    if (nThreads == 1)
      Timings::takeTime("(ij|ab) integrals");
    calculate_ijab_integrals(pair, auxBasisController, pairDomainToK, k_redToFullMapsPAO, abK, invV_ijK);
    if (nThreads == 1)
      Timings::timeTaken("(ij|ab) integrals");
    if (nThreads == 1)
      Timings::takeTime("(ik|ja) 3 occ. integrals");
    calculate_ikjaANDikjlANDkilaANDkjla_integrals(pair, auxBasisController, llt, pairDomainToK, k_redToFullMapsPAO,
                                                  reducedOccIndices, iaK, klK);
    if (nThreads == 1)
      Timings::timeTaken("(ik|ja) 3 occ. integrals");
    if (nThreads == 1)
      Timings::takeTime("(ij|ka) integrals");
    calculate_ijka_integrals(pair, auxBasisController, pairDomainToK, k_redToFullMapsPAO, reducedOccIndices, iaK, invV_ijK);
    if (nThreads == 1)
      Timings::timeTaken("(ij|ka) integrals");
    if (nThreads == 1)
      Timings::takeTime("Mixed PNO-basis integrals");
    calculateMixedIntegrals(pair, auxBasisController, M, k_redToFullMapsPAO, reducedOccIndices, occToK, iaK, abK, klK);
    if (nThreads == 1)
      Timings::timeTaken("Mixed PNO-basis integrals");
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
                                                             std::vector<OrbitalPairSet>& orbitalPairSets, bool dumpIntegrals) {
  for (unsigned int iSet = 0; iSet < orbitalPairSets.size(); ++iSet) {
    auto& orbitalPairSet = orbitalPairSets[iSet];
    // Transform integrals
    transformAllIntegrals(auxBasisController, mo3CenterIntegralController, orbitalPairSet);
    // Write to file
    if (dumpIntegrals) {
      takeTime("Writing Integral Sets to File");
      for (unsigned int iPair = 0; iPair < orbitalPairSet.size(); ++iPair) {
        auto& pair = orbitalPairSet[iPair];
        pair->writeIntegralsToFile();
        // Flush integrals if there are more than one set of orbital pairs.
        if (iSet != orbitalPairSets.size() - 1)
          pair->flushIntegrals();
      }
      timeTaken(2, "Writing Integral Sets to File");
    }
  } // for orbitalPairSets
}

} /* namespace Serenity */
