/**
 * @file OrbitalTriple.cpp
 *
 * @author Moritz Bensberg
 * @date Oct 31, 2019
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
#include "data/OrbitalTriple.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalPair.h"                                       //Orbital pair definition.
#include "data/SingleSubstitution.h"                                //Definition of a single substitution.
#include "integrals/MO3CenterIntegralController.h"                  //MO 3 center integrals.
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h" //Getter for integral lists.
#include "misc/SystemSplittingTools.h"                              //Extraction of coulomb metric block.

namespace Serenity {

OrbitalTriple::OrbitalTriple(std::shared_ptr<OrbitalPair> ikPair, std::shared_ptr<OrbitalPair> jkPair,
                             std::shared_ptr<OrbitalPair> ijPair, std::vector<std::shared_ptr<OrbitalPair>> ilPairs,
                             std::vector<std::shared_ptr<OrbitalPair>> jlPairs, std::vector<std::shared_ptr<OrbitalPair>> klPairs)
  : _ikPair(ikPair), _jkPair(jkPair), _ijPair(ijPair), _ilPairs(ilPairs), _jlPairs(jlPairs), _klPairs(klPairs) {
  _weakTriple = (_ikPair->type != OrbitalPairTypes::CLOSE || _jkPair->type != OrbitalPairTypes::CLOSE ||
                 _ijPair->type != OrbitalPairTypes::CLOSE);
  _i = (_ijPair->i == _ikPair->i) ? _ijPair->i : _ijPair->j;
  _k = (_jkPair->j == _ikPair->j) ? _jkPair->j : _jkPair->i;
  _j = (_ijPair->i == _i) ? _ijPair->j : _ijPair->i;
  if (_i == _j && _j == _k)
    throw SerenityError("ERROR: An artificial orbital triplet i==j==k was constructed.");
  _paoDomain = _ijPair->paoDomain + _ikPair->paoDomain + _jkPair->paoDomain;
  _iSingle = (_ijPair->i == _i) ? _ijPair->singles_i : _ijPair->singles_j;
  _jSingle = (_ijPair->j == _j) ? _ijPair->singles_j : _ijPair->singles_i;
  _kSingle = (_jkPair->j == _k) ? _jkPair->singles_j : _jkPair->singles_i;
  assert(_iSingle->i == _i);
  assert(_jSingle->i == _j);
  assert(_kSingle->i == _k);
  if (_iSingle->i != _i || _jSingle->i != _j || _kSingle->i != _k)
    throw SerenityError("Single extraction or index assignment failed in triple construction");
  _f_ii = _iSingle->getDiagonalPair()->f_ij;
  _f_jj = _jSingle->getDiagonalPair()->f_ij;
  _f_kk = _kSingle->getDiagonalPair()->f_ij;
  for (const auto& ilPair : ilPairs)
    _ilIndices.push_back((ilPair->i == _i) ? ilPair->j : ilPair->i);
  for (const auto& jlPair : jlPairs)
    _jlIndices.push_back((jlPair->i == _j) ? jlPair->j : jlPair->i);
  for (const auto& klPair : klPairs)
    _klIndices.push_back((klPair->i == _k) ? klPair->j : klPair->i);
}

void OrbitalTriple::cleanUp() {
  for (auto ints : _ia_bc)
    ints.resize(0, 0);
  _ia_bc = {};
  for (auto ints : _ja_bc)
    ints.resize(0, 0);
  _ja_bc = {};
  for (auto ints : _ka_bc)
    ints.resize(0, 0);
  _ka_bc = {};
  _ja_il.resize(0, 0);
  _ia_jl.resize(0, 0);
  _ka_il.resize(0, 0);
  _ia_kl.resize(0, 0);
  _ka_jl.resize(0, 0);
  _ja_kl.resize(0, 0);
  _ia_jb.resize(0, 0);
  _ia_kb.resize(0, 0);
  _ja_kb.resize(0, 0);

  _s_ij_ijk.resize(0, 0);
  _s_ik_ijk.resize(0, 0);
  _s_jk_ijk.resize(0, 0);
  _s_i_ijk.resize(0, 0);
  _s_j_ijk.resize(0, 0);
  _s_k_ijk.resize(0, 0);
  for (auto ints : _s_il_ijk)
    ints.resize(0, 0);
  _s_il_ijk = {};
  for (auto ints : _s_kl_ijk)
    ints.resize(0, 0);
  _s_kl_ijk = {};
  for (auto ints : _s_jl_ijk)
    ints.resize(0, 0);
  _s_jl_ijk = {};

  _toPAODomain.resize(0, 0);
  _tnoEigenvalues.resize(0);
  _tnoInit = false;
  _integralsCalc = false;
}

double OrbitalTriple::calculateTripleEnergy() {
  assert(_integralsCalc);
  const unsigned int nTNOs = _ia_bc.size();
  // Note: Capital letters are used in order to be as close to the reference publication as possible.
  std::vector<Eigen::MatrixXd> W = calculateW();
  std::vector<Eigen::MatrixXd> V = calculateV(W);
  std::vector<Eigen::MatrixXd> X = calculateX(W, V);
  std::vector<Eigen::MatrixXd> Y = calculateY(V);
  std::vector<Eigen::MatrixXd> Z = calculateZ(V);
  double delta_ij = (_i == _j) ? 1.0 : 0.0;
  double delta_jk = (_j == _k) ? 1.0 : 0.0;
  double P_ijk = 2 - delta_ij - delta_jk;
  assert(P_ijk > 1e-1 && "This triplet does not exist!");
  double f_ijk = _f_ii + _f_jj + _f_kk;
  double triplesEnergy = 0.0;
  for (unsigned int a = 0; a < nTNOs; ++a) {
    const double eps_a = _tnoEigenvalues(a);
    for (unsigned int b = 0; b <= a; ++b) {
      const double eps_ab = eps_a + _tnoEigenvalues(b);
      for (unsigned int c = 0; c <= b; ++c) {
        const double eps_abc = eps_ab + _tnoEigenvalues(c);
        double tmp = (Y[c](a, b) - 2.0 * Z[c](a, b)) * (W[c](a, b) + W[a](b, c) + W[b](c, a));
        tmp += (Z[c](a, b) - 2 * Y[c](a, b)) * (W[b](a, c) + W[c](b, a) + W[a](c, b));
        tmp += 3.0 * X[c](a, b);
        triplesEnergy += tmp / (f_ijk - eps_abc);
      } // for c <= b
    }   // for b <= a
  }     // for a
  triplesEnergy *= P_ijk;
  cleanUp();
  return triplesEnergy;
}

std::vector<Eigen::MatrixXd> OrbitalTriple::calculateZ(const std::vector<Eigen::MatrixXd>& V) {
  const unsigned int nTNOs = _ia_bc.size();
  std::vector<Eigen::MatrixXd> Z_abc(nTNOs, Eigen::MatrixXd::Zero(nTNOs, nTNOs));
  for (unsigned int c = 0; c < nTNOs; ++c) {
    Z_abc[c] = V[c].transpose(); // V_bac
    for (unsigned int b = 0; b < nTNOs; ++b) {
      Z_abc[c].col(b) += V[b].col(c); // V_acb
      for (unsigned int a = 0; a < nTNOs; ++a) {
        Z_abc[c](a, b) += V[a](c, b); // Inefficient!
      }                               // for a
    }                                 // for b
  }                                   // for c
  return Z_abc;
}

std::vector<Eigen::MatrixXd> OrbitalTriple::calculateY(const std::vector<Eigen::MatrixXd>& V) {
  const unsigned int nTNOs = _ia_bc.size();
  std::vector<Eigen::MatrixXd> Y_abc = V; // Copy! --> V_abc term
  for (unsigned int c = 0; c < nTNOs; ++c) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(nTNOs, nTNOs);
    for (unsigned int a = 0; a < nTNOs; ++a) {
      tmp.col(a) = V[a].col(c); // V_bca
    }                           // for a
    Y_abc[c] += tmp.transpose();
    for (unsigned int b = 0; b < nTNOs; ++b) {
      Y_abc[c].col(b) += V[b].row(c).transpose(); // V_cab
    }                                             // for b
  }                                               // for c
  return Y_abc;
}

std::vector<Eigen::MatrixXd> OrbitalTriple::calculateX(const std::vector<Eigen::MatrixXd>& W,
                                                       const std::vector<Eigen::MatrixXd>& V) {
  const unsigned int nTNOs = _ia_bc.size();
  std::vector<Eigen::MatrixXd> X_abc;
  for (unsigned int c = 0; c < nTNOs; ++c) {
    X_abc.push_back((W[c].array() * V[c].array()).eval());                 // first
    X_abc[c] += (W[c].array() * V[c].array()).matrix().transpose().eval(); // third
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(nTNOs, nTNOs);
    for (unsigned int b = 0; b < nTNOs; ++b) {
      tmp.col(b).array() += (W[b].col(c).array() * V[b].col(c).array()).eval();              // second and fourth
      tmp.col(b) += (W[b].row(c).array() * V[b].row(c).array()).matrix().transpose().eval(); // last two
    }                                                                                        // for b
    X_abc[c] += tmp + tmp.transpose();
  } // for c
  return X_abc;
}

std::vector<Eigen::MatrixXd> OrbitalTriple::calculateV(const std::vector<Eigen::MatrixXd>& W) {
  const unsigned int nTNOs = _ia_bc.size();
  std::vector<Eigen::MatrixXd> V_abc = W; // Copy! W_abc term
  const Eigen::VectorXd t_i = _s_i_ijk.transpose() * _iSingle->t_i;
  const Eigen::VectorXd t_j = _s_j_ijk.transpose() * _jSingle->t_i;
  const Eigen::VectorXd t_k = _s_k_ijk.transpose() * _kSingle->t_i;
  Eigen::MatrixXd P_abc = Eigen::MatrixXd::Constant(nTNOs, nTNOs, 1.0);
  P_abc.diagonal().array() += 1.0;
  for (unsigned int c = 0; c < nTNOs; ++c) {
    V_abc[c] += t_i * _ja_kb.col(c).transpose();
    V_abc[c] += _ia_kb.col(c) * t_j.transpose();
    V_abc[c] += t_k(c) * _ia_jb;
    Eigen::MatrixXd finalP_abc = P_abc;
    finalP_abc.col(c).array() += 1.0;
    V_abc[c].array() /= finalP_abc.array();
  } // for c
  return V_abc;
}

std::vector<Eigen::MatrixXd> OrbitalTriple::calculateW() {
  const unsigned int nTNOs = _ia_bc.size();
  std::vector<Eigen::MatrixXd> W_abc(nTNOs, Eigen::MatrixXd::Zero(nTNOs, nTNOs));
  const Eigen::MatrixXd t_kj = (_jkPair->i == _k) ? (_s_jk_ijk.transpose() * _jkPair->t_ij * _s_jk_ijk).eval()
                                                  : (_s_jk_ijk.transpose() * _jkPair->t_ij.transpose() * _s_jk_ijk).eval();
  const Eigen::MatrixXd t_jk = t_kj.transpose();
  const Eigen::MatrixXd t_ij = (_ijPair->i == _i) ? (_s_ij_ijk.transpose() * _ijPair->t_ij * _s_ij_ijk).eval()
                                                  : (_s_ij_ijk.transpose() * _ijPair->t_ij.transpose() * _s_ij_ijk).eval();
  const Eigen::MatrixXd t_ji = t_ij.transpose();
  const Eigen::MatrixXd t_ik = (_ikPair->i == _i) ? (_s_ik_ijk.transpose() * _ikPair->t_ij * _s_ik_ijk).eval()
                                                  : (_s_ik_ijk.transpose() * _ikPair->t_ij.transpose() * _s_ik_ijk).eval();
  const Eigen::MatrixXd t_ki = t_ik.transpose();

  // Sum_d part. All permutations.
  for (unsigned int d = 0; d < nTNOs; ++d) {
    const Eigen::MatrixXd ad_kb = _ka_bc[d].transpose();
    for (unsigned int c = 0; c < nTNOs; ++c) {
      // sum_d t^kj_cd (ia|bd)
      W_abc[c] += t_kj(c, d) * _ia_bc[d];                     //(abc--ijk) permutation.
      W_abc[c] += t_ki(c, d) * _ja_bc[d].transpose();         //(bac--jik)
      W_abc[c] += t_ij.col(d) * ad_kb.col(c).transpose();     //(cba--kji)
      W_abc[c] += t_ik.col(d) * _ja_bc[d].col(c).transpose(); //(bca-jki)
      W_abc[c] += _ia_bc[d].col(c) * t_jk.col(d).transpose(); //(acb--ikj)
      W_abc[c] += ad_kb.col(c) * t_ji.col(d).transpose();     //(cab--kij)
    }                                                         // for c
  }                                                           // for d

  // Sum_l part. All permutations.
  // Contractions with il-amplitudes.
  for (unsigned int ilIndex = 0; ilIndex < _ilPairs.size(); ++ilIndex) {
    const auto& ilPair = _ilPairs[ilIndex];
    const Eigen::MatrixXd t_il =
        (ilPair->i == _i) ? (_s_il_ijk[ilIndex].transpose() * ilPair->t_ij * _s_il_ijk[ilIndex]).eval()
                          : (_s_il_ijk[ilIndex].transpose() * ilPair->t_ij.transpose() * _s_il_ijk[ilIndex]).eval();
    for (unsigned int c = 0; c < nTNOs; ++c) {
      W_abc[c] -= _ka_jl(c, ilIndex) * t_il;                     //(abc--ijk)
      W_abc[c] -= t_il.col(c) * _ja_kl.col(ilIndex).transpose(); //(acb--ikj)
    }                                                            // for c
  }                                                              // for ilIndex

  // Contractions with kl-amplitudes
  for (unsigned int klIndex = 0; klIndex < _klPairs.size(); ++klIndex) {
    const auto& klPair = _klPairs[klIndex];
    const Eigen::MatrixXd t_lk =
        (klPair->j == _k) ? (_s_kl_ijk[klIndex].transpose() * klPair->t_ij * _s_kl_ijk[klIndex]).eval()
                          : (_s_kl_ijk[klIndex].transpose() * klPair->t_ij.transpose() * _s_kl_ijk[klIndex]).eval();
    for (unsigned int c = 0; c < nTNOs; ++c) {
      W_abc[c] -= t_lk.col(c) * _ja_il.col(klIndex).transpose(); //(cab--kij)
      W_abc[c] -= _ia_jl.col(klIndex) * t_lk.col(c).transpose(); //(cba--kji)
    }                                                            // for c
  }                                                              // for klIndex

  // Contractions with jl-amplitudes
  for (unsigned int jlIndex = 0; jlIndex < _jlPairs.size(); ++jlIndex) {
    const auto& jlPair = _jlPairs[jlIndex];
    const Eigen::MatrixXd t_jl =
        (jlPair->i == _j) ? (_s_jl_ijk[jlIndex].transpose() * jlPair->t_ij * _s_jl_ijk[jlIndex]).eval()
                          : (_s_jl_ijk[jlIndex].transpose() * jlPair->t_ij.transpose() * _s_jl_ijk[jlIndex]).eval();
    const Eigen::MatrixXd t_lj = t_jl.transpose();
    for (unsigned int c = 0; c < nTNOs; ++c) {
      W_abc[c] -= _ia_kl.col(jlIndex) * t_jl.col(c).transpose(); //(bca--jki)
      W_abc[c] -= t_lj * _ka_il(c, jlIndex);                     //(bac--jik)
    }                                                            // for c
  }                                                              // for jlIndex

  return W_abc;
}

void OrbitalTriple::calculateIntegrals(std::shared_ptr<BasisController> auxBasisController,
                                       std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                                       const MatrixInBasis<RESTRICTED>& coulombMetric) {
  if (!_integralsCalc) {
    const auto& sparseMaps = mo3CenterIntegralController->getSparseMapsController();
    const SparseMap& occToK = sparseMaps->getTripletOccToAuxShellMap(_weakTriple);
    const Eigen::SparseVector<int> fittingDomain = (occToK.col(_i) + occToK.col(_j) + occToK.col(_k)).pruned().eval();
    const Eigen::MatrixXd localBlock =
        SystemSplittingTools<RESTRICTED>::getMatrixBlockShellWise(coulombMetric, fittingDomain, fittingDomain);
    const unsigned int nLocalAux = localBlock.cols();
    // Calculate the Cholesky decomposition of the local metric.
    const Eigen::LLT<Eigen::MatrixXd> llt = localBlock.llt();

    std::vector<Eigen::Triplet<double>> projectionTriplets;
    unsigned int row = 0;
    for (Eigen::SparseVector<int>::InnerIterator itPAO(_paoDomain); itPAO; ++itPAO) {
      projectionTriplets.push_back(Eigen::Triplet<double>(row, itPAO.row(), 1.0));
      row++;
    } // for itMu
    Eigen::SparseMatrix<double> projectionMatrix(_paoDomain.nonZeros(), _paoDomain.rows());
    projectionMatrix.setFromTriplets(projectionTriplets.begin(), projectionTriplets.end());

    const auto& paoMapAndIndices = mo3CenterIntegralController->getProjectionAndIndices(ORBITAL_TYPE::VIRTUAL);
    const auto& occMapAndIndices = mo3CenterIntegralController->getProjectionAndIndices(ORBITAL_TYPE::OCCUPIED);
    const auto& k_redToFullMapsPAO = paoMapAndIndices.first;
    const auto& reducedOccIndices = occMapAndIndices.second;

    Eigen::SparseVector<int> auxSuperDomain =
        Eigen::VectorXi::Constant(auxBasisController->getBasis().size(), 1).sparseView();
    const MO3CenterIntegrals& iaK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ia_K, auxSuperDomain);
    const MO3CenterIntegrals& abK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ab_K, auxSuperDomain);
    const MO3CenterIntegrals& klK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::kl_K, auxSuperDomain);

    // TODO write get_iaK and get_ikK functions for more than one occupied index.
    // Note: This is currently not a bottleneck or close to one. Thus, I will ignore it for now.
    const Eigen::MatrixXd ibQ =
        Ao2MoExchangeIntegralTransformer::get_iaK(auxBasisController, nLocalAux, fittingDomain, projectionMatrix,
                                                  k_redToFullMapsPAO, reducedOccIndices, _toPAODomain, _i, iaK);
    const Eigen::MatrixXd jbQ =
        Ao2MoExchangeIntegralTransformer::get_iaK(auxBasisController, nLocalAux, fittingDomain, projectionMatrix,
                                                  k_redToFullMapsPAO, reducedOccIndices, _toPAODomain, _j, iaK);
    const Eigen::MatrixXd kbQ =
        Ao2MoExchangeIntegralTransformer::get_iaK(auxBasisController, nLocalAux, fittingDomain, projectionMatrix,
                                                  k_redToFullMapsPAO, reducedOccIndices, _toPAODomain, _k, iaK);

    const Eigen::MatrixXd ibQ_llt = llt.solve(ibQ.transpose()).eval().transpose();
    const Eigen::MatrixXd jbQ_llt = llt.solve(jbQ.transpose()).eval().transpose();
    const Eigen::MatrixXd kbQ_llt = llt.solve(kbQ.transpose()).eval().transpose();

    const std::vector<Eigen::MatrixXd> acQ =
        Ao2MoExchangeIntegralTransformer::get_abK(auxBasisController, nLocalAux, fittingDomain, projectionMatrix,
                                                  projectionMatrix, k_redToFullMapsPAO, _toPAODomain, _toPAODomain, abK);
    const Eigen::MatrixXd ilQ = Ao2MoExchangeIntegralTransformer::get_ikK(auxBasisController, nLocalAux, fittingDomain,
                                                                          reducedOccIndices, _ilIndices, _i, klK);
    //(ib|ac), (jb|ac), (kb|ac)
    for (const auto& m_acQ : acQ) {
      _ia_bc.push_back((ibQ_llt * m_acQ.transpose()).eval());
      _ja_bc.push_back((jbQ_llt * m_acQ.transpose()).eval());
      _ka_bc.push_back((kbQ_llt * m_acQ.transpose()).eval());
    }
    //(il|ja), (jl|ia),(il|ka),(kl|ia),(jl|ka),(kl|ja)
    {
      const Eigen::MatrixXd ilQ_kl = Ao2MoExchangeIntegralTransformer::get_ikK(auxBasisController, nLocalAux, fittingDomain,
                                                                               reducedOccIndices, _klIndices, _i, klK);
      const Eigen::MatrixXd jlQ_kl = Ao2MoExchangeIntegralTransformer::get_ikK(auxBasisController, nLocalAux, fittingDomain,
                                                                               reducedOccIndices, _klIndices, _j, klK);
      _ja_il = (jbQ_llt * ilQ_kl.transpose()).eval(); // contracted with kl
      _ia_jl = (ibQ_llt * jlQ_kl.transpose()).eval(); // contracted with kl
    }
    {
      const Eigen::MatrixXd ilQ_jl = Ao2MoExchangeIntegralTransformer::get_ikK(auxBasisController, nLocalAux, fittingDomain,
                                                                               reducedOccIndices, _jlIndices, _i, klK);
      const Eigen::MatrixXd klQ_jl = Ao2MoExchangeIntegralTransformer::get_ikK(auxBasisController, nLocalAux, fittingDomain,
                                                                               reducedOccIndices, _jlIndices, _k, klK);
      _ka_il = (kbQ_llt * ilQ_jl.transpose()).eval(); // contracted with jl
      _ia_kl = (ibQ_llt * klQ_jl.transpose()).eval(); // contracted with jl
    }
    {
      const Eigen::MatrixXd jlQ_il = Ao2MoExchangeIntegralTransformer::get_ikK(auxBasisController, nLocalAux, fittingDomain,
                                                                               reducedOccIndices, _ilIndices, _j, klK);
      const Eigen::MatrixXd klQ_il = Ao2MoExchangeIntegralTransformer::get_ikK(auxBasisController, nLocalAux, fittingDomain,
                                                                               reducedOccIndices, _ilIndices, _k, klK);
      _ka_jl = (kbQ_llt * jlQ_il.transpose()).eval(); // contracted with il
      _ja_kl = (jbQ_llt * klQ_il.transpose()).eval(); // contracted with il
    }
    //(ia|jb),(ia|kb),(ja|kb)
    _ia_jb = ibQ * jbQ_llt.transpose();
    _ia_kb = ibQ * kbQ_llt.transpose();
    _ja_kb = jbQ * kbQ_llt.transpose();
    _integralsCalc = true;
  } // if !_integralsCalc
}

} /* namespace Serenity */
