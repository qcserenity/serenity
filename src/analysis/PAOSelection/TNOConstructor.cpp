/**
 * @file TNOConstructor.cpp
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
#include "analysis/PAOSelection/TNOConstructor.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"       //of the environment systems.
#include "data/OrbitalPair.h"               //Orbital Pair definition
#include "data/OrbitalTriple.h"             //Orbital Triple definition
#include "data/PAOController.h"             //PAOController
#include "data/SingleSubstitution.h"        //Singles--triples overlap.
#include "misc/SystemSplittingTools.h"      //Orthogonalization in non redundant basis.
#include "potentials/LevelshiftPotential.h" //Levelshift of environment orbitals
#include "system/SystemController.h"        //Get the electronic structure

namespace Serenity {

TNOConstructor::TNOConstructor(std::shared_ptr<const MatrixInBasis<Options::SCF_MODES::RESTRICTED>> overlapMatrix,
                               std::shared_ptr<const FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                               std::shared_ptr<PAOController> paoController, double paoOrthogonalizationThreshold,
                               double tnoThreshold, double coreScaling,
                               std::vector<std::shared_ptr<SystemController>> environmentSystems, double levelShiftParameter)
  : _S(overlapMatrix),
    _paoController(paoController),
    _paoOrthogonalizationThreshold(paoOrthogonalizationThreshold),
    _tnoThreshold(tnoThreshold),
    _tnoCoreThreshold(tnoThreshold * coreScaling) {
  _f = *f;
  for (const auto& envSys : environmentSystems) {
    assert(envSys->getLastSCFMode() == Options::SCF_MODES::RESTRICTED);
    LevelshiftPotential<RESTRICTED> levelshift(f->getBasisController(),
                                               envSys->getElectronicStructure<RESTRICTED>()->getDensityMatrixController(),
                                               levelShiftParameter);
    _f += levelshift.getMatrix();
  } // for envSys
}

Eigen::MatrixXd TNOConstructor::constructOverlapMatrix(std::shared_ptr<OrbitalPair> pair, const Eigen::MatrixXd& s_p_ijk) {
  const Eigen::MatrixXd p_ij_T = (_paoController->getPAOsFromDomain(pair->paoDomain) * pair->toPAODomain).transpose();
  return p_ij_T * s_p_ijk;
}
Eigen::MatrixXd TNOConstructor::constructOverlapMatrix(std::shared_ptr<SingleSubstitution> single,
                                                       const Eigen::MatrixXd& s_p_ijk) {
  const Eigen::MatrixXd p_i_T =
      (_paoController->getPAOsFromDomain(single->getDiagonalPair()->paoDomain) * single->toPAODomain).transpose();
  return p_i_T * s_p_ijk;
}
Eigen::MatrixXd TNOConstructor::constructDensityMatrix(std::shared_ptr<OrbitalPair> pair) {
  double prefactor = (pair->i == pair->j) ? 0.5 : 1.0;
  const Eigen::MatrixXd t_ji = pair->t_ij.transpose();
  const Eigen::MatrixXd t_ij_tilde = prefactor * (4.0 * pair->t_ij - 2.0 * t_ji).eval();
  return (t_ij_tilde * t_ji + t_ij_tilde.transpose() * pair->t_ij).eval();
}

void TNOConstructor::transformToTNOBasis(std::shared_ptr<OrbitalTriple> triple) {
  // Transform to non-redundant PAO basis.
  const Eigen::MatrixXd& R_ijk = _paoController->getPAOsFromDomain(triple->getPAODomain());
  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd transformation;
  SystemSplittingTools<RESTRICTED>::diagonalizationInNonRedundantPAOBasis(R_ijk, *_S, _f, _paoOrthogonalizationThreshold,
                                                                          eigenvalues, transformation);

  const auto& ijPair = triple->getIJPair();
  const auto& ikPair = triple->getIKPair();
  const auto& jkPair = triple->getJKPair();
  // Construct quasi-canonical PAO overlap matrices.
  Eigen::MatrixXd p_ijk = R_ijk * transformation;
  Eigen::MatrixXd s_p_ijk = ((*_S) * p_ijk).eval();
  const Eigen::MatrixXd qcS_ij_ijk = constructOverlapMatrix(ijPair, s_p_ijk);
  const Eigen::MatrixXd qcS_ik_ijk = constructOverlapMatrix(ikPair, s_p_ijk);
  const Eigen::MatrixXd qcS_jk_ijk = constructOverlapMatrix(jkPair, s_p_ijk);
  // Construct density matrices
  const Eigen::MatrixXd d_ij = constructDensityMatrix(ijPair);
  const Eigen::MatrixXd d_ik = constructDensityMatrix(ikPair);
  const Eigen::MatrixXd d_jk = constructDensityMatrix(jkPair);
  // Construct triples density matrix.
  const Eigen::MatrixXd D_ijk = 1.0 / 3.0 *
                                (qcS_ij_ijk.transpose() * d_ij * qcS_ij_ijk + qcS_ik_ijk.transpose() * d_ik * qcS_ik_ijk +
                                 qcS_jk_ijk.transpose() * d_jk * qcS_jk_ijk);
  // Diagonalize triples density matrix.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(D_ijk);
  auto& d_ijk = es.eigenvectors();
  auto& n_ijk = es.eigenvalues();
  // Truncate by occupation number.
  int nTNOs = (n_ijk.array() >= _tnoThreshold).count();
  if (triple->getISingle()->coreLikeOrbital || triple->getJSingle()->coreLikeOrbital || triple->getKSingle()->coreLikeOrbital)
    nTNOs = (n_ijk.array() >= _tnoCoreThreshold).count();
  // Transform to quasi-canonical basis.
  const Eigen::MatrixXd d_ij_red = d_ijk.rightCols(nTNOs);
  p_ijk *= d_ij_red;
  Eigen::MatrixXd f_pao = (p_ijk.transpose() * _f * p_ijk).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ef(f_pao);

  // Set the transformation, eigenvalues and TNO-PNO overlap matrices.
  s_p_ijk = ((*_S) * p_ijk * ef.eigenvectors()).eval();
  std::vector<Eigen::MatrixXd> s_il_ijk;
  std::vector<Eigen::MatrixXd> s_jl_ijk;
  std::vector<Eigen::MatrixXd> s_kl_ijk;
  for (const auto& ilPair : triple->getILPairs())
    s_il_ijk.push_back(constructOverlapMatrix(ilPair, s_p_ijk));
  for (const auto& jlPair : triple->getJLPairs())
    s_jl_ijk.push_back(constructOverlapMatrix(jlPair, s_p_ijk));
  for (const auto& klPair : triple->getKLPairs())
    s_kl_ijk.push_back(constructOverlapMatrix(klPair, s_p_ijk));
  triple->setTNOTransformation(
      (transformation * d_ij_red * ef.eigenvectors()).eval(), ef.eigenvalues(), constructOverlapMatrix(ijPair, s_p_ijk),
      constructOverlapMatrix(ikPair, s_p_ijk), constructOverlapMatrix(jkPair, s_p_ijk), s_il_ijk, s_kl_ijk, s_jl_ijk,
      constructOverlapMatrix(triple->getISingle(), s_p_ijk), constructOverlapMatrix(triple->getJSingle(), s_p_ijk),
      constructOverlapMatrix(triple->getKSingle(), s_p_ijk));
}

} /* namespace Serenity */
