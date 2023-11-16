/**
 * @file PNOConstructor.cpp
 *
 * @date Apr 11, 2019
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
#include "analysis/PAOSelection/PNOConstructor.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalPair.h"        //OrbitalPair definition
#include "data/PAOController.h"      //PAOController
#include "data/SingleSubstitution.h" //Singles definition
#include "misc/WarningTracker.h"     //Warnings
/* Include Std and External Headers */
#include "Eigen/Dense" //Dense-Matrix definition

namespace Serenity {

PNOConstructor::PNOConstructor(const CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& coefficients,
                               std::shared_ptr<const FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                               std::shared_ptr<PAOController> paoController, double paoOrthogonalizationThreshold,
                               std::vector<std::shared_ptr<SystemController>> environmentSystems,
                               double levelShiftParameter, bool setFaiZero, double ssScaling, double osScaling)
  : QuasiCanonicalPAODomainConstructor(coefficients, f, paoController, paoOrthogonalizationThreshold,
                                       environmentSystems, levelShiftParameter, ssScaling, osScaling),
    _setFaiZero(setFaiZero),
    _ssScaling(ssScaling),
    _osScaling(osScaling){};

void PNOConstructor::transformExternalBasis(std::shared_ptr<OrbitalPair> pair) {
  transformToQuasiCanonicalPAOBasis(pair);
}

void PNOConstructor::transformToPNOBasis(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs) {
  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
    auto& pair = orbitalPairs[iPair];
    transformToPNOBasis(pair);
  } // for pair
  Eigen::setNbThreads(0);
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> PNOConstructor::orthogonalizeFockMatrix(const Eigen::MatrixXd& d_ij_red,
                                                                                    std::shared_ptr<OrbitalPair> pair) {
  const Eigen::MatrixXd p_ij = pair->toPAODomain * d_ij_red;
  const Eigen::MatrixXd tmp = pair->domainProjection * this->_f_pao * pair->domainProjection.transpose();
  const Eigen::MatrixXd f_pao = (p_ij.transpose() * tmp * p_ij).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ef(f_pao);
  Eigen::MatrixXd transformation = (d_ij_red * ef.eigenvectors()).eval();
  return std::pair<Eigen::VectorXd, Eigen::MatrixXd>(ef.eigenvalues(), transformation);
}

void PNOConstructor::postProcessing(std::shared_ptr<OrbitalPair> pair) {
  initializeAmplitudes(pair);
  transformToPNOBasis(pair);
}

void PNOConstructor::transformToPNOBasis(std::shared_ptr<OrbitalPair> pair) {
  unsigned int i = pair->i;
  unsigned int j = pair->j;
  assert(pair->toPAODomain.array().size() > 0);
  assert(pair->t_ij.array().size() > 0 && "No initial amplitudes given for PNO construction!");

  double prefactor = (pair->i == pair->j) ? 0.5 : 1.0;
  Eigen::MatrixXd t_ji = pair->t_ij.transpose();
  Eigen::MatrixXd t_ij_tilde = 4.0 * pair->t_ij - 2.0 * t_ji;
  const double fullPNOPairEnergy = pair->scMP2PairEnergy;
  Eigen::MatrixXd D_ij = prefactor * (t_ij_tilde * t_ji + t_ij_tilde.transpose() * pair->t_ij);
  double initialT2Norm = D_ij.trace();

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(D_ij);
  auto& d_ij = es.eigenvectors();
  auto& n_ij = es.eigenvalues();
  int nPNOs = (n_ij.array() >= pair->getPNOThreshold()).count();
  if (pair->singles_i && pair->singles_j) {
    if ((pair->singles_i->coreLikeOrbital && not pair->singles_j->coreLikeOrbital) ||
        (pair->singles_j->coreLikeOrbital && not pair->singles_i->coreLikeOrbital)) {
      if (std::fabs(pair->f_ij) > 1e-6) {
        WarningTracker::printWarning(
            (std::string) "\n Small WARNING: Non-zero Fock matrix elements between core and non-core orbitals!\n" +
                "                  Results can change significantly if core and valence orbitals are\n" +
                "                  localized independently! F_ij= " + pair->f_ij,
            true);
      }
    }
  }
  if (nPNOs == 0) {
    if (pair->i != pair->j) {
      pair->type = OrbitalPairTypes::VERY_DISTANT;
      pair->deltaPNO = fullPNOPairEnergy;
      pair->t_ij.resize(0, 0);
      pair->k_ij.resize(0, 0);
      pair->uncoupledTerm.resize(0, 0);
      pair->f_ab.resize(0);
      return;
    }
    else {
      nPNOs = 1;
    }
  }

  auto epsAndTrafo = orthogonalizeFockMatrix(d_ij.rightCols(nPNOs).eval(), pair);
  const Eigen::VectorXd& eigenvalues = epsAndTrafo.first;
  const Eigen::MatrixXd& transformation = epsAndTrafo.second;

  Eigen::MatrixXd uncoupledTerm = Eigen::MatrixXd::Zero(eigenvalues.size(), eigenvalues.size());
  uncoupledTerm.colwise() += eigenvalues.eval();
  uncoupledTerm.rowwise() += eigenvalues.transpose().eval();
  uncoupledTerm.array() -= this->_f_MO(i, i) + this->_f_MO(j, j);

  pair->k_ij = (transformation.transpose() * pair->k_ij * transformation).eval();
  pair->uncoupledTerm = uncoupledTerm;
  pair->t_ij = -(pair->k_ij.array() / pair->uncoupledTerm.array()).matrix();
  pair->f_ab = eigenvalues;

  if (pair->i == pair->j && pair->singles_i) {
    auto single = pair->singles_i;
    int nSinglesPNOs = (n_ij.array() >= single->getPNOThreshold()).count();
    nSinglesPNOs = (nSinglesPNOs > 0) ? nSinglesPNOs : 1;
    auto epsAndTrafo_singles = orthogonalizeFockMatrix(d_ij.rightCols(nSinglesPNOs).eval(), pair);
    const Eigen::VectorXd& eigenvalues_singles = epsAndTrafo_singles.first;
    const Eigen::MatrixXd& transformation_singles = epsAndTrafo_singles.second;
    single->toPAODomain = (pair->toPAODomain * transformation_singles).eval();
    single->t_i = Eigen::VectorXd::Zero(eigenvalues_singles.rows());
    single->f_ai = (single->toPAODomain.transpose() * pair->domainProjection * _f_PAO_MO.col(i)).eval();
    if (_setFaiZero)
      single->f_ai.setZero();
    single->epsMinusF = eigenvalues_singles.array() - _f_MO(i, i);
    single->f_ab = eigenvalues_singles;
  }
  pair->toPAODomain = (pair->toPAODomain * transformation).eval();

  t_ij_tilde = prefactor * (4.0 * pair->t_ij - 2.0 * pair->t_ij.transpose());
  Eigen::MatrixXd finalD = t_ij_tilde * pair->t_ij.transpose() + t_ij_tilde.transpose() * pair->t_ij;
  double finalT2Norm = finalD.trace();
  pair->pnoNormError = 1.0 - (initialT2Norm - finalT2Norm) / initialT2Norm;
  double ssEnergy =
      (pair->i == pair->j ? 1.0 : 2.0) * ((pair->t_ij - pair->t_ij.transpose()).array() * pair->k_ij.array()).sum();
  double osEnergy = (pair->i == pair->j ? 1.0 : 2.0) * (pair->t_ij.array() * pair->k_ij.array()).sum();
  double truncPNOPairEnergy = _ssScaling * ssEnergy + _osScaling * osEnergy;
  pair->deltaPNO = fullPNOPairEnergy - truncPNOPairEnergy;
  pair->lMP2PairEnergy = truncPNOPairEnergy;
}

} /* namespace Serenity */
