/**
 * @file QuasiCanonicalPAODomainConstructor.cpp
 *
 * @date May 13, 2019
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
#include "analysis/PAOSelection/QuasiCanonicalPAODomainConstructor.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"                   //of the environment systems.
#include "data/OrbitalPair.h"                           //OrbitalPair definition
#include "data/PAOController.h"                         //PAOController
#include "misc/SystemSplittingTools.h"                  //Orthogonalization in non redundant basis.
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h" //Definition of CouplingOrbitalSet.
#include "potentials/LevelshiftPotential.h"             //Levelshift of environment orbitals
#include "system/SystemController.h"                    //Get the electronic structure

namespace Serenity {

QuasiCanonicalPAODomainConstructor::QuasiCanonicalPAODomainConstructor(
    const CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& coefficients,
    std::shared_ptr<const MatrixInBasis<Options::SCF_MODES::RESTRICTED>> overlapMatrix,
    std::shared_ptr<const FockMatrix<Options::SCF_MODES::RESTRICTED>> f, std::shared_ptr<PAOController> paoController,
    double paoOrthogonalizationThreshold, std::vector<std::shared_ptr<SystemController>> environmentSystems,
    double levelShiftParameter, double ssScaling, double osScaling)
  : _paoController(paoController),
    _S(overlapMatrix),
    _paoOrthogonalizationThreshold(paoOrthogonalizationThreshold),
    _coeff(coefficients),
    _ssScaling(ssScaling),
    _osScaling(osScaling) {
  assert(_paoController);
  assert(f);
  _f = *f;
  _f_MO = (coefficients.transpose() * _f * coefficients).eval();
  for (const auto& envSys : environmentSystems) {
    assert(envSys->getLastSCFMode() == Options::SCF_MODES::RESTRICTED);
    LevelshiftPotential<RESTRICTED> levelshift(f->getBasisController(),
                                               envSys->getElectronicStructure<RESTRICTED>()->getDensityMatrixController(),
                                               levelShiftParameter);
    _f += levelshift.getMatrix();
  } // for envSys
  _f_AO_MO = (_f * coefficients).eval();
}

void QuasiCanonicalPAODomainConstructor::initializeAmplitudes(std::shared_ptr<OrbitalPair> pair) {
  pair->t_ij = -(pair->k_ij.array() / pair->uncoupledTerm.array()).matrix();
  double ssEnergy =
      (pair->i == pair->j ? 1.0 : 2.0) * ((pair->t_ij - pair->t_ij.transpose()).array() * pair->k_ij.array()).sum();
  double osEnergy = (pair->i == pair->j ? 1.0 : 2.0) * (pair->t_ij.array() * pair->k_ij.array()).sum();
  double scMP2PairEnergy = _ssScaling * ssEnergy + _osScaling * osEnergy;
  pair->scMP2PairEnergy = scMP2PairEnergy;
}

void QuasiCanonicalPAODomainConstructor::transformToQuasiCanonicalPAOBasis(std::shared_ptr<OrbitalPair> pair) {
  // Get the PAO coefficients of the pair domain.
  const Eigen::MatrixXd& R_ij = _paoController->getPAOsFromDomain(pair->paoDomain);
  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd transformation;
  SystemSplittingTools<RESTRICTED>::diagonalizationInNonRedundantPAOBasis(R_ij, *_S, _f, _paoOrthogonalizationThreshold,
                                                                          eigenvalues, transformation);
  // Save transformation to linear independent PAO-pair basis that diagonalize the fock matrix (in this block).
  pair->toPAODomain = transformation;
  // Save matrix of eps_r + eps_s - f_ii - f_jj. This is needed in every amplitude opt. iteration.
  unsigned int i = pair->i;
  unsigned int j = pair->j;
  pair->uncoupledTerm = Eigen::MatrixXd::Zero(eigenvalues.size(), eigenvalues.size());
  pair->uncoupledTerm.colwise() += eigenvalues;
  pair->uncoupledTerm.rowwise() += eigenvalues.transpose();
  pair->uncoupledTerm.array() -= _f_MO(i, i) + _f_MO(j, j);
  pair->f_ij = _f_MO(i, j);
  pair->f_ab = eigenvalues;
}

} /* namespace Serenity */
