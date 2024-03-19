/**
 * @file DipoleApproximationToPairEnergies.cpp
 *
 * @date Feb 1, 2019
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
#include "postHF/MPn/DipoleApproximationToPairEnergies.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h" //Loop shells.
#include "basis/BasisController.h"
#include "data/OrbitalPair.h"
#include "data/PAOController.h"
#include "integrals/wrappers/Libint.h" //Dipole integrals
#include "io/FormattedOutputStream.h"  //Filtered output streams.
#include "misc/SystemSplittingTools.h" //diagonalizationInNonRedundantPAOBasis
/* Include Std and External Headers */
#include <algorithm> //std::find

namespace Serenity {

DipoleApproximationToPairEnergies::DipoleApproximationToPairEnergies(
    std::shared_ptr<BasisController> basisController, std::shared_ptr<MatrixInBasis<RESTRICTED>> overlapMatrix,
    std::shared_ptr<Eigen::MatrixXd> occCoefficients, std::shared_ptr<PAOController> paoController,
    std::shared_ptr<FockMatrix<RESTRICTED>> fockMatrix,
    std::shared_ptr<Eigen::SparseMatrix<int>> paoToOccupiedOrbitalMap, double paoOrthogonalizationThreshold)
  : _basisController(basisController),
    _overlapMatrix(overlapMatrix),
    _occCoefficients(occCoefficients),
    _paoController(paoController),
    _fockMatrix(fockMatrix),
    _paoToOccupiedOrbitalMap(paoToOccupiedOrbitalMap),
    _paoOrthogonalizationThreshold(paoOrthogonalizationThreshold) {
  _pao_r_i.resize(occCoefficients->cols(), nullptr);
  _eigenvalues_ii.resize(occCoefficients->cols(), nullptr);
}

std::vector<Eigen::MatrixXd> DipoleApproximationToPairEnergies::calculateDipoleIntegrals() {
  unsigned int nBasisFunctions = _basisController->getNBasisFunctions();
  std::vector<Eigen::MatrixXd> results(3, Eigen::MatrixXd::Zero(nBasisFunctions, nBasisFunctions));
  auto& basis = _basisController->getBasis();
  auto& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::emultipole1, 0, 2);
  // Loop over shells i and i <= j.
  for (unsigned int iShell = 0; iShell < basis.size(); iShell++) {
    auto nContShellI = basis[iShell]->getNContracted();
    for (unsigned int jShell = 0; jShell <= iShell; jShell++) {
      auto nContShellJ = basis[jShell]->getNContracted();
      Eigen::MatrixXd multiPoleInts;
      if (libint.compute(LIBINT_OPERATOR::emultipole1, 0, *basis[jShell], *basis[iShell], multiPoleInts)) {
        /*
         * set vector size: libint will return a matrix containing
         * <j|i>,<j|x|i>,<j|y|i>,<j|z|i>
         * Order in raw buffer: (j_1|i_1),(j_1|i_2)...,(j_2|i_1) ... for each column.
         */
        unsigned int iExtended = _basisController->extendedIndex(iShell);
        unsigned int jExtended = _basisController->extendedIndex(jShell);
        for (unsigned int dim = 0; dim < 3; ++dim) {
          Eigen::MatrixXd resultBlock =
              Eigen::Map<Eigen::MatrixXd>(multiPoleInts.col(dim + 1).data(), nContShellI, nContShellJ);
          results[dim].block(iExtended, jExtended, nContShellI, nContShellJ) = resultBlock;
          // Symmetry of i and j.
          if (iShell != jShell)
            results[dim].block(jExtended, iExtended, nContShellJ, nContShellI) = resultBlock.transpose();
        } // for dim
      }   // if compute
    }     // for jShell
  }       // for iShell
  libint.finalize(LIBINT_OPERATOR::emultipole1, 0, 2);
  // Transform right side to the occupied orbital basis.
  const auto& occOrbCoeff = *_occCoefficients;
  for (auto& ints : results) {
    ints = ints * occOrbCoeff;
  }
  // Construct the (i|r|i) integral matrix.
  _iri.resize(3, occOrbCoeff.cols());
  _iri.row(0) = (occOrbCoeff.array() * results[0].array()).colwise().sum();
  _iri.row(1) = (occOrbCoeff.array() * results[1].array()).colwise().sum();
  _iri.row(2) = (occOrbCoeff.array() * results[2].array()).colwise().sum();
  return results;
}

inline void DipoleApproximationToPairEnergies::calculateTransformationAndEigenvalues(unsigned int i) {
  const Eigen::MatrixXd R_ii = _paoController->getPAOsFromDomain(_paoToOccupiedOrbitalMap->col(i));
  Eigen::MatrixXd transformation;
  _eigenvalues_ii[i] = std::make_shared<Eigen::VectorXd>();
  SystemSplittingTools<RESTRICTED>::diagonalizationInNonRedundantPAOBasis(
      R_ii, *_overlapMatrix, *_fockMatrix, _paoOrthogonalizationThreshold, *_eigenvalues_ii[i], transformation);
  transformation = R_ii * transformation;
  // Calculate integrals (pao_ii|r|i)
  _pao_r_i[i] = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd(transformation.cols(), 3));
  _pao_r_i[i]->col(0).noalias() = transformation.transpose() * _dipoleInts[0].col(i);
  _pao_r_i[i]->col(1).noalias() = transformation.transpose() * _dipoleInts[1].col(i);
  _pao_r_i[i]->col(2).noalias() = transformation.transpose() * _dipoleInts[2].col(i);
}

inline void DipoleApproximationToPairEnergies::preCalculateTransformations(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs) {
  std::vector<unsigned int> orbitalIndices;
  for (const auto& pair : orbitalPairs) {
    if (std::find(orbitalIndices.begin(), orbitalIndices.end(), pair->i) == orbitalIndices.end())
      orbitalIndices.push_back(pair->i);
    if (pair->i == pair->j)
      continue;
    if (std::find(orbitalIndices.begin(), orbitalIndices.end(), pair->j) == orbitalIndices.end())
      orbitalIndices.push_back(pair->j);
  }
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iOrb = 0; iOrb < orbitalIndices.size(); ++iOrb) {
    unsigned int i = orbitalIndices[iOrb];
    if (!_pao_r_i[i] || !_eigenvalues_ii[i])
      calculateTransformationAndEigenvalues(i);
  }
}

void DipoleApproximationToPairEnergies::calculateDipoleApproximation(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs) {
  takeTime("Dipole Approximation");
  /*
   * 1. Calculate all (mu|r|nu) integrals in AO basis.
   * For each pair
   *   2. Build Fock matrix in non-redundant PAO basis for i and j.
   *   3. Calculate all:
   *        (i|r|mu_pao),
   *        (j|r|mu_pao)
   *      and the two integrals (i|r|i) and (j|r|j)
   *   4. Calculate R_ij = (i|r|i) - (j|r|j)
   * 5. Sum all up according to equation 17 J.Chem.Phys. 143, 034108 (2015)
   */
  // Calculate all (mu|r|nu) integrals in AO basis if not already done.
  if (_dipoleInts.size() < 1)
    _dipoleInts = calculateDipoleIntegrals();
  preCalculateTransformations(orbitalPairs);
  // Transform the Fock matrix.
  Eigen::MatrixXd f_MO = _occCoefficients->transpose() * *_fockMatrix * *_occCoefficients;
  // Loop over pairs:
  unsigned int pairIndex = 0;
  unsigned int nThreads = omp_get_max_threads();
  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
    auto& pair = orbitalPairs[iPair];
    unsigned int i = pair->i;
    unsigned int j = pair->j;
    // Extract precalculated values.
    const Eigen::MatrixXd& muPAO_r_i = *_pao_r_i[i];
    const Eigen::MatrixXd& muPAO_r_j = *_pao_r_i[j];
    const Eigen::Vector3d r_ij = _iri.col(i) - _iri.col(j);
    Eigen::MatrixXd denominator = Eigen::MatrixXd::Zero(_eigenvalues_ii[i]->size(), _eigenvalues_ii[j]->size());
    denominator.colwise() += *_eigenvalues_ii[i];
    denominator.rowwise() += _eigenvalues_ii[j]->transpose().eval();
    denominator.array() -= f_MO(i, i) + f_MO(j, j);
    /*
     * Step 2 to 4.
     */
    /*
     * According to equation 17 J.Chem.Phys. 143, 034108 (2015)
     */
    // First term: (i|r|mu_PAO)(j|r|nu_PAO)
    // n_ix3 * 3xn_j = n_ixn_j
    Eigen::MatrixXd epsilon_mu_nu_i_j = muPAO_r_i * muPAO_r_j.transpose();
    const Eigen::Vector3d unit_r_ij = r_ij.array() / r_ij.norm();
    // Second term: -3((i|r|mu_PAO)r_ij)((j|r|nu_PAO)r_ij
    //                        (  n_ix3 * 3x1     *     (n_jx3 * 3x1)^T            )= n_ixn_j
    epsilon_mu_nu_i_j -= 3 * ((muPAO_r_i * unit_r_ij) * (muPAO_r_j * unit_r_ij).transpose());
    // Square it and divide by eps_mu+eps_nu-F_ii-F_jj
    epsilon_mu_nu_i_j.array() *= epsilon_mu_nu_i_j.array() / denominator.array();
    // The pair energy is given as the sum of all matrix elements multiplied with -4/r_ij^6
    double pairEnergy = epsilon_mu_nu_i_j.sum();
    const double r_ijSquared = r_ij.squaredNorm();
    pairEnergy *= -4.0 / (r_ijSquared * r_ijSquared * r_ijSquared);
    pair->dipolePairEnergy = pairEnergy;
    ++pairIndex;
  } // for pair
  Eigen::setNbThreads(nThreads);
  timeTaken(2, "Dipole Approximation");
}
void DipoleApproximationToPairEnergies::calculateDipoleApproximationCollinear(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs) {
  takeTime("Collinear Dipole Approximation");
  if (_dipoleInts.size() < 1)
    _dipoleInts = calculateDipoleIntegrals();
  preCalculateTransformations(orbitalPairs);
  // Transform the fock matrix.
  Eigen::MatrixXd f_MO = _occCoefficients->transpose() * *_fockMatrix * *_occCoefficients;
  // Loop over pairs:
  unsigned int pairIndex = 0;
  unsigned int nThreads = omp_get_max_threads();
  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(static)
  for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
    auto& pair = orbitalPairs[iPair];
    unsigned int i = pair->i;
    unsigned int j = pair->j;
    // Extract precalculated values.
    Eigen::MatrixXd muPAO_r_i = *_pao_r_i[i];
    Eigen::MatrixXd muPAO_r_j = *_pao_r_i[j];
    const Eigen::Vector3d r_ij = _iri.col(i) - _iri.col(j);
    Eigen::MatrixXd denominator = Eigen::MatrixXd::Zero(_eigenvalues_ii[i]->size(), _eigenvalues_ii[j]->size());
    denominator.colwise() += *_eigenvalues_ii[i];
    denominator.rowwise() += _eigenvalues_ii[j]->transpose().eval();
    denominator.array() -= f_MO(i, i) + f_MO(j, j);
    /*
     * Sum all up according to equation 17 J.Chem.Phys. 143, 034108 (2015)
     * Assume collinear orientation of the dipoles --> take the norm of all
     * (l|r|k) integral sets.
     *
     * The paper is not very clear here... I am actual not sure what they mean with
     * collinear dipoles in the equation.
     */
    muPAO_r_i = muPAO_r_i.rowwise().norm().eval();
    muPAO_r_j = muPAO_r_j.rowwise().norm().eval();
    // First term: (i|r|mu_PAO)(j|r|nu_PAO)
    // nx1 * 1xn = nxn
    Eigen::MatrixXd epsilon_mu_nu_i_j = muPAO_r_i * muPAO_r_j.transpose();
    Eigen::RowVector3d unit_r_ij_T = (r_ij.array() / r_ij.norm()).transpose().eval();
    // Second term: -3((i|r|mu_PAO)r_ij)((j|r|nu_PAO)r_ij
    //                               (  nx1 * 1x3       *     (nx1 * 1x3)^T              )= nxn
    epsilon_mu_nu_i_j -= 3 * ((muPAO_r_i * unit_r_ij_T) * (muPAO_r_j * unit_r_ij_T).transpose());
    // Square it and divide by eps_mu+eps_nu-F_ii-F_jj
    epsilon_mu_nu_i_j.array() *= epsilon_mu_nu_i_j.array() / denominator.array();
    // Pair energy is given as the sum of all matrix elements multiplied with -4/r_ij^6
    double pairEnergy = epsilon_mu_nu_i_j.sum();
    double r_ijSquared = r_ij.squaredNorm();
    pairEnergy *= -4 / (r_ijSquared * r_ijSquared * r_ijSquared);
    pair->dipoleCollinearPairEnergy = pairEnergy;
    ++pairIndex;
  } // for pair
  Eigen::setNbThreads(nThreads);
  timeTaken(2, "Collinear Dipole Approximation");
}

} /* namespace Serenity */
