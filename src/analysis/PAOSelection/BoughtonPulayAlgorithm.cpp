/**
 * @file BoughtonPulayAlgorithm.cpp
 *
 * @date Dec 10, 2018
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
#include "analysis/PAOSelection/BoughtonPulayAlgorithm.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"       //Atom<-->Shell map etc.
#include "data/matrices/MatrixInBasis.h"             //Coefficient matrix operations.
#include "geometry/Atom.h"                           //Definition of an atom.
#include "integrals/OneElectronIntegralController.h" //Overlap integrals.
#include "io/FormattedOutputStream.h"                //Filtered output streams.
#include "misc/SystemSplittingTools.h"               //Matrix block selection.

namespace Serenity {

BoughtonPulayAlgorithm::BoughtonPulayAlgorithm(std::shared_ptr<OneElectronIntegralController> oneElectronIntegralController,
                                               std::shared_ptr<AtomCenteredBasisController> atomCenteredBasisController,
                                               std::shared_ptr<Eigen::MatrixXd> mullikenGrossCharges,
                                               std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients, double completeness)
  : _oneElectronIntegralController(oneElectronIntegralController),
    _atomCenteredBasisController(atomCenteredBasisController),
    _mullikenGrossCharges(mullikenGrossCharges), // nAtoms x nOcc
    _occupiedCoefficients(occupiedCoefficients), // nBasis x nOcc
    _completeness(completeness) {
  assert(_mullikenGrossCharges->cols() == _occupiedCoefficients->cols());
}

std::shared_ptr<Eigen::SparseMatrix<int>> BoughtonPulayAlgorithm::selectPAOs() {
  // Loop over all orbitals
  unsigned int nPAOs = _atomCenteredBasisController->getNBasisFunctions();
  auto paoToOccMap = std::make_shared<Eigen::SparseMatrix<int>>(nPAOs, _occupiedCoefficients->cols());
  std::vector<Eigen::Triplet<int>> tripletList;
  auto idx = _atomCenteredBasisController->getBasisIndices();
  for (unsigned int iOrb = 0; iOrb < _occupiedCoefficients->cols(); ++iOrb) {
    // Assign atoms to each orbital
    auto atomIndices = this->assignAtoms(_mullikenGrossCharges->col(iOrb), _occupiedCoefficients->col(iOrb));
    // Get the basis function indices of the atoms and add them to the triplet list
    for (const auto& atomIndex : atomIndices) {
      for (unsigned int mu = idx[atomIndex].first; mu < idx[atomIndex].second; ++mu) {
        tripletList.push_back(Eigen::Triplet<int>(mu, iOrb, 1));
      } // for mu
    }   // for atomIndex
  }     // for iOrb
  paoToOccMap->setFromTriplets(tripletList.begin(), tripletList.end());
  return paoToOccMap;
}

std::vector<unsigned int> BoughtonPulayAlgorithm::assignAtoms(Eigen::VectorXd mullikenGrossCharges, const Eigen::VectorXd& c_i) {
  /*
   * Select the atoms with the highest Mulliken charges up to
   * an accumulative charge of 0.9. (JCP 111, 5691 (1999))
   * Then test the "completeness" of the AO space for the given
   * LMO.
   */
  double accuMullikenCharge = 0.0;
  double minimumCharge = 0.9;
  std::vector<unsigned int> allIndices;
  for (unsigned int index = 0; index < mullikenGrossCharges.size(); ++index)
    allIndices.push_back(index);
  std::vector<unsigned int> atomIndices;
  double completeness = 0.0;

  while (true) {
    unsigned int largestIndex = 0;
    double largestCharge = mullikenGrossCharges.maxCoeff(&largestIndex);
    atomIndices.push_back(largestIndex);
    accuMullikenCharge += largestCharge;
    mullikenGrossCharges[largestIndex] = -1000.0;
    // Check the "completeness"
    if (accuMullikenCharge >= minimumCharge || atomIndices.size() == allIndices.size()) {
      auto s_A = SystemSplittingTools<Options::SCF_MODES::RESTRICTED>::getMatrixBlock(
          _oneElectronIntegralController->getOverlapIntegrals(), atomIndices, atomIndices,
          _atomCenteredBasisController->getBasisIndices());
      auto s_nu = SystemSplittingTools<Options::SCF_MODES::RESTRICTED>::getMatrixBlock(
          _oneElectronIntegralController->getOverlapIntegrals(), atomIndices, allIndices,
          _atomCenteredBasisController->getBasisIndices());
      // calculate invers overlap matrix in the selected basis set.
      auto s_Ainv = s_A.completeOrthogonalDecomposition().pseudoInverse();
      // calculate new coefficients a_i
      auto a_i = s_Ainv * s_nu * c_i;
      completeness = 1.0 - a_i.transpose() * s_A * a_i;
      if (std::fabs(completeness) <= _completeness || atomIndices.size() == allIndices.size())
        break;
    }
  } // while accuMullikenCharge < minimumCharge
  OutputControl::vOut << "AO selection completness " << completeness << " Atoms: ";
  for (const auto& atom : atomIndices)
    OutputControl::vOut << atom << " ";
  OutputControl::vOut << std::endl;
  return atomIndices;
}

} /* namespace Serenity */
