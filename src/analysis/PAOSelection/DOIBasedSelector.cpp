/**
 * @file DOIBasedSelector.cpp
 *
 * @date Apr 3, 2019
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
#include "analysis/PAOSelection/DOIBasedSelector.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"               //Atom<-->Shell map etc.
#include "data/PAOController.h"                              //PAOController definition.
#include "data/grid/DifferentialOverlapIntegralCalculator.h" //DOI calculator.
#include "misc/SystemSplittingTools.h"                       //Presceening matrices.

namespace Serenity {

std::shared_ptr<Eigen::SparseMatrix<int>> DOIBasedSelector::selectPAOs() {
  const auto& occC = *_occupiedCoefficients;
  const auto& paos = _paoController->getAllPAOs();
  unsigned int nPAOs = _paoController->getNPAOs();

  Eigen::SparseMatrix<int> basisFunctionsToPAOMap =
      SystemSplittingTools<Options::SCF_MODES::RESTRICTED>::reduceMatrixToMullikenNetPopulationMap(paos, _mnpPreThreshold);
  Eigen::SparseMatrix<int> basisFunctionsToOccMap =
      SystemSplittingTools<Options::SCF_MODES::RESTRICTED>::reduceMatrixToMullikenNetPopulationMap(occC, _mnpPreThreshold);
  Eigen::MatrixXd dois;
  DifferentialOverlapIntegralCalculator::calculateDOI(paos, occC, basisFunctionsToPAOMap, basisFunctionsToOccMap,
                                                      _basOnGridController, dois);
  auto idx = _atomCenteredBasisController->getBasisIndices();
  std::vector<Eigen::Triplet<int>> tripletList;
  for (unsigned int iOrb = 0; iOrb < occC.cols(); ++iOrb) {
    for (const auto& index : idx) {
      unsigned int nAtomPAOs = index.second - index.first;
      unsigned int nImportant = (dois.block(index.first, iOrb, nAtomPAOs, 1).array().abs() >= _doiThresholds(iOrb)).count();
      if (nImportant != 0) {
        for (unsigned int mu = index.first; mu < index.second; ++mu) {
          tripletList.push_back(Eigen::Triplet<int>(mu, iOrb, 1));
        } // for mu
      }   // if nImportant!=0
    }     // for index (atom)
  }       // for iOrb
  auto paoToOccMap = std::make_shared<Eigen::SparseMatrix<int>>(nPAOs, occC.cols());
  paoToOccMap->setFromTriplets(tripletList.begin(), tripletList.end());
  return paoToOccMap;
}

} /* namespace Serenity */
