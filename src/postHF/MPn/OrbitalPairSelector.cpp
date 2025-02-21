/**
 * @file OrbitalPairSelector.cpp
 *
 * @date Apr 2, 2019
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
#include "postHF/MPn/OrbitalPairSelector.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"                           //Basis function values on the grid
#include "data/OrbitalController.h"                          //Coefficient matrix
#include "data/OrbitalPair.h"                                //Orbital pair definition
#include "data/PAOController.h"                              //PAOs
#include "data/grid/BasisFunctionOnGridControllerFactory.h"  //Basis function values on the grid
#include "data/grid/DifferentialOverlapIntegralCalculator.h" //DOI calculation
#include "data/matrices/CoefficientMatrix.h"                 //Coefficient matrix
#include "integrals/OneElectronIntegralController.h"         //Overlap matrix
#include "io/FormattedOutputStream.h"                        //Filtered output streams.
#include "misc/SystemSplittingTools.h"                       //MnP prescreening for DOI
#include "misc/Timing.h"                                     //Timings
#include "postHF/MPn/DipoleApproximationToPairEnergies.h"    //Dipole approximation
#include "system/SystemController.h"                         //SystemController
/* Include Std and External Headers */
#include <Eigen/SparseCore> //Sparse maps for prescreening.

namespace Serenity {

void OrbitalPairSelector::calculateDOIs(std::vector<std::shared_ptr<OrbitalPair>> initialPairs, double mnpPreThreshold) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  unsigned int nOcc = _systemController->getNOccupiedOrbitals<scfMode>();
  unsigned int nBasisFunctions = _systemController->getBasisController()->getNBasisFunctions();
  const CoefficientMatrix<scfMode>& coefficients =
      _systemController->getActiveOrbitalController<scfMode>()->getCoefficients();
  std::vector<bool> map_i(nOcc, false);
  std::vector<bool> map_j(nOcc, false);
  Eigen::MatrixXd c_x = Eigen::MatrixXd::Zero(nBasisFunctions, nOcc);
  Eigen::MatrixXd c_y = Eigen::MatrixXd::Zero(nBasisFunctions, nOcc);
  for (const auto& pair : initialPairs) {
    unsigned int i = pair->i;
    unsigned int j = pair->j;
    if (!map_i[i])
      c_x.col(i) = coefficients.col(i).eval();
    if (!map_j[j])
      c_y.col(j) = coefficients.col(j).eval();
    map_i[i] = true;
    map_j[j] = true;
  }
  Eigen::SparseMatrix<int> basisFunctionsToXMap =
      SystemSplittingTools<scfMode>::reduceMatrixToMullikenNetPopulationMap(c_x, mnpPreThreshold);
  Eigen::SparseMatrix<int> basisFunctionsToYMap =
      SystemSplittingTools<scfMode>::reduceMatrixToMullikenNetPopulationMap(c_y, mnpPreThreshold);
  auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _systemController->getSettings(), _systemController->getBasisController(), _systemController->getGridController());
  Eigen::MatrixXd dois;
  DifferentialOverlapIntegralCalculator::calculateDOI(c_x, c_y, basisFunctionsToXMap, basisFunctionsToYMap,
                                                      basFuncOnGridController, dois);
  // Pass DOIs to pairs.
  for (const auto& pair : initialPairs) {
    pair->doi = dois(pair->i, pair->j);
  }
}

void OrbitalPairSelector::calculateDipoleApproximation(std::vector<std::shared_ptr<OrbitalPair>> initialPairs,
                                                       double paoOrthoThreshold,
                                                       std::shared_ptr<FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                                                       std::shared_ptr<Eigen::SparseMatrix<int>> paoToOccupiedOrbitalMap) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  unsigned int nOcc = _systemController->getNOccupiedOrbitals<scfMode>();
  auto S = _systemController->getOneElectronIntegralController()->getOverlapIntegrals();
  auto s_ptr = std::make_shared<MatrixInBasis<scfMode>>(S);
  const CoefficientMatrix<scfMode>& coefficients =
      _systemController->getActiveOrbitalController<scfMode>()->getCoefficients();
  auto cOcc = std::make_shared<Eigen::MatrixXd>(coefficients.leftCols(nOcc).eval());
  DipoleApproximationToPairEnergies approxCalculator(_systemController->getBasisController(), s_ptr, cOcc,
                                                     _paoController, f, paoToOccupiedOrbitalMap, paoOrthoThreshold);
  approxCalculator.calculateDipoleApproximationCollinear(initialPairs);
  std::vector<std::shared_ptr<OrbitalPair>> veryDistantPairs;
  for (auto& pair : initialPairs) {
    if (std::fabs(pair->dipoleCollinearPairEnergy) < pair->getCollinearDipolePairThreshold())
      veryDistantPairs.push_back(pair);
  } // for pair
  approxCalculator.calculateDipoleApproximation(veryDistantPairs);
}

std::pair<std::vector<std::shared_ptr<OrbitalPair>>, std::vector<std::shared_ptr<OrbitalPair>>>
OrbitalPairSelector::selectOrbitalPairs(std::vector<std::shared_ptr<OrbitalPair>> initialPairs, double mnpPreThreshold,
                                        double paoOrthoThreshold,
                                        std::shared_ptr<FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                                        std::shared_ptr<Eigen::SparseMatrix<int>> paoToOccupiedOrbitalMap,
                                        double doiThreshold) {
  takeTime("Pair Prescreening");
  assert(_systemController->getLastSCFMode() == Options::SCF_MODES::RESTRICTED &&
         "Local MP2 is only implemented for RESTRICTED calculations.");
  // Results are stored directly in the orbital pairs.
  calculateDOIs(initialPairs, mnpPreThreshold);
  std::vector<std::shared_ptr<OrbitalPair>> dipoleCheckPairs;
  for (auto& pair : initialPairs) {
    if (std::fabs(pair->doi) < doiThreshold)
      dipoleCheckPairs.push_back(pair);
  } // for pair
  calculateDipoleApproximation(dipoleCheckPairs, paoOrthoThreshold, f, paoToOccupiedOrbitalMap);

  std::vector<std::shared_ptr<OrbitalPair>> significantPairs;
  std::vector<std::shared_ptr<OrbitalPair>> inSignificantPairs;
  for (const auto& pair : initialPairs) {
    bool significantDoi = false;
    bool significantCollinearPair = false;
    if (std::fabs(pair->doi >= doiThreshold))
      significantDoi = true;
    if (std::fabs(pair->dipoleCollinearPairEnergy) >= pair->getCollinearDipolePairThreshold())
      significantCollinearPair = true;
    if (!significantDoi && !significantCollinearPair) {
      pair->type = OrbitalPairTypes::VERY_DISTANT;
      inSignificantPairs.push_back(pair);
    } // if !significantDoi && !significantCollinearPair
    else {
      significantPairs.push_back(pair);
    } // if else !significantDoi && !significantCollinearPair
  }   // for pair

  OutputControl::nOut << " Orbital Pair Prescreening:" << std::endl;
  OutputControl::nOut << "  Significant pairs: " << std::to_string(significantPairs.size()) << std::endl;
  OutputControl::nOut << "  Distant pairs:     " + std::to_string(inSignificantPairs.size()) << std::endl;
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
  timeTaken(2, "Pair Prescreening");
  return std::make_pair(significantPairs, inSignificantPairs);
}

} /* namespace Serenity */
