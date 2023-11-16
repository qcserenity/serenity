/**
 * @file BasisFunctionMapper.cpp
 *
 * @date 10 Aug 2019
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
#include "basis/BasisFunctionMapper.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"                 //Access shell info..
#include "basis/BasisController.h"       //Basis construction from shell lists.
#include "basis/CustomBasisController.h" //Basis construction from shell lists.
#include "basis/Shell.h"                 //Basis construction from shell lists and shell comparison.
#include "io/FormattedOutputStream.h"    //Filtered output.
#include "io/IOOptions.h"                //Print debug info.

namespace Serenity {

BasisFunctionMapper::BasisFunctionMapper(std::shared_ptr<BasisController> basisControllerA)
  : _basisControllerA(basisControllerA) {
}

std::shared_ptr<BasisController> BasisFunctionMapper::getDifferentialBasis(std::shared_ptr<BasisController> basisControllerB) {
  std::vector<std::shared_ptr<const Shell>> differentialBasis;
  unsigned int nShellsA = _basisControllerA->getBasis().size();

  for (const auto& shellB : basisControllerB->getBasis()) {
    unsigned int shellIndexInA = getShellIndex(*shellB);
    if (shellIndexInA >= nShellsA)
      differentialBasis.push_back(shellB);
  }
  if (differentialBasis.size() == 0) {
    if (iOOptions.printDebugInfos)
      OutputControl::dOut << "Remark: Basis Function Mapper -- Basis sets are identical!" << std::endl;
    return nullptr;
  }
  auto customBasisController = std::make_shared<CustomBasisController>(differentialBasis, "DifferentialBasis");
  return customBasisController;
}

std::shared_ptr<BasisController> BasisFunctionMapper::getCombinedBasis(std::shared_ptr<BasisController> basisControllerB) {
  auto differentialBasisController = this->getDifferentialBasis(basisControllerB);
  if (!differentialBasisController)
    return _basisControllerA;
  auto differentialBasis = differentialBasisController->getBasis();
  Basis combinedBas = _basisControllerA->getBasis();
  combinedBas.insert(combinedBas.end(), differentialBasis.begin(), differentialBasis.end());
  auto customBasisController = std::make_shared<CustomBasisController>(combinedBas, "CombinedBasis");
  return customBasisController;
}

std::shared_ptr<Eigen::SparseMatrix<double>>
BasisFunctionMapper::getSparseProjection(std::shared_ptr<BasisController> basisControllerB) {
  const auto basisB = basisControllerB->getBasis();
  unsigned int nShellsA = _basisControllerA->getBasis().size();
  std::vector<Eigen::Triplet<double>> tripletList;
  unsigned int row = 0;
  /*
   * The resulting matrix should have the property that multiplying it
   * to a matrix of type basisA x something extracts the rows mapped to
   * the basis B.
   *                      1 2 3  4  5 6
   * Example for basis B: B B A1 A3 B A4 shells
   *                      A1 A2 A3 A4 A5
   *                  B   0  0  0  0  0
   *                  B   0  0  0  0  0
   *                  A1  1  0  0  0  0
   *                  A3  0  0  1  0  0
   *                  B   0  0  0  0  0
   *                  A4  0  0  0  1  0
   * 1 corresponds to a unit matrix over the contracted functions and zero
   * is a matrix containing only zeros.
   */
  for (const auto& shellB : basisB) {
    unsigned int shellStart = getShellIndex(*shellB);
    unsigned int nContractedB = shellB->getNContracted();
    if (shellStart < nShellsA) {
      unsigned int colStart = _basisControllerA->extendedIndex(shellStart);
      for (unsigned int col = colStart; col < colStart + nContractedB; ++col) {
        tripletList.push_back(Eigen::Triplet<double>(row, col, 1.0));
        ++row;
      }
    }
    else {
      row += nContractedB;
    }
  } // for shellB
  assert(row == basisControllerB->getNBasisFunctions());
  Eigen::SparseMatrix<double> sparseProjection(basisControllerB->getNBasisFunctions(),
                                               _basisControllerA->getNBasisFunctions());
  sparseProjection.setFromTriplets(tripletList.begin(), tripletList.end());
  for (unsigned int iCol = 0; iCol < sparseProjection.cols(); ++iCol)
    assert(sparseProjection.col(iCol).sum() <= 1);
  return std::make_shared<Eigen::SparseMatrix<double>>(sparseProjection);
}

unsigned int BasisFunctionMapper::getShellIndex(const Shell& other) {
  auto basisA = _basisControllerA->getBasis();
  unsigned int iShellA;
  for (iShellA = 0; iShellA < basisA.size(); ++iShellA) {
    if (*basisA[iShellA] == other)
      return iShellA;
  }
  return iShellA + 1;
}

} /* namespace Serenity */
