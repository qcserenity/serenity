/**
 * @file   BasisController__TEST_SUPPLY.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   23. August 2015, 18:25
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
#include "testsupply/BasisController__TEST_SUPPLY.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/Shell.h"

namespace Serenity {

std::map<TEST_BASIS_CONTROLLERS, std::shared_ptr<AtomCenteredBasisController>> BasisController__TEST_SUPPLY::_testBasisControllers;

void BasisController__TEST_SUPPLY::prepare(TEST_BASIS_CONTROLLERS kind) {
  auto thisBasisController = std::shared_ptr<AtomCenteredBasisController>(new AtomCenteredBasisController());
  switch (kind) {
      /***********/
      /* MINIMAL */
      /***********/
    case TEST_BASIS_CONTROLLERS::MINIMAL: {
      thisBasisController->_basis = std::unique_ptr<Basis>(new Basis());
      thisBasisController->_basis->push_back(std::shared_ptr<Shell>(new Shell({2.0}, {0.5}, 0, false, {{0.0, 0.0, -1.0}})));
      thisBasisController->_basis->push_back(std::shared_ptr<Shell>(new Shell({1.0}, {0.8}, 0, false, {{0.0, 0.0, +1.0}})));
      /*
       * Fill all the other member variables
       */
      thisBasisController->_nBasisFunctions = 2;
      thisBasisController->_nBasisFunctionsCart = thisBasisController->_nBasisFunctions;
      thisBasisController->_nBasisFunctionsSpherical = thisBasisController->_nBasisFunctions;
      thisBasisController->_reducedIndex = {0, 1};
      thisBasisController->_extendedIndex = thisBasisController->_reducedIndex;
      thisBasisController->_extendedIndexCart = thisBasisController->_extendedIndex;
      thisBasisController->_extendedIndexSpherical = thisBasisController->_extendedIndex;
      thisBasisController->_maxAngularMomentum = 0;
      thisBasisController->_geometry.reset();
      thisBasisController->_basisIndicesOfAtom = {{0, 1}, {1, 2}};
      thisBasisController->_basisIndicesRedOfAtom = thisBasisController->_basisIndicesOfAtom;
      break;
    }
      /***************/
      /* SMALL_MIXED */
      /***************/
    case TEST_BASIS_CONTROLLERS::SMALL_MIXED: {
      thisBasisController->_basis = std::unique_ptr<Basis>(new Basis());
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({1.0, 2.0}, {1.0, 0.1}, 0, false, {{0.0, 0.0, 0.0}})));
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({1.0, 2.0}, {1.0, 0.1}, 1, false, {{1.0, 0.0, 0.0}})));
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({0.5, 1.2}, {0.8, 0.4}, 2, false, {{-1.0, 2.5, 1.0}})));
      /*
       * Fill all the other member variables
       */
      thisBasisController->_nBasisFunctions = 10;
      thisBasisController->_nBasisFunctionsCart = thisBasisController->_nBasisFunctions;
      thisBasisController->_nBasisFunctionsSpherical = 9;
      thisBasisController->_reducedIndex = {0, 1, 1, 1, 2, 2, 2, 2, 2, 2};
      thisBasisController->_extendedIndex = {0, 1, 4};
      thisBasisController->_extendedIndexCart = thisBasisController->_extendedIndex;
      thisBasisController->_extendedIndexSpherical = thisBasisController->_extendedIndex;
      thisBasisController->_maxAngularMomentum = 2;
      thisBasisController->_geometry.reset();
      thisBasisController->_basisIndicesOfAtom = {{0, 1}, {1, 4}, {4, 10}};
      thisBasisController->_basisIndicesRedOfAtom = {{0, 1}, {1, 2}, {2, 3}};
      break;
    }
      /*****************/
      /* ASYMMETRIC_CH3 */
      /*****************/
    case TEST_BASIS_CONTROLLERS::ASYMMETRIC_CH3: {
      /*
       * Basis
       */
      thisBasisController->_basis = std::unique_ptr<Basis>(new Basis());
      // H1
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({18.7311370, 2.8253937, 0.6401217}, {0.03349460, 0.23472695, 0.81375733}, 0,
                                           false, {{-1.01199, 0.28861, 1.33493}})));
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({0.1612778}, {1.0000000}, 0, false, {{-1.01199, 0.28861, 1.33493}})));
      // H2
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({18.7311370, 2.8253937, 0.6401217}, {0.03349460, 0.23472695, 0.81375733}, 0,
                                           false, {{-1.05807, 0.62810, 0.08722}})));
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({0.1612778}, {1.0000000}, 0, false, {{-1.05807, 0.62810, 0.08722}})));
      // H3
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({18.7311370, 2.8253937, 0.6401217}, {0.03349460, 0.23472695, 0.81375733}, 0,
                                           false, {{-1.15616, 2.73228, 0.30975}})));
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({0.1612778}, {1.0000000}, 0, false, {{-1.15616, 2.73228, 0.30975}})));
      // C1
      thisBasisController->_basis->push_back(std::shared_ptr<Shell>(new Shell(
          {3047.5249000, 457.3695100, 103.9486900, 29.2101550, 9.2866630, 3.1639270},
          {0.0018347, 0.0140373, 0.0688426, 0.2321844, 0.4679413, 0.3623120}, 0, false, {{-0.57388, 1.74108, 0.00000}})));
      thisBasisController->_basis->push_back(std::shared_ptr<Shell>(new Shell(
          {7.8682724, 1.8812885, 0.5442493}, {-0.1193324, -0.1608542, 1.1434564}, 0, false, {{-0.57388, 1.74108, 0.00000}})));
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({0.1687144}, {1.0000000}, 0, false, {{-0.57388, 1.74108, 0.00000}})));
      thisBasisController->_basis->push_back(std::shared_ptr<Shell>(new Shell(
          {7.8682724, 1.8812885, 0.5442493}, {0.0689991, 0.3164240, 0.7443083}, 1, false, {{-0.57388, 1.74108, 0.00000}})));
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({0.1687144}, {1.0000000}, 1, false, {{-0.57388, 1.74108, 0.00000}})));
      thisBasisController->_basis->push_back(
          std::shared_ptr<Shell>(new Shell({0.8000000}, {1.0000000}, 2, false, {{-0.57388, 1.74108, 0.00000}})));
      /*
       * Fill all the other member variables
       */
      thisBasisController->_nBasisFunctions = 21;
      thisBasisController->_nBasisFunctionsCart = thisBasisController->_nBasisFunctions;
      thisBasisController->_nBasisFunctionsSpherical = 20;
      thisBasisController->_reducedIndex = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 11, 11};
      thisBasisController->_extendedIndex = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15};
      thisBasisController->_extendedIndexCart = thisBasisController->_extendedIndex;
      thisBasisController->_extendedIndexSpherical = thisBasisController->_extendedIndex;
      thisBasisController->_maxAngularMomentum = 2;
      thisBasisController->_geometry.reset();
      thisBasisController->_basisIndicesOfAtom = {{0, 2}, {2, 4}, {4, 21}};
      thisBasisController->_basisIndicesRedOfAtom = {{0, 2}, {2, 4}, {4, 12}};
      break;
    }
  }
  _testBasisControllers[kind] = thisBasisController;
}

} /* namespace Serenity */
