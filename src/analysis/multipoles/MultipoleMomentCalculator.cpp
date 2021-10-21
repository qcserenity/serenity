/**
 * @file MultipoleMomentCalculator.cpp
 *
 * @date Dec 8, 2015
 * @author David Schnieders
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
#include "analysis/multipoles/MultipoleMomentCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"                         //Loop shells.
#include "basis/BasisController.h"               //getBasis().
#include "data/ElectronicStructure.h"            //getDensityMatrix().
#include "data/matrices/DensityMatrix.h"         //.total().
#include "geometry/Atom.h"                       //getX() etc.
#include "geometry/Geometry.h"                   //getAtoms()
#include "integrals/wrappers/Libint.h"           //compute integrals.
#include "math/IntegerMaths.h"                   //factorial().
#include "settings/ElectronicStructureOptions.h" //RESTRICTED/UNRESTRICTED
#include "system/SystemController.h"             //getGeometry()
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
std::vector<std::vector<double>>
MultipoleMomentCalculator::calculateMultipoleMoment(std::shared_ptr<SystemController> system, unsigned int highestOrder) {
  // Geometry
  auto geometry = system->getGeometry();

  // Density Matrix
  auto densMatrix(system->getElectronicStructure<SCFMode>()->getDensityMatrix().total());

  // Basis
  auto basisController = system->getBasisController();
  auto& basis = basisController->getBasis();

  // prepare vector
  std::vector<std::vector<double>> multipoleMoment;
  for (unsigned int multipole = 0; multipole < highestOrder; multipole++) {
    std::vector<double> subVec(factorial(3 + multipole) / (factorial(2) * factorial(multipole + 1)), 0.0);
    multipoleMoment.push_back(subVec);
  }

  auto& libint = Libint::getInstance();
  if (highestOrder == 1) {
    libint.initialize(LIBINT_OPERATOR::emultipole1, 0, 2);
  }
  else if (highestOrder == 2) {
    libint.initialize(LIBINT_OPERATOR::emultipole2, 0, 2);
  }

  assert(highestOrder <= 2 && highestOrder >= 1 && "Only orders 1-2 supported for calculation of multipole moments!");

  // Subtract electronic part

  for (unsigned int i = 0; i < basis.size(); i++) {
    for (unsigned int j = 0; j < basis.size(); j++) {
      auto shellA = basis[i]->getNContracted();
      auto shellB = basis[j]->getNContracted();
      Eigen::MatrixXd multiPoleInts;
      if (highestOrder == 1) {
        libint.compute(LIBINT_OPERATOR::emultipole1, 0, *basis[i], *basis[j], multiPoleInts);
      }
      else if (highestOrder == 2) {
        libint.compute(LIBINT_OPERATOR::emultipole2, 0, *basis[i], *basis[j], multiPoleInts);
      }
      /*
       * set vector size: libint will return a matrix containing
       * <mu|nu>,<mu|x|nu>,<mu|y|nu>,<mu|z|nu>,<mu|xx|nu>,<mu|xy|nu>,<mu|xz|nu>,
       *                                       <mu|yy|nu>,<mu|yz|nu>,<mu|zz|nu>
       * with nBasTot integrals each. Therefore, depending on the requested order,
       * we will therefore calculate (multiPoleInts.size()/nBasTot) - 1(overlap)
       * multipole elements.
       */
      for (unsigned int k = 0; k < shellA; k++) {
        auto mu = basisController->extendedIndex(i) + k;
        for (unsigned int l = 0; l < shellB; l++) {
          auto nu = basisController->extendedIndex(j) + l;
          for (unsigned int multipole = 0; multipole < highestOrder; multipole++) {
            for (unsigned int element = 0; element < (multipole + 1) * 3; element++) {
              multipoleMoment[multipole][element] -=
                  densMatrix(mu, nu) * multiPoleInts((shellB * k + l), (multipole * 3 + 1 + element));
            }
          }
        }
      }
    }
  };
  if (highestOrder == 1) {
    libint.finalize(LIBINT_OPERATOR::emultipole1, 0, 2);
  }
  else if (highestOrder == 2) {
    libint.finalize(LIBINT_OPERATOR::emultipole2, 0, 2);
  }
  // Add nuclear part
  for (const auto& atom : geometry->getAtoms()) {
    const double charge = atom->getEffectiveCharge();
    /*
     * since libint dictates a specific order, I was not able to find a clever way to calculate
     * the current multipole element from coord1 and coord2. Therefore, an additional increment
     * variable is introduced. Initialized to 2 in order to skip the dipole part of the vector.
     */
    unsigned int quadrupoleElement = 0;
    std::vector<double> coords(3, 0.0);
    coords[0] = atom->getX();
    coords[1] = atom->getY();
    coords[2] = atom->getZ();
    /*
     * calculate dipole part, i.e. q*x, q*y, q*z
     */
    for (unsigned int coord1 = 0; coord1 < 3; coord1++) {
      multipoleMoment[0][coord1] += charge * coords[coord1];
      if (highestOrder > 1) {
        /*
         * calculate quadrupole part, i.e. q*xx, q*xy, q*xz, ...
         */
        for (unsigned int coord2 = coord1; coord2 < 3; coord2++) {
          multipoleMoment[1][quadrupoleElement] += charge * coords[coord1] * coords[coord2];
          quadrupoleElement++;
        }
      }
    }
  };

  return multipoleMoment;
}

template std::vector<std::vector<double>>
MultipoleMomentCalculator::calculateMultipoleMoment<Options::SCF_MODES::RESTRICTED>(std::shared_ptr<SystemController> system,
                                                                                    unsigned int highestOrder);
template std::vector<std::vector<double>>
MultipoleMomentCalculator::calculateMultipoleMoment<Options::SCF_MODES::UNRESTRICTED>(std::shared_ptr<SystemController> system,
                                                                                      unsigned int highestOrder);

} /* namespace Serenity */
