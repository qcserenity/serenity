/**
 * @file   CombinedShellPair.cpp
 *
 * @date   Jul 30, 2020
 * @author Lars Hellmann
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "basis/CombinedShellPair.h"
#include "math/IntegerMaths.h"
#include "misc/SerenityError.h"

namespace Serenity {

CombinedShellPair::CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB, bool spherical)
  : Shell(generateExponents(shellA, shellB),
          generateContractions(shellA, shellB, shellA->getAngularMomentum() + shellB->getAngularMomentum()),
          shellA->getAngularMomentum() + shellB->getAngularMomentum(), spherical, checkCoords(shellA->O, shellB->O),
          checkElement(shellA->getElement(), shellB->getElement())),
    _shellA(shellA),
    _shellB(shellB){};

CombinedShellPair::CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB,
                                     unsigned int angularMomentum, bool spherical)
  : Shell(generateExponents(shellA, shellB), generateContractions(shellA, shellB, angularMomentum), angularMomentum,
          spherical, checkCoords(shellA->O, shellB->O), checkElement(shellA->getElement(), shellB->getElement())),
    _shellA(shellA),
    _shellB(shellB){};

CombinedShellPair::CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB,
                                     libint2::svector<double> expo, libint2::svector<double> contr,
                                     unsigned int angularMomentum, bool spherical)
  : Shell(expo, reverseNormalization(expo, contr, angularMomentum), angularMomentum, spherical,
          checkCoords(shellA->O, shellB->O), checkElement(shellA->getElement(), shellB->getElement())),
    _shellA(shellA),
    _shellB(shellB){};

libint2::svector<double> CombinedShellPair::generateExponents(std::shared_ptr<const Shell> shellA,
                                                              std::shared_ptr<const Shell> shellB) {
  libint2::svector<double> exponentsA = shellA->getExponents();
  libint2::svector<double> exponentsB = shellB->getExponents();

  libint2::svector<double> combinedExponents;
  if ((*shellA) != (*shellB)) {
    for (double i : exponentsA) {
      for (double j : exponentsB) {
        combinedExponents.push_back(i + j);
      }
    }
  }
  else {
    unsigned int nI = exponentsA.size();
    for (unsigned int i = 0; i < nI; i++) {
      for (unsigned int j = 0; j <= i; j++) {
        combinedExponents.push_back(exponentsA[i] + exponentsA[j]);
      }
    }
  }

  return combinedExponents;
}

libint2::svector<double> CombinedShellPair::generateContractions(std::shared_ptr<const Shell> shellA,
                                                                 std::shared_ptr<const Shell> shellB,
                                                                 unsigned int angularMomentum) {
  libint2::svector<double> contractionsA = shellA->getNormContractions();
  libint2::svector<double> contractionsB = shellB->getNormContractions();
  libint2::svector<double> exponentsA = shellA->getExponents();
  libint2::svector<double> exponentsB = shellB->getExponents();

  // The following reverses the normalization always performed by libint::shell() if it is constructed from scratch
  // and not by the copy constructor. [see libint::Shell::renorm()]
  double two_to_am = pow(2, angularMomentum);
  double sqrt_pi_cubed = sqrt(pow(PI, 3));
  unsigned int doubleFactorial = double_factorial(2 * angularMomentum - 1);
  unsigned int nI = contractionsA.size();
  unsigned int nJ = contractionsB.size();

  libint2::svector<double> newExponents = generateExponents(shellA, shellB);

  libint2::svector<double> newContractions;

  if ((*shellA) != (*shellB)) {
    for (unsigned int i = 0, ij = 0; i < nI; i++) {
      for (unsigned int j = 0; j < nJ; j++, ij++) {
        double two_alpha_to_am32 = pow((2 * newExponents[ij]), (angularMomentum + 1)) * sqrt(2 * newExponents[ij]);
        double normFactor = sqrt(two_to_am * two_alpha_to_am32 / (sqrt_pi_cubed * doubleFactorial));
        newContractions.push_back((contractionsA[i] * contractionsB[j]) / normFactor);
      }
    }
  }
  else {
    for (unsigned int i = 0, ij = 0; i < nI; i++) {
      for (unsigned int j = 0; j <= i; j++, ij++) {
        double two_alpha_to_am32 = pow((2 * newExponents[ij]), (angularMomentum + 1)) * sqrt(2 * newExponents[ij]);
        double normFactor = sqrt(two_to_am * two_alpha_to_am32 / (sqrt_pi_cubed * doubleFactorial));
        if (i == j) {
          newContractions.push_back((contractionsA[i] * contractionsB[j]) / normFactor);
        }
        else {
          newContractions.push_back((2 * contractionsA[i] * contractionsB[j]) / normFactor);
        }
      }
    }
  }

  return newContractions;
}

libint2::svector<double> CombinedShellPair::reverseNormalization(libint2::svector<double> expo, libint2::svector<double> contr,
                                                                 unsigned int angularMomentum) {
  // The following reverses the normalization always performed by libint::shell() if it is constructed from scratch
  // and not by the copy constructor. [see libint::Shell::renorm()]
  double two_to_am = pow(2, angularMomentum);
  double sqrt_pi_cubed = sqrt(pow(PI, 3));
  unsigned int doubleFactorial = double_factorial(2 * angularMomentum - 1);
  unsigned int nIJ = contr.size();

  libint2::svector<double> newContractions;

  for (unsigned int ij = 0; ij < nIJ; ij++) {
    double two_alpha_to_am32 = pow((2 * expo[ij]), (angularMomentum + 1)) * sqrt(2 * expo[ij]);
    double normFactor = sqrt(two_to_am * two_alpha_to_am32 / (sqrt_pi_cubed * doubleFactorial));
    newContractions.push_back(contr[ij] / normFactor);
  }
  return newContractions;
}

std::array<double, 3> CombinedShellPair::checkCoords(std::array<double, 3> a, std::array<double, 3> b) {
  if (a[0] != b[0] || a[1] != b[1] || a[2] != b[2])
    throw SerenityError("CombinedShellPair: Base shells need to have the same origin");
  return a;
}

std::string CombinedShellPair::checkElement(std::string a, std::string b) {
  if (a != b)
    throw SerenityError("CombinedShellPair: Base shells have to correspond to the same element");
  return a;
}

} // namespace Serenity
