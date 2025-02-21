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

/* Include Class Header*/
#include "basis/CombinedShellPair.h"
/* Include Serenity Internal Headers */
#include "math/IntegerMaths.h"
#include "misc/SerenityError.h"

namespace Serenity {

CombinedShellPair::CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB, bool spherical)
  : CombinedShellPair(shellA, shellB, generateExponents(shellA, shellB), generateContractionsNoNorm(shellA, shellB),
                      (shellA->getAngularMomentum() + shellB->getAngularMomentum()), spherical){};

CombinedShellPair::CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB,
                                     unsigned int angularMomentum, bool spherical)
  : CombinedShellPair(shellA, shellB, generateExponents(shellA, shellB), generateContractionsNoNorm(shellA, shellB),
                      angularMomentum, spherical){};

CombinedShellPair::CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB,
                                     std::vector<double> expo, std::vector<double> contr, unsigned int angularMomentum,
                                     bool spherical)
  : CombinedShellPair(shellA, shellB, expo, reverseNormalization(expo, contr, angularMomentum), angularMomentum, spherical,
                      checkCoords(shellA->O, shellB->O), checkElement(shellA->getElement(), shellB->getElement())){};

CombinedShellPair::CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB,
                                     std::vector<double> expo, std::vector<double> contr, unsigned int angularMomentum,
                                     bool spherical, std::array<double, 3> coords, std::string element)
  : Shell(libint2::svector<double>(expo.begin(), expo.end()), libint2::svector<double>(expo.begin(), expo.end()),
          libint2::svector<double>(contr.begin(), contr.end()), angularMomentum, spherical, coords, element),
    _shellA(shellA),
    _shellB(shellB){};

std::vector<double> CombinedShellPair::generateExponents(std::shared_ptr<const Shell> shellA,
                                                         std::shared_ptr<const Shell> shellB) {
  auto libExpoA = shellA->getExponents();
  auto libExpoB = shellB->getExponents();
  std::vector<double> exponentsA(libExpoA.begin(), libExpoA.end());
  std::vector<double> exponentsB(libExpoB.begin(), libExpoB.end());

  std::vector<double> tmpVec;

  if ((*shellA) != (*shellB)) {
    for (double i : exponentsA) {
      for (double j : exponentsB) {
        tmpVec.push_back(i + j);
      }
    }
  }
  else {
    unsigned int nI = exponentsA.size();
    for (unsigned int i = 0; i < nI; i++) {
      for (unsigned int j = 0; j <= i; j++) {
        tmpVec.push_back(exponentsA[i] + exponentsB[j]);
      }
    }
  }

  return tmpVec;
}

std::vector<double> CombinedShellPair::generateContractions(std::shared_ptr<const Shell> shellA,
                                                            std::shared_ptr<const Shell> shellB, unsigned int angularMomentum) {
  auto libContrA = shellA->getNormContractions();
  auto libContrB = shellB->getNormContractions();
  std::vector<double> contractionsA(libContrA.begin(), libContrA.end());
  std::vector<double> contractionsB(libContrB.begin(), libContrB.end());
  auto libExpoA = shellA->getExponents();
  auto libExpoB = shellB->getExponents();
  std::vector<double> exponentsA(libExpoA.begin(), libExpoA.end());
  std::vector<double> exponentsB(libExpoB.begin(), libExpoB.end());

  // The following reverses the normalization always performed by libint::shell() if it is constructed from scratch
  // and not by the copy constructor. [see libint::Shell::renorm()]
  double two_to_am = pow(2, angularMomentum);
  double sqrt_pi_cubed = sqrt(pow(PI, 3));
  unsigned int doubleFactorial = double_factorial(2 * angularMomentum - 1);
  unsigned int nI = contractionsA.size();
  unsigned int nJ = contractionsB.size();

  std::vector<double> newExponents = generateExponents(shellA, shellB);

  std::vector<double> newContractions;

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

std::vector<double> CombinedShellPair::generateContractionsNoNorm(std::shared_ptr<const Shell> shellA,
                                                                  std::shared_ptr<const Shell> shellB) {
  auto libContrA = shellA->getNormContractions();
  auto libContrB = shellB->getNormContractions();
  std::vector<double> contractionsA(libContrA.begin(), libContrA.end());
  std::vector<double> contractionsB(libContrB.begin(), libContrB.end());

  unsigned int nI = contractionsA.size();
  unsigned int nJ = contractionsB.size();

  unsigned int size = ((*shellA) != (*shellB)) ? nI * nJ : nI * (nI + 1) / 2;
  std::vector<double> newContractions(size);

  if ((*shellA) != (*shellB)) {
    for (unsigned int i = 0, ij = 0; i < nI; i++) {
      for (unsigned int j = 0; j < nJ; j++, ij++) {
        newContractions[ij] = contractionsA[i] * contractionsB[j];
      }
    }
  }
  else {
    for (unsigned int i = 0, ij = 0; i < nI; i++) {
      for (unsigned int j = 0; j <= i; j++, ij++) {
        if (i == j) {
          newContractions[ij] = contractionsA[i] * contractionsB[j];
        }
        else {
          newContractions[ij] = 2 * contractionsA[i] * contractionsB[j];
        }
      }
    }
  }
  return newContractions;
}

std::vector<double> CombinedShellPair::reverseNormalization(std::vector<double> expo, std::vector<double> contr,
                                                            unsigned int angularMomentum) {
  // The following reverses the normalization always performed by libint::shell() if it is constructed from scratch
  // and not by the copy constructor. [see libint::Shell::renorm()]
  double two_to_am = pow(2, angularMomentum);
  double sqrt_pi_cubed = sqrt(pow(PI, 3));
  unsigned int doubleFactorial = double_factorial(2 * angularMomentum - 1);
  unsigned int nIJ = contr.size();

  std::vector<double> newContractions;

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
