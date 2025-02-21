/**
 * @file   BasisFunctionProvider.cpp
 *
 * @date   24.03.2013
 * @author Thomas Dresselhaus
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
#include "basis/BasisFunctionProvider.h"
/* Include Serenity Internal Headers */
#include "basis/Shell.h"
#include "geometry/Atom.h"
/* Include Std and External Headers */
#include <algorithm>
#include <fstream>
#include <libecpint/ecp.hpp>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Serenity {

void BasisFunctionProvider::provideAtomWithBasisFunction(Atom& atom, const std::string libraryPath,
                                                         const std::string basisType, const bool isSpherical,
                                                         const bool isPrimary, const int firstECP) {
  std::string basisPath = libraryPath + basisType;

  std::ifstream basisFile(&basisPath[0]);

  if (!basisFile.good()) {
    throw SerenityError((std::string) "Error while parsing basis file " + basisPath +
                        "\n Make sure the directory exists, change the path in the input, or set $SERENITY_RESOURCES.");
  }
  /*
   * The actual content of the loaded basis
   */
  std::string loadedBasisFile;
  basisFile.seekg(0, std::ios::end);
  loadedBasisFile.reserve(basisFile.tellg());
  basisFile.seekg(0, std::ios::beg);
  loadedBasisFile.assign(std::istreambuf_iterator<char>(basisFile), std::istreambuf_iterator<char>());
  basisFile.close();
  /*
   * Basis file loaded into string. Getting Element from atom
   * and determining the search string.
   */
  std::string element = atom.getAtomType()->getElementSymbol();
  if (element.length() == 1)
    element = "\n" + element;
  transform(element.begin(), element.end(), element.begin(), ::tolower);
  // Escape special characters + and * in basisType
  const std::regex specialChars("[+*]");
  std::string escapedBasisType = std::regex_replace(basisType, specialChars, R"(\$&)");
  const std::string searchString = element + R"(\s+)" + escapedBasisType;
  std::regex regexPattern(searchString, std::regex_constants::icase);

  const std::string errorMessage = (std::string) "Error while parsing basis set file " + basisPath + " for element " +
                                   element + ". (Basis type: " + basisType + ")\n";

  /* Searching for entry. */
  std::smatch match;
  unsigned matchsize = 0;
  if (std::regex_search(loadedBasisFile, match, regexPattern)) {
    std::string matches = match[0];
    matchsize = matches.size();
  }
  else {
    throw SerenityError(errorMessage + "The used basis (file) is not defined for this element.");
  }
  std::string::size_type startPosition = match.position() + matchsize;
  /* Turbomole-format basis files sometimes have a comment along the lines of # o     (7s4p1d) / [3s2p1d]     {511/31/1}
  If the next line begins with a hash, skip it */
  if (loadedBasisFile[startPosition + 1] == '#') {
    startPosition += loadedBasisFile.find("\n", startPosition + 1);
  }
  /* The basis information for one element is terminated by a star */
  std::string::size_type endPosition = loadedBasisFile.find("*", startPosition + 3);

  /* Extract the basis data for this element */
  std::string substring = loadedBasisFile.substr(startPosition + 3, endPosition - startPosition - 3);
  /* Replace D+ and D- by E+ and E-, respectively */
  std::string modifiedstring = std::regex_replace(substring, std::regex(R"(D(\+|-))"), "E$1");
  std::stringstream workStream(modifiedstring);

  std::vector<std::shared_ptr<Shell>> basisFunctions;
  int nPrimitives;
  char type;
  /* loop over basis functions for this atom */
  while (workStream >> nPrimitives) {
    if (!(workStream >> type))
      throw SerenityError(errorMessage + "Number of primitives for a contraction could not be parsed.");
    const unsigned int angularMomentum = resolveAngularMomentumChar(type, errorMessage);

    libint2::svector<double> exponents;
    libint2::svector<double> contractions;
    /* loop over primitives for this basis function */
    for (int i = 0; i < nPrimitives; i++) {
      double exponent;
      double contraction;
      if (!(workStream >> exponent))
        throw SerenityError(errorMessage + "A possible source of error is a mismatch between the number of primitives "
                                           "and the number of exponent-contraction entries.");
      if (!(workStream >> contraction))
        throw SerenityError(errorMessage + "A possible source of error is a mismatch between the number of primitives "
                                           "and the number of exponent-contraction entries.");

      exponents.emplace_back(exponent);
      contractions.emplace_back(contraction);
    }
    std::array<double, 3> coords;
    coords[0] = atom.getX();
    coords[1] = atom.getY();
    coords[2] = atom.getZ();
    /* Create the basis function, centered on atom */
    basisFunctions.push_back(std::make_shared<Shell>(exponents, contractions, angularMomentum, isSpherical, coords, element));
  }
  /*
   * Search for effective core potentials
   */
  const std::string::size_type ecpStart = loadedBasisFile.find("$ecp");
  if (ecpStart == loadedBasisFile.npos || !isPrimary) {
    // No ECPs found in the file -> return now
    atom.addBasis(make_pair(basisType, basisFunctions), isPrimary);
    return;
  }
  // Search for ECP for this element
  const std::string searchStringECP = element + "  "; // The basis name may be different for the ECP
  const std::string::size_type startPositionECP = loadedBasisFile.find(searchStringECP, ecpStart);
  if (startPositionECP == loadedBasisFile.npos) {
    // No ECPs found in the file for this element -> return now
    atom.addBasis(make_pair(basisType, basisFunctions), isPrimary);
    return;
  }
  // Create empty ECP object to be filled below
  auto ecp = std::make_shared<libecpint::ECP>();
  unsigned int nCoreElectrons = 0;
  if (atom.getNuclearCharge() >= firstECP and isPrimary and !atom.isDummy()) {
    // The basis set contains an ECP for this element. It must be used as the primary basis.
    const std::string::size_type nCorePosition = loadedBasisFile.find("ncore", startPositionECP);
    const std::string::size_type endPositionECP = loadedBasisFile.find("*", nCorePosition);
    /* Extract the basis data for this element */
    std::stringstream workStreamECP(loadedBasisFile.substr(nCorePosition + 8, endPositionECP - (nCorePosition + 8)));
    if (!(workStreamECP >> nCoreElectrons))
      throw SerenityError(errorMessage + "Reading in ECP failed, ncore could not be parsed.");
    std::string unneededContent;
    workStreamECP >> unneededContent;
    workStreamECP >> unneededContent;
    unsigned int lmax;
    workStreamECP >> lmax;
    ecp->setPos(atom.getX(), atom.getY(), atom.getZ());
    unsigned int angularMomentumECP = 999999;
    std::string testString;
    // The angular momenta are given in a specific order: l_max, 0, 1,... l_max-1
    bool firstAngularMomentum = true;
    unsigned int angularMomentumCounter = 0;
    while (workStreamECP >> testString) {
      // Test whether a primitive follows
      std::stringstream maybeDouble(testString);
      double c;
      if (maybeDouble >> c) {
        if (angularMomentumECP == 999999)
          throw SerenityError(errorMessage +
                              "Reading in ECP failed, no angular momentum has been read while parsing a primitive.");
        // We are within the definition of a primitive
        unsigned int n;
        if (!(workStreamECP >> n))
          throw SerenityError(errorMessage + "Reading in ECP failed, n in primitive could not be parsed.");
        double exponentECP;
        if (!(workStreamECP >> exponentECP))
          throw SerenityError(errorMessage + "Reading in ECP failed, coefficient in primitive could not be parsed.");
        // The exponent is always greater or equal to 0
        assert(exponentECP >= 0.0);
        ecp->addPrimitive(n, angularMomentumECP, exponentECP, c, false);
      }
      else {
        // We arrived at a new angular momentum
        angularMomentumECP = resolveAngularMomentumChar(testString[0], errorMessage);
        // Check input consistency
        if ((firstAngularMomentum && lmax != angularMomentumECP) ||
            (!firstAngularMomentum && angularMomentumECP != angularMomentumCounter))
          throw SerenityError(errorMessage + "Reading in ECP failed, angular momentum is in a wrong order.");
        if (!firstAngularMomentum)
          angularMomentumCounter++;
        firstAngularMomentum = false;
      }
    }
  }
  ecp->sort();
  atom.addBasis(make_pair(basisType, basisFunctions), ecp, nCoreElectrons);
}

unsigned int BasisFunctionProvider::resolveAngularMomentumChar(char type, std::string errorMessage) {
  unsigned int angularMomentum = 999999;
  if (type == 's') {
    angularMomentum = 0;
  }
  else if (type == 'p') {
    angularMomentum = 1;
  }
  else if (type == 'd') {
    angularMomentum = 2;
  }
  else if (type == 'f') {
    angularMomentum = 3;
  }
  else if (type == 'g') {
    angularMomentum = 4;
  }
  else if (type == 'h') {
    angularMomentum = 5;
  }
  else if (type == 'i') {
    angularMomentum = 6;
  }
  else if (type == 'k') {
    angularMomentum = 7;
  }
  else if (type == 'm') {
    angularMomentum = 8;
  }
  else if (type == 'n') {
    angularMomentum = 9;
  }
  else if (type == 'o') {
    angularMomentum = 10;
  }
  if (angularMomentum == 999999) {
    throw SerenityError(errorMessage +
                        "The BasisFunctionProvider tries to read in a basis with an unknown symbol "
                        "for the angular momentum. Either you try to use a basis function with a ridiculously "
                        "high angular momentum, or the basis file is corrupted.");
  }
  return angularMomentum;
}

} /* namespace Serenity */
