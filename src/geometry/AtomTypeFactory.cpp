/**
 * @file   AtomTypeFactory.cpp
 *
 * @date   Mar 19, 2013
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
#include "geometry/AtomTypeFactory.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
#include "parameters/AtomicParameters.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <ctype.h>
#include <utility>
#include <vector>

namespace Serenity {

std::map<std::string, std::shared_ptr<const AtomType>> AtomTypeFactory::_atomTypes;

std::shared_ptr<const AtomType> AtomTypeFactory::getAtomType(const std::string name) {
  std::string str = name;
  for (auto& c : str)
    c = toupper(c);
  auto it = _atomTypes.find(str);
  if (it != _atomTypes.end()) {
    return _atomTypes.at(str);
  }
  else {
    AtomTypeFactory::generateAtomType(str);
    it = _atomTypes.find(str);
    if (it != _atomTypes.end()) {
      return _atomTypes.at(str);
    }
    else {
      throw SerenityError((std::string) "There is no atom type defined with name " + name);
    }
  }
}

/*
 * This macro helps to define element and isotope atom types. Some data is taken from tabulated values
 * for the elements.
 * Usage: see below. For "H" the macro has not yet been used (technical reason), but the code for the
 * elements below that will look similar after expansion of the macro.
 */
// clang-format off
#define NEW_ATOM_TYPE(NAME, ATOMIC_NUMBER, MASS)                                        \
        } else if (name == NAME ) {                                                     \
          _atomTypes.insert(std::pair<std::string, std::shared_ptr<const AtomType> >(   \
          NAME, std::shared_ptr<const AtomType>(new AtomType(                           \
          NAME,                                                                         \
          ATOMIC_NUMBER,                                                                \
          MASS,                                                                         \
          ELEMENTAL_BRAGG_SLATER_RADII[ATOMIC_NUMBER] * ANGSTROM_TO_BOHR,               \
          ELEMENTAL_VAN_DER_WAALS_RADII[ATOMIC_NUMBER] * ANGSTROM_TO_BOHR,              \
          ELEMENTAL_UFF_RADII[ATOMIC_NUMBER] * ANGSTROM_TO_BOHR,                        \
          NUMBER_OF_CORE_ELECTRONS[ATOMIC_NUMBER],                                      \
          ELEMENT_OCCUPATIONS[ATOMIC_NUMBER],                                           \
          CHEMICAL_HARDNESS[ATOMIC_NUMBER]))));
// clang-format on

void AtomTypeFactory::generateAtomType(std::string name) {
  /*
   * To define your own atom type look at the definition of the hydrogen atom type and make use
   * of tabulated data for the elements if appropriate (as is done in the macro above).
   */
  /*
   * Elements:
   *
   * (data of the respectively most abundant isotope are used, other isotopes: see below)
   * The masses of elements and isotopes have been taken from
   * http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl (column: "Relative Atomic Mass").
   * The data (although not listed in that paper) belong to:
   * Peter J. Mohr, David B. Newell, Barry N. Taylor:
   * "CODATA Recommended Values of the Fundamental Physical Constants: 2014",
   * arXiv:1507.07956.
   * As a comment behind the element/isotope we give the natural abundance from the same source.
   * TODO move the masses and natural abundances to the AtomicParameters.h
   */
  for (auto& c : name)
    c = toupper(c);

  if (name.substr(name.length() - 1) == ":") {
    auto atomType = AtomTypeFactory::getAtomType(name.substr(0, name.length() - 1));
    _atomTypes.insert(std::pair<std::string, std::shared_ptr<const AtomType>>(
        name, std::make_shared<AtomType>(name, atomType->getNuclearCharge(), 0.0, atomType->getBraggSlaterRadius(),
                                         atomType->getVanDerWaalsRadius(), atomType->getUFFRadius(), atomType->getNCoreElectrons(),
                                         atomType->getOccupations(), atomType->getChemicalHardness(), true)));
  }
  else if (name == "H") {
    _atomTypes.insert(std::pair<std::string, std::shared_ptr<const AtomType>>(
        "H", std::shared_ptr<const AtomType>(new AtomType(
                 "H", 1, 1.00782503223, 0.25 * ANGSTROM_TO_BOHR, 1.20 * ANGSTROM_TO_BOHR, 1.4430 * ANGSTROM_TO_BOHR, 0,
                 {{{ANGULAR_QUANTUM_NUMBER::s, 1}}}, CHEMICAL_HARDNESS[1], false))));
    // 0.999 885
    NEW_ATOM_TYPE("HE", 2, 4.00260325413) // 0.999 998 66
    // 2nd row
    NEW_ATOM_TYPE("LI", 3, 7.0160034366)   // 0.9241
    NEW_ATOM_TYPE("BE", 4, 9.012183065)    // 1
    NEW_ATOM_TYPE("B", 5, 11.00930536)     // 0.801
    NEW_ATOM_TYPE("C", 6, 12.0)            // 0.9893, Mass is exactly 12 by definition
    NEW_ATOM_TYPE("N", 7, 14.00307400443)  // 0.996 36
    NEW_ATOM_TYPE("O", 8, 15.99491461957)  // 0.997 57
    NEW_ATOM_TYPE("F", 9, 18.99840316273)  // 1
    NEW_ATOM_TYPE("NE", 10, 19.9924401762) // 0.9048
    // 3rd row
    NEW_ATOM_TYPE("NA", 11, 22.9897692820)  // 1
    NEW_ATOM_TYPE("MG", 12, 23.985041697)   // 0.7899
    NEW_ATOM_TYPE("AL", 13, 26.98153853)    // 1
    NEW_ATOM_TYPE("SI", 14, 27.97692653465) // 0.922 23
    NEW_ATOM_TYPE("P", 15, 30.97376199842)  // 1
    NEW_ATOM_TYPE("S", 16, 31.9720711744)   // 0.9499
    NEW_ATOM_TYPE("CL", 17, 34.968852682)   // 0.7576
    NEW_ATOM_TYPE("AR", 18, 39.9623831237)  // 0.996 035
    // 4th row (incl. transition metals)
    NEW_ATOM_TYPE("K", 19, 38.9637064864)  // 0.932 581
    NEW_ATOM_TYPE("CA", 20, 39.962590863)  // 0.969 41
    NEW_ATOM_TYPE("SC", 21, 44.95590828)   // 1
    NEW_ATOM_TYPE("TI", 22, 47.94794198)   // 0.7372
    NEW_ATOM_TYPE("V", 23, 50.94395704)    // 0.997 50
    NEW_ATOM_TYPE("CR", 24, 51.94050623)   // 0.837 89
    NEW_ATOM_TYPE("MN", 25, 54.93804391)   // 1
    NEW_ATOM_TYPE("FE", 26, 55.93493633)   // 0.917 54
    NEW_ATOM_TYPE("CO", 27, 58.93319429)   // 1
    NEW_ATOM_TYPE("NI", 28, 57.93534241)   // 0.680 77
    NEW_ATOM_TYPE("CU", 29, 62.92959772)   // 0.6915
    NEW_ATOM_TYPE("ZN", 30, 63.92914201)   // 0.4917
    NEW_ATOM_TYPE("GA", 31, 68.9255735)    // 0.601 08
    NEW_ATOM_TYPE("GE", 32, 73.921177761)  // 0.3650
    NEW_ATOM_TYPE("AS", 33, 74.92159457)   // 1
    NEW_ATOM_TYPE("SE", 34, 79.9165218)    // 0.4961
    NEW_ATOM_TYPE("BR", 35, 78.9183376)    // 0.5069
    NEW_ATOM_TYPE("KR", 36, 83.9114977282) // 0.569 87
    // 5th row (incl. transition metals)
    NEW_ATOM_TYPE("RB", 37, 84.9117897379) // 0.7217
    NEW_ATOM_TYPE("SR", 38, 87.9056125)    // 0.8258
    NEW_ATOM_TYPE("Y", 39, 88.9058403)     // 1
    NEW_ATOM_TYPE("ZR", 40, 89.9046977)    // 0.5145
    NEW_ATOM_TYPE("NB", 41, 92.9063730)    // 1
    NEW_ATOM_TYPE("MO", 42, 97.90540482)   // 0.2439
    /*
     * Technetium is the first element of which no isotopes occur in nature (except for traces).
     * In these cases, we simply use the most stable isotope here to have an entry at all.
     * No natural abundance is given, for obvious reasons.
     */
    NEW_ATOM_TYPE("TC", 43, 97.9072124)
    NEW_ATOM_TYPE("RU", 44, 101.9043441)   // 0.3155
    NEW_ATOM_TYPE("RH", 45, 102.9054980)   // 1
    NEW_ATOM_TYPE("PD", 46, 105.9034804)   // 0.2733
    NEW_ATOM_TYPE("AG", 47, 106.9050916)   // 0.518 39
    NEW_ATOM_TYPE("CD", 48, 113.90336509)  // 0.2873
    NEW_ATOM_TYPE("IN", 49, 114.903878776) // 0.9571
    NEW_ATOM_TYPE("SN", 50, 119.90220163)  // 0.3258
    NEW_ATOM_TYPE("SB", 51, 120.9038120)   // 0.5721
    NEW_ATOM_TYPE("TE", 52, 129.906222748) // 0.3408
    NEW_ATOM_TYPE("I", 53, 126.9044719)    // 1
    // 6th row (incl. lantahnoids)
    NEW_ATOM_TYPE("XE", 54, 131.9041550856) // 0.269086
    NEW_ATOM_TYPE("CS", 55, 132.9054519610) // 1
    NEW_ATOM_TYPE("BA", 56, 137.90524700)   // 0.71698
    NEW_ATOM_TYPE("LA", 57, 138.9063563)    // 0.9991119
    NEW_ATOM_TYPE("CE", 58, 139.9054431)    // 0.88450
    NEW_ATOM_TYPE("PR", 59, 140.9076576)    // 1
    NEW_ATOM_TYPE("ND", 60, 141.9077290)    // 0.27152
    NEW_ATOM_TYPE("PM", 61, 144.9127559)
    NEW_ATOM_TYPE("SM", 62, 151.9197397)  // 0.2675
    NEW_ATOM_TYPE("EU", 63, 152.9212380)  // 0.5219
    NEW_ATOM_TYPE("GD", 64, 157.9241123)  // 0.2484
    NEW_ATOM_TYPE("TB", 65, 158.9253547)  // 1
    NEW_ATOM_TYPE("DY", 66, 163.9291819)  // 0.28260
    NEW_ATOM_TYPE("HO", 67, 164.9303288)  // 1
    NEW_ATOM_TYPE("ER", 68, 165.9302995)  // 0.33503
    NEW_ATOM_TYPE("TM", 69, 168.9342179)  // 1
    NEW_ATOM_TYPE("YB", 70, 173.9388664)  // 0.32026
    NEW_ATOM_TYPE("LU", 71, 174.9407752)  // 0.97401
    NEW_ATOM_TYPE("HF", 72, 179.9465570)  // 0.3508
    NEW_ATOM_TYPE("TA", 73, 180.9479958)  // 0.9998799
    NEW_ATOM_TYPE("W", 74, 183.95093092)  // 0.3064
    NEW_ATOM_TYPE("RE", 75, 186.9557501)  // 0.6260
    NEW_ATOM_TYPE("OS", 76, 191.9614770)  // 0.4078
    NEW_ATOM_TYPE("IR", 77, 192.9629216)  // 0.627
    NEW_ATOM_TYPE("PT", 78, 194.9647917)  // 0.3378
    NEW_ATOM_TYPE("AU", 79, 196.96656879) // 1
    NEW_ATOM_TYPE("HG", 80, 201.97064340) // 0.2986
    NEW_ATOM_TYPE("TL", 81, 204.9744278)  // 0.7048
    NEW_ATOM_TYPE("PB", 82, 207.9766525)  // 0.524
    NEW_ATOM_TYPE("BI", 83, 208.9803991)  // 1
    NEW_ATOM_TYPE("PO", 84, 208.9824308)
    NEW_ATOM_TYPE("AT", 85, 209.9871479)
    NEW_ATOM_TYPE("RN", 86, 222.0175782)
    // 7th row (incl. actinoids)
    NEW_ATOM_TYPE("FR", 87, 223.0197360)
    NEW_ATOM_TYPE("RA", 88, 226.0254103)
    NEW_ATOM_TYPE("AC", 89, 227.0277523)
    NEW_ATOM_TYPE("TH", 90, 232.0380558) // 1
    NEW_ATOM_TYPE("PA", 91, 231.0358842) // 1
    NEW_ATOM_TYPE("U", 92, 238.0507884)  // 0.992742
    NEW_ATOM_TYPE("NP", 93, 237.0481736)
    NEW_ATOM_TYPE("PU", 94, 244.0642053)
    NEW_ATOM_TYPE("AM", 95, 243.0613813)
    NEW_ATOM_TYPE("CM", 96, 247.0703541)
    NEW_ATOM_TYPE("BK", 97, 247.0703073)
    NEW_ATOM_TYPE("CF", 98, 251.0795886)
    NEW_ATOM_TYPE("ES", 99, 252.082980)
    NEW_ATOM_TYPE("FM", 100, 257.0951061)
    NEW_ATOM_TYPE("MD", 101, 258.0984315)
    NEW_ATOM_TYPE("NO", 102, 259.10103)
    NEW_ATOM_TYPE("LR", 103, 262.10961)
    NEW_ATOM_TYPE("RF", 104, 267.12179)
    NEW_ATOM_TYPE("DB", 105, 268.12567)
    NEW_ATOM_TYPE("SG", 106, 271.13393)
    NEW_ATOM_TYPE("BH", 107, 272.13826)
    NEW_ATOM_TYPE("HS", 108, 270.13429)
    NEW_ATOM_TYPE("MT", 109, 276.15159)
    NEW_ATOM_TYPE("DS", 110, 281.16451)
    NEW_ATOM_TYPE("RG", 111, 280.16514)
    NEW_ATOM_TYPE("CN", 112, 285.17712)
    NEW_ATOM_TYPE("NH", 113, 284.17873)
    NEW_ATOM_TYPE("FL", 114, 289.19042)
    NEW_ATOM_TYPE("MC", 115, 288.19274)
    NEW_ATOM_TYPE("LV", 116, 293.20449)
    NEW_ATOM_TYPE("TS", 117, 292.20746)
    NEW_ATOM_TYPE("OG", 118, 294.21392)
    /*
     * Isotopes
     */
    NEW_ATOM_TYPE("D", 1, 2.01410177812) // 0.000 115
    NEW_ATOM_TYPE("T", 1, 3.0160492779)  // small
  }
  else {
    throw SerenityError((std::string) "There is no atom type defined with name " + name);
  }
}

#undef NEW_ATOM_TYPE

AtomTypeFactory& AtomTypeFactory::getInstance() {
  // The only instance
  // Guaranteed to be lazy initialized
  // Guaranteed that it will be destroyed correctly
  static AtomTypeFactory instance;
  return instance;
}

} /* namespace Serenity */
