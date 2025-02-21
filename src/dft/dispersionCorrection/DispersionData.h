/**
 * @file   DispersionData.h
 *
 * @date   Nov 26, 2015
 * @author Jan Unsleber
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

/* Include Serenity Internal Headers */
#include "dft/dispersionCorrection/DispersionRawData.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <array>
#include <cassert>
#include <cmath>

#ifndef DISPERSIONDATA_H_
#define DISPERSIONDATA_H_

namespace Serenity {
/* Forward Declarations */
namespace CompositeFunctionals {
enum class XCFUNCTIONALS;
}
namespace Options {
enum class DFT_DISPERSION_CORRECTIONS;
}

/**
 * @class DispersionData DispersionData.h
 * @brief Class providing usable versions of all the parameters for D3(0)/D3(BJ)
 *
 */
class DispersionData {
 public:
  DispersionData() = default;
  virtual ~DispersionData() = default;

  /**
   * @brief Returns the parameters for a given functional
   *
   * @param functional The functional the parameters shall be given for.
   * @param s6         Variable for the s6 parameter.
   * @param rs6        Variable for the rs6 parameter.
   * @param s18        Variable for the s18 parameter.
   * @param rs18       Variable for the rs18 parameter.
   * @param alp        Variable for the alp parameter.
   */
  template<Options::DFT_DISPERSION_CORRECTIONS dispLevel>
  static void getFunctionalParameters(CompositeFunctionals::XCFUNCTIONALS functional, double& s6, double& rs6,
                                      double& s18, double& rs18, double& alp);

  /**
   * @brief Wraps the raw r0ab data.
   *
   * @param Nuclear charge of the first atom.
   * @param Nuclear charge of the second atom.
   * @return r0ab element (i,j)
   */
  static double r0ab(unsigned int i, unsigned int j) {
    // in the case of a ghost atom the nuclear charge is 0 and no dispersion contribution
    // should be calculated.
    assert(i * j > 0 && "no r0ab data for atoms with charge 0 or lower!");
    return (j <= i) ? DispersionRawData::r0abRaw[(j - 1) + (i - 1) * i / 2] * ANGSTROM_TO_BOHR
                    : DispersionRawData::r0abRaw[(i - 1) + (j - 1) * j / 2] * ANGSTROM_TO_BOHR;
  }

  /**
   * @brief Wraps the raw rcov data.
   *
   * @param Nuclear charge of the atom.
   * @return Covalent radius of the element with the given nuclear charge.
   */
  static double rcov(unsigned int i) {
    assert(i > 0 && "A ghost atom or atom of charge 0 does not have a covalent radius");
    return DispersionRawData::K2 * DispersionRawData::rcovRaw[i - 1] * ANGSTROM_TO_BOHR;
  }

  /**
   * @brief Wraps the raw r2r4 data.
   *
   * @param i
   * @return r2r4 element (i)
   */
  static double r2r4(unsigned int i) {
    /*   Original comment:
     scale r4/r2 values of the atoms by sqrt(Z)
     sqrt is also globally close to optimum
     together with the factor 1/2 this yield reasonable
     c8 for he, ne and ar. for larger Z, C8 becomes too large
     which effectively mimics higher R^n terms neglected due
     to stability reasons
    */
    assert(i > 0 && "r2r4 not defined for nuclear charge < 0. Did you use ghost atoms?");
    return sqrt(0.5 * DispersionRawData::r2r4Raw[i - 1] * sqrt(double(i)));
  }

  /**
   * @brief Wraps the raw rcov data.
   *
   * @param Nuclear charge of the atom.
   * @return Maximum number of C6 parameters of the element with the given nuclear charge.
   */
  static unsigned int maxNC6(unsigned int i) {
    return DispersionRawData::maxNC6Raw[i - 1];
  }

  /**
   * @brief Wraps the raw c6ab data.
   *
   * @param i Atom number/nuclear charge of the first atom.
   * @param j Atom number/nuclear charge of the second atom.
   * @param k Range: [0,maxNC6(i)]
   * @param l Range: [0,maxNC6(j)]
   * @param m 0=y, 1=cn1, 2=cn2
   * @return c6ab element (i,j,k,l,m)
   */
  static double c6ab(unsigned int i, unsigned int j, unsigned int k, unsigned int l, unsigned int m) {
    return DispersionRawData::c6abRaw[i - 1][j - 1][k][l][m];
  }

  /*
   * Forwarded parameters
   */
  static constexpr int MAX_ELEMENT_NUMBER = DispersionRawData::MAX_ELEMENT_NUMBER;
  static constexpr double K1 = DispersionRawData::K1;
  static constexpr double K2 = DispersionRawData::K2;
  static constexpr double K3 = DispersionRawData::K3;
  static constexpr double RANGE_THRESHOLD_N2 = DispersionRawData::RANGE_THRESHOLD_N2;
  static constexpr double RANGE_THRESHOLD_N3 = DispersionRawData::RANGE_THRESHOLD_N3;
};

} /* namespace Serenity */

#endif /* DISPERSIONDATA_H_ */
