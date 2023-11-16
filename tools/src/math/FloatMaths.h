/**
 * @file   FloatMaths.h
 *
 * @date   Nov 12, 2013
 * @author Thomas Dresselhaus
 *
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
#ifndef FLOATMATHS_H_
#define FLOATMATHS_H_
/* Include Std and External Headers */
#include <cmath>
#include <iostream>
#include <limits>

namespace Serenity {

/**
 * These are accuracy thresholds to be used globally, e.g. for float comparisons.
 */
constexpr double SLOPPY_D = 1e-3;
constexpr double LOOSE_D = 1e-6;
constexpr double NORMAL_D = 1e-9;
constexpr double TIGHT_D = 1e-12;

/**
 * @brief   Check if two numbers are equal.
 * @param   lhs, rhs the numbers to be compared.
 * @param   epsilon the precision to which the numbers shall be equal, default: 100 \* machine precision.
 * @returns true iff lhs and rhs are equal to up to the precision defined by epsilon
 */
inline bool isEqual(const double& lhs, const double& rhs, const double epsilon = 100 * std::numeric_limits<double>::epsilon()) {
  return (fabs((lhs - rhs) / (rhs + lhs)) < epsilon || fabs(lhs - rhs) < epsilon);
}

/**
 * @brief own pow with unsigned int as exponent, i.e. base**exp
 *
 * @param base basis for an exponential expression
 * @param exp exponent
 * @return The double to the power of exp.
 */
constexpr double fipow(const double& base, const unsigned int& exp) {
  return (exp >= 1) ? base * fipow(base, exp - 1) : 1;
}

} /* namespace Serenity */
#endif /* FLOATMATHS_H_ */
