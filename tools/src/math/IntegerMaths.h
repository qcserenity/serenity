/**
 * @file   IntegerMaths.h
 *
 * @date   Aug 13, 2013
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
#ifndef INTEGERMATHS_H_
#define INTEGERMATHS_H_
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {

/**
 * @param   arg
 * @returns arg\!
 */
inline static constexpr unsigned int factorial(const unsigned int& arg) {
  return (arg == 0 || arg == 1) ? 1 : arg * factorial(arg - 1);
}

/**
 * @param   arg
 * @returns arg\!\!
 */
inline static constexpr unsigned int double_factorial(const int& arg) {
  return (arg < 2) ? 1 : arg * double_factorial(arg - 2);
}

/**
 * @param arg
 * @returns true if arg is an even integer number
 */
inline static constexpr bool isEven(const unsigned int& arg) {
  return (arg % 2 == 0);
}

/**
 * @param base The base
 * @param exp  The exponent.
 * @return The base to the power of exp.
 */
inline static constexpr unsigned int intPow(const unsigned int base, const unsigned int exp) {
  return (exp == 0) ? 1 : (exp == 1) ? base : base * intPow(base, exp - 1);
}

// Taken from: https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
// inline static unsigned int binomial(unsigned int n, unsigned int k) {
//  assert(n > k && "For the calculation of a binomial coefficient n > k has to be true!");
//  int res = 1;
//  if (k > n-k) k=n-k;
//  for (unsigned int i = 0; i < k; ++i) {
//    res *= (n-i);
//    res /= (i+1);
//  }
//  return res;
//}

} /* namespace Serenity */
#endif /* INTEGERMATHS_H_ */
