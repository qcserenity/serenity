/**
 * @file   VectorMaths.h
 *
 * @date   Sep 28, 2014
 * @author Thomas Dresselhaus
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
#ifndef VECTORMATHS_H_
#define VECTORMATHS_H_
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <cassert>
#include <vector>

namespace Serenity {
/**
 * @param[in] lhs will be summed with rhs
 * @param[in] rhs will be summed with lhs
 * @returns a vector with elements which are the sum of the corresponding elements of rhs and lhs
 */
template<typename T>
inline std::vector<T> operator+(const std::vector<T>& lhs, const std::vector<T>& rhs) {
  if (lhs.size() != rhs.size())
    throw SerenityError("VectorMath: Vectors are of different sizes.");
  std::vector<T> result;
  result.resize(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); ++i) {
    result[i] = lhs[i] + rhs[i];
  }
  return result;
}
/**
 * @param[inout] lhs rhs will be added to this
 * @param[in] rhs is added to lhs
 * @returns the mutated lhs (as a reference; the actually incoming lhs is mutated)
 */
template<typename T>
inline std::vector<T>& addToVector(std::vector<T>& lhs, const std::vector<T>& rhs) {
  if (lhs.size() != rhs.size())
    throw SerenityError("VectorMath: Vectors are of different sizes.");
  for (unsigned int i = 0; i < lhs.size(); ++i) {
    lhs[i] += rhs[i];
  }
  return lhs;
}
/**
 * @param[in] lhs from this rhs will be subtracted
 * @param[in] rhs will be subtracted from lhs
 * @returns a vector with elements which are the difference between the corresponding elements of
 *          rhs and lhs
 */
template<typename T>
inline std::vector<T> operator-(const std::vector<T>& lhs, const std::vector<T>& rhs) {
  if (lhs.size() != rhs.size())
    throw SerenityError("VectorMath: Vectors are of different sizes.");
  std::vector<T> result;
  result.resize(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); ++i) {
    result[i] = lhs[i] - rhs[i];
  }
  return result;
}
/**
 * @param[inout] lhs rhs will be subtracted from this
 * @param[in] rhs is subtracted from lhs
 * @returns the mutated lhs (as a reference; the actually incoming lhs is mutated)
 */
template<typename T>
inline std::vector<T>& subtractFromVector(std::vector<T>& lhs, const std::vector<T>& rhs) {
  if (lhs.size() != rhs.size())
    throw SerenityError("VectorMath: Vectors are of different sizes.");
  for (unsigned int i = 0; i < lhs.size(); ++i) {
    lhs[i] -= rhs[i];
  }
  return lhs;
}
/**
 * @param[in] lhs will be summed with rhs
 * @param[in] rhs will be summed with lhs
 * @returns the scalar product of the two vectors which is the sum of the corresponding elements of rhs and lhs
 */
template<typename T>
inline T scalarProduct(const std::vector<T>& lhs, const std::vector<T>& rhs) {
  if (lhs.size() != rhs.size())
    throw SerenityError("VectorMath: Vectors are of different sizes.");
  T result = 0;
  for (unsigned int i = 0; i < lhs.size(); ++i) {
    result += lhs[i] * rhs[i];
  }
  return result;
}

} /* namespace Serenity */

#endif /* VECTORMATHS_H_ */
