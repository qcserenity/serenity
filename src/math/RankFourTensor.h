/**
 * @file   RankFourTensor.h
 *
 * @date   28. Juli 2013, 20:47
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
#ifndef RANKFOURTENSOR_H
#define RANKFOURTENSOR_H
/* Include Std and External Headers */
#include <cassert>
#include <iostream>

namespace Serenity {
/**
 * @class RankFourTensor RankFourTensor.h
 *
 * @brief an A x B x C x D type of object. Can be accessed via foo(i,j,k,l).
 */
template<typename T>
class RankFourTensor {
 public:
  /**
   * Creates a filled tensor.
   *
   * @param nDimOne, nDimTwo, nDimThree, nDimFour the length in each dimension
   * @param data the contents of the produced tensor are taken from it
   */
  RankFourTensor(const int& nDimOne, const int& nDimTwo, const int& nDimThree, const int& nDimFour, T* data)
    : _nDimOne(nDimOne),
      _nDimTwo(nDimTwo),
      _nDimThree(nDimThree),
      _nDimFour(nDimFour),
      _kOffset(_nDimFour),
      _jOffset(_kOffset * _nDimThree),
      _iOffset(_jOffset * nDimTwo),
      _nElements(_iOffset * nDimOne) {
    _data = new T[_nElements];
    for (unsigned int i = 0; i < _nElements; i++) {
      _data[i] = data[i];
    }
  }
  /**
   * Creates an empty tensor, initialized to zero.
   *
   * @param nDimOne, nDimTwo, nDimThree, nDimFour the length in each dimension
   */
  RankFourTensor(const int& nDimOne, const int& nDimTwo, const int& nDimThree, const int& nDimFour)
    : // Tensor<T>(nDimOne*nDimTwo*nDimThree*nDimFour),
      _nDimOne(nDimOne),
      _nDimTwo(nDimTwo),
      _nDimThree(nDimThree),
      _nDimFour(nDimFour),
      _kOffset(_nDimFour),
      _jOffset(_kOffset * _nDimThree),
      _iOffset(_jOffset * nDimTwo),
      _nElements(_iOffset * nDimOne) {
    _data = new T[_nElements];
    for (unsigned int i = 0; i < _nElements; i++) {
      _data[i] = 0.0;
    }
  }

  virtual ~RankFourTensor() {
    delete[] _data;
  }

  /**
   * @brief The length in each dimension
   */
  const unsigned int _nDimOne, _nDimTwo, _nDimThree, _nDimFour;

  /**
   * @brief Straightforward and rather safe access to the tensor. Slow!
   * @param i, j, k, l the index in each of the four dimensions.
   * @return Returns the value stored at this position in the tensor.
   *
   * For a faster access in time critical parts use the data() function. But be careful with
   * that, using that can more easily cause errors.
   */
  inline const T& operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const {
    assert(i < this->_nDimOne && j < this->_nDimTwo && k < this->_nDimThree && l < this->_nDimFour);
    return this->_data[i * _iOffset + j * _jOffset + k * _kOffset + l];
  }
  /**
   * @brief Straightforward and rather safe access to the tensor. Slow!
   * @param i, j, k, l the index in each of the four dimensions.
   * @return Returns the value stored at this position in the tensor.
   *
   * For a faster access in time critical parts use the data() function. But be careful with
   * that, using that can more easily cause errors.
   */
  inline T& operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
    assert(i < this->_nDimOne && j < this->_nDimTwo && k < this->_nDimThree && l < this->_nDimFour);
    return this->_data[i * _iOffset + j * _jOffset + k * _kOffset + l];
  }
  /**
   * @brief Prints the tensor to stdout. Careful! Will produce a lot of output!
   */
  void print() const {
    for (unsigned int i = 0; i < this->_nDimOne; ++i) {
      for (unsigned int j = 0; j < this->_nDimTwo; ++j) {
        for (unsigned int k = 0; k < this->_nDimThree; ++k) {
          for (unsigned int l = 0; l < this->_nDimFour; ++l) {
            std::cout << std::scientific << i << " ";
            std::cout << std::scientific << j << " ";
            std::cout << std::scientific << k << " ";
            std::cout << std::scientific << l << " ";
            std::cout << std::scientific << (*this)(i, j, k, l) << " ";
            std::cout << std::endl;
          }
        }
      }
    }
  }

 private:
  int _kOffset, _jOffset, _iOffset;

  unsigned int _nElements;

  T* _data;
};

} /* namespace Serenity */
#endif /* RANKFOURTENSOR_H */
