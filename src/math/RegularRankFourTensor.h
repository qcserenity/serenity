/**
 * @file   RegularRankFourTensor.h
 *
 * @date   Sep 1, 2013
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
#ifndef REGULARRANKFOURTENSOR_H_
#define REGULARRANKFOURTENSOR_H_

/* Include Std and External Headers */
#include <cassert>
#include <iostream>

namespace Serenity {
/**
 * @class RegularRankFourTensor RegularRankFourTensor.h
 * @brief An object of type A x A x A x A. Can be accessed via foo(i,j,k,l)
 *
 * A four-dimensional object with the same length in each dimension. Used to store electron
 * repulsion integrals.
 */
template<typename T>
class RegularRankFourTensor {
 public:
  /**
   * @brief Copy constructor.
   * @param old tensor
   */
  RegularRankFourTensor(const RegularRankFourTensor<T>& tensor)
    : _lengthPerDim(tensor._lengthPerDim),
      _kOffset(tensor._lengthPerDim),
      _jOffset(_kOffset * tensor._lengthPerDim),
      _iOffset(_jOffset * tensor._lengthPerDim),
      _nElements(_iOffset * tensor._lengthPerDim) {
    _data = new T[tensor._nElements];
    for (unsigned int i = 0; i < tensor._nElements; i++) {
      _data[i] = tensor._data[i];
    }
  }
  /**
   * @brief Produces a filled tensor.
   * @param lengthPerDim number of entries into each dimension.
   * @param data the content of the produced tensor is taken from it
   */
  RegularRankFourTensor(const unsigned int& lengthPerDim, T* data)
    : _lengthPerDim(lengthPerDim),
      _kOffset(lengthPerDim),
      _jOffset(_kOffset * lengthPerDim),
      _iOffset(_jOffset * lengthPerDim),
      _nElements(_iOffset * lengthPerDim) {
    _data = new T[_nElements];
    for (unsigned int i = 0; i < _nElements; i++) {
      _data[i] = data[i];
    }
  }
  /**
   * @brief Produces an empty tensor. Caution! Values are NOT initialized!
   * @param lengthPerDim number of entries into each dimension.
   */
  RegularRankFourTensor(const unsigned int& lengthPerDim)
    : _lengthPerDim(lengthPerDim),
      _kOffset(lengthPerDim),
      _jOffset(_kOffset * lengthPerDim),
      _iOffset(_jOffset * lengthPerDim),
      _nElements(_iOffset * lengthPerDim) {
    _data = new T[_nElements];
  }
  /**
   * @brief Produces a tensor with all members set to initValue
   * @param lengthPerDim number of entries into each dimension.
   * @param initValue is assigned to each field in the tensor.
   */
  RegularRankFourTensor(const unsigned int& lengthPerDim, T initValue)
    : _lengthPerDim(lengthPerDim),
      _kOffset(lengthPerDim),
      _jOffset(_kOffset * lengthPerDim),
      _iOffset(_jOffset * lengthPerDim),
      _nElements(_iOffset * lengthPerDim) {
    _data = new T[_nElements];
    for (unsigned int i = 0; i < _nElements; ++i) {
      _data[i] = initValue;
    }
  }

  virtual ~RegularRankFourTensor() {
    delete[] _data;
  }
  /// number of entries into each dimension.
  const unsigned int _lengthPerDim;
  /**
   * @brief Straightforward and rather safe access to the tensor. Slow!
   * @param i, j, k, l the index in each of the four dimensions.
   * @returns the value stored at this position in the tensor.
   *
   * For a faster access in time critical parts use the data() function. But be careful with
   * that, using that can more easily cause errors.
   */
  inline const T& operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const {
    assert(i < this->_lengthPerDim && j < this->_lengthPerDim && k < this->_lengthPerDim && l < this->_lengthPerDim);
    return this->_data[i * _iOffset + j * _jOffset + k * _kOffset + l];
  }
  /**
   * @brief Straightforward and rather safe access to the tensor. Slow!
   * @param i, j, k, l the index in each of the four dimensions.
   * @returns the value stored at this position in the tensor.
   *
   * For a faster access in time critical parts use the data() function. But be careful with
   * that, using that can more easily cause errors.
   */
  inline T& operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
    assert(i < this->_lengthPerDim && j < this->_lengthPerDim && k < this->_lengthPerDim && l < this->_lengthPerDim);
    return this->_data[i * _iOffset + j * _jOffset + k * _kOffset + l];
  }
  /**
   * @brief Prints the tensor to stdout. Careful! Will produce a lot of output!
   * TODO switch to managed internal output system
   */
  void print() const {
    for (unsigned int i = 0; i < this->_lengthPerDim; ++i) {
      for (unsigned int j = 0; j < this->_lengthPerDim; ++j) {
        for (unsigned int k = 0; k < this->_lengthPerDim; ++k) {
          for (unsigned int l = 0; l < this->_lengthPerDim; ++l) {
            std::cout << std::scientific << (*this)(i, j, k, l) << " ";
          }
          std::cout << std::endl;
        }
        std::cout << "------------------ j was " << j << std::endl;
      }
      std::cout << "------------------        i was " << i << std::endl;
    }
  }
  /**
   * @brief Caution! Only use if you know what you are doing!
   * @returns a pointer to the first data element.
   *
   * This provides the possibility for a fast but less robust access to the underlying data.
   * If you use this, you can orient by looking at the standard data access.
   */
  inline T* data() {
    return _data;
  }
  /**
   * @brief Caution! Only use if you know what you are doing!
   * @returns a pointer to the first data element.
   *
   * This provides the possibility for a fast but less robust access to the underlying data.
   * If you use this, you can orient by looking at the standard data access.
   */
  inline const T* data() const {
    return _data;
  }

  inline void cleardata() {
    delete[] _data;
  }
  /// @returns the offset in the data array for each increase in the index of the first dimension.
  inline unsigned int getIOffset() const {
    return _iOffset;
  }
  /// @returns the offset in the data array for each increase in the index of the second dimension.
  inline unsigned int getJOffset() const {
    return _jOffset;
  }
  /// @returns the offset in the data array for each increase in the index of the third dimension.
  inline unsigned int getKOffset() const {
    return _kOffset;
  }

 private:
  unsigned int _kOffset, _jOffset, _iOffset, _nElements;

  T* _data;
};

} /* namespace Serenity */
#endif /* REGULARRANKFOURTENSOR_H_ */
