/**
 * @file   Matrix.h
 *
 * @date   28. Juli 2013, 20:22
 * @author Thomas Dresselhaus
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef MATRIX_H
#define MATRIX_H
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <vector>

namespace Serenity {
/**
 * @class Matrix Matrix.h
 * @brief an object holding numbers in rows and columns. Access via foo(i,j)
 *
 * with i being the row index and j the column index.
 */
template<typename T>
class Matrix {
 public:
  /**
   * @brief Provides a filled nRows x nCols matrix.
   * @param nRows the number of rows
   * @param nCols the number of columns
   * @param init initial value
   */
  Matrix(const int& nRows, const int& nCols, T init)
    : _nRows(nRows), _nColumns(nCols), _size(nRows * nCols), _data(nRows * nCols, init) {
  }

  /**
   * @brief Provides a filled nRows x nCols matrix.
   * @param nRows the number of rows
   * @param nCols the number of columns
   * @param data  the matrix will be constructed using this data
   */
  Matrix(const int& nRows, const int& nCols, std::vector<T>&& data)
    : _nRows(nRows), _nColumns(nCols), _size(nRows * nCols), _data(std::move(data)) {
  }

  /**
   * @brief Straightforward and safe access to the matrix. Slow!
   * @param i the row index
   * @param j the column index
   * @returns the value stored at this position in the matrix.
   *
   * For a faster access in time critical parts use the data() function. But be careful with
   * that, using that can more easily cause errors.
   */
  inline const T& operator()(unsigned int i, unsigned int j) const {
    assert(i < this->_nRows && j < this->_nColumns);
    return this->_data[j * this->_nRows + i];
  }
  /**
   * @brief Straightforward and safe access to the matrix. Slow!
   * @param i the row index
   * @param j the column index
   * @returns the value stored at this position in the matrix.
   *
   * For a faster access in time critical parts use the data() function. But be careful with
   * that, using that can more easily cause errors.
   */
  inline T& operator()(unsigned int i, unsigned int j) {
    assert(i < this->_nRows && j < this->_nColumns);
    return this->_data[j * this->_nRows + i];
  }
  /**
   * @brief Prints the matrix to stdout.
   * TODO switch to managed internal output system
   */
  void print() const {
    for (unsigned int i = 0; i < this->_nRows; i++) {
      for (unsigned int j = 0; j < this->_nColumns; j++) {
        std::cout << std::scientific << (*this)(i, j) << " ";
      }
      std::cout << std::endl;
    }
  }
  /**
   * @returns the underlying data structure (array)
   *
   * Caution! Only use this if you are totally sure what you
   * are doing. The structure of the array might not be the
   * way you expect! Suggestion: use the () operators.
   */
  inline T* data() {
    return _data.data();
  }
  /**
   * @returns the underlying data structure (array)
   *
   * Caution! Only use this if you are totally sure what you
   * are doing. The structure of the array might not be the
   * way you expect! Suggestion: use the () operators.
   */
  inline const T* data() const {
    return _data.data();
  }

  /**
   * @brief returns a specific row of the matrix
   * @param iRow the number of the row that should be returned
   * @return a vector containing all elements of the row
   */
  std::vector<T> row(unsigned int iRow) const {
    assert(iRow < _nRows);
    std::vector<T> returnRow(_nColumns);
    for (unsigned int i = 0; i < _nColumns; ++i) {
      returnRow[i] = _data[i * _nRows + iRow];
    }
    return returnRow;
  }

  /**
   * @brief returns a specific column of the matrix
   * @param iCol the number of the column that should be returned
   * @return a vector containing all elements of the column
   */
  std::vector<T> col(unsigned int iCol) const {
    assert(iCol < _nColumns);
    std::vector<T> returnCol(_nRows);
    for (unsigned int i = 0; i < _nRows; i++) {
      returnCol[i] = _data[iCol * _nRows + i];
    }
    return returnCol;
  }

  /**
   * @return Returns the number of elements in the matrix.
   */
  inline unsigned int size() const {
    return _size;
  }
  /**
   * @return Returns the number of columns.
   */
  inline unsigned int cols() const {
    return _nColumns;
  }
  /**
   * @return Returns the number of rows.
   */
  inline unsigned int rows() const {
    return _nRows;
  }

 private:
  ///@brief The rank, if viewed as a tensor
  const static int RANK = 2;
  ///@{Size indicators for the matrix
  unsigned int _nRows, _nColumns, _size;
  ///}
 protected:
  /// the underlying data object
  std::vector<T> _data;
};

/*
 * Eigen3 versions for Matrix<double> and Matrix<int>
 *
 * For the complete documentation check the Eigen3 website:
 *  http://eigen.tuxfamily.org/dox/
 */
template<>
class Matrix<double> : public Eigen::MatrixXd {
  using Eigen::MatrixXd::MatrixXd;

 public:
  /**
   * @brief Prints the matrix to stdout.
   * TODO switch to managed internal output system
   */
  void print() const {
    for (unsigned int i = 0; i < this->rows(); i++) {
      for (unsigned int j = 0; j < this->cols(); j++) {
        std::cout << std::scientific << (*this)(i, j) << " ";
      }
      std::cout << std::endl;
    }
  }

 private:
  ///@brief The rank, if viewed as a tensor
  const static int RANK = 2;
};

template<>
class Matrix<int> : public Eigen::MatrixXi {
  using Eigen::MatrixXi::MatrixXi;

 public:
  /**
   * @brief Prints the matrix to stdout.
   * TODO switch to managed internal output system
   */
  void print() const {
    for (unsigned int i = 0; i < this->rows(); i++) {
      for (unsigned int j = 0; j < this->cols(); j++) {
        std::cout << std::scientific << (*this)(i, j) << " ";
      }
      std::cout << std::endl;
    }
  }

 private:
  ///@brief The rank, if viewed as a tensor
  const static int RANK = 2;
};

} /* namespace Serenity */
#endif /* MATRIX_H */
