/**
 * @file   ScalarProductMatrix.h
 *
 * @date   Jun 11, 2014
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
#ifndef SCALARPRODUCTMATRIX_H_
#define SCALARPRODUCTMATRIX_H_
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
template<class T>
class Matrix;
/**
 * @class ScalarProductMatrix ScalarProductMatrix.h
 * @brief Holds and manages a matrix of varying size consisting of scalar products
 *
 */
class ScalarProductMatrix {
 public:
  /**
   * @param dataLength of each data set which will be attached. Memory will be reserved accordingly.
   */
  ScalarProductMatrix(unsigned int rows, unsigned int cols);
  virtual ~ScalarProductMatrix();
  /**
   * @param newData of size dataLength (see constructor)
   */
  void putNewData(Matrix<double> newData);
  /**
   * Erases all data stored by this object
   */
  void eraseAllData();
  /**
   * @returns the current matrix which holds the scalar products of all currently attached data sets
   */
  inline const Matrix<double>& getScalarProductMatrix() const {
    return *_scalarProductMatrix;
  }
  /**
   * @returns a vector holding all stored data sets
   */
  inline const std::vector<Matrix<double>>& getRawData() const {
    return _rawData;
  }

 private:
  const unsigned int _rows;
  const unsigned int _cols;

  std::unique_ptr<Matrix<double>> _scalarProductMatrix;
  /**
   * Scalar products are constructed from the underlying data of pairs of held pointers
   */
  std::vector<Matrix<double>> _rawData;
  unsigned int _nStored;
};

} /* namespace Serenity */

#endif /* SCALARPRODUCTMATRIX_H_ */
