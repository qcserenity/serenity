/**
 * @file   SPMatrix.h
 * @author M. Boeckers
 * @date   April 20. 2018
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
#ifndef SPMATRIX_H
#define SPMATRIX_H
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class SPMatrix;
template<>
class SPMatrix<Options::SCF_MODES::RESTRICTED> : public Eigen::MatrixXd {
 public:
  /**
   * @brief Constructor.
   * @param nRows Number of rows of the matrix.
   * @param nCols Number of columns of rows of the matrix.
   */
  inline SPMatrix(unsigned int nRows, unsigned int nCols) : Eigen::MatrixXd(nRows, nCols) {
    this->Base::setZero();
  }

  virtual ~SPMatrix<RESTRICTED>() = default;
  /**
   * @brief Explicit copy constructor to avoid unintentional copying
   */
  inline SPMatrix<RESTRICTED>(const SPMatrix<RESTRICTED>& orig) : Eigen::MatrixXd(Eigen::MatrixXd(orig)) {
  }
  /**
   * @brief Constructor for the use in vector initializations etc.
   */
  inline SPMatrix<RESTRICTED>() : Eigen::MatrixXd(Eigen::MatrixXd()) {
  }
  /**
   * @brief Move constructor
   */
  inline SPMatrix<RESTRICTED>(SPMatrix<RESTRICTED>&& orig) : Eigen::MatrixXd(Eigen::MatrixXd(orig)) {
  }
  /**
   * @brief Stores the SPMatrix into a HDF5 file.
   * @param fBaseName  The base name of the HDF5 file
   * @param matrixName The name of the dataset inside the file, optional.
   *
   */
  void toHDF5(std::string fBaseName, std::string matrixName);
  /**
   * @brief Total data.
   * @return Returns the sum of alpha an beta..
   */
  inline SPMatrix<RESTRICTED> total() const {
    return *this;
  }
  /**
   * @brief Difference between alpha and beta.
   * @return Returns the difference of between alpha and beta (alpha-beta).
   */
  inline SPMatrix<RESTRICTED> difference() const {
    SPMatrix<RESTRICTED> ret(*this);
    ret.setZero();
    return ret;
  }

  /*
   * Operators explicitly needed in order not to loose
   * the BasisController while doing basic operations
   * with the underlying Eigen3 matrices.
   *
   * If you want to read up on this search for:
   * 'inheritance and slicing'
   *
   * - JU
   */
 public:
  template<typename OtherDerived>
  __attribute__((always_inline)) inline SPMatrix<RESTRICTED>(const Eigen::MatrixBase<OtherDerived>& other)
    : Eigen::MatrixXd(other) {
  }

  // Operator overloads for other SPMatrix objects
  /// @brief Assignment operator.
  __attribute__((always_inline)) inline SPMatrix<RESTRICTED>& operator=(const SPMatrix<RESTRICTED>& other) {
    this->Base::operator=(other);
    return *this;
  }
  /// @brief Move assignment operator.
  __attribute__((always_inline)) inline SPMatrix<RESTRICTED>& operator=(SPMatrix<RESTRICTED>&& other) {
    this->Base::operator=(other);
    return *this;
  }
  __attribute__((always_inline)) inline void operator+=(const SPMatrix<RESTRICTED>& other) {
    this->Base::operator+=(other);
  }
  __attribute__((always_inline)) inline SPMatrix<RESTRICTED> operator+(const SPMatrix<RESTRICTED>& y) const {
    SPMatrix<RESTRICTED> result(*this);
    result += y;
    return result;
  }
  using Eigen::MatrixXd::operator-=;
  using Eigen::MatrixXd::operator+=;
  __attribute__((always_inline)) inline void operator-=(const SPMatrix<RESTRICTED>& other) {
    this->Base::operator-=(other);
  }
  __attribute__((always_inline)) inline SPMatrix<RESTRICTED> operator-(const SPMatrix<RESTRICTED>& y) const {
    SPMatrix<RESTRICTED> result(*this);
    result -= y;
    return result;
  }
  __attribute__((always_inline)) inline SPMatrix<RESTRICTED> operator*(const SPMatrix<RESTRICTED>& y) const {
    SPMatrix<RESTRICTED> result(this->rows(), this->cols());
    result.Base::operator=(this->Base::operator*(y));
    return result;
  }
  __attribute__((always_inline)) inline Eigen::MatrixXd operator+(const Eigen::MatrixXd& y) const {
    Eigen::MatrixXd result((Eigen::MatrixXd)(*this));
    result += y;
    return result;
  }
  __attribute__((always_inline)) inline Eigen::MatrixXd operator-(const Eigen::MatrixXd& y) const {
    Eigen::MatrixXd result((Eigen::MatrixXd)(*this));
    result -= y;
    return result;
  }
  __attribute__((always_inline)) inline Eigen::MatrixXd operator*(const Eigen::MatrixXd& y) const {
    Eigen::MatrixXd result((Eigen::MatrixXd)(*this));
    result = this->Base::operator*(y);
    return result;
  }

  // Operator overloads for other Eigen3 objects
  /// @brief Assignment operator for Eigen3 objects.
  template<typename OtherDerived>
  __attribute__((always_inline)) inline SPMatrix<RESTRICTED>& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }
  /// @brief Move assignment operator for Eigen3 objects.
  template<typename OtherDerived>
  __attribute__((always_inline)) inline SPMatrix<RESTRICTED>& operator=(Eigen::MatrixBase<OtherDerived>&& other) {
    this->Base::operator=(other);
    return *this;
  }

  /**
   * @brief Function to print SPMatrix<RESTRICTED> to e.g. std::cout .
   * @param stream The stream to print to.
   * @param matrix The matrix to be printed.
   * @return Returns the stream after usage.
   */
  friend std::ostream& operator<<(std::ostream& stream, const SPMatrix<RESTRICTED>& matrix) {
    stream << Eigen::MatrixXd(matrix);
    return stream;
  }
  /// @brief The SCF_MODE of this matrix
  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::RESTRICTED;
  /// @brief The type
  typedef Eigen::MatrixXd type;
  /// @brief Base type for Eigen3.
  typedef Eigen::MatrixXd Base;
  /// @brief The underlying type
  typedef Eigen::MatrixXd& spinlesstype;
  /// @brief The underlying type
  typedef const Eigen::MatrixXd& constspinlesstype;
};

template<>
class SPMatrix<Options::SCF_MODES::UNRESTRICTED> {
 public:
  /**
   * @brief Constructor
   * @param basisController in its basis the matrix will be defined
   */
  inline SPMatrix(unsigned int nRows, unsigned int nCols) : alpha(nRows, nCols), beta(nRows, nCols) {
    alpha.setZero();
    beta.setZero();
  }
  /**
   * @brief Constructor for the use in vector initializations etc.
   */
  inline SPMatrix<UNRESTRICTED>() : alpha(Eigen::MatrixXd()), beta(Eigen::MatrixXd()) {
  }

  virtual ~SPMatrix<UNRESTRICTED>() = default;
  /**
   * @brief Copy constructor.
   */
  inline SPMatrix<UNRESTRICTED>(const SPMatrix<UNRESTRICTED>& orig) : alpha(orig.alpha), beta(orig.beta) {
  }
  inline SPMatrix<UNRESTRICTED>(const SPMatrix<RESTRICTED>& orig) : alpha(orig), beta(orig) {
  }
  /**
   * @brief Move constructor
   */
  inline SPMatrix<UNRESTRICTED>(SPMatrix<UNRESTRICTED>&& orig)
    : alpha((Eigen::MatrixXd)orig.alpha), beta((Eigen::MatrixXd)orig.beta) {
  }
  /// @brief assignment operator
  SPMatrix<UNRESTRICTED>& operator=(const SPMatrix<UNRESTRICTED>&) = default;
  /// @brief move assignment operator
  SPMatrix<UNRESTRICTED>& operator=(SPMatrix<UNRESTRICTED>&&) = default;
  /**
   * @brief Stores the SPMatrix into a HDF5 file.
   * @param fBaseName  The base name of the HDF5 file
   * @param matrixName The name of the dataset inside the file, optional.
   *
   */
  void toHDF5(std::string fBaseName, std::string matrixName);
  /**
   * @brief Total data.
   * @return Returns the sum of alpha an beta..
   */
  inline SPMatrix<RESTRICTED> total() const {
    SPMatrix<RESTRICTED> ret(this->alpha.rows(), this->alpha.cols());
    ret = this->alpha + this->beta;
    return ret;
  }
  /**
   * @brief Difference between alpha and beta.
   * @return Returns the difference of between alpha and beta (alpha-beta).
   */
  inline SPMatrix<RESTRICTED> difference() const {
    SPMatrix<RESTRICTED> ret(this->alpha.rows(), this->alpha.cols());
    ret = this->alpha - this->beta;
    return ret;
  }

  /**
   * @brief Returns the number of rows of the matrix.
   */
  inline unsigned int rows() const {
    assert(alpha.rows() == beta.rows());
    return alpha.rows();
  }

  /**
   * @brief Returns the number of columns of the matrix.
   */
  inline unsigned int cols() const {
    assert(alpha.cols() == beta.cols());
    return alpha.cols();
  }

  /**
   * @brief Returns the number of entries in the alpha matrix.
   */
  inline unsigned int size() const {
    return alpha.size();
  }

  /// @brief The alpha part.
  Eigen::MatrixXd alpha;
  /// @brief The beta part.
  Eigen::MatrixXd beta;

  /*
   * Operators explicitly needed in order not to loose
   * the BasisController while doing basic operations
   * with the underlying Eigen3 matrices.
   *
   * If you want to read up on this search for:
   * 'inheritance and slicing'
   *
   * - JU
   */
 public:
  inline void operator+=(const SPMatrix<UNRESTRICTED>& y) {
    (*this).alpha += y.alpha;
    (*this).beta += y.beta;
  }
  inline void operator+=(const SPMatrix<RESTRICTED>& y) {
    (*this).alpha += y;
    (*this).beta += y;
  }
  inline SPMatrix<UNRESTRICTED> operator+(const SPMatrix<UNRESTRICTED>& y) const {
    SPMatrix<UNRESTRICTED> result(*this);
    result.alpha += y.alpha;
    result.beta += y.beta;
    return result;
  }
  inline SPMatrix<UNRESTRICTED> operator+(const SPMatrix<RESTRICTED>& y) const {
    SPMatrix<UNRESTRICTED> result(*this);
    result.alpha += y;
    result.beta += y;
    return result;
  }
  inline void operator-=(const SPMatrix<UNRESTRICTED>& y) {
    (*this).alpha -= y.alpha;
    (*this).beta -= y.beta;
  }
  inline void operator-=(const SPMatrix<RESTRICTED>& y) {
    (*this).alpha -= y;
    (*this).beta -= y;
  }
  inline SPMatrix<UNRESTRICTED> operator-(const SPMatrix<UNRESTRICTED>& y) const {
    SPMatrix<UNRESTRICTED> result(*this);
    result.alpha -= y.alpha;
    result.beta -= y.beta;
    return result;
  }
  inline SPMatrix<UNRESTRICTED> operator-(const SPMatrix<RESTRICTED>& y) const {
    SPMatrix<UNRESTRICTED> result(*this);
    result.alpha -= y;
    result.beta -= y;
    return result;
  }
  inline SPMatrix<UNRESTRICTED> operator*(const SPMatrix<UNRESTRICTED>& y) const {
    SPMatrix<UNRESTRICTED> result(this->alpha.rows(), this->alpha.cols());
    result.alpha = (*this).alpha * y.alpha;
    result.beta = (*this).beta * y.beta;
    return result;
  }
  inline SPMatrix<UNRESTRICTED> operator*(const SPMatrix<RESTRICTED>& y) const {
    SPMatrix<UNRESTRICTED> result(this->alpha.rows(), this->alpha.cols());
    result.alpha = (*this).alpha * (Eigen::MatrixXd)y;
    result.beta = (*this).beta * (Eigen::MatrixXd)y;
    return result;
  }

  friend std::ostream& operator<<(std::ostream& stream, const SPMatrix<UNRESTRICTED>& matrix) {
    stream << matrix.alpha << "\n \n" << matrix.beta;
    return stream;
  }
  /// @brief The SCF_MODE of this matrix
  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::UNRESTRICTED;
  /// @brief The type
  typedef SPMatrix<Options::SCF_MODES::UNRESTRICTED> type;
  /// @brief The underlying type
  typedef Eigen::MatrixXd& spinlesstype;
  /// @brief The underlying type
  typedef const Eigen::MatrixXd& constspinlesstype;
};

} /* namespace Serenity */
#endif /* SPMATRIX_H */
