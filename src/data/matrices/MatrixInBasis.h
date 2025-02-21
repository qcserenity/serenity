/**
 * @file   MatrixInBasis.h
 * @author Jan Unsleber, Thomas Dresselhaus
 * @date   rework on April 09. 2017
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
#ifndef MATRIXINBASIS_H
#define MATRIXINBASIS_H
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/SpinPolarizedData.h"
#include "data/matrices/SPMatrix.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward Declarations */
class BasisController;
/**
 * @brief Little helper function to easily assert that two objects are defined in the same basis
 * @param t
 * @param u
 * @returns True if t and u use the same BasisController, based on corresponding getter functions.
 */
template<class T, class U>
inline bool isDefinedInSameBasis(const T& t, const U& u) {
  return (t.getBasisController() == u.getBasisController());
}
/**
 * @class MatrixInBasis MatrixInBasis.h
 * @brief Data class for matrices defined within a basis set.
 *
 * Many operations in a quantum chemistry program are performed within a basis set which typically
 * consists of Gaussian (or Slater) functions. Most often, combinations of two basis functions
 * are used and data for each combination is stored in a matrix. Often, such data is obtained by
 * 'sandwiching' an operator between the basis functions, i.e.
 * \f${\rm Matrix}_{\mu, \nu} = \left< \chi_\mu | Op | \chi_\nu \right>\f$.
 * That kind of matrices are derived from this basic class. Take a look at e.g. the FockMatrix or
 * DensityMatrix.
 */
template<Options::SCF_MODES SCFMode>
class MatrixInBasis;
template<>
class MatrixInBasis<Options::SCF_MODES::RESTRICTED> : public SPMatrix<Options::SCF_MODES::RESTRICTED> {
 public:
  /**
   * @brief Constructor
   * @param basisController in its basis the matrix will be defined
   */
  inline MatrixInBasis<RESTRICTED>(std::shared_ptr<BasisController> basisController)
    : SPMatrix<RESTRICTED>(basisController->getNBasisFunctions(), basisController->getNBasisFunctions()),
      _basisController(basisController) {
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
  }
  inline MatrixInBasis<RESTRICTED>() = delete;
  /**
   * @brief Default destructor
   */
  virtual ~MatrixInBasis<RESTRICTED>() = default;
  /**
   * @brief Explicit copy constructor to avoid unintentional copying
   */
  inline MatrixInBasis<RESTRICTED>(const MatrixInBasis<RESTRICTED>& orig)
    : SPMatrix<RESTRICTED>(orig), _basisController(orig.getBasisController()) {
  }
  /**
   * @brief Move constructor
   */
  inline MatrixInBasis<RESTRICTED>(MatrixInBasis<RESTRICTED>&& orig)
    : SPMatrix<RESTRICTED>(orig), _basisController(orig.getBasisController()) {
  }
  /**
   * @brief Stores the MatrixInBasis into a HDF5 file.
   * @param fBaseName  The base name of the HDF5 file
   * @param matrixName The name of the dataset inside the file, optional.
   */
  void toHDF5(std::string fBaseName, std::string matrixName);
  /**
   * @returns The controller for the basis in which this matrix is defined.
   */
  inline std::shared_ptr<BasisController> getBasisController() const {
    return _basisController;
  }
  /**
   * @returns The dimension of the matrix, which equals the number of basis functions of the
   *          basisController.
   */
  inline unsigned int getNBasisFunctions() const {
    return _basisController->getNBasisFunctions();
  }
  /**
   * @brief Total data.
   * @return Returns the sum of alpha an beta..
   */
  inline MatrixInBasis<RESTRICTED> total() const {
    return *this;
  }
  /**
   * @brief Difference between alpha and beta.
   * @return Returns the difference of between alpha and beta (alpha-beta).
   */
  inline MatrixInBasis<RESTRICTED> difference() const {
    MatrixInBasis<RESTRICTED> ret(*this);
    ret.setZero();
    return ret;
  }
  /**
   * @brief Calculate the absolute max. coefficient for each shell-wise block and return the shell-wise matrix.
   * @return The shell-wise matrix.
   */
  SPMatrix<RESTRICTED> shellWiseAbsMax() const;

 private:
  /**
   * The controller of the basis the matrix is defined in
   */
  std::shared_ptr<BasisController> _basisController;

  /*
   * Operators explicitly needed in order not to lose
   * the BasisController while doing basic operations
   * with the underlying Eigen3 matrices.
   *
   * If you want to read up on this search for:
   * 'inheritance and slicing'
   *
   * - JU
   */
 public:
  using SPMatrix<RESTRICTED>::operator+=;
  using SPMatrix<RESTRICTED>::operator-=;

  // Operator overloads for other MatrixInBasis objects
  /// @brief Assignment operator.
  inline MatrixInBasis<RESTRICTED>& operator=(const MatrixInBasis<RESTRICTED>& other) {
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    if (other.getBasisController()) {
      if (_basisController != other.getBasisController())
        throw SerenityError("MatrixInBasis: BasisController do not match");
    }
    this->Base::operator=(other);
    return *this;
  }
  /// @brief Move assignment operator.
  inline MatrixInBasis<RESTRICTED>& operator=(MatrixInBasis<RESTRICTED>&& other) {
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    if (other.getBasisController()) {
      if (_basisController != other.getBasisController())
        throw SerenityError("MatrixInBasis: BasisController do not match");
    }
    this->Base::operator=(other);
    return *this;
  }

  inline void operator+=(const MatrixInBasis<RESTRICTED>& other) {
    if (_basisController != other.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    this->Base::operator+=(other);
  }
  inline MatrixInBasis<RESTRICTED> operator+(const MatrixInBasis<RESTRICTED>& y) const {
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    MatrixInBasis<RESTRICTED> result(*this);
    result += y;
    return result;
  }

  inline void operator-=(const MatrixInBasis<RESTRICTED>& other) {
    if (_basisController != other.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    this->Base::operator-=(other);
  }
  inline MatrixInBasis<RESTRICTED> operator-(const MatrixInBasis<RESTRICTED>& y) const {
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    MatrixInBasis<RESTRICTED> result(*this);
    result -= y;
    return result;
  }
  inline MatrixInBasis<RESTRICTED> operator*(const MatrixInBasis<RESTRICTED>& y) const {
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    MatrixInBasis<RESTRICTED> result(_basisController);
    result.Base::operator=(this->Base::operator*(y));
    return result;
  }
  inline MatrixInBasis<RESTRICTED> operator*(const double& y) const {
    MatrixInBasis<RESTRICTED> result(_basisController);
    result.Base::operator=(this->Base::operator*(y));
    return result;
  }

  template<typename OtherDerived>
  inline Eigen::MatrixXd operator+(const Eigen::MatrixBase<OtherDerived>& y) const {
    Eigen::MatrixXd result((Eigen::MatrixXd)(*this));
    result += y;
    return result;
  }
  template<typename OtherDerived>
  inline Eigen::MatrixXd operator-(const Eigen::MatrixBase<OtherDerived>& y) const {
    Eigen::MatrixXd result((Eigen::MatrixXd)(*this));
    result -= y;
    return result;
  }
  inline Eigen::MatrixXd operator*(const Eigen::MatrixXd& y) const {
    Eigen::MatrixXd result((Eigen::MatrixXd)(*this));
    result = this->Base::operator*(y);
    return result;
  }

  // Operator overloads for other Eigen3 objects
  /// @brief Assignment operator for Eigen3 objects.
  template<typename OtherDerived>
  inline MatrixInBasis<RESTRICTED>& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    this->Base::operator=(other);
    return *this;
  }
  /// @brief Move assignment operator for Eigen3 objects.
  template<typename OtherDerived>
  inline MatrixInBasis<RESTRICTED>& operator=(Eigen::MatrixBase<OtherDerived>&& other) {
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    this->Base::operator=(other);
    return *this;
  }

  /**
   * @brief Function to print MatrixInBasis<RESTRICTED> to e.g. std::cout .
   * @param stream The stream to print to.
   * @param matrix The matrix to be printed.
   * @return Returns the stream after usage.
   */
  friend std::ostream& operator<<(std::ostream& stream, const MatrixInBasis<RESTRICTED>& matrix) {
    stream << Eigen::MatrixXd(matrix);
    return stream;
  }
  /// @brief The SCF_MODE of this matrix
  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::RESTRICTED;
  /// @brief The type
  typedef MatrixInBasis<Options::SCF_MODES::RESTRICTED> type;
  /// @brief Base type for Eigen3.
  typedef Eigen::MatrixXd Base;
  /// @brief The underlying type
  typedef Eigen::MatrixXd& spinlesstype;
  /// @brief The underlying type
  typedef const Eigen::MatrixXd& constspinlesstype;

 protected:
  template<typename OtherDerived>
  inline MatrixInBasis<RESTRICTED>(const Eigen::MatrixBase<OtherDerived>& other)
    : SPMatrix<RESTRICTED>(other), _basisController(nullptr) {
  }
};

template<>
class MatrixInBasis<Options::SCF_MODES::UNRESTRICTED> : public SPMatrix<Options::SCF_MODES::UNRESTRICTED> {
 public:
  /**
   * @brief Constructor
   * @param basisController in its basis the matrix will be defined
   */
  inline MatrixInBasis<UNRESTRICTED>(std::shared_ptr<BasisController> basisController)
    : SPMatrix<UNRESTRICTED>(basisController->getNBasisFunctions(), basisController->getNBasisFunctions()),
      _basisController(basisController) {
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
  }

  inline MatrixInBasis<UNRESTRICTED>() = delete;
  /**
   * Default destructor
   */
  virtual ~MatrixInBasis<UNRESTRICTED>() = default;
  /**
   * @brief Copy constructor.
   */
  inline MatrixInBasis<UNRESTRICTED>(const MatrixInBasis<UNRESTRICTED>& orig)
    : SPMatrix<UNRESTRICTED>(orig), _basisController(orig.getBasisController()) {
  }
  inline MatrixInBasis<UNRESTRICTED>(const MatrixInBasis<RESTRICTED>& orig)
    : SPMatrix<UNRESTRICTED>(orig), _basisController(orig.getBasisController()) {
  }
  /**
   * @brief Move constructor
   */
  inline MatrixInBasis<UNRESTRICTED>(MatrixInBasis<UNRESTRICTED>&& orig)
    : SPMatrix<UNRESTRICTED>(orig), _basisController(orig.getBasisController()) {
    if (!orig.getBasisController())
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
  }
  /// @brief Assignment operator
  MatrixInBasis<UNRESTRICTED>& operator=(const MatrixInBasis<UNRESTRICTED>&) = default;
  /// @brief Move assignment operator
  MatrixInBasis<UNRESTRICTED>& operator=(MatrixInBasis<UNRESTRICTED>&&) = default;
  /**
   * @brief Stores the MatrixInBasis into a HDF5 file.
   * @param fBaseName  The base name of the HDF5 file
   * @param matrixName The name of the dataset inside the file, optional.
   */
  void toHDF5(std::string fBaseName, std::string matrixName);
  /**
   * @returns The controller for the basis in which this matrix is defined.
   */
  inline std::shared_ptr<BasisController> getBasisController() const {
    return _basisController;
  }
  /**
   * @returns The dimension of the matrix, which equals the number of basis functions of the
   *          basisController.
   */
  inline unsigned int getNBasisFunctions() const {
    return _basisController->getNBasisFunctions();
  }
  /**
   * @brief Total data.
   * @return Returns the sum of alpha an beta.
   */
  inline MatrixInBasis<RESTRICTED> total() const {
    MatrixInBasis<RESTRICTED> ret(_basisController);
    ret = this->alpha + this->beta;
    return ret;
  }
  /**
   * @brief Calculate the absolute max. coefficient for each shell-wise block and return the shell-wise matrix.
   * @return The shell-wise matrix.
   */
  SPMatrix<UNRESTRICTED> shellWiseAbsMax() const;
  /**
   * @brief Difference between alpha and beta.
   * @return Returns the difference of between alpha and beta (alpha-beta).
   */
  inline MatrixInBasis<RESTRICTED> difference() const {
    MatrixInBasis<RESTRICTED> ret(_basisController);
    ret = this->alpha - this->beta;
    return ret;
  }

 private:
  /**
   * The controller of the basis the matrix is defined in
   */
  std::shared_ptr<BasisController> _basisController;

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
  inline void operator+=(const MatrixInBasis<UNRESTRICTED>& y) {
    (*this).alpha += y.alpha;
    (*this).beta += y.beta;
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
  }
  inline void operator+=(const MatrixInBasis<RESTRICTED>& y) {
    (*this).alpha += y;
    (*this).beta += y;
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
  }
  inline MatrixInBasis<UNRESTRICTED> operator+(const MatrixInBasis<UNRESTRICTED>& y) const {
    MatrixInBasis<UNRESTRICTED> result(*this);
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    result.alpha += y.alpha;
    result.beta += y.beta;
    if (!result.getBasisController())
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    return result;
  }
  inline MatrixInBasis<UNRESTRICTED> operator+(const MatrixInBasis<RESTRICTED>& y) const {
    MatrixInBasis<UNRESTRICTED> result(*this);
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    result.alpha += y;
    result.beta += y;
    if (!result.getBasisController())
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    return result;
  }
  inline void operator-=(const MatrixInBasis<UNRESTRICTED>& y) {
    (*this).alpha -= y.alpha;
    (*this).beta -= y.beta;
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
  }
  inline void operator-=(const MatrixInBasis<RESTRICTED>& y) {
    (*this).alpha -= y;
    (*this).beta -= y;
    if (!_basisController)
      throw SerenityError("MatrixInBasis: Missing basis controller.");
  }
  inline MatrixInBasis<UNRESTRICTED> operator-(const MatrixInBasis<UNRESTRICTED>& y) const {
    MatrixInBasis<UNRESTRICTED> result(*this);
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    result.alpha -= y.alpha;
    result.beta -= y.beta;
    if (!result.getBasisController())
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    return result;
  }
  inline MatrixInBasis<UNRESTRICTED> operator-(const MatrixInBasis<RESTRICTED>& y) const {
    MatrixInBasis<UNRESTRICTED> result(*this);
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    result.alpha -= y;
    result.beta -= y;
    if (!result.getBasisController())
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    return result;
  }
  inline MatrixInBasis<UNRESTRICTED> operator*(const MatrixInBasis<UNRESTRICTED>& y) const {
    MatrixInBasis<UNRESTRICTED> result(_basisController);
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    result.alpha = (*this).alpha * y.alpha;
    result.beta = (*this).beta * y.beta;
    if (!result.getBasisController())
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    return result;
  }
  inline MatrixInBasis<UNRESTRICTED> operator*(const MatrixInBasis<RESTRICTED>& y) const {
    MatrixInBasis<UNRESTRICTED> result(_basisController);
    if (_basisController != y.getBasisController())
      throw SerenityError("MatrixInBasis: BasisController do not match");
    result.alpha = (*this).alpha * (Eigen::MatrixXd)y;
    result.beta = (*this).beta * (Eigen::MatrixXd)y;
    if (!result.getBasisController())
      throw SerenityError("MatrixInBasis: Missing basis controller.");
    return result;
  }

  friend std::ostream& operator<<(std::ostream& stream, const MatrixInBasis<UNRESTRICTED>& matrix) {
    stream << matrix.alpha << "\n \n" << matrix.beta;
    return stream;
  }
  /// @brief The SCF_MODE of this matrix
  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::UNRESTRICTED;
  /// @brief The type
  typedef MatrixInBasis<Options::SCF_MODES::UNRESTRICTED> type;
  /// @brief The underlying type
  typedef Eigen::MatrixXd& spinlesstype;
  /// @brief The underlying type
  typedef const Eigen::MatrixXd& constspinlesstype;
};

} /* namespace Serenity */
#endif /* MATRIXINBASIS_H */
