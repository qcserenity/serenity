/**
 * @file   GridData.h
 * @author Thomas Dresselhaus, Jan Unsleber
 *
 * @date   05.07.2015, last rework May 7, 2017 (JU)
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
#ifndef GRIDDATA_H
#define GRIDDATA_H
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "grid/GridController.h"
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declarations */
class Grid;
/**
 * @brief little helper function to easily assert that two objects are defined on the same integration grid
 *
 * This function should not only be used for data which are defined on an integration grid, but
 * also for functionalities which work with a certain grid.
 *
 * @param t
 * @param u
 * @returns true iff t and u use the same GridController, based on corresponding getter functions.
 */
template<class T, class U>
inline bool isDefinedOnSameGrid(const T& t, const U& u) {
  return (t.getGridController() == u.getGridController());
}

/**
 * @class GridData GridData.h
 * @brief Class to hold data that is defined on an integration grid. Usage: like std::vector
 *
 * Caution! Before any use you should make sure that the object is still in a valid state
 * (isValid()). An init() will put the object back into a valid state, but will throw away
 * all old data.\n
 *
 * Example usage (assuming you have a 'shared_ptr<GridController> gridController' available):\n
 *
 * // Construct\n
 * GridData myGridData(gridController);\n
 * // Access particular elements\n
 * myGridData[0] = 1.23;\n
 * myGridData[1] = 4.56;\n
 * // Loop over all elements\n
 * for (auto& dataPoint : myGridData) dataPoint *= 2;\n
 *
 * // At a point in time which may be long after the construction\n
 * // Either (e.g. if you only want to read data)\n
 * assert(myGridData.isValid());\n
 * // Or (if you recalculate all the data)\n
 * myGridData.init();
 */
template<Options::SCF_MODES SCFMode>
class GridData;

template<>
class GridData<Options::SCF_MODES::RESTRICTED> : public Eigen::VectorXd, public ObjectSensitiveClass<Grid> {
 public:
  /**
   * @brief Constructor
   * @param gridController The grid this data is defined on.
   */
  GridData(std::shared_ptr<GridController> gridController)
    : Eigen::VectorXd(gridController->getNGridPoints()), _gridController(gridController), _valid(true) {
    unsigned int nPoints = gridController->getNGridPoints();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->segment(start, n).setZero();
    }
    assert(_gridController);
  }

  virtual ~GridData<RESTRICTED>() = default;
  /**
   * @brief Explicit copy constructor to avoid unintentional copying
   */
  inline GridData<RESTRICTED>(const GridData<RESTRICTED>& orig)
    : Eigen::VectorXd(Eigen::VectorXd(orig)), _gridController(orig.getGridController()), _valid(true) {
  }
  /**
   * @brief Move constructor
   */
  inline GridData<RESTRICTED>(GridData<RESTRICTED>&& orig)
    : Eigen::VectorXd(Eigen::VectorXd(orig)), _gridController(orig.getGridController()), _valid(true) {
  }
  /**
   * @returns the controller for the basis in which this matrix is defined.
   */
  inline std::shared_ptr<GridController> getGridController() const {
    return _gridController;
  }
  /**
   * @returns the dimension of the matrix, which equals the number of basis functions of the
   *          GridController.
   */
  inline unsigned int getNGridPoints() const {
    return _gridController->getNGridPoints();
  }
  /**
   * @brief Total data.
   * @return Returns the sum of alpha an beta..
   */
  inline GridData<RESTRICTED> total() const {
    return *this;
  }
  /**
   * @brief Difference between alpha and beta.
   * @return Returns the difference between alpha and beta (alpha-beta).
   */
  inline GridData<RESTRICTED> difference() const {
    GridData<RESTRICTED> ret(*this);
    ret.setZero();
    return ret;
  }
  /**
   * @returns true iff the data still matches the grid on which it is defined.
   */
  inline bool isValid() const {
    return _valid;
  }
  /// @brief See ObjectSensitveClass.
  virtual void notify() {
    _valid = false;
  }

 private:
  std::shared_ptr<GridController> _gridController;

  /*
   * Operators explicitly needed in order not to loose
   * the GridController while doing basic operations
   * with the underlying Eigen3 matrices.
   *
   * If you want to read up on this search for:
   * 'inheritance and slicing'
   *
   * - JU
   */
 public:
  // Operator overloads for other GridData objects
  /// @brief Assignment operator.
  inline GridData<RESTRICTED>& operator=(const GridData<RESTRICTED>& other) {
    if (other.getGridController())
      assert(other.getGridController() == _gridController);
    this->Base::operator=(other);
    return *this;
  }
  /// @brief Move assignment operator.
  inline GridData<RESTRICTED>& operator=(GridData<RESTRICTED>&& other) {
    if (other.getGridController())
      assert(other.getGridController() == _gridController);
    this->Base::operator=(other);
    return *this;
  }
  inline void operator+=(const GridData<RESTRICTED>& other) {
    assert(_gridController == other.getGridController());
    unsigned int nPoints = this->size();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->segment(start, n) += other.segment(start, n);
    }
  }
  inline void operator+=(const Eigen::VectorXd& other) {
    unsigned int nPoints = this->size();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->segment(start, n) += other.segment(start, n);
    }
  }
  inline GridData<RESTRICTED> operator+(const GridData<RESTRICTED>& y) const {
    assert(_gridController == y.getGridController());
    GridData<RESTRICTED> result(*this);
    result += y;
    return result;
  }
  inline void operator-=(const Eigen::VectorXd& other) {
    unsigned int nPoints = this->size();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->segment(start, n) -= other.segment(start, n);
    }
  }
  inline void operator-=(const GridData<RESTRICTED>& other) {
    assert(_gridController == other.getGridController());
    unsigned int nPoints = this->size();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->segment(start, n) -= other.segment(start, n);
    }
  }
  inline GridData<RESTRICTED> operator-(const GridData<RESTRICTED>& y) const {
    assert(_gridController == y.getGridController());
    GridData<RESTRICTED> result(*this);
    result -= y;
    return result;
  }
  inline Eigen::VectorXd operator+(const Eigen::VectorXd& y) const {
    Eigen::VectorXd result((Eigen::VectorXd)(*this));
    result += y;
    return result;
  }
  inline Eigen::VectorXd operator-(const Eigen::VectorXd& y) const {
    Eigen::VectorXd result((Eigen::VectorXd)(*this));
    result -= y;
    return result;
  }

  // Operator overloads for other Eigen3 objects
  /// @brief Assignment operator for Eigen3 objects.
  template<typename OtherDerived>
  __attribute__((always_inline)) inline GridData<RESTRICTED>& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }
  /// @brief Move assignment operator for Eigen3 objects.
  template<typename OtherDerived>
  __attribute__((always_inline)) inline GridData<RESTRICTED>& operator=(Eigen::MatrixBase<OtherDerived>&& other) {
    this->Base::operator=(other);
    return *this;
  }

  /**
   * @brief Function to print GridData<RESTRICTED> to e.g. std::cout .
   * @param stream The stream to print to.
   * @param matrix The matrix to be printed.
   * @return Returns the stream after usage.
   */
  friend std::ostream& operator<<(std::ostream& stream, const GridData<RESTRICTED>& vector) {
    stream << Eigen::VectorXd(vector);
    return stream;
  }
  /// @brief The SCF_MODE of this matrix
  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::RESTRICTED;
  /// @brief The type
  typedef Eigen::VectorXd type;
  /// @brief Base type for Eigen3.
  typedef Eigen::VectorXd Base;
  /// @brief The underlying type
  typedef Eigen::VectorXd& spinlesstype;
  /// @brief The underlying type
  typedef const Eigen::VectorXd& constspinlesstype;
  /// @brief A switch to determine if the data is still valid, or if the Grid changed.
  bool _valid;

 protected:
  template<typename OtherDerived>
  __attribute__((always_inline)) inline GridData<RESTRICTED>(const Eigen::MatrixBase<OtherDerived>& other)
    : Eigen::VectorXd(other), _gridController(nullptr), _valid(true) {
  }
};

template<>
class GridData<Options::SCF_MODES::UNRESTRICTED> : public ObjectSensitiveClass<Grid> {
 public:
  /**
   * @brief Constructor
   * @param GridController The grid this data is defined on.
   */
  GridData(std::shared_ptr<GridController> gridController)
    : alpha(gridController->getNGridPoints()),
      beta(gridController->getNGridPoints()),
      _gridController(gridController),
      _valid(true) {
    unsigned int nPoints = gridController->getNGridPoints();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->alpha.segment(start, n).setZero();
      this->beta.segment(start, n).setZero();
    }
    assert(_gridController);
  }
  virtual ~GridData<UNRESTRICTED>() = default;
  /**
   * @brief Explicit copy constructor to avoid unintentional copying
   */
  explicit inline GridData<UNRESTRICTED>(const GridData<UNRESTRICTED>& orig)
    : alpha(orig.alpha), beta(orig.beta), _gridController(orig.getGridController()), _valid(true) {
  }
  inline GridData<UNRESTRICTED>(const GridData<RESTRICTED>& orig)
    : alpha(orig), beta(orig), _gridController(orig.getGridController()), _valid(true) {
  }
  /**
   * @brief Move constructor
   */
  inline GridData<UNRESTRICTED>(GridData<UNRESTRICTED>&& orig)
    : alpha((Eigen::VectorXd)orig.alpha),
      beta((Eigen::VectorXd)orig.beta),
      _gridController(orig.getGridController()),
      _valid(true) {
    assert(orig.getGridController());
    assert(_gridController);
  }
  /// @brief assignment operator
  inline GridData<UNRESTRICTED>& operator=(const GridData<UNRESTRICTED>& orig) {
    this->alpha = orig.alpha;
    this->beta = orig.beta;
    this->_gridController = orig.getGridController();
    this->_valid = orig.isValid();
    return *this;
  };
  /// @brief move assignment operator
  inline GridData<UNRESTRICTED>& operator=(GridData<UNRESTRICTED>&& orig) {
    this->alpha = std::move(orig.alpha);
    this->beta = std::move(orig.beta);
    this->_gridController = orig.getGridController();
    this->_valid = orig.isValid();
    return *this;
  }

  /**
   * @returns the controller for the basis in which this matrix is defined.
   */
  inline std::shared_ptr<GridController> getGridController() const {
    return _gridController;
  }
  /**
   * @returns the dimension of the matrix, which equals the number of basis functions of the
   *          GridController.
   */
  inline unsigned int getNGridPoints() const {
    return _gridController->getNGridPoints();
  }
  /**
   * @brief Total data.
   * @return Returns the sum of alpha an beta..
   */
  inline GridData<RESTRICTED> total() const {
    GridData<RESTRICTED> ret(_gridController);
    ret = this->alpha + this->beta;
    return ret;
  }
  /**
   * @brief Difference between alpha and beta.
   * @return Returns the difference of between alpha and beta (alpha-beta).
   */
  inline GridData<RESTRICTED> difference() const {
    GridData<RESTRICTED> ret(_gridController);
    ret = this->alpha - this->beta;
    return ret;
  }
  /**
   * @returns true iff the data still matches the grid on which it is defined.
   */
  inline bool isValid() const {
    return _valid;
  }
  /// @brief See ObjectSensitveClass.
  virtual void notify() {
    _valid = false;
  }
  /// @brief The alpha part.
  Eigen::VectorXd alpha;
  /// @brief The beta part.
  Eigen::VectorXd beta;

 private:
  std::shared_ptr<GridController> _gridController;

  /*
   * Operators explicitly needed in order not to loose
   * the GridController while doing basic operations
   * with the underlying Eigen3 matrices.
   *
   * If you want to read up on this search for:
   * 'inheritance and slicing'
   *
   * - JU
   */
 public:
  inline void operator+=(const GridData<UNRESTRICTED>& other) {
    unsigned int nPoints = this->alpha.size();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->alpha.segment(start, n) += other.alpha.segment(start, n);
      this->beta.segment(start, n) += other.beta.segment(start, n);
    }
    assert(_gridController);
  }
  inline void operator+=(const GridData<RESTRICTED>& other) {
    unsigned int nPoints = this->alpha.size();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->alpha.segment(start, n) += other.segment(start, n);
      this->beta.segment(start, n) += other.segment(start, n);
    }
    assert(_gridController);
  }
  inline GridData<UNRESTRICTED> operator+(const GridData<UNRESTRICTED>& y) const {
    GridData<UNRESTRICTED> result(*this);
    assert(_gridController == y.getGridController());
    result.alpha += y.alpha;
    result.beta += y.beta;
    assert(result.getGridController());
    return result;
  }
  inline GridData<UNRESTRICTED> operator+(const GridData<RESTRICTED>& y) const {
    GridData<UNRESTRICTED> result(*this);
    assert(_gridController == y.getGridController());
    result.alpha += y;
    result.beta += y;
    assert(result.getGridController());
    return result;
  }
  inline void operator-=(const GridData<UNRESTRICTED>& other) {
    unsigned int nPoints = this->alpha.size();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->alpha.segment(start, n) -= other.alpha.segment(start, n);
      this->beta.segment(start, n) -= other.beta.segment(start, n);
    }
    assert(_gridController);
  }
  inline void operator-=(const GridData<RESTRICTED>& other) {
    unsigned int nPoints = this->alpha.size();
    unsigned int nBlocks = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
      unsigned int n = (unsigned int)(nPoints / nBlocks);
      const unsigned int start = iBlock * n;
      if (iBlock == nBlocks - 1)
        n += nPoints % nBlocks;
      this->alpha.segment(start, n) -= other.segment(start, n);
      this->beta.segment(start, n) -= other.segment(start, n);
    }
    assert(_gridController);
  }
  inline GridData<UNRESTRICTED> operator-(const GridData<UNRESTRICTED>& y) const {
    GridData<UNRESTRICTED> result(*this);
    assert(_gridController == y.getGridController());
    result.alpha -= y.alpha;
    result.beta -= y.beta;
    assert(result.getGridController());
    return result;
  }
  inline GridData<UNRESTRICTED> operator-(const GridData<RESTRICTED>& y) const {
    GridData<UNRESTRICTED> result(*this);
    assert(_gridController == y.getGridController());
    result.alpha -= y;
    result.beta -= y;
    assert(result.getGridController());
    return result;
  }

  friend std::ostream& operator<<(std::ostream& stream, const GridData<UNRESTRICTED>& matrix) {
    stream << matrix.alpha << "\n \n" << matrix.beta;
    return stream;
  }
  /// @brief The SCF_MODE of this matrix
  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::UNRESTRICTED;
  /// @brief The type
  typedef GridData<Options::SCF_MODES::UNRESTRICTED> type;
  /// @brief The underlying type
  typedef Eigen::VectorXd& spinlesstype;
  /// @brief The underlying type
  typedef const Eigen::VectorXd& constspinlesstype;
  /// @brief A switch to determine if the data is still valid, or if the Grid changed.
  bool _valid;
};

} /* namespace Serenity */

#endif /* GRIDDATA_H */
