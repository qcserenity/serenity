/**
 * @file   SpinPolarizedData.h
 * @author Jan Unsleber, Thomas Dresselhaus
 * @date   rework on April 09. 2017
 *
 * April 09. 2017: defined the for_spin macro a bit more broadly. -JU
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
#ifndef SPINPOLARIZEDDATA_H_
#define SPINPOLARIZEDDATA_H_
/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <functional>
#include <type_traits>
#include <utility>

namespace Serenity {
/* Forward Declarations */
/**
 * @class SpinPolarizedData SpinPolarizedData.h
 * @brief General version of objects which may contain one or two sets of data.
 *
 * The intention is that often a restricted and an unrestricted version of an object is needed,
 * where the unrestricted version is typically just twice the restricted version, once for alpha
 * and once for beta spin (access via yourData.alpha and yourData.beta).
 *
 * Technical information (Sept 11, 2017 - JU):
 * SpinPolarizedData has three template parameters.
 * The first one is the SCF_MODE determining if there is one set of data stored in the
 *   restricted case or if there are two sets (alpha and beta) in case of an unrestricted
 *   mode.
 * The second one is the class of the underlying data this might be an integer or double for
 *   something like a number of electrons or it might be an Eigen::VectorXd in case of orbital
 *   populations.
 * The third template parameter is a purely technical one, in case of a restricted data set,
 *   the class needs different specializations for the second parameter being a class and
 *   for it being a primitive.
 *   The last template parameter is helping to ensure that these cases are handled as expected.
 *   For further information see e.g.:
 *   https://stackoverflow.com/questions/11287043/is-there-a-way-to-specialize-a-template-to-target-primitives
 */
template<Options::SCF_MODES SCFMode, class U, typename E = void>
class SpinPolarizedData;

/* =================================================
 *   Specialization for RESTRICTED, non-primitives
 * =================================================*/
template<class U>
class SpinPolarizedData<RESTRICTED, U, typename std::enable_if<std::is_class<U>::value>::type> : public U {
 public:
  using U::U;

  /// @brief Default destructor.
  virtual ~SpinPolarizedData() = default;

  /// @brief Default constructor.
  inline SpinPolarizedData() : U() {
  }
  /**
   * @brief Copy constructor using base object.
   * @param init Initial value.
   */
  inline SpinPolarizedData(const U& init) : U(init) {
  }
  /**
   * @brief Move constructor using base object.
   * @param init Initial value.
   */
  inline SpinPolarizedData(U&& init) : U(std::move(init)) {
  }
  /**
   * @brief Copy constructor using SpinPolarizedData.
   * @param orig Initial/original value.
   */
  inline SpinPolarizedData(const SpinPolarizedData<RESTRICTED, U, void>& orig) : U(orig) {
  }
  /**
   * @brief Copy constructor using SpinPolarizedData.
   * @param orig Initial/original value.
   */
  inline SpinPolarizedData(SpinPolarizedData<RESTRICTED, U, void>& orig) : U(orig) {
  }
  /**
   * @brief Move constructor using SpinPolarizedData.
   * @param orig Initial/original value.
   */
  inline SpinPolarizedData(SpinPolarizedData<RESTRICTED, U, void>&& orig) : U(std::move(orig)) {
  }
  /**
   * @brief Assignment operator.
   * @param rhs Other object.
   */
  SpinPolarizedData<RESTRICTED, U, void>& operator=(const SpinPolarizedData<RESTRICTED, U, void>& rhs) {
    if (this != &rhs)
      this->U::operator=(rhs);
    return *this;
  }
  /**
   * @brief Move assignment operator.
   * @param rhs Other object.
   */
  SpinPolarizedData<RESTRICTED, U, void>& operator=(SpinPolarizedData<RESTRICTED, U, void>&& rhs) {
    if (this != &rhs)
      this->U::operator=(std::move(rhs));
    return *this;
  }

  /**
   * @brief Function to get total value (see unrestricted implementations).
   * @return Returns the total value.
   */
  inline U total() {
    return *this;
  }

  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::RESTRICTED;
  /// @brief The type
  typedef U type;
  /// @brief The underlying type
  typedef U& spinlesstype;
  /// @brief The underlying type
  typedef const U& constspinlesstype;
};

/* =============================================
 *   Specialization for Restricted, primitives
 * =============================================*/
template<class U>
class SpinPolarizedData<RESTRICTED, U, typename std::enable_if<!std::is_class<U>::value>::type> {
 public:
  // Check whether the constructors have to be explicit (if that is possible at all in a sensible way)
  inline SpinPolarizedData(){};
  inline SpinPolarizedData(const U& initValue) : _m(initValue){};
  inline operator const U&() const {
    return _m;
  };
  inline operator U&() {
    return _m;
  };

  inline U total() {
    return _m;
  }

  U _m;
  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::RESTRICTED;
  /// @brief The type
  typedef U type;
  /// @brief The underlying type
  typedef U& spinlesstype;
  /// @brief The underlying type
  typedef const U& constspinlesstype;
};

/* ================================================
 *   Specialization for UNRESTRICTED, (all cases)
 * ================================================*/
template<class U>
class SpinPolarizedData<UNRESTRICTED, U, void> : public std::pair<U, U> {
 public:
  using std::pair<U, U>::operator=;

  inline SpinPolarizedData(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& orig)
    : std::pair<U, U>((const std::pair<U, U>&)orig), alpha(this->first), beta(this->second) {
  }
  inline SpinPolarizedData(SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& orig)
    : std::pair<U, U>((std::pair<U, U>&)orig), alpha(this->first), beta(this->second) {
  }
  inline SpinPolarizedData(SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&& orig)
    : std::pair<U, U>(std::move(orig)), alpha(this->first), beta(this->second) {
  }
  inline SpinPolarizedData(std::piecewise_construct_t, U&& alpha, U&& beta)
    : std::pair<U, U>(std::forward<U>(alpha), std::forward<U>(beta)) {
  }
  inline SpinPolarizedData(std::piecewise_construct_t, const U& alpha, const U& beta) : std::pair<U, U>(alpha, beta) {
  }

  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&
  operator=(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& orig) {
    this->alpha = orig.alpha;
    this->beta = orig.beta;
    return *this;
  }
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&
  operator=(SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&& orig) {
    this->alpha = std::move(orig.alpha);
    this->beta = std::move(orig.beta);
    return *this;
  }
  /*
   * Forward constructors of a single U-type object. I looked at std::make_shared for this.
   */
  template<typename... _Args>
  inline SpinPolarizedData(_Args&&... __args)
    : std::pair<U, U>(std::piecewise_construct, std::forward_as_tuple(__args...), std::forward_as_tuple(__args...)) {
  }
  U& alpha = this->first;
  U& beta = this->second;

  inline U total() {
    return this->alpha + this->beta;
  }

  static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::UNRESTRICTED;
  /// @brief The type
  typedef SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void> type;
  /// @brief The underlying type
  typedef U& spinlesstype;
  /// @brief The underlying type
  typedef const U& constspinlesstype;
  /*
   * Arithmetic operators
   */

  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>&
  operator+=(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) {
    this->alpha += y.alpha;
    this->beta += y.beta;
    return *this;
  }
  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>
  operator+(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) const {
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> result(*this);
    result += y;
    return result;
  }

  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>&
  operator*=(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) {
    this->alpha *= y.alpha;
    this->beta *= y.beta;
    return *this;
  }
  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>
  operator*(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) const {
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> result(*this);
    result *= y;
    return result;
  }

  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>&
  operator-=(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) {
    this->alpha -= y.alpha;
    this->beta -= y.beta;
    return *this;
  }
  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>
  operator-(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) const {
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> result(*this);
    result -= y;
    return result;
  }

  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>&
  operator/=(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) {
    this->alpha /= y.alpha;
    this->beta /= y.beta;
    return *this;
  }
  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>
  operator/(const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) const {
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> result(*this);
    result /= y;
    return result;
  }
};

template<class U>
inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> makeUnrestrictedFromPieces(U&& alpha, U&& beta) {
  return SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>(std::piecewise_construct, std::forward<U>(alpha),
                                                                std::forward<U>(beta));
}
template<class U>
inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> makeUnrestrictedFromPieces(const U& alpha, const U& beta) {
  return SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>(std::piecewise_construct, alpha, beta);
}

/**************************************************************************************************
 *
 * In the following, the for_spin loop structure is defined for up to five arguments. Unfortunately
 * this became quite lengthy, because I did not find a way to generalize a few things (e.g. the
 * number of arguments). - TD
 *
 * Any class that should be usable in a for_spin macro needs the following things:
 * - it needs to be templated as: template<Options::SCF_MODES SCFModes>
 * - it needs the following member variables for the macro:
 *   - a member for restricted/unrestricted:
 *     static constexpr Options::SCF_MODES scf_mode = Options::SCF_MODES::RESTRICTED;
 *   - a member detailing the class type:
 *     typedef MatrixInBasis<UNRESTRICTED> type;
 *   - a member detailing the polarized bas type:
 *     typedef Eigen::MatrixXd spinlesstype;
 *   - a member detailing a const version of the latter:
 *     typedef const Eigen::MatrixXd constspinlesstype;
 * - it needs a member for alpha and beta version of the polarized/unrestricted code:
 *   Eigen::MatrixXd alpha;
 *   Eigen::MatrixXd beta;
 * - the restricted version is expected to inherit from 'spinlesstype'
 * While this description does not exactly explain what the code below does,it should be enough
 * to help write code using the mess below. -JU
 **************************************************************************************************/
/*
 * Single argument
 */
template<Options::SCF_MODES SCFMode, class U, bool isConstU>
class ForSpinHelperOne;
template<class U, bool isConstU>
class ForSpinHelperOne<Options::SCF_MODES::RESTRICTED, U, isConstU> {
 public:
  inline ForSpinHelperOne(typename std::conditional<isConstU, const U&, U&>::type u) : _u(u){};
  void operator<<(const std::function<void(typename std::conditional<isConstU, const U&, U&>::type&)>& f) {
    f(_u);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
};
template<class U, bool isConstU>
class ForSpinHelperOne<Options::SCF_MODES::UNRESTRICTED, U, isConstU> {
 public:
  inline ForSpinHelperOne(typename std::conditional<isConstU, const U&, U&>::type u) : _u(u){};
  void operator<<(
      const std::function<void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&)>& f) {
    f(_u.alpha);
    f(_u.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
};
#define _FOR_SPIN_1(a)                                                                                                          \
  ForSpinHelperOne<std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type, \
                   std::is_const<typename std::remove_reference<decltype(a)>::type>::value>(a)                                  \
      << [&](typename ForSpinHelperOne<std::remove_reference<decltype(a)>::type::scf_mode,                                      \
                                       typename std::remove_reference<decltype(a)>::type::type,                                 \
                                       std::is_const<typename std::remove_reference<decltype(a)>::type>::value>::typeU &        \
             a##_spin)

/*
 * Two arguments
 */
template<Options::SCF_MODES SCFMode, class U, class V, bool isConstU, bool isConstV>
class ForSpinHelperTwo;
template<class U, class V, bool isConstU, bool isConstV>
class ForSpinHelperTwo<Options::SCF_MODES::RESTRICTED, U, V, isConstU, isConstV> {
 public:
  inline ForSpinHelperTwo(typename std::conditional<isConstU, const U&, U&>::type u,
                          typename std::conditional<isConstV, const V&, V&>::type v)
    : _u(u), _v(v){};
  void operator<<(const std::function<void(typename std::conditional<isConstU, const U&, U&>::type&,
                                           typename std::conditional<isConstV, const V&, V&>::type&)>& f) {
    f(_u, _v);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
};
template<class U, class V, bool isConstU, bool isConstV>
class ForSpinHelperTwo<Options::SCF_MODES::UNRESTRICTED, U, V, isConstU, isConstV> {
 public:
  inline ForSpinHelperTwo(typename std::conditional<isConstU, const U&, U&>::type u,
                          typename std::conditional<isConstV, const V&, V&>::type v)
    : _u(u), _v(v){};
  void operator<<(
      const std::function<void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                               typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha);
    f(_u.beta, _v.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
};
#define _FOR_SPIN_2(a, b)                                                                                                       \
  ForSpinHelperTwo<std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type, \
                   typename std::remove_reference<decltype(b)>::type::type,                                                     \
                   std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                     \
                   std::is_const<typename std::remove_reference<decltype(b)>::type>::value>(a, b)                               \
      << [&](typename ForSpinHelperTwo<std::remove_reference<decltype(a)>::type::scf_mode,                                      \
                                       typename std::remove_reference<decltype(a)>::type::type,                                 \
                                       typename std::remove_reference<decltype(b)>::type::type,                                 \
                                       std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                 \
                                       std::is_const<typename std::remove_reference<decltype(b)>::type>::value>::typeU &        \
                 a##_spin,                                                                                                      \
             typename ForSpinHelperTwo<std::remove_reference<decltype(a)>::type::scf_mode,                                      \
                                       typename std::remove_reference<decltype(a)>::type::type,                                 \
                                       typename std::remove_reference<decltype(b)>::type::type,                                 \
                                       std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                 \
                                       std::is_const<typename std::remove_reference<decltype(b)>::type>::value>::typeV &        \
                 b##_spin)

/*
 * Three arguments
 */
template<Options::SCF_MODES SCFMode, class U, class V, class W, bool isConstU, bool isConstV, bool isConstW>
class ForSpinHelperThree;
template<class U, class V, class W, bool isConstU, bool isConstV, bool isConstW>
class ForSpinHelperThree<Options::SCF_MODES::RESTRICTED, U, V, W, isConstU, isConstV, isConstW> {
 public:
  inline ForSpinHelperThree(typename std::conditional<isConstU, const U&, U&>::type u,
                            typename std::conditional<isConstV, const V&, V&>::type v,
                            typename std::conditional<isConstW, const W&, W&>::type w)
    : _u(u), _v(v), _w(w){};
  void operator<<(const std::function<void(typename std::conditional<isConstU, const U&, U&>::type&,
                                           typename std::conditional<isConstV, const V&, V&>::type&,
                                           typename std::conditional<isConstW, const W&, W&>::type&)>& f) {
    f(_u, _v, _w);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
};
template<class U, class V, class W, bool isConstU, bool isConstV, bool isConstW>
class ForSpinHelperThree<Options::SCF_MODES::UNRESTRICTED, U, V, W, isConstU, isConstV, isConstW> {
 public:
  inline ForSpinHelperThree(typename std::conditional<isConstU, const U&, U&>::type u,
                            typename std::conditional<isConstV, const V&, V&>::type v,
                            typename std::conditional<isConstW, const W&, W&>::type w)
    : _u(u), _v(v), _w(w){};
  void
  operator<<(const std::function<
             void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                  typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&,
                  typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha);
    f(_u.beta, _v.beta, _w.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type typeW;
};
#define _FOR_SPIN_3(a, b, c)                                                                                                           \
  ForSpinHelperThree<std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
                     typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
                     std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
                     std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
                     std::is_const<typename std::remove_reference<decltype(c)>::type>::value>(a, b, c)                                 \
      << [&](typename ForSpinHelperThree<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,          \
                 typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,     \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                              \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                              \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value>::typeU &                                     \
                 a##_spin,                                                                                                             \
             typename ForSpinHelperThree<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,          \
                 typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,     \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                              \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                              \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value>::typeV &                                     \
                 b##_spin,                                                                                                             \
             typename ForSpinHelperThree<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,          \
                 typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,     \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                              \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                              \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value>::typeW &                                     \
                 c##_spin)

/*
 * Four arguments
 */
template<Options::SCF_MODES SCFMode, class U, class V, class W, class X, bool isConstU, bool isConstV, bool isConstW, bool isConstX>
class ForSpinHelperFour;
template<class U, class V, class W, class X, bool isConstU, bool isConstV, bool isConstW, bool isConstX>
class ForSpinHelperFour<Options::SCF_MODES::RESTRICTED, U, V, W, X, isConstU, isConstV, isConstW, isConstX> {
 public:
  inline ForSpinHelperFour(typename std::conditional<isConstU, const U&, U&>::type u,
                           typename std::conditional<isConstV, const V&, V&>::type v,
                           typename std::conditional<isConstW, const W&, W&>::type w,
                           typename std::conditional<isConstX, const X&, X&>::type x)
    : _u(u), _v(v), _w(w), _x(x){};
  void operator<<(const std::function<void(typename std::conditional<isConstU, const U&, U&>::type&,
                                           typename std::conditional<isConstV, const V&, V&>::type&,
                                           typename std::conditional<isConstW, const W&, W&>::type&,
                                           typename std::conditional<isConstX, const X&, X&>::type&)>& f) {
    f(_u, _v, _w, _x);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
};
template<class U, class V, class W, class X, bool isConstU, bool isConstV, bool isConstW, bool isConstX>
class ForSpinHelperFour<Options::SCF_MODES::UNRESTRICTED, U, V, W, X, isConstU, isConstV, isConstW, isConstX> {
 public:
  inline ForSpinHelperFour(typename std::conditional<isConstU, const U&, U&>::type u,
                           typename std::conditional<isConstV, const V&, V&>::type v,
                           typename std::conditional<isConstW, const W&, W&>::type w,
                           typename std::conditional<isConstX, const X&, X&>::type x)
    : _u(u), _v(v), _w(w), _x(x){};
  void
  operator<<(const std::function<
             void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                  typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&,
                  typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type&,
                  typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type typeX;
};
#define _FOR_SPIN_4(a, b, c, d)                                                                                                       \
  ForSpinHelperFour<std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
                    typename std::remove_reference<decltype(b)>::type::type,                                                          \
                    typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type, \
                    std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
                    std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
                    std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
                    std::is_const<typename std::remove_reference<decltype(d)>::type>::value>(a, b, c, d)                              \
      << [&](typename ForSpinHelperFour<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type,                                                             \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value>::typeU &                                    \
                 a##_spin,                                                                                                            \
             typename ForSpinHelperFour<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type,                                                             \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value>::typeV &                                    \
                 b##_spin,                                                                                                            \
             typename ForSpinHelperFour<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type,                                                             \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value>::typeW &                                    \
                 c##_spin,                                                                                                            \
             typename ForSpinHelperFour<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type,                                                             \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value>::typeX &                                    \
                 d##_spin)

/*
 * Five arguments
 */
template<Options::SCF_MODES SCFMode, class U, class V, class W, class X, class Y, bool isConstU, bool isConstV,
         bool isConstW, bool isConstX, bool isConstY>
class ForSpinHelperFive;
template<class U, class V, class W, class X, class Y, bool isConstU, bool isConstV, bool isConstW, bool isConstX, bool isConstY>
class ForSpinHelperFive<Options::SCF_MODES::RESTRICTED, U, V, W, X, Y, isConstU, isConstV, isConstW, isConstX, isConstY> {
 public:
  inline ForSpinHelperFive(typename std::conditional<isConstU, const U&, U&>::type u,
                           typename std::conditional<isConstV, const V&, V&>::type v,
                           typename std::conditional<isConstW, const W&, W&>::type w,
                           typename std::conditional<isConstX, const X&, X&>::type x,
                           typename std::conditional<isConstY, const Y&, Y&>::type y)
    : _u(u), _v(v), _w(w), _x(x), _y(y){};
  void operator<<(const std::function<void(typename std::conditional<isConstU, const U&, U&>::type&,
                                           typename std::conditional<isConstV, const V&, V&>::type&,
                                           typename std::conditional<isConstW, const W&, W&>::type&,
                                           typename std::conditional<isConstX, const X&, X&>::type&,
                                           typename std::conditional<isConstY, const Y&, Y&>::type&)>& f) {
    f(_u, _v, _w, _x, _y);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
  typedef typename std::conditional<isConstY, const Y&, Y&>::type typeY;
};
template<class U, class V, class W, class X, class Y, bool isConstU, bool isConstV, bool isConstW, bool isConstX, bool isConstY>
class ForSpinHelperFive<Options::SCF_MODES::UNRESTRICTED, U, V, W, X, Y, isConstU, isConstV, isConstW, isConstX, isConstY> {
 public:
  inline ForSpinHelperFive(typename std::conditional<isConstU, const U&, U&>::type u,
                           typename std::conditional<isConstV, const V&, V&>::type v,
                           typename std::conditional<isConstW, const W&, W&>::type w,
                           typename std::conditional<isConstX, const X&, X&>::type x,
                           typename std::conditional<isConstY, const Y&, Y&>::type y)
    : _u(u), _v(v), _w(w), _x(x), _y(y){};
  void
  operator<<(const std::function<
             void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                  typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&,
                  typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type&,
                  typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type&,
                  typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha, _y.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta, _y.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type typeX;
  typedef typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type typeY;
};
#define _FOR_SPIN_5(a, b, c, d, e)                                                                                                    \
  ForSpinHelperFive<std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
                    typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
                    typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
                    std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
                    std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
                    std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
                    std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
                    std::is_const<typename std::remove_reference<decltype(e)>::type>::value>(a, b, c, d, e)                           \
      << [&](typename ForSpinHelperFive<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,    \
                 typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeU &                                    \
                 a##_spin,                                                                                                            \
             typename ForSpinHelperFive<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,    \
                 typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeV &                                    \
                 b##_spin,                                                                                                            \
             typename ForSpinHelperFive<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,    \
                 typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeW &                                    \
                 c##_spin,                                                                                                            \
             typename ForSpinHelperFive<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,    \
                 typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeX &                                    \
                 d##_spin,                                                                                                            \
             typename ForSpinHelperFive<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,         \
                 typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,    \
                 typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,    \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                             \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeY &                                    \
                 e##_spin)

/*
 * Six arguments
 */

template<Options::SCF_MODES SCFMode, class U, class V, class W, class X, class Y, class Z, bool isConstU, bool isConstV,
         bool isConstW, bool isConstX, bool isConstY, bool isConstZ>
class ForSpinHelperSix;
template<class U, class V, class W, class X, class Y, class Z, bool isConstU, bool isConstV, bool isConstW, bool isConstX, bool isConstY, bool isConstZ>
class ForSpinHelperSix<Options::SCF_MODES::RESTRICTED, U, V, W, X, Y, Z, isConstU, isConstV, isConstW, isConstX, isConstY, isConstZ> {
 public:
  inline ForSpinHelperSix(typename std::conditional<isConstU, const U&, U&>::type u,
                          typename std::conditional<isConstV, const V&, V&>::type v,
                          typename std::conditional<isConstW, const W&, W&>::type w,
                          typename std::conditional<isConstX, const X&, X&>::type x,
                          typename std::conditional<isConstY, const Y&, Y&>::type y,
                          typename std::conditional<isConstZ, const Z&, Z&>::type z)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z){};
  void operator<<(
      const std::function<void(
          typename std::conditional<isConstU, const U&, U&>::type&, typename std::conditional<isConstV, const V&, V&>::type&,
          typename std::conditional<isConstW, const W&, W&>::type&, typename std::conditional<isConstX, const X&, X&>::type&,
          typename std::conditional<isConstY, const Y&, Y&>::type&, typename std::conditional<isConstZ, const Z&, Z&>::type&)>& f) {
    f(_u, _v, _w, _x, _y, _z);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
  typedef typename std::conditional<isConstY, const Y&, Y&>::type typeY;
  typedef typename std::conditional<isConstZ, const Z&, Z&>::type typeZ;
};

template<class U, class V, class W, class X, class Y, class Z, bool isConstU, bool isConstV, bool isConstW, bool isConstX, bool isConstY, bool isConstZ>
class ForSpinHelperSix<Options::SCF_MODES::UNRESTRICTED, U, V, W, X, Y, Z, isConstU, isConstV, isConstW, isConstX, isConstY, isConstZ> {
 public:
  inline ForSpinHelperSix(typename std::conditional<isConstU, const U&, U&>::type u,
                          typename std::conditional<isConstV, const V&, V&>::type v,
                          typename std::conditional<isConstW, const W&, W&>::type w,
                          typename std::conditional<isConstX, const X&, X&>::type x,
                          typename std::conditional<isConstY, const Y&, Y&>::type y,
                          typename std::conditional<isConstZ, const Z&, Z&>::type z)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z){};
  void
  operator<<(const std::function<
             void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                  typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&,
                  typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type&,
                  typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type&,
                  typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type&,
                  typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha, _y.alpha, _z.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta, _y.beta, _z.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type typeX;
  typedef typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type typeY;
  typedef typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type typeZ;
};

#define _FOR_SPIN_6(a, b, c, d, e, f)                                                                                                \
  ForSpinHelperSix<std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
                   typename std::remove_reference<decltype(b)>::type::type,                                                          \
                   typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type, \
                   typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type, \
                   std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
                   std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
                   std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
                   std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
                   std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
                   std::is_const<typename std::remove_reference<decltype(f)>::type>::value>(a, b, c, d, e, f)                        \
      << [&](typename ForSpinHelperSix<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,        \
                 typename std::remove_reference<decltype(b)>::type::type,                                                            \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,   \
                 typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,   \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(f)>::type>::value>::typeU &                                   \
                 a##_spin,                                                                                                           \
             typename ForSpinHelperSix<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,        \
                 typename std::remove_reference<decltype(b)>::type::type,                                                            \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,   \
                 typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,   \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(f)>::type>::value>::typeV &                                   \
                 b##_spin,                                                                                                           \
             typename ForSpinHelperSix<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,        \
                 typename std::remove_reference<decltype(b)>::type::type,                                                            \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,   \
                 typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,   \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(f)>::type>::value>::typeW &                                   \
                 c##_spin,                                                                                                           \
             typename ForSpinHelperSix<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,        \
                 typename std::remove_reference<decltype(b)>::type::type,                                                            \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,   \
                 typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,   \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(f)>::type>::value>::typeX &                                   \
                 d##_spin,                                                                                                           \
             typename ForSpinHelperSix<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,        \
                 typename std::remove_reference<decltype(b)>::type::type,                                                            \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,   \
                 typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,   \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(f)>::type>::value>::typeY &                                   \
                 e##_spin,                                                                                                           \
             typename ForSpinHelperSix<                                                                                              \
                 std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,        \
                 typename std::remove_reference<decltype(b)>::type::type,                                                            \
                 typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,   \
                 typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,   \
                 std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                            \
                 std::is_const<typename std::remove_reference<decltype(f)>::type>::value>::typeZ &                                   \
                 f##_spin)

/*
 * Seven arguments
 */

template<Options::SCF_MODES SCFMode, class U, class V, class W, class X, class Y, class Z, class Q, bool isConstU,
         bool isConstV, bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ>
class ForSpinHelperSeven;
template<class U, class V, class W, class X, class Y, class Z, class Q, bool isConstU, bool isConstV, bool isConstW,
         bool isConstX, bool isConstY, bool isConstZ, bool isConstQ>
class ForSpinHelperSeven<Options::SCF_MODES::RESTRICTED, U, V, W, X, Y, Z, Q, isConstU, isConstV, isConstW, isConstX, isConstY, isConstZ, isConstQ> {
 public:
  inline ForSpinHelperSeven(typename std::conditional<isConstU, const U&, U&>::type u,
                            typename std::conditional<isConstV, const V&, V&>::type v,
                            typename std::conditional<isConstW, const W&, W&>::type w,
                            typename std::conditional<isConstX, const X&, X&>::type x,
                            typename std::conditional<isConstY, const Y&, Y&>::type y,
                            typename std::conditional<isConstZ, const Z&, Z&>::type z,
                            typename std::conditional<isConstQ, const Q&, Q&>::type q)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z), _q(q){};
  void operator<<(
      const std::function<void(
          typename std::conditional<isConstU, const U&, U&>::type&, typename std::conditional<isConstV, const V&, V&>::type&,
          typename std::conditional<isConstW, const W&, W&>::type&, typename std::conditional<isConstX, const X&, X&>::type&,
          typename std::conditional<isConstY, const Y&, Y&>::type&, typename std::conditional<isConstZ, const Z&, Z&>::type&,
          typename std::conditional<isConstQ, const Q&, Q&>::type&)>& f) {
    f(_u, _v, _w, _x, _y, _z, _q);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typename std::conditional<isConstQ, const Q&, Q&>::type _q;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
  typedef typename std::conditional<isConstY, const Y&, Y&>::type typeY;
  typedef typename std::conditional<isConstZ, const Z&, Z&>::type typeZ;
  typedef typename std::conditional<isConstQ, const Q&, Q&>::type typeQ;
};

template<class U, class V, class W, class X, class Y, class Z, class Q, bool isConstU, bool isConstV, bool isConstW,
         bool isConstX, bool isConstY, bool isConstZ, bool isConstQ>
class ForSpinHelperSeven<Options::SCF_MODES::UNRESTRICTED, U, V, W, X, Y, Z, Q, isConstU, isConstV, isConstW, isConstX,
                         isConstY, isConstZ, isConstQ> {
 public:
  inline ForSpinHelperSeven(typename std::conditional<isConstU, const U&, U&>::type u,
                            typename std::conditional<isConstV, const V&, V&>::type v,
                            typename std::conditional<isConstW, const W&, W&>::type w,
                            typename std::conditional<isConstX, const X&, X&>::type x,
                            typename std::conditional<isConstY, const Y&, Y&>::type y,
                            typename std::conditional<isConstZ, const Z&, Z&>::type z,
                            typename std::conditional<isConstQ, const Q&, Q&>::type q)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z), _q(q){};
  void
  operator<<(const std::function<
             void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                  typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&,
                  typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type&,
                  typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type&,
                  typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type&,
                  typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type&,
                  typename std::conditional<isConstQ, typename Q::constspinlesstype, typename Q::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha, _y.alpha, _z.alpha, _q.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta, _y.beta, _z.beta, _q.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typename std::conditional<isConstQ, const Q&, Q&>::type _q;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type typeX;
  typedef typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type typeY;
  typedef typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type typeZ;
  typedef typename std::conditional<isConstQ, typename Q::constspinlesstype, typename Q::spinlesstype>::type typeQ;
};

#define _FOR_SPIN_7(a, b, c, d, e, f, g)                                                                                        \
  ForSpinHelperSeven<                                                                                                           \
      std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
      typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,         \
      typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,         \
      typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type,         \
      std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(g)>::type>::value>(a, b, c, d, e, f, g)                             \
      <<                                                                                                                        \
      [&](typename ForSpinHelperSeven<                                                                                          \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value>::typeU &                                 \
              a##_spin,                                                                                                         \
          typename ForSpinHelperSeven<                                                                                          \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value>::typeV &                                 \
              b##_spin,                                                                                                         \
          typename ForSpinHelperSeven<                                                                                          \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value>::typeW &                                 \
              c##_spin,                                                                                                         \
          typename ForSpinHelperSeven<                                                                                          \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value>::typeX &                                 \
              d##_spin,                                                                                                         \
          typename ForSpinHelperSeven<                                                                                          \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value>::typeY &                                 \
              e##_spin,                                                                                                         \
          typename ForSpinHelperSeven<                                                                                          \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value>::typeZ &                                 \
              f##_spin,                                                                                                         \
          typename ForSpinHelperSeven<                                                                                          \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value>::typeQ &                                 \
              g##_spin)

/*
 * Eight arguments
 */

template<Options::SCF_MODES SCFMode, class U, class V, class W, class X, class Y, class Z, class Q, class R, bool isConstU,
         bool isConstV, bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ, bool isConstR>
class ForSpinHelperEight;
template<class U, class V, class W, class X, class Y, class Z, class Q, class R, bool isConstU, bool isConstV,
         bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ, bool isConstR>
class ForSpinHelperEight<Options::SCF_MODES::RESTRICTED, U, V, W, X, Y, Z, Q, R, isConstU, isConstV, isConstW, isConstX,
                         isConstY, isConstZ, isConstQ, isConstR> {
 public:
  inline ForSpinHelperEight(typename std::conditional<isConstU, const U&, U&>::type u,
                            typename std::conditional<isConstV, const V&, V&>::type v,
                            typename std::conditional<isConstW, const W&, W&>::type w,
                            typename std::conditional<isConstX, const X&, X&>::type x,
                            typename std::conditional<isConstY, const Y&, Y&>::type y,
                            typename std::conditional<isConstZ, const Z&, Z&>::type z,
                            typename std::conditional<isConstQ, const Q&, Q&>::type q,
                            typename std::conditional<isConstR, const R&, R&>::type r)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z), _q(q), _r(r){};
  void operator<<(
      const std::function<void(
          typename std::conditional<isConstU, const U&, U&>::type&, typename std::conditional<isConstV, const V&, V&>::type&,
          typename std::conditional<isConstW, const W&, W&>::type&, typename std::conditional<isConstX, const X&, X&>::type&,
          typename std::conditional<isConstY, const Y&, Y&>::type&, typename std::conditional<isConstZ, const Z&, Z&>::type&,
          typename std::conditional<isConstQ, const Q&, Q&>::type&, typename std::conditional<isConstR, const R&, R&>::type&)>& f) {
    f(_u, _v, _w, _x, _y, _z, _q, _r);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typename std::conditional<isConstQ, const Q&, Q&>::type _q;
  typename std::conditional<isConstR, const R&, R&>::type _r;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
  typedef typename std::conditional<isConstY, const Y&, Y&>::type typeY;
  typedef typename std::conditional<isConstZ, const Z&, Z&>::type typeZ;
  typedef typename std::conditional<isConstQ, const Q&, Q&>::type typeQ;
  typedef typename std::conditional<isConstR, const R&, R&>::type typeR;
};

template<class U, class V, class W, class X, class Y, class Z, class Q, class R, bool isConstU, bool isConstV,
         bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ, bool isConstR>
class ForSpinHelperEight<Options::SCF_MODES::UNRESTRICTED, U, V, W, X, Y, Z, Q, R, isConstU, isConstV, isConstW,
                         isConstX, isConstY, isConstZ, isConstQ, isConstR> {
 public:
  inline ForSpinHelperEight(typename std::conditional<isConstU, const U&, U&>::type u,
                            typename std::conditional<isConstV, const V&, V&>::type v,
                            typename std::conditional<isConstW, const W&, W&>::type w,
                            typename std::conditional<isConstX, const X&, X&>::type x,
                            typename std::conditional<isConstY, const Y&, Y&>::type y,
                            typename std::conditional<isConstZ, const Z&, Z&>::type z,
                            typename std::conditional<isConstQ, const Q&, Q&>::type q,
                            typename std::conditional<isConstR, const R&, R&>::type r)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z), _q(q), _r(r){};
  void
  operator<<(const std::function<
             void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                  typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&,
                  typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type&,
                  typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type&,
                  typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type&,
                  typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type&,
                  typename std::conditional<isConstQ, typename Q::constspinlesstype, typename Q::spinlesstype>::type&,
                  typename std::conditional<isConstR, typename R::constspinlesstype, typename R::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha, _y.alpha, _z.alpha, _q.alpha, _r.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta, _y.beta, _z.beta, _q.beta, _r.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typename std::conditional<isConstQ, const Q&, Q&>::type _q;
  typename std::conditional<isConstR, const R&, R&>::type _r;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type typeX;
  typedef typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type typeY;
  typedef typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type typeZ;
  typedef typename std::conditional<isConstQ, typename Q::constspinlesstype, typename Q::spinlesstype>::type typeQ;
  typedef typename std::conditional<isConstR, typename R::constspinlesstype, typename R::spinlesstype>::type typeR;
};

#define _FOR_SPIN_8(a, b, c, d, e, f, g, h)                                                                                             \
  ForSpinHelperEight<                                                                                                                   \
      std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,                      \
      typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,                 \
      typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,                 \
      typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type,                 \
      typename std::remove_reference<decltype(h)>::type::type, std::is_const<typename std::remove_reference<decltype(a)>::type>::value, \
      std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(h)>::type>::value>(a, b, c, d, e, f, g, h)                                  \
      <<                                                                                                                                \
      [&](typename ForSpinHelperEight<                                                                                                  \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value>::typeU &                                         \
              a##_spin,                                                                                                                 \
          typename ForSpinHelperEight<                                                                                                  \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value>::typeV &                                         \
              b##_spin,                                                                                                                 \
          typename ForSpinHelperEight<                                                                                                  \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value>::typeW &                                         \
              c##_spin,                                                                                                                 \
          typename ForSpinHelperEight<                                                                                                  \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value>::typeX &                                         \
              d##_spin,                                                                                                                 \
          typename ForSpinHelperEight<                                                                                                  \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value>::typeY &                                         \
              e##_spin,                                                                                                                 \
          typename ForSpinHelperEight<                                                                                                  \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value>::typeZ &                                         \
              f##_spin,                                                                                                                 \
          typename ForSpinHelperEight<                                                                                                  \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value>::typeQ &                                         \
              g##_spin,                                                                                                                 \
          typename ForSpinHelperEight<                                                                                                  \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value>::typeR &                                         \
              h##_spin)

/*
 * Nine arguments
 */

template<Options::SCF_MODES SCFMode, class U, class V, class W, class X, class Y, class Z, class Q, class R, class S, bool isConstU,
         bool isConstV, bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ, bool isConstR, bool isConstS>
class ForSpinHelperNine;
template<class U, class V, class W, class X, class Y, class Z, class Q, class R, class S, bool isConstU, bool isConstV,
         bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ, bool isConstR, bool isConstS>
class ForSpinHelperNine<Options::SCF_MODES::RESTRICTED, U, V, W, X, Y, Z, Q, R, S, isConstU, isConstV, isConstW,
                        isConstX, isConstY, isConstZ, isConstQ, isConstR, isConstS> {
 public:
  inline ForSpinHelperNine(typename std::conditional<isConstU, const U&, U&>::type u,
                           typename std::conditional<isConstV, const V&, V&>::type v,
                           typename std::conditional<isConstW, const W&, W&>::type w,
                           typename std::conditional<isConstX, const X&, X&>::type x,
                           typename std::conditional<isConstY, const Y&, Y&>::type y,
                           typename std::conditional<isConstZ, const Z&, Z&>::type z,
                           typename std::conditional<isConstQ, const Q&, Q&>::type q,
                           typename std::conditional<isConstR, const R&, R&>::type r,
                           typename std::conditional<isConstS, const S&, S&>::type s)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z), _q(q), _r(r), _s(s){};
  void operator<<(
      const std::function<void(
          typename std::conditional<isConstU, const U&, U&>::type&, typename std::conditional<isConstV, const V&, V&>::type&,
          typename std::conditional<isConstW, const W&, W&>::type&, typename std::conditional<isConstX, const X&, X&>::type&,
          typename std::conditional<isConstY, const Y&, Y&>::type&, typename std::conditional<isConstZ, const Z&, Z&>::type&,
          typename std::conditional<isConstQ, const Q&, Q&>::type&, typename std::conditional<isConstR, const R&, R&>::type&,
          typename std::conditional<isConstS, const S&, S&>::type&)>& f) {
    f(_u, _v, _w, _x, _y, _z, _q, _r, _s);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typename std::conditional<isConstQ, const Q&, Q&>::type _q;
  typename std::conditional<isConstR, const R&, R&>::type _r;
  typename std::conditional<isConstS, const S&, S&>::type _s;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
  typedef typename std::conditional<isConstY, const Y&, Y&>::type typeY;
  typedef typename std::conditional<isConstZ, const Z&, Z&>::type typeZ;
  typedef typename std::conditional<isConstQ, const Q&, Q&>::type typeQ;
  typedef typename std::conditional<isConstR, const R&, R&>::type typeR;
  typedef typename std::conditional<isConstS, const S&, S&>::type typeS;
};

template<class U, class V, class W, class X, class Y, class Z, class Q, class R, class S, bool isConstU, bool isConstV,
         bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ, bool isConstR, bool isConstS>
class ForSpinHelperNine<Options::SCF_MODES::UNRESTRICTED, U, V, W, X, Y, Z, Q, R, S, isConstU, isConstV, isConstW,
                        isConstX, isConstY, isConstZ, isConstQ, isConstR, isConstS> {
 public:
  inline ForSpinHelperNine(typename std::conditional<isConstU, const U&, U&>::type u,
                           typename std::conditional<isConstV, const V&, V&>::type v,
                           typename std::conditional<isConstW, const W&, W&>::type w,
                           typename std::conditional<isConstX, const X&, X&>::type x,
                           typename std::conditional<isConstY, const Y&, Y&>::type y,
                           typename std::conditional<isConstZ, const Z&, Z&>::type z,
                           typename std::conditional<isConstQ, const Q&, Q&>::type q,
                           typename std::conditional<isConstR, const R&, R&>::type r,
                           typename std::conditional<isConstS, const S&, S&>::type s)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z), _q(q), _r(r), _s(s){};
  void
  operator<<(const std::function<
             void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                  typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&,
                  typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type&,
                  typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type&,
                  typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type&,
                  typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type&,
                  typename std::conditional<isConstQ, typename Q::constspinlesstype, typename Q::spinlesstype>::type&,
                  typename std::conditional<isConstR, typename R::constspinlesstype, typename R::spinlesstype>::type&,
                  typename std::conditional<isConstS, typename S::constspinlesstype, typename S::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha, _y.alpha, _z.alpha, _q.alpha, _r.alpha, _s.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta, _y.beta, _z.beta, _q.beta, _r.beta, _s.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typename std::conditional<isConstQ, const Q&, Q&>::type _q;
  typename std::conditional<isConstR, const R&, R&>::type _r;
  typename std::conditional<isConstS, const S&, S&>::type _s;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type typeX;
  typedef typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type typeY;
  typedef typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type typeZ;
  typedef typename std::conditional<isConstQ, typename Q::constspinlesstype, typename Q::spinlesstype>::type typeQ;
  typedef typename std::conditional<isConstR, typename R::constspinlesstype, typename R::spinlesstype>::type typeR;
  typedef typename std::conditional<isConstS, typename S::constspinlesstype, typename S::spinlesstype>::type typeS;
};

#define _FOR_SPIN_9(a, b, c, d, e, f, g, h, k)                                                                                  \
  ForSpinHelperNine<                                                                                                            \
      std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
      typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,         \
      typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,         \
      typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type,         \
      typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type,         \
      std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
      std::is_const<typename std::remove_reference<decltype(k)>::type>::value>(a, b, c, d, e, f, g, h, k)                       \
      <<                                                                                                                        \
      [&](typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeU &                                 \
              a##_spin,                                                                                                         \
          typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeV &                                 \
              b##_spin,                                                                                                         \
          typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeW &                                 \
              c##_spin,                                                                                                         \
          typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeX &                                 \
              d##_spin,                                                                                                         \
          typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeY &                                 \
              e##_spin,                                                                                                         \
          typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeZ &                                 \
              f##_spin,                                                                                                         \
          typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeQ &                                 \
              g##_spin,                                                                                                         \
          typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeR &                                 \
              h##_spin,                                                                                                         \
          typename ForSpinHelperNine<                                                                                           \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,      \
              typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type, \
              typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type, \
              typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type, \
              typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type, \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                          \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value>::typeS &                                 \
              k##_spin)

/*
 * Ten arguments
 */

template<Options::SCF_MODES SCFMode, class U, class V, class W, class X, class Y, class Z, class Q, class R, class S,
         class P, bool isConstU, bool isConstV, bool isConstW, bool isConstX, bool isConstY, bool isConstZ,
         bool isConstQ, bool isConstR, bool isConstS, bool isConstP>
class ForSpinHelperTen;
template<class U, class V, class W, class X, class Y, class Z, class Q, class R, class S, class P, bool isConstU, bool isConstV,
         bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ, bool isConstR, bool isConstS, bool isConstP>
class ForSpinHelperTen<Options::SCF_MODES::RESTRICTED, U, V, W, X, Y, Z, Q, R, S, P, isConstU, isConstV, isConstW,
                       isConstX, isConstY, isConstZ, isConstQ, isConstR, isConstS, isConstP> {
 public:
  inline ForSpinHelperTen(
      typename std::conditional<isConstU, const U&, U&>::type u, typename std::conditional<isConstV, const V&, V&>::type v,
      typename std::conditional<isConstW, const W&, W&>::type w, typename std::conditional<isConstX, const X&, X&>::type x,
      typename std::conditional<isConstY, const Y&, Y&>::type y, typename std::conditional<isConstZ, const Z&, Z&>::type z,
      typename std::conditional<isConstQ, const Q&, Q&>::type q, typename std::conditional<isConstR, const R&, R&>::type r,
      typename std::conditional<isConstS, const S&, S&>::type s, typename std::conditional<isConstP, const P&, P&>::type p)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z), _q(q), _r(r), _s(s), _p(p){};
  void operator<<(
      const std::function<void(
          typename std::conditional<isConstU, const U&, U&>::type&, typename std::conditional<isConstV, const V&, V&>::type&,
          typename std::conditional<isConstW, const W&, W&>::type&, typename std::conditional<isConstX, const X&, X&>::type&,
          typename std::conditional<isConstY, const Y&, Y&>::type&, typename std::conditional<isConstZ, const Z&, Z&>::type&,
          typename std::conditional<isConstQ, const Q&, Q&>::type&, typename std::conditional<isConstR, const R&, R&>::type&,
          typename std::conditional<isConstS, const S&, S&>::type&, typename std::conditional<isConstP, const P&, P&>::type&)>& f) {
    f(_u, _v, _w, _x, _y, _z, _q, _r, _s, _p);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typename std::conditional<isConstQ, const Q&, Q&>::type _q;
  typename std::conditional<isConstR, const R&, R&>::type _r;
  typename std::conditional<isConstS, const S&, S&>::type _s;
  typename std::conditional<isConstP, const P&, P&>::type _p;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
  typedef typename std::conditional<isConstY, const Y&, Y&>::type typeY;
  typedef typename std::conditional<isConstZ, const Z&, Z&>::type typeZ;
  typedef typename std::conditional<isConstQ, const Q&, Q&>::type typeQ;
  typedef typename std::conditional<isConstR, const R&, R&>::type typeR;
  typedef typename std::conditional<isConstS, const S&, S&>::type typeS;
  typedef typename std::conditional<isConstP, const P&, P&>::type typeP;
};

template<class U, class V, class W, class X, class Y, class Z, class Q, class R, class S, class P, bool isConstU, bool isConstV,
         bool isConstW, bool isConstX, bool isConstY, bool isConstZ, bool isConstQ, bool isConstR, bool isConstS, bool isConstP>
class ForSpinHelperTen<Options::SCF_MODES::UNRESTRICTED, U, V, W, X, Y, Z, Q, R, S, P, isConstU, isConstV, isConstW,
                       isConstX, isConstY, isConstZ, isConstQ, isConstR, isConstS, isConstP> {
 public:
  inline ForSpinHelperTen(
      typename std::conditional<isConstU, const U&, U&>::type u, typename std::conditional<isConstV, const V&, V&>::type v,
      typename std::conditional<isConstW, const W&, W&>::type w, typename std::conditional<isConstX, const X&, X&>::type x,
      typename std::conditional<isConstY, const Y&, Y&>::type y, typename std::conditional<isConstZ, const Z&, Z&>::type z,
      typename std::conditional<isConstQ, const Q&, Q&>::type q, typename std::conditional<isConstR, const R&, R&>::type r,
      typename std::conditional<isConstS, const S&, S&>::type s, typename std::conditional<isConstP, const P&, P&>::type p)
    : _u(u), _v(v), _w(w), _x(x), _y(y), _z(z), _q(q), _r(r), _s(s), _p(p){};
  void
  operator<<(const std::function<
             void(typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type&,
                  typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type&,
                  typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type&,
                  typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type&,
                  typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type&,
                  typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type&,
                  typename std::conditional<isConstQ, typename Q::constspinlesstype, typename Q::spinlesstype>::type&,
                  typename std::conditional<isConstR, typename R::constspinlesstype, typename R::spinlesstype>::type&,
                  typename std::conditional<isConstS, typename S::constspinlesstype, typename S::spinlesstype>::type&,
                  typename std::conditional<isConstP, typename P::constspinlesstype, typename P::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha, _y.alpha, _z.alpha, _q.alpha, _r.alpha, _s.alpha, _p.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta, _y.beta, _z.beta, _q.beta, _r.beta, _s.beta, _p.beta);
  }
  typename std::conditional<isConstU, const U&, U&>::type _u;
  typename std::conditional<isConstV, const V&, V&>::type _v;
  typename std::conditional<isConstW, const W&, W&>::type _w;
  typename std::conditional<isConstX, const X&, X&>::type _x;
  typename std::conditional<isConstY, const Y&, Y&>::type _y;
  typename std::conditional<isConstZ, const Z&, Z&>::type _z;
  typename std::conditional<isConstQ, const Q&, Q&>::type _q;
  typename std::conditional<isConstR, const R&, R&>::type _r;
  typename std::conditional<isConstS, const S&, S&>::type _s;
  typename std::conditional<isConstP, const P&, P&>::type _p;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype, typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype, typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype, typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype, typename X::spinlesstype>::type typeX;
  typedef typename std::conditional<isConstY, typename Y::constspinlesstype, typename Y::spinlesstype>::type typeY;
  typedef typename std::conditional<isConstZ, typename Z::constspinlesstype, typename Z::spinlesstype>::type typeZ;
  typedef typename std::conditional<isConstQ, typename Q::constspinlesstype, typename Q::spinlesstype>::type typeQ;
  typedef typename std::conditional<isConstR, typename R::constspinlesstype, typename R::spinlesstype>::type typeR;
  typedef typename std::conditional<isConstS, typename S::constspinlesstype, typename S::spinlesstype>::type typeS;
  typedef typename std::conditional<isConstP, typename P::constspinlesstype, typename P::spinlesstype>::type typeP;
};

#define _FOR_SPIN_10(a, b, c, d, e, f, g, h, k, l)                                                                                      \
  ForSpinHelperTen<                                                                                                                     \
      std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,                      \
      typename std::remove_reference<decltype(b)>::type::type, typename std::remove_reference<decltype(c)>::type::type,                 \
      typename std::remove_reference<decltype(d)>::type::type, typename std::remove_reference<decltype(e)>::type::type,                 \
      typename std::remove_reference<decltype(f)>::type::type, typename std::remove_reference<decltype(g)>::type::type,                 \
      typename std::remove_reference<decltype(h)>::type::type, typename std::remove_reference<decltype(k)>::type::type,                 \
      typename std::remove_reference<decltype(l)>::type::type, std::is_const<typename std::remove_reference<decltype(a)>::type>::value, \
      std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                          \
      std::is_const<typename std::remove_reference<decltype(l)>::type>::value>(a, b, c, d, e, f, g, h, k, l)                            \
      <<                                                                                                                                \
      [&](typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeU &                                         \
              a##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeV &                                         \
              b##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeW &                                         \
              c##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeX &                                         \
              d##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeY &                                         \
              e##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeZ &                                         \
              f##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeQ &                                         \
              g##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeR &                                         \
              h##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeS &                                         \
              k##_spin,                                                                                                                 \
          typename ForSpinHelperTen<                                                                                                    \
              std::remove_reference<decltype(a)>::type::scf_mode, typename std::remove_reference<decltype(a)>::type::type,              \
              typename std::remove_reference<decltype(b)>::type::type,                                                                  \
              typename std::remove_reference<decltype(c)>::type::type, typename std::remove_reference<decltype(d)>::type::type,         \
              typename std::remove_reference<decltype(e)>::type::type, typename std::remove_reference<decltype(f)>::type::type,         \
              typename std::remove_reference<decltype(g)>::type::type, typename std::remove_reference<decltype(h)>::type::type,         \
              typename std::remove_reference<decltype(k)>::type::type, typename std::remove_reference<decltype(l)>::type::type,         \
              std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(e)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(f)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(g)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(h)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(k)>::type>::value,                                                  \
              std::is_const<typename std::remove_reference<decltype(l)>::type>::value>::typeP &                                         \
              l##_spin)

#define _GET_FOR_SPIN_MACRO(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, NAME, ...) NAME

/**
 * TODO Check whether this is the correct place for a comment -> check whether/how doxygen can pick it up!
 *
 * @fn for_spin(Args... args)
 * @brief Attempt to resemble a for-loop syntax to loop over spins.
 *
 * @param args all objects used inside which may exist in a restricted or unrestricted version
 *             (atm only up to four arguments supported, but can be trivially extended).
 *
 * In our project the restricted and unrestricted calculations are seperated via a template argument.
 * That means, that in the end the compiler generates all code which differs for these cases twice.
 * This indroduces a lot of safety, e.g. w.r.t. the used data types. As a drawback, one cannot simply
 * differ between restricted and unrestricted via a simple if (at least not up to a point where it
 * is useful). However, in most places an unrestricted calculation is basically identical to a restricted
 * one, except that things are done once for the alpha data and once for the beta data. This concept
 * is abstracted into this 'function'.\n \n
 *
 * \b Usage:\n
 * At any point of the code where in the unrestricted case everything should be done once for alpha
 * spin and once for beta spin, which is done in a restricted calculation on the plain object, this
 * function can be used similar to a for loop. Here is an example of the possible usage:\n
 * \code{.cpp}
 * int five=5;
 * for_spin(myData) {
 *   myData_spin.member;
 *   doFoo(myData_spin);
 *   myData_spin *= five;
 * };
 * \endcode
 * Note that the variable which may be restricted or unrestricted is accessed with a suffix '_spin'
 * inside the structure, and note that it has to be ended with a semicolon!
 */
#define for_spin(...)                                                                                             \
  _GET_FOR_SPIN_MACRO(__VA_ARGS__, _FOR_SPIN_10, _FOR_SPIN_9, _FOR_SPIN_8, _FOR_SPIN_7, _FOR_SPIN_6, _FOR_SPIN_5, \
                      _FOR_SPIN_4, _FOR_SPIN_3, _FOR_SPIN_2, _FOR_SPIN_1, DUMMY)                                  \
  (__VA_ARGS__)
// Note: the final argument DUMMY is only there to avoid a compiler warning.
} /* namespace Serenity */
#endif /* SpinPolarizedData_H */
