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
#ifndef SPINPOLARIZEDDATA_H_
#define	SPINPOLARIZEDDATA_H_
/* Include Serenity Internal Headers */
#include "settings/Options.h"
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
template<Options::SCF_MODES T, class U, typename E=void>class SpinPolarizedData;

/* =================================================
 *   Specialization for RESTRICTED, non-primitives
 * =================================================*/
template<class U>
class SpinPolarizedData<RESTRICTED,U,
    typename std::enable_if<std::is_class<U>::value>::type> : public U {
public:
  using U::U;
  
  /// @brief Default destructor.
  virtual ~SpinPolarizedData() = default;
  
  /// @brief Default constructor.
  inline SpinPolarizedData() : U() {}
  /**
   * @brief Copy constructor using base object.
   * @param init Initial value.
   */
  inline SpinPolarizedData(const U& init) : U(init) {}
  /**
   * @brief Move constructor using base object.
   * @param init Initial value.
   */
  inline SpinPolarizedData(U&& init) : U(std::move(init)) {}
  /**
   * @brief Copy constructor using SpinPolarizedData.
   * @param orig Initial/original value.
   */
  inline SpinPolarizedData(const SpinPolarizedData<RESTRICTED, U, void>&  orig) : U(orig){}
  /**
   * @brief Copy constructor using SpinPolarizedData.
   * @param orig Initial/original value.
   */
  inline SpinPolarizedData(      SpinPolarizedData<RESTRICTED, U, void>&  orig) : U(orig){}
  /**
   * @brief Move constructor using SpinPolarizedData.
   * @param orig Initial/original value.
   */
  inline SpinPolarizedData(      SpinPolarizedData<RESTRICTED, U, void>&& orig) : U(std::move(orig)){}
  /**
   * @brief Assignment operator.
   * @param rhs Other object.
   */
  SpinPolarizedData<RESTRICTED, U, void>& operator=(const SpinPolarizedData<RESTRICTED, U, void>& rhs){
    if (this != &rhs)  this->U::operator=(rhs);
    return *this;
  }
  /**
   * @brief Move assignment operator.
   * @param rhs Other object.
   */
  SpinPolarizedData<RESTRICTED, U, void>& operator=(SpinPolarizedData<RESTRICTED, U, void>&& rhs){
    if (this != &rhs)  this->U::operator=(std::move(rhs));
    return *this;
  }

  /**
   * @brief Function to get total value (see unrestricted implementations).
   * @return Returns the total value.
   */
  inline U total() { return *this; }

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
class SpinPolarizedData<RESTRICTED,U,
    typename std::enable_if<!std::is_class<U>::value>::type> {
public:
  // Check whether the constructors have to be explicit (if that is possible at all in a sensible way)
  inline SpinPolarizedData() {};
  inline SpinPolarizedData(const U& initValue) : _m(initValue) {};
  inline operator const U&() const {return _m;};
  inline operator U&() {return _m;};

  inline U total() { return _m; }

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

  inline SpinPolarizedData(
    const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& orig) :
      std::pair<U, U>((const std::pair<U, U>&)orig), alpha(this->first), beta(this->second) {}
  inline SpinPolarizedData(
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& orig) :
      std::pair<U, U>((std::pair<U, U>&)orig), alpha(this->first), beta(this->second) {}
  inline SpinPolarizedData(
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&& orig) :
      std::pair<U, U>(std::move(orig)), alpha(this->first), beta(this->second) {}
  inline SpinPolarizedData(std::piecewise_construct_t, U&& alpha, U&& beta) :
      std::pair<U, U>(std::forward<U>(alpha), std::forward<U>(beta)) {}
  inline SpinPolarizedData(std::piecewise_construct_t, const U& alpha, const U& beta) :
      std::pair<U, U>(alpha, beta) {}

  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& operator=(
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& orig) {
    this->alpha = orig.alpha;
    this->beta = orig.beta;
    return *this;
  }
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& operator=(
      SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&& orig) {
    this->alpha = std::move(orig.alpha);
    this->beta = std::move(orig.beta);
    return *this;
  }
  /*
   * Forward constructors of a single U-type object. I looked at std::make_shared for this.
   */
  template<typename... _Args>inline SpinPolarizedData(_Args&&... __args) :
    std::pair<U, U>(std::piecewise_construct,
          std::forward_as_tuple(__args...),
          std::forward_as_tuple(__args...)) {}
  U& alpha = this->first;
  U& beta = this->second;

  inline U total() { return this->alpha + this->beta; }

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
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>& operator+= (
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) {
    this->alpha += y.alpha;
    this->beta += y.beta;
    return *this;
  }
  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> operator+ (
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) const {
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> result(*this);
    result += y;
    return result;
  }

  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>& operator*= (
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) {
    this->alpha *= y.alpha;
    this->beta *= y.beta;
    return *this;
  }
  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> operator* (
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) const {
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> result(*this);
    result *= y;
    return result;
  }

  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>& operator-= (
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) {
    this->alpha -= y.alpha;
    this->beta -= y.beta;
    return *this;
  }
  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> operator- (
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) const {
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> result(*this);
    result -= y;
    return result;
  }

  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>& operator/= (
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) {
    this->alpha /= y.alpha;
    this->beta /= y.beta;
    return *this;
  }
  template<class Y>
  inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> operator/ (
      const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Y>& y) const {
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> result(*this);
    result /= y;
    return result;
  }

};

template<class U>
inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> makeUnrestrictedFromPieces(
          U&& alpha, U&& beta) {
  return SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>(
      std::piecewise_construct, std::forward<U>(alpha), std::forward<U>(beta));
}
template<class U>
inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U> makeUnrestrictedFromPieces(
          const U& alpha, const U& beta) {
  return SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U>(
      std::piecewise_construct, alpha, beta);
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
template<Options::SCF_MODES T, class U, bool isConstU>class ForSpinHelperOne;
template<class U, bool isConstU>class ForSpinHelperOne<Options::SCF_MODES::RESTRICTED, U, isConstU> {
  public:
  inline ForSpinHelperOne(
      typename std::conditional<isConstU,const U&,U&>::type u) :
    _u(u) {};
  void operator<<(const std::function <void (
         typename std::conditional<isConstU, const U&, U&>::type&)>& f) {
    f(_u);
  }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
};
template<class U, bool isConstU>class ForSpinHelperOne<Options::SCF_MODES::UNRESTRICTED, U, isConstU> {
  public:
  inline ForSpinHelperOne(typename std::conditional<isConstU,const U&,U&>::type u) :
    _u(u) {};
  void operator<<(const std::function <void (
         typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type&)>& f) {
    f(_u.alpha);
    f(_u.beta);
  }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type typeU;
};
#define _FOR_SPIN_1(a) ForSpinHelperOne<                                                           \
    std::remove_reference<decltype(a)>::type::scf_mode,                                            \
    typename std::remove_reference<decltype(a)>::type::type,                                       \
    std::is_const<typename std::remove_reference<decltype(a)>::type>::value>(a) << [&](            \
      typename ForSpinHelperOne<                                                                   \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value>::typeU& a ## _spin)

/*
 * Two arguments
 */
template<Options::SCF_MODES T, class U, class V, bool isConstU, bool isConstV>class ForSpinHelperTwo;
template<class U, class V, bool isConstU, bool isConstV>
class ForSpinHelperTwo<Options::SCF_MODES::RESTRICTED, U, V, isConstU, isConstV> {
  public:
  inline ForSpinHelperTwo(
      typename std::conditional<isConstU,const U&,U&>::type u,
      typename std::conditional<isConstV,const V&,V&>::type v) :
    _u(u), _v(v) {};
  void operator<<(const std::function <void (
          typename std::conditional<isConstU, const U&, U&>::type&,
          typename std::conditional<isConstV, const V&, V&>::type&)>& f){
    f(_u, _v);
   }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typename std::conditional<isConstV,const V&,V&>::type _v;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
};
template<class U, class V, bool isConstU, bool isConstV>
class ForSpinHelperTwo<Options::SCF_MODES::UNRESTRICTED, U, V, isConstU, isConstV> {
  public:
  inline ForSpinHelperTwo(
      typename std::conditional<isConstU,const U&,U&>::type u,
      typename std::conditional<isConstV,const V&,V&>::type v) :
    _u(u), _v(v) {};
  void operator<<(const std::function <void (
         typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type& ,
         typename std::conditional<isConstV, typename V::constspinlesstype,typename V::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha);
    f(_u.beta, _v.beta);
  }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typename std::conditional<isConstV,const V&,V&>::type _v;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype,typename V::spinlesstype>::type typeV;
};
#define _FOR_SPIN_2(a, b) ForSpinHelperTwo<                                                        \
    std::remove_reference<decltype(a)>::type::scf_mode,                                            \
    typename std::remove_reference<decltype(a)>::type::type,                                       \
    typename std::remove_reference<decltype(b)>::type::type,                                       \
    std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(b)>::type>::value>(a, b) << [&](         \
      typename ForSpinHelperTwo<                                                                   \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value>::typeU& a ## _spin,\
      typename ForSpinHelperTwo<                                                                   \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value>::typeV& b ## _spin)

/*
 * Three arguments
 */
template<Options::SCF_MODES T, class U, class V, class W, bool isConstU, bool isConstV, bool isConstW>
class ForSpinHelperThree;
template<class U, class V, class W, bool isConstU, bool isConstV, bool isConstW>
class ForSpinHelperThree<Options::SCF_MODES::RESTRICTED, U, V, W, isConstU, isConstV, isConstW> {
  public:
  inline ForSpinHelperThree(
      typename std::conditional<isConstU,const U&,U&>::type u,
      typename std::conditional<isConstV,const V&,V&>::type v,
      typename std::conditional<isConstW,const W&,W&>::type w) :
    _u(u), _v(v), _w(w) {};
  void operator<<(const std::function <void (
          typename std::conditional<isConstU, const U&, U&>::type&,
          typename std::conditional<isConstV, const V&, V&>::type&,
          typename std::conditional<isConstW, const W&, W&>::type&)>& f){
    f(_u, _v, _w);
   }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typename std::conditional<isConstV,const V&,V&>::type _v;
  typename std::conditional<isConstW,const W&,W&>::type _w;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
};
template<class U, class V, class W, bool isConstU, bool isConstV, bool isConstW>
class ForSpinHelperThree<Options::SCF_MODES::UNRESTRICTED, U, V, W, isConstU, isConstV, isConstW> {
  public:
  inline ForSpinHelperThree(
      typename std::conditional<isConstU,const U&,U&>::type u,
      typename std::conditional<isConstV,const V&,V&>::type v,
      typename std::conditional<isConstW,const W&,W&>::type w) :
    _u(u), _v(v), _w(w) {};
  void operator<<(const std::function <void (
         typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type& ,
         typename std::conditional<isConstV, typename V::constspinlesstype,typename V::spinlesstype>::type& ,
         typename std::conditional<isConstW, typename W::constspinlesstype,typename W::spinlesstype>::type&)>& f){
    f(_u.alpha, _v.alpha, _w.alpha);
    f(_u.beta, _v.beta, _w.beta);
  }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typename std::conditional<isConstV,const V&,V&>::type _v;
  typename std::conditional<isConstW,const W&,W&>::type _w;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype,typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype,typename W::spinlesstype>::type typeW;
};
#define _FOR_SPIN_3(a, b, c) ForSpinHelperThree<                                                   \
    std::remove_reference<decltype(a)>::type::scf_mode,                                            \
    typename std::remove_reference<decltype(a)>::type::type,                                       \
    typename std::remove_reference<decltype(b)>::type::type,                                       \
    typename std::remove_reference<decltype(c)>::type::type,                                       \
    std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(c)>::type>::value>(a, b, c) << [&](      \
      typename ForSpinHelperThree<                                                                 \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value>::typeU& a ## _spin,\
      typename ForSpinHelperThree<                                                                 \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value>::typeV& b ## _spin,\
      typename ForSpinHelperThree<                                                                 \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value>::typeW& c ## _spin)


/*
 * Four arguments
 */
template<Options::SCF_MODES T, class U, class V, class W, class X,
        bool isConstU, bool isConstV, bool isConstW, bool isConstX>class ForSpinHelperFour;
template<class U, class V, class W, class X,
    bool isConstU, bool isConstV, bool isConstW, bool isConstX>
class ForSpinHelperFour<Options::SCF_MODES::RESTRICTED, U, V, W, X,
        isConstU, isConstV, isConstW, isConstX> {
  public:
  inline ForSpinHelperFour(
      typename std::conditional<isConstU,const U&,U&>::type u,
      typename std::conditional<isConstV,const V&,V&>::type v,
      typename std::conditional<isConstW,const W&,W&>::type w,
      typename std::conditional<isConstX,const X&,X&>::type x):
    _u(u), _v(v), _w(w), _x(x) {};
  void operator<<(const std::function <void (
          typename std::conditional<isConstU, const U&, U&>::type&,
          typename std::conditional<isConstV, const V&, V&>::type&,
          typename std::conditional<isConstW, const W&, W&>::type&,
          typename std::conditional<isConstX, const X&, X&>::type&)>& f) {
    f(_u, _v, _w, _x);
   }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typename std::conditional<isConstV,const V&,V&>::type _v;
  typename std::conditional<isConstW,const W&,W&>::type _w;
  typename std::conditional<isConstX,const X&,X&>::type _x;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
};
template<class U, class V, class W, class X,
    bool isConstU, bool isConstV, bool isConstW, bool isConstX>
class ForSpinHelperFour<Options::SCF_MODES::UNRESTRICTED, U, V, W, X,
        isConstU, isConstV, isConstW, isConstX> {
  public:
  inline ForSpinHelperFour(
      typename std::conditional<isConstU,const U&,U&>::type u,
      typename std::conditional<isConstV,const V&,V&>::type v,
      typename std::conditional<isConstW,const W&,W&>::type w,
      typename std::conditional<isConstX,const X&,X&>::type x) :
    _u(u), _v(v), _w(w), _x(x) {};
  void operator<<(const std::function <void (
         typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type& ,
         typename std::conditional<isConstV, typename V::constspinlesstype,typename V::spinlesstype>::type& ,
         typename std::conditional<isConstW, typename W::constspinlesstype,typename W::spinlesstype>::type& ,
         typename std::conditional<isConstX, typename X::constspinlesstype,typename X::spinlesstype>::type&)>& f){
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta);
  }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typename std::conditional<isConstV,const V&,V&>::type _v;
  typename std::conditional<isConstW,const W&,W&>::type _w;
  typename std::conditional<isConstX,const X&,X&>::type _x;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype,typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype,typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype,typename X::spinlesstype>::type typeX;
};
#define _FOR_SPIN_4(a, b, c, d) ForSpinHelperFour<                                                 \
    std::remove_reference<decltype(a)>::type::scf_mode,                                            \
    typename std::remove_reference<decltype(a)>::type::type,                                       \
    typename std::remove_reference<decltype(b)>::type::type,                                       \
    typename std::remove_reference<decltype(c)>::type::type,                                       \
    typename std::remove_reference<decltype(d)>::type::type,                                       \
    std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(d)>::type>::value>(a, b, c, d) << [&](   \
      typename ForSpinHelperFour<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value>::typeU& a ## _spin,\
      typename ForSpinHelperFour<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value>::typeV& b ## _spin,\
      typename ForSpinHelperFour<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value>::typeW& c ## _spin,\
      typename ForSpinHelperFour<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value>::typeX& d ## _spin)

template<Options::SCF_MODES T, class U, class V, class W, class X, class Y,
        bool isConstU, bool isConstV, bool isConstW, bool isConstX, bool isConstY>class ForSpinHelperFive;
template<class U, class V, class W, class X, class Y,
    bool isConstU, bool isConstV, bool isConstW, bool isConstX, bool isConstY>
class ForSpinHelperFive<Options::SCF_MODES::RESTRICTED, U, V, W, X, Y,
        isConstU, isConstV, isConstW, isConstX, isConstY> {
  public:
  inline ForSpinHelperFive(
      typename std::conditional<isConstU,const U&,U&>::type u,
      typename std::conditional<isConstV,const V&,V&>::type v,
      typename std::conditional<isConstW,const W&,W&>::type w,
      typename std::conditional<isConstX,const X&,X&>::type x,
      typename std::conditional<isConstY,const Y&,Y&>::type y) :
    _u(u), _v(v), _w(w), _x(x), _y(y) {};
  void operator<<(const std::function <void (
          typename std::conditional<isConstU, const U&, U&>::type&,
          typename std::conditional<isConstV, const V&, V&>::type&,
          typename std::conditional<isConstW, const W&, W&>::type&,
          typename std::conditional<isConstX, const X&, X&>::type&,
          typename std::conditional<isConstY, const Y&, Y&>::type&)>& f) {
    f(_u, _v, _w, _x, _y);
   }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typename std::conditional<isConstV,const V&,V&>::type _v;
  typename std::conditional<isConstW,const W&,W&>::type _w;
  typename std::conditional<isConstX,const X&,X&>::type _x;
  typename std::conditional<isConstY,const Y&,Y&>::type _y;
  typedef typename std::conditional<isConstU, const U&, U&>::type typeU;
  typedef typename std::conditional<isConstV, const V&, V&>::type typeV;
  typedef typename std::conditional<isConstW, const W&, W&>::type typeW;
  typedef typename std::conditional<isConstX, const X&, X&>::type typeX;
  typedef typename std::conditional<isConstY, const Y&, Y&>::type typeY;
};
template<class U, class V, class W, class X, class Y,
    bool isConstU, bool isConstV, bool isConstW, bool isConstX, bool isConstY>
class ForSpinHelperFive<Options::SCF_MODES::UNRESTRICTED, U, V, W, X, Y,
        isConstU, isConstV, isConstW, isConstX, isConstY> {
  public:
  inline ForSpinHelperFive(
      typename std::conditional<isConstU,const U&,U&>::type u,
      typename std::conditional<isConstV,const V&,V&>::type v,
      typename std::conditional<isConstW,const W&,W&>::type w,
      typename std::conditional<isConstX,const X&,X&>::type x,
      typename std::conditional<isConstY,const Y&,Y&>::type y) :
    _u(u), _v(v), _w(w), _x(x), _y(y) {};
  void operator<<(const std::function <void (
         typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type& ,
         typename std::conditional<isConstV, typename V::constspinlesstype,typename V::spinlesstype>::type& ,
         typename std::conditional<isConstW, typename W::constspinlesstype,typename W::spinlesstype>::type& ,
         typename std::conditional<isConstX, typename X::constspinlesstype,typename X::spinlesstype>::type& ,
         typename std::conditional<isConstY, typename Y::constspinlesstype,typename Y::spinlesstype>::type&)>& f) {
    f(_u.alpha, _v.alpha, _w.alpha, _x.alpha, _y.alpha);
    f(_u.beta, _v.beta, _w.beta, _x.beta, _y.beta);
  }
  typename std::conditional<isConstU,const U&,U&>::type _u;
  typename std::conditional<isConstV,const V&,V&>::type _v;
  typename std::conditional<isConstW,const W&,W&>::type _w;
  typename std::conditional<isConstX,const X&,X&>::type _x;
  typename std::conditional<isConstY,const Y&,Y&>::type _y;
  typedef typename std::conditional<isConstU, typename U::constspinlesstype,typename U::spinlesstype>::type typeU;
  typedef typename std::conditional<isConstV, typename V::constspinlesstype,typename V::spinlesstype>::type typeV;
  typedef typename std::conditional<isConstW, typename W::constspinlesstype,typename W::spinlesstype>::type typeW;
  typedef typename std::conditional<isConstX, typename X::constspinlesstype,typename X::spinlesstype>::type typeX;
  typedef typename std::conditional<isConstY, typename Y::constspinlesstype,typename Y::spinlesstype>::type typeY;
};
#define _FOR_SPIN_5(a, b, c, d, e) ForSpinHelperFive<                                              \
    std::remove_reference<decltype(a)>::type::scf_mode,                                            \
    typename std::remove_reference<decltype(a)>::type::type,                                       \
    typename std::remove_reference<decltype(b)>::type::type,                                       \
    typename std::remove_reference<decltype(c)>::type::type,                                       \
    typename std::remove_reference<decltype(d)>::type::type,                                       \
    typename std::remove_reference<decltype(e)>::type::type,                                       \
    std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                       \
    std::is_const<typename std::remove_reference<decltype(e)>::type>::value>(a, b, c, d, e) << [&](\
      typename ForSpinHelperFive<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        typename std::remove_reference<decltype(e)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeU& a ## _spin,\
      typename ForSpinHelperFive<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        typename std::remove_reference<decltype(e)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeV& b ## _spin,\
      typename ForSpinHelperFive<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        typename std::remove_reference<decltype(e)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeW& c ## _spin,\
      typename ForSpinHelperFive<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        typename std::remove_reference<decltype(e)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeX& d ## _spin,\
      typename ForSpinHelperFive<                                                                  \
        std::remove_reference<decltype(a)>::type::scf_mode,                                        \
        typename std::remove_reference<decltype(a)>::type::type,                                   \
        typename std::remove_reference<decltype(b)>::type::type,                                   \
        typename std::remove_reference<decltype(c)>::type::type,                                   \
        typename std::remove_reference<decltype(d)>::type::type,                                   \
        typename std::remove_reference<decltype(e)>::type::type,                                   \
        std::is_const<typename std::remove_reference<decltype(a)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(b)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(c)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(d)>::type>::value,                   \
        std::is_const<typename std::remove_reference<decltype(e)>::type>::value>::typeY& e ## _spin)

#define _GET_FOR_SPIN_MACRO(_1,_2,_3,_4,_5,NAME,...) NAME
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
#define for_spin(...) _GET_FOR_SPIN_MACRO(__VA_ARGS__, _FOR_SPIN_5, _FOR_SPIN_4, _FOR_SPIN_3, _FOR_SPIN_2, _FOR_SPIN_1, DUMMY)(__VA_ARGS__)
// Note: the final argument DUMMY is only there to avoid a compiler warning.
} /* namespace Serenity */
#endif	/* SpinPolarizedData_H */
