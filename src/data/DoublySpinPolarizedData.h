/**
 * @file   DoublySpinPolarizedData.h
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   4. September 2015, 15:13
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
#ifndef DOUBLYSPINPOLARIZEDDATA_H
#define DOUBLYSPINPOLARIZEDDATA_H
/* Include Serenity Internal Headers */
#include "data/grid/GridData.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>
#include <type_traits>
#include <utility>

namespace Serenity {
/* Forward declarations */
class GridController;
/**
 * @class DoublySpinPolarizedData DoublySpinPolarizedData.h
 * @brief Analog of SpinPolarizedData for a double spin polarization
 *
 * Especially when dealing with higher derivatives of data which may be spin-polarized, additional
 * combinations of spin polarizations arise. Let's take as an example a function f of the electron
 * density rho. In the spin-unpolarized (restricted) case f is only dependent on one variable, so
 * for each derivative level one only gets one new function: f'(rho), f''(rho) and so on. With spin-
 * polarization f is actually dependent on two variables, rho_alpha and rho_beta. So the derivatives
 * can also be taken w.r.t. these two different variables and one obtaines for the first derivative:
 * df/drho_alpha and df/drho_beta. Taking higher derivatives more combinations arise:
 * d^2f/drho_alpha^2, d^2f/drho_beta^2, d^2f/(drho_alpha drho_beta) = d^2f/(drho_beta drho_alpha).
 * Objects of this class gather e.g. all these three derivatives together.
 */
template<Options::SCF_MODES SCFMode, class U, typename E = void>
class DoublySpinPolarizedData;
/**
 * Restricted, non-primitives
 */
template<class U>
class DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, U, typename std::enable_if<std::is_class<U>::value>::type>
  : public U {
 public:
  using U::U;
  using U::operator=;
};
/**
 * Restricted, primitives (requires special treatment)
 */
template<class U>
class DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, U, typename std::enable_if<!std::is_class<U>::value>::type> {
 public:
  inline DoublySpinPolarizedData(){};
  inline DoublySpinPolarizedData(U initValue) : _m(initValue){};
  inline operator const U&() const {
    return _m;
  };
  inline operator U&() {
    return _m;
  };

  U _m;
};
/**
 * Unrestricted
 */
template<class U>
class DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void> {
 public:
  inline explicit DoublySpinPolarizedData(DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&) = default;
  inline explicit DoublySpinPolarizedData(const DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&) = default;
  inline DoublySpinPolarizedData(DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&&) = default;
  inline DoublySpinPolarizedData(std::piecewise_construct_t, U&& aa_, U&& ab_, U&& ba_, U&& bb_)
    : aa(std::move(aa_)), ab(std::move(ab_)), ba(std::move(ba_)), bb(std::move(bb_)) {
  }

  inline DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&
  operator=(const DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>& orig) = default;
  inline DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&
  operator=(DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, U, void>&& orig) = default;
  /*
   * Forward constructors of a single U-type object. I looked at std::make_shared for this.
   */
  template<typename... _Args>
  inline DoublySpinPolarizedData(_Args&&... __args)
    : DoublySpinPolarizedData(std::piecewise_construct, U(std::forward<_Args>(__args)...), U(std::forward<_Args>(__args)...),
                              U(std::forward<_Args>(__args)...), U(std::forward<_Args>(__args)...)) {
  }
  U aa;
  U ab;
  U ba;
  U bb;
};
/**
 * Unrestricted, special case: U=GridData, it is more convenient to have a usage analogous to the
 * SpinPolarizedGridData.
 */
template<>
class DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<RESTRICTED>, void> {
 public:
  inline explicit DoublySpinPolarizedData(
      DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<RESTRICTED>, void>&) = default;
  inline explicit DoublySpinPolarizedData(
      const DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<RESTRICTED>, void>&) = default;
  inline DoublySpinPolarizedData(DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<RESTRICTED>, void>&&) = default;
  inline DoublySpinPolarizedData(std::piecewise_construct_t, GridData<RESTRICTED>&& aa_, GridData<RESTRICTED>&& ab_,
                                 GridData<RESTRICTED>&& ba_, GridData<RESTRICTED>&& bb_)
    : aa(std::move(aa_)), ab(std::move(ab_)), ba(std::move(ba_)), bb(std::move(bb_)) {
  }

  inline DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<RESTRICTED>, void>&
  operator=(const DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<RESTRICTED>, void>& orig) = default;
  inline DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<RESTRICTED>, void>&
  operator=(DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<RESTRICTED>, void>&& orig) = default;
  /*
   * Forward constructors of a single U-type object. I looked at std::make_shared for this.
   */
  template<typename... _Args>
  inline DoublySpinPolarizedData(_Args&&... __args)
    : DoublySpinPolarizedData(std::piecewise_construct, GridData<RESTRICTED>(std::forward<_Args>(__args)...),
                              GridData<RESTRICTED>(std::forward<_Args>(__args)...),
                              GridData<RESTRICTED>(std::forward<_Args>(__args)...),
                              GridData<RESTRICTED>(std::forward<_Args>(__args)...)) {
  }
  GridData<RESTRICTED> aa;
  GridData<RESTRICTED> ab;
  GridData<RESTRICTED> ba;
  GridData<RESTRICTED> bb;

  inline std::shared_ptr<GridController> getGridController() const {
    assert(aa.getGridController() == ab.getGridController());
    assert(aa.getGridController() == ba.getGridController());
    assert(aa.getGridController() == bb.getGridController());
    return aa.getGridController();
  }

  inline bool isValid() const {
    return (aa.isValid() && ab.isValid() && ba.isValid() && bb.isValid());
  }
};
template<class T, class X>
void apply(DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, T>& target, X func) {
  func(target);
}
template<class T, class X>
void apply(const DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, T>& target, X func) {
  func(target);
}
template<class T, class X>
void apply(DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, T>& target, X func) {
  func(target.aa);
  func(target.ab);
  func(target.ba);
  func(target.bb);
}
template<class T, class X>
void apply(const DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, T>& target, X func) {
  func(target.aa);
  func(target.ab);
  func(target.ba);
  func(target.bb);
}

} /* namespace Serenity */

#endif /* DOUBLYSPINPOLARIZEDDATA_H */
