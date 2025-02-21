/**
 * @file PartialDerivatives.h
 *
 * @date Apr 4, 2017
 * @author Michael Boeckers
 * @date Jul 12, 2015
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

#ifndef DFT_FUNCTIONALS_WRAPPERS_PARTIALDERIVATIVES_H_
#define DFT_FUNCTIONALS_WRAPPERS_PARTIALDERIVATIVES_H_

/* Include Serenity Internal Headers */
#include "data/DoublySpinPolarizedData.h"
#include "data/grid/GridData.h"
#include "data/grid/GridPotential.h"
#include "dft/Functional.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "math/Derivatives.h"
#include "settings/Options.h"

namespace Serenity {
/**
 * @class  dF_dRho PartialDerivatives.h
 * @brief Object to store the derivative of the enhancement factor F w.r.t. the density, i.e.
 *        \f$ \frac{\partial F}{\partial \rho}\f$. In the unrestricted case, this object has
 *        two elements which can be accessed as usual (See SpinPolarizedData).
 *        \n
 *        Previously, this derivative was stored as GridPotential. However, this may be confusing
 *        as one may think this is the complete XC-potential which is only true when using pure LDA
 *        functionals.
 */
template<Options::SCF_MODES SCFMode>
class dF_dRho : public GridPotential<SCFMode> {
 public:
  using GridPotential<SCFMode>::GridPotential;
  explicit dF_dRho(dF_dRho<SCFMode>&) = default;
  explicit dF_dRho(const dF_dRho<SCFMode>&) = default;
  dF_dRho(dF_dRho<SCFMode>&&) = default;
  dF_dRho& operator=(const dF_dRho<SCFMode>&) = default;
  dF_dRho& operator=(dF_dRho<SCFMode>&&) = default;
  using GridPotential<SCFMode>::operator+=;
  using GridPotential<SCFMode>::operator-=;
};

/**
 * @class   d2F_dRho2 PartialDerivatives.h
 * @brief   Object to store the second derivative of the enhancement factor F w.r.t. the density, i.e.
 *          \f$ \frac{\partial^2 F}{\partial \rho^2}\f$. In the unrestricted case, this object holds
 *          three elements, \f$ \frac{\partial^2 F}{\partial \rho_\alpha \rho_\alpha}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \rho_\alpha \partial\rho_\beta} \f$ and
 *           \f$ \frac{\partial^2 F}{\partial \rho_\beta \partial\rho_\beta} \f$.
 *          Every distinct derivative is stored as GridData.
 */
template<Options::SCF_MODES SCFMode>
class d2F_dRho2;

template<>
class d2F_dRho2<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
 public:
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<>
class d2F_dRho2<Options::SCF_MODES::UNRESTRICTED> {
 public:
  /**
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  d2F_dRho2(std::shared_ptr<GridController> gridController)
    : aa(gridController), ab(gridController), bb(gridController) {
  }

  GridData<RESTRICTED> aa;
  GridData<RESTRICTED> ab;
  GridData<RESTRICTED> bb;
};

/**
 * @class   dF_dSigma PartialDerivatives.h
 * @brief   Object to store the derivative of the enhancement factor F w.r.t. \f$ \sigma \f$, i.e.
 *          \f$ \frac{\partial F}{\partial \sigma}\f$ where \f$ \sigma =  (\nabla \rho)^2 \f$. In the
 *          unrestricted case this objects holds three elements, i.e. \f$ \frac{\partial F}{\partial
 * \sigma_{\alpha\alpha}}\f$, \f$ \frac{\partial F}{\partial \sigma_{\alpha\beta}}\f$ and \f$ \frac{\partial F}{\partial
 * \sigma_{\beta\beta}}\f$ . Every distinct derivative is stored as GridData.
 *
 */
template<Options::SCF_MODES SCFMode>
class dF_dSigma;

template<>
class dF_dSigma<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<>
class dF_dSigma<Options::SCF_MODES::UNRESTRICTED> {
 public:
  /**
   *
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  dF_dSigma(std::shared_ptr<GridController> gridController)
    : aa(gridController), ab(gridController), bb(gridController) {
  }

  GridData<RESTRICTED> aa;
  GridData<RESTRICTED> ab;
  GridData<RESTRICTED> bb;
};

/**
 * @class   d2F_dSigma2 PartialDerivatives.h
 * @brief   Object to store the second derivative of the enhancement factor F w.r.t. \f$ \sigma \f$, i.e.
 *          \f$ \frac{\partial^2 F}{\partial \sigma^2}\f$ where \f$ \sigma =  (\nabla \rho)^2 \f$. In the
 *          unrestricted case this objects holds six elements, i.e.
 *          \f$ \frac{\partial^2 F}{\partial \sigma_{\alpha\alpha} \partial\sigma_{\alpha\alpha}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \sigma_{\alpha\alpha} \partial\sigma_{\alpha\beta}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \sigma_{\alpha\alpha} \partial\sigma_{\beta\beta}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \sigma_{\alpha\beta} \partial\sigma_{\alpha\beta}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \sigma_{\alpha\beta} \partial\sigma_{\beta\beta}}\f$ and
 *          \f$ \frac{\partial^2 F}{\partial \sigma_{\beta \beta} \partial\sigma_{\beta\beta}}\f$.
 *          Every distinct derivative is stored as GridData.
 *          Note that the first two letters of the GridData object xxxx within this class always refer
 *          to the first derivative w.r.t. \f$ \sigma \f$, while the latter two refer to the second derivative.
 */
template<Options::SCF_MODES SCFMode>
class d2F_dSigma2;

template<>
class d2F_dSigma2<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<>
class d2F_dSigma2<Options::SCF_MODES::UNRESTRICTED> {
 public:
  /**
   *
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  d2F_dSigma2(std::shared_ptr<GridController> gridController)
    : aaaa(gridController), aaab(gridController), aabb(gridController), abab(gridController), abbb(gridController), bbbb(gridController) {
  }

  GridData<RESTRICTED> aaaa;
  GridData<RESTRICTED> aaab;
  GridData<RESTRICTED> aabb;
  GridData<RESTRICTED> abab;
  GridData<RESTRICTED> abbb;
  GridData<RESTRICTED> bbbb;
};

/**
 * @class   d2F_dRhodSigma PartialDerivatives.h
 * @brief   Object to store the second derivative of the enhancement factor F w.r.t. \f$ \rho \f$ and \f$ \sigma \f$,
 *           i.e. \f$ \frac{\partial^2 F}{\partial \rho \partial \sigma}\f$ where \f$ \sigma =  (\nabla \rho)^2 \f$. In
 * the unrestricted case this objects holds six elements, i.e. \f$ \frac{\partial^2 F}{\partial \rho_\alpha
 * \partial\sigma_{\alpha\alpha}}\f$, \f$ \frac{\partial^2 F}{\partial \rho_\alpha \partial\sigma_{\alpha\beta}}\f$, \f$
 * \frac{\partial^2 F}{\partial \rho_\alpha \partial\sigma_{\beta\beta}}\f$, \f$ \frac{\partial^2 F}{\partial \rho_\beta
 * \partial\sigma_{\alpha\alpha}}\f$, \f$ \frac{\partial^2 F}{\partial \rho_\beta \partial\sigma_{\alpha\beta}}\f$, \f$
 * \frac{\partial^2 F}{\partial \rho_\beta \partial\sigma_{\beta\beta}}\f$. Every distinct derivative is stored as
 * GridData. Note that the first letter of the GridData object xxx within this class always refers to the derivative
 * w.r.t. the density while the latter two determine the spin-indices of \f$ \sigma \f$.
 */
template<Options::SCF_MODES SCFMode>
class d2F_dRhodSigma;

template<>
class d2F_dRhodSigma<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<>
class d2F_dRhodSigma<Options::SCF_MODES::UNRESTRICTED> {
 public:
  /**
   *
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  d2F_dRhodSigma(std::shared_ptr<GridController> gridController)
    : aaa(gridController), aab(gridController), abb(gridController), baa(gridController), bab(gridController), bbb(gridController) {
  }

  GridData<RESTRICTED> aaa;
  GridData<RESTRICTED> aab;
  GridData<RESTRICTED> abb;
  GridData<RESTRICTED> baa;
  GridData<RESTRICTED> bab;
  GridData<RESTRICTED> bbb;
};

/**
 * @class   d3F_dRho2dSigma PartialDerivatives.h
 * @brief   Object to store the third derivative of the enhancement factor F w.r.t. the density twice and \f$ \sigma \f$
 * once, i.e. \f$ \frac{\partial^3 F}{\partial \rho^2 \partial\sigma}\f$. In the unrestricted case, this object holds 9
 * elements, which are each stored as GridData. The first indices always correspond to Rho, the remaining to sigma.
 *
 */
template<Options::SCF_MODES SCFMode>
class d3F_dRho2dSigma;

template<>
class d3F_dRho2dSigma<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
 public:
  using GridData<RESTRICTED>::GridData;
};

template<>
class d3F_dRho2dSigma<Options::SCF_MODES::UNRESTRICTED> {
 public:
  d3F_dRho2dSigma(std::shared_ptr<GridController> gridController)
    : aaaa(gridController),
      aaab(gridController),
      aabb(gridController),
      abaa(gridController),
      abab(gridController),
      abbb(gridController),
      bbaa(gridController),
      bbab(gridController),
      bbbb(gridController) {
  }

  GridData<RESTRICTED> aaaa;
  GridData<RESTRICTED> aaab;
  GridData<RESTRICTED> aabb;
  GridData<RESTRICTED> abaa;
  GridData<RESTRICTED> abab;
  GridData<RESTRICTED> abbb;
  GridData<RESTRICTED> bbaa;
  GridData<RESTRICTED> bbab;
  GridData<RESTRICTED> bbbb;
};

/**
 * @class   d3F_dRhodSigma2 PartialDerivatives.h
 * @brief   Object to store the third derivative of the enhancement factor F w.r.t. the density once and \f$ \sigma \f$
 * twice, i.e. \f$ \frac{\partial^3 F}{\partial \rho \partial\sigma^2}\f$. In the unrestricted case, this object holds
 *          12 elements, which are each stored as GridData. The first indices always correspond to Rho, the remaining to
 * sigma.
 *
 */
template<Options::SCF_MODES SCFMode>
class d3F_dRhodSigma2;

template<>
class d3F_dRhodSigma2<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
 public:
  using GridData<RESTRICTED>::GridData;
};

template<>
class d3F_dRhodSigma2<Options::SCF_MODES::UNRESTRICTED> {
 public:
  d3F_dRhodSigma2(std::shared_ptr<GridController> gridController)
    : aaaaa(gridController),
      aaaab(gridController),
      aaabb(gridController),
      aabab(gridController),
      aabbb(gridController),
      abbbb(gridController),
      baaaa(gridController),
      baaab(gridController),
      baabb(gridController),
      babab(gridController),
      babbb(gridController),
      bbbbb(gridController) {
  }

  GridData<RESTRICTED> aaaaa;
  GridData<RESTRICTED> aaaab;
  GridData<RESTRICTED> aaabb;
  GridData<RESTRICTED> aabab;
  GridData<RESTRICTED> aabbb;
  GridData<RESTRICTED> abbbb;
  GridData<RESTRICTED> baaaa;
  GridData<RESTRICTED> baaab;
  GridData<RESTRICTED> baabb;
  GridData<RESTRICTED> babab;
  GridData<RESTRICTED> babbb;
  GridData<RESTRICTED> bbbbb;
};

/**
 * @class   d3F_dSigma3 PartialDerivatives.h
 * @brief   Object to store the third derivative of the enhancement factor F w.r.t.  \f$ \sigma \f$ thrice, i.e.
 *          \f$ \frac{\partial^3 F}{\partial \partial\sigma^3}\f$. In the unrestricted case, this object holds
 *          10 elements, which are each stored as GridData.
 *
 */
template<Options::SCF_MODES SCFMode>
class d3F_dSigma3;

template<>
class d3F_dSigma3<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
 public:
  using GridData<RESTRICTED>::GridData;
};

template<>
class d3F_dSigma3<Options::SCF_MODES::UNRESTRICTED> {
 public:
  d3F_dSigma3(std::shared_ptr<GridController> gridController)
    : aaaaaa(gridController),
      aaaaab(gridController),
      aaaabb(gridController),
      aaabab(gridController),
      aaabbb(gridController),
      aabbbb(gridController),
      ababab(gridController),
      ababbb(gridController),
      abbbbb(gridController),
      bbbbbb(gridController) {
  }

  GridData<RESTRICTED> aaaaaa;
  GridData<RESTRICTED> aaaaab;
  GridData<RESTRICTED> aaaabb;
  GridData<RESTRICTED> aaabab;
  GridData<RESTRICTED> aaabbb;
  GridData<RESTRICTED> aabbbb;
  GridData<RESTRICTED> ababab;
  GridData<RESTRICTED> ababbb;
  GridData<RESTRICTED> abbbbb;
  GridData<RESTRICTED> bbbbbb;
};

/**
 * @class   TriplySpinPolarizedGridData PartialDerivatives.h
 * @brief   Helper class to organize member variables for triply spin-polarized objects with three distinct spin
 * polarizations, e.g. d3F_dGradRho3.xyz.
 */
class TriplySpinPolarizedGridData {
 public:
  TriplySpinPolarizedGridData(std::shared_ptr<GridController> gridController)
    : aaa(gridController),
      aab(gridController),
      aba(gridController),
      abb(gridController),
      baa(gridController),
      bab(gridController),
      bba(gridController),
      bbb(gridController) {
  }

  GridData<RESTRICTED> aaa;
  GridData<RESTRICTED> aab;
  GridData<RESTRICTED> aba;
  GridData<RESTRICTED> abb;
  GridData<RESTRICTED> baa;
  GridData<RESTRICTED> bab;
  GridData<RESTRICTED> bba;
  GridData<RESTRICTED> bbb;

  /**
   * @brief += operator for other TriplySpinPolarizedGridData.
   */
  inline void operator+=(const TriplySpinPolarizedGridData& rhs) {
    this->aaa += rhs.aaa;
    this->aab += rhs.aab;
    this->aba += rhs.aba;
    this->abb += rhs.abb;
    this->baa += rhs.baa;
    this->bab += rhs.bab;
    this->bba += rhs.bba;
    this->bbb += rhs.bbb;
  }

  /**
   * @brief Sets all elements to zero at iPoint that depend on the alpha electron density or its gradient.
   */
  inline void setAlphaZero(unsigned iPoint) {
    this->aaa(iPoint) = 0.0;
    this->aab(iPoint) = 0.0;
    this->aba(iPoint) = 0.0;
    this->abb(iPoint) = 0.0;
    this->baa(iPoint) = 0.0;
    this->bab(iPoint) = 0.0;
    this->bba(iPoint) = 0.0;
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the beta electron density or its gradient.
   */
  inline void setBetaZero(unsigned iPoint) {
    this->aab(iPoint) = 0.0;
    this->aba(iPoint) = 0.0;
    this->abb(iPoint) = 0.0;
    this->baa(iPoint) = 0.0;
    this->bab(iPoint) = 0.0;
    this->bba(iPoint) = 0.0;
    this->bbb(iPoint) = 0.0;
  }
};

/**
 * @class   TriplySpinPolarizedGridDataTwoOne PartialDerivatives.h
 * @brief   Helper class to organize member variables for triply spin-polarized objects where the first two spin
 * polarizations are indistinguishable, e.g. d3F_dRho2dGradRho.
 */
class TriplySpinPolarizedGridDataTwoOne {
 public:
  TriplySpinPolarizedGridDataTwoOne(std::shared_ptr<GridController> gridController)
    : aaa(gridController), aab(gridController), aba(gridController), abb(gridController), bba(gridController), bbb(gridController) {
  }

  GridData<RESTRICTED> aaa;
  GridData<RESTRICTED> aab;
  GridData<RESTRICTED> aba;
  GridData<RESTRICTED> abb;
  GridData<RESTRICTED> bba;
  GridData<RESTRICTED> bbb;

  /**
   * @brief += operator for other TriplySpinPolarizedGridDataTwoOne.
   */
  inline void operator+=(const TriplySpinPolarizedGridDataTwoOne& rhs) {
    this->aaa += rhs.aaa;
    this->aab += rhs.aab;
    this->aba += rhs.aba;
    this->abb += rhs.abb;
    this->bba += rhs.bba;
    this->bbb += rhs.bbb;
  }

  /**
   * @brief Sets all elements to zero at iPoint that depend on the alpha electron density or its gradient.
   */
  inline void setAlphaZero(unsigned iPoint) {
    this->aaa(iPoint) = 0.0;
    this->aab(iPoint) = 0.0;
    this->aba(iPoint) = 0.0;
    this->abb(iPoint) = 0.0;
    this->bba(iPoint) = 0.0;
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the beta electron density or its gradient.
   */
  inline void setBetaZero(unsigned iPoint) {
    this->aab(iPoint) = 0.0;
    this->aba(iPoint) = 0.0;
    this->abb(iPoint) = 0.0;
    this->bba(iPoint) = 0.0;
    this->bbb(iPoint) = 0.0;
  }
};

/**
 * @class   TriplySpinPolarizedGridDataOneTwo PartialDerivatives.h
 * @brief   Helper class to organize member variables for triply spin-polarized objects where the last two spin
 * polarizations are indistinguishable, e.g. d3F_dRhodGradRho2.xx
 */
class TriplySpinPolarizedGridDataOneTwo {
 public:
  TriplySpinPolarizedGridDataOneTwo(std::shared_ptr<GridController> gridController)
    : aaa(gridController), aab(gridController), abb(gridController), baa(gridController), bab(gridController), bbb(gridController) {
  }

  GridData<RESTRICTED> aaa;
  GridData<RESTRICTED> aab;
  GridData<RESTRICTED> abb;
  GridData<RESTRICTED> baa;
  GridData<RESTRICTED> bab;
  GridData<RESTRICTED> bbb;

  /**
   * @brief += operator for other TriplySpinPolarizedGridDataOneTwo.
   */
  inline void operator+=(const TriplySpinPolarizedGridDataOneTwo& rhs) {
    this->aaa += rhs.aaa;
    this->aab += rhs.aab;
    this->abb += rhs.abb;
    this->baa += rhs.baa;
    this->bab += rhs.bab;
    this->bbb += rhs.bbb;
  }

  /**
   * @brief Sets all elements to zero at iPoint that depend on the alpha electron density or its gradient.
   */
  inline void setAlphaZero(unsigned iPoint) {
    this->aaa(iPoint) = 0.0;
    this->aab(iPoint) = 0.0;
    this->abb(iPoint) = 0.0;
    this->baa(iPoint) = 0.0;
    this->bab(iPoint) = 0.0;
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the beta electron density or its gradient.
   */
  inline void setBetaZero(unsigned iPoint) {
    this->aab(iPoint) = 0.0;
    this->abb(iPoint) = 0.0;
    this->baa(iPoint) = 0.0;
    this->bab(iPoint) = 0.0;
    this->bbb(iPoint) = 0.0;
  }
};

/**
 * @class   TriplySpinPolarizedGridDataThree PartialDerivatives.h
 * @brief   Helper class to organize member variables for triply spin-polarized objects where the three spin
 * polarizations are all indistinguishable, e.g. d3F_dRho3. While something like this could be stored inside
 * DoublySpinPolarizedData by mappin aa -> aaa, ba -> aab, ab -> abb, bb -> bbb, this class more clearly conveys the
 * meaning.
 */
class TriplySpinPolarizedGridDataThree {
 public:
  TriplySpinPolarizedGridDataThree(std::shared_ptr<GridController> gridController)
    : aaa(gridController), aab(gridController), abb(gridController), bbb(gridController) {
  }

  GridData<RESTRICTED> aaa;
  GridData<RESTRICTED> aab;
  GridData<RESTRICTED> abb;
  GridData<RESTRICTED> bbb;

  /**
   * @brief += operator for other TriplySpinPolarizedGridDataThree.
   */
  inline void operator+=(const TriplySpinPolarizedGridDataThree& rhs) {
    this->aaa += rhs.aaa;
    this->aab += rhs.aab;
    this->abb += rhs.abb;
    this->bbb += rhs.bbb;
  }

  /**
   * @brief Sets all elements to zero at iPoint that depend on the alpha electron density or its gradient.
   */
  inline void setAlphaZero(unsigned iPoint) {
    this->aaa(iPoint) = 0.0;
    this->aab(iPoint) = 0.0;
    this->abb(iPoint) = 0.0;
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the beta electron density or its gradient.
   */
  inline void setBetaZero(unsigned iPoint) {
    this->aab(iPoint) = 0.0;
    this->abb(iPoint) = 0.0;
    this->bbb(iPoint) = 0.0;
  }
};

/**
 * @class   d3F_dRho3 PartialDerivatives.h
 * @brief   Object to store the third derivative of the enhancement factor F w.r.t. the density, i.e.
 *          \f$ \frac{\partial^3 F}{\partial \rho^3}\f$. In the unrestricted case, this object holds
 *          four elements, \f$ \frac{\partial^3 F}{\partial \rho_\alpha \rho_\alpha \rho_\alpha}\f$,
 *          \f$ \frac{\partial^3 F}{\partial \rho_\alpha \partial\rho_\alpha \partial\rho_\beta}\f$,
 *          \f$ \frac{\partial^3 F}{\partial \rho_\alpha \partial\rho_\beta \partial\rho_\beta}\f$, and
 *          \f$ \frac{\partial^3 F}{\partial \rho_\beta \partial\rho_\beta \partial\rho_\beta}\f$.
 *          Every distinct derivative is stored as GridData.
 *
 */
template<Options::SCF_MODES SCFMode>
class d3F_dRho3;

template<>
class d3F_dRho3<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
 public:
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<>
class d3F_dRho3<Options::SCF_MODES::UNRESTRICTED> : public TriplySpinPolarizedGridDataThree {
 public:
  using TriplySpinPolarizedGridDataThree::TriplySpinPolarizedGridDataThree;
};

/**
 * @class   d3F_dRho2dGradRho PartialDerivatives.h
 * @brief   Object to store the third derivative of the enhancement factor F: twice w.r.t. the density and once w.r.t.
 * the gradient of the density, i.e. \f$ \frac{\partial^3 F}{\partial \rho^2 \partial\nabla\rho}\f$. Contains three
 * spatial variables. In the unrestricted case, this object holds six elements (see TriplySpinPolarizedGridDataTwoOne).
 * The first two (spin) indices refer to rho, the third refers to GradRho. Example: .x.aba denotes \f$ \frac{\partial^3
 * F}{\partial \rho_\alpha \partial \rho_\beta \partial \mathrm d_x \rho_\alpha}\f$.
 */
template<Options::SCF_MODES SCFMode>
class d3F_dRho2dGradRho;

template<>
class d3F_dRho2dGradRho<Options::SCF_MODES::RESTRICTED> {
 public:
  d3F_dRho2dGradRho(std::shared_ptr<GridController> gridController)
    : x(gridController), y(gridController), z(gridController) {
  }

  GridData<RESTRICTED> x;
  GridData<RESTRICTED> y;
  GridData<RESTRICTED> z;

  /**
   * @brief += operator for other d3F_dRho2dGradRho, adds them componentwise (no mixing of restricted/unrestricted
   * possible).
   */
  inline void operator+=(const d3F_dRho2dGradRho<Options::SCF_MODES::RESTRICTED>& rhs) {
    this->x += rhs.x;
    this->y += rhs.y;
    this->z += rhs.z;
  }
  /**
   * @brief Sets all elements to zero at iPoint.
   */
  inline void setZero(unsigned iPoint) {
    this->x(iPoint) = 0.0;
    this->y(iPoint) = 0.0;
    this->z(iPoint) = 0.0;
  }
};

template<>
class d3F_dRho2dGradRho<Options::SCF_MODES::UNRESTRICTED> {
 public:
  d3F_dRho2dGradRho(std::shared_ptr<GridController> gridController)
    : x(gridController), y(gridController), z(gridController) {
  }

  TriplySpinPolarizedGridDataTwoOne x;
  TriplySpinPolarizedGridDataTwoOne y;
  TriplySpinPolarizedGridDataTwoOne z;

  /**
   * @brief += operator for other d3F_dRho2dGradRho, adds them componentwise (no mixing of restricted/unrestricted
   * possible).
   */
  inline void operator+=(const d3F_dRho2dGradRho<Options::SCF_MODES::UNRESTRICTED>& rhs) {
    this->x += rhs.x;
    this->y += rhs.y;
    this->z += rhs.z;
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the alpha electron density or its gradient.
   */
  inline void setAlphaZero(unsigned iPoint) {
    this->x.setAlphaZero(iPoint);
    this->y.setAlphaZero(iPoint);
    this->z.setAlphaZero(iPoint);
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the beta electron density or its gradient.
   */
  inline void setBetaZero(unsigned iPoint) {
    this->x.setBetaZero(iPoint);
    this->y.setBetaZero(iPoint);
    this->z.setBetaZero(iPoint);
  }
};

/**
 * @class   d3F_dRhodGradRho2 PartialDerivatives.h
 * @brief   Object to store the third derivative of the enhancement factor F: once w.r.t. the density and twice w.r.t.
 * the gradient of the density, i.e. \f$ \frac{\partial^3 F}{\partial \rho \partial\left(\nabla\rho \right)^2}\f$.
 * Contains six spatial variables. In the unrestricted case, this object holds six (see
 * TriplySpinPolarizedGridDataOneTwo, e.g. xx) or eight elements (see TriplySpinPolarizedGridDataThree, e.g. xy). The
 * first (spin) index refers to rho, the second and third refer to GradRho. Example: .xy.aba denotes \f$
 * \frac{\partial^3 F}{\partial \rho_\alpha \partial  \mathrm{d}_x \rho_\beta \partial \mathrm{d}_y \rho_\alpha}\f$.
 */
template<Options::SCF_MODES SCFMode>
class d3F_dRhodGradRho2;

template<>
class d3F_dRhodGradRho2<Options::SCF_MODES::RESTRICTED> {
 public:
  d3F_dRhodGradRho2(std::shared_ptr<GridController> gridController)
    : xx(gridController), xy(gridController), xz(gridController), yy(gridController), yz(gridController), zz(gridController) {
  }

  GridData<RESTRICTED> xx;
  GridData<RESTRICTED> xy;
  GridData<RESTRICTED> xz;
  GridData<RESTRICTED> yy;
  GridData<RESTRICTED> yz;
  GridData<RESTRICTED> zz;

  /**
   * @brief += operator for other d3F_dRhodGradRho2, adds them componentwise (no mixing of restricted/unrestricted
   * possible).
   */
  inline void operator+=(const d3F_dRhodGradRho2<Options::SCF_MODES::RESTRICTED>& rhs) {
    this->xx += rhs.xx;
    this->xy += rhs.xy;
    this->xz += rhs.xz;
    this->yy += rhs.yy;
    this->yz += rhs.yz;
    this->zz += rhs.zz;
  }
  /**
   * @brief Sets all elements to zero at iPoint.
   */
  inline void setZero(unsigned iPoint) {
    this->xx(iPoint) = 0.0;
    this->xy(iPoint) = 0.0;
    this->xz(iPoint) = 0.0;
    this->yy(iPoint) = 0.0;
    this->yz(iPoint) = 0.0;
    this->zz(iPoint) = 0.0;
  }
};

template<>
class d3F_dRhodGradRho2<Options::SCF_MODES::UNRESTRICTED> {
 public:
  d3F_dRhodGradRho2(std::shared_ptr<GridController> gridController)
    : xx(gridController), xy(gridController), xz(gridController), yy(gridController), yz(gridController), zz(gridController) {
  }

  TriplySpinPolarizedGridDataOneTwo xx;
  TriplySpinPolarizedGridData xy;
  TriplySpinPolarizedGridData xz;
  TriplySpinPolarizedGridDataOneTwo yy;
  TriplySpinPolarizedGridData yz;
  TriplySpinPolarizedGridDataOneTwo zz;

  /**
   * @brief += operator for other d3F_dRhodGradRho2, adds them componentwise (no mixing of restricted/unrestricted
   * possible).
   */
  inline void operator+=(const d3F_dRhodGradRho2<Options::SCF_MODES::UNRESTRICTED>& rhs) {
    this->xx += rhs.xx;
    this->xy += rhs.xy;
    this->xz += rhs.xz;
    this->yy += rhs.yy;
    this->yz += rhs.yz;
    this->zz += rhs.zz;
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the alpha electron density or its gradient.
   */
  inline void setAlphaZero(unsigned iPoint) {
    this->xx.setAlphaZero(iPoint);
    this->xy.setAlphaZero(iPoint);
    this->xz.setAlphaZero(iPoint);
    this->yy.setAlphaZero(iPoint);
    this->yz.setAlphaZero(iPoint);
    this->zz.setAlphaZero(iPoint);
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the beta electron density or its gradient.
   */
  inline void setBetaZero(unsigned iPoint) {
    this->xx.setBetaZero(iPoint);
    this->xy.setBetaZero(iPoint);
    this->xz.setBetaZero(iPoint);
    this->yy.setBetaZero(iPoint);
    this->yz.setBetaZero(iPoint);
    this->zz.setBetaZero(iPoint);
  }
};

/**
 * @class   d3F_dGradRho3 PartialDerivatives.h
 * @brief   Object to store the third derivative of the enhancement factor F w.r.t. the gradient of the density, i.e.
 *          \f$ \frac{\partial^3 F}{\partial\left(\nabla\rho \right)^3}\f$. Contains ten spatial variables.
 *          In the unrestricted case, each of them holds
 *          four, six or eight elements, depending on the combination of spatial coordinates.
 *          Example: .xyz.aba denotes \f$ \frac{\partial^3 F}{\partial \mathrm{d}_x \rho_\alpha \partial  \mathrm{d}_y
 * \rho_\beta \partial \mathrm{d}_z \rho_\alpha}\f$. But: .xyy.aba does not exist, it needs to be .xyy.aab (since xyy is
 * a TriplySpinPolarizedGridDataOneTwo)
 */
template<Options::SCF_MODES SCFMode>
class d3F_dGradRho3;

template<>
class d3F_dGradRho3<Options::SCF_MODES::RESTRICTED> {
 public:
  d3F_dGradRho3(std::shared_ptr<GridController> gridController)
    : xxx(gridController),
      xxy(gridController),
      xxz(gridController),
      xyy(gridController),
      xyz(gridController),
      xzz(gridController),
      yyy(gridController),
      yyz(gridController),
      yzz(gridController),
      zzz(gridController) {
  }
  GridData<RESTRICTED> xxx;
  GridData<RESTRICTED> xxy;
  GridData<RESTRICTED> xxz;
  GridData<RESTRICTED> xyy;
  GridData<RESTRICTED> xyz;
  GridData<RESTRICTED> xzz;
  GridData<RESTRICTED> yyy;
  GridData<RESTRICTED> yyz;
  GridData<RESTRICTED> yzz;
  GridData<RESTRICTED> zzz;

  /**
   * @brief += operator for other d3F_dGradRho3, adds them componentwise (no mixing of restricted/unrestricted possible).
   */
  inline void operator+=(const d3F_dGradRho3<Options::SCF_MODES::RESTRICTED>& rhs) {
    this->xxx += rhs.xxx;
    this->xxy += rhs.xxy;
    this->xxz += rhs.xxz;
    this->xyy += rhs.xyy;
    this->xyz += rhs.xyz;
    this->xzz += rhs.xzz;
    this->yyy += rhs.yyy;
    this->yyz += rhs.yyz;
    this->yzz += rhs.yzz;
    this->zzz += rhs.zzz;
  }
  /**
   * @brief Sets all elements to zero at iPoint.
   */
  inline void setZero(unsigned iPoint) {
    this->xxx(iPoint) = 0.0;
    this->xxy(iPoint) = 0.0;
    this->xxz(iPoint) = 0.0;
    this->xyy(iPoint) = 0.0;
    this->xyz(iPoint) = 0.0;
    this->xzz(iPoint) = 0.0;
    this->yyy(iPoint) = 0.0;
    this->yyz(iPoint) = 0.0;
    this->yzz(iPoint) = 0.0;
    this->zzz(iPoint) = 0.0;
  }
};

template<>
class d3F_dGradRho3<Options::SCF_MODES::UNRESTRICTED> {
 public:
  d3F_dGradRho3(std::shared_ptr<GridController> gridController)
    : xxx(gridController),
      xxy(gridController),
      xxz(gridController),
      xyy(gridController),
      xyz(gridController),
      xzz(gridController),
      yyy(gridController),
      yyz(gridController),
      yzz(gridController),
      zzz(gridController) {
  }
  TriplySpinPolarizedGridDataThree xxx;  // 4
  TriplySpinPolarizedGridDataTwoOne xxy; // 6
  TriplySpinPolarizedGridDataTwoOne xxz; // 6
  TriplySpinPolarizedGridDataOneTwo xyy; // 6
  TriplySpinPolarizedGridData xyz;       // 8
  TriplySpinPolarizedGridDataOneTwo xzz; // 6
  TriplySpinPolarizedGridDataThree yyy;  // 4
  TriplySpinPolarizedGridDataTwoOne yyz; // 6
  TriplySpinPolarizedGridDataOneTwo yzz; // 6
  TriplySpinPolarizedGridDataThree zzz;  // 4

  /**
   * @brief += operator for other d3F_dGradRho3, adds them componentwise (no mixing of restricted/unrestricted possible).
   */
  inline void operator+=(const d3F_dGradRho3<Options::SCF_MODES::UNRESTRICTED>& rhs) {
    this->xxx += rhs.xxx;
    this->xxy += rhs.xxy;
    this->xxz += rhs.xxz;
    this->xyy += rhs.xyy;
    this->xyz += rhs.xyz;
    this->xzz += rhs.xzz;
    this->yyy += rhs.yyy;
    this->yyz += rhs.yyz;
    this->yzz += rhs.yzz;
    this->zzz += rhs.zzz;
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the alpha electron density or its gradient.
   */
  inline void setAlphaZero(unsigned iPoint) {
    this->xxx.setAlphaZero(iPoint);
    this->xxy.setAlphaZero(iPoint);
    this->xxz.setAlphaZero(iPoint);
    this->xyy.setAlphaZero(iPoint);
    this->xyz.setAlphaZero(iPoint);
    this->xzz.setAlphaZero(iPoint);
    this->yyy.setAlphaZero(iPoint);
    this->yyz.setAlphaZero(iPoint);
    this->yzz.setAlphaZero(iPoint);
    this->zzz.setAlphaZero(iPoint);
  }
  /**
   * @brief Sets all elements to zero at iPoint that depend on the beta electron density or its gradient.
   */
  inline void setBetaZero(unsigned iPoint) {
    this->xxx.setBetaZero(iPoint);
    this->xxy.setBetaZero(iPoint);
    this->xxz.setBetaZero(iPoint);
    this->xyy.setBetaZero(iPoint);
    this->xyz.setBetaZero(iPoint);
    this->xzz.setBetaZero(iPoint);
    this->yyy.setBetaZero(iPoint);
    this->yyz.setBetaZero(iPoint);
    this->yzz.setBetaZero(iPoint);
    this->zzz.setBetaZero(iPoint);
  }
};

/**
 * @brief Derivatives of the exchange--correlation functional can be formulated in different ways. Expressions involving
 * gradient invariants (\f$ \sigma_{\alpha\beta} = \nabla\rho_\alpha \nabla\rho_\beta \f$) can be written in terms of
 * gradients by applying the chain rule. The POTENTIAL type can be used when only the evaluation of the potential (i.e.
 * the functional derivative of the xc functional) is required. In the LDA case, the potential is simply the derivative
 * of the enhancement factor w.r.t. the density. For GGAs, however, the potential also contains a term depending on the
 * density gradient. Libxc and XCfun both offer the possibility of grouping these terms together to simplify the xc
 * potential evaluation.
 */
enum class FUNCTIONAL_DATA_TYPE { GRADIENT_INVARIANTS = 0, GRADIENTS = 1, POTENTIAL = 2 };

/**
 * @class FunctionalData PartialDerivatives.h
 * @brief A class to group the various derivatives of the xc functional. The exchange--correlation energy is written as
 * an integral over all space of an enhancement factor f, which depends on the density in the case of LDAs and on the
 * density and its gradients in the case of GGAs.
 */
template<Options::SCF_MODES SCFMode>
class FunctionalData {
 public:
  FunctionalData(const unsigned int order, const FUNCTIONAL_DATA_TYPE type, Functional functional,
                 const std::shared_ptr<GridController> gridController)
    : _order(order),
      _type(type),
      _functional(functional),
      _gridController(gridController),
      energy(0.0),
      epuv(nullptr),
      potential(nullptr),
      dFdRho(nullptr),
      d2FdRho2(nullptr),
      d3FdRho3(nullptr),
      // Gradient invariants
      dFdSigma(nullptr),
      d2FdSigma2(nullptr),
      d2FdRhodSigma(nullptr),
      d3FdRho2dSigma(nullptr),
      d3FdRhodSigma2(nullptr),
      d3FdSigma3(nullptr),
      // Gradients
      dFdGradRho(nullptr),
      d2FdRhodGradRho(nullptr),
      d2FdGradRho2(nullptr),
      d3FdRho2dGradRho(nullptr),
      d3FdRhodGradRho2(nullptr),
      d3FdGradRho3(nullptr) {
    if (order > 3)
      throw SerenityError("Partial derivatives only implemented up to third order!");
    // Initialize objects
    epuv = std::make_shared<GridData<RESTRICTED>>(_gridController);
    if (_type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
      potential = std::make_shared<GridPotential<SCFMode>>(_gridController);
    }
    else {
      if (_order >= 1) {
        dFdRho = std::make_shared<dF_dRho<SCFMode>>(_gridController);
      }
      if (_order >= 2) {
        d2FdRho2 = std::make_shared<d2F_dRho2<SCFMode>>(_gridController);
      }
      if (_order >= 3) {
        d3FdRho3 = std::make_shared<d3F_dRho3<SCFMode>>(_gridController);
      }
      if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA && _order > 0) {
        switch (_type) {
          case FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS:
            if (_order >= 1) {
              dFdSigma = std::make_shared<dF_dSigma<SCFMode>>(_gridController);
            }
            if (_order >= 2) {
              d2FdSigma2 = std::make_shared<d2F_dSigma2<SCFMode>>(_gridController);
              d2FdRhodSigma = std::make_shared<d2F_dRhodSigma<SCFMode>>(_gridController);
            }
            if (_order >= 3) {
              d3FdRho2dSigma = std::make_shared<d3F_dRho2dSigma<SCFMode>>(_gridController);
              d3FdRhodSigma2 = std::make_shared<d3F_dRhodSigma2<SCFMode>>(_gridController);
              d3FdSigma3 = std::make_shared<d3F_dSigma3<SCFMode>>(_gridController);
            }
            break;
          case FUNCTIONAL_DATA_TYPE::GRADIENTS:
            if (_order >= 1) {
              dFdGradRho = makeGradientPtr<GridPotential<SCFMode>>(_gridController);
            }
            if (_order >= 2) {
              d2FdRhodGradRho = makeGradientPtr<DoublySpinPolarizedData<SCFMode, GridData<RESTRICTED>>>(_gridController);
              d2FdGradRho2 = makeHessianPtr<DoublySpinPolarizedData<SCFMode, GridData<RESTRICTED>>>(_gridController);
            }
            if (_order >= 3) {
              d3FdRho2dGradRho = std::make_shared<d3F_dRho2dGradRho<SCFMode>>(_gridController);
              d3FdRhodGradRho2 = std::make_shared<d3F_dRhodGradRho2<SCFMode>>(_gridController);
              d3FdGradRho3 = std::make_shared<d3F_dGradRho3<SCFMode>>(_gridController);
            }
            break;
          default:
            // case POTENTIAL can not appear here since it was already checked above
            throw SerenityError("PartialDerivatives: This code should never be reached");
            break;
        }
      }
    }
  };

  /**
   * @brief Getter for the GridController associated with this FunctionalData.
   */
  const std::shared_ptr<GridController> getGridController() {
    return _gridController;
  }

  /**
   * @brief Getter for the Functional associated with this FunctionalData.
   */
  Functional getFunctional() {
    return _functional;
  }

  /**
   * @brief Getter for the derivative order associated with this FunctionalData.
   */
  unsigned int getOrder() {
    return _order;
  }

  /**
   * @brief Getter for the FUNCTIONAL_DATA_TYPE (i.e. potential, gradients or gradient invariants) associated with this
   * FunctionalData.
   */
  FUNCTIONAL_DATA_TYPE getType() {
    return _type;
  }

 private:
  const unsigned int _order;
  const FUNCTIONAL_DATA_TYPE _type;
  const Functional _functional;
  const std::shared_ptr<GridController> _gridController;

 public:
  double energy;
  /**
   * @brief The energy density  \f$ \epsilon(r) \f$
   * defined as
   *  \f$ E_{xc}=\int \epsilon (r) \mathrm{d}r \f$.\n
   *  NOT as \f$ E_{xc}=\int \epsilon (r) \rho(r)\mathrm{d}r \f$ !
   */
  std::shared_ptr<GridData<RESTRICTED>> epuv;
  std::shared_ptr<GridPotential<SCFMode>> potential;
  std::shared_ptr<dF_dRho<SCFMode>> dFdRho;
  std::shared_ptr<d2F_dRho2<SCFMode>> d2FdRho2;
  std::shared_ptr<d3F_dRho3<SCFMode>> d3FdRho3;
  std::shared_ptr<dF_dSigma<SCFMode>> dFdSigma;
  std::shared_ptr<d2F_dSigma2<SCFMode>> d2FdSigma2;
  std::shared_ptr<d2F_dRhodSigma<SCFMode>> d2FdRhodSigma;
  std::shared_ptr<d3F_dRho2dSigma<SCFMode>> d3FdRho2dSigma;
  std::shared_ptr<d3F_dRhodSigma2<SCFMode>> d3FdRhodSigma2;
  std::shared_ptr<d3F_dSigma3<SCFMode>> d3FdSigma3;
  std::shared_ptr<Gradient<GridPotential<SCFMode>>> dFdGradRho;
  // First index of DoublySpinPolarizedData corresponds to spin index of Rho.
  // Second index of DoublySpinPolarizedData corresponds to spin index of GradRho.
  std::shared_ptr<Gradient<DoublySpinPolarizedData<SCFMode, GridData<RESTRICTED>>>> d2FdRhodGradRho;
  std::shared_ptr<Hessian<DoublySpinPolarizedData<SCFMode, GridData<RESTRICTED>>>> d2FdGradRho2;
  std::shared_ptr<d3F_dRho2dGradRho<SCFMode>> d3FdRho2dGradRho;
  std::shared_ptr<d3F_dRhodGradRho2<SCFMode>> d3FdRhodGradRho2;
  std::shared_ptr<d3F_dGradRho3<SCFMode>> d3FdGradRho3;
};

} /* namespace Serenity */

#endif /* DFT_FUNCTIONALS_WRAPPERS_PARTIALDERIVATIVES_H_ */
