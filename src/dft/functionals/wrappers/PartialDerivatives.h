/**
 * @file PartialDerivatives.h
 *
 * @date Apr 4, 2017
 * @author Michael Boeckers
 * @date Jul 12, 2015
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

#ifndef DFT_FUNCTIONALS_WRAPPERS_PARTIALDERIVATIVES_H_
#define DFT_FUNCTIONALS_WRAPPERS_PARTIALDERIVATIVES_H_

/* Include Serenity Internal Headers */
#include "math/Derivatives.h"
#include "dft/Functional.h"
#include "data/grid/GridData.h"
#include "data/grid/GridPotential.h"
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
template<Options::SCF_MODES T> class dF_dRho : public GridPotential<T>{
public:
  using GridPotential<T>::GridPotential;
  explicit dF_dRho(dF_dRho<T>&) = default;
  explicit dF_dRho(const dF_dRho<T>&) = default;
  dF_dRho(dF_dRho<T> &&) = default;
  dF_dRho& operator=(const dF_dRho<T>&) = default;
  dF_dRho& operator=(dF_dRho<T> &&) = default;
  using GridPotential<T>::operator+=;
  using GridPotential<T>::operator-=;
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
template<Options::SCF_MODES T> class d2F_dRho2;

template<> class d2F_dRho2<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
public:
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<> class d2F_dRho2<Options::SCF_MODES::UNRESTRICTED>  {
public:
  /**
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  d2F_dRho2(std::shared_ptr<GridController> gridController):
    aa(gridController),
    ab(gridController),
    bb(gridController){}

  GridData<RESTRICTED> aa;
  GridData<RESTRICTED> ab;
  GridData<RESTRICTED> bb;
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
template<Options::SCF_MODES T> class d3F_dRho3;

template<> class d3F_dRho3<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
public:
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<> class d3F_dRho3<Options::SCF_MODES::UNRESTRICTED> {
public:
  /**
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  d3F_dRho3(std::shared_ptr<GridController> gridController):
    aaa(gridController),
    aab(gridController),
    abb(gridController),
    bbb(gridController){}

  GridData<RESTRICTED> aaa;
  GridData<RESTRICTED> aab;
  GridData<RESTRICTED> abb;
  GridData<RESTRICTED> bbb;
};

/**
 * @class   dF_dSigma PartialDerivatives.h
 * @brief   Object to store the derivative of the enhancement factor F w.r.t. \f$ \sigma \f$, i.e.
 *          \f$ \frac{\partial F}{\partial \sigma}\f$ where \f$ \sigma =  (\nabla \rho)^2 \f$. In the
 *          unrestricted case this objects holds three elements, i.e. \f$ \frac{\partial F}{\partial \sigma_{\alpha\alpha}}\f$,
 *          \f$ \frac{\partial F}{\partial \sigma_{\alpha\beta}}\f$ and
 *          \f$ \frac{\partial F}{\partial \sigma_{\beta\beta}}\f$ .
 *          Every distinct derivative is stored as GridData.
 *
 */
template<Options::SCF_MODES T> class dF_dSigma;

template<> class dF_dSigma<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<> class dF_dSigma<Options::SCF_MODES::UNRESTRICTED>  {
public:
  /**
   *
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  dF_dSigma(std::shared_ptr<GridController> gridController):
    aa(gridController),
    ab(gridController),
    bb(gridController){}

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
template<Options::SCF_MODES T> class d2F_dSigma2;

template<> class d2F_dSigma2<Options::SCF_MODES::RESTRICTED>: public GridData<RESTRICTED> {
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<> class d2F_dSigma2<Options::SCF_MODES::UNRESTRICTED> {
public:
  /**
   *
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  d2F_dSigma2(std::shared_ptr<GridController> gridController):
    aaaa(gridController),
    aaab(gridController),
    aabb(gridController),
    abab(gridController),
    abbb(gridController),
    bbbb(gridController){}

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
 *           i.e. \f$ \frac{\partial^2 F}{\partial \rho \partial \sigma}\f$ where \f$ \sigma =  (\nabla \rho)^2 \f$. In the
 *          unrestricted case this objects holds six elements, i.e.
 *          \f$ \frac{\partial^2 F}{\partial \rho_\alpha \partial\sigma_{\alpha\alpha}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \rho_\alpha \partial\sigma_{\alpha\beta}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \rho_\alpha \partial\sigma_{\beta\beta}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \rho_\beta \partial\sigma_{\alpha\alpha}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \rho_\beta \partial\sigma_{\alpha\beta}}\f$,
 *          \f$ \frac{\partial^2 F}{\partial \rho_\beta \partial\sigma_{\beta\beta}}\f$.
 *          Every distinct derivative is stored as GridData.
 *          Note that the first letter of the GridData object xxx within this class always refers to the derivative w.r.t.
 *          the density while the latter two determine the spin-indices of \f$ \sigma \f$.
 */
template<Options::SCF_MODES T> class d2F_dRhodSigma;

template<> class d2F_dRhodSigma<Options::SCF_MODES::RESTRICTED> : public GridData<RESTRICTED> {
  /**
   * @brief See GridData for more information.
   */
  using GridData<RESTRICTED>::GridData;
};

template<> class d2F_dRhodSigma<Options::SCF_MODES::UNRESTRICTED> {
public:
  /**
   *
   * @param gridController The gridController controlling the grid on which the derivatives are calculated.
   */
  d2F_dRhodSigma(std::shared_ptr<GridController> gridController):
    aaa(gridController),
    aab(gridController),
    abb(gridController),
    baa(gridController),
    bab(gridController),
    bbb(gridController){}

  GridData<RESTRICTED> aaa;
  GridData<RESTRICTED> aab;
  GridData<RESTRICTED> abb;
  GridData<RESTRICTED> baa;
  GridData<RESTRICTED> bab;
  GridData<RESTRICTED> bbb;
};

enum class FUNCTIONAL_DATA_TYPE {
  GRADIENT_INVARIANTS = 0,
  GRADIENTS = 1,
  POTENTIAL = 2
};

template<Options::SCF_MODES T> class FunctionalData {
public:
  FunctionalData(
      const unsigned int order,
      const FUNCTIONAL_DATA_TYPE type,
      Functional functional,
      const std::shared_ptr<GridController> gridController):
    _order(order),
    _type(type),
    _functional(functional),
    _gridController(gridController),
    energy(0.0),
    epuv(nullptr),
    potential(nullptr),
    dFdRho(nullptr),
    d2FdRho2(nullptr),
    dFdSigma(nullptr),
    d2FdSigma2(nullptr),
    d2FdRhodSigma(nullptr),
    dFdGradRho(nullptr){
    assert(order <= 2 && "Partial derivatives only implemented up to second order");
    assert(_gridController);
    //Initialize objects
    epuv = std::make_shared<GridData<RESTRICTED> >(_gridController);
    if (_type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
      potential = std::make_shared<GridPotential<T> >(_gridController);
    } else {
      if (_order >= 1) {
        dFdRho = std::make_shared<dF_dRho<T> >(_gridController);
      }
      if (_order >= 2) {
        d2FdRho2 = std::make_shared<d2F_dRho2<T> >(_gridController);
      }
      if (_functional.getFunctionalClass() == FUNCTIONAL_CLASSES::GGA && _order > 0) {
        switch(_type) {
          case FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS :
            if (_order >= 1) {
              dFdSigma = std::make_shared<dF_dSigma<T> >(_gridController);
            }
            if (_order >= 2) {
              d2FdSigma2 = std::make_shared<d2F_dSigma2<T> >(_gridController);
              d2FdRhodSigma = std::make_shared<d2F_dRhodSigma<T> >(_gridController);
            }
            break;
          case FUNCTIONAL_DATA_TYPE::GRADIENTS :
            assert (order < 2 && "GRADIENTS type only implemented up to first order");
            if (_order >= 1) {
              dFdGradRho = makeGradientPtr<GridPotential<T> >(_gridController);
            }
            break;
          case FUNCTIONAL_DATA_TYPE::POTENTIAL :
            assert(false);
            break;
        }
      }
    }
  };

  const std::shared_ptr<GridController> getGridController() {
    return _gridController;
  }

  Functional getFunctional() {
    return _functional;
  }

  unsigned int getOrder() {
    return _order;
  }

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
  std::shared_ptr<GridData<RESTRICTED> > epuv;
  std::shared_ptr<GridPotential<T> > potential;
  std::shared_ptr<dF_dRho<T> > dFdRho;
  std::shared_ptr<d2F_dRho2<T> > d2FdRho2;
  std::shared_ptr<dF_dSigma<T> > dFdSigma;
  std::shared_ptr<d2F_dSigma2<T> > d2FdSigma2;
  std::shared_ptr<d2F_dRhodSigma<T> > d2FdRhodSigma;
  std::shared_ptr<Gradient<GridPotential<T> > > dFdGradRho;
};

} /* namespace Serenity */

#endif /* DFT_FUNCTIONALS_WRAPPERS_PARTIALDERIVATIVES_H_ */
