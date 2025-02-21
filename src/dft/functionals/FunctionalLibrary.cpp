/**
 * @file   FunctionalLibrary.cpp
 *
 * @date   Sep 3, 2020
 * @author Jan Unsleber
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

/* Include Class Header*/
#include "dft/functionals/FunctionalLibrary.h"
/* Include Serenity Internal Headers */
#include "dft/functionals/wrappers/LibXC.h"
#include "dft/functionals/wrappers/XCFun.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FunctionalLibrary<SCFMode>::FunctionalLibrary(unsigned int maxBlockSize) {
#ifdef SERENITY_USE_LIBXC
  _libxc = std::make_unique<LibXC<SCFMode>>(maxBlockSize);
#endif /* SERENITY_USE_LIBXC */
#ifdef SERENITY_USE_XCFUN
  _xcfun = std::make_unique<XCFun<SCFMode>>(maxBlockSize);
#endif /* SERENITY_USE_XCFUN */
}
template<Options::SCF_MODES SCFMode>
FunctionalData<SCFMode>
FunctionalLibrary<SCFMode>::calcData(FUNCTIONAL_DATA_TYPE type, const Functional functional,
                                     const std::shared_ptr<DensityOnGridController<SCFMode>> densityOnGridController,
                                     unsigned int order) {
  if (functional.implementation() == CompositeFunctionals::IMPLEMENTATIONS::BOTH) {
    throw SerenityError("Composite functionals mixing basic functionals from LibXC & XCFun and not yet supported.");
  }
  else if (functional.implementation() == CompositeFunctionals::IMPLEMENTATIONS::XCFUN) {
#ifdef SERENITY_USE_XCFUN
    return _xcfun->calcData(type, functional, densityOnGridController, order);
#else  /* SERENITY_USE_XCFUN */
    throw SerenityError("XCFun functional requested but XCFun was not compiled.");
#endif /* SERENITY_USE_XCFUN */
  }
  else if (functional.implementation() == CompositeFunctionals::IMPLEMENTATIONS::LIBXC) {
#ifdef SERENITY_USE_LIBXC
    return _libxc->calcData(type, functional, densityOnGridController, order);
#else  /* SERENITY_USE_LIBXC */
    throw SerenityError("LibXC functional requested but LibXC was not compiled.");
#endif /* SERENITY_USE_LIBXC */
  }
  else {
#if defined SERENITY_PREFER_XCFUN && defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC
    return _xcfun->calcData(type, functional, densityOnGridController, order);
#elif defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC
    return _libxc->calcData(type, functional, densityOnGridController, order);
#elif defined SERENITY_USE_XCFUN
    return _xcfun->calcData(type, functional, densityOnGridController, order);
#else
    return _libxc->calcData(type, functional, densityOnGridController, order);
#endif
  }
}

template class FunctionalLibrary<RESTRICTED>;
template class FunctionalLibrary<UNRESTRICTED>;

} /* namespace Serenity */
