/**
 * @file Functional.cpp
 * @author Thomas Dresselhaus
 *
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
/* Include Class Header*/
#include "dft/Functional.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <stdexcept>

namespace Serenity {
Functional::Functional(CompositeFunctionals::IMPLEMENTATIONS impl,
                       const std::vector<BasicFunctionals::BASIC_FUNCTIONALS>& basicFunctionals,
                       const std::vector<double>& mixingFactors, const double hfExchangeRatio, const double hfCorrelRatio,
                       const double lrExchangeRatio, const double mu, const double ssScaling, const double osScaling)
  : _impl(impl),
    _basicFunctionals(basicFunctionals),
    _mixingFactors(mixingFactors),
    _functionalClass(CompositeFunctionals::CLASSES::NONE),
    _hfExchangeRatio(hfExchangeRatio),
    _hfCorrelRatio(hfCorrelRatio),
    _lrExchangeRatio(lrExchangeRatio),
    _mu(mu),
    _ssScaling(ssScaling),
    _osScaling(osScaling) {
  for (const auto& func : _basicFunctionals) {
    switch (BasicFunctionals::getClass[(int)func]) {
      case BasicFunctionals::CLASSES::NONE:
        // nothing to be done here, NONE is default
        break;
      case BasicFunctionals::CLASSES::LDA:
        if (_functionalClass == CompositeFunctionals::CLASSES::NONE)
          _functionalClass = CompositeFunctionals::CLASSES::LDA;
        break;
      case BasicFunctionals::CLASSES::GGA:
        if (_functionalClass == CompositeFunctionals::CLASSES::NONE || _functionalClass == CompositeFunctionals::CLASSES::LDA)
          _functionalClass = CompositeFunctionals::CLASSES::GGA;
        break;
      case BasicFunctionals::CLASSES::META_GGA:
        if (_functionalClass == CompositeFunctionals::CLASSES::NONE ||
            _functionalClass == CompositeFunctionals::CLASSES::LDA || _functionalClass == CompositeFunctionals::CLASSES::GGA)
          _functionalClass = CompositeFunctionals::CLASSES::META_GGA;
        break;
      case BasicFunctionals::CLASSES::MODELL:
        _functionalClass = CompositeFunctionals::CLASSES::MODELL;
        break;
      default:
        // Should not be reached. The BASIC_FUNCTIONALS should only work with the density and derivatives.
        throw SerenityError("A BASIC_FUNCTIONAL with class that is not LDA, GGA or META_GGA is created.");
        break;
    }
  }
}

} /* namespace Serenity */
