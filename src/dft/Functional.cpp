/**
 * @file Functional.cpp
 * @author Thomas Dresselhaus
 * 
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
/* Include Class Header*/
#include "dft/Functional.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <stdexcept>


namespace Serenity {

Functional::Functional(const std::vector<BASIC_FUNCTIONALS>& basicFunctionals) :
    Functional(basicFunctionals, std::vector<double>(basicFunctionals.size(), 1.0)) {
}

Functional::Functional(
    const std::vector<BASIC_FUNCTIONALS>& basicFunctionals,
    const std::vector<double>& mixingFactors,
    const double hfExchangeRatio,
    const double hfCorrelRatio,
    const double lrExchangeRatio,
    const double mu,
    const double ssScaling,
    const double osScaling) :
      _basicFunctionals(basicFunctionals),
      _mixingFactors(mixingFactors),
      _functionalClass(FUNCTIONAL_CLASSES::LDA),
      _hfExchangeRatio(hfExchangeRatio),
      _hfCorrelRatio(hfCorrelRatio),
      _lrExchangeRatio(lrExchangeRatio),
      _mu(mu),
      _ssScaling(ssScaling),
      _osScaling(osScaling){
  for (const auto& func : _basicFunctionals) {
    switch (FunctionalClassResolver::resolveFunctionalClass(func)) {
      case FUNCTIONAL_CLASSES::NONE:
        _functionalClass = FUNCTIONAL_CLASSES::NONE;
        break;
      case FUNCTIONAL_CLASSES::LDA:
      // do nothing, this is the initial setting
        break;
      case FUNCTIONAL_CLASSES::GGA:
        if (_functionalClass == FUNCTIONAL_CLASSES::LDA) _functionalClass = FUNCTIONAL_CLASSES::GGA;
        break;
      case FUNCTIONAL_CLASSES::META_GGA:
        if (_functionalClass == FUNCTIONAL_CLASSES::LDA || _functionalClass == FUNCTIONAL_CLASSES::GGA)
          _functionalClass = FUNCTIONAL_CLASSES::META_GGA;
        break;
      case FUNCTIONAL_CLASSES::MODELL:
        _functionalClass = FUNCTIONAL_CLASSES::MODELL;
        break;
      default:
      // Should not be reached. The BASIC_FUNCTIONALS should only work with the density and derivatives.
        throw SerenityError(
                "A BASIC_FUNCTIONAL with class that is not LDA, GGA or META_GGA is created.");
        break;
    }
  }
}

} /* namespace Serenity */
