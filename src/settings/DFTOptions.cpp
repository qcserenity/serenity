/**
 * @file DFTOptions.cpp
 *
 * @author Moritz Bensberg
 * @date May 11, 2020
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
#include "settings/DFTOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<DFT_DISPERSION_CORRECTIONS>(std::string& value, DFT_DISPERSION_CORRECTIONS& field) {
  static const std::map<std::string, DFT_DISPERSION_CORRECTIONS> m = {{"NONE", DFT_DISPERSION_CORRECTIONS::NONE},
                                                                      {"D3", DFT_DISPERSION_CORRECTIONS::D3},
                                                                      {"D3ABC", DFT_DISPERSION_CORRECTIONS::D3ABC},
                                                                      {"D3BJ", DFT_DISPERSION_CORRECTIONS::D3BJ},
                                                                      {"D3BJABC", DFT_DISPERSION_CORRECTIONS::D3BJABC}};
  check(m, value, field);
}

template<>
void resolve<std::vector<CompositeFunctionals::XCFUNCTIONALS>>(std::string& value,
                                                               std::vector<CompositeFunctionals::XCFUNCTIONALS>& field) {
  if (value.empty()) {
  }
  else {
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      while (iss >> word) {
        CompositeFunctionals::XCFUNCTIONALS toResolve;
        resolve<CompositeFunctionals::XCFUNCTIONALS>(word, toResolve);
        field.push_back(toResolve);
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of XCFunctionals.");
    }
  }
}

template<>
void resolve<std::vector<CompositeFunctionals::KINFUNCTIONALS>>(std::string& value,
                                                                std::vector<CompositeFunctionals::KINFUNCTIONALS>& field) {
  if (value.empty()) {
  }
  else {
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      while (iss >> word) {
        CompositeFunctionals::KINFUNCTIONALS toResolve;
        resolve<CompositeFunctionals::KINFUNCTIONALS>(word, toResolve);
        field.push_back(toResolve);
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of KINFUNCTIONALS.");
    }
  }
}

template<>
void resolve<CompositeFunctionals::IMPLEMENTATIONS>(std::string& value, CompositeFunctionals::IMPLEMENTATIONS& field) {
  // the lowercase strings are not supposed to act as keys, but are intended for the reverse resolve call
  static const std::map<std::string, CompositeFunctionals::IMPLEMENTATIONS> m1 = {
      {"XCFUN", CompositeFunctionals::IMPLEMENTATIONS::XCFUN},
      {"LIBXC", CompositeFunctionals::IMPLEMENTATIONS::LIBXC},
      {"xcfun", CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR}};
  static const std::map<std::string, CompositeFunctionals::IMPLEMENTATIONS> m2 = {
      {"XCFUN", CompositeFunctionals::IMPLEMENTATIONS::XCFUN},
      {"LIBXC", CompositeFunctionals::IMPLEMENTATIONS::LIBXC},
      {"libxc", CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR}};
#if defined SERENITY_PREFER_XCFUN && defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC
  check(m1, value, field);
#elif defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC
  check(m2, value, field);
#elif defined SERENITY_USE_XCFUN
  check(m1, value, field);
#else
  check(m2, value, field);
#endif
}

template<>
void resolve<std::vector<BasicFunctionals::BASIC_FUNCTIONALS>>(std::string& value,
                                                               std::vector<BasicFunctionals::BASIC_FUNCTIONALS>& field) {
  if (value.empty()) {
    if (field.size()) {
      value = "{ ";
      for (BasicFunctionals::BASIC_FUNCTIONALS val : field) {
        std::string varAsString;
        resolve<BasicFunctionals::BASIC_FUNCTIONALS>(varAsString, val);
        value += (varAsString + " ");
      }
      value += "}";
    }
  }
  else {
    field.clear();
    std::istringstream iss(value);
    std::string word;
    while (iss >> word) {
      BasicFunctionals::BASIC_FUNCTIONALS toResolve;
      resolve<BasicFunctionals::BASIC_FUNCTIONALS>(word, toResolve);
      field.push_back(toResolve);
    }
  }
}

} /* namespace Options */
} /* namespace Serenity */
