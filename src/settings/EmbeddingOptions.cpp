/**
 * @file EmbeddingOptions.cpp
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
#include "settings/EmbeddingOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<KIN_EMBEDDING_MODES>(std::string& value, KIN_EMBEDDING_MODES& field) {
  static const std::map<std::string, KIN_EMBEDDING_MODES> m = {
      {"NONE", KIN_EMBEDDING_MODES::NONE},
      {"NADDFUNC", KIN_EMBEDDING_MODES::NADD_FUNC},
      {"NADD_FUNC", KIN_EMBEDDING_MODES::NADD_FUNC},
      {"LEVELSHIFT", KIN_EMBEDDING_MODES::LEVELSHIFT},
      {"HUZINAGA", KIN_EMBEDDING_MODES::HUZINAGA},
      {"HOFFMANN", KIN_EMBEDDING_MODES::HOFFMANN},
      {"RECONSTRUCTION", KIN_EMBEDDING_MODES::RECONSTRUCTION},
      {"FERMI_SHIFTED_HUZINAGA", KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA},
      {"FERMI", KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA},
      {"LOEWDIN", KIN_EMBEDDING_MODES::LOEWDIN},
      {"ALMO", KIN_EMBEDDING_MODES::ALMO}};
  check(m, value, field);
}
template<>
void resolve<std::vector<KIN_EMBEDDING_MODES>>(std::string& value, std::vector<KIN_EMBEDDING_MODES>& field) {
  if (value.empty()) {
  }
  else {
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      while (iss >> word) {
        KIN_EMBEDDING_MODES toResolve;
        resolve<KIN_EMBEDDING_MODES>(word, toResolve);
        field.push_back(toResolve);
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of KIN_EMBEDDING_MODES.");
    }
  }
}

template<>
void resolve<EMBEDDING_SCHEME>(std::string& value, EMBEDDING_SCHEME& field) {
  static const std::map<std::string, EMBEDDING_SCHEME> m = {{"NONE", EMBEDDING_SCHEME::NONE},
                                                            {"ISOLATED", EMBEDDING_SCHEME::ISOLATED},
                                                            {"FDE", EMBEDDING_SCHEME::FDE},
                                                            {"FAT", EMBEDDING_SCHEME::FAT}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
