/**
 * @file BasisOptions.cpp
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
#include "settings/BasisOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"
/* Include Std and External Headers */
#include <map>
#include <utility>

namespace Serenity {
namespace Options {

template<>
void resolve<BASIS_PURPOSES>(std::string& value, BASIS_PURPOSES& field) {
  static const std::map<std::string, BASIS_PURPOSES> m = {
      {"DEFAULT", BASIS_PURPOSES::DEFAULT},
      {"AUX_COULOMB", BASIS_PURPOSES::AUX_COULOMB},
      {"MINBAS", BASIS_PURPOSES::MINBAS},
      {"HUECKEL", BASIS_PURPOSES::HUECKEL},
      {"IAO_LOCALIZATION", BASIS_PURPOSES::IAO_LOCALIZATION},
      {"SCF_DENS_GUESS", BASIS_PURPOSES::SCF_DENS_GUESS},
      {"AUX_CORREL", BASIS_PURPOSES::AUX_CORREL},
      {"ATOMIC_CHOLESKY", BASIS_PURPOSES::ATOMIC_CHOLESKY},
      {"ATOMIC_COMPACT_CHOLESKY", BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY},
      {"ATOMIC_CHOLESKY_ERF", BASIS_PURPOSES::ERF_ATOMIC_CHOLESKY},
      {"ATOMIC_COMPACT_CHOLESKY_ERF", BASIS_PURPOSES::ERF_ATOMIC_COMPACT_CHOLESKY},
      {"AUX_JK", BASIS_PURPOSES::AUX_JK}};
  check(m, value, field);
}

template<>
void resolve<AUX_BASIS_PURPOSES>(std::string& value, AUX_BASIS_PURPOSES& field) {
  static const std::map<std::string, AUX_BASIS_PURPOSES> m = {{"COULOMB", AUX_BASIS_PURPOSES::COULOMB},
                                                              {"EXCHANGE", AUX_BASIS_PURPOSES::EXCHANGE},
                                                              {"LREXCHANGE", AUX_BASIS_PURPOSES::LREXCHANGE},
                                                              {"CORRELATION", AUX_BASIS_PURPOSES::CORRELATION}};
  check(m, value, field);
}

template<>
void resolve<DENS_FITS>(std::string& value, DENS_FITS& field) {
  static const std::map<std::string, DENS_FITS> m = {{"RI", DENS_FITS::RI},     {"NORI", DENS_FITS::NONE},
                                                     {"NONE", DENS_FITS::NONE}, {"CD", DENS_FITS::CD},
                                                     {"ACD", DENS_FITS::ACD},   {"ACCD", DENS_FITS::ACCD}};
  check(m, value, field);
}

template<>
void resolve<EXTEND_ACD>(std::string& value, EXTEND_ACD& field) {
  static const std::map<std::string, EXTEND_ACD> m = {{"NONE", EXTEND_ACD::NONE},
                                                      {"SIMPLE", EXTEND_ACD::SIMPLE},
                                                      {"FIRST", EXTEND_ACD::FIRST},
                                                      {"COMPLETE", EXTEND_ACD::COMPLETE}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
