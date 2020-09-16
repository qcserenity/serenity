/**
 * @file CorrelatedMethodsOptions.h
 *
 * @author Moritz Bensberg
 * @date May 12, 2020
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

#ifndef SETTINGS_CORRELATEDMETHODSOPTIONS_H_
#define SETTINGS_CORRELATEDMETHODSOPTIONS_H_

#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                      Coupled-Cluster                                           */
/**************************************************************************************************/
/**
 * Coupled-Cluster type methods.
 *  CCSD    CCSD
 *  CCSD_T  CCSD(T)
 */
enum class CC_LEVEL { CCSD, CCSD_T, DLPNO_CCSD, DLPNO_CCSD_T0 };
template<>
void resolve<CC_LEVEL>(std::string& value, CC_LEVEL& field);
/**************************************************************************************************/
/*                                           MP2                                                  */
/**************************************************************************************************/
/**
 * Various MP2 types:
 *   AO      Employ full four center integrals.
 *   RI      Use RI approximation for four center integrals.
 *   Local   Use orbital invariant formulation of MP2 using PNO, RI and pair approximation.
 */
enum class MP2_TYPES { AO = 0, RI = 1, LOCAL = 2 };
template<>
void resolve<MP2_TYPES>(std::string& value, MP2_TYPES& field);
/**************************************************************************************************/
/*                                  PNO Macro Settings                                            */
/**************************************************************************************************/
/*
 * These flags are used to adjust multiple threshold in PNO based
 * calculation at once. Note that they may have different effects
 * based on the calculation that is performed. NORMAL-PNO for DLPNO-MP2
 * means something different than NORMAL-PNO for DLPNO-CCSD.
 */
enum class PNO_SETTINGS { LOOSE, NORMAL, TIGHT };
template<>
void resolve<PNO_SETTINGS>(std::string& value, PNO_SETTINGS& field);
} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_CORRELATEDMETHODSOPTIONS_H_ */
