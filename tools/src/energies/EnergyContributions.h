/**
 * @file   EnergyContributions.h
 * @author Thomas Dresselhaus
 *
 * @date März 12, 2015
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
#ifndef ENERGYCONTRIBUTIONS_H
#define ENERGYCONTRIBUTIONS_H
/* Includes */

namespace Serenity {
/* Forward declarations */
/**
 * These types of energy can appear in Serenity. They are related to each other via the
 * ENERGY_CONTRIBUTIONS_CHILDREN_MAP and ENERGY_CONTRIBUTIONS_PARENTS_MAP. For each entry these
 * relations are clearly defined. If you need a different combination of energy values, you may
 * think about extending this list (and the maps accordingly).
 */
enum class ENERGY_CONTRIBUTIONS {
  // Primitive values
  ELECTRON_KINETIC = 1,
  ELECTRON_NUCLEUS_COULOMB = 2,
  ELECTRON_ELECTRON_COULOMB = 3,
  ELECTRON_ELECTRON_CORRELATION = 4,
  NUCLEUS_NUCLEUS_COULOMB = 5,
  ONE_ELECTRON_ENERGY = 6,
  ELECTRON_ECP = 7,

  // HF Values
  HF_ELECTRON_ELECTRON_EXCHANGE = 101,
  HF_TWO_ELECTRON_ENERGY = 102,
  HF_ENERGY = 103,

  // KS-DFT Values
  KS_DFT_ELECTRON_ELECTRON_EXCHANGE = 201,
  KS_DFT_EXCHANGE_CORRELATION_NO_HF = 202,
  KS_DFT_EXACT_EXCHANGE = 203,
  KS_DFT_EXCHANGE_CORRELATION_ENERGY = 204,
  KS_DFT_TWO_ELECTRON_ENERGY = 205,
  KS_DFT_DISPERSION_CORRECTION = 206,
  KS_DFT_ENERGY = 207,
  KS_DFT_PERTURBATIVE_CORRELATION = 208,

  // Energies related to FDE
  FDE_NAD_EXCHANGE = 301,
  FDE_NAD_CORRELATION = 302,
  FDE_NAD_XC = 303,
  FDE_NAD_EXACT_EXCHANGE = 304,
  FDE_NAD_KINETIC = 305,
  FDE_NAD_DISP = 306,
  FDE_NUCLEI_ENV_ELECTRONS_COULOMB = 307,
  FDE_NUCLEI_ENV_NUCLEI_COULOMB = 308,
  FDE_ELECTRONS_ENV_ELECTRONS_COULOMB = 309,
  FDE_ELECTRONS_ENV_NUCLEI_COULOMB = 310,
  FDE_ELECTROSTATIC_INTERACTIONS = 311,
  FDE_INTERACTION_ENERGY = 312,
  FDE_EMBEDDED_KS_DFT_ENERGY = 313,
  FDE_EMBEDDED_HF_ENERGY = 314,
  FDE_SUPERSYSTEM_ENERGY_DFT_DFT = 315,
  FDE_FROZEN_SUBSYSTEM_ENERGIES = 316,
  FDE_SUPERSYSTEM_ENERGY_WF_DFT = 317,
  FDE_ECP_ENV_ENERGY = 318,
  FDE_MP2_INT_ENERGY = 319,

  // FDE solvation energies
  FDE_SOLV_SCALED_NAD_EXCHANGE = 320,
  FDE_SOLV_SCALED_NAD_CORRELATION = 321,
  FDE_SOLV_SCALED_NAD_XC = 322,
  FDE_SOLV_SCALED_NAD_EXACT_EXCHANGE = 323,
  FDE_SOLV_SCALED_NAD_KINETIC = 324,
  FDE_SOLV_SCALED_NAD_DISP = 325,
  FDE_SOLV_SCALED_NUCLEI_ENV_ELECTRONS_COULOMB = 326,
  FDE_SOLV_SCALED_NUCLEI_ENV_NUCLEI_COULOMB = 327,
  FDE_SOLV_SCALED_ELECTRONS_ENV_ELECTRONS_COULOMB = 328,
  FDE_SOLV_SCALED_ELECTRONS_ENV_NUCLEI_COULOMB = 329,
  FDE_SOLV_SCALED_ELECTROSTATIC_INTERACTIONS = 330,
  FDE_SOLV_SCALED_INTERACTION_ENERGY = 331,
  FDE_SOLV_SCALED_EMBEDDED_KS_DFT_ENERGY = 332,
  FDE_SOLV_SCALED_EMBEDDED_HF_ENERGY = 333,

  // Energies related to Projection-Based-Embedding (PBE)
  PBE_LINEAR_CORRECTION = 401,
  PBE_SUPERSYSTEM_ENERGY_WFT_DFT = 402,
  PBE_SUPERSYSTEM_ENERGY_DFT_DFT = 403,

  // Post HF Corrections and combined energies
  MP2_CORRECTION = 501,
  MP2_TOTAL_ENERGY = 502,
  CCSD_CORRECTION = 503,
  CCSD_TOTAL_ENERGY = 504,
  TRIPLES_CORRECTION = 505,
  CCSD_T_TOTAL_ENERGY = 506,

  // solvation energy
  SOLVATION_ENERGY = 601,
};

} /* namespace Serenity */
#endif /* ENERGYCONTRIBUTIONS_H */
