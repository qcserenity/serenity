/**
 * @file ROHF.h
 *
 * @date Jun 21, 2023
 * @author Niklas Niemeyer
 *
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

#ifndef SCF_ROHF
#define SCF_ROHF

/* Include Serenity Internal Headers */
#include "data/matrices/FockMatrix.h"
#include "settings/ElectronicStructureOptions.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class ElectronicStructure;

/**
 * @class ROHF
 *
 * Implements two flavors of ROHF as constrained variants of UHF.
 *
 * CUHF: J. Chem. Phys. 133, 141102 (2010)
 * SUHF: Chem. Phys. Lett. 183, 423 (1991)
 */
template<Options::SCF_MODES SCFMode>
class ROHF {
 public:
  /**
   * @brief Adds the part to the UHF Fock matrices (alpha and beta) representing the constraint.
   * @param F The Fock matrix. Modified in place.
   * @param es The electronic structure (needed for the density matrix).
   * @param rohf The ROHF type (CUHF or SUHF with corresponding lambda scaling parameter).
   * @param suhfLambda Lambda scaling parameter in case of the SUHF method.
   */
  static void addConstraint(FockMatrix<SCFMode>& F, std::shared_ptr<ElectronicStructure<SCFMode>> es,
                            Options::ROHF_TYPES rohf, double suhfLambda);
};

} /* namespace Serenity */
#endif /* SCF_ROHF */
