/**
 * @file DeltaScf.h
 *
 * @date Oct. 11, 2022
 * @author Niklas Goellmann
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
#ifndef SCF_DELTASCF
#define SCF_DELTASCF

/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class ElectronicStructure;
template<Options::SCF_MODES SCFMode>
class SPMatrix;

template<Options::SCF_MODES SCFMode>
/**
 * @class DeltaScf DeltaScf.h
 * @brief Helps in setting up a DeltaScf calculation.
 */
class DeltaScf {
 public:
  /**
   * @brief Prepares the DeltaScf calculation.
   * @param exca Excitation of alpha electrons. {0 0} means an excitation from HOMO to LUMO,
   *             {1 1} from HOMO - 1 to Lumo + 1. {0} gives a momMatrix without excitations.
   * @param excb Excitation of beta electrons. {0 0} means an excitation from HOMO to LUMO,
   *             {1 1} from HOMO - 1 to Lumo + 1. {0} gives a momMatrix without excitations.
   * @param es Electronic structure.
   */
  static std::shared_ptr<SPMatrix<SCFMode>> prepareMOMMatrix(const std::vector<int>& exca, const std::vector<int>& excb,
                                                             std::shared_ptr<ElectronicStructure<SCFMode>> es);
};

} /* namespace Serenity */
#endif /* SCF_DELTASCF */