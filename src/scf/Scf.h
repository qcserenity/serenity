/**
 * @file   Scf.h
 *
 * @date   last rework Nov 29. 2016
 * @author Thomas Dresselhaus, Jan Unsleber
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
#ifndef SCF_H_
#define SCF_H_
/* Include Serenity Internal Headers */
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>


namespace Serenity {
/* Forward Declarations */
struct Settings;
template <Options::SCF_MODES T>class PotentialBundle;
template <Options::SCF_MODES T>class ElectronicStructure;
/**
 * @class Scf Scf.h
 * @brief The algorithm of a self consistent field calculation.
 *
 * This class doesn't really hold data.
 * It just puts the different steps of the algorithm after each other.
 *
 * Purely static, cannot be instantiated.
 */
template <Options::SCF_MODES SCFMode>
class Scf {
  Scf() = delete;

public:
  /**
   * @brief Runs the algorithm of a self consistent field procedure.
   *
   * This method does not manipulate any external data directly. It just puts together the
   * different steps in the right order. For the function arguments see the corresponding classes
   * of those objects.
   *
   * @param es
   * @param potentials
   * @param maxNumberOfCycles
   * @param writeOrbitals
   */
  static void perform(
      const Settings& settings,
      std::shared_ptr<ElectronicStructure<SCFMode> > es,
      std::shared_ptr<PotentialBundle<SCFMode> > potentials,
      bool allowNotConverged = false) ;


};

} /* namespace Serenity */

#endif /* SCF_H_ */
