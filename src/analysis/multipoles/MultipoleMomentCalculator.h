/**
 * @file MultipoleMomentCalculator.h
 *
 * @date Dec 8, 2015
 * @author David Schnieders
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

#ifndef MULTIPOLEMOMENTCALCULATOR_H_
#define MULTIPOLEMOMENTCALCULATOR_H_

/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {

namespace Options {
enum class SCF_MODES;
}

class SystemController;

/**
 * @class MultipoleMomentCalculator MultipoleMomentCalculator.h
 *
 * @brief Calculates multipole moments analytically from one electron integrals
 *
 * Equations and Definitions taken from: (28.11.2018)
 * http://crm2.univ-lorraine.fr/pages_perso/Angyan/Documents/IMF/pdf/fi_07.pdf
 *
 */

class MultipoleMomentCalculator {
 public:
  /**
   * @brief Constructor
   */
  MultipoleMomentCalculator() = default;
  /**
   * @brief Destructor
   */
  virtual ~MultipoleMomentCalculator() = default;
  /**
   * @brief Calculates the multipole moments
   *
   * @param system The system of which the quadrupole moment should be calculated
   * @param highestOrder The highest order of multipoles to be calculated
   * @return multipoleMoments vector{x,y,z} or
   *                          vector{x,y,z,xx,xy,xz,yy,yz,zz}(if requested)
   */
  template<Options::SCF_MODES SCFMode>
  static std::vector<std::vector<double>> calculateMultipoleMoment(std::shared_ptr<SystemController> system,
                                                                   unsigned int highestOrder);
};

} /* namespace Serenity */

#endif /* MULTIPOLEMOMENTCALCULATOR_H_ */
