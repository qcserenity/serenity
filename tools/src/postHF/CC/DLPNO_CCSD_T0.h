/**
 * @file DLPNO_CCSD_T0.h
 *
 * @author Moritz Bensberg
 * @date Nov 5, 2019
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

#ifndef POSTHF_CC_DLPNO_CCSD_T0_H_
#define POSTHF_CC_DLPNO_CCSD_T0_H_

/* Include Std and External Headers */
#include <memory> //smrt_ptr

namespace Serenity {
/* Forward Declarations */
class LocalCorrelationController;
class OrbitalTriple;
class OrbitalTripleSet;

/**
 * @class DLPNO_CCSD_T0 DLPNO_CCSD_T0.h
 * @brief A static class that performs the calculation of the semi-canonical triples
 *        correction for DLNP-CCSD(T0).\n
 *        See J. Chem. Phys. 139, 134101 (2013) for details.\n
 *        The actual triples correction is calculated in @see data/OrbitalTriple.h
 */
class DLPNO_CCSD_T0 {
 private:
  DLPNO_CCSD_T0() = default;

 public:
  virtual ~DLPNO_CCSD_T0() = default;
  /**
   * @brief Calculates the semi-canonical triples correction for the given correlation set up.
   * @param localCorrelationController The local correlation controller.
   * @return The semi-canonical triples correction.
   */
  static double calculateEnergyCorrection(std::shared_ptr<LocalCorrelationController> localCorrelationController);
};

} /* namespace Serenity */

#endif /* POSTHF_CC_DLPNO_CCSD_T0_H_ */
