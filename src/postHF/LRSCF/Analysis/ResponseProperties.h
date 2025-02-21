/**
 * @file ResponseProperties.h
 * @author Niklas Niemeyer
 *
 * @date May 12, 2018
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

#ifndef LRSCF_RESPONSEPROPERTIES
#define LRSCF_RESPONSEPROPERTIES

/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h" //Options::SCF_MODES
#include "settings/LRSCFOptions.h"               //Options::GAUGE
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices.
#include <memory>      //shared_ptr

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class DipoleIntegrals;

template<Options::SCF_MODES SCFMode>
/**
 * @class ResponseProperties
 *
 * Prints the dynamic polarizabilities, optical rotation for the frequencies given in the input.
 * If a complex frequency was employed, also prints and their imaginary parts
 * (linear-absorption cross section and electronic circular dichroism).
 */
class ResponseProperties {
 public:
  /**
   * @brief Calculates response properties from solution vectors (X+Y,X-Y).
   * @param isNotCC2 Determines if CC2/ADC(2) is used.
   * @param dipoles All property integrals needed (electric length/velocity and magnetic).
   * @param perturbeddensities Solutionvectors sorted frequency-wise (3 Cartesian components each).
   * @param frequencies The corresponding frequencies.
   * @param rotFactor Conversion factor from optical rotation parameter to specific rotation.
   * @param damping Damping factor.
   * @param gauge The gauge of the properties.
   * @param results The results of the property calculation.
   */
  static void
  printProperties(bool isNotCC2, const std::shared_ptr<DipoleIntegrals<SCFMode>> dipoles,
                  const std::vector<Eigen::MatrixXd>& perturbeddensities, const std::vector<double>& frequencies,
                  const double rotFactor, const double damping, Options::GAUGE gauge,
                  std::vector<Eigen::Matrix3d> Fdipdip, std::vector<Eigen::Matrix3d> Fdipmag,
                  std::vector<std::tuple<double, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d>>& results);
};

} /* namespace Serenity */
#endif /* LRSCF_RESPONSEPROPERTIES */
