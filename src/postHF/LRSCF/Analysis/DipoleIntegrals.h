/**
 * @file DipoleIntegrals.h
 *
 * @date Dec. 17, 2018
 * @author Niklas Niemeyer
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

#ifndef LRSCF_DIPOLEINTEGRALS
#define LRSCF_DIPOLEINTEGRALS

/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h"
#include "geometry/Point.h"                      //Default constructor of Point.
#include "settings/ElectronicStructureOptions.h" //RESTRICTED
/* Include Std and External Headers */
#include <Eigen/Dense> //MatrixXd
#include <memory>      //shared_ptr

namespace Serenity {
template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class DipoleIntegrals DipoleIntegrals.h
 * @brief Calculates integrals of the form <i|O|a>, where O is one of the following
 * one-particle operators: \n
 *  - electric-dipole operator in length representation   (mu = -r)\n
 *  - electric-dipole operator in velocity representation (p  = -i nabla)\n
 *  - magnetic-dipole operator in length representation   (m  = -0.5 r x p \n
 *                                                            =  0.5 i r x nabla)\n
 * For all these operators, there are three spatial components, such that three vectors
 * will be obtained. Furthermore, the fact that the latter two are imaginary is ignored.\n\n
 *
 * These integrals are stored and organized for a given linear-response problem over
 * occupied--virtual orbital pairs. They can be used to obtain\n
 *  - oscillator strengths\n
 *  - rotatory strengths\n
 *  - polarizabilities\n
 *  - optical rotation\n
 *  - ...
 */

template<Options::SCF_MODES SCFMode>
class DipoleIntegrals {
 public:
  /**
   * @brief Constructor
   * @param lrscf A vector of the LRSCF controllers of the response problem.
   * @param gaugeOrigin A real-space point that will be used as the gauge-origin. Normally
   *        chosen to be the center of mass of a system.
   */
  DipoleIntegrals(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, Point gaugeOrigin);

  /**
   * @brief Default destructor.
   */
  virtual ~DipoleIntegrals() = default;

  /**
   * @brief Returns the electric-dipole integrals in length representation.
   * One vector for each spatial component.
   * Sorted as:
   *                     int-Sys1-x   int-Sys1-y   int-Sys1-z
   *                     int-Sys1-x   int-Sys1-y   int-Sys1-z
   *                     int-Sys2-x   int-Sys2-y   int-Sys2-z
   *                     int-Sys2-x   int-Sys2-y   int-Sys2-z
   *                         ...           ...         ...
   * Integrals are given over pairs of occupied and virtual MOs.
   */
  std::shared_ptr<const Eigen::MatrixXd> getLengths();

  /**
   * @brief Returns the electric-dipole integrals in velocity representation.
   * One vector for each spatial component.
   */
  std::shared_ptr<const Eigen::MatrixXd> getVelocities();

  /**
   * @brief Returns the magnetic-dipole integrals.
   * One vector for each spatial component.
   */
  std::shared_ptr<const Eigen::MatrixXd> getMagnetics();

  /**
   * @brief Returns the gauge origin that was used to set up these integrals
   */
  Point getGaugeOrigin();

  void setFullSpace(bool fullSpace) {
    // Reset the integrals.
    if (_fullSpace != fullSpace) {
      _lengths = nullptr;
      _velocities = nullptr;
      _magnetics = nullptr;
    }
    _fullSpace = fullSpace;
  }

 private:
  ///@brief Contains all LRSCF controller.
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;

  ///@brief The gauge-origin for dipole integrals.
  Point _gaugeOrigin;

  ///@brief Compute all property integrals.
  void computeIntegrals();

  ///@brief Transforms the AO integrals to MO integrals, sorted in the response fashion (ia-wise).
  Eigen::MatrixXd ao2mo(std::vector<std::vector<MatrixInBasis<RESTRICTED>>>& ao_xyz);

  ///@brief Stores the electric length integrals computed analytically.
  std::shared_ptr<const Eigen::MatrixXd> _lengths;

  ///@brief Stores the electric velocity integrals computed analytically.
  std::shared_ptr<const Eigen::MatrixXd> _velocities;

  ///@brief Stores the magnetic integrals computed analytically.
  std::shared_ptr<const Eigen::MatrixXd> _magnetics;

  bool _fullSpace = false;
};

} /* namespace Serenity */
#endif /* LRSCF_DIPOLEINTEGRALS */
