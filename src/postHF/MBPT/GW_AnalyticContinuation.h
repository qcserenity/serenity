/**
 * @file GW_AnalyticContinuation.h
 *
 * @date Nov 10, 2020
 * @author Johannes Toelle
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

#ifndef POSTHF_MBPT_GW_ANALYTICCONTINUATION_H_
#define POSTHF_MBPT_GW_ANALYTICCONTINUATION_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "geometry/Geometry.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/MBPT/MBPT.h"
#include "tasks/GWTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/**
 * @class GW_AnalyticContinuation GW_AnalyticContinuation.h
 * @brief The GW_AnalyticContinuation class calculates GW quasi-particle energies via the Analytic Continuation
 * technique. \n Coded in reference to: https://publikationen.bibliothek.kit.edu/1000095752,
 *                           https://authors.library.caltech.edu/104904/1/2007.03148.pdf \n
 */
template<Options::SCF_MODES SCFMode>
class GW_AnalyticContinuation : protected MBPT<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param systemController
   * @param settings The GW task settings
   * @param envSystemController
   * @param riInts The RI integrals
   * @param startOrb The first orbital index included in a GW calculation
   * @param endOrb The last orbital index included in a GW calculation
   */
  GW_AnalyticContinuation(std::shared_ptr<LRSCFController<SCFMode>> lrscf, GWTaskSettings settings,
                          std::vector<std::shared_ptr<SystemController>> envSystemController,
                          std::shared_ptr<RIIntegrals<SCFMode>> riInts, int startOrb = 0, int endOrb = 0);
  /**
   * @brief Default destructor;
   */
  virtual ~GW_AnalyticContinuation() = default;
  /**
   * @brief Calculates the GW quasi-particle energies
   * @param qpEnergy The quasi-particle energies
   * @param vxc_energies The xc orbital energies
   * @param x_energies The exchange orbital energies
   * @param correlation The correlation orbital energies
   * @param dsigma_de The derivative of the self-energy with respect to the frequency
   * @param z The linearization factor
   */
  void calculateGWOrbitalenergies(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& vxc_energies,
                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& x_energies,
                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation,
                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& dsigma_de,
                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& z);
  /**
   * @brief Performs the analytic continuation
   * @param qpEnergy The quasi-particle energies
   * @param correlation The correlation orbital energies
   * @param wnm The screened Coulomb interaction for complex frequencies
   * @param z The linearization factor
   */
  void analyticContinuation(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                            SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation,
                            SpinPolarizedData<SCFMode, Eigen::MatrixXd>& wnm,
                            SpinPolarizedData<SCFMode, Eigen::VectorXd>& dsigma_de,
                            SpinPolarizedData<SCFMode, Eigen::VectorXd>& z);

 private:
  /// @brief The imaginary frequencies for the analytic continuation
  Eigen::VectorXcd _nodesPade;
};

} /* namespace Serenity */

#endif /* POSTHF_MBPT_GW_ANALYTICCONTINUATION_H_ */
