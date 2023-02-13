/**
 * @file GW_Analytic.h
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

#ifndef POSTHF_MBPT_GW_ANALYTIC_H_
#define POSTHF_MBPT_GW_ANALYTIC_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "geometry/Geometry.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/MBPT/MBPT.h"
#include "tasks/GWTask.h"
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class GW_Analytic GW_Analytic.h
 * @brief The GW_Analytic class calculates GW quasi-particle energies analytically. \n
 *        Coded in reference to: J. Chem. Theory Comput. 2013, 9, 1, 232–246 \n
 */
template<Options::SCF_MODES SCFMode>
class GW_Analytic : public MBPT<SCFMode> {
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
  GW_Analytic(std::shared_ptr<LRSCFController<SCFMode>> lrscf, GWTaskSettings settings,
              std::vector<std::shared_ptr<SystemController>> envSystemController,
              std::shared_ptr<RIIntegrals<SCFMode>> riInts = nullptr, int startOrb = 0, int endOrb = 0);
  /**
   * @brief Default destructor;
   */
  virtual ~GW_Analytic() = default;
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
   * @brief Calculates orbital energy-excitation energy difference vector
   * @param freq The frequency for which the dVector is evaluated
   * @param orbEigenValues The orbital eigenvalues
   * @param excitationEnergy The excitation energy
   * @param nOcc The number of occupied orbitals
   * @param nVirt The number of virtual orbitals
   * @return orbital energy-excitation energy difference vector
   */
  Eigen::VectorXd calculatedVector(double freq, Eigen::VectorXd orbEigenValues, double excitationEnergy,
                                   unsigned int nOcc, unsigned int nVirt);
};

} /* namespace Serenity */

#endif /* POSTHF_MBPT_GW_ANALYTIC_H_ */