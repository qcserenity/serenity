/**
 * @file GW_ContourDeformation.h
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

#ifndef POSTHF_MBPT_GW_CONTOURDEFORMATION_H_
#define POSTHF_MBPT_GW_CONTOURDEFORMATION_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/MBPT/MBPT.h"
#include "tasks/GWTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/**
 * @class GW_ContourDeformation GW_ContourDeformation.h
 * @brief The GW_ContourDeformation class calculates GW quasi-particle energies via the Contour-Deformation technique.
 * \n Coded in reference to: JCTC 9 4856–4869 2018, https://publikationen.bibliothek.kit.edu/1000095752 \n
 */
template<Options::SCF_MODES SCFMode>
class GW_ContourDeformation : protected MBPT<SCFMode> {
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
  GW_ContourDeformation(std::shared_ptr<LRSCFController<SCFMode>> lrscf, GWTaskSettings settings,
                        std::vector<std::shared_ptr<SystemController>> envSystemController,
                        std::shared_ptr<RIIntegrals<SCFMode>> riInts, int startOrb = 0, int endOrb = 0);
  /**
   * @brief Default destructor;
   */
  virtual ~GW_ContourDeformation() = default;
  /**
   * @brief Calculates the GW quasi-particle energies
   * @param qpEnergy The quasi-particle energies
   * @param vxc_energies The xc orbital energies
   * @param x_energies The exchange orbital energies
   * @param correlation The correlation orbital energies
   */
  void calculateGWOrbitalenergies(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& vxc_energies,
                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& x_energies,
                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation);

  /**
   * @brief Calculates the contour residues
   * @param qpEnergy The quasi-particle energies
   * @param fermiLevel The fermi-level
   * @returns contour residue quasi-particle correlation energy
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> calculateContourResidues(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                                                       SpinPolarizedData<SCFMode, double>& fermiLevel);
  /**
   * @brief Calculates the contour residues
   * @param wnm the screened Coulomb interaction for imaginary frquencies
   * @param qpEnergy The quasi-particle energies
   * @returns imaginary frquency quasi-particle correlation energy
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> calculateImagCorr(SpinPolarizedData<SCFMode, Eigen::MatrixXd>& wnm,
                                                                SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy);

 private:
  bool double_equals(double a, double b, double epsilon = 1e-8) {
    return std::abs(a - b) < epsilon;
  }
  bool double_smaller(double a, double b, double epsilon = 1e-8) {
    if (std::abs(a - b) < epsilon)
      return false;
    else
      return a < b;
  }
};

} /* namespace Serenity */

#endif /* POSTHF_MBPT_GW_CONTOURDEFORMATION_H_ */