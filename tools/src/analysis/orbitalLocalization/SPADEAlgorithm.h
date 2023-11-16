/**
 * @file SPADEAlgorithm.h
 *
 * @date 8 Mar 2020
 * @author Moritz Bensberg
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

#ifndef ANALYSIS_ORBITALLOCALIZATION_SPADEALGORITHM_H_
#define ANALYSIS_ORBITALLOCALIZATION_SPADEALGORITHM_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h" //Return value of run.
#include "settings/Options.h"       //SCF_MODES
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices.
#include <memory>      //smrt_ptr.

namespace Serenity {

/* Forward Declarations */
class SystemController;
/**
 * @class SPADEAlgorithm SPADEAlgorithm.h
 * @brief Performs the SPADE algorithm as documented in:\n
 *      J. Chem. Theory Comput. 2019, 15, 6085–6096.\n\n
 *
 *  The SPADE algorithm effectively localizes the supersystem orbitals such that
 *  a subset is mostly localized on the active system. The orbitals may be
 *  delocalized over the complete, respective fragment. The singular values
 *  associated to the localized orbitals are commonly used for the orbital
 *  partitioning. The algorithm suggests a number of orbitals that should be
 *  assigned to the active system.
 *
 */
template<Options::SCF_MODES SCFMode>
class SPADEAlgorithm {
 public:
  /**
   * @brief Constructor.
   * @param supersystem  The supersystem.
   * @param activeSystem The active system (has to be a fragment of the supersystem).
   */
  SPADEAlgorithm(std::shared_ptr<SystemController> supersystem, std::shared_ptr<SystemController> activeSystem);
  /**
   * @brief Run the SPADE algorithm.
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXi> run();

 private:
  ///@brief The supersystem.
  std::shared_ptr<SystemController> _supersystem;
  ///@brief The active system.
  std::shared_ptr<SystemController> _activeSystem;
};

} /* namespace Serenity */

#endif /* ANALYSIS_ORBITALLOCALIZATION_SPADEALGORITHM_H_ */
