/**
 * @file   CDHFPotential.h
 *
 * @date   Oct 26, 2018
 * @author Lars Hellmann
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

#ifndef POTENTIALS_CDHFPOTENTIAL_H_
#define POTENTIALS_CDHFPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/HFPotential.h"
#include "settings/Options.h"

namespace Serenity {
/* Forward Declarations */
class CDIntegralController;

/**
 * @class CDHFPotential CDHFPotential.h
 *
 * @brief A class for a single System Hartree Fock potential
 * (Coulomb + exact Exchange) calculated using the Cholesky Decomposition.
 * The amount of exchange in the potential can be changed.
 *
 */
template<Options::SCF_MODES SCFMode>
class CDHFPotential : public HFPotential<SCFMode> {
 public:
  /**
   * @brief Constructor
   *
   * @param systemController Controller which provides you with all needed information
   *                         about your system and configuration.
   * @param dMAt The density matrix (controller) for this Coulomb potential.
   * @param exchangeRatio The amount of exchange added.
   * @param prescreeningThreshold Threshold parameter for integral screening (using Cauchy-Schwarz)
   * @param prescreeningIncrementStart The start integrals prescreening threshold for the incremental Fock-matrix build
   * @param prescreeningIncrementEnd   The end integrals prescreening threshold for the incremental Fock-matrix build
   */
  CDHFPotential(std::shared_ptr<SystemController> systemController,
                std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double exchangeRatio,
                const double prescreeningThreshold, double prescreeningIncrementStart, double prescreeningIncrementEnd);

  /**
   * @brief Getter for an incremental potential.
   *
   * @param F      A reference to the potential to be added to.
   * @param deltaP An increment of the density matrix.
   */
  void addToMatrix(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& densityMatrix) override final;

 private:
  // The integral Controller in the Cholesky basis
  std::shared_ptr<CDIntegralController> _cdIntController;
  // Vector  that contains all prescreening factors
  std::vector<std::vector<double>> _shellPairFactors;
  // Label for the CDStorageController to be used in this potential
  std::string _storageLabel;
};

} /* namespace Serenity */

#endif /* POTENTIALS_HFPOTENTIAL_H_ */
