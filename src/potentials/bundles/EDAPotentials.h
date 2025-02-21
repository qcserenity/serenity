/**
 * @file   EDAPotentials.h
 *
 * @date   Aug 18, 2017
 * @author Moritz Bensberg, Jan Unsleber
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

#ifndef POTENTIALS_EDAPOTENTIALS_H_
#define POTENTIALS_EDAPOTENTIALS_H_

/* Include Serenity Internal Headers */
#include "potentials/bundles/PotentialBundle.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/* forward declarations */
class SystemController;
template<Options::SCF_MODES SCFMode>
class HCorePotential;
template<Options::SCF_MODES SCFMode>
class ERIPotential;

/**
 * @class EDAEnergyContributions EDAPotentials.h
 *
 * @brief The different energy contributions used in EDA calculations.
 *
 * For a detailed describtion of the energy contributions see the original puplication
 * by Kitaura and Morokuma:\n
 * K. Kitaura and K. Morokuma, Int. J. Quantum Chem. 10, 325 (1976)
 */
enum class EDAEnergyContributions { ES, ESX, ESXCT, ESXPLX, ESPL, ESXEX };
/**
 * @class EDAPotentials EDAPotentials.h
 * @brief A class containing only a single potential. Currently used for EDA calculations.
 *
 * In EDA calculations the fock which is used in the calculation of the combined system (A+B)
 * is manipulated. This is done in the MO-(Monomer) basis of the fock matrix by setting specific blocks
 * of the matrix to 0. Then the SCF calculation is performed with the constrained fock matrix.
 *
 * In this class the fock matrix is constrained according to the block selection. Thus, the transformation matrices for
 * the transformation between MO-(Monomer) and AO basis are prepared in the constructor. In the actual fock matrix
 * calculation, the unconstrained fock matrix is transformed into the MO-(Monomer) basis, the blocks are selected and
 * than the matrix is transformed back into the AO basis. The actual SCF calculation is performed in the EDATask
 * (tasks/EDATask.h/.cpp)
 *
 * Ref: Kitaura-Morokuma : K. Kitaura and K. Morokuma, Int. J. Quantum Chem. 10, 325 (1976)
 */
template<Options::SCF_MODES SCFMode>
class EDAPotentials : public PotentialBundle<SCFMode> {
 public:
  /**
   * @brief The constructor.
   * @param systemA The first system.
   * @param systemB The second system.
   * @param super The supersystem (expected to be unoptimized).
   */
  EDAPotentials(std::shared_ptr<SystemController> systemA, std::shared_ptr<SystemController> systemB,
                std::shared_ptr<SystemController> super, EDAEnergyContributions contribution);
  /**
   * @brief Default destructor.
   */
  virtual ~EDAPotentials() = default;

  /**
   * @brief A function to get the entire fock matrix.
   * @param P The density matrix.
   * @param energies The controller to add all the energy contributions to.
   * @return The current fock matrix (rebuilt on every call).
   *         (Note that this does not imply that each Potential is rebuilt, they
   *          may very well be cached. But all potentials contributing are added
   *          together in every call)
   */
  FockMatrix<SCFMode> getFockMatrix(const DensityMatrix<SCFMode>& P,
                                    std::shared_ptr<EnergyComponentController> energies) override final;

  /**
   * @brief Dummy function, which leads to an error.
   *
   */
  Eigen::MatrixXd getGradients() override final;

  /**
   * @brief Get the energy contribution calculated.
   * @return Returns the calculated energy contribution.
   */
  double getEDAEnergyContribution() {
    return _EDAEnergy;
  }

 private:
  /// @brief The first system.
  std::shared_ptr<SystemController> _a;
  /// @brief The second system.
  std::shared_ptr<SystemController> _b;
  /// @brief The supersystem.
  std::shared_ptr<SystemController> _super;
  /// @brief The potential.
  FockMatrix<SCFMode> _pot;
  /// @brief Transformation matrix AO->MO
  MatrixInBasis<SCFMode> _AO2MO;
  /// @brief Transformation matrix AO->MO
  MatrixInBasis<SCFMode> _MO2AO;
  /// @brief The EDA energy contribution.
  double _EDAEnergy = 0.0;

  /**
   * @brief Sets appropriate blocks of the input matrix to zero.
   *
   * @param matrix the matrix in which blocks are zeroed (for each spin)
   * @param x The choice in x, can be: "ESX + EX", "ESX + CT", "ESX + PLX", "ESX".
   */
  void setBlocksZero(MatrixInBasis<SCFMode>& matrix, EDAEnergyContributions x);

  ///@brief The H core potential contribtution.
  std::shared_ptr<HCorePotential<SCFMode>> _hcorepot;
  ///@brief The Hartree-Fock potential contribution.
  std::shared_ptr<ERIPotential<SCFMode>> _hfpot;
  ///@brief The EDA energy contribution calculated with this potential.
  EDAEnergyContributions _energyContr;
};

} /* namespace Serenity */

#endif /* POTENTIALS_EDAPOTENTIALS_H_ */
