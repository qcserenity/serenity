/**
 * @file ABCoulombInteractionPotential.h
 *
 * @date May 8, 2018
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

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABCOULOMBINTERACTIONPOTENTIAL_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABCOULOMBINTERACTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "integrals/RI_J_IntegralController.h"
#include "notification/ObjectSensitiveClass.h"
#include "potentials/ABFockMatrixConstruction/ABPotential.h"
#include "settings/BasisOptions.h"

namespace Serenity {

/* Forward Declarations */
class SystemController;
class Libint;
/**
 * @class ABCoulombInteractionPotential ABCoulombInteractionPotential.h
 *
 * A class that calculates the Coulomb interaction contribution for an outer diagonal
 * block of the fock matrix.\n\n
 *
 * The matrix entries will have the form\n
 * \f$ J^{AB}_{\nu\mu}=\langle \chi^A_\nu| \hat{J}_C | \chi^B_\mu \rangle \f$ ,\n
 * where \f$ \chi^A_\nu \f$ and \f$ \chi^B_\mu \f$ are basis functions of the systems A
 * and B, respectively and \f$ \hat{J}_C \f$ is the Coulomb repulsion operator of system C.
 *
 */
template<Options::SCF_MODES SCFMode>
class ABCoulombInteractionPotential : public ABPotential<SCFMode>,
                                      public ObjectSensitiveClass<Basis>,
                                      public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param actSystem The system controller which's settings are used. (Could be replaced by Settings).
   * @param basisA The basis controller A.
   * @param basisB The basis controller B.
   * @param envDensityMatrixController The density matrix controllers which represent the interacting densities (system
   * C).
   * @param topDown Flag whether the calculation is a top down calculation. Reduces the cost of some evaluations.
   * @param densFitJ Flag for the density fitting which is used.
   * @param auxBasisAB The auxillary basis controller for A+B.
   * @param envAuxBasisController The auxiallary basis controller for the environment systems.
   */
  ABCoulombInteractionPotential(std::shared_ptr<SystemController> actSystem, std::shared_ptr<BasisController> basisA,
                                std::shared_ptr<BasisController> basisB,
                                std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensityMatrixController,
                                bool topDown = false, Options::DENS_FITS densFitJ = Options::DENS_FITS::NONE,
                                std::shared_ptr<BasisController> auxBasisAB = nullptr,
                                std::vector<std::shared_ptr<BasisController>> envAuxBasisController = {});

  /**
   * @brief Default destructor.
   */
  ~ABCoulombInteractionPotential() = default;
  /**
   * @brief Getter for the AB fock matrix contribution.
   * @return The AB fock matrix contribution.
   */
  SPMatrix<SCFMode>& getMatrix() override final;
  /**
   * @brief Deletes the AB fock matrix contribution if it is out of date.
   */
  void notify() override final {
    _abPotential = nullptr;
  };

 private:
  ///@brief The active system controller for the settings.
  std::weak_ptr<SystemController> _actSystem;
  ///@brief Environment densities.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _envDMatController;
  ///@brief A Libint instance.
  const std::shared_ptr<Libint> _libint;
  ///@brief The outer diagonal block of the fock matrix.
  std::unique_ptr<SPMatrix<SCFMode>> _abPotential;
  ///@brief The AB shell pair data.
  std::shared_ptr<std::vector<ShellPairData>> _shellPairsAB;
  bool _topDown;
  ///@brief The density fitting mode.
  Options::DENS_FITS _mode;
  ///@brief The auxillary basis controller A.
  std::shared_ptr<BasisController> _auxBasisControllerAB;
  ///@brief The auxillary basis controllers for the systems C.
  std::vector<std::shared_ptr<BasisController>> _envAuxBasisController;
  ///@brief The ri_j_Integral controller of BasisA+auxBasisAB.
  std::shared_ptr<RI_J_IntegralController> _ri_j_IntegralController_A_AB = nullptr;
  ///@brief The combined ABBasis.
  std::shared_ptr<BasisController> _superABBasis;
  ///@brief The RI_J_Integral controller for C == A. This has to be evaluated in most calls.
  std::shared_ptr<RI_J_IntegralController> _ri_j_IntegralController_A_B_AB = nullptr;

  /**
   * @brief Calculates the AB Fock matrix with RI-Fitting.
   */
  void calculateFockMatrixRI();
  /**
   * @brief Calculates the AB Fock matrix without any density fitting.
   */
  void calculateFockMatrixNoFitting();
};

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABCOULOMBINTERACTIONPOTENTIAL_H_ */
