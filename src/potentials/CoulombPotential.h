/**
 * @file CoulombPotential.h
 *
 * @date Nov 22, 2016
 * @author Jan Unsleber
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

#ifndef POTENTIALS_COULOMBPOTENTIAL_H_
#define POTENTIALS_COULOMBPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
//#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/ElectronicStructure.h"
#include "integrals/wrappers/Libint.h"
#include "settings/Options.h"
#include "data/OrbitalController.h"
#include "potentials/Potential.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "system/SystemController.h"

namespace Serenity {
class RI_J_IntegralController;

/**
 * @class CoulombPotential CoulombPotential.h
 *
 * A class for a single System Coulomb potential.
 * The potential is calculated using the RI approximation.
 *
 * Implementation according to:
 * [1] Weigend, F.; Kattannek, M.; Ahlrichs, R.; J.Chem.Phys (2009), 130, 164106
 *
 * Gradients according to:
 * [2] Aquilante, F.; Lindh, R.; Pedersen, T. B.; J.Chem.Phys. (2008), 129, 034106
 *
 */
template <Options::SCF_MODES SCFMode>
class CoulombPotential : public Potential<SCFMode>,
                         public ObjectSensitiveClass<Basis>,
                         public ObjectSensitiveClass<DensityMatrix<SCFMode> >{
public:
  /**
   * @brief Constructor
   * @param actSystem The active system.
   * @param dMAt The density matrix (controller) for this Coulomb potential.
   * @param ri_j_IntController The controller for rij integrals.
   * @param prescreeningThreshold The Schwartz prescreening threshold.
   */
  CoulombPotential(std::shared_ptr<SystemController> actSystem,
                   std::shared_ptr<DensityMatrixController<SCFMode> > dMat,
                   std::shared_ptr<RI_J_IntegralController> ri_j_IntController,
                   const double prescreeningThreshold):
    Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _actSystem(actSystem),
    _prescreeningThreshold(prescreeningThreshold),
    _dMatController(dMat),
    _ri_j_IntController(ri_j_IntController),
    _potential(nullptr){
    this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
    this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode> >::_self);
  };
  /// @brief Default destructor.
  virtual ~CoulombPotential() =default;
  /**
   * @brief Getter for the actual potential.
   *
   * This function makes use of the RI approximation.
   *
   * @return Returns the potential in matrix form.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Adds increment to an existing potential in matrix form.
   *
   * Working function, that calculates the Fock matrix.
   * This function will prescreen based on the delta-density matrix.
   * This function makes use of the RI approximation.
   *
   * @param F      A reference to the potential to be added to.
   * @param deltaP An increment of the density matrix.
   */
  void addToMatrix(FockMatrix<SCFMode>& F,
                   const DensityMatrix<SCFMode>& deltaP);

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * This potential gradient part still contains the two electron interaction potential contribution!!!
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Matches each basis function shell to its respective atom center.
   *
   * @param basisIndicesRed see AtomCenteredBasisController
   * @param nBasisFunctionRed the (reduced) number of basis functions
   */
  static std::vector<unsigned int> createBasisToAtomIndexMapping(
      const std::vector<std::pair<unsigned int, unsigned int> >& basisIndicesRed,
      unsigned int nBasisFunctionsRed);

  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final{
    _potential = nullptr;
  };


private:
  ///@brief active system controller
  std::shared_ptr<SystemController> _actSystem;
  ///@brief Threshold for the integral prescreening.
  const double _prescreeningThreshold;
  ///@brief The basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode> > _dMatController;
  ///@brief A Libint instance.
  const std::shared_ptr<Libint> _libint = Libint::getSharedPtr();
  ///@brief The controller for ri integrals
  const std::shared_ptr<RI_J_IntegralController> _ri_j_IntController;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode> >_potential;
};

} /* namespace Serenity */

#endif /* POTENTIALS_COULOMBPOTENTIAL_H_ */
