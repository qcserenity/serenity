/**
 * @file RICoulombPotential.h
 *
 * @date Nov 22, 2016
 * @author Jan Unsleber
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

#ifndef POTENTIALS_RICOULOMBPOTENTIAL_H_
#define POTENTIALS_RICOULOMBPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "misc/WarningTracker.h"
#include "potentials/IncrementalFockMatrix.h"
#include "potentials/Potential.h"

namespace Serenity {

class RI_J_IntegralController;
class SystemController;
template<Options::SCF_MODES SCFMode>
class Potential;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;

/**
 * @class RICoulombPotential RICoulombPotential.h
 *
 * A class for a single System Coulomb potential.
 * The potential is calculated using the RI approximation.
 *
 * Implementation according to:
 * [1] Weigend, F.; Kattannek, M.; Ahlrichs, R.; J.Chem.Phys. (2009), 130, 164106
 *
 * Gradients according to:
 * [2] Aquilante, F.; Lindh, R.; Pedersen, T. B.; J.Chem.Phys. (2008), 129, 034106
 *
 */
template<Options::SCF_MODES SCFMode>
class RICoulombPotential : public Potential<SCFMode>,
                           public ObjectSensitiveClass<Basis>,
                           public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor
   * @param actSystem The active system.
   * @param dMAt The density matrix (controller) for this Coulomb potential.
   * @param ri_j_IntController The controller for rij integrals.
   * @param prescreeningThreshold The Schwartz prescreening threshold.
   * @param prescreeningIncrementStart The start integrals prescreening thresold for the incremental Fock-matrix build
   * @param prescreeningIncrementEnd The end integrals prescreening thresold for the incremental Fock-matrix build
   * @param incrementSteps The number of steps of an incremental Fock-matrix build until it gets rebuild
   */
  RICoulombPotential(std::shared_ptr<SystemController> actSystem, std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                     std::shared_ptr<RI_J_IntegralController> ri_j_IntController, const double prescreeningThreshold,
                     double prescreeningIncrementStart, double prescreeningIncrementEnd, unsigned int incrementSteps);
  /// @brief Default destructor.
  virtual ~RICoulombPotential() = default;
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
  void addToMatrix(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& deltaP);

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * This potential gradient part still contains the two electron interaction potential contribution!!!
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _outOfDate = true;
  };

 private:
  ///@brief active system controller
  std::weak_ptr<SystemController> _actSystem;
  ///@brief The basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode>> _dMatController;
  ///@brief The controller for ri integrals
  const std::shared_ptr<RI_J_IntegralController> _ri_j_IntController;
  ///@brief The entire potential
  std::shared_ptr<FockMatrix<SCFMode>> _fullpotential;
  ///@brief Checks if the data is up to date
  bool _outOfDate;
  ///@brief Screening Threshold for the current iteration
  double _screening;
  ///@brief Internal iteration counter
  unsigned int _counter = 0;
  // Helper for the incremental fock matrix construction.
  std::shared_ptr<IncrementalFockMatrix<SCFMode>> _incrementHelper;
};

} /* namespace Serenity */

#endif /* POTENTIALS_COULOMBPOTENTIAL_H_ */
