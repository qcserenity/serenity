/**
 * @file NAddFuncPotential.h
 *
 * @date Nov 24, 2016
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

#ifndef POTENTIALS_NADDFUNCPOTENTIAL_H_
#define POTENTIALS_NADDFUNCPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGrid.h"
#include "data/grid/GridPotential.h"
#include "dft/Functional.h"
#include "potentials/Potential.h"
#include "settings/Options.h"

namespace Serenity {
/*forward declarations*/

template<Options::SCF_MODES SCFMode>
class ScalarOperatorToMatrixAdder;
template<Options::SCF_MODES SCFMode>
class SupersystemDensityOnGridController;
template<Options::SCF_MODES SCFMode>
class DensityOnGridController;
template<class T>
struct Gradient;
class BasisFunctionOnGridController;
class EnergyComponentController;
class SystemController;
template<Options::SCF_MODES SCFMode>
class ExchangeInteractionPotential;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;

/**
 * @class NAddEnergyHelper NAddFuncPotential.h
 * @brief A helper class to track density changes in the environment
 *          seperate of the changes in the active system
 */
template<Options::SCF_MODES SCFMode>
class NAddEnergyHelper : public ObjectSensitiveClass<DensityOnGrid<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param
   * @param
   */
  NAddEnergyHelper(Functional func, std::shared_ptr<DensityOnGridController<SCFMode>> ctrs);

  /// @brief Default destructor.
  virtual ~NAddEnergyHelper() = default;

  /**
   * @brief The main function, tracking the environment system XC energies.
   * @returns Returns the current negative sum of subsystems XC energies.
   */
  double getEnergy();

  /// @brief @see notification.
  virtual void notify() override {
    _old = true;
  };

 private:
  double _energy = 0.0;
  bool _old = true;
  ///@brief The functional.
  Functional _functional;
  ///@brief The controller for the densities on the grid.
  std::shared_ptr<DensityOnGridController<SCFMode>> _densOnGridControllers;
};

/**
 * @class NAddFuncPotential NAddFuncPotential.h
 * @brief
 */
template<Options::SCF_MODES SCFMode>
class NAddFuncPotential : public Potential<SCFMode>,
                          public ObjectSensitiveClass<Grid>,
                          public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param system The system needed for its config.
   * @param activeDMat The active density (controller) and basis this
   *                   potential is defined in.
   * @param envDMats The environment densities.
   * @param grid The grid for the intermediate calculation on the grid.
   * @param functional The functional to be used.
   * @param eCon A pair of a boolean and a list of EnergyComponentControllers.
   *             The boolean describes if a XC (true) or kinetic functional
   *             is tracked. The EnergyComponentControllers are used to store
   *             the respective XC or kinetic energies of each subsystem when
   *             they are evaluated via a functional on the supersystem grid.
   *             This allows for the reuse of these subtracted energies across
   *             multiple calculations in a freeze-and-thaw run.
   *             If this optional argument is not given, the subtracted energies
   *             of environment systems are reevaluated at the start of each
   *             SCF run. This becomes costly for more than 50 subsystems (water),
   *             or even earlier for bigger systems.
   *             The active system controller has to be at position 0 in the vector.
   */
  NAddFuncPotential(std::shared_ptr<SystemController> system, std::shared_ptr<DensityMatrixController<SCFMode>> activeDMat,
                    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDMats,
                    std::shared_ptr<GridController> grid, Functional functional,
                    std::pair<bool, std::vector<std::shared_ptr<EnergyComponentController>>> eCon =
                        {true, std::vector<std::shared_ptr<EnergyComponentController>>(0)},
                    bool evaluateEnergy = true, bool evaluateExactX = true, bool calculateSolvationEnergy = false);

  /**
   * @brief Constructor.
   * @param system The system needed for its config.
   * @param activeDMat The active density (controller) and basis this
   *                   potential is defined in.
   * @param otherExactDmats The Density matrices of other exactly treated subsystems in mixed exact/approx embedding.
   * @param BtoAProjections The projector to sort the exactly treated env density matrices
   * @param envDMats The Density matrices of approximated treated subsystems
   * @param grid The for the intermediate calculation on the grid.
   * @param functional The functional to be used.
   * @param eCon A pair of a boolean and a list of EnergyComponentControllers.
   *             The boolean describes if a XC (true) or kinetic functional
   *             is tracked. The EnergyComponentControllers are used to store
   *             the respective XC or kinetic energies of each subsystem when
   *             they are evaluated via a functional on the supersystem grid.
   *             This allows for the reuse of these subtracted energies across
   *             multiple calculations in a freeze-and-thaw run.
   *             If this optional argument is not given, the subtracted energies
   *             of environment systems are reevaluated at the start of each
   *             SCF run. This becomes costly for more than 50 subsystems (water),
   *             or even earlier for bigger systems.
   *             The active system controller has to be at position 0 in the vector.
   */
  NAddFuncPotential(std::shared_ptr<SystemController> system, std::shared_ptr<DensityMatrixController<SCFMode>> activeDMat,
                    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> otherExactDmats,
                    std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> BtoAProjections,
                    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDMats,
                    std::shared_ptr<GridController> grid, Functional functional,
                    std::pair<bool, std::vector<std::shared_ptr<EnergyComponentController>>> eCon =
                        {true, std::vector<std::shared_ptr<EnergyComponentController>>(0)},
                    bool evaluateEnergy = true, bool evaluateExactX = true, bool calculateSolvationEnergy = false);

  /// @brief Default destructor.
  ~NAddFuncPotential();

  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override;

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override;

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  virtual double getEnergy(const DensityMatrix<SCFMode>& P) override;
  /**
   * @brief Calculates the energy as scaling * (P dot F).sum().
   * @param P        The density matrix.
   * @param scaling  The scaling
   * @return The linearized and scaled energy.
   */
  double getLinearizedEnergy(const DensityMatrix<SCFMode>& P, double scaling);

  /**
   * @brief Getter.
   * @return Returns the functional.
   */
  Functional getFunctional() {
    return _functional;
  }
  /**
   * @brief Getter.
   * @return Returns the grid (controller).
   */
  std::shared_ptr<GridController> getGridController() {
    return _grid;
  }

  /**
   * @brief Getter.
   * @return Returns the grid potential derivative w.r.t. the density.
   */
  std::unique_ptr<Gradient<GridPotential<SCFMode>>> getGridPotentialDerivative();

  /**
   * @brief This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  virtual void notify() override {
    _potential = nullptr;
  };

 protected:
  ///@brief The system only needed for its config.
  std::weak_ptr<SystemController> _system;
  ///@brief The active density and basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode>> _actDMatController;
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _otherExactDmats;
  ///@brief The environment densities
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _envDMatController;
  ///@brief The grid of the supersystem.
  std::shared_ptr<GridController> _grid;
  ///@brief The functional.
  Functional _functional;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  ///@brief The controller for the densities on the supersystem grid.
  std::vector<std::shared_ptr<DensityOnGridController<SCFMode>>> _densOnGridControllers;
  ///@brief A controller for the total density (act+env) on the supersystem grid.
  std::shared_ptr<SupersystemDensityOnGridController<SCFMode>> _supersysDensOnGridController;
  ///@brief The controller of the total environment density.
  std::shared_ptr<SupersystemDensityOnGridController<SCFMode>> _environmentDensOnGridController = nullptr;
  ///@brief Conversion too from supersystem grid to matrix.
  std::shared_ptr<ScalarOperatorToMatrixAdder<SCFMode>> _gridToMatrix;
  ///@brief Basis functions on grid of the supersystem, need to be kept for gradient calculations
  std::shared_ptr<BasisFunctionOnGridController> _basisFunctionOnGridController;
  ///@brief Exact XC pot for hybrid functionals (experimental).
  std::unique_ptr<ExchangeInteractionPotential<SCFMode>> _excPot;
  ///@brief NaddFunc energy
  double _energy;
  ///@brief helper
  std::vector<std::unique_ptr<NAddEnergyHelper<SCFMode>>> _helper;
  ///@brief EnergyController to avoid recalculation of data in FaT runs
  std::vector<std::shared_ptr<EnergyComponentController>> _energyController;
  ///@brief boolean to track which energy to take from the _energyController
  bool _isXC;
  ///@brief boolean to check if energy needs to be evaluated
  bool _evaluateEnergy;
  ///@brief boolean to check if exact exchange interaction needs to be evaluated
  bool _evaluateExactX;
  bool _calculateSolvationEnergy;
};

} /* namespace Serenity */

#endif /* POTENTIALS_NADDFUNCPOTENTIAL_H_ */
