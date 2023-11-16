/**
 * @file BUReconstructionPotential.h
 *
 * @date Oct 13, 2016
 * @author David Schnieders
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

#ifndef POTENTIALS_BURECONSTRUCTIONPOTENTIAL_H_
#define POTENTIALS_BURECONSTRUCTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/NAddFuncPotential.h"
#include "settings/Options.h"

namespace Serenity {

class SystemController;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;

template<Options::SCF_MODES SCFMode>
class BUReconstructionPotential : public NAddFuncPotential<SCFMode> {
 public:
  // Todo add documentation for all params.
  /**
   * @brief Constructor
   * @param actSys The active system.
   * @param activeDMat The densitymatrixcontroller for the active system.
   * @param envDMats A vector of all densitymatrixcontrollers for the environment systems.
   * @param superSystemGrid The supersystem grid.
   * @param functional The functional.
   * @param envSystems A vector of all environment systems.
   * @param smoothFactor A factor used to smoothen the potential.
   * @param potBasisLabel
   * @param singValThreshold
   * @param lbDamping
   * @param lbCycles
   * @param carterCycles
   */
  BUReconstructionPotential(std::shared_ptr<SystemController> actSys,
                            std::shared_ptr<DensityMatrixController<SCFMode>> activeDMat,
                            std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDMats,
                            std::shared_ptr<GridController> superSystemGrid, Functional functional,
                            std::vector<std::shared_ptr<SystemController>> envSystems, double smoothFactor = 0.0,
                            std::string potBasisLabel = "", const double singValThreshold = 0.0,
                            double lbDamping = 0.995, unsigned int lbCycles = 0, unsigned int carterCycles = 0);

  /**
   * @brief Getter for the potential matrix.
   * @return The potential in its matrix representation.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The density matrix for which the energy should be
   *                   calculated.
   * @return The energy associated with the potential and P.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Geometry gradient contribution from this Potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override;

  /**
   * @brief Overrides the notifying function of NAddFuncPotential
   *        to do nothing (we do not want to reperform the potential
   *        reconstructions on every SCF step).
   */
  void notify() override final{};

 private:
  /**
   * @brief Calculates the non-additive kinetic potential as described
   * in J. Chem. Phys. 133, 084103 (2010)
   */
  void calculatePotential();

  // Todo Add documentation for all member variables
  ///@brief The basis the system is defined in.
  std::shared_ptr<BasisController> _basis;
  ///@brief A vector holding all environment systems.
  std::vector<std::weak_ptr<SystemController>> _envSystems;
  ///@brief A factor used to smoothen the potential.
  double _smoothFactor;
  ///@brief
  std::string _potBasisLabel;
  ///@brief Threshold for the singular value decomposition.
  const double _singValThreshold;
  ///@brief
  double _lbDamping;
  ///@brief
  unsigned int _lbCycles;
  ///@brief
  unsigned int _carterCycles;
  ///@brief The grid representation of the potential.
  std::unique_ptr<GridPotential<SCFMode>> _potentialOnGrid;
  ///@brief The supersystem energy.
  double _supEnergy;
};

} /* namespace Serenity */

#endif /* BASICS_GRID_DATAONGRID_EXACTNADDKINPOTCALCULATOR_H_ */
