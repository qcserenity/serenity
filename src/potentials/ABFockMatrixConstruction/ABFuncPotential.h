/**
 * @file ABFuncPotential.h
 *
 * @date May 16, 2018
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

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABFUNCPOTENTIAL_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABFUNCPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGridController.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/DensityMatrixController.h"
#include "dft/Functional.h"
#include "notification/ObjectSensitiveClass.h"
#include "potentials/ABFockMatrixConstruction/ABPotential.h"

namespace Serenity {

/* Forward declarations */
class SystemController;

template<Options::SCF_MODES SCFMode>
class ABFuncPotential : public ABPotential<SCFMode>,
                        public ObjectSensitiveClass<Basis>,
                        public ObjectSensitiveClass<Grid>,
                        public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param activeSystem The active system controller needed for the settings.
   * @param basisA The basis controller of basis A.
   * @param basisB The basis controller of basis B.
   * @param grid The grid on which is integrated.
   * @param dMats The density matrices which contribute to the potential.
   * @param functional The exchange--correlation functional.
   */
  ABFuncPotential(std::shared_ptr<SystemController> activeSystem, std::shared_ptr<BasisController> basisA,
                  std::shared_ptr<BasisController> basisB, std::shared_ptr<GridController> grid,
                  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> dMats, Functional functional);
  /**
   * @brief Default destructor.
   */
  ~ABFuncPotential() = default;
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
  std::weak_ptr<SystemController> _actSystem;
  ///@brief The outer diagonal block of the fock matrix.
  std::unique_ptr<SPMatrix<SCFMode>> _abPotential;
  ///@brief The density matrices which contribute to the potential.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _densityMatrices;
  ///@brief The grid on which is integrated.
  std::shared_ptr<GridController> _grid;
  ///@brief The functional.
  Functional _functional;
  ///@brief The integration instance.
  std::shared_ptr<ScalarOperatorToMatrixAdder<SCFMode>> _gridToMatrix_AB;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABFUNCPOTENTIAL_H_ */
