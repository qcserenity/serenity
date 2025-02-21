/**
 * @file   LevelshiftHybridPotential.h
 *
 * @date Nov 10, 2018
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

#ifndef POTENTIALS_LEVELSHIFTHYBRIDPOTENTIAL_H_
#define POTENTIALS_LEVELSHIFTHYBRIDPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "dft/Functional.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "notification/ObjectSensitiveClass.h"
#include "potentials/Potential.h"
#include "settings/DFTOptions.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

/* Forward Declarations */
class SystemController;
class Atom;
class GridController;
class Functional;
class EnergyComponentController;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
template<Options::SCF_MODES SCFMode>
class SPMatrix;
template<Options::SCF_MODES SCFMode>
class LevelshiftPotential;
template<Options::SCF_MODES SCFMode>
class NAddFuncPotential;

/**
 * @class LevelshiftHybridPotential LevelshiftHybridPotential.h
 * @brief A generalized version of the LevelshiftPotential (see potentials/LevelshiftPotential.h) which accepts an
 *        arbitrary number of environment systems and is able to use a hybrid method of level-shift and non-additive
 * kinetic energy functional.
 */
template<Options::SCF_MODES SCFMode>
class LevelshiftHybridPotential : public Potential<SCFMode>,
                                  public ObjectSensitiveClass<Basis>,
                                  public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param activeSystemController The active system.
   * @param envSystems             The environment systems.
   * @param levelShiftParameter    The level-shift parameter.
   * @param isHybrid (optional)    Flag for using a hybrid method of level-shift + NAddFunc
   * @param basisFunctionRatio     (optional) The basis function ratio of kept basis functions. Determination of distant
   * atoms. (needs isHybrid = true)
   * @param borderAtomThreshold    (optional) The population threshold of the environment orbitals on non-distant atoms.
   *                               (needs isHybrid = true)
   * @param grid                   (optional) The grid used for the non-additive functional part.
   *                               (needs isHybrid = true)
   * @param functional             (optional) The non-additive functional.
   *                               (needs isHybrid = true)
   * @param eCon                   (optional) The energy component controller for the non-additive functional.
   *                               (needs isHybrid = true)
   */
  LevelshiftHybridPotential(
      std::shared_ptr<SystemController> activeSystemController,
      std::vector<std::shared_ptr<SystemController>> envSystems, const double levelShiftParameter, bool isHybrid = false,
      double basisFunctionRatio = 0.0, double borderAtomThreshold = 0.02, std::shared_ptr<GridController> grid = nullptr,
      Functional functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::KINFUNCTIONALS::LLP91K),
      bool localizedEnv = false,
      std::pair<bool, std::vector<std::shared_ptr<EnergyComponentController>>> eCon = {
          true, std::vector<std::shared_ptr<EnergyComponentController>>(0)});
  /**
   * @brief Default destructor.
   */
  virtual ~LevelshiftHybridPotential() = default;

  /**
   * @brief Getter for the actual potential.
   *
   * This function makes use of the RI approximation.
   *
   * @return Returns the active Systems potential in matrix form.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Getter for the energy associated with this potential.
   *
   * In this case the energy is a linear correction to the projection
   * operator plus possible non-additive contributions.
   *
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };

 private:
  /// @brief The density controllers of the not-projected densities.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _notProjectedEnvDensities;
  /// @brief The level-shift potentials.
  std::vector<std::shared_ptr<LevelshiftPotential<SCFMode>>> _levelshiftPotentials;
  /// @brief The optional non-additive kinetic energy potential.
  std::shared_ptr<NAddFuncPotential<SCFMode>> _naddFuncPotential = nullptr;
  /// @brief Flag for the use of a hybrid method of level-shift and non-additive kin. potential.
  bool _isHybrid;

  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
};

} /* namespace Serenity */

#endif /* POTENTIALS_LEVELSHIFTHYBRIDPOTENTIAL_H_ */
