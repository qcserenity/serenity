/**
 * @file NEInteractionPotential.h
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

#ifndef POTENTIALS_NEINTERACTIONPOTENTIAL_H_
#define POTENTIALS_NEINTERACTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "geometry/Geometry.h"
#include "potentials/Potential.h"

namespace Serenity {
/* Forward declaration */
class SystemController;
/**
 * @class NEInteractionPotential NEInteractionPotential.h
 * @brief The attraction potential between electrons in a given
 *        basis and sets of nucleii.
 */
template<Options::SCF_MODES SCFMode>
class NEInteractionPotential : public Potential<SCFMode>, public ObjectSensitiveClass<Atom> {
 public:
  /**
   * @brief Constructor.
   * @param basis      The basis of the potential.
   * @param geometries The nucleii in form of multiple geometries.
   */
  NEInteractionPotential(std::shared_ptr<SystemController> activeSystem,
                         std::vector<std::shared_ptr<SystemController>> environmentSystems,
                         std::shared_ptr<BasisController> basis, std::vector<std::shared_ptr<const Geometry>> geometries);
  /// @brief Default destructor.
  virtual ~NEInteractionPotential() = default;

  /**
   * @brief Getter for the actual potential.
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
   * @brief Geometry gradient contribution from this Potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Matches each basis function shell to its respective atom center.
   *
   * @param basisIndicesRed see AtomCenteredBasisController
   * @param nBasisFunctionRed the (reduced) number of basis functions
   */
  static std::vector<unsigned int>
  createBasisToAtomIndexMapping(const std::vector<std::pair<unsigned int, unsigned int>>& basisIndicesRed,
                                unsigned int nBasisFunctionsRed);

  /**
   * @brief Potential is linked to the grid and density.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };

 private:
  ///@brief The active systems controller.
  std::weak_ptr<SystemController> _actSystem;
  ///@brief The environment systems controllers.
  std::vector<std::weak_ptr<SystemController>> _envSystems;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  ///@brief The geometries.
  std::vector<std::shared_ptr<const Geometry>> _geometries;
};

} /* namespace Serenity */

#endif /* POTENTIALS_NEINTERACTIONPOTENTIAL_H_ */
