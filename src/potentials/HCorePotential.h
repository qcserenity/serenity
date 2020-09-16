/**
 * @file HCorePotential.h
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

#ifndef POTENTIALS_HCOREPOTENTIAL_H_
#define POTENTIALS_HCOREPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/Potential.h"
#include "settings/Options.h"
#include "system/SystemController.h"

namespace Serenity {
/**
 * @class HCorePotential HCorePotential.h
 */
template<Options::SCF_MODES SCFMode>
class HCorePotential : public Potential<SCFMode>, public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @brief Constructor
   * @param system The system this potential is defined for.
   */
  HCorePotential(std::shared_ptr<SystemController> system);
  /// @brief Default destructor.
  virtual ~HCorePotential() = default;
  /**
   * @brief Getter for the actual potential.
   *
   * This function uses references to a MatrixInBasis
   * because this way the unrestricted version can still hold
   * two matrices for potentials that are equal for alpha and beta
   * spin without two copies of the same matrix.
   *
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override final;

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
  std::vector<unsigned int> createBasisToAtomIndexMapping(const std::vector<std::pair<unsigned int, unsigned int>>& basisIndicesRed,
                                                          unsigned int nBasisFunctionsRed);

  /**
   * @brief Calculates the energy weighted density matrix.
   * @param systemController The system.
   * @param orbitalSet The chosen orbitals.
   * @return Returns the energy weighted density matrix.
   */
  DensityMatrix<SCFMode> calcEnergyWeightedDensityMatrix(std::shared_ptr<SystemController> systemController,
                                                         const std::shared_ptr<OrbitalController<SCFMode>>& orbitalSet);

  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };

 private:
  ///@brief The system potential is defined for.
  std::weak_ptr<SystemController> _system;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
};

} /* namespace Serenity */

#endif /* POTENTIALS_HCOREPOTENTIAL_H_ */
