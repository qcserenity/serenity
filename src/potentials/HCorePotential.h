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

namespace Serenity {
/* Forward declaration */
class SystemController;
class OneElectronIntegralController;
template<Options::SCF_MODES>
class OrbitalController;
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
  virtual ~HCorePotential();

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
   * @brief Getter for the gradient contribution with respect to the point charge positions.
   * @return The point charge gradients.
   */
  const Eigen::MatrixXd& getPointChargeGradients();

  /**
   * @brief Calculates the energy weighted density matrix.
   *
   * Energy weighted density matrix \f$ W_{\alpha\beta\sigma} = \sum_{\mu\nu} D_{\alpha\mu\sigma} F_{\mu\nu\sigma}
   * D_{\nu\beta\sigma} = \sum_i^{occ} \epsilon_{i\sigma} c_{\alpha i\sigma} c_{\beta i\sigma} \f$, with the latter
   * equality only holding for canonical orbitals. For the restricted case, it has an overall factor of 2, but since two
   * factors of 2 are pulled into the density matrices, it becomes a factor of one half in front of the density matrix -
   * fock matrix expression.
   * @return Returns the energy weighted density matrix.
   */
  DensityMatrix<SCFMode> calcEnergyWeightedDensityMatrix();

  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
    _pointChargeGradients = nullptr;
  };

 private:
  ///@brief The system potential is defined for.
  std::weak_ptr<SystemController> _system;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  ///@brief The one electron integral controller. We keep this as a shared pointer to
  /// ensure that it is not deleted by the remembering factory and all cached integrals
  /// are deleted.
  std::shared_ptr<OneElectronIntegralController> _oneElectronIntegrals;
  ///@brief Read external potential defined on a grid from file.
  void importExternalGridPotential(std::string inputFile);
  ///@brief If true, the potential includes a contribution from external charges.
  bool _hasExternalCharges = false;
  ///@brief Getter for all charges (external and nuclei). Charge ordering: First all atoms followed by all external
  /// charges.
  std::vector<std::pair<double, Point>> getAllCharges();
  ///@brief Gradient contribution with respect to the point charge positions.
  std::unique_ptr<Eigen::MatrixXd> _pointChargeGradients;
};

} /* namespace Serenity */

#endif /* POTENTIALS_HCOREPOTENTIAL_H_ */
