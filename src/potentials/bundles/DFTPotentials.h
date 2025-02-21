/**
 * @file DFTPotentials.h
 *
 * @date Nov 24, 2016
 * @author: Jan Unsleber
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

#ifndef POTENTIALS_BUNDLES_DFTPOTENTIALS_H_
#define POTENTIALS_BUNDLES_DFTPOTENTIALS_H_

/* Include Serenity Internal Headers */
#include "potentials/bundles/PotentialBundle.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward Declarations */
class Geometry;
template<Options::SCF_MODES SCFMode>
class Potential;
template<Options::SCF_MODES SCFMode>
class FuncPotential;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
template<Options::SCF_MODES SCFMode>
class HCorePotential;

/**
 * @class DFTPotentials DFTPotentials.h
 * @brief A class containing all the potentials relevant for a KS-DFT-SCF.
 */
template<Options::SCF_MODES SCFMode>
class DFTPotentials : public PotentialBundle<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param hcore                 The one electron potential.
   * @param J                     The Coulomb potential.
   * @param Vxc                   The exchange-correlation potential coming from a functional.
   * @param pcm                   The implicit solvent model.
   * @param geom                  The geometry.
   * @param dMatController        The density matrix controller.
   * @param prescreeningThreshold The integral prescreening threshold.
   */
  DFTPotentials(std::shared_ptr<HCorePotential<SCFMode>> hcore, std::shared_ptr<Potential<SCFMode>> J,
                std::shared_ptr<FuncPotential<SCFMode>> Vxc, std::shared_ptr<Potential<SCFMode>> pcm,
                std::shared_ptr<const Geometry> geom, std::shared_ptr<DensityMatrixController<SCFMode>> dMatController,
                const double prescreeningThreshold);
  /// @brief Default destructor.
  virtual ~DFTPotentials() = default;

  /**
   * @brief A function to get the entire Fock matrix.
   * @param P The density matrix.
   * @param energies The controller to add all the energy contributions to.
   * @return The current Fock matrix (rebuilt on every call).
   *         (Note that this does not imply that each Potential is rebuilt, they
   *          may very well be cached. But all potentials contributing are added
   *          together in every call)
   */
  FockMatrix<SCFMode> getFockMatrix(const DensityMatrix<SCFMode>& P,
                                    std::shared_ptr<EnergyComponentController> energies) override final;

  /**
   * @brief Returns gradients of all underlying potentials
   */
  Eigen::MatrixXd getGradients() override final;
  /**
   * @brief Replace the functional potential.
   * @param newVxc The new funcational potential.
   */
  void replaceFunctionalPotential(std::shared_ptr<FuncPotential<SCFMode>> newVxc);
  /**
   * @brief Getter for the gradients of external point charges.
   * @return The point charge gradients.
   */
  Eigen::MatrixXd getPointChargeGradients() override;

 private:
  ///@brief The one electron potential.
  std::shared_ptr<HCorePotential<SCFMode>> _h;
  ///@brief The Coulomb potential.
  std::shared_ptr<Potential<SCFMode>> _J;
  ///@brief The exchange-correlation potential.
  std::shared_ptr<FuncPotential<SCFMode>> _Vxc;
  ///@brief The implcit solvation model.
  std::shared_ptr<Potential<SCFMode>> _pcm;
  ///@brief The geometry.
  std::shared_ptr<const Geometry> _geom;
  ///@brief The density matrix controller for this potential.
  std::shared_ptr<DensityMatrixController<SCFMode>> _dMatController;
  ///@brief Threshold for the integral prescreening.
  const double _prescreeningThreshold;
};

} /* namespace Serenity */

#endif /* POTENTIALS_BUNDLES_DFTPOTENTIALS_H_ */
