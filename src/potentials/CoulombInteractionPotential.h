/**
 * @file CoulombInteractionPotential.h
 * @author: Kevin Klahr
 *
 * @date 29. November 2016
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

#ifndef POTENTIALS_COULOMBINTERACTIONPOTENTIAL_H_
#define POTENTIALS_COULOMBINTERACTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/Potential.h" //Base class.
#include "settings/Options.h"     //SCF_MODES

namespace Serenity {

/* Forward Declarations */
class AtomCenteredBasisController;
class SystemController;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
class Libint;
namespace Options {
enum class DENS_FITS;
}

/**
 * @class CoulombInteractionPotential CoulombInteractionPotential.h
 *
 * A class for the Coulomb interaction potential of an active system with several environment systems.
 *
 *
 * Implementation according to:
 * [1] Weigend, F.; Kattannek, M.; Ahlrichs, R.; J.Chem.Phys (2009), 130, 164106 (eq. 4)
 *
 * See also:
 * [2] Neese, F.; J.Comput.Chem. (2003), 24, 1740 (Scheme 3)
 *
 */
template<Options::SCF_MODES SCFMode>
class CoulombInteractionPotential : public Potential<SCFMode>,
                                    public ObjectSensitiveClass<Basis>,
                                    public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor employing the RI approximation.
   * @param actSystem The active system.
   * @param envSystems The environment systems.
   * @param actBasis The basis controller of the active system.
   * @param envDensityMatrixController The density matrix controllers of the environment systems.
   * @param actAuxBasis The auxiliary basis sets of the active systems.
   * @param envAuxBasis The auxiliary basis sets of the environment systems.
   * @param isPassive If true the Fock matrix block is written to disk in order
   *                  to avoid recalculating it in future freeze-and-thaw iterations.
   */
  CoulombInteractionPotential(std::shared_ptr<SystemController> actSystem,
                              std::vector<std::shared_ptr<SystemController>> envSystems,
                              const std::shared_ptr<BasisController> actBasis,
                              std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensityMatrixController,
                              const std::shared_ptr<BasisController> actAuxBasis,
                              std::vector<std::shared_ptr<BasisController>> envAuxBasis, bool isPassive = false);
  /**
   * @brief Constructor without density fitting.
   * @param actSystem The active system.
   * @param envSystems The environment systems.
   * @param actBasis The basis controller of the active system.
   * @param envDensityMatrixController The density matrix controllers of the environment systems.
   * @param isPassive If true the Fock matrix block is written to disk in order
   *                  to avoid recalculating it in future freeze-and-thaw iterations.
   */
  CoulombInteractionPotential(std::shared_ptr<SystemController> actSystem,
                              std::vector<std::shared_ptr<SystemController>> envSystems,
                              const std::shared_ptr<BasisController> actBasis,
                              std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensityMatrixController,
                              bool isPassive = false);
  /**
   * @brief Default destructor.
   */
  virtual ~CoulombInteractionPotential() = default;

  /**
   * @brief Getter for the actual potential.
   * @return Returns the active system's potential in matrix form.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

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
  static std::vector<unsigned int>
  createBasisToAtomIndexMapping(const std::vector<std::pair<unsigned int, unsigned int>>& basisIndicesRed,
                                unsigned int nBasisFunctionsRed);
  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };

 private:
  ///@brief active system controller
  std::weak_ptr<SystemController> _actSystem;
  ///@brief environment system controller
  std::vector<std::weak_ptr<SystemController>> _envSystems;
  ///@brief Environment densities.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _envDMatController;
  ///@brief Active system  aux. basis set.
  const std::shared_ptr<BasisController> _actAuxBasis;
  ///@brief Environment aux. basis sets.
  std::vector<std::shared_ptr<BasisController>> _envAuxBasis;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  ///@brief switch for RI.
  Options::DENS_FITS _mode;
  ///@brief Flag for a passive potential.
  bool _isPassive = false;
  ///@brief Try to get the Fock matrix from disk.
  void getFromDisk();
  ///@brief Write the Fock matrix contribution to disk.
  void toHDF5();
  ///@brief Load the Fock matrix contribution from disk.
  bool fromHDF5();
  ///@brief Base name of the Fock matrix contribution file.
  std::string _fBaseName;
  ///@brief Calculate the Fock matrix contribution.
  void calculateFockMatrix();
};

} /* namespace Serenity */

#endif /* POTENTIALS_COULOMBINTERACTIONPOTENTIAL_H_ */
