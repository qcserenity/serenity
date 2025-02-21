/**
 * @file PCMPotential.h
 *
 * @author Moritz Bensberg
 * @date Apr 28, 2020
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

#ifndef POTENTIALS_PCMPOTENTIAL_H_
#define POTENTIALS_PCMPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/grid/GridPotential.h"           //Can not use forward declaration due to this class being a typedef.
#include "notification/ObjectSensitiveClass.h" //Notify.
#include "potentials/Potential.h"              //Base class.
#include "settings/Options.h"                  //SCF_MODES
/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {
/* Forward Declarations */
class BasisController;
class Geometry;
struct PCMSettings;
template<Options::SCF_MODES SCFMode>
class ContinuumModel;
template<Options::SCF_MODES SCFMode>
class ElectrostaticPotentialOnGridController;
class MolecularSurfaceController;
template<Options::SCF_MODES SCFMode>
class OrbitalController;

/**
 * @class PCMPotential PCMPotential.h
 * @brief A class that handles Fock matrix, energy and gradient contributions that
 *        originate from an implicit solvent model.
 */
template<Options::SCF_MODES SCFMode>
class PCMPotential : public Potential<SCFMode>,
                     public ObjectSensitiveClass<Basis>,
                     public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param pcmSettings                   The settings associated with the PCM
   * @param basisController               The basis controller to express the Fock matrix in.
   * @param geometry                      The geometry of the active system.
   * @param activeDensityMatrixController The active density matrix controller.
   * @param densityMatrixControllers      The density matrix controller of the envrionments.
   * @param environmentGeometries         The geometries of the environemt.
   */
  PCMPotential(const PCMSettings& pcmSettings, std::shared_ptr<BasisController> basisController,
               std::shared_ptr<const Geometry> geometry, std::shared_ptr<MolecularSurfaceController> molecularSurface = nullptr,
               std::shared_ptr<MolecularSurfaceController> vdWmolecularSurface = nullptr,
               std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>> activePotential = nullptr,
               std::vector<std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>> environmentPotentials = {});
  /**
   * @brief Getter for the Fock matrix.
   * @return Returns the Fock matrix.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override;
  /**
   * @brief Getter for the energy. Uses the active electrostatic potential.
   * @return The energy.
   */
  virtual double getEnergy(const DensityMatrix<SCFMode>& P) override;
  /**
   * @brief Getter for the total dielectic energy. Uses the total electrostatic potential.
   * @return The energy.
   */
  double getTotalEnergy();
  /**
   * @brief Destructor. Frees libint engines.
   */
  ~PCMPotential();
  /**
   * @brief Potential is linked to the basis etc. it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };
  /**
   * @brief Geometry gradient contribution from this Potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override final;

 private:
  ///@brief Calculate the cavity formation energy.
  bool _calculateG_cav;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  ///@brief Integrate the PCM charges to the final Fock matrix.
  Eigen::MatrixXd integrateToFockMatrix(const GridPotential<RESTRICTED>& pcmCharges);
  ///@brief The wrapper for the external library.
  std::shared_ptr<ContinuumModel<SCFMode>> _continuumModel;
  ///@brief The geometry (needed for the gradients).
  std::shared_ptr<const Geometry> _geometry;
  ///@brief The vdW surface used for the calculation of the cavity formation energy.
  std::shared_ptr<MolecularSurfaceController> _vdWmolecularSurface;
  ///@brief The controller of the electrostatic potential of the active system. This
  /// is used to access the two center integrals used by it.
  std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>> _activePotential;
  ///@brief The electrostatic potential controller for the environment systems.
  std::vector<std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>> _environmentPotentials;
};

} /* namespace Serenity */

#endif /* POTENTIALS_PCMPOTENTIAL_H_ */
