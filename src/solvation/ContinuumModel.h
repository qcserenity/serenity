/**
 * @file ContinuumModel.h
 *
 * @author Moritz Bensberg
 * @date Mai 18, 2020
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

#ifndef SOLVATION_CONTINUUMMODEL_H_
#define SOLVATION_CONTINUUMMODEL_H_

/* Include Serenity Internal Headers */
#include "data/grid/GridPotential.h"           //Can not use forward declaration due to this class being a typedef.
#include "data/matrices/DensityMatrix.h"       //Can not use forward declaration due to this class being a typedef.
#include "notification/ObjectSensitiveClass.h" //Notify.
#include "settings/PCMSettings.h"              //Default constructor of PCMSettings.
/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {
/* Forward Declarations */
class GridController;
class MolecularSurfaceController;
template<Options::SCF_MODES SCFMode>
class ElectrostaticPotentialOnGridController;
template<Options::SCF_MODES SCFMode>
class MatrixInBasis;

/**
 * @class ContinuumModel ContinuumModel.h
 * @brief Base class for any continuum model.
 */
template<Options::SCF_MODES SCFMode>
class ContinuumModel : public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor. Creates the input for PCMSolver.
   * @param pcmSettings                      The settings associated with the PCM.
   * @param molecularSurface                 The molecular surface.
   * @param activePotential                  The electrostatic potential controller of the active system.
   * @param environmentPotentials            The elecstrostatic potentials of the environment systems.
   */
  ContinuumModel(const PCMSettings& pcmSettings, std::shared_ptr<MolecularSurfaceController> molecularSurface,
                 std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>> activePotential,
                 std::vector<std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>> environmentPotentials = {});
  /**
   * @brief Default destructor.
   */
  ~ContinuumModel();
  /**
   * @brief Getter for the PCM charges on a grid.
   * @brief The PCM charges.
   */
  const GridPotential<RESTRICTED>& getPCMCharges();
  /**
   * @brief Getter for the interaction of the active system with the PCM charges.
   */
  double getActivePCMEnergy();
  /**
   * @brief Getter for the interaction of the total system (active + all environment)
   *        with the PCM charges.
   */
  double getTotalPCMEnergy();
  /**
   * @brief Force reinit of the partial charges, given that the charges are not loaded from a file.
   */
  void notify() override final {
    if (!_settings.loadedPCM)
      _pcmCharges = nullptr;
  };
  /**
   * @brief Getter for the underlying PCMSettings.
   */
  const PCMSettings& getPCMSettings();
  /**
   * @brief Getter for the molecular surface controller.
   */
  std::shared_ptr<MolecularSurfaceController> getMolecularSurfaceController();
  /**
   * @brief Getter for the CPCM scaling factor.
   */
  double getCPCMScaling();

 private:
  // The PCMSettings.
  PCMSettings _settings;
  // The molecular surface.
  std::shared_ptr<MolecularSurfaceController> _molecularSurface;
  // The electrostatic potential of the active system.
  std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>> _activePotential;
  // The electrostatic potential of the environment systems.
  std::vector<std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>> _environmentPotentials;
  // The PCM charges defined on the cavity grid.
  std::shared_ptr<GridPotential<RESTRICTED>> _pcmCharges;
  // The matrix K. The PCM charges (P) are given
  // with the electrostatic potential V as -K V = P.
  std::shared_ptr<Eigen::MatrixXd> _K;
  // Calculate K.
  void decomposeCavityMatrix();
  // Calculate the diielectric energy.
  double calculateEnergy(const GridPotential<RESTRICTED>& pcmCharges, const GridPotential<RESTRICTED>& potential);
  // The static permittivity.
  double _eps;
};

} /* namespace Serenity */

#endif /*  SOLVATION_CONTINUUMMODEL_H_ */
