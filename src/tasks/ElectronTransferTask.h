/**
 * @file   ElectronTransferTask.h
 *
 * @date   Oct 25, 2019
 * @author Patrick Eschenbach, Niklas Niemeyer
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef ELECTRONTRANSFERTASK_H_
#define ELECTRONTRANSFERTASK_H_

/* Include Serenity Internal Headers */
#include "settings/DFTOptions.h"
#include "settings/Reflection.h"
#include "settings/Settings.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {
/* Forward declarations */
class SystemController;
class FDEETController;
class FDEDiabController;
class FDEETCalculator;

using namespace Serenity::Reflection;

struct ElectronTransferTaskSettings {
  ElectronTransferTaskSettings()
    : couple({}),
      states({}),
      population({0}),
      spindensity(false),
      spinpopulation(false),
      disjoint({}),
      diskMode(false),
      printContributions(false),
      useHFCoupling(false),
      coupleAdiabaticStates(false),
      configurationWeights({}){};
  REFLECTABLE((std::vector<std::vector<unsigned int>>)couple, (std::vector<unsigned int>)states,
              (std::vector<unsigned int>)population, (bool)spindensity, (bool)spinpopulation,
              (std::vector<unsigned int>)disjoint, (bool)diskMode, (bool)printContributions, (bool)useHFCoupling,
              (bool)coupleAdiabaticStates, (std::vector<std::vector<double>>)configurationWeights)
};

/**
 * @class  ElectronTransferTask ElectronTransferTask.h
 * @brief  Perform a FDE-ET or ALMO-ET calculation for a given set of (sub)systems and states.
 *         All (sub)systems have to be set to active, environment systems will be ignored.
 *         The number of subsystems must be the same for all states. \n
 * References:\n
 *
 * FDE-ET procedure:\n
 * [1] M. Pavanello, J. Neugebauer. Modelling charge transfer reactions with the frozen
 *     density embedding formalism. J. Chem. Phys., 135 (2011) 234103.\n
 *
 * [2] M. Pavanello, T. Van Voorhis, L. Visscher, J. Neugebauer. An accurate and linear-
 *     scaling method for calculating charge-transfer excitation energies and diabatic
 *     couplings. J. Chem. Phys., 138 (2013) 054101.\n
 *
 * [3] P. Ramos, M. Papadakis, M. Pavanello. Performance of Frozen Density Embedding
 *     for Modeling Hole Transfer Reactions. J. Phys. Chem. B, 119 (2015) 7541–7557.\n
 *
 * FDE-diab procedure:
 * [4] D. G. Artiukhin, J. Neugebauer. Frozen-density embedding as a quasi-diabatization
 *     tool: Charge-localized states for spin-density calculations. J. Chem. Phys., 148 (2018) 214104.\n
 *
 * ALMO MS-DFT procedure:
 * [5] Y. Mao, A. Motoya-Castillo, T. E. Markland. Accurate and efficient DFT-based diabatization
 *     for hole and electron transfer using absolutely localized molecular orbitals. J. Chem. Phys., 151 (2018) 164114.\n
 *
 * [6] A. Cembran, L. Song, Y. Mo, J. Gao. Block-Localized Density Functional Theory (BLDFT), Diabatic Coupling, and
 *     Their Use in Valence Bond Theory for Representing Reactive Potential Energy Surfaces.
 *     J. Chem. Theory Comput., 5 (2009) 2702-2716.\n
 *
 */
class ElectronTransferTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param activeSystems  The active systems.
   */
  ElectronTransferTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems);
  /**
   * @brief Default destructor.
   */
  virtual ~ElectronTransferTask() = default;
  /**
   * @see Task.h
   */
  void run();

  /**
   * @brief The settings/keywords for the ElectronTransferTask: \n
   *        @param couple Specifies which diabatic states are be coupled. E.g. "couple {1 2}" should be understood as
   *                    "couple diabatic states 1 and 2", whereas "couple {1 2 ; 1 3}" means "do two separate runs:
   * first couple 1 and 2 second couple 1 and 3." (If none is given: All states are coupled.)
   *        @param states Specifies the number of subsystems contained in one diabatic state. states{2 2} can be,
   *                         thus understood that we have two diabatic states, both containing two subsystems. The first
   * state contains the first two systems specified in the active systems, whereas the second state contains the next
   * two systems. In all states the same subsystems must occur!! (If none is given: Each state consist of one
   * (sub)system.)
   *        @param spindensity Specifies if adiabatic spin density matrices can be calculated and written to
   *                            cube-files (Default: false).
   *        @param spinpopulation Specifies if adiabatic spin populations can be calculated and printed (Default:
   *                               false).
   *        @param population Specifies in a list of which adiabatic states spin populations can be calculated and
   * printed (Default: 0).
   *        @param diskMode Specifies whether the transition-density matrices are written to HDF5 file or kept in
   * memory. (Default: false).
   *        @param printContributions Specifies if the alpha and beta contributions of the diabatic density matrices
   *                                   can be written to cube-files (Default: false).
   *        @param disjoint Specifies which subsystems are joined in case of using the disjoint approximation. E.g.
   * disjoint{1 2} joins systems 1 and 2 in the list of systems for each state. See Ref.~[2] for more details about the
   * disjoint approximation
   *
   *        @param useHFCoupling Specifies whether the off-diagonal elements of the Hamilton matrix are calculated
   * according to Ref.~[6].
   *        @param coupleAdiabaticStates Specifies whether the adiabatic states of different runs (see keyword "couple")
   * should be coupled again.
   *        @param configurationWeights Specifies weights of the adiabatic states of that are coupled. Must have the
   * same dimensions as keyword "couple". Requires "coupleAdiabaticStates".
   *
   */
  ElectronTransferTaskSettings settings;

  /**
   * @brief Set up the supersystem from the subsystems contained in each diabatic state.
   * @param coupleSystems The subsystems that are contained in each diabatic state.
   * @return The constructed supersystem
   */
  std::shared_ptr<SystemController> setUpSupersystem(std::vector<std::shared_ptr<SystemController>> coupleSystems,
                                                     std::string name);

  /**
   * @brief Getter for the FDEETController.
   * @return The FDEETController
   */
  std::shared_ptr<FDEETController> getFDEETController();
  /**
   * @brief Getter for the FDEETCalculator.
   * @return The FDEETCalculator
   */
  std::shared_ptr<FDEETCalculator> getFDEETCalculator();
  /**
   * @brief Getter for the FDEDiabController
   * @return The FDEDiabController
   */
  std::shared_ptr<FDEDiabController> getFDEDiabController();

 private:
  ///@brief The active systems
  std::vector<std::shared_ptr<SystemController>> _act;
  ///@brief The FDEETController
  std::shared_ptr<FDEETController> _fdeetController = nullptr;
  ///@brief The FDEETCalculator
  std::shared_ptr<FDEETCalculator> _fdeetCalculator = nullptr;
  ///@brief The FDEDiabController
  std::shared_ptr<FDEDiabController> _fdeDiabController = nullptr;
};

} /* namespace Serenity */

#endif /* ELECTRONTRANSFERTASK_H_ */
