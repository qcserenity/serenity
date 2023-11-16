/**
 * @file   FDEETTask.h
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
#ifndef FDEETTASK_H_
#define FDEETTASK_H_

/* Include Serenity Internal Headers */
#include "settings/DFTOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"

namespace Serenity {
/* Forward declarations */
class SystemController;
class FDEETController;
class FDEDiabController;
class FDEETCalculator;

using namespace Serenity::Reflection;

struct FDEETTaskSettings {
  FDEETTaskSettings()
    : functional(CompositeFunctionals::XCFUNCTIONALS::PW91),
      couple({}),
      states({}),
      population({0}),
      spindensity(false),
      spinpopulation(false),
      disjoint({}),
      invThreshold(1e-3),
      diskMode(false),
      printContributions(false){};
  REFLECTABLE((CompositeFunctionals::XCFUNCTIONALS)functional, (std::vector<std::vector<unsigned int>>)couple,
              (std::vector<unsigned int>)states, (std::vector<unsigned int>)population, (bool)spindensity, (bool)spinpopulation,
              (std::vector<unsigned int>)disjoint, (double)invThreshold, (bool)diskMode, (bool)printContributions)
};

/**
 * @class  FDEETTask FDEETTask.h
 * @brief  Perform a FDE-ET calculation for a given set of subsystems and states.
 *         All subsystems have to be set to active, environment systems will be ignored.
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
 * FDE-Diab procedure:
 * [4] D. G. Artiukhin, J. Neugebauer. Frozen-density embedding as a quasi-diabatization
 *     tool: Charge-localized states for spin-density calculations. J. Chem. Phys., 148 (2018) 214104.\n
 *
 */
class FDEETTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param activeSystems  The active systems.
   */
  FDEETTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems);
  /**
   * @brief Default destructor.
   */
  virtual ~FDEETTask() = default;
  /**
   * @see Task.h
   */
  void run();

  /**
   * @brief The settings/keywords for the FDEETTask: \n
   *        @param couple Specifies which diabatic states are be coupled. E.g. couple{12} should be understood as
   *                    "couple diabatic states 1 and 2", whereas couple{12 13} means "do two separate runs: first
   *                     couple 1 and 2 second couple 1 and 3."
   *        @param states Specifies the number of subsystems contained in one diabatic state. states{2 2} can be,
   *                         thus understood that we have two diabatic states, both containing two subsystems. The first
   * state contains the first two systems specified in the active systems, whereas the second state contains the next
   * two systems. In all states the same subsystems must occur!!
   *        @param spindensity Specifies if adiabatic spin density matrices can be calculated and written to
   *                            cube-files (Default: false).
   *        @param spinpopulation Specifies if adiabatic spin populations can be calculated and printed (Default:
   *                               false).
   *        @param population Specifies in a list of which adiabatic states spin populations can be calculated and
   * printed (Default: 0).
   *        @param invThreshold Specifies the threshold for the Moore-Penrose inversion (Default: 1e-3).
   *        @param printContributions Specifies if the alpha and beta contributions of the diabatic density matrices
   *                                   can be written to cube-files (Default: false).
   *        @param disjoint Specifies which subsystems are joined in case of using the disjoint approximation. E.g.
   * disjoint{1 2} joins systems 1 and 2 in the list of systems for each state. See Ref.~[2] for more details about the
   * disjoint approximation
   *
   */
  FDEETTaskSettings settings;

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
  std::shared_ptr<FDEETController> _fdeetController;
  ///@brief The FDEETCalculator
  std::shared_ptr<FDEETCalculator> _fdeetCalculator;
  ///@brief The FDEDiabController
  std::shared_ptr<FDEDiabController> _fdeDiabController;
};

} /* namespace Serenity */

#endif /* FDEETTASK_H_ */
