/**
 * @file EvaluateEnergyTask.h
 *
 * @author Moritz Bensberg, Anja Massolle
 * @date Mar 10, 2020
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

#ifndef TASKS_EVALUATEENERGYTASK_H_
#define TASKS_EVALUATEENERGYTASK_H_

/* Include Serenity Internal Headers */
#include "dft/functionals/CompositeFunctionals.h"
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //Local correlation settings for local MP2.
#include "settings/EmbeddingSettings.h"
#include "settings/OrthogonalizationOptions.h"
#include "settings/Reflection.h"
#include "settings/Settings.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/* Forward Declarations */
class SystemController;

using namespace Serenity::Reflection;
struct EvaluateEnergyTaskSettings {
  EvaluateEnergyTaskSettings()
    : evalTsOrtho(false),
      evalAllOrtho(false),
      orthogonalizationScheme(Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN),
      useDifferentXCFunc(false),
      XCfunctional(CompositeFunctionals::XCFUNCTIONALS::BP86),
      mp2Type(Options::MP2_TYPES::DF),
      maxResidual(1e-5),
      maxCycles(100) {
    lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
    lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  }
  REFLECTABLE((bool)evalTsOrtho, (bool)evalAllOrtho, (Options::ORTHOGONALIZATION_ALGORITHMS)orthogonalizationScheme,
              (bool)useDifferentXCFunc, (CompositeFunctionals::XCFUNCTIONALS)XCfunctional, (Options::MP2_TYPES)mp2Type,
              (double)maxResidual, (int)maxCycles)
 public:
  LocalCorrelationSettings lcSettings;
  EmbeddingSettings embedding;
  CUSTOMFUNCTIONAL customFunc;
};
/**
 * @class EvaluateEnergyTask EvaluateEnergyTask.h
 * @brief Evaluates the energy for the given system, density and settings.
 */
template<Options::SCF_MODES SCFMode>
class EvaluateEnergyTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param systemController The system for which the energy should be evaluated.
   */
  EvaluateEnergyTask(std::vector<std::shared_ptr<SystemController>> systemController,
                     std::shared_ptr<SystemController> superSystem = nullptr);
  /**
   * @brief Default destructor.
   */
  virtual ~EvaluateEnergyTask() = default;

  /**
   * @see Task
   */
  void run();
  /**
   * @brief Settings.
   *        - mp2Type:               Type of MP2 used for double-hybrid correlation part.
   *        - maxResidual:           Maximum residual threshold for local MP2.
   *        - maxCycles:             Maximum number of iterations before cancelling the amplitude optimization
   *                                in local MP2.
   *        - lcSettings:            Local correlation settings for local MP2.
   *        - evalTsOrtho:           if true, the non-additive kinetic energy is evaluated from orthogonalized subsystem
   * orbitals.
   *        - evalAllOrtho:          if true, all energy contributions are evaluated from orthogonalized subsystem
   * orbitals.
   *        - orthogonalizationScheme: the orthogonalization procedure used for evalTsOrtho and evalAllOrtho
   *        - useDifferentXCFunc:    if true a different XC functional than the one defined in the system is used for
   * the energy evaluation. - XCfunctional:          The XC functional which shall be used for the energy evaluation
   */
  EvaluateEnergyTaskSettings settings;
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(EvaluateEnergyTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.customFunc.visitAsBlockSettings(v, blockname)) {
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname)) {
      return;
    }
    if (c.lcSettings.visitAsBlockSettings(v, blockname)) {
      return;
    }
    throw SerenityError((std::string) "Unknown block in EvaluateEnergyTaskSettings: " + blockname);
  }

 private:
  // The system controller.
  std::vector<std::shared_ptr<SystemController>> _systemController;
  // The system controller for the supersystem used for the orthogonalization
  std::shared_ptr<SystemController> _superSystem;
  // calculates the non-additive kinetic energy from orthogonalized supersystem orbitals
  double calcNaddKin(std::shared_ptr<SystemController> supersystem, std::vector<std::shared_ptr<SystemController>> subsystems);
};

} /* namespace Serenity */

#endif /* TASKS_EVALUATEENERGYTASK_H_ */
