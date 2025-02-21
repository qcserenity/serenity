/**
 * @file BrokenSymmetryTask.h
 *
 * @date Feb 24, 2020
 * @author Anja Massolle
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

#ifndef TASKS_BROKENSYMMETRYTASK_H_
#define TASKS_BROKENSYMMETRYTASK_H_

/* Include Serenity Internal Headers */
#include "dft/functionals/CompositeFunctionals.h"
#include "settings/EmbeddingOptions.h"
#include "settings/EmbeddingSettings.h"
#include "settings/LocalizationOptions.h"
#include "settings/OrthogonalizationOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace Serenity::Reflection;

struct BrokenSymmetryTaskSettings {
  BrokenSymmetryTaskSettings()
    : nA(1),
      nB(1),
      embeddingScheme(Options::EMBEDDING_SCHEME::NONE),
      maxCycles(50),
      convThresh(1.0e-6),
      noThreshold(0.2),
      printLevel(2),
      evalTsOrtho(false),
      evalAllOrtho(false),
      orthogonalizationScheme(Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN),
      locType(Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  }

  REFLECTABLE((unsigned int)nA, (unsigned int)nB, (Options::EMBEDDING_SCHEME)embeddingScheme,
              (bool)calculateSolvationEnergy, (Options::KIN_EMBEDDING_MODES)embeddingMode, (unsigned int)maxCycles,
              (double)convThresh, (double)noThreshold, (unsigned int)printLevel, (bool)evalTsOrtho, (bool)evalAllOrtho,
              (Options::ORTHOGONALIZATION_ALGORITHMS)orthogonalizationScheme, (Options::ORBITAL_LOCALIZATION_ALGORITHMS)locType)
 public:
  EmbeddingSettings embedding;
};
/**
 * @class BrokenSymmetryTask BrokenSymmetryTask.h
 * @brief A task for Broken-Symmetry calculations. It calculates magnetic exchange couplings using DFT and the
 * broken-symmetry approximation. When the embedding scheme: NONE is used the LNO guess (According to: Shoji M., Chem.
 * Phys. Lett. 608, 50-54) is performed to obtain the start orbitals for the broken-symmetry calculation. Otherwise each
 * subsystem is taken as a spin site and the orbitals of one subsystem are flipped to converge to the broken-symmetry
 * solution.
 */
class BrokenSymmetryTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param systemController the high-spin system as one supersystem or both spin sites as subsystems.
   * @param bsSystemController if specified these systems are loaded from disk as broken-symmetry systems. E.g.
   *                           for a different energy evaluation.
   */
  BrokenSymmetryTask(std::vector<std::shared_ptr<SystemController>> hsSystemController,
                     std::vector<std::shared_ptr<SystemController>> bsSystemController = {});
  /**
   * @brief Default destructor.
   */
  virtual ~BrokenSymmetryTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(BrokenSymmetryTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the blockname is unknown.
    throw SerenityError((std::string) "Unknown block in BrokenSymmetryTaskSettings: " + blockname);
  }

  /**
   * @brief The settings/keywords for the Broken-Symmetry Task:\n
   *  - nA: Number of unpaired electrons on the first spin site.
   *  - nB: Number of unpaired electrons on the second spin site.
   *  - embeddingScheme The kind of sDFT embedding: FDETask, FreezeAndThawTask or only Isolated subsystem densities.\n
   *  If the embeddingScheme is None, a standard BS-DFT calculation is performed.
   *  - maxCycles Max FreezeAndThawTask cycles.
   *  - convThresh FreezeAndThawTask convergence threshold.
   *  - noThreshold Threshold for the assignment of the NO orbitals as SONO, UONO and DONO.
   *  - printLevel The print level of the Task.
   *  - smallSupersystemGrid see FreezeAndThawTask.
   *  - evalTsOrtho if true, the non-additive kinetic energy is evaluated from orthogonalized subsystem orbitals.
   *  - evalAllOrtho if true, all energy contributions are evaluated from orthogonalized subsystem orbitals.
   *  - orthogonalizationScheme the orthogonalization procedure used for evalTsOrtho and evalAllOrtho.
   *  - locType the localization procedure applied for localizing the SONOs.
   */
  BrokenSymmetryTaskSettings settings;

 private:
  // The system Controller holding the HS state of the spin sites
  std::vector<std::shared_ptr<SystemController>> _hsSystemController;
  // The system Controller holding the BS states of the spin sites (if loaded from disk)
  std::vector<std::shared_ptr<SystemController>> _bsSystemController;
};

} /* namespace Serenity */

#endif /* TASKS_BROKENSYMMETRYTASK_H_ */
