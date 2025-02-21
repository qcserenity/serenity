/**
 * @file QuasiRestrictedOrbitalsTask.h
 *
 * @date Mai 10, 2021
 * @author Moritz Bensberg
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

#ifndef SERENITY_QUASIRESTRICTEDORBITALSTASK_H
#define SERENITY_QUASIRESTRICTEDORBITALSTASK_H
/* Include Serenity Internal Headers */
#include "settings/EmbeddingSettings.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {
class SystemController;
using namespace Serenity::Reflection;

struct QuasiRestrictedOrbitalsTaskSettings {
  QuasiRestrictedOrbitalsTaskSettings() : canonicalize(true) {
  }
  REFLECTABLE((bool)canonicalize)
 public:
  EmbeddingSettings embedding;
};
/**
 * @class
 * @brief Calculates natural UHF/UKS orbitals and quasi-restricted orbitals.
 *        Note that this task will only do something for SCFMode = UNRESTRICTED.
 *        Otherwise it will simply skip any calculations.
 *
 * Quasi-restricted orbitals:\n
 * Diagonalize UHF/UKS-density matrix:\n
 * A) MOs with occ. number n= 1.0: SOMOs.\n
 * B) MOs with occ. number n~2.0:  DOMOs (double occ. orbitals).\n
 * C) MOs with occ. number n~0.0:  virtual.\n
 * Canonicalization.\n
 * DOMOs: diag. F-beta.\n
 * virtual: diag. F-alpha.\n
 * SOMOs: diag. (F-alpha + F-beta)/2, two orbital energies:\n
 *        e-alpha = <SOMO|f-alpha|SOMO> and e-beta = <SOMO|f-beta|SOMO>\n
 * According to J. Am. Chem. Soc., 128, 10213-10222 (2006)
 */
template<Options::SCF_MODES SCFMode>
class QuasiRestrictedOrbitalsTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param activeSystem The active system.
   * @param environmentSystems The environment systems (energy calculation only).
   */
  QuasiRestrictedOrbitalsTask(std::shared_ptr<SystemController> activeSystem,
                              std::vector<std::shared_ptr<SystemController>> environmentSystems);
  /**
   * @brief Default destructor.
   */
  ~QuasiRestrictedOrbitalsTask() = default;
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
  void visit(QuasiRestrictedOrbitalsTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the blockname is unknown.
    throw SerenityError((std::string) "Unknown block in QuasiRestrictedOrbitalsTaskSettings: " + blockname);
  }
  /**
   * @brief Settings.
   * - canonicalize: If true, the orbitals are re-canonicalized as described above.
   * - embedding:    Embedding settings for the energy calculation.
   */
  QuasiRestrictedOrbitalsTaskSettings settings;
  /**
   * @brief Getter for the number of DOMOs.
   * @return The number of DOMOs.
   */
  unsigned int getNDOMOs();
  /**
   * @brief Getter for the number of SOMOs.
   * @return The number of SOMOs.
   */
  unsigned int getNSOMOs();
  /**
   * @brief Getter for the number of purely virtual orbitals.
   * @return The number of purely virtual orbitals.
   */
  unsigned int getNVirtuals();
  /**
   * @brief Getter for the UHF/UKS occupation numbers.
   * @return The occupation numbers.
   */
  Eigen::VectorXd getOccupationNumbers();
  /**
   * @brief Calculate the UHF/UKS natural orbitals for the active system.
   * @return The natural orbital coefficients.
   */
  Eigen::MatrixXd calculateNaturalOrbitals();

 private:
  std::shared_ptr<SystemController> _activeSystem;
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  void updateToNonCanonicalNaturalOrbitals(const Eigen::MatrixXd& naturalOrbitals);
  void canonicalizeOrbitals(unsigned int nDOMOs, unsigned int nSOMOs, unsigned int nVirt);
  void updateEnergy();
  /*
   * Initialized upon natural orbital calculation.
   */
  bool _initialized = false;
  unsigned int _nDOMOs = 0;
  unsigned int _nSOMOs = 0;
  unsigned int _nVirtuals = 0;
  Eigen::VectorXd _occupationNumbers;
};

} /* namespace Serenity */
#endif // SERENITY_QUASIRESTRICTEDORBITALSTASK_H
