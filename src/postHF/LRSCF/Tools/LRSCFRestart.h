/**
 * @file LRSCFRestart.h
 *
 * @date Jan 6, 2020
 * @author Niklas Niemeyer
 *
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

#ifndef LRSCF_LRSCFRESTART
#define LRSCF_LRSCFRESTART

/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
#include "settings/LRSCFOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

struct LRSCFTaskSettings;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class LRSCFRestart
 * A static class supposed to keep the task neat and clean.\n
 *
 * Fetches and writes eigenpairs from and to disk.
 */
template<Options::SCF_MODES SCFMode>
class LRSCFRestart {
 public:
  /**
   * @brief Constructor.
   * @param lrscf LRSCF Controller.
   * @param settings The LRSCF settings of the LRSCFTask.
   * @param type What type of LRSCFTask is being run.
   */
  LRSCFRestart(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, const LRSCFTaskSettings& settings,
               Options::LRSCF_TYPE type);

  /**
   * @brief Destructor.
   */
  virtual ~LRSCFRestart() = default;

  /**
   * @brief Tries to find eigenpairs on disk.
   * @param settings The LRSCF settings of the LRSCFTask.
   * @param eigenvalues The eigenvalues.
   * @return The eigenvectors.
   */
  std::shared_ptr<std::vector<Eigen::MatrixXd>> fetchEigenpairs(Eigen::VectorXd& eigenvalues);

  /**
   * @brief Stores converged solution.
   * @param eigenvectors The eigenvectors.
   * @param eigenvalues The eigenvalues.
   */
  void storeConvergedSolution(std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors, Eigen::VectorXd eigenvalues);

  /**
   * @brief Stores converged response solution.
   * @param solutionvectors The solutionvectors.
   * @param frequencies The perturbing frequencies.
   */
  void storeConvergedResponse(std::shared_ptr<std::vector<Eigen::MatrixXd>> solutionvectors, Eigen::VectorXd frequencies);

  /**
   * @brief A lambda to print temporary iteration data to disk.
   * @return A lambda to print temporary iteration data to disk.
   */
  std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> getWriteToDisk();

 private:
  ///@brief The LRSCFController.
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;

  ///@brief The LRSCFTaskSettings.
  const LRSCFTaskSettings& _settings;

  ///@brief The type of the current LRSCFTask run.
  Options::LRSCF_TYPE _type;
};

} /* namespace Serenity */

#endif /* LRSCF_LRSCFRestart */
