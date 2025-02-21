/**
 * @file   FiniteFieldTask.h
 *
 * @date   Apr 29, 2021
 * @author Patrick Eschenbach, Niklas Niemeyer
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
#ifndef FINITEFIELDTASK_H_
#define FINITEFIELDTASK_H_
/* Include Serenity Internal Headers */
#include "dft/functionals/CompositeFunctionals.h" //functionals in the embedding block
#include "settings/EmbeddingOptions.h"            //settings in the embedding block
#include "settings/EmbeddingSettings.h"           //settings in the embedding block
#include "settings/Reflection.h"                  //keywords for the task
#include "tasks/Task.h"                           //task base class
/* Include Std and External Headers */
#include <Eigen/Dense> //Eigen::Matrix objects
#include <vector>      //std::vector objects

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace Serenity::Reflection;

enum class BASE_PROPERTY { DIPOLE_MOMENT = 1, DYNAMIC_POLARIZABILITY = 2 };

struct FiniteFieldTaskSettings {
  FiniteFieldTaskSettings() : finiteFieldStrength(1.0e-2), frequency(0), hyperPolarizability(false) {
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::NONE;
    embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  }
  REFLECTABLE((double)finiteFieldStrength, (double)frequency, (bool)hyperPolarizability)
  EmbeddingSettings embedding;
};

/**
 * @class  FiniteFieldTask FiniteFieldTask.h
 * @brief  A Task that can perform finite field calculations on molecules to obtain dipole moments, polarizabilities and
 * hyperpolarizabilities
 */
template<Options::SCF_MODES SCFMode>
class FiniteFieldTask : public Task {
 public:
  /**
   * @brief Default constructor
   * @param activeSystem The active system
   * @param environmentSystems A vector of all environment SystemControllers.
   */
  FiniteFieldTask(const std::shared_ptr<SystemController>& activeSystem,
                  const std::vector<std::shared_ptr<SystemController>>& environmentSystems = {});
  /**
   * @brief Default destructor
   */
  virtual ~FiniteFieldTask() = default;
  /**
   * @brief see Task.h
   */
  void run();

  /**
   * @brief The settings for this task.
   * @param finiteFieldStrength The difference of negative and positive electric field strength applied.
   * @param frequency Frequency chosen for the calculation of static/dynamic polarizabilities
   * @param hyperPolarizability Boolean that decides whehter hyper polarizabilities are calculated
   *                            numerically based on the polarizability.
   */
  FiniteFieldTaskSettings settings;

  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(FiniteFieldTaskSettings& c, set_visitor v, std::string blockname);

  /**
   * @brief Return the polarizability tensor
   * @return Eigen::MatrixXd containing the polarizability tensor
   */
  Eigen::Matrix3d getPolarizability();
  /**
   * @brief Return the hyper polarizability tensor
   * @return std::vector<Eigen::MatrixXd> containing the hyper polarizability tensor for each cartesian coordinate
   */
  std::vector<Eigen::Matrix3d> getHyperPolarizability();
  /**
   * @brief Return the isotropic polarizability
   * @return double containing the isotropic polarizability
   */
  double getIsotropicPolarizability();

 private:
  ///@brief performs a scf under the influence of the electric field and returns the desired property as a matrix
  Eigen::MatrixXd perturbedSCF(unsigned direction, double fStrength, BASE_PROPERTY propertyType, double frequency);
  ///@brief The active system.
  std::shared_ptr<SystemController> _activeSystem;
  ///@brief The vector of environment systemControllers.
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  ///@brief The embedding scheme used for the calculation.
  Options::EMBEDDING_SCHEME _embeddingScheme;
  ///@brief The polarizability
  Eigen::Matrix3d _polarizability;
  ///@brief The hyper polarizability vector
  std::vector<Eigen::Matrix3d> _hyperPolarizability;
  ///@brief the isotropic polarizability
  double _isotropicPolarizability;
};

} /* namespace Serenity */

#endif /* FINITEFIELDTASK_H_ */
