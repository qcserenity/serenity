/**
 * @file   FXDTask.h
 *
 * @date   Mar 3, 2020
 * @author Johannes Toelle
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
#ifndef FXDTASK_H_
#define FXDTASK_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "settings/LRSCFOptions.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

using namespace Serenity::Reflection;

struct FXDTaskSettings {
  FXDTaskSettings()
    : loadType(Options::LRSCF_TYPE::ISOLATED),
      donoratoms({}),
      acceptoratoms({}),
      FED(false),
      FCD(false),
      multistateFXD(false),
      states(100),
      loewdinpopulation(true),
      writeTransformedExcitationVectors(false){};
  REFLECTABLE((Options::LRSCF_TYPE)loadType, (std::vector<unsigned int>)donoratoms,
              (std::vector<unsigned int>)acceptoratoms, (bool)FED, (bool)FCD, (bool)multistateFXD, (unsigned int)states,
              (bool)loewdinpopulation, (bool)writeTransformedExcitationVectors)
};
/**
 * @class FXDTask FXDTask.h
 *
 * Performs fragment charge difference (FCD), fragment excitation difference (FED) and multistate FED-FCD calculations
 * Literature:
 * FCD: J. Chem. Phys, 117, 5607 (2002)
 * FED: J. Phys. Chem. C, 112, 1204 (2008)
 * Multistate FED-FCD: Photosynth. Res., 137, 215 (2018)
 */
template<Options::SCF_MODES SCFMode>
class FXDTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param systemController
   */
  FXDTask(std::shared_ptr<SystemController> systemController);

  /**
   * @brief Default destructor.
   */
  virtual ~FXDTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief A function to calculate the FCD/FED matrices
   * @param excVector The excitation vector
   * @param bool Whether FED or matrix or FCD matrix is evaluated
   * @return Returns the FCD/FED matrix
   */
  std::vector<Eigen::MatrixXd> calculateFXDMatrix(std::shared_ptr<Eigen::MatrixXd> excVector, bool FED);

  // To access the data from the tests
  std::vector<Eigen::MatrixXd> _fcd_matrix;
  std::vector<Eigen::MatrixXd> _fed_matrix;
  Eigen::MatrixXd _fxd_matrix;

  /**
   * @brief The settings/keywords for the FXDTask:
   *       - loadType:           Loading type of the excitation vectors (isolated,uncoupled,coupled)
   *       - donoratoms:         Atom indices defining the donor fragment(starting from 0!)
   *       - acceptoratoms:      Atom indices defininf the acceptor atoms (starting from donoratoms+1)
   *       - FED:                Decides if a FED calculation should be performed
   *       - FCD:                Decides if a FCD calculation should be performed
   *       - multistateFXD:      Decides if a multistate FED-FCD calculation should be performed
   *       - states:             The number of states used for the multistate FED-FCD calculation
   *       - loewdinpopulation:  The population analyis used to construct the FED and FCD matrices (Loewdin/Mulliken)
   *       - writeTransformedExcitationVectors Transforms the excitation vectors in the diabatic basis and overwrites
   * the old excitation vectors
   */
  FXDTaskSettings settings;

 private:
  /// @brief The system controller.
  std::shared_ptr<SystemController> _systemController;
  ///@brief The lrscfController
  std::shared_ptr<LRSCFController<SCFMode>> _lrscf;
  ///@brief The excitation energies
  std::shared_ptr<Eigen::VectorXd> _excEner;
  ///@brief The excitation vectors
  std::shared_ptr<Eigen::MatrixXd> _excVecs;
};

} /* namespace Serenity */

#endif /* FXDTASK_H_ */