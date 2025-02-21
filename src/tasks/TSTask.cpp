/**
 * @file TSTask.cpp
 *
 * @date May 22, 2017
 * @author Jan Unsleber
 * @copyright \n
 * This file is part of the program Serenity.\n\n
 * Serenity is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.\n\n
 * Serenity is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.\n\n
 * You should have received a copy of the GNU Lesser General
 * Public License along with Serenity.
 * If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "tasks/TSTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutput.h"
#include "io/HDF5.h"
#include "math/saddlepoint/Bofill.h"
#include "math/saddlepoint/QST.h"
#include "system/SystemController.h"
#include "tasks/GradientTask.h"
#include "tasks/HessianTask.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {

TSTask::TSTask(std::shared_ptr<SystemController> ts, const std::vector<std::shared_ptr<SystemController>> env)
  : _ts(ts), _env(env) {
  assert(_ts && "Something is wrong with the TS guess and its shared_ptr.");
}

void TSTask::run() {
  /* ================================================
   *   Evaluation Lambda for Energies and Gradients
   * ================================================ */
  assert(_ts);
  auto system = _ts;
  auto getData = [&system](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients, bool print) {
    unsigned int nAtoms = system->getGeometry()->getNAtoms();
    assert(parameters.size() == 3 * nAtoms);
    assert(gradients.size() == 3 * nAtoms);
    Eigen::MatrixXd coordinates = Eigen::Map<const Eigen::MatrixXd>(parameters.data(), nAtoms, 3);
    system->getGeometry()->setCoordinates(coordinates);
    if (print) {
      system->getGeometry()->printToFile(system->getHDF5BaseName(), system->getSystemIdentifier());
      system->getGeometry()->updateTrajFile(system->getHDF5BaseName(), value, gradients.norm());
    }
    ScfTask<Options::SCF_MODES::RESTRICTED> scfTask(system);
    scfTask.run();

    GradientTask<Options::SCF_MODES::RESTRICTED> gradientTask({system});
    gradientTask.settings.print = false;
    gradientTask.run();

    auto tmp = system->getGeometry()->getGradients();
    gradients = Eigen::Map<const Eigen::VectorXd>(tmp.data(), nAtoms * 3);
    value = system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  };
  // prepare output and save engines
  const bool io1 = iOOptions.printSCFCycleInfo;
  const bool io2 = iOOptions.printSCFResults;
  iOOptions.printSCFCycleInfo = false;
  iOOptions.printSCFResults = false;
  iOOptions.printGridInfo = false;
  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 4);
  /* o======================o
   * |  OPTIMIZATION MODES  |
   * o======================o */

  /* ========================
   *   Bofill Using Hessian
   * ======================== */
  if (_env.size() == 0) {
    std::string hessianpath(_ts->getSystemPath() + _ts->getSystemName() + ".hess.h5");
    struct stat buffer;
    if (stat(hessianpath.c_str(), &buffer) != 0) {
      if (_ts->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        HessianTask<Options::SCF_MODES::RESTRICTED> hessianTask({_ts});
        hessianTask.run();
      }
      else {
        HessianTask<Options::SCF_MODES::UNRESTRICTED> hessianTask({_ts});
        hessianTask.run();
      }
      HDF5::Filepath path(hessianpath);
    }

    /*
     * Load normal modes
     */
    HDF5::H5File file(hessianpath.c_str(), H5F_ACC_RDONLY);
    HDF5::dataset_exists(file, "normalModes");

    unsigned int nAtoms = _ts->getGeometry()->getNAtoms();
    unsigned int nModes = 3 * nAtoms - (_ts->getGeometry()->isLinear() ? 5 : 6);
    Eigen::MatrixXd normalModes(3 * nAtoms, nModes);
    HDF5::load(file, "normalModes", normalModes);
    file.close();
    std::unique_ptr<Eigen::VectorXd> guess;
    std::unique_ptr<Eigen::VectorXd> mode(new Eigen::VectorXd(normalModes.col(settings.normalmode - 1)));
    normalModes.resize(0, 0);

    /*
     * Run TS search
     */
    auto ts = _ts->getGeometry()->getCoordinates();
    guess.reset(new Eigen::VectorXd(Eigen::Map<const Eigen::VectorXd>(ts.data(), 3 * nAtoms)));
    Bofill bofill(*guess, 0.3, std::move(mode));
    bofill.optimize(getData);
    /* =======================
     *   LST/QST then Bofill
     * ======================= */
  }
  else if (_env.size() == 2) {
    // sanity checks
    unsigned int nAtoms = _ts->getGeometry()->getNAtoms();
    assert(nAtoms == _env[0]->getGeometry()->getNAtoms() && "The systems have different numbers of atoms.");
    assert(nAtoms == _env[1]->getGeometry()->getNAtoms() && "The systems have different numbers of atoms.");
    // get coordinates
    auto ts = _ts->getGeometry()->getAlignedCoordinates();
    auto coords1 = _env[0]->getGeometry()->getAlignedCoordinates();
    auto coords2 = _env[1]->getGeometry()->getAlignedCoordinates();

    // each environment is one minimum
    Eigen::VectorXd min1 = Eigen::Map<const Eigen::VectorXd>(coords1.data(), 3 * nAtoms);
    Eigen::VectorXd min2 = Eigen::Map<const Eigen::VectorXd>(coords2.data(), 3 * nAtoms);

    std::unique_ptr<Eigen::VectorXd> guess;
    if (!settings.lst) {
      guess.reset(new Eigen::VectorXd(Eigen::Map<const Eigen::VectorXd>(ts.data(), 3 * nAtoms)));
    }
    // Run LST/QST and gather final tangent
    std::unique_ptr<Eigen::VectorXd> tangent;
    if (!settings.lst) {
      printSubSectionTitle("QST");
      QST lstqst(min1, min2, !settings.lstqstonly, std::move(guess));
      lstqst.optimize(getData);
      tangent.reset(new Eigen::VectorXd(lstqst.getTangent()));
    }
    else {
      printSubSectionTitle("LST/QST");
      QST qst(min1, min2, !settings.lstqstonly);
      qst.optimize(getData);
      tangent.reset(new Eigen::VectorXd(qst.getTangent()));
    }
    // Run dimer afterwards
    if (!settings.lstqstonly) {
      printSubSectionTitle("Bofill");
      ts = _ts->getGeometry()->getCoordinates();
      guess.reset(new Eigen::VectorXd(Eigen::Map<const Eigen::VectorXd>(ts.data(), 3 * nAtoms)));
      Bofill bofill(*guess, 0.3, std::move(tangent));
      bofill.optimize(getData);
    }
  }
  else {
    assert(false && "The amount of systems given does not fit any TS optimization mode.");
  }
  // Free engines and reset output
  iOOptions.printSCFCycleInfo = io1;
  iOOptions.printSCFResults = io2;
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 4);
};
} /* namespace Serenity */
