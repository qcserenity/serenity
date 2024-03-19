/**
 * @file   GeometryOptimizationTask.cpp
 *
 * @date   May 21, 2015
 * @author Jan Unsleber
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
/* Include Class Header*/
#include "tasks/GeometryOptimizationTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "energies/EnergyContributions.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "geometry/gradients/GeometryGradientCalculator.h"
#include "geometry/gradients/NumericalGeomGradCalc.h"
#include "io/FormattedOutput.h"
#include "io/IOOptions.h"
#include "math/optimizer/BFGS.h"
#include "math/optimizer/LBFGS.h"
#include "math/optimizer/Optimizer.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/GradientTask.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <cmath>
#include <iostream>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GeometryOptimizationTask<SCFMode>::GeometryOptimizationTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                                                            const std::vector<std::shared_ptr<SystemController>>& passiveSystems)
  : _activeSystems(activeSystems), _passiveSystems(passiveSystems) {
  assert(_activeSystems[0]);
}

template<Options::SCF_MODES SCFMode>
void GeometryOptimizationTask<SCFMode>::run() {
  /*
   * Output Options
   */
  bool info(iOOptions.printSCFCycleInfo);
  bool results(iOOptions.printSCFResults);
  bool gridAcc(iOOptions.gridAccuracyCheck);
  switch (GLOBAL_PRINT_LEVEL) {
    case Options::GLOBAL_PRINT_LEVELS::MINIMUM:
      iOOptions.printSCFCycleInfo = false;
      iOOptions.printSCFResults = false;
      iOOptions.printGridInfo = false;
      break;
    default:
      iOOptions.printSCFCycleInfo = true;
      iOOptions.printSCFResults = true;
      break;
  }
  iOOptions.gridAccuracyCheck = false;

  if (_passiveSystems.size() == 0 and _activeSystems.size() == 1) {
    auto nAtoms = _activeSystems[0]->getGeometry()->getNAtoms();
    Eigen::VectorXd coords(nAtoms * 3);
    coords = Eigen::Map<const Eigen::VectorXd>(_activeSystems[0]->getGeometry()->getCoordinates().data(), nAtoms * 3);
    /*
     * Get the optimizer
     */
    std::shared_ptr<Optimizer> optimizer = std::make_shared<BFGS>(coords, 1.0, true);
    /*
     * Optimization Cycle
     */
    unsigned int cycle = 1;
    Eigen::VectorXd oldParams = coords;
    double oldEnergy = _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy();
    /*
     * Lambda Function for the Optimizer
     */
    auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                    std::shared_ptr<Eigen::MatrixXd> hessian, bool printInfo) {
      (void)hessian;
      // Atom Coordinates (the parameters to optimize)
      Matrix<double> coordinates = Eigen::Map<const Eigen::MatrixXd>(parameters.data(), nAtoms, 3);
      // Set Coordinates of all atoms to new paramater data from previous optimization cycle
      _activeSystems[0]->getGeometry()->setCoordinates(coordinates);

      // Calculate Energy for new Coordinate Set
      auto scf = ScfTask<SCFMode>(_activeSystems[0]);
      scf.run();
      value = _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy();

      // Calc/Print Gradients

      if (printInfo)
        printSubSectionTitle("Cycle: " + std::to_string(cycle));
      GradientTask<SCFMode> task(_activeSystems);
      task.settings.transInvar = settings.transInvar;
      task.settings.gradType = settings.gradType;
      task.settings.numGradStepSize = settings.numGradStepSize;
      task.run();

      gradients = Eigen::Map<const Eigen::VectorXd>(_activeSystems[0]->getGeometry()->getGradients().data(), nAtoms * 3);

      // Calculate convergence parameters
      double RMSgrad = std::sqrt(gradients.squaredNorm() / gradients.size());
      double gradNorm = gradients.norm();
      auto step = (parameters - oldParams);
      /*
       * Output and increment
       */
      if (printInfo)
        printSmallCaption("Geometry Relaxation");
      printf("%4s Total Energy  %15.10f\n", "", value);
      printf("%4s Energy Change %15.10f\n", "", value - oldEnergy);
      printf("%4s Gradient Norm %15.10f\n", "", gradNorm);
      printf("%4s RMS Gradient  %15.10f\n", "", RMSgrad);
      printf("%4s Max Gradient  %15.10f\n", "", gradients.array().abs().maxCoeff());
      printf("%4s RMS Step      %15.10f\n", "", std::sqrt(step.squaredNorm() / step.size()));
      printf("%4s Max Step      %15.10f\n\n", "", step.array().abs().maxCoeff());
      _activeSystems[0]->getGeometry()->print();
      _activeSystems[0]->getGeometry()->printToFile(_activeSystems[0]->getHDF5BaseName(),
                                                    _activeSystems[0]->getSystemIdentifier());
      _activeSystems[0]->getGeometry()->updateTrajFile(_activeSystems[0]->getHDF5BaseName(), value, gradNorm);

      /*
       * Convergence check
       */
      bool converged = false;
      unsigned int convCriteriaMet = 0;
      if (abs(RMSgrad) < settings.rmsgradThresh)
        convCriteriaMet++;
      if (abs(value - oldEnergy) < settings.energyChangeThresh)
        convCriteriaMet++;
      if (gradients.array().abs().maxCoeff() < settings.maxGradThresh)
        convCriteriaMet++;
      if (std::sqrt(step.squaredNorm() / step.size()) < settings.stepThresh)
        convCriteriaMet++;
      if (step.array().abs().maxCoeff() < settings.maxStepThresh)
        convCriteriaMet++;
      if (convCriteriaMet >= 3 and cycle > 1) {
        print("Convergence reached. Exiting...");
        converged = true;
      }
      if (cycle == settings.maxCycles) {
        if (printInfo)
          WarningTracker::printWarning("WARNING: Geometry Optimization not yet converged!", true);
        converged = true;
      }

      oldEnergy = value;
      oldParams = parameters;

      ++cycle;

      return converged;
    };

    optimizer->optimize(updateFunction);
  }
  else {
    /*
     * Initial FaT
     */
    FreezeAndThawTask<SCFMode> task(_activeSystems, _passiveSystems);
    task.settings.embedding = settings.embedding;
    if (settings.embedding.embeddingMode != Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
      throw SerenityError("Only non-additive kinetic embedding is supported for embedded Gradients!");
    }
    task.settings.gridCutOff = settings.FaTgridCutOff;
    task.run();

    Geometry SupersystemGeometry;

    for (auto sys : _activeSystems) {
      SupersystemGeometry += *(sys->getGeometry());
    }
    for (auto sys : _passiveSystems) {
      SupersystemGeometry += *(sys->getGeometry());
    }

    const unsigned int nAtoms = SupersystemGeometry.getNAtoms();
    Eigen::VectorXd tmpCoords = Eigen::Map<const Eigen::VectorXd>(SupersystemGeometry.getCoordinates().data(), nAtoms * 3);
    /*
     * Get the optimizer
     */
    std::shared_ptr<Optimizer> optimizer = std::make_shared<BFGS>(tmpCoords, 1.0, true);

    auto oldParams = tmpCoords;
    assert(_activeSystems[0]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT &&
           "Check logic in this task for its compatible with HF.");
    double oldEnergy = _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy(
        ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT);
    unsigned int optCycle = 1;

    auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                    std::shared_ptr<Eigen::MatrixXd> hessian, bool printInfo) {
      (void)hessian;
      if (printInfo)
        printSubSectionTitle("Cycle: " + std::to_string(optCycle));

      Matrix<double> coordinates = Eigen::Map<const Eigen::MatrixXd>(parameters.data(), nAtoms, 3);
      SupersystemGeometry.setCoordinates(coordinates);

      // Calc/Print Gradients

      GradientTask<SCFMode> task(_activeSystems, _passiveSystems);
      task.settings.gradType = settings.gradType;
      task.settings.embedding = settings.embedding;
      if (settings.embedding.embeddingMode != Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
        throw SerenityError("Only non-additive kinetic embedding is supported for embedded Gradients!");
      }
      task.settings.FDEgridCutOff = settings.FaTgridCutOff;
      task.settings.transInvar = settings.transInvar;
      task.settings.numGradStepSize = settings.numGradStepSize;
      task.run();

      /* Variables for optimization. */
      value = _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy(
          ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT);

      gradients =
          Eigen::Map<const Eigen::VectorXd>(SupersystemGeometry.getGradients().data(), SupersystemGeometry.getNAtoms() * 3);

      double gradNorm = gradients.norm();

      double RMSgrad = std::sqrt(gradients.squaredNorm() / gradients.size());

      auto step = parameters - oldParams;

      /*
       * Output and increment
       */
      printSmallCaption("Geometry Relaxation");
      printf("%4s Total Energy  %15.10f\n", "", value);
      printf("%4s Energy Change %15.10f\n", "", value - oldEnergy);
      printf("%4s Gradient Norm %15.10f\n", "", gradNorm);
      printf("%4s RMS Gradient  %15.10f\n", "", RMSgrad);
      printf("%4s Max Gradient  %15.10f\n", "", gradients.array().abs().maxCoeff());
      printf("%4s RMS Step      %15.10f\n", "", std::sqrt(step.squaredNorm() / step.size()));
      printf("%4s Max Step      %15.10f\n\n", "", step.array().abs().maxCoeff());
      SupersystemGeometry.print();
      SupersystemGeometry.printToFile("./opt", "");
      SupersystemGeometry.updateTrajFile("./opt", value, gradNorm);
      for (auto sys : _activeSystems) {
        sys->getGeometry()->printToFile(sys->getHDF5BaseName(), sys->getSystemIdentifier());
      }

      /*
       * Convergence check
       */
      bool converged = false;
      unsigned int convCriteriaMet = 0;
      if (abs(RMSgrad) < settings.rmsgradThresh)
        convCriteriaMet++;
      if (abs(value - oldEnergy) < settings.energyChangeThresh)
        convCriteriaMet++;
      if (gradients.array().abs().maxCoeff() < settings.maxGradThresh)
        convCriteriaMet++;
      if (std::sqrt(step.squaredNorm() / step.size()) < settings.stepThresh)
        convCriteriaMet++;
      if (step.array().abs().maxCoeff() < settings.maxStepThresh)
        convCriteriaMet++;
      if (convCriteriaMet >= 3 and optCycle > 1) {
        print("Convergence reached. Exiting...");
        converged = true;
      }
      if (optCycle == settings.maxCycles) {
        printSubSectionTitle("WARNING: Geometry Optimization not yet converged!");
        converged = true;
      }

      oldEnergy = value;
      oldParams = parameters;

      ++optCycle;

      return converged;
    };

    optimizer->optimize(updateFunction);
  }

  iOOptions.printSCFCycleInfo = info;
  iOOptions.printSCFResults = results;
  iOOptions.gridAccuracyCheck = gridAcc;
}

template class GeometryOptimizationTask<Options::SCF_MODES::RESTRICTED>;
template class GeometryOptimizationTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
