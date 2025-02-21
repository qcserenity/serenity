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
#include "energies/EnergyContributions.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h"
#include "io/IOOptions.h"
#include "math/optimizer/BFGS.h"
#include "math/optimizer/Optimizer.h"
#include "math/optimizer/SQNM.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/GradientTask.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <ctime>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GeometryOptimizationTask<SCFMode>::GeometryOptimizationTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                                                            const std::vector<std::shared_ptr<SystemController>>& passiveSystems)
  : _activeSystems(activeSystems), _passiveSystems(passiveSystems) {
}

template<Options::SCF_MODES SCFMode>
void GeometryOptimizationTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _activeSystems);
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
      _minimumPrint = true;
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
    std::shared_ptr<Optimizer> optimizer;
    switch (settings.optAlgorithm) {
      case Options::OPTIMIZATION_ALGORITHMS::BFGS:
        optimizer = std::make_shared<BFGS>(coords, 1.0, true);
        break;
      case Options::OPTIMIZATION_ALGORITHMS::SQNM:
        optimizer = std::make_shared<SQNM>(coords, settings.sqnm.historyLength, settings.sqnm.epsilon,
                                           settings.sqnm.alpha, settings.sqnm.energyThreshold, settings.sqnm.trustRadius);
        break;
    }
    /*
     * Optimization Cycle
     */
    unsigned int cycle = 1;
    Eigen::VectorXd oldParams = coords;
    double oldEnergy = std::numeric_limits<double>::infinity();
    /*
     * Lambda Function for the Optimizer
     */
    auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                    std::shared_ptr<Eigen::MatrixXd> hessian, bool printInfo) {
      // Print header before the first cycle if minimal output is requested
      if (_minimumPrint && cycle == 1) {
        printf("%7s%17s%17s%17s%17s%17s%17s%17s%17s\n", "cycle", "Total Energy", "Energy Change", "Gradient Norm",
               "RMS Gradient", "Max Gradient", "RMS Step", "Max step", "s/cycle");
      }
      // get timings for minimal output
      if (_minimumPrint) {
        clock_gettime(CLOCK_REALTIME, &_time);
      }
      if (printInfo && !_minimumPrint) {
        printSubSectionTitle("Cycle: " + std::to_string(cycle));
      }
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
      GradientTask<SCFMode> task(_activeSystems);
      task.settings.transInvar = settings.transInvar;
      task.settings.gradType = settings.gradType;
      task.settings.numGradStepSize = settings.numGradStepSize;
      task.settings.print = !_minimumPrint;
      task.run();

      gradients = Eigen::Map<const Eigen::VectorXd>(_activeSystems[0]->getGeometry()->getGradients().data(), nAtoms * 3);

      // Calculate convergence parameters
      double RMSgrad = std::sqrt(gradients.squaredNorm() / gradients.size());
      double gradNorm = gradients.norm();
      auto step = (parameters - oldParams);
      /*
       * Output and increment
       */
      if (printInfo && !_minimumPrint) {
        printSmallCaption("Geometry Relaxation");
        OutputControl::n.printf("%4s Total Energy  %15.10f\n", "", value);
        OutputControl::n.printf("%4s Energy Change %15.10f\n", "", value - oldEnergy);
        OutputControl::n.printf("%4s Gradient Norm %15.10f\n", "", gradNorm);
        OutputControl::n.printf("%4s RMS Gradient  %15.10f\n", "", RMSgrad);
        OutputControl::n.printf("%4s Max Gradient  %15.10f\n", "", gradients.array().abs().maxCoeff());
        OutputControl::n.printf("%4s RMS Step      %15.10f\n", "", std::sqrt(step.squaredNorm() / step.size()));
        OutputControl::n.printf("%4s Max Step      %15.10f\n\n", "", step.array().abs().maxCoeff());
        _activeSystems[0]->getGeometry()->print();
      }

      if (_minimumPrint) {
        timespec now;
        clock_gettime(CLOCK_REALTIME, &now);
        double sec = (double)(now.tv_sec - _time.tv_sec) + (now.tv_nsec - _time.tv_nsec) * 0.000000001;
        printf("%7d%17f%17f%17f%17f%17f%17f%17f%17f\n", cycle, value, value - oldEnergy, gradNorm, RMSgrad,
               gradients.array().abs().maxCoeff(), std::sqrt(step.squaredNorm() / step.size()),
               step.array().abs().maxCoeff(), sec);
      }

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
        converged = true;
      }
      if (cycle == settings.maxCycles) {
        if (printInfo)
          WarningTracker::printWarning("WARNING: Geometry Optimization not yet converged!", true);
        converged = true;
      }

      oldEnergy = value;
      oldParams = parameters;

      if (converged) {
        OutputControl::m.printf("Convergence reached after %3u cycles.\n", cycle);
      }

      ++cycle;

      return converged;
    };
    std::shared_ptr<unsigned int> nRejected = std::make_shared<unsigned int>(0);
    optimizer->optimize(updateFunction, nRejected);

    /*
     * Print Final Results
     */
    printSubSectionTitle("Final Results");
    _activeSystems[0]->getGeometry()->print();
    OutputControl::m.printf("%4s Final Energy: %15.10f\n\n\n", "",
                            _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy());
    OutputControl::m.printf("Convergence reached after %3u cycles.\n", cycle - 1);
    if (settings.optAlgorithm == Options::OPTIMIZATION_ALGORITHMS::SQNM) {
      OutputControl::m.printf("A total of %3u cycles have been rejected.\n", *(nRejected));
    }
  }

  else {
    Geometry supersystemGeometry;

    for (auto sys : _activeSystems) {
      supersystemGeometry += *(sys->getGeometry());
    }
    for (auto sys : _passiveSystems) {
      supersystemGeometry += *(sys->getGeometry());
    }
    const unsigned int nAtoms = supersystemGeometry.getNAtoms();
    Eigen::VectorXd tmpCoords = Eigen::Map<const Eigen::VectorXd>(supersystemGeometry.getCoordinates().data(), nAtoms * 3);
    /*
     * Get the optimizer
     */
    std::shared_ptr<Optimizer> optimizer;
    switch (settings.optAlgorithm) {
      case Options::OPTIMIZATION_ALGORITHMS::BFGS:
        optimizer = std::make_shared<BFGS>(tmpCoords, 1.0, true);
        break;
      case Options::OPTIMIZATION_ALGORITHMS::SQNM:
        optimizer = std::make_shared<SQNM>(tmpCoords, settings.sqnm.historyLength, settings.sqnm.epsilon,
                                           settings.sqnm.alpha, settings.sqnm.energyThreshold, settings.sqnm.trustRadius);
        break;
    }
    /*
     * Optimization cycle
     */
    unsigned int optCycle = 1;
    Eigen::VectorXd oldParams = tmpCoords;
    double oldEnergy = std::numeric_limits<double>::infinity();
    /*
     * Lambda Function for the Optimizer
     */
    auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                    std::shared_ptr<Eigen::MatrixXd> hessian, bool printInfo) {
      // Print header before the first cycle if minimal output is requested
      if (_minimumPrint && optCycle == 1) {
        printf("%7s%17s%17s%17s%17s%17s%17s%17s%17s\n", "cycle", "Total Energy", "Energy Change", "Gradient Norm",
               "RMS Gradient", "Max Gradient", "RMS Step", "Max step", "s/cycle");
      }
      // get timings for minimal output
      if (_minimumPrint) {
        clock_gettime(CLOCK_REALTIME, &_time);
      }
      if (printInfo && !_minimumPrint) {
        printSubSectionTitle("Cycle: " + std::to_string(optCycle));
      }
      (void)hessian;
      // Atom Coordinates (the parameters to optimize)
      Matrix<double> coordinates = Eigen::Map<const Eigen::MatrixXd>(parameters.data(), nAtoms, 3);
      // Set coordinates of all atoms to new parameter data from previous optimization cycle
      supersystemGeometry.setCoordinates(coordinates);

      // Calc/Print Gradients - in a subsystem case the gradienttask also always performes a FaTTask, thus the energy
      // calculation is outsourced

      GradientTask<SCFMode> gradientTask(_activeSystems, _passiveSystems);
      gradientTask.settings.gradType = settings.gradType;
      gradientTask.settings.embedding = settings.embedding;
      gradientTask.settings.FDEgridCutOff = settings.FaTgridCutOff;
      gradientTask.settings.transInvar = settings.transInvar;
      gradientTask.settings.numGradStepSize = settings.numGradStepSize;
      gradientTask.settings.print = !_minimumPrint;
      gradientTask.run();

      gradients =
          Eigen::Map<const Eigen::VectorXd>(supersystemGeometry.getGradients().data(), supersystemGeometry.getNAtoms() * 3);

      switch (_activeSystems[0]->getSettings().method) {
        case Options::ELECTRONIC_STRUCTURE_THEORIES::DFT:
          value = _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy(
              ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT);
          break;
        case Options::ELECTRONIC_STRUCTURE_THEORIES::HF:
          value = _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy(
              ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT);
          break;
      }

      // Calculate convergence parameters
      double gradNorm = gradients.norm();
      double RMSgrad = std::sqrt(gradients.squaredNorm() / gradients.size());
      auto step = parameters - oldParams;
      /*
       * Output and increment
       */
      if (printInfo && !_minimumPrint) {
        printSmallCaption("Geometry Relaxation");
        OutputControl::n.printf("%4s Total Energy  %15.10f\n", "", value);
        OutputControl::n.printf("%4s Energy Change %15.10f\n", "", value - oldEnergy);
        OutputControl::n.printf("%4s Gradient Norm %15.10f\n", "", gradNorm);
        OutputControl::n.printf("%4s RMS Gradient  %15.10f\n", "", RMSgrad);
        OutputControl::n.printf("%4s Max Gradient  %15.10f\n", "", gradients.array().abs().maxCoeff());
        OutputControl::n.printf("%4s RMS Step      %15.10f\n", "", std::sqrt(step.squaredNorm() / step.size()));
        OutputControl::n.printf("%4s Max Step      %15.10f\n\n", "", step.array().abs().maxCoeff());
        supersystemGeometry.print();
      }
      if (_minimumPrint) {
        timespec now;
        clock_gettime(CLOCK_REALTIME, &now);
        double sec = (double)(now.tv_sec - _time.tv_sec) + (now.tv_nsec - _time.tv_nsec) * 0.000000001;
        printf("%7d%17f%17f%17f%17f%17f%17f%17f%17f\n", optCycle, value, value - oldEnergy, gradNorm, RMSgrad,
               gradients.array().abs().maxCoeff(), std::sqrt(step.squaredNorm() / step.size()),
               step.array().abs().maxCoeff(), sec);
      }
      supersystemGeometry.printToFile("./opt", "");
      supersystemGeometry.updateTrajFile("./opt", value, gradNorm);
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
        converged = true;
      }
      if (optCycle == settings.maxCycles) {
        WarningTracker::printWarning("WARNING: Geometry Optimization not yet converged!", true);
        converged = true;
      }

      oldEnergy = value;
      oldParams = parameters;

      if (converged) {
        printf("Convergence reached after %3u cycles.\n", optCycle);
      }

      ++optCycle;

      return converged;
    };

    std::shared_ptr<unsigned int> nRejected = std::make_shared<unsigned int>(0);
    optimizer->optimize(updateFunction, nRejected);

    printSubSectionTitle("Final Results");
    supersystemGeometry.print();
    switch (_activeSystems[0]->getSettings().method) {
      case Options::ELECTRONIC_STRUCTURE_THEORIES::DFT:
        OutputControl::m.printf("%4s Final Energy: %15.10f\n\n\n", "",
                                _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy(
                                    ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT));
        break;
      case Options::ELECTRONIC_STRUCTURE_THEORIES::HF:
        OutputControl::m.printf("%4s Final Energy: %15.10f\n\n\n", "",
                                _activeSystems[0]->template getElectronicStructure<SCFMode>()->getEnergy(
                                    ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT));
        break;
    }
    OutputControl::m.printf("Convergence reached after %3u cycles.\n", optCycle - 1);
    if (settings.optAlgorithm == Options::OPTIMIZATION_ALGORITHMS::SQNM) {
      OutputControl::m.printf("A total of %3u cycles have been rejected.\n", *(nRejected));
    }
  }

  iOOptions.printSCFCycleInfo = info;
  iOOptions.printSCFResults = results;
  iOOptions.gridAccuracyCheck = gridAcc;
}

template class GeometryOptimizationTask<Options::SCF_MODES::RESTRICTED>;
template class GeometryOptimizationTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
