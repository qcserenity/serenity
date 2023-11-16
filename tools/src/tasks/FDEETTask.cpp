/**
 * @file   FDEETTask.cpp
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

/* Include Class Header*/
#include "tasks/FDEETTask.h" // class header
/* Include Serenity Internal Headers */
#include "dft/Functional.h"              // check the functional
#include "geometry/Geometry.h"           //
#include "io/FormattedOutputStream.h"    // formatted output
#include "misc/WarningTracker.h"         // warnings
#include "postHF/ET/FDEDiabController.h" // acces fde-diab functionality
#include "postHF/ET/FDEETCalculator.h"   // calcualte and store objects
#include "postHF/ET/FDEETController.h"   // calculate and store diabatic state info
#include "settings/Settings.h"           // settings
#include "system/SystemController.h"     // setting up supersystem
#include "tasks/SystemAdditionTask.h"    // setting up supersystem

namespace Serenity {

FDEETTask::FDEETTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems) : _act(activeSystems) {
}

void FDEETTask::run() {
  printSectionTitle("FDE-ET/diab");
  if (settings.disjoint.size() > 0) {
    printSubSectionTitle("Disjoint Approximation");
    OutputControl::nOut << "  Systems";
    for (auto index : settings.disjoint)
      OutputControl::nOut << "  " << index << "  ";
    OutputControl::nOut << "are disjoint" << std::endl << std::endl;
  }
  // Sanity checks
  if (resolveFunctional(settings.functional).isHybrid())
    throw SerenityError("Hybrid functionals cannot be used so far!");
  // check for active systems and states consisitency
  unsigned sum = 0;
  for (auto& item : settings.states)
    sum += item;
  if (sum > _act.size())
    throw SerenityError("You provided less active systems than specified in the states keyword!");
  if (sum < _act.size())
    WarningTracker::printWarning(
        "WARNING: You provided more active systems than specified in the states keyword. They will be ignored!", true);
  // check if all systems are unrestricted
  for (auto sys : _act) {
    if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) // keep in mind: if system carries odd spin
                                                                      // serenity uses unrestricted anyway
      throw SerenityError("Cannot perform FDE-ET for spin restricted SCFMode!");
  }
  // Get all information of systems from input and parse to controller
  _fdeetController = std::make_shared<FDEETController>(_act, settings.states.size(), settings.states[0]);
  // Large loop over the whole formalism coupling different states from the input
  // Perform calculations for each coupling step
  for (unsigned iCoupling = 0; iCoupling < settings.couple.size(); ++iCoupling) {
    unsigned nStates = settings.couple[iCoupling].size();
    if (nStates < 2)
      throw SerenityError("You must at least couple two states in FDE-ET/diab!");
    // check if the the populations vector is set up accordingly
    // sort entries in vector and get rid of duplicates
    std::sort(settings.population.begin(), settings.population.end());
    settings.population.erase(std::unique(settings.population.begin(), settings.population.end()), settings.population.end());
    for (auto& item : settings.population) {
      if (item > nStates - 1)
        throw SerenityError(
            "You want to calculate spin populations for an adiabatic state which will not be calculated!");
    }
    // get from controller which states shall be coupled in this step -> info stored in stateIndex
    _fdeetController->setNStatesCouple(nStates);
    _fdeetController->setupIndexVector(settings.couple[iCoupling]);
    auto stateIndex = _fdeetController->getIndexVector();
    OutputControl::nOut << "  Coupling following states from your input:" << std::endl;
    std::string name = "superSystem";
    if (settings.disjoint.size() > 0)
      name += "-disjoint-";
    for (auto index : (*stateIndex)) {
      OutputControl::nOut << "  " << index + 1 << "  ";
      name += std::to_string(index + 1);
    }
    OutputControl::nOut << std::endl;
    OutputControl::nOut << std::endl;
    // check if the same subsystems are contained in each state, and if electronic structure is already there
    for (unsigned iState = 1; iState < nStates; ++iState) {
      for (unsigned iSys = 0; iSys < _fdeetController->getNSysPerState(); ++iSys) {
        std::shared_ptr<Geometry> lhsGeo = _fdeetController->getSubsystems()[iSys]->getGeometry();
        std::shared_ptr<Geometry> rhsGeo =
            _fdeetController->getSubsystems()[iSys + _fdeetController->getNSysPerState() * iState]->getGeometry();
        if (*lhsGeo != *rhsGeo) {
          OutputControl::nOut << "Mismatch of the system's geometries in states  1 and " << iState + 1 << std::endl;
          OutputControl::nOut << "System " << _fdeetController->getSubsystems()[iSys]->getSystemName() << std::endl;
          lhsGeo->print();
          OutputControl::nOut
              << "System "
              << _fdeetController->getSubsystems()[iSys + _fdeetController->getNSysPerState() * iState]->getSystemName()
              << std::endl;
          rhsGeo->print();
          throw SerenityError("All quasi-diabatic states should contain the same subsystems within the same order!");
        }
        if (!_fdeetController->getSubsystems()[iSys]->template hasElectronicStructure<Options::SCF_MODES::UNRESTRICTED>())
          WarningTracker::printWarning(
              "WARNING: System " + _fdeetController->getSubsystems()[iSys]->getSystemName() +
                  " has no initial electronic structure!!! Usually a prior FaT calculation needs to be performed.",
              true);
      } /* iSys */
    }   /* iState */
    // get the dimensions of occupied orbitals and basis functions
    _fdeetController->setupOrbitalInfo();
    // check for empty spin orbitals, if needed crash
    if ((_fdeetController->getNOccState().array() == 0).any())
      throw SerenityError("The alpha and beta orbitals should have at least one electron, each!");
    // info gathering for coupling step is complete, now run the actual calculation
    // and then stores all necessary objects
    _fdeetCalculator = std::make_shared<FDEETCalculator>();
    // Calc MO Transition Overlap Matrix
    OutputControl::nOut << "Calculating Transition Overlap Matrix...\n" << std::endl;
    _fdeetCalculator->calculateMoOverlap(nStates, _fdeetController->getNSysPerState(), settings.disjoint,
                                         _fdeetController->getStateCoeffs(), _fdeetController->getStateVector(),
                                         _fdeetController->getIndexVector(), _fdeetController->getNOccState());
    // Calc Moore Penrose Pseudo Inverse
    OutputControl::nOut << "Calculating Pseudo Inverse of Transition Overlap Matrix...\n\n" << std::endl;
    _fdeetCalculator->calcBlockedMoorePenrose(_fdeetCalculator->getMoOverlap(), _fdeetController->getNOccState(),
                                              nStates, settings.invThreshold);
    // Building an artifical supersystem we construct adiabatic wave functions for
    OutputControl::nOut << "Constructing Artificial Supersystem...\n" << std::endl;
    std::vector<std::shared_ptr<SystemController>> coupleSystems;
    for (unsigned i = 0; i < _fdeetController->getNSysPerState(); ++i)
      coupleSystems.push_back(_act[i]);
    std::shared_ptr<SystemController> superSystem = setUpSupersystem(coupleSystems, name);
    superSystem->getGeometry()->print();
    // Calc Transition Density Matrices
    OutputControl::nOut << "Calculating Transition Density Matrices...\n" << std::endl;
    if (settings.diskMode) {
      _fdeetCalculator->calcTransDensMatsDisk(_fdeetCalculator->getPseudoInverse(), superSystem,
                                              _fdeetController->getStateCoeffs(), _fdeetController->getNOccState(), nStates);
    }
    else {
      _fdeetCalculator->calcTransDensMats(_fdeetCalculator->getPseudoInverse(), superSystem,
                                          _fdeetController->getStateCoeffs(), _fdeetController->getNOccState(), nStates);
    }

    if (settings.printContributions)
      _fdeetCalculator->printTransDensContributions(superSystem, nStates);
    // Calc Hamilton Matrix
    if (settings.diskMode) {
      OutputControl::nOut << "Calculating Hamilton Matrix...\n" << std::endl;
      _fdeetCalculator->calcHamiltonianDisk(superSystem, nStates);
    }
    else {
      OutputControl::nOut << "Calculating Hamilton Matrix...\n" << std::endl;
      _fdeetCalculator->calcHamiltonian(_fdeetCalculator->getTransDensMats(), superSystem, nStates);
    }
    // Solve Generalized Eigenvalue Problem HC = SCE
    OutputControl::nOut << "\nSolving Generalized Eigenvalue Problem...\n" << std::endl;
    _fdeetCalculator->solveEigenValueProblem(_fdeetCalculator->getHamiltonian(), _fdeetCalculator->getDeterminants(), nStates);
    // Print FDE-ET Results
    _fdeetCalculator->printResults(nStates);
    _fdeDiabController = nullptr;
    // Calculate and print adiabatic spin densities to cube file
    if (settings.spindensity == true or settings.spinpopulation == true) {
      if (settings.diskMode) {
        _fdeDiabController =
            std::make_shared<FDEDiabController>(_fdeetCalculator->getDensityMatrixFiles(), _fdeetCalculator->getLinCoeffs(),
                                                _fdeetCalculator->getDeterminants(), superSystem, nStates);
      }
      else {
        _fdeDiabController =
            std::make_shared<FDEDiabController>(_fdeetCalculator->getTransDensMats(), _fdeetCalculator->getLinCoeffs(),
                                                _fdeetCalculator->getDeterminants(), superSystem, nStates);
      }
      if (settings.spindensity == true) {
        _fdeDiabController->calcAdiabDensities();
        _fdeDiabController->printAdiabDensities();
      }
      if (settings.spinpopulation == true) {
        _fdeDiabController->calcPopulations(settings.population);
        _fdeDiabController->printPopulations(settings.population);
        _fdeDiabController->adiabaticWavefunctionAnalysis(settings.population);
      }
    }
  } /* End of coupling */
} /* End of task run */

std::shared_ptr<SystemController> FDEETTask::setUpSupersystem(std::vector<std::shared_ptr<SystemController>> coupleSystems,
                                                              std::string name) {
  Settings superSysSettings = _act[0]->getSettings();
  superSysSettings.spin = 0;
  superSysSettings.charge = 0;
  superSysSettings.name = name.c_str();
  superSysSettings.path = "sys-";
  superSysSettings.dft.functional = settings.functional;
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), superSysSettings);
  SystemAdditionTask<Options::SCF_MODES::UNRESTRICTED> additionTask(supersystem, coupleSystems);
  additionTask.settings.addOccupiedOrbitals = false;
  additionTask.run();
  return supersystem;
}

std::shared_ptr<FDEETController> FDEETTask::getFDEETController() {
  return _fdeetController;
}

std::shared_ptr<FDEETCalculator> FDEETTask::getFDEETCalculator() {
  return _fdeetCalculator;
}

std::shared_ptr<FDEDiabController> FDEETTask::getFDEDiabController() {
  return _fdeDiabController;
}

} /* namespace Serenity */
