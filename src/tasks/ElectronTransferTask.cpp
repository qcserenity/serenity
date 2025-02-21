/**
 * @file   ElectronTransferTask.cpp
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
#include "tasks/ElectronTransferTask.h" // class header
/* Include Serenity Internal Headers */
#include "dft/Functional.h"              // check the functional
#include "geometry/Geometry.h"           // geometry
#include "io/FormattedOutputStream.h"    // formatted output
#include "misc/WarningTracker.h"         // warnings
#include "postHF/ET/FDEDiabController.h" // acces fde-diab functionality
#include "postHF/ET/FDEETCalculator.h"   // calcualte and store objects
#include "postHF/ET/FDEETController.h"   // calculate and store diabatic state info
#include "settings/Settings.h"           // settings
#include "system/SystemController.h"     // setting up supersystem
#include "tasks/SystemAdditionTask.h"    // setting up supersystem
/* Include Std and External Headers */
#include <numeric>

namespace Serenity {

ElectronTransferTask::ElectronTransferTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems)
  : _act(activeSystems) {
}

void ElectronTransferTask::run() {
  printSectionTitle("Diabatization for Electron Transfer");
  if (settings.disjoint.size() > 0) {
    printSubSectionTitle("Disjoint Approximation");
    OutputControl::nOut << "  Systems";
    for (auto index : settings.disjoint)
      OutputControl::nOut << "  " << index << "  ";
    OutputControl::nOut << "are disjoint" << std::endl << std::endl;
  }
  // Sanity checks
  if (_act[0]->getSettings().customFunc.basicFunctionals.size()
          ? Functional(_act[0]->getSettings().customFunc).isDoubleHybrid()
          : resolveFunctional(_act[0]->getSettings().dft.functional).isDoubleHybrid())
    throw SerenityError("Double hybrid functionals cannot be used so far!");
  if (_act[0]->getSettings().pcm.use)
    throw SerenityError("Solvent models cannot be used so far!");
  // check input
  if (settings.states.size() == 0)
    settings.states = std::vector<unsigned int>(_act.size(), 1);
  if (settings.couple.size() == 0) {
    std::vector<unsigned int> couple(settings.states.size());
    std::iota(couple.begin(), couple.end(), 1);
    settings.couple = {couple};
  }
  if (settings.coupleAdiabaticStates) {
    if (settings.couple.size() < 2) {
      throw SerenityError("The keyword coupledAdiabaticStates requires at least two runs.");
    }
    else {
      // set up coupling of adiabatic states
      std::vector<unsigned int> allStatesCoupled = {};
      for (auto coupledStates : settings.couple) {
        allStatesCoupled.insert(allStatesCoupled.end(), coupledStates.begin(), coupledStates.end());
      }
      std::sort(allStatesCoupled.begin(), allStatesCoupled.end());
      allStatesCoupled.erase(std::unique(allStatesCoupled.begin(), allStatesCoupled.end()), allStatesCoupled.end());
      settings.couple.push_back(allStatesCoupled);
    }
  }
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
                                                                      // Serenity uses unrestricted anyway
      throw SerenityError("Cannot perform electron transfer for spin restricted SCFMode!");
  }
  // Get all information of systems from input and parse to controller
  _fdeetController = std::make_shared<FDEETController>(_act, settings.states.size(), settings.states[0]);
  // Initialize linear coefficients used for coupling of adiabatic states
  Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(settings.couple.back().size(), settings.couple.size() - 1);
  // If configuration weights are given, the corresponding coupling steps are skipped
  if (settings.configurationWeights.size() != 0) {
    if (!settings.coupleAdiabaticStates)
      throw SerenityError("Cannot use configurationWeights if coupleAdiabaticStates is turned off!");
    if (settings.configurationWeights.size() != settings.couple.size() - 1)
      throw SerenityError("There is mismatch in the dimensions of keywords configurationWeights and couple!");
    for (unsigned iCoupling = 0; iCoupling < settings.configurationWeights.size(); ++iCoupling) {
      unsigned nStates = settings.couple[iCoupling].size();
      if (settings.configurationWeights[iCoupling].size() != nStates)
        throw SerenityError("There is mismatch in the dimensions of keywords configurationWeights and couple!");
      for (unsigned iState = 0; iState < settings.couple.back().size(); ++iState) {
        for (unsigned jState = 0; jState < nStates; ++jState) {
          if (settings.couple.back()[iState] == settings.couple[iCoupling][jState])
            coefficients(iState, iCoupling) = settings.configurationWeights[iCoupling][jState];
        }
      }
    }
    settings.couple = {settings.couple.back()};
  }
  // Large loop over the whole formalism coupling different states from the input
  // Perform calculations for each coupling step
  for (unsigned iCoupling = 0; iCoupling < settings.couple.size(); ++iCoupling) {
    unsigned nStates = settings.couple[iCoupling].size();
    if (nStates < 2)
      WarningTracker::printWarning(
          "WARNING: You provided only one state to couple. This results in an energy evaluation of this state.\n"
          "         Your results may not be meaningful!\n",
          true);
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
    if (settings.coupleAdiabaticStates && iCoupling == settings.couple.size() - 1)
      OutputControl::nOut << "  Coupling previous adiabatic states based on following states from your input:" << std::endl;
    else
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
    // Building an artifical supersystem we construct adiabatic wave functions for
    OutputControl::nOut << "Constructing Artificial Supersystem...\n" << std::endl;
    std::vector<std::shared_ptr<SystemController>> coupleSystems;
    for (unsigned i = 0; i < _fdeetController->getNSysPerState(); ++i)
      coupleSystems.push_back(_act[i]);
    std::shared_ptr<SystemController> superSystem = setUpSupersystem(coupleSystems, name);
    superSystem->getGeometry()->print();
    // Calc MO Transition Overlap Matrix
    OutputControl::nOut << "Calculating Transition Overlap Matrix...\n" << std::endl;
    _fdeetCalculator->calculateMoOverlap(nStates, _fdeetController->getNSysPerState(), settings.disjoint,
                                         _fdeetController->getStateCoeffs(), _fdeetController->getStateVector(),
                                         _fdeetController->getIndexVector(), _fdeetController->getNOccState());
    // Calc Moore Penrose Pseudo Inverse
    OutputControl::nOut << "Calculating Pseudo Inverse of Transition Overlap Matrix...\n\n" << std::endl;
    _fdeetCalculator->calcBlockedMoorePenrose(_fdeetCalculator->getMoOverlap(), _fdeetController->getNOccState(),
                                              nStates, superSystem->getBasisController()->getPrescreeningThreshold());
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
      _fdeetCalculator->calcHamiltonianDisk(superSystem, nStates, settings.useHFCoupling);
    }
    else {
      OutputControl::nOut << "Calculating Hamilton Matrix...\n" << std::endl;
      _fdeetCalculator->calcHamiltonian(_fdeetCalculator->getTransDensMats(), superSystem, nStates, settings.useHFCoupling);
    }
    // Solve Generalized Eigenvalue Problem HC = SCE and Print Results
    OutputControl::nOut << "\nSolving Generalized Eigenvalue Problem..." << std::endl;
    unsigned nAdiab = nStates;
    if (settings.coupleAdiabaticStates) {
      // Normalize coefficients from configuration weights
      if (settings.configurationWeights.size() != 0) {
        for (unsigned iCol = 0; iCol < coefficients.cols(); ++iCol) {
          double norm =
              sqrt(coefficients.col(iCol).transpose() * _fdeetCalculator->getDeterminants() * coefficients.col(iCol));
          coefficients.col(iCol) = coefficients.col(iCol) / norm;
        }
      }
      // Couple adiabatic states
      if (iCoupling == settings.couple.size() - 1) {
        nAdiab = coefficients.cols();
        Eigen::MatrixXd HPrime = coefficients.transpose() *
                                 _fdeetCalculator->getHamiltonian().cwiseProduct(_fdeetCalculator->getDeterminants()) *
                                 coefficients;
        Eigen::MatrixXd S = coefficients.transpose() * _fdeetCalculator->getDeterminants() * coefficients;
        _fdeetCalculator->solveEigenValueProblem(HPrime, S, nAdiab, true, true);
        Eigen::MatrixXd newCoeff = coefficients * _fdeetCalculator->getLinCoeffs();
        _fdeetCalculator->printTransformedMatrices(HPrime, S, _fdeetCalculator->getLinCoeffs());
        _fdeetCalculator->setLinCoeffs(newCoeff);
        _fdeetCalculator->printResults(nAdiab);
      }
      else {
        // Couple diabatic states for subsequent coupling of adiabatic states
        _fdeetCalculator->solveEigenValueProblem(_fdeetCalculator->getHamiltonian(),
                                                 _fdeetCalculator->getDeterminants(), nStates, false);
        _fdeetCalculator->printResults(nStates, false);
      }
    }
    else {
      // Couple diabatic states
      _fdeetCalculator->solveEigenValueProblem(_fdeetCalculator->getHamiltonian(), _fdeetCalculator->getDeterminants(), nStates);
      _fdeetCalculator->printResults(nStates);
    }
    // Store linear coefficients for coupling of adiabatic states
    if (settings.coupleAdiabaticStates && iCoupling < settings.couple.size() - 1) {
      Eigen::VectorXd coeff = _fdeetCalculator->getLinCoeffs().col(0);
      for (unsigned iState = 0; iState < settings.couple.back().size(); ++iState) {
        for (unsigned jState = 0; jState < nStates; ++jState) {
          if (settings.couple.back()[iState] == settings.couple[iCoupling][jState])
            coefficients(iState, iCoupling) = coeff(jState);
        }
      }
    }
    // Calculate and print adiabatic spin densities to cube file
    _fdeDiabController = nullptr;
    if (settings.spindensity == true or settings.spinpopulation == true) {
      if (settings.diskMode) {
        _fdeDiabController =
            std::make_shared<FDEDiabController>(_fdeetCalculator->getDensityMatrixFiles(), _fdeetCalculator->getLinCoeffs(),
                                                _fdeetCalculator->getDeterminants(), superSystem, nStates, nAdiab);
      }
      else {
        _fdeDiabController =
            std::make_shared<FDEDiabController>(_fdeetCalculator->getTransDensMats(), _fdeetCalculator->getLinCoeffs(),
                                                _fdeetCalculator->getDeterminants(), superSystem, nStates, nAdiab);
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

std::shared_ptr<SystemController>
ElectronTransferTask::setUpSupersystem(std::vector<std::shared_ptr<SystemController>> coupleSystems, std::string name) {
  Settings superSysSettings = _act[0]->getSettings();
  superSysSettings.spin = 0;
  superSysSettings.charge = 0;
  superSysSettings.name = name.c_str();
  superSysSettings.path = "sys-";
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), superSysSettings);
  SystemAdditionTask<Options::SCF_MODES::UNRESTRICTED> additionTask(supersystem, coupleSystems);
  additionTask.settings.addOccupiedOrbitals = false;
  additionTask.run();
  return supersystem;
}

std::shared_ptr<FDEETController> ElectronTransferTask::getFDEETController() {
  if (!_fdeetController)
    this->run();
  return _fdeetController;
}

std::shared_ptr<FDEETCalculator> ElectronTransferTask::getFDEETCalculator() {
  if (!_fdeetCalculator)
    this->run();
  return _fdeetCalculator;
}

std::shared_ptr<FDEDiabController> ElectronTransferTask::getFDEDiabController() {
  return _fdeDiabController;
}

} /* namespace Serenity */
