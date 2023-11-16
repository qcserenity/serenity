/**
 * @file   FDEETController.cpp
 *
 * @date   Apr. 20, 2020
 * @author Patrick Eschenbach
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

/* Include Class Header*/
#include "postHF/ET/FDEETController.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/OrbitalController.h"
#include "io/FormattedOutputStream.h" // formatted output
#include "system/SystemController.h"

namespace Serenity {

FDEETController::FDEETController(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                                 unsigned nStatesInput, unsigned nSysPerState)
  : _activeSystems(activeSystems), _nStatesInput(nStatesInput), _nSysPerState(nSysPerState) {
  this->setupStateVector();
}

unsigned FDEETController::getNStatesInput() {
  return _nStatesInput;
}

unsigned FDEETController::getNSysPerState() {
  return _nSysPerState;
}

const std::vector<std::shared_ptr<SystemController>>& FDEETController::getActiveSystems() {
  return _activeSystems;
}

void FDEETController::setStateVector(std::vector<std::vector<std::shared_ptr<SystemController>>> newV) {
  _stateVector = newV;
}

std::vector<std::vector<std::shared_ptr<SystemController>>> FDEETController::getStateVector() {
  return _stateVector;
}

unsigned FDEETController::getNStatesCouple() {
  return _nStatesCouple;
}

void FDEETController::setNStatesCouple(unsigned nStates) {
  _nStatesCouple = nStates;
}

std::shared_ptr<std::vector<unsigned>> FDEETController::getIndexVector() {
  return _indexVector;
}

void FDEETController::setIndexVector(std::shared_ptr<std::vector<unsigned>> states) {
  _indexVector = states;
}

void FDEETController::setupIndexVector(std::vector<unsigned> couples) {
  auto stateIndex = std::make_shared<std::vector<unsigned>>(_nStatesCouple, 0);
  for (unsigned item = 0; item < couples.size(); ++item) {
    (*stateIndex)[item] = couples[item] - 1;
  }
  std::sort((*stateIndex).begin(), (*stateIndex).end());
  this->setIndexVector(stateIndex);
}

Eigen::MatrixXi FDEETController::getNOccState() {
  return _nOccState;
}

void FDEETController::setNOccState(Eigen::MatrixXi newOcc) {
  _nOccState = newOcc;
}

Eigen::VectorXi FDEETController::getNBfState() {
  return _nBfState;
}

void FDEETController::setNBfState(Eigen::VectorXi newBf) {
  _nBfState = newBf;
}

void FDEETController::setupOrbitalInfo() {
  // vectors containing numbers of occupied orbs and basis functions for each state
  Eigen::MatrixXi nOccState = Eigen::MatrixXi::Zero(_nStatesCouple, 2);
  Eigen::VectorXi nBFState = Eigen::VectorXi::Zero(_nStatesCouple);
  // Get all necessary dimensions right
  for (unsigned iState = 0; iState < _nStatesCouple; ++iState) {
    for (unsigned iSys = 0; iSys < this->getNSysPerState(); ++iSys) {
      nBFState(iState) += _stateVector[(*_indexVector)[iState]][iSys]->getBasisController()->getNBasisFunctions();
      auto nOcc = _stateVector[(*_indexVector)[iState]][iSys]->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
      unsigned iSpin = 0;
      for_spin(nOcc) {
        nOccState(iState, iSpin) += nOcc_spin;
        ++iSpin;
      };
    }
  }
  this->setNOccState(nOccState);
  this->setNBfState(nBFState);
  printTableHead(" Orbital Info ");
  OutputControl::nOut << std::string(1, 'o') << std::string(46, '-') << std::string(1, 'o') << std::endl;
  printf("| State |   nBF   |  nOcc Alpha  |  nOcc Beta  |\n");
  OutputControl::nOut << std::string(1, 'o') << std::string(7, '-') << std::string(1, '+') << std::string(9, '-')
                      << std::string(1, '+') << std::string(14, '-') << std::string(1, '+') << std::string(13, '-')
                      << std::string(1, 'o') << std::endl;
  for (unsigned iState = 0; iState < nBFState.size(); ++iState) {
    printf("|  %3i  |  %5i  |  %6i      |  %6i     |\n", iState, nBFState[iState], nOccState(iState, 0), nOccState(iState, 1));
  }
  OutputControl::nOut << std::string(1, 'o') << std::string(46, '-') << std::string(1, 'o') << std::endl;
  OutputControl::nOut << std::endl;
  OutputControl::nOut << std::endl;
}

std::vector<std::shared_ptr<SystemController>> FDEETController::getSubsystems() {
  this->collectSubsystems();
  return _allSystems;
}

void FDEETController::setSubsystems(std::vector<std::shared_ptr<SystemController>> systems) {
  _allSystems = systems;
}

std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>> FDEETController::getStateCoeffs() {
  this->setupStateCoeffs();
  return _stateCoeffs;
}

void FDEETController::setupStateCoeffs() {
  // Build state specific coefficient matrices (alpha and beta)
  // this object contains a coefficient matrix for each state, where each coefficient matrix
  // contains the coefficient blocks of the subsystems within this state on its diagonal
  unsigned iSpin = 0;
  std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>> tmp(_nStatesCouple,
                                                                 std::vector<std::shared_ptr<Eigen::MatrixXd>>(2));
  for (unsigned iState = 0; iState < _nStatesCouple; ++iState) {
    std::vector<unsigned> iCount = {0, 0};
    std::vector<unsigned> jCount = {0, 0};
    for (unsigned iSpin = 0; iSpin < 2; ++iSpin) {
      // initialize block of coeff matrix
      tmp[iState][iSpin] =
          std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(_nBfState(iState), _nOccState(iState, iSpin)));
    }
    for (unsigned iSys = 0; iSys < _nSysPerState; ++iSys) {
      // Get coefficients and nOcc of this subsystem
      auto coeff =
          _stateVector[(*_indexVector)[iState]][iSys]->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
      auto nOcc = _stateVector[(*_indexVector)[iState]][iSys]->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
      unsigned iSize = _stateVector[(*_indexVector)[iState]][iSys]->getBasisController()->getNBasisFunctions();
      iSpin = 0;
      for_spin(nOcc, coeff) {
        unsigned jSize = nOcc_spin;
        // get coefficient blocks
        tmp[iState][iSpin]->block(iCount[0], jCount[iSpin], iSize, jSize) = coeff_spin.leftCols(nOcc_spin);
        jCount[iSpin] += jSize;
        ++iSpin;
      };
      // Refers to basis functions, so it's the same for alpha and beta, hence only [0] considered.
      iCount[0] += iSize;
    } /* iSys */
  }   /* iState */
  _stateCoeffs = tmp;
}

void FDEETController::collectSubsystems() {
  // Construct vector containing all subsystems for this coupling problem
  std::vector<std::shared_ptr<SystemController>> systems;
  for (unsigned iState = 0; iState < _nStatesCouple; ++iState) {
    for (const auto& sys : _stateVector[(*_indexVector)[iState]]) {
      systems.emplace_back(sys);
    }
  }
  this->setSubsystems(systems);
}

void FDEETController::setupStateVector() {
  _stateVector = std::vector<std::vector<std::shared_ptr<SystemController>>>(this->getNStatesInput());
  auto act = this->getActiveSystems();
  unsigned iSystem = 0;
  for (unsigned iState = 0; iState < this->getNStatesInput(); ++iState) {
    for (unsigned iSys = 0; iSys < this->getNSysPerState(); ++iSys) {
      _stateVector[iState].emplace_back(act[iSystem + iSys]);
    }
    iSystem += this->getNSysPerState();
  }
}

} /* namespace Serenity */