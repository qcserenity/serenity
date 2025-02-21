/**
 * @file   FDEDiabController.cpp
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
#include "postHF/ET/FDEDiabController.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/BeckePopulationCalculator.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/matrices/DensityMatrixController.h"
#include "geometry/Geometry.h"
#include "io/CubeFileWriter.h"
#include "io/Filesystem.h"
#include "io/FormattedOutputStream.h"
#include "io/HDF5.h"
#include "parameters/Constants.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/PlotTask.h"

namespace Serenity {

FDEDiabController::FDEDiabController(const std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>>& transDensMats,
                                     const Eigen::MatrixXd& linCoeffs, const Eigen::MatrixXd& determinants,
                                     std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple, unsigned nStatesAdiab)
  : _transDensMats(transDensMats),
    _linCoeffs(linCoeffs),
    _determinants(determinants),
    _superSystem(superSystem),
    _nStates(nStatesCouple),
    _nAdiab(nStatesAdiab) {
  if (_superSystem->hasElectronicStructure<UNRESTRICTED>())
    throw SerenityError("The supersystem should not posses an electronic structure!");
}

FDEDiabController::FDEDiabController(std::shared_ptr<std::vector<std::string>> densityMatrixFiles,
                                     const Eigen::MatrixXd& linCoeffs, const Eigen::MatrixXd& determinants,
                                     std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple, unsigned nStatesAdiab)
  : _transDensMats({}),
    _linCoeffs(linCoeffs),
    _determinants(determinants),
    _superSystem(superSystem),
    _nStates(nStatesCouple),
    _nAdiab(nStatesAdiab),
    _densityMatrixFiles(densityMatrixFiles) {
  if (_superSystem->hasElectronicStructure<UNRESTRICTED>())
    throw SerenityError("The supersystem should not posses an electronic structure!");
}

std::vector<std::shared_ptr<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> FDEDiabController::getAdiabDensities() {
  if (_adiabDensities.size() == 0) {
    this->calcAdiabDensities();
  }
  return _adiabDensities;
}

void FDEDiabController::setAdiabDensities(std::vector<std::shared_ptr<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> newD) {
  _adiabDensities = newD;
}

void FDEDiabController::calcAdiabDensities() {
  auto adiabDensities = std::vector<std::shared_ptr<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>>(
      _nAdiab, std::make_shared<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>(_superSystem->getBasisController()));
  for (unsigned adiaState = 0; adiaState < adiabDensities.size(); ++adiaState) {
    auto adiaDens = std::make_shared<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>(_superSystem->getBasisController());
    unsigned ij = 0;
    for (unsigned iState = 0; iState < _nStates; ++iState) {
      for (unsigned jState = 0; jState < _nStates; ++jState) {
        auto diabDens = DensityMatrix<Options::SCF_MODES::UNRESTRICTED>(_superSystem->getBasisController());
        if (_transDensMats.size() > 0) {
          diabDens = DensityMatrix<Options::SCF_MODES::UNRESTRICTED>(_transDensMats[iState][jState]);
        }
        else {
          std::string pathUnres = (*_densityMatrixFiles)[ij];
          ++ij;
          std::vector<Eigen::MatrixXd> temp(
              2, Eigen::MatrixXd::Zero(_superSystem->getBasisController()->getNBasisFunctions(),
                                       _superSystem->getBasisController()->getNBasisFunctions()));
          HDF5::Filepath name(pathUnres);
          HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
          HDF5::dataset_exists(file, "densityMatrix_alpha");
          HDF5::load(file, "densityMatrix_alpha", temp[0]);
          HDF5::dataset_exists(file, "densityMatrix_beta");
          HDF5::load(file, "densityMatrix_beta", temp[1]);
          file.close();
          unsigned iSpin = 0;
          for_spin(diabDens) {
            diabDens_spin = temp[iSpin];
            ++iSpin;
          };
        }
        (*adiaDens).alpha +=
            _linCoeffs(iState, adiaState) * _linCoeffs(jState, adiaState) * _determinants(iState, jState) * diabDens.alpha;
        (*adiaDens).beta +=
            _linCoeffs(iState, adiaState) * _linCoeffs(jState, adiaState) * _determinants(iState, jState) * diabDens.beta;
      }
    }
    adiabDensities[adiaState] = adiaDens;
  }
  this->setAdiabDensities(adiabDensities);
}

void FDEDiabController::printAdiabDensities() {
  printSubSectionTitle("Writing Adiabatic Spin Densities to Cube");
  struct PlotTaskSettings plotSettings;
  plotSettings.gridSpacing = {0.12, 0.12, 0.12};
  auto cubeWriter = CubeFileWriter(_superSystem->getSettings(), plotSettings);
  for (unsigned adiaState = 0; adiaState < _adiabDensities.size(); ++adiaState) {
    print((std::string) "Printing density matrix to file: " + _superSystem->getSystemPath() + "adiabState" + adiaState +
          "_SpinDensity");
    cubeWriter.writeMatrixToGrid(_superSystem->getSystemPath() + "adiabState" + adiaState + "_SpinDensity",
                                 _superSystem->getGeometry(), _superSystem->getBasisController(),
                                 (*_adiabDensities[adiaState]).difference());
  }
  OutputControl::nOut << std::endl;
}

std::vector<Eigen::VectorXd> FDEDiabController::getSortedAtomPopulations() {
  return _sortedAtomPopulations;
}

void FDEDiabController::setSortedAtomPopulations(std::vector<Eigen::VectorXd> newPop) {
  _sortedAtomPopulations = newPop;
}

void FDEDiabController::calcPopulations(std::vector<unsigned int> list) {
  printSubSectionTitle("Performing Becke Spin Population Analysis");
  auto densPtrVec = this->getAdiabDensities();
  std::vector<Eigen::VectorXd> sortedAtomPopulations;
  for (auto index : list) {
    auto beckeCalculator = BeckePopulationCalculator<Options::SCF_MODES::UNRESTRICTED>(_superSystem, densPtrVec[index]);
    sortedAtomPopulations.push_back(beckeCalculator.getSpinPopulations());
  }
  this->setSortedAtomPopulations(sortedAtomPopulations);
}

void FDEDiabController::printPopulations(std::vector<unsigned int> list) {
  printSubSectionTitle("Spin Populations of Adiabatic Spin Densities");
  auto allPops = this->getSortedAtomPopulations();
  auto atoms = _superSystem->getGeometry()->getAtoms();
  auto atomSymbols = _superSystem->getGeometry()->getAtomSymbols();
  for (auto index : list) {
    printf(" Adiabatic State %4i\n", index);
    const auto& pops = allPops[index];
    OutputControl::nOut << std::string(1, 'o') << std::string(65, '-') << std::string(1, 'o') << std::endl;
    printf("|  Atoms  |      X      |      Y      |      Z      |      P      |\n");
    OutputControl::nOut << std::string(1, 'o') << std::string(9, '-') << std::string(1, '+') << std::string(13, '-')
                        << std::string(1, '+') << std::string(13, '-') << std::string(1, '+') << std::string(13, '-')
                        << std::string(1, '+') << std::string(13, '-') << std::string(1, 'o') << std::endl;
    for (unsigned i = 0; i < atoms.size(); ++i) {
      printf("|%5s    |%+12.5f |%+12.5f |%+12.5f |%+12.5f |\n", atomSymbols[i].c_str(), atoms[i]->getX() * BOHR_TO_ANGSTROM,
             atoms[i]->getY() * BOHR_TO_ANGSTROM, atoms[i]->getZ() * BOHR_TO_ANGSTROM, pops[i]);
    }
    OutputControl::nOut << std::string(1, 'o') << std::string(65, '-') << std::string(1, 'o') << std::endl;
    printf("|%-50s  %12.5f |\n", " Total Integral: ", pops.sum());
    OutputControl::nOut << std::string(1, 'o') << std::string(65, '-') << std::string(1, 'o') << std::endl;
    OutputControl::nOut << std::endl;
    OutputControl::nOut << std::endl;
  }
}

double FDEDiabController::getS2(unsigned iState) {
  auto densPtrVec = this->getAdiabDensities();
  double S = 0.5 * _superSystem->getSpin();
  // Calculate <S*S> as functional of density. Calculate spin density
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _superSystem->getSettings(), _superSystem->getBasisController(), _superSystem->getGridController());
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>>(
      basisFunctionOnGridController, _superSystem->getSettings().grid.blockAveThreshold);
  auto densityMatrixController =
      std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(*densPtrVec[iState]);
  auto densityOnGridController = std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>>(
      densityOnGridCalculator, densityMatrixController);
  auto& densityOnGrid = densityOnGridController->getDensityOnGrid();
  Eigen::VectorXd spinDensity = densityOnGrid.alpha - densityOnGrid.beta;
  // Integrage over all regions where spin density is negative
  for (unsigned int i = 0; i < spinDensity.rows(); ++i) {
    if (spinDensity(i) > 0)
      spinDensity(i) = 0.0;
  }
  double integral = spinDensity.dot(_superSystem->getGridController()->getWeights());
  return S * (S + 1) - integral;
}

void FDEDiabController::adiabaticWavefunctionAnalysis(std::vector<unsigned int> list) {
  printSubSectionTitle("Analysis of Adiabatic Spin Densities");
  printf(" State   tot. Int.    <S*S>\n");
  printf(" -----   ---------   -------\n");
  auto pops = this->getSortedAtomPopulations();
  for (auto index : list) {
    double sSquared = this->getS2(index);
    auto totInt = pops[index].sum();
    printf("%6i %10.5f %10.5f\n", index, totInt, sSquared);
  }
}

} /* namespace Serenity */
