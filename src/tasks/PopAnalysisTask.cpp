/**
 * @file   PopAnalysisTask.cpp
 *
 * @date   last rework Jun 30, 2017
 * @author Thomas Dresselhaus, last rework Jan Unsleber
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
#include "tasks/PopAnalysisTask.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/BeckePopulationCalculator.h"
#include "analysis/populationAnalysis/CHELPGPopulationCalculator.h"
#include "analysis/populationAnalysis/CM5PopulationCalculator.h"
#include "analysis/populationAnalysis/HirshfeldPopulationCalculator.h"
#include "analysis/populationAnalysis/IAOPopulationCalculator.h"
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/Basis.h" //Shell wise populations.
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "geometry/Atom.h"
#include "geometry/AtomType.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutput.h"
#include "system/SystemController.h"

namespace Serenity {
template<>
PopAnalysisTask<RESTRICTED>::PopAnalysisTask(std::shared_ptr<SystemController> systemController)
  : _systemController(systemController), _modestring(" ") {
  assert(systemController);
}

template<>
PopAnalysisTask<UNRESTRICTED>::PopAnalysisTask(std::shared_ptr<SystemController> systemController)
  : _systemController(systemController), _modestring() {
  assert(systemController);
  _modestring.alpha = "a";
  _modestring.beta = "b";
}

template<Options::SCF_MODES SCFMode>
void PopAnalysisTask<SCFMode>::run() {
  if (settings.algorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN) {
    auto es = _systemController->getElectronicStructure<SCFMode>();
    MullikenPopulationCalculator<SCFMode> calculator;
    auto pop = calculator.calculateMullikenPopulations(_systemController);
    this->print("Mulliken", pop);

    /*
     * Also print orbital occupations in this special case
     */
    std::shared_ptr<AtomCenteredBasisController> bfCont = nullptr;
    bfCont = std::dynamic_pointer_cast<AtomCenteredBasisController>(es->getDensityMatrix().getBasisController());
    // Only try this if the cast to an AtomCenteredBasisController worked
    if (bfCont) {
      auto indices = bfCont->getBasisIndices();
      auto dMat = es->getDensityMatrix();
      printSubSectionTitle("Mulliken Orbital Population Analysis");
      auto bfpop =
          calculator.calculateBasisFunctionPopulations(dMat, es->getOneElectronIntegralController()->getOverlapIntegrals());
      std::vector<char> symbols;
      const auto& basis = bfCont->getBasis();
      for (auto& s : basis) {
        for (unsigned int i = 0; i < s->getNContracted(); i++) {
          symbols.push_back(ANGMOM_TO_LABEL[s->getAngularMomentum()]);
        }
      }
      const auto& atoms = _systemController->getAtoms();
      printf("%4s %5s %2s %5s %3s %23s\n", "", "Atom#", "AT", " BF# ", "BFT", " Population (tot or a/b)");
      printf("%4s %5s %2s %5s %3s %23s\n", "", "-----", "--", "-----", "---", "-------------------------");
      for (unsigned int idx = 0; idx < indices.size(); idx++) {
        for (unsigned int i = indices[idx].first; i < indices[idx].second; i++) {
          if (i == indices[idx].first) {
            printf("%4s %5d %2s %5d  %c   ", "", (idx + 1), atoms[idx]->getAtomType()->getName().c_str(), i + 1, symbols[i]);
          }
          else {
            printf("%4s %5s %2s %5d  %c   ", "", "", "", i + 1, symbols[i]);
          }
          for_spin(bfpop) {
            printf("%11f ", bfpop_spin(i));
          };
          printf("\n");
        }
      }
      printf("\n");
    }
  }

  if (settings.algorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::HIRSHFELD) {
    auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
        128, 0.0, 2, _systemController->getBasisController(), _systemController->getGridController());
    auto densOnGridCalc = std::make_shared<DensityOnGridCalculator<SCFMode>>(basFuncOnGridController, 0.0);
    auto densMatController = _systemController->getElectronicStructure<SCFMode>()->getDensityMatrixController();
    auto densOnGridController =
        std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densOnGridCalc, densMatController);

    HirshfeldPopulationCalculator<SCFMode> hirshFeld(_systemController, densOnGridController);
    this->print("Hirshfeld", hirshFeld.getAtomPopulations());
  }
  if (settings.algorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::BECKE) {
    auto densPtr =
        std::make_shared<DensityMatrix<SCFMode>>(_systemController->getElectronicStructure<SCFMode>()->getDensityMatrix());
    auto beckeCalculator = BeckePopulationCalculator<SCFMode>(_systemController, densPtr);
    auto atomPops = beckeCalculator.getAtomPopulations();
    auto spinPops = beckeCalculator.getSpinPopulations();
    const auto& atoms = _systemController->getAtoms();
    unsigned int nAtoms = atoms.size();
    printSubSectionTitle("Becke Population Analysis");
    printf("%4s %5s %9s %11s %11s\n", "", "No.", "Atomtype", "Electrons", "Spin Pop.");
    printf("%4s %5s %9s %11s %11s\n", "", "---", "--------", "---------", "---------");
    for (unsigned int i = 0; i < nAtoms; ++i) {
      printf("%4s %5d %9s %11f %11f\n", "", (i + 1), atoms[i]->getAtomType()->getName().c_str(), atomPops[i], spinPops[i]);
    }
    printf("%4s %5s %9s %11s %11s\n", "", "---", "--------", "---------", "---------");
    printf("%4s %5s %9s %11f %11f\n\n\n", "", "sum", "", atomPops.sum(), spinPops.sum());
  }

  if (settings.algorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::IAO ||
      settings.algorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell) {
    // Ensure that an electronic structure is present.
    _systemController->getElectronicStructure<SCFMode>()->getEnergy();
    auto pop = IAOPopulationCalculator<SCFMode>::calculateIAOPopulations(_systemController);
    this->print("IAO", pop);

    if (settings.algorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell) {
      std::shared_ptr<AtomCenteredBasisController> bfCont =
          _systemController->getAtomCenteredBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION);
      auto indices = bfCont->getBasisIndicesRed();
      printSubSectionTitle("IAO-Shell Orbital Population Analysis");
      auto iaoOrbitalWiseShell = IAOPopulationCalculator<SCFMode>::calculateShellwiseOrbitalPopulations(_systemController);

      SPMatrix<SCFMode> bfpop;
      for_spin(bfpop, iaoOrbitalWiseShell) {
        bfpop_spin = (SCFMode == RESTRICTED ? 2.0 : 1.0) * iaoOrbitalWiseShell_spin.rowwise().sum();
      };
      std::vector<char> symbols;
      const auto& basis = bfCont->getBasis();
      for (auto& s : basis) {
        symbols.push_back(ANGMOM_TO_LABEL[s->getAngularMomentum()]);
      }
      const auto& atoms = _systemController->getAtoms();
      printf("%4s %5s %2s %5s %3s %23s\n", "", "Atom#", "AT", " BF# ", "BFT", " Population (tot or a/b)");
      printf("%4s %5s %2s %5s %3s %23s\n", "", "-----", "--", "-----", "---", "-------------------------");
      for (unsigned int idx = 0; idx < indices.size(); idx++) {
        for (unsigned int i = indices[idx].first; i < indices[idx].second; i++) {
          if (i == indices[idx].first) {
            printf("%4s %5d %2s %5d  %c   ", "", (idx + 1), atoms[idx]->getAtomType()->getName().c_str(), i + 1, symbols[i]);
          }
          else {
            printf("%4s %5s %2s %5d  %c   ", "", "", "", i + 1, symbols[i]);
          }
          for_spin(bfpop) {
            printf("%11f ", bfpop_spin(i));
          };
          printf("\n");
        }
      }
      printf("\n");
    }
  }

  if (settings.algorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::CM5) {
    auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
        128, 0.0, 2, _systemController->getBasisController(), _systemController->getGridController());
    auto densOnGridCalc = std::make_shared<DensityOnGridCalculator<SCFMode>>(basFuncOnGridController, 0.0);
    auto densMatController = _systemController->getElectronicStructure<SCFMode>()->getDensityMatrixController();
    auto densOnGridController =
        std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densOnGridCalc, densMatController);

    std::shared_ptr<HirshfeldPopulationCalculator<SCFMode>> hirshFeld =
        std::make_shared<HirshfeldPopulationCalculator<SCFMode>>(_systemController, densOnGridController);
    CM5PopulationCalculator<SCFMode> calculator(_systemController, hirshFeld);
    this->print("CM5", calculator.getAtomPopulations());
  }

  if (settings.algorithm == Options::POPULATION_ANALYSIS_ALGORITHMS::CHELPG) {
    CHELPGPopulationCalculator<SCFMode> calculator(_systemController);
    auto atoms = _systemController->getAtoms();
    auto atomPops = calculator.getAtomPopulations();
    Eigen::VectorXd electrons = Eigen::VectorXd(atoms.size()).setZero();
    for (unsigned int i = 0; i < atoms.size(); i++) {
      electrons[i] = atoms[i]->getEffectiveCharge() - atomPops[i];
    }
    printSubSectionTitle("CHELPG Population Analysis");
    printf("%4s %5s %9s %11s %11s\n", "", "No.", "Atomtype", "Electrons", "Charge");
    printf("%4s %5s %9s %11s %11s\n", "", "---", "--------", "---------", "------");
    for (unsigned int i = 0; i < atoms.size(); i++) {
      printf("%4s %5d %9s %11f %11f\n", "", (i + 1), atoms[i]->getAtomType()->getName().c_str(), electrons[i], atomPops[i]);
    }
    printf("\n");
  }
}

template<>
void PopAnalysisTask<RESTRICTED>::print(std::string type, const SpinPolarizedData<RESTRICTED, Eigen::VectorXd>& populations) {
  // getting the atoms of the system and the number of atoms.
  const auto& atoms = _systemController->getAtoms();
  unsigned int nAtoms = atoms.size();
  printSubSectionTitle(type + " Population Analysis");
  printf("%4s %5s %9s %11s %11s\n", "", "No.", "Atomtype", "Electrons", "Charge");
  printf("%4s %5s %9s %11s %11s\n", "", "---", "--------", "---------", "------");
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto charge = atoms[i]->getEffectiveCharge() - populations[i];
    if (atoms[i]->usesECP()) {
      printf("%4s %5d %9s %11f %11f %2s %3i %24s\n", "", (i + 1), atoms[i]->getAtomType()->getName().c_str(),
             populations[i], charge, " (", atoms[i]->getNECPElectrons(), " core electrons present)");
    }
    else {
      printf("%4s %5d %9s %11f %11f\n", "", (i + 1), atoms[i]->getAtomType()->getName().c_str(), populations[i], charge);
    }
  }
  printf("\n");
}

template<>
void PopAnalysisTask<UNRESTRICTED>::print(std::string type, const SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>& populations) {
  // getting the atoms of the system and the number of atoms.
  const auto& atoms = _systemController->getAtoms();
  unsigned int nAtoms = atoms.size();
  printSubSectionTitle(type + " Population Analysis");
  printf("%4s %5s %9s %11s %11s %11s\n", "", "No.", "Atomype", "Electrons", "Spin Pop.", "Charge");
  printf("%4s %5s %9s %11s %11s %11s\n", "", "---", "-------", "---------", "---------", "------");
  for (unsigned int i = 0; i < nAtoms; i++) {
    auto charge = atoms[i]->getEffectiveCharge() - (populations.alpha[i] + populations.beta[i]);
    if (atoms[i]->usesECP()) {
      printf("%4s %5d %9s %11f %11f %11f %2s %3i %24s\n", "", (i + 1), atoms[i]->getAtomType()->getName().c_str(),
             populations.alpha[i] + populations.beta[i], populations.alpha[i] - populations.beta[i], charge, " (",
             atoms[i]->getNECPElectrons(), " core electrons present)");
    }
    else {
      printf("%4s %5d %9s %11f %11f %11f\n", "", (i + 1), atoms[i]->getAtomType()->getName().c_str(),
             populations.alpha[i] + populations.beta[i], populations.alpha[i] - populations.beta[i], charge);
    }
  }
  printf("\n");
}

template class PopAnalysisTask<Options::SCF_MODES::RESTRICTED>;
template class PopAnalysisTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */