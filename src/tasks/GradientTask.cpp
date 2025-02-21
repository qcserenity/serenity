/**
 * @file   GradientTask.cpp
 *
 * @date   Mar 23, 2015
 * @author Kevin Klahr
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
#include "tasks/GradientTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "geometry/gradients/GeometryGradientCalculator.h"
#include "geometry/gradients/NumericalGeomGradCalc.h"
#include "io/FormattedOutput.h"
#include "io/IOOptions.h"
#include "misc/Timing.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/bundles/FDEPotentialBundleFactory.h"
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/ScfTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GradientTask<SCFMode>::GradientTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                                    const std::vector<std::shared_ptr<SystemController>>& passiveSystems)
  : _activeSystems(activeSystems), _passiveSystems(passiveSystems) {
}

template<Options::SCF_MODES SCFMode>
GradientTask<SCFMode>::~GradientTask() = default;

template<Options::SCF_MODES SCFMode>
void GradientTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _activeSystems);
  takeTime("Gradient Calculation");
  bool info(iOOptions.printSCFCycleInfo);
  bool results(iOOptions.printSCFResults);
  bool check(iOOptions.gridAccuracyCheck = false);
  int timings(iOOptions.timingsPrintLevel);
  iOOptions.gridAccuracyCheck = false;
  if (!(settings.print)) {
    iOOptions.printSCFCycleInfo = false;
    iOOptions.printSCFResults = false;
  }

  if (_passiveSystems.size() == 0 && _activeSystems.size() == 1) {
    if (settings.gradType == Options::GRADIENT_TYPES::NUMERICAL) {
      NumericalGeomGradCalc<SCFMode> numGradCalc(settings.numGradStepSize);
      if (settings.print)
        printSubSectionTitle("Numerical Gradient Calculation");

      numGradCalc.calcGradients(_activeSystems[0]);

      if (settings.transInvar) {
        _activeSystems[0]->getGeometry()->makeGradientsTranslationallyInvariant();
      }
    }
    else {
      auto es = _activeSystems[0]->template getElectronicStructure<SCFMode>();
      if (!es->potentialsAvailable()) {
        ScfTask<SCFMode> scf(_activeSystems[0]);
        scf.run();
      }
      Timings::takeTime("Tech. -          Gradient Step");
      auto potBundle = es->getPotentialBundle();
      auto potentialGradients = potBundle->getGradients();
      Matrix<double> dispCorr(_activeSystems[0]->getGeometry()->getAtoms().size(), 3);
      dispCorr.setZero();
      if (_activeSystems[0]->getSettings().dft.dispersion != Options::DFT_DISPERSION_CORRECTIONS::NONE) {
        // Dispersion Correction components
        dispCorr += DispersionCorrectionCalculator::calcDispersionGradientCorrection(
            _activeSystems[0]->getSettings().dft.dispersion, _activeSystems[0]->getGeometry(),
            _activeSystems[0]->getSettings().dft.functional);
      }
      Matrix<double> gradient = potentialGradients + dispCorr;
      _activeSystems[0]->getGeometry()->setGradients(gradient);
      if (_activeSystems[0]->hasExternalCharges()) {
        _activeSystems[0]->setPointChargeGradients(potBundle->getPointChargeGradients());
      }
      if (settings.transInvar) {
        _activeSystems[0]->getGeometry()->makeGradientsTranslationallyInvariant();
      }
      Timings::timeTaken("Tech. -          Gradient Step");
    }
  }
  else {
    // Change the SCF mode of the system. The FAT will check for
    // the system SCF mode and adjust the FDE-SCF mode accordingly.
    // Hence, we need to ensure that we actually calculate the system
    // electronic structure using the correct SCF mode. Otherwise the
    // FAT does not make any sense.
    std::vector<Options::SCF_MODES> scfModes;
    for (auto& sys : _activeSystems) {
      scfModes.push_back(sys->getSCFMode());
      sys->setSCFMode(SCFMode);
    }
    /*
     * Initial FaT
     */
    FreezeAndThawTask<SCFMode> task(_activeSystems, _passiveSystems);
    task.settings.embedding = settings.embedding;
    if (settings.embedding.embeddingMode != Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
      throw SerenityError("Only non-additive kinetic embedding is supported for embedded Gradients!");
    }
    task.settings.gridCutOff = settings.FDEgridCutOff;
    task.settings.maxCycles = settings.FaTmaxCycles;
    task.settings.convThresh = settings.FaTenergyConvThresh;
    task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
    task.settings.printResults = false;
    task.run();
    if (settings.gradType == Options::GRADIENT_TYPES::NUMERICAL) {
      NumericalGeomGradCalc<SCFMode> numGradCalc(settings.numGradStepSize);
      if (settings.print)
        printSubSectionTitle("Numerical FDE Gradient Calculation");

      numGradCalc.calcFDEGradients(_activeSystems, _passiveSystems, settings.embedding.naddKinFunc,
                                   settings.embedding.naddXCFunc, settings.FDEgridCutOff, settings.FaTmaxCycles,
                                   settings.FaTenergyConvThresh, settings.embedding.dispersion);

      for (unsigned int i = 0; i < _activeSystems.size(); i++) {
        printBigCaption((std::string) "Active System: " + (i + 1));
        if (settings.transInvar) {
          _activeSystems[i]->getGeometry()->makeGradientsTranslationallyInvariant();
        }
      }
    }
    else {
      // Active system cycles
      for (unsigned int nSystem = 0; nSystem < _activeSystems.size(); nSystem++) {
        // Output and init
        if (settings.print)
          printBigCaption((std::string) "Active System: " + (nSystem + 1));
        auto activeSystem = _activeSystems[nSystem];

        // Set up systems
        auto passiveSystems = _activeSystems;
        passiveSystems.erase(std::remove(passiveSystems.begin(), passiveSystems.end(), activeSystem), passiveSystems.end());

        passiveSystems.insert(passiveSystems.end(), _passiveSystems.begin(), _passiveSystems.end());

        auto supersystemGeometry = std::make_shared<Geometry>();
        *supersystemGeometry += *activeSystem->getGeometry();
        for (auto sys : passiveSystems) {
          *supersystemGeometry += *sys->getGeometry();
        }

        FDETask<SCFMode> task(activeSystem, passiveSystems);
        task.settings.embedding = settings.embedding;
        if (settings.embedding.embeddingMode != Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
          throw SerenityError("Only non-additive kinetic embedding is supported for embedded Gradients!");
        }
        task.settings.gridCutOff = settings.FDEgridCutOff;
        // There is no point in converging the FDE-SCF again. We only need its potential bundle.
        task.settings.skipSCF = true;
        task.run();

        Timings::takeTime("Tech. -          Gradient Step");
        auto es = activeSystem->template getElectronicStructure<SCFMode>();
        auto potBundle = es->getPotentialBundle();
        auto potentialGradients = potBundle->getGradients();
        Matrix<double> ccRepDerivative(activeSystem->getGeometry()->getAtoms().size(), 3);
        ccRepDerivative.setZero();
        for (auto sys : passiveSystems) {
          ccRepDerivative += CoreCoreRepulsionDerivative::calculateDerivative(activeSystem->getAtoms(), sys->getAtoms());
        }

        Matrix<double> dispCorr(activeSystem->getGeometry()->getAtoms().size(), 3);
        dispCorr.setZero();
        if (!(settings.embedding.dispersion == Options::DFT_DISPERSION_CORRECTIONS::NONE)) {
          dispCorr += DispersionCorrectionCalculator::calcDispersionGradientCorrection(
                          settings.embedding.dispersion, supersystemGeometry, settings.embedding.naddXCFunc)
                          .topRows(activeSystem->getGeometry()->getNAtoms());
        }

        // Add interaction parts into atoms
        Matrix<double> gradientSum = potentialGradients + ccRepDerivative + dispCorr;

        activeSystem->getGeometry()->setGradients(gradientSum);
        if (activeSystem->hasExternalCharges()) {
          activeSystem->setPointChargeGradients(potBundle->getPointChargeGradients());
        }

        if (settings.transInvar) {
          activeSystem->getGeometry()->makeGradientsTranslationallyInvariant();
        }
        Timings::timeTaken("Tech. -          Gradient Step");
      }
    }
    // Reset the SCFMode of the system.
    for (unsigned int iSys = 0; iSys < _activeSystems.size(); ++iSys) {
      _activeSystems[iSys]->setSCFMode(scfModes[iSys]);
    }
  }
  /*
   * Output
   */
  if (settings.print) {
    for (auto sys : _activeSystems) {
      printBigCaption((std::string) "Active System: " + sys->getSystemName());
      sys->getGeometry()->printGradients();
      std::cout << "\n" << std::endl;
    }
  }
  if (settings.printTotal)
    this->printTotalGradient();

  printPointChargeGradients();

  iOOptions.printSCFCycleInfo = info;
  iOOptions.printSCFResults = results;
  iOOptions.timingsPrintLevel = timings;
  iOOptions.gridAccuracyCheck = check;
  timeTaken(3, "Gradient Calculation");
}

template<Options::SCF_MODES SCFMode>
void GradientTask<SCFMode>::printPointChargeGradients() {
  std::vector<std::shared_ptr<SystemController>> allSystems = _activeSystems;
  allSystems.insert(allSystems.end(), _passiveSystems.begin(), _passiveSystems.end());
  unsigned int nSystemsWithGradients = 0;
  Eigen::MatrixXd totalPointChargeGradients;
  for (const auto& sys : allSystems) {
    if (sys->hasExternalCharges()) {
      const auto& pointChargeGradients = sys->getPointChargeGradients();
      if (nSystemsWithGradients == 0) {
        totalPointChargeGradients = pointChargeGradients;
      }
      else {
        if (pointChargeGradients.rows() != totalPointChargeGradients.rows()) {
          throw SerenityError(
              "The number of point charges must be identical for every system! Otherwise, forces and energies"
              " are calculated inconstent.");
        }
        totalPointChargeGradients += pointChargeGradients;
      }
      nSystemsWithGradients++;
    }
  }
  if (nSystemsWithGradients > 0 && nSystemsWithGradients != allSystems.size()) {
    throw SerenityError("Only some systems are calculated with external charges! This is inconsistent.");
  }
  if (nSystemsWithGradients > 0) {
    printSmallCaption("Point Charge Gradients (a.u.)");
    for (unsigned int iRow = 0; iRow < totalPointChargeGradients.rows(); ++iRow) {
      printf("%4s %4d %2s %+15.10f %+15.10f %+15.10f\n", "", iRow, "Q", totalPointChargeGradients(iRow, 0),
             totalPointChargeGradients(iRow, 1), totalPointChargeGradients(iRow, 2));
    }
  }
}

template<Options::SCF_MODES SCFMode>
void GradientTask<SCFMode>::printTotalGradient() {
  /*
   * Note: Serenity is compatible with the program SNF, which
   * searches for the string 'Total Geometry Gradients' in
   * Serenity's output. Thus, please make sure that at least this
   * substring remains intact!
   */
  printSmallCaption("Total Geometry Gradients (a.u.)");
  int counter = 0;
  for (auto sys : _activeSystems) {
    const auto atoms = sys->getGeometry()->getAtoms();
    for (const auto& atom : atoms) {
      ++counter;
      printf("%4s %4d %2s %+15.10f %+15.10f %+15.10f\n", "", counter, atom->getAtomType()->getElementSymbol().c_str(),
             atom->getGradient()[0], atom->getGradient()[1], atom->getGradient()[2]);
    }
  }
}

template class GradientTask<Options::SCF_MODES::RESTRICTED>;
template class GradientTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
