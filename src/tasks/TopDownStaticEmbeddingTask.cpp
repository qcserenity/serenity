/**
 * @file TopDownStaticEmbeddingTask.cpp
 *
 * @date Mar 21, 2024
 * @author Moritz Bensberg
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
#include "tasks/TopDownStaticEmbeddingTask.h"
/* Include Serenity Internal Headers */
#include "basis/BasisFunctionMapper.h"    // Cache one electron integrals for QM/MM.
#include "data/ElectronicStructure.h"     // Get energy component controller.
#include "energies/EnergyContributions.h" // Get energy of the uncorrelated system.
#include "geometry/Geometry.h"            // Supersystem construction
#include "grid/AtomCenteredGridControllerFactory.h"
#include "integrals/OneElectronIntegralController.h" // Cache one electron integrals for QM/MM.
#include "io/FormattedOutputStream.h"                //Filtered output streams.
#include "math/Derivatives.h"
#include "settings/Settings.h" // Supersystem construction
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/LocalCorrelationTask.h"
#include "tasks/QuasiRestrictedOrbitalsTask.h"
#include "tasks/ScfTask.h" // Supersystem construction
/* Include Std and External Headers */
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
TopDownStaticEmbeddingTask<SCFMode>::TopDownStaticEmbeddingTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                                                                const std::vector<std::shared_ptr<SystemController>>& passiveSystems,
                                                                std::shared_ptr<SystemController> supersystem)
  : _activeSystems(activeSystems), _passiveSystems(passiveSystems), _supersystem(supersystem) {
}

template<Options::SCF_MODES SCFMode>
void TopDownStaticEmbeddingTask<SCFMode>::run() {
  printSectionTitle("Static DFT-Embedding Calculation");
  std::string embeddingMode, naddXCFunc, naddKinFunc;
  Options::resolve<Options::KIN_EMBEDDING_MODES>(embeddingMode, this->settings.lcSettings.embeddingSettings.embeddingMode);
  Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(naddXCFunc, this->settings.lcSettings.embeddingSettings.naddXCFunc);
  Options::resolve<CompositeFunctionals::KINFUNCTIONALS>(naddKinFunc, this->settings.lcSettings.embeddingSettings.naddKinFunc);
  OutputControl::nOut << "  embedding mode:  " << embeddingMode << std::endl;
  OutputControl::nOut << "  naddXCFunc:      " << naddXCFunc << std::endl;
  if (this->settings.lcSettings.embeddingSettings.embeddingMode == Options::KIN_EMBEDDING_MODES::NADD_FUNC)
    OutputControl::nOut << "  naddKinFunc:     " << naddKinFunc << std::endl;

  // Basis set truncation is impossible if the orbitals are not relaxed!
  BasisSetTruncationTaskSettings trunc;
  trunc.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  setUpSubsystems(_supersystem, _activeSystems, _passiveSystems, this->settings.loc, this->settings.split,
                  this->settings.add, trunc, settings.useQuasiRestrictedOrbitals);

  std::vector<std::shared_ptr<SystemController>> allSystemsButTheFirst = _passiveSystems;
  for (unsigned int iAct = 0; iAct < _activeSystems.size(); ++iAct) {
    if (iAct == 0) {
      continue;
    }
    allSystemsButTheFirst.insert(allSystemsButTheFirst.begin(), _activeSystems[iAct]);
  }
  printSubSectionTitle("Reference Energy Calculation");
  this->calculateReferenceEnergy();

  Eigen::VectorXd correlationEnergyContributions = Eigen::VectorXd::Zero(_activeSystems.size());
  if (this->settings.lcSettings.method != Options::PNO_METHOD::NONE) {
    correlationEnergyContributions = runCorrelationCalculations();
  }
  this->printResults(correlationEnergyContributions);
}

template<Options::SCF_MODES SCFMode>
void TopDownStaticEmbeddingTask<SCFMode>::updateSubsystemEnergy(unsigned int subsystemIndex,
                                                                std::vector<std::shared_ptr<SystemController>> allSubsystems) {
  std::vector<std::shared_ptr<SystemController>> allSubsystemsButI;
  for (unsigned int jSub = 0; jSub < allSubsystems.size(); ++jSub) {
    if (subsystemIndex == jSub) {
      continue;
    }
    allSubsystemsButI.push_back(allSubsystems[jSub]);
  }
  FDETask<SCFMode> fdeTask(allSubsystems[subsystemIndex], allSubsystemsButI);
  fdeTask.settings.embedding = this->settings.lcSettings.embeddingSettings;
  fdeTask.settings.lcSettings = this->settings.lcSettings;
  fdeTask.settings.loc = this->settings.loc;
  fdeTask.settings.skipSCF = true;
  fdeTask.run();
  allSubsystems[subsystemIndex]->setDiskMode(true);
}

template<Options::SCF_MODES SCFMode>
void TopDownStaticEmbeddingTask<SCFMode>::calculateReferenceEnergy() {
  _allSubsystems = _passiveSystems;
  _allSubsystems.insert(_allSubsystems.begin(), _activeSystems.begin(), _activeSystems.end());
  for (unsigned int iSub = 0; iSub < _allSubsystems.size(); ++iSub) {
    this->updateSubsystemEnergy(iSub, _allSubsystems);
  }
  this->updateFrozenEnvironmentEnergies();
}

template<Options::SCF_MODES SCFMode>
void TopDownStaticEmbeddingTask<SCFMode>::updateFrozenEnvironmentEnergies() {
  for (auto referenceSystem : _allSubsystems) {
    double totEnvEnergy = 0.0;
    for (auto sys : _allSubsystems) {
      if (sys == referenceSystem) {
        continue;
      }
      auto energyContribution = (sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT)
                                    ? ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY
                                    : ENERGY_CONTRIBUTIONS::HF_ENERGY;
      totEnvEnergy += sys->template getElectronicStructure<SCFMode>()->getEnergy(energyContribution);
      totEnvEnergy += 0.5 * sys->template getElectronicStructure<SCFMode>()->getEnergy(
                                ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
    }
    totEnvEnergy -= 0.5 * referenceSystem->template getElectronicStructure<SCFMode>()->getEnergy(
                              ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
    auto eCont = referenceSystem->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES, totEnvEnergy);
  }
}

template<Options::SCF_MODES SCFMode>
void TopDownStaticEmbeddingTask<SCFMode>::printResults(const Eigen::VectorXd& correlationEnergyContributions) {
  auto systemController = _allSubsystems[_referenceIndex];
  auto energyComponentController =
      systemController->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
  double referenceEnergy = NAN;
  if (energyComponentController->checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT)) {
    referenceEnergy = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT);
  }
  else if (energyComponentController->checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT)) {
    referenceEnergy = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT);
  }
  else {
    throw SerenityError("No embedded energy is available after FDE energy evaluation. Something went wrong!");
  }
  unsigned int lineLength = 80;
  unsigned int indent = 2;
  unsigned int numberLength = lineLength - 8 - 34 - 2 * indent;
  OutputControl::mOut << std::fixed << std::right;
  printSectionTitle("Embedding Results");
  if (abs(correlationEnergyContributions.sum()) > 0.0) {
    OutputControl::mOut << std::string(indent, ' ') << "Reference Energy                  " << std::setw(numberLength)
                        << referenceEnergy << " Hartree" << std::endl;
    OutputControl::mOut << std::string(lineLength, '-') << std::endl;
    for (unsigned int iAct = 0; iAct < _activeSystems.size(); ++iAct) {
      OutputControl::mOut << std::string(indent, ' ') << "Subsystem " << iAct << std::setw(4) << " total correlation energy"
                          << std::setw(numberLength - 2) << correlationEnergyContributions(iAct) << " Hartree" << std::endl;
    }
    OutputControl::mOut << std::string(lineLength, '-') << std::endl;
  }
  _finalTotalEnergy = std::make_unique<double>(correlationEnergyContributions.sum() + referenceEnergy);
  OutputControl::mOut << std::string(indent, ' ') << "Total Supersystem Energy          " << std::setw(numberLength)
                      << *_finalTotalEnergy << " Hartree\n"
                      << std::endl;
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd TopDownStaticEmbeddingTask<SCFMode>::runCorrelationCalculations() {
  printSubSectionTitle("Correlation Energy Calculation");
  Eigen::VectorXd correlationEnergyContributions = Eigen::VectorXd::Zero(_activeSystems.size());
  for (unsigned int iAct = 0; iAct < _activeSystems.size(); ++iAct) {
    auto activeSystemController = _activeSystems[iAct];
    if (activeSystemController->getLastSCFMode() != RESTRICTED) {
      throw SerenityError("Local correlation calculations are only available for spin RESTRICTED subsystems.");
    }

    if (activeSystemController->getSettings().method != Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
      continue;
    }
    std::vector<std::shared_ptr<SystemController>> allSystems = _passiveSystems;
    for (unsigned int jAct = 0; jAct < _activeSystems.size(); ++jAct) {
      if (iAct == jAct) {
        continue;
      }
      allSystems.insert(allSystems.begin(), _activeSystems[jAct]);
    }
    // Run localCorrelation task
    LocalCorrelationTask localCorrelationTask(activeSystemController, allSystems);
    localCorrelationTask.settings.lcSettings = this->settings.lcSettings;
    localCorrelationTask.settings.normThreshold = this->settings.normThreshold;
    localCorrelationTask.settings.maxCycles = this->settings.maxCycles;
    localCorrelationTask.settings.writePairEnergies = this->settings.writePairEnergies;
    localCorrelationTask.settings.loc = this->settings.loc;
    localCorrelationTask.run();
    correlationEnergyContributions(iAct) =
        localCorrelationTask.getCorrelationEnergy(activeSystemController, this->settings.lcSettings.method);
  }
  return correlationEnergyContributions;
}

template<Options::SCF_MODES SCFMode>
void TopDownStaticEmbeddingTask<SCFMode>::setUpSubsystems(
    std::shared_ptr<SystemController> supersystem, std::vector<std::shared_ptr<SystemController>> activeSystems,
    std::vector<std::shared_ptr<SystemController>> environmentSystems, const LocalizationTaskSettings& loc,
    const SystemSplittingTaskSettings& split, const SystemAdditionTaskSettings& add,
    const BasisSetTruncationTaskSettings& trunc, bool useQuasiRestrictedOrbitals) {
  std::vector<std::shared_ptr<SystemController>> allSystems = environmentSystems;
  allSystems.insert(allSystems.end(), activeSystems.begin(), activeSystems.end());
  if (!supersystem) {
    if (allSystems.size() < 1)
      throw SerenityError("ERROR: At least one system must be available in order to built a supersystem.\n");
    Settings supersystemSettings = allSystems[0]->getSettings();
    supersystemSettings.name = "TMP_Supersystem";
    supersystemSettings.path = supersystemSettings.path + supersystemSettings.name + "/";
    supersystemSettings.charge = 0;
    supersystemSettings.spin = 0;
    supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersystemSettings);
    // Addition
    SystemAdditionTask<SCFMode> additionTask(supersystem, allSystems);
    additionTask.settings = add;
    additionTask.run();
    // Scf
    printSubSectionTitle("Supersystem SCF Calculation");
    ScfTask<SCFMode> scfTask(supersystem);
    scfTask.run();
  }
  else {
    printSubSectionTitle("Supersystem Energy Evaluation");
    ScfTask<SCFMode> scfTask(supersystem);
    scfTask.settings.skipSCF = true;
    scfTask.run();
  }
  // Quasi restricted orbitals
  if (useQuasiRestrictedOrbitals) {
    QuasiRestrictedOrbitalsTask<SCFMode> quasiRestrictedOrbitalsTask(supersystem, {});
    quasiRestrictedOrbitalsTask.run();
  }
  // Localization
  LocalizationTask localizationTask(supersystem);
  localizationTask.settings = loc;
  localizationTask.settings.separateSOMOs = useQuasiRestrictedOrbitals || loc.separateSOMOs;
  localizationTask.run();

  // Splitting
  SystemSplittingTask<SCFMode> splittingTask(supersystem, allSystems);
  splittingTask.settings = split;
  splittingTask.run();
  // Basis set truncation
  if (trunc.truncAlgorithm != Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE) {
    for (auto& activeSystem : activeSystems) {
      BasisSetTruncationTask<SCFMode> truncationTask(activeSystem);
      truncationTask.settings = trunc;
      truncationTask.run();
    }
  }
  cacheHCoreIntegrals(supersystem, activeSystems, environmentSystems);
  supersystem->setDiskMode(true);
}

template<Options::SCF_MODES SCFMode>
double TopDownStaticEmbeddingTask<SCFMode>::getFinalEnergy() {
  if (!_finalTotalEnergy) {
    this->run();
  }
  return *_finalTotalEnergy;
}

template<Options::SCF_MODES SCFMode>
void TopDownStaticEmbeddingTask<SCFMode>::cacheHCoreIntegrals(std::shared_ptr<SystemController> supersystem,
                                                              std::vector<std::shared_ptr<SystemController>> activeSystems,
                                                              std::vector<std::shared_ptr<SystemController>> environmentSystems) {
  const auto& extCharges = supersystem->getOneElectronIntegralController()->getExternalCharges();
  if (!extCharges.empty()) {
    auto superBasis = supersystem->getBasisController();
    const auto& extChargeIntegrals = supersystem->getOneElectronIntegralController()->getExtChargeIntegrals();
    std::vector<std::shared_ptr<SystemController>> allSubsystems = activeSystems;
    allSubsystems.insert(allSubsystems.end(), environmentSystems.begin(), environmentSystems.end());
    BasisFunctionMapper basisFunctionMapper(superBasis);
    for (const auto& subsystem : allSubsystems) {
      auto subBasis = subsystem->getBasisController();
      const auto projectionMatrix = *basisFunctionMapper.getSparseProjection(subBasis);
      MatrixInBasis<RESTRICTED> reshuffledIntegrals(subBasis);
      unsigned int nBasisFunctions = reshuffledIntegrals.rows();
      reshuffledIntegrals.block(0, 0, nBasisFunctions, nBasisFunctions) =
          projectionMatrix * extChargeIntegrals * projectionMatrix.transpose();
      subsystem->getOneElectronIntegralController()->cacheExtChargeIntegrals(reshuffledIntegrals);
    }
  }
}

template class TopDownStaticEmbeddingTask<Options::SCF_MODES::RESTRICTED>;
template class TopDownStaticEmbeddingTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
