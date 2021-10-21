/**
 * @file DOSCCTask.cpp
 *
 * @author Moritz Bensberg
 * @date Jul 8, 2021
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
#include "tasks/DOSCCTask.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"        //Dummy constructor for the geometry.
#include "io/FormattedOutput.h"       //Print nice section titles.
#include "io/FormattedOutputStream.h" //Filtered output.
#include "settings/Settings.h"        //New system settings.
#include "system/SystemController.h"  //Construct dummy systems.
/* Include Std and External Headers */
#include <iomanip> //setw(...)

namespace Serenity {

DOSCCTask::DOSCCTask(std::vector<std::shared_ptr<SystemController>> supersystems) : _supersystems(supersystems) {
  settings.wfemb.lcSettings = std::vector<LocalCorrelationSettings>(_maxFragments);
}

DOSCCTask::~DOSCCTask() = default;

void DOSCCTask::run() {
  // Resolve DOS settings.
  settings.resolveDOSSettings();
  _nFragments = settings.gdos.similarityLocThreshold.size() + 1;
  // Resize the settings vector the number of desired fragments.
  settings.wfemb.lcSettings.resize(_nFragments);
  // Build subsystems
  setUpSubsystems();
  // Run the orbital localization for each supersystem.
  localizeOrbitals();
  // Run generalized DOS.
  runGDOS();
  // Run multi-level CC
  _totalEnergies = std::make_unique<Eigen::VectorXd>(runCC());
  // Print the final result / energy differences with respect to the 0-th energy.
  printResults(*_totalEnergies);
}

Eigen::VectorXd DOSCCTask::runCC() {
  Eigen::VectorXd energies(_supersystems.size());
  for (unsigned int i = 0; i < _supersystems.size(); ++i) {
    auto supersystem = _supersystems[i];
    auto subsystemSet = _subSystems[i];
    WavefunctionEmbeddingTask wfemb(supersystem, subsystemSet);
    wfemb.settings = settings.wfemb;
    wfemb.settings.fromFragments = true;
    wfemb.run();
    energies(i) = wfemb.getTotalEnergy();
  } // for i
  return energies;
}

void DOSCCTask::runGDOS() {
  std::vector<std::shared_ptr<SystemController>> allSubsystems;
  for (const auto& subsystemSet : _subSystems) {
    allSubsystems.insert(allSubsystems.end(), subsystemSet.begin(), subsystemSet.end());
  }
  GeneralizedDOSTask<RESTRICTED> gdos(_supersystems, allSubsystems);
  gdos.settings = settings.gdos;
  gdos.run();
}

void DOSCCTask::localizeOrbitals() {
  // Run localization for the first point
  LocalizationTask locOne(_supersystems[0]);
  locOne.settings = settings.wfemb.loc;
  locOne.run();
  // Run localization for all other points + optional alignment
  // with respect to the first point.
  for (unsigned int i = 1; i < _supersystems.size(); ++i) {
    auto& system = _supersystems[i];
    if (settings.alignOrbitals) {
      auto reference = _supersystems[0];
      LocalizationTask align(system, {reference});
      align.settings = settings.wfemb.loc;
      align.settings.maxSweeps = 10;
      align.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN;
      align.run();
    }
    LocalizationTask locI(system);
    locI.settings = settings.wfemb.loc;
    locI.run();
  } // for i
}

void DOSCCTask::setUpSubsystems() {
  for (const auto& sys : _supersystems) {
    std::vector<std::shared_ptr<SystemController>> newSubsystemSet;
    for (unsigned int i = 0; i < _nFragments; ++i) {
      Settings newSettings = sys->getSettings();
      newSettings.name = newSettings.name + std::to_string(i);
      newSettings.identifier = newSettings.identifier + std::to_string(i);
      newSubsystemSet.push_back(std::make_shared<SystemController>(std::make_shared<Geometry>(), newSettings));
    } // for i
    _subSystems.push_back(newSubsystemSet);
  } // for sys
}

void DOSCCTask::printResults(const Eigen::VectorXd& energies) {
  double zeroPoint = energies(0);
  _relativeEnergies = std::make_unique<Eigen::VectorXd>(Eigen::VectorXd::Zero(_supersystems.size()));
  printSubSectionTitle("DOS-Coupled-Cluster Results");
  OutputControl::mOut << "  Name                          Total Energy     Relative Energy" << std::endl;
  OutputControl::mOut << "----------------------------------------------------------------" << std::endl;
  for (unsigned int i = 0; i < _supersystems.size(); ++i) {
    double relativeEnergy = energies[i] - zeroPoint;
    (*_relativeEnergies)(i) = relativeEnergy;
    OutputControl::mOut << "  " << std::left << std::setw(22) << std::setprecision(10)
                        << _supersystems[i]->getSystemName() << std::right << std::setw(20) << energies[i]
                        << std::setw(20) << std::setprecision(6) << relativeEnergy << std::endl;
  }
}

const Eigen::VectorXd& DOSCCTask::getRelativeEnergies() {
  if (!_relativeEnergies)
    this->run();
  return *_relativeEnergies;
}

const Eigen::VectorXd& DOSCCTask::getTotalEnergies() {
  if (!_totalEnergies)
    this->run();
  return *_totalEnergies;
}

void DOSCCTaskSettings::resolveDOSSettings() {
  const std::map<Options::DOS_SETTINGS, std::vector<double>> m = {{Options::DOS_SETTINGS::LOOSE, {1e-1, 1e-2, 1e-3}},
                                                                  {Options::DOS_SETTINGS::NORMAL, {5e-2, 5e-3, 1e-3}},
                                                                  {Options::DOS_SETTINGS::TIGHT, {2e-2, 2e-3, 8e-4}},
                                                                  {Options::DOS_SETTINGS::VERY_TIGHT, {8e-3, 1e-3, 1e-4}},
                                                                  {Options::DOS_SETTINGS::EXTREME, {5e-3, 5e-4, 1e-4}}};

  if (gdos.similarityKinEnergyThreshold.size() == 0) {
    gdos.similarityKinEnergyThreshold = m.at(dosSettings);
  }
  if (gdos.similarityLocThreshold.size() == 0) {
    gdos.similarityLocThreshold = m.at(dosSettings);
  }
}

} /* namespace Serenity */
