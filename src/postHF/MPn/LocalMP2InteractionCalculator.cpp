/**
 * @file LocalMP2InteractionCalculator.cpp
 *
 * @date 18 Apr 2020
 * @author Moritz Bensberg
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
#include "postHF/MPn/LocalMP2InteractionCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/BasisFunctionMapper.h"
#include "data/OrbitalController.h"
#include "dft/Functional.h"                       //Functional definition (getHfCorrelRatio() etc.).
#include "dft/functionals/CompositeFunctionals.h" //Resolve functional
#include "geometry/Geometry.h"                    //Supersystem construction.
#include "io/FormattedOutputStream.h"
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h" //K-Set definition.
#include "postHF/MPn/LocalMP2.h"
#include "settings/Settings.h" //Check electronic structure method of the system.
#include "system/SystemController.h"
#include "tasks/SystemAdditionTask.h"

namespace Serenity {

LocalMP2InteractionCalculator::LocalMP2InteractionCalculator(std::shared_ptr<SystemController> activeSystem,
                                                             std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                                             LocalCorrelationSettings lcSettings, unsigned int maxCycles,
                                                             double maxResidual, bool fullCoupling,
                                                             std::shared_ptr<SystemController> supersystem,
                                                             std::vector<unsigned int> environmentOrbitalIndices)
  : _activeSystem(activeSystem),
    _environmentSystems(environmentSystems),
    _lcSettings(lcSettings),
    _maxCycles(maxCycles),
    _maxResidual(maxResidual),
    _fullCoupling(fullCoupling),
    _supersystem(supersystem),
    _environmentOrbitalIndices(environmentOrbitalIndices) {
  _lcSettings.resolvePNOSettings();
}

double LocalMP2InteractionCalculator::getInteractionEnergy() {
  if (not _energiesAvailable)
    calculateLocalMP2Energies();
  return _interactionEnergy;
}

double LocalMP2InteractionCalculator::getEnvironmentEnergy() {
  if (not _energiesAvailable)
    calculateLocalMP2Energies();
  return _environmentEnergy;
}

double LocalMP2InteractionCalculator::getActiveEnergy() {
  if (not _energiesAvailable)
    calculateLocalMP2Energies();
  return _activeEnergy;
}

double LocalMP2InteractionCalculator::getCouplingEnergyCorrection() {
  if (not _energiesAvailable)
    calculateLocalMP2Energies();
  return _couplingCorrectionEnergy;
}

void LocalMP2InteractionCalculator::setUpSupersystem() {
  Settings supersystemSettings = _environmentSystems[0]->getSettings();
  supersystemSettings.dft.functional = _lcSettings.embeddingSettings.naddXCFunc;
  supersystemSettings.name = "TMP_Supersystem";
  supersystemSettings.charge = 0;
  supersystemSettings.spin = 0;
  _supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersystemSettings);
  std::vector<std::shared_ptr<SystemController>> allFragments = _environmentSystems;
  allFragments.insert(allFragments.begin(), _activeSystem);
  SystemAdditionTask<RESTRICTED> additionTask(_supersystem, allFragments);
  additionTask.settings.addOccupiedOrbitals = true;
  additionTask.run();
  _environmentOrbitalIndices = {};
  unsigned int nOccSup = _supersystem->getNOccupiedOrbitals<RESTRICTED>();
  unsigned int nOccAct = _activeSystem->getNOccupiedOrbitals<RESTRICTED>();
  for (unsigned int i = nOccAct; i < nOccSup; ++i)
    _environmentOrbitalIndices.push_back(i);
}

void LocalMP2InteractionCalculator::buildOrbitalPairs() {
  takeTime("Initial orbital pairs");
  unsigned int nOccSup = _supersystem->getNOccupiedOrbitals<RESTRICTED>();
  unsigned int nOccEnv = 0;
  for (auto env : _environmentSystems)
    nOccEnv += env->getNOccupiedOrbitals<RESTRICTED>();
  if (_environmentOrbitalIndices.size() != nOccEnv)
    throw SerenityError("Not all environment orbitals are contained in the supersystem. Something went wrong.");
  if (_lcSettings.useFrozenCore) {
    OutputControl::nOut << "  The frozen core approximation is used!\n";
    OutputControl::nOut << "  Please ensure that core orbitals have been localized independently!\n";
    OutputControl::nOut << "  This may lead to artifacts in the core-orbital selection otherwise.\n";
  }
  const auto& coreOrbitals = _supersystem->getActiveOrbitalController<RESTRICTED>()->getOrbitalFlags();
  const double pnoCoreThreshold = _lcSettings.pnoThreshold * _lcSettings.pnoCoreScaling;
  for (unsigned int iOcc = 0; iOcc < nOccSup; ++iOcc) {
    bool iIsCore = coreOrbitals(iOcc);
    if (_lcSettings.useFrozenCore && iIsCore)
      continue;
    bool iOccIsEnv = std::find(_environmentOrbitalIndices.begin(), _environmentOrbitalIndices.end(), iOcc) !=
                     _environmentOrbitalIndices.end();
    for (unsigned int jOcc = 0; jOcc <= iOcc; ++jOcc) {
      bool jIsCore = coreOrbitals(jOcc);
      if (_lcSettings.useFrozenCore && jIsCore)
        continue;
      auto newPair =
          std::make_shared<OrbitalPair>(iOcc, jOcc, (iIsCore || jIsCore) ? pnoCoreThreshold : _lcSettings.pnoThreshold,
                                        _lcSettings.ccsdPairThreshold, _lcSettings.collinearDipoleScaling);
      bool jOccIsEnv = std::find(_environmentOrbitalIndices.begin(), _environmentOrbitalIndices.end(), jOcc) !=
                       _environmentOrbitalIndices.end();
      if (iOccIsEnv and jOccIsEnv) {
        // Env--Env pair.
        _env_envPairs.push_back(newPair);
      }
      else if (iOccIsEnv or jOccIsEnv) {
        // Env--Act pair
        _env_actPairs.push_back(newPair);
      }
      else if (_fullCoupling) {
        // Act--Act pair
        _act_actPairs.push_back(newPair);
      }
    } // for jOcc
  }   // for iOcc
  timeTaken(2, "Initial orbital pairs");
}

void LocalMP2InteractionCalculator::restrictCouplings(std::vector<std::shared_ptr<OrbitalPair>> pairs) {
  for (auto& pair : pairs) {
    std::vector<std::shared_ptr<CouplingOrbitalSet>> kSets = pair->coupledPairs;
    std::vector<std::shared_ptr<CouplingOrbitalSet>> newKSets;
    for (auto kSet : kSets) {
      unsigned int k = kSet->getK();
      bool kIsEnv = std::find(_environmentOrbitalIndices.begin(), _environmentOrbitalIndices.end(), k) !=
                    _environmentOrbitalIndices.end();
      if (not kIsEnv)
        newKSets.push_back(kSet);
    } // for kSet
    pair->coupledPairs = newKSets;
  } // for pair
}

double LocalMP2InteractionCalculator::calculateActiveOnlyLocalMP2Energy(std::shared_ptr<LocalCorrelationController> lcController) {
  if (_lcSettings.embeddingSettings.embeddingMode == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
    LocalMP2 localMP2(lcController);
    localMP2.settings.maxResidual = _maxResidual;
    localMP2.settings.maxCycles = _maxCycles;
    if (_supersystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      auto functional = _supersystem->getSettings().customFunc.basicFunctionals.size()
                            ? Functional(_supersystem->getSettings().customFunc)
                            : resolveFunctional(_supersystem->getSettings().dft.functional);
      localMP2.settings.osScaling = functional.getosScaling();
      localMP2.settings.ssScaling = functional.getssScaling();
    }
    restrictCouplings(_act_actPairs);
    localMP2.calculateEnergyCorrection(_act_actPairs);
    double activeOnlyLocalMP2Energy = 0.0;
    for (auto pair : _act_actPairs) {
      activeOnlyLocalMP2Energy += getPairEnergy(pair);
    } // for pair
    return activeOnlyLocalMP2Energy;
  }
  else {
    Settings tmpActiveSettings = _activeSystem->getSettings();
    tmpActiveSettings.dft.functional = _lcSettings.embeddingSettings.naddXCFunc;
    tmpActiveSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
    tmpActiveSettings.name = "TMP_Active";
    tmpActiveSettings.load = "";
    tmpActiveSettings.charge = 0;
    tmpActiveSettings.spin = 0;
    auto tmpActiveSystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), tmpActiveSettings);
    SystemAdditionTask<RESTRICTED> additionTask(tmpActiveSystem, {_activeSystem});
    additionTask.settings.addOccupiedOrbitals = true;
    additionTask.run();
    auto localCorrelationController =
        std::make_shared<LocalCorrelationController>(tmpActiveSystem, lcController->getSettings(), _environmentSystems);
    LocalMP2 localMP2(localCorrelationController);
    localMP2.settings.maxResidual = _maxResidual;
    localMP2.settings.maxCycles = _maxCycles;

    auto functional = tmpActiveSystem->getSettings().customFunc.basicFunctionals.size()
                          ? Functional(tmpActiveSystem->getSettings().customFunc)
                          : resolveFunctional(tmpActiveSystem->getSettings().dft.functional);
    localMP2.settings.osScaling = functional.getosScaling();
    localMP2.settings.ssScaling = functional.getssScaling();
    return localMP2.calculateEnergyCorrection().sum();
  }
}

std::shared_ptr<LocalCorrelationController> LocalMP2InteractionCalculator::runLocalMP2() {
  std::vector<std::shared_ptr<OrbitalPair>> allPairs = _env_actPairs;
  allPairs.insert(allPairs.begin(), _env_envPairs.begin(), _env_envPairs.end());
  allPairs.insert(allPairs.begin(), _act_actPairs.begin(), _act_actPairs.end());
  std::vector<std::shared_ptr<SystemController>> dummy = {};
  auto localCorrelationController =
      std::make_shared<LocalCorrelationController>(_supersystem, _lcSettings, dummy, nullptr, allPairs);
  LocalMP2 localMP2(localCorrelationController);
  localMP2.settings.maxResidual = _maxResidual;
  localMP2.settings.maxCycles = _maxCycles;
  if (_supersystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    auto functional = _supersystem->getSettings().customFunc.basicFunctionals.size()
                          ? Functional(_supersystem->getSettings().customFunc)
                          : resolveFunctional(_supersystem->getSettings().dft.functional);
    localMP2.settings.osScaling = functional.getosScaling();
    localMP2.settings.ssScaling = functional.getssScaling();
  }
  localMP2.calculateEnergyCorrection();
  return localCorrelationController;
}

double LocalMP2InteractionCalculator::getPairEnergy(std::shared_ptr<OrbitalPair> pair) {
  double pairEnergy = 0.0;
  if (pair->type == OrbitalPairTypes::CLOSE) {
    pairEnergy = pair->lMP2PairEnergy;
    pairEnergy += pair->deltaPNO;
  }
  else {
    if (pair->scMP2PairEnergy != 0.0) {
      pairEnergy = pair->scMP2PairEnergy;
    }
    else {
      pairEnergy = pair->dipolePairEnergy;
    }
  }
  return pairEnergy;
}

void LocalMP2InteractionCalculator::calculateLocalMP2Energies() {
  // Build supersystem. The orbital ordering will be active ... environment.
  if (not _supersystem)
    setUpSupersystem();
  // Construct orbital pairs for the supersystem.
  buildOrbitalPairs();
  OutputControl::vOut << "Total number of pairs: " << _env_actPairs.size() + _env_envPairs.size() + _act_actPairs.size()
                      << std::endl;
  OutputControl::vOut << "# act-env " << _env_actPairs.size() << std::endl;
  OutputControl::vOut << "# env-env " << _env_envPairs.size() << std::endl;
  OutputControl::vOut << "# act-act " << _act_actPairs.size() << std::endl;
  // Calculate MP2 amplitudes and pair energies which will be stored in the pairs.
  auto lcController = runLocalMP2();
  // Sum pair energies up.
  _interactionEnergy = 0.0;
  for (auto pair : _env_actPairs) {
    _interactionEnergy += getPairEnergy(pair);
    pair->cleanUp();
  }
  _environmentEnergy = 0.0;
  for (auto pair : _env_envPairs) {
    _environmentEnergy += getPairEnergy(pair);
    pair->cleanUp();
  }
  _activeEnergy = 0.0;
  for (auto pair : _act_actPairs) {
    _activeEnergy += getPairEnergy(pair);
  }

  double activeOnlyMP2Energy = 0.0;
  if (_fullCoupling)
    activeOnlyMP2Energy = calculateActiveOnlyLocalMP2Energy(lcController);
  for (auto pair : _act_actPairs) {
    pair->cleanUp();
  }
  _couplingCorrectionEnergy = _activeEnergy - activeOnlyMP2Energy;
  // Scaling for DFT/double hybrid functionals.
  if (_supersystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    auto functional = _supersystem->getSettings().customFunc.basicFunctionals.size()
                          ? Functional(_supersystem->getSettings().customFunc)
                          : resolveFunctional(_supersystem->getSettings().dft.functional);
    _interactionEnergy *= functional.getHfCorrelRatio();
    _environmentEnergy *= functional.getHfCorrelRatio();
    _activeEnergy *= functional.getHfCorrelRatio();
    _couplingCorrectionEnergy *= functional.getHfCorrelRatio();
  }
  _energiesAvailable = true;
}

LocalMP2InteractionCalculator::~LocalMP2InteractionCalculator() {
  removeSystemFiles(_activeSystem->getSystemPath() + "TMP_Active/", "TMP_Active");
  removeSystemFiles(_environmentSystems[0]->getSystemPath() + "TMP_Supersystem/", "TMP_Supersystem");
}

} /* namespace Serenity */
