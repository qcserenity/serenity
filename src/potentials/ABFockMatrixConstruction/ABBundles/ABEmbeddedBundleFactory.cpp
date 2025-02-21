/**
 * @file ABEmbeddedBundleFactory.cpp
 *
 * @date 18 Aug 2019
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
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundleFactory.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h" //Mulliken populations/Hybrid approach.
#include "basis/AtomCenteredBasisController.h"                        //Combined AB aux. basis sets.
#include "basis/AtomCenteredBasisControllerFactory.h"                 //Combined AB aux. basis sets.
#include "data/ElectronicStructure.h"                                 //Check for electronic structure.
#include "data/OrbitalController.h"                                   //Hybrid approach.
#include "dft/functionals/CompositeFunctionals.h"                     //Resolving of settings.
#include "geometry/Geometry.h"                                        //Combined AB geometries.
#include "grid/AtomCenteredGridControllerFactory.h"                   //Combinded AB grids.
#include "integrals/CDIntegralController.h"                           //Generate AB aCD-auxiliary basis
#include "integrals/OneElectronIntegralController.h"                  //Mulliken populations/Hybrid approach.
#include "misc/SystemSplittingTools.h" //Environment density matrix controller construction.
#include "potentials/ABFockMatrixConstruction/ABCoreHamiltonian.h"
#include "potentials/ABFockMatrixConstruction/ABCoulombInteractionPotential.h"
#include "potentials/ABFockMatrixConstruction/ABERIPotential.h"
#include "potentials/ABFockMatrixConstruction/ABExchangePotential.h"
#include "potentials/ABFockMatrixConstruction/ABFuncPotential.h"
#include "potentials/ABFockMatrixConstruction/ABNAddFuncPotential.h"
#include "potentials/ABFockMatrixConstruction/ABZeroPotential.h"
#include "settings/Settings.h"       //Settings.
#include "system/SystemController.h" //Definition of a SystemController.

namespace Serenity {

template<>
std::shared_ptr<ABEmbeddedBundle<RESTRICTED>>
ABEmbeddedBundleFactory<RESTRICTED>::produce(std::shared_ptr<SystemController> activeSystem,
                                             std::shared_ptr<BasisController> basisControllerB,
                                             std::shared_ptr<Geometry> geometryB,
                                             std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                             const std::shared_ptr<EmbeddingSettings> embeddingSettings, bool topDown) {
  if (!_restrictedFactory)
    _restrictedFactory.reset(new ABEmbeddedBundleFactory<RESTRICTED>);
  return _restrictedFactory->getOrProduce(activeSystem, basisControllerB, geometryB, environmentSystems,
                                          embeddingSettings, topDown);
}

template<>
std::shared_ptr<ABEmbeddedBundle<UNRESTRICTED>>
ABEmbeddedBundleFactory<UNRESTRICTED>::produce(std::shared_ptr<SystemController> activeSystem,
                                               std::shared_ptr<BasisController> basisControllerB,
                                               std::shared_ptr<Geometry> geometryB,
                                               std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                               const std::shared_ptr<EmbeddingSettings> embeddingSettings, bool topDown) {
  if (!_unrestrictedFactory)
    _unrestrictedFactory.reset(new ABEmbeddedBundleFactory<UNRESTRICTED>);
  return _unrestrictedFactory->getOrProduce(activeSystem, basisControllerB, geometryB, environmentSystems,
                                            embeddingSettings, topDown);
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>>
ABEmbeddedBundleFactory<SCFMode>::getNotProjectedEnvironmentDensityMatrixControllers(
    std::shared_ptr<SystemController> activeSystem, std::vector<std::shared_ptr<SystemController>> environmentSystems,
    double basisFunctionRatio, double borderAtomThreshold) {
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> notProjectedEnvDensities;
  for (auto sys : environmentSystems) {
    const auto& envCoeff = sys->getActiveOrbitalController<SCFMode>()->getCoefficients();
    auto orbitalPopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
        envCoeff, sys->getOneElectronIntegralController()->getOverlapIntegrals(),
        sys->getAtomCenteredBasisController()->getBasisIndices());
    auto distantOrbitals = SystemSplittingTools<SCFMode>::selectDistantOrbitals(orbitalPopulations, activeSystem, sys,
                                                                                basisFunctionRatio, borderAtomThreshold);
    auto nonOrthogonalDensityMatrix = SystemSplittingTools<SCFMode>::buildNonOrthogonalDensityMatrix(sys, distantOrbitals);
    notProjectedEnvDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(*nonOrthogonalDensityMatrix));
  }
  return notProjectedEnvDensities;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<BasisController>
ABEmbeddedBundleFactory<SCFMode>::getABAuxBasisController(std::shared_ptr<SystemController> activeSystem,
                                                          std::shared_ptr<Geometry> geometryB) {
  auto densFitJ = activeSystem->getSettings().basis.densFitJ;
  auto combinedGeometry = std::make_shared<Geometry>();
  auto actGeom = activeSystem->getGeometry();
  *combinedGeometry += *actGeom;
  *combinedGeometry += *geometryB;
  combinedGeometry->deleteIdenticalAtoms();
  std::shared_ptr<BasisController> abAuxBasisController;
  if (densFitJ == Options::DENS_FITS::NONE) {
    throw SerenityError("An AB auxiliary basis controller was requested but no densityfitting was specified.");
  }
  else if (densFitJ == Options::DENS_FITS::RI) {
    abAuxBasisController = AtomCenteredBasisControllerFactory::produce(
        combinedGeometry, activeSystem->getSettings().basis.basisLibPath, activeSystem->getSettings().basis.makeSphericalBasis,
        false, activeSystem->getSettings().basis.firstECP, activeSystem->getSettings().basis.auxJLabel);
  }
  else if (densFitJ == Options::DENS_FITS::ACD) {
    activeSystem->getCDIntegralController()->generateACDBasis(combinedGeometry, activeSystem->getSettings().basis.label, "-AB");
    std::string label = "ACD-" + activeSystem->getSettings().basis.label + "-AB";
    abAuxBasisController = AtomCenteredBasisControllerFactory::produce(
        combinedGeometry, activeSystem->getSettings().path, activeSystem->getSettings().basis.makeSphericalBasis, false,
        activeSystem->getSettings().basis.firstECP, label);
  }
  else if (densFitJ == Options::DENS_FITS::ACCD) {
    activeSystem->getCDIntegralController()->generateACCDBasis(combinedGeometry, activeSystem->getSettings().basis.label, "-AB");
    std::string label = "ACCD-" + activeSystem->getSettings().basis.label + "-AB";
    abAuxBasisController = AtomCenteredBasisControllerFactory::produce(
        combinedGeometry, activeSystem->getSettings().path, activeSystem->getSettings().basis.makeSphericalBasis, false,
        activeSystem->getSettings().basis.firstECP, label);
  }

  return abAuxBasisController;
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<ABEmbeddedBundle<SCFMode>>
ABEmbeddedBundleFactory<SCFMode>::produceNew(std::shared_ptr<SystemController> activeSystem,
                                             std::shared_ptr<BasisController> basisControllerB,
                                             std::shared_ptr<Geometry> geometryB,
                                             std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                             const std::shared_ptr<EmbeddingSettings> embeddingSettings, bool topDown) {
  auto supersystemGeometry = std::make_shared<Geometry>();
  auto actGeom = activeSystem->getGeometry();
  *supersystemGeometry += *actGeom;
  for (unsigned int i = 0; i < environmentSystems.size(); ++i) {
    *supersystemGeometry += *environmentSystems[i]->getGeometry();
  }
  supersystemGeometry->deleteIdenticalAtoms();

  auto environmentDensityControllers = SystemSplittingTools<SCFMode>::getEnvironmentDensityControllers(environmentSystems);
  std::vector<std::shared_ptr<BasisController>> envAuxBasis;
  if (activeSystem->getSettings().basis.densFitJ != Options::DENS_FITS::NONE) {
    for (const auto& envSys : environmentSystems) {
      envAuxBasis.push_back(envSys->getAuxBasisController(Options::AUX_BASIS_PURPOSES::COULOMB,
                                                          activeSystem->getSettings().basis.densFitJ));
    }
  }
  else {
    for (unsigned int i = 0; i < environmentSystems.size(); ++i) {
      envAuxBasis.push_back(nullptr);
    }
  }

  auto basisContA = activeSystem->getBasisController();
  auto hcore = std::make_shared<ABCoreHamiltonian<SCFMode>>(basisContA, basisControllerB, supersystemGeometry);

  auto combinedGeometry = std::make_shared<Geometry>();
  *combinedGeometry += *activeSystem->getGeometry();
  *combinedGeometry += *geometryB;
  combinedGeometry->deleteIdenticalAtoms();
  std::shared_ptr<Serenity::BasisController> abAuxBasisController = nullptr;
  if (activeSystem->getSettings().basis.densFitJ != Options::DENS_FITS::NONE) {
    abAuxBasisController = getABAuxBasisController(activeSystem, geometryB);
  }
  const auto activeFunctional = activeSystem->getSettings().customFunc.basicFunctionals.size()
                                    ? Functional(activeSystem->getSettings().customFunc)
                                    : resolveFunctional(activeSystem->getSettings().dft.functional);
  double exchangeRatioActive = (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF)
                                   ? 1.0
                                   : activeFunctional.getHfExchangeRatio();
  double rangeSeperationParameterActive = activeFunctional.getRangeSeparationParameter();
  double lrExchangeRatioActive = activeFunctional.getLRExchangeRatio();
  std::vector<std::shared_ptr<BasisController>> actAuxVec = {abAuxBasisController};
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> actDensityMatrixController = {
      activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()};
  auto activeCoulomb = std::make_shared<ABERIPotential<SCFMode>>(
      activeSystem, basisContA, basisControllerB, actDensityMatrixController, exchangeRatioActive, lrExchangeRatioActive,
      rangeSeperationParameterActive, topDown, activeSystem->getSettings().basis.densFitJ, abAuxBasisController, actAuxVec);
  const auto naddXCFunc = embeddingSettings->customNaddXCFunc.basicFunctionals.size()
                              ? Functional(embeddingSettings->customNaddXCFunc)
                              : resolveFunctional(embeddingSettings->naddXCFunc);
  double exchangeRatioNAdd = naddXCFunc.getHfExchangeRatio();
  double rangeSeperationParameterNAdd = naddXCFunc.getRangeSeparationParameter();
  double lrExchangeRatioNAdd = naddXCFunc.getLRExchangeRatio();
  auto environmentCoulomb = std::make_shared<ABERIPotential<SCFMode>>(
      activeSystem, basisContA, basisControllerB, environmentDensityControllers, exchangeRatioNAdd, lrExchangeRatioNAdd,
      rangeSeperationParameterNAdd, topDown, activeSystem->getSettings().basis.densFitJ, abAuxBasisController, envAuxBasis);
  // Build A+B grid controller
  std::shared_ptr<GridController> grid_AB =
      AtomCenteredGridControllerFactory::produce(combinedGeometry, activeSystem->getSettings().grid);
  std::shared_ptr<ABPotential<SCFMode>> activeExchangeCorrelation;
  if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    auto actDensityMatrixController = {activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()};
    activeExchangeCorrelation = std::make_shared<ABFuncPotential<SCFMode>>(
        activeSystem, basisContA, basisControllerB, grid_AB, actDensityMatrixController, activeFunctional);
  }
  else if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    activeExchangeCorrelation = std::make_shared<ABZeroPotential<SCFMode>>(basisContA, basisControllerB);
  }
  else {
    assert(false && "Unknwon electronic structure theory!");
  }
  auto naddExchangeCorrelation = std::make_shared<ABNAddFuncPotential<SCFMode>>(
      activeSystem, basisContA, basisControllerB, environmentDensityControllers, grid_AB, naddXCFunc);
  // build non-additive kinetic energy contribution in case of an hybrid method
  std::shared_ptr<ABPotential<SCFMode>> naddKinetic = nullptr;
  if (embeddingSettings->longRangeNaddKinFunc != CompositeFunctionals::KINFUNCTIONALS::NONE) {
    auto notProjectedEnvDensMatCont = getNotProjectedEnvironmentDensityMatrixControllers(
        activeSystem, environmentSystems, embeddingSettings->basisFunctionRatio, embeddingSettings->borderAtomThreshold);
    naddKinetic = std::make_shared<ABNAddFuncPotential<SCFMode>>(
        activeSystem, basisContA, basisControllerB, notProjectedEnvDensMatCont, grid_AB,
        embeddingSettings->customNaddKinFunc.basicFunctionals.size() ? Functional(embeddingSettings->customNaddKinFunc)
                                                                     : resolveFunctional(embeddingSettings->naddKinFunc));
  }
  return std::unique_ptr<ABEmbeddedBundle<SCFMode>>(new ABEmbeddedBundle<SCFMode>(
      hcore, activeCoulomb, environmentCoulomb, activeExchangeCorrelation, naddExchangeCorrelation, naddKinetic));
}
template class ABEmbeddedBundleFactory<Options::SCF_MODES::RESTRICTED>;
template class ABEmbeddedBundleFactory<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
