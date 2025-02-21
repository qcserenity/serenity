/**
 * @file HuzinagaProjectionPotential.cpp
 *
 * @date Nov 23, 2017
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
#include "potentials/HuzinagaProjectionPotential.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "basis/BasisFunctionMapper.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/SPMatrix.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/IOOptions.h"
#include "misc/SystemSplittingTools.h"
#include "misc/Timing.h"
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundleFactory.h"
#include "settings/EmbeddingSettings.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HuzinagaProjectionPotential<SCFMode>::HuzinagaProjectionPotential(
    std::shared_ptr<SystemController> activeSystem, std::vector<std::shared_ptr<SystemController>> environmentSystems,
    const EmbeddingSettings& settings, std::shared_ptr<PotentialBundle<SCFMode>> activeFockMatrix, bool topDown,
    std::shared_ptr<GridController> supersystemgrid, double gridCutOff,
    std::vector<std::shared_ptr<EnergyComponentController>> allEConts, double fermiShift)
  : Potential<SCFMode>(activeSystem->getBasisController()),
    _activeSystem(activeSystem),
    _activeFockMatrix(activeFockMatrix),
    _fermiShift(fermiShift) {
  takeTime("Huzinaga -- pre-calculations");
  for (auto e : environmentSystems)
    _environmentSystems.push_back(e);
  // Expects that an electronic structure is already in place for the active system!
  activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->addSensitiveObject(
      ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  // A projection does only make sense if there is something to project
  assert(environmentSystems.size() > 0 && "There is no environment system given!");
  _s_ABs = SystemSplittingTools<SCFMode>::getProjectedSubsystems(
      activeSystem, environmentSystems, (settings.truncateProjector) ? settings.projecTruncThresh : -10.0);
  _envDensityCont = SystemSplittingTools<SCFMode>::getEnvironmentDensityControllers(environmentSystems, topDown);
  for (auto sys : environmentSystems) {
    // If it is a hybrid within each environment system. Construct the not projected
    // density matrix.
    if (settings.longRangeNaddKinFunc != CompositeFunctionals::KINFUNCTIONALS::NONE) {
      // Localization (already done for top-down calculations)
      if (!topDown) {
        LocalizationTask locTask(sys);
        locTask.run();
      }
      const auto& envCoeff = sys->getActiveOrbitalController<SCFMode>()->getCoefficients();
      auto orbitalPopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
          envCoeff, sys->getOneElectronIntegralController()->getOverlapIntegrals(),
          sys->getAtomCenteredBasisController()->getBasisIndices());
      auto distantOrbitals = SystemSplittingTools<SCFMode>::selectDistantOrbitals(
          orbitalPopulations, activeSystem, sys, settings.basisFunctionRatio, settings.borderAtomThreshold);
      auto nonOrthogonalDensityMatrix = SystemSplittingTools<SCFMode>::buildNonOrthogonalDensityMatrix(sys, distantOrbitals);
      _notProjectedEnvDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(*nonOrthogonalDensityMatrix));
    } // if settings->longRangeNaddKinFunc != CompositeFunctionals::KINFUNCTIONALS::NONE
  }   // for sys
  if (settings.longRangeNaddKinFunc != CompositeFunctionals::KINFUNCTIONALS::NONE) {
    assert(supersystemgrid);
    // Build non-additive kinetic energy potential for not projected systems
    _naddKinPot = std::make_shared<NAddFuncPotential<SCFMode>>(
        activeSystem, activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        _notProjectedEnvDensities, supersystemgrid,
        settings.customLongRangeNaddKinFunc.basicFunctionals.size() ? Functional(settings.customLongRangeNaddKinFunc)
                                                                    : resolveFunctional(settings.longRangeNaddKinFunc),
        std::make_pair(false, (gridCutOff < 0) ? allEConts : std::vector<std::shared_ptr<EnergyComponentController>>(0)),
        (gridCutOff < 0.0) ? true : false);
  }
  BasisFunctionMapper basisFunctionMapper(activeSystem->getBasisController());
  //  ABEmbeddedBundleFactory<SCFMode> abBundleFactory;
  for (unsigned int iEnv = 0; iEnv < environmentSystems.size(); ++iEnv) {
    if (_s_ABs[iEnv]) {
      auto basisControllerB = environmentSystems[iEnv]->getBasisController();
      if (activeFockMatrix) {
        _AtoBProjections.push_back(basisFunctionMapper.getSparseProjection(basisControllerB));
        auto differentialBasis = basisFunctionMapper.getDifferentialBasis(environmentSystems[iEnv]->getBasisController());
        if (differentialBasis) {
          BasisFunctionMapper toBMapper(differentialBasis);
          _difftoBProjections.push_back(toBMapper.getSparseProjection(basisControllerB));
        }
        basisControllerB = differentialBasis;
      }
      else {
        _AtoBProjections.push_back(nullptr);
        Eigen::MatrixXd I =
            Eigen::MatrixXd::Identity(basisControllerB->getNBasisFunctions(), basisControllerB->getNBasisFunctions());
        _difftoBProjections.push_back(std::make_shared<Eigen::SparseMatrix<double>>(I.sparseView()));
      }
      if (basisControllerB) {
        EmbeddingSettings abBundleSettings = settings;
        auto adjustedSettings_ptr = adjustSettings(abBundleSettings);
        _abEmbeddedBundles.push_back(ABEmbeddedBundleFactory<SCFMode>::produce(
            activeSystem, basisControllerB, environmentSystems[iEnv]->getGeometry(), environmentSystems,
            adjustedSettings_ptr, topDown));
      } // if basisControllerB
      else {
        _abEmbeddedBundles.push_back(nullptr);
        _difftoBProjections.push_back(nullptr);
      }
    }
    else {
      _abEmbeddedBundles.push_back(nullptr);
      _AtoBProjections.push_back(nullptr);
      _difftoBProjections.push_back(nullptr);
    }
  } // for iEnv
  if (iOOptions.printDebugInfos) {
    writeInterSubsystemOccOverlap();
  }
  timeTaken(2, "Huzinaga -- pre-calculations");
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& HuzinagaProjectionPotential<SCFMode>::getMatrix() {
  if (!_potential) {
    if (iOOptions.printDebugInfos) {
      writeInterSubsystemOccOverlap();
    }
    // reset the fock matrix
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot) {
      pot_spin.setZero();
    };
    // Loop over all environmental systems.
    for (unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
      if (_s_ABs[iEnv]) {
        auto environmentSystem = _environmentSystems[iEnv].lock();
        const SPMatrix<SCFMode> f_AB = buildOuterDiagonalFockMatrix(iEnv); // N(act) x N(env)
        // get the outer diagonal overlap matrix s_AB
        const Eigen::MatrixXd& s_AB = *_s_ABs[iEnv]; // N(act) x N(env)
        DensityMatrix<SCFMode> projectedMatrix = _envDensityCont[iEnv]->getDensityMatrix();
        if (_naddKinPot) {
          const auto& nonOrthogonalDensityMatrix = _notProjectedEnvDensities[iEnv]->getDensityMatrix();
          projectedMatrix -= nonOrthogonalDensityMatrix;
        }
        // A factor of one half for restricted to account for the factor of two in the density matrix.
        double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 0.5 : 1.0;
        for_spin(f_AB, pot, projectedMatrix) {
          auto matrix = (f_AB_spin - _fermiShift * s_AB) * projectedMatrix_spin * s_AB.transpose();
          pot_spin += -scfFactor * (matrix + matrix.transpose());
        };
      } // if _projectedEnvSystems[iEnv]
    }   // for iEnv
    if (_naddKinPot) {
      *_potential += _naddKinPot->getMatrix();
    }
    if (iOOptions.printDebugInfos) {
      writeInterSubsystemOccOverlap();
    }
  } // if !_potential
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
double HuzinagaProjectionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  double energy = 0.0;
  if (!_potential)
    getMatrix();
  if (_naddKinPot)
    energy += _naddKinPot->getEnergy(P);
  return energy;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd HuzinagaProjectionPotential<SCFMode>::getGeomGradients() {
  Eigen::MatrixXd gradientContr(1, 3);
  gradientContr.setZero();

  return gradientContr;
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode> HuzinagaProjectionPotential<SCFMode>::buildOuterDiagonalFockMatrix(unsigned int iEnv) {
  takeTime("Huzinaga -- Build core potential");
  auto activeSystem = _activeSystem.lock();
  const unsigned int nBasisAct = activeSystem->getBasisController()->getNBasisFunctions();
  auto environmentSystem = _environmentSystems[iEnv].lock();
  const unsigned int nBasisEnv = environmentSystem->getBasisController()->getNBasisFunctions();
  // initialize the Fock matrix
  SPMatrix<SCFMode> f_AB(Eigen::MatrixXd::Zero(nBasisAct, nBasisEnv));
  if (_activeFockMatrix) {
    SPMatrix<SCFMode> f_AA =
        _activeFockMatrix->getFockMatrix(activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix(),
                                         std::make_shared<EnergyComponentController>());
    if (_naddKinPot)
      f_AA += _naddKinPot->getMatrix();
    const auto& projectionAtoB = *_AtoBProjections[iEnv];
    for_spin(f_AB, f_AA) {
      f_AB_spin += f_AA_spin * projectionAtoB.transpose();
    };
  }
  if (_abEmbeddedBundles[iEnv]) {
    const auto& projectionB = *_difftoBProjections[iEnv];
    SPMatrix<SCFMode> f_ABRemains = _abEmbeddedBundles[iEnv]->getABMatrix();
    for_spin(f_AB, f_ABRemains) {
      f_AB_spin += f_ABRemains_spin * projectionB.transpose();
    };
  }
  timeTaken(2, "Huzinaga -- Build core potential");
  return f_AB;
}

template<Options::SCF_MODES SCFMode>
void HuzinagaProjectionPotential<SCFMode>::writeInterSubsystemOccOverlap() {
  auto activeSystem = _activeSystem.lock();
  auto nOccAct = activeSystem->template getNOccupiedOrbitals<SCFMode>();
  const auto& coeffAct = activeSystem->template getActiveOrbitalController<SCFMode>()->getCoefficients();
  for (unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
    if (_s_ABs[iEnv]) {
      auto environmentSystem = _environmentSystems[iEnv].lock();
      const auto& coeffEnv = environmentSystem->template getActiveOrbitalController<SCFMode>()->getCoefficients();
      auto nOccEnv = environmentSystem->template getNOccupiedOrbitals<SCFMode>();
      assert(_s_ABs[iEnv]);
      const auto& s_AB = *_s_ABs[iEnv];
      // Calculate the average overlap between the occupied orbitals of the subsystems.
      for_spin(coeffAct, coeffEnv, nOccAct, nOccEnv) {
        double totalInterSystemOccOverlap = 0.0;
        Eigen::MatrixXd totalOverlap =
            (coeffAct_spin.leftCols(nOccAct_spin).transpose() * s_AB * coeffEnv_spin.leftCols(nOccEnv_spin))
                .array()
                .abs()
                .matrix();
        assert(totalOverlap.rows() == nOccAct_spin && totalOverlap.cols() == nOccEnv_spin);
        for (unsigned int i = 0; i < nOccAct_spin; ++i) {
          for (unsigned int j = 0; j < nOccEnv_spin; ++j) {
            double overlap = totalOverlap(i, j);
            totalInterSystemOccOverlap += overlap;
          }
        }
        unsigned int loopEnd = (nOccAct_spin * nOccEnv_spin >= 10) ? 10 : nOccAct_spin * nOccEnv_spin;
        std::cout << loopEnd << " largest overlaps " << activeSystem->getSystemName() << " "
                  << environmentSystem->getSystemName() << std::endl;
        for (unsigned int i = 0; i < loopEnd; ++i) {
          int indexAct = 0;
          int indexEnv = 0;
          totalOverlap.maxCoeff(&indexAct, &indexEnv);
          std::cout << "Act " << indexAct << "<->" << indexEnv << "  " << totalOverlap(indexAct, indexEnv) << std::endl;
          totalOverlap(indexAct, indexEnv) = 0.0;
        }
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "Total intersystem overlap of occ. orbitals(" << activeSystem->getSystemName() << " and "
                  << environmentSystem->getSystemName() << "): " << totalInterSystemOccOverlap << std::endl;
      };
    }
  }
}
template<Options::SCF_MODES SCFMode>
std::shared_ptr<EmbeddingSettings> HuzinagaProjectionPotential<SCFMode>::adjustSettings(EmbeddingSettings& settings) {
  settings.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
  if (settings.truncateProjector) {
    settings.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  }
  return std::make_shared<EmbeddingSettings>(settings);
}

template class HuzinagaProjectionPotential<Options::SCF_MODES::RESTRICTED>;
template class HuzinagaProjectionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
