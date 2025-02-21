/**
 * @file ALMOPotentials.cpp
 *
 * @author Lukas Lampe
 * @date Nov 6, 2024
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
#include "potentials/bundles/ALMOPotentials.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "basis/BasisFunctionMapper.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h"
#include "misc/SystemSplittingTools.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/SystemAdditionTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ALMOPotentials<SCFMode>::ALMOPotentials(std::shared_ptr<SystemController> activeSystem,
                                        std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _activeSystem(activeSystem) {
  std::vector<std::shared_ptr<SystemController>> allSystems = {_activeSystem};
  Settings settings = _activeSystem->getSettings();
  settings.spin = 0;
  settings.charge = 0;
  settings.name += "_tmp";
  for (auto e : environmentSystems) {
    _environmentSystems.push_back(e);
    allSystems.push_back(e);
  }

  _supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settings);
  SystemAdditionTask<SCFMode> additionTask(_supersystem, allSystems);
  additionTask.run();

  OutputControl::nOut << std::endl << "Construct ALMO potential from supersystem..." << std::endl << std::endl;

  if (settings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT)
    _superPotentialBundle = _supersystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  else
    _superPotentialBundle = _supersystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();

  auto es = _supersystem->getElectronicStructure<SCFMode>();
  es->getDensityMatrixController()->setALMO(es->getOneElectronIntegralController());

  _s_ABs = SystemSplittingTools<SCFMode>::getProjectedSubsystems(_activeSystem, _environmentSystems);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> ALMOPotentials<SCFMode>::getFockMatrix(const DensityMatrix<SCFMode>& P,
                                                           std::shared_ptr<EnergyComponentController> energies) {
  // The density matrix P is not used since the Fock matrix is constructed from the supersystem
  // Update supersystem coefficients and density
  auto es = _supersystem->getElectronicStructure<SCFMode>();
  auto basisController = _supersystem->getBasisController();
  unsigned int nBasis = basisController->getNBasisFunctions();
  CoefficientMatrix<SCFMode> coefficients(basisController);
  SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenvalues = Eigen::VectorXd::Zero(nBasis);
  SpinPolarizedData<SCFMode, Eigen::VectorXi> coreOrbitals = Eigen::VectorXi::Zero(nBasis);
  SpinPolarizedData<SCFMode, unsigned int> nOccSuper(0);
  BasisFunctionMapper mapper(basisController);
  std::vector<std::shared_ptr<SystemController>> allSystems = {_activeSystem};
  for (auto e : _environmentSystems)
    allSystems.push_back(e);
  for (unsigned int iSub = 0; iSub < allSystems.size(); ++iSub) {
    CoefficientMatrix<SCFMode> subCoefficients = allSystems[iSub]->getActiveOrbitalController<SCFMode>()->getCoefficients();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> subEigenvalues =
        allSystems[iSub]->getActiveOrbitalController<SCFMode>()->getEigenvalues();
    SpinPolarizedData<SCFMode, Eigen::VectorXi> subCoreOrbitals =
        allSystems[iSub]->getActiveOrbitalController<SCFMode>()->getOrbitalFlags();
    auto nOccSub = allSystems[iSub]->getNOccupiedOrbitals<SCFMode>();
    Eigen::MatrixXd projection = *mapper.getSparseProjection(allSystems[iSub]->getBasisController());
    for_spin(subCoefficients, coefficients, subEigenvalues, eigenvalues, subCoreOrbitals, coreOrbitals, nOccSuper, nOccSub) {
      coefficients_spin.block(0, nOccSuper_spin, nBasis, nOccSub_spin) =
          projection.transpose() * subCoefficients_spin.leftCols(nOccSub_spin);
      eigenvalues_spin.segment(nOccSuper_spin, nOccSub_spin) = subEigenvalues_spin.head(nOccSub_spin);
      coreOrbitals_spin.segment(nOccSuper_spin, nOccSub_spin) = subCoreOrbitals_spin.head(nOccSub_spin);
      nOccSuper_spin += nOccSub_spin;
    };
  }
  _supersystem->getActiveOrbitalController<SCFMode>()->updateOrbitals(coefficients, eigenvalues, coreOrbitals);
  auto superP = es->getDensityMatrix();
  // Construct activeSystem Fock matrix from supersystem Fock matrix
  double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 0.5 : 1.0;
  auto F = _superPotentialBundle->getFockMatrix(superP, energies);
  energies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES, 0.0);
  unsigned int nActBasis = _activeSystem->getBasisController()->getNBasisFunctions();
  _potential.reset(new FockMatrix<SCFMode>(_activeSystem->getBasisController()));
  FockMatrix<SCFMode>& pot = *_potential;
  for_spin(pot, F) {
    pot_spin = F_spin.block(0, 0, nActBasis, nActBasis);
  };
  unsigned int iBlockStart = nActBasis;
  for (unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
    unsigned int nEnvBasisI = _environmentSystems[iEnv]->getBasisController()->getNBasisFunctions();
    auto D_BB = _environmentSystems[iEnv]->getElectronicStructure<SCFMode>()->getDensityMatrix();
    Eigen::MatrixXd S_AB = *_s_ABs[iEnv];
    Eigen::MatrixXd S_BA = (*_s_ABs[iEnv]).transpose();
    for_spin(pot, F, D_BB) {
      Eigen::MatrixXd F_AB = F_spin.block(0, iBlockStart, nActBasis, nEnvBasisI);
      Eigen::MatrixXd F_BA = F_spin.block(iBlockStart, 0, nEnvBasisI, nActBasis);
      pot_spin -= scfFactor * S_AB * D_BB_spin * F_BA;
      pot_spin -= scfFactor * F_AB * D_BB_spin * S_BA;
    };
    unsigned int jBlockStart = nActBasis;
    for (unsigned int jEnv = 0; jEnv < _environmentSystems.size(); ++jEnv) {
      unsigned int nEnvBasisJ = _environmentSystems[jEnv]->getBasisController()->getNBasisFunctions();
      auto D_CC = _environmentSystems[jEnv]->getElectronicStructure<SCFMode>()->getDensityMatrix();
      Eigen::MatrixXd S_CA = (*_s_ABs[jEnv]).transpose();
      for_spin(pot, F, D_BB, D_CC) {
        Eigen::MatrixXd F_BC = F_spin.block(iBlockStart, jBlockStart, nEnvBasisI, nEnvBasisJ);
        pot_spin += scfFactor * scfFactor * S_AB * D_BB_spin * F_BC * D_CC_spin * S_CA;
      };
      jBlockStart += nEnvBasisJ;
    }
    iBlockStart += nEnvBasisI;
  }
  return *_potential;
}
#pragma GCC diagnostic pop

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd ALMOPotentials<SCFMode>::getGradients() {
  throw SerenityError("No geometry gradients available for the ALMO potentials.");
  Eigen::MatrixXd gradientContr(1, 3);
  return gradientContr;
}

template class ALMOPotentials<Options::SCF_MODES::RESTRICTED>;
template class ALMOPotentials<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
