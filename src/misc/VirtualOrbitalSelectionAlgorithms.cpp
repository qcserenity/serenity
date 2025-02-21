/**
 * @file VirtualOrbitalSelectionAlgorithms.cpp
 *
 * @date May 28, 2020
 * @author Johannes Toelle
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
#include "misc/VirtualOrbitalSelectionAlgorithms.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutputStream.h"
#include "math/linearAlgebra/MatrixFunctions.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "potentials/bundles/FDEPotentialBundleFactory.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <iomanip>

namespace Serenity {
template<Options::SCF_MODES SCFMode>
VirtualOrbitalSelectionAlgorithms<SCFMode>::VirtualOrbitalSelectionAlgorithms(
    std::shared_ptr<SystemController> activeSystem, std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _activeSystem(activeSystem), _environmentSystems(environmentSystems) {
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSelectionAlgorithms<SCFMode>::excludeProjection(CoefficientMatrix<SCFMode>& coefs,
                                                                   SpinPolarizedData<SCFMode, unsigned int>& nOcc,
                                                                   SpinPolarizedData<SCFMode, std::vector<unsigned int>>& indices) {
  auto indices_old = indices;
  // calculate overlap
  auto& libint = Libint::getInstance();
  for (unsigned int iSysEnv = 0; iSysEnv < _environmentSystems.size(); iSysEnv++) {
    auto envSys = _environmentSystems[iSysEnv];
    // First Basiscontroller --> Col entries in Overlap
    // Second Basiscontroller --> Row entries in Overlap
    // Order therefore here important
    auto overlapAB =
        libint.compute1eInts(LIBINT_OPERATOR::overlap, _activeSystem->getBasisController(), envSys->getBasisController());
    SpinPolarizedData<SCFMode, Eigen::MatrixXd> overlapMO;
    auto nOccEnv = envSys->template getNOccupiedOrbitals<SCFMode>();
    auto coefEnv = envSys->template getActiveOrbitalController<SCFMode>()->getCoefficients();
    // Calculate overlap Sia between occ env and virt act
    for_spin(nOccEnv, nOcc, coefEnv, coefs, overlapMO) {
      auto nBasisEnv = envSys->getBasisController()->getNBasisFunctions();
      auto nBasisAct = _activeSystem->getBasisController()->getNBasisFunctions();
      overlapMO_spin = coefEnv_spin.block(0, 0, nBasisEnv, nOccEnv_spin).transpose() * overlapAB *
                       coefs_spin.block(0, nOcc_spin, nBasisAct, nBasisAct - nOcc_spin);
    };
    // Find orbitals of env system through overlap criteria
    for_spin(nOcc, indices, indices_old, overlapMO) {
      for (unsigned int col = 0; col < overlapMO_spin.cols(); col++) {
        double sum = 0.0;
        for (unsigned int row = 0; row < overlapMO_spin.rows(); row++) {
          sum += std::fabs(overlapMO_spin(row, col));
        }
        // remove Orbital from Index
        if (sum > 0.95) {
          indices_spin.erase(std::remove(indices_spin.begin(), indices_spin.end(), indices_old_spin[col + nOcc_spin]),
                             indices_spin.end());
        }
      }
    };
  }
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSelectionAlgorithms<SCFMode>::virtualCanonicalOrbitalSpaceSelection(
    CoefficientMatrix<SCFMode>& coefs, SpinPolarizedData<SCFMode, unsigned int>& nOcc,
    SpinPolarizedData<SCFMode, unsigned int>& nVirt, SpinPolarizedData<SCFMode, std::vector<unsigned int>>& indices,
    double& localThresh, double& envThres, bool onlyOne) {
  // Get atom indices for non-dummy atoms
  auto atoms = _activeSystem->getAtoms();
  // Calculate Overlap between Virtual MOs in basis of atoms of active subsystem and entire basis
  const auto nBasisFunc = _activeSystem->getBasisController()->getNBasisFunctions();
  const auto basisIndices = _activeSystem->getAtomCenteredBasisController()->getBasisIndices();
  const auto& oneIntController = _activeSystem->getOneElectronIntegralController();
  const auto& overlaps = oneIntController->getOverlapIntegrals();

  unsigned int counter = 0;
  std::vector<unsigned int> index;
  for (auto atom : atoms) {
    if (atom->isDummy() == false) {
      index.push_back(counter);
    }
    counter += 1;
  }
  unsigned int start = basisIndices[index[0]].first;
  unsigned int end = basisIndices[index[index.size() - 1]].second;
  for_spin(coefs, indices, nOcc, nVirt) {
    unsigned int one = 0;
    Eigen::MatrixXd overlapMOModifief = coefs_spin.block(start, nOcc_spin, end - start, nVirt_spin).transpose() *
                                        overlaps.block(start, 0, end - start, nBasisFunc) *
                                        coefs_spin.block(0, nOcc_spin, nBasisFunc, nVirt_spin);
    for (unsigned int col = 0; col < overlapMOModifief.cols(); col++) {
      if (localThresh != 0 && one == 0) {
        if (overlapMOModifief(col, col) > localThresh) {
          indices_spin.push_back(col + nOcc_spin);
          if (onlyOne)
            one++;
        }
      }
      if (envThres != 0 && one == 0) {
        if ((1 - overlapMOModifief(col, col)) >= envThres) {
          indices_spin.push_back(col + nOcc_spin);
          if (onlyOne)
            one++;
        }
      }
    }
  };
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSelectionAlgorithms<SCFMode>::virtualOrbitalSpaceLocalization(
    CoefficientMatrix<SCFMode>& coefs, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenValues,
    SpinPolarizedData<SCFMode, unsigned int>& nOcc, SpinPolarizedData<SCFMode, unsigned int>& nVirt,
    FockMatrix<SCFMode>& embeddedFock, bool local) {
  auto atoms = _activeSystem->getAtoms();
  std::vector<unsigned int> index;
  unsigned int counter = 0;
  for (auto atom : atoms) {
    if (atom->isDummy() == false) {
      index.push_back(counter);
    }
    counter += 1;
  }
  auto basisIndices = _activeSystem->getAtomCenteredBasisController()->getBasisIndices();
  unsigned int start = basisIndices[index[0]].first;
  unsigned int end = basisIndices[index[index.size() - 1]].second;
  // Set up Matrix for SVD and perform SVD
  Eigen::MatrixXd overlapMatrix = _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals();
  Eigen::MatrixXd overlapSqrt = mSqrt_Sym(overlapMatrix);
  // new Orbital Indices
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> indices;
  SpinPolarizedData<SCFMode, Eigen::VectorXd> newEigenvaluesTemp;

  for_spin(coefs, embeddedFock, nOcc, nVirt, indices, eigenValues, newEigenvaluesTemp) {
    Eigen::MatrixXd coefNEW = overlapSqrt * coefs_spin;
    Eigen::MatrixXd forSVd = coefNEW.block(start, nOcc_spin, end - start, nVirt_spin).transpose() *
                             coefNEW.block(start, nOcc_spin, end - start, nVirt_spin);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(forSVd, Eigen::ComputeThinU | Eigen::ComputeThinV);
    coefNEW = coefs_spin.block(0, nOcc_spin, overlapMatrix.rows(), nVirt_spin) * svd.matrixV();

    std::cout << std::fixed << std::setprecision(4);
    if (local)
      OutputControl::dOut << "Singular Values:\n" << svd.singularValues() << std::endl;
    if (!local)
      OutputControl::dOut << "Singular Values:\n" << (svd.singularValues().array() - 1).abs() << std::endl;

    for (unsigned int i = 0; i < svd.singularValues().size(); i++) {
      if (local && (svd.singularValues())[i] > 0.5) {
        indices_spin.push_back(i);
      }
      else if (!local && (svd.singularValues())[i] < 0.5) {
        indices_spin.push_back(i);
      }
    }
    // Reshape coefficients based on singular values
    Eigen::MatrixXd temp = coefNEW;
    coefNEW.conservativeResize(overlapMatrix.rows(), indices_spin.size());
    coefNEW.setZero();
    for (unsigned int i = 0; i < indices_spin.size(); i++) {
      coefNEW.col(i) = temp.col(indices_spin[i]);
    }
    // Transform Fock Matrix
    Eigen::MatrixXd f_transformed = coefNEW.transpose() * embeddedFock_spin * coefNEW;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> subspaceSolver(f_transformed);
    Eigen::MatrixXd eigenvectors = subspaceSolver.eigenvectors();
    Eigen::VectorXd eigenvalues = subspaceSolver.eigenvalues();
    newEigenvaluesTemp_spin = eigenvalues;
    // Transform new eigenvectors (coefficients) back in old basis
    coefs_spin.block(0, nOcc_spin, coefNEW.rows(), coefs_spin.cols() - nOcc_spin).setZero();
    coefs_spin.block(0, nOcc_spin, coefNEW.rows(), coefNEW.cols()) = coefNEW * eigenvectors;
    eigenValues_spin.conservativeResize(nOcc_spin + newEigenvaluesTemp_spin.size());
    eigenValues_spin.segment(nOcc_spin, newEigenvaluesTemp_spin.size()) = newEigenvaluesTemp_spin;
  };
}

template<Options::SCF_MODES SCFMode>
void VirtualOrbitalSelectionAlgorithms<SCFMode>::occupiedVirtualMixing(bool relaxation, EmbeddingSettings settings) {
  // Combined Supersystem
  CoefficientMatrix<SCFMode> coefs(_activeSystem->getBasisController());
  SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenValues;
  // Env system 1
  CoefficientMatrix<SCFMode> coefsEnv1 =
      _environmentSystems[0]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
  SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenValuesEnv1 =
      _environmentSystems[0]->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
  SpinPolarizedData<SCFMode, Eigen::VectorXi> coreOrbitalsEnv1 =
      _environmentSystems[0]->template getActiveOrbitalController<SCFMode>()->getOrbitalFlags();
  auto nBasisFuncEnv1 = _environmentSystems[0]->getBasisController()->getNBasisFunctions();
  auto nOccEnv1 = _environmentSystems[0]->template getNOccupiedOrbitals<SCFMode>();
  // Env system 2
  CoefficientMatrix<SCFMode> coefsEnv2 =
      _environmentSystems[1]->template getActiveOrbitalController<SCFMode>()->getCoefficients();
  SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenValuesEnv2 =
      _environmentSystems[1]->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
  SpinPolarizedData<SCFMode, Eigen::VectorXi> coreOrbitalsEnv2 =
      _environmentSystems[1]->template getActiveOrbitalController<SCFMode>()->getOrbitalFlags();
  auto nBasisFuncEnv2 = _environmentSystems[1]->getBasisController()->getNBasisFunctions();
  auto nOccEnv2 = _environmentSystems[1]->template getNOccupiedOrbitals<SCFMode>();
  // The charge of the new system needs to be adjusted with respect to subsystem 1 from which the occupied orbitals are
  // taken. Otherwise the wrong number of occ orbitals in that system is set
  const int charge = _environmentSystems[0]->getCharge();
  for_spin(coefs, coefsEnv1, coefsEnv2, nOccEnv1, nOccEnv2) {
    coefs_spin.conservativeResize(nBasisFuncEnv1 + nBasisFuncEnv2, nBasisFuncEnv1 + nBasisFuncEnv2);
    // ToDo assumes that the active and evironment systems do not have ghost atoms !
    coefs_spin.block(0, 0, nBasisFuncEnv1, nOccEnv1_spin) = coefsEnv1_spin.block(0, 0, nBasisFuncEnv1, nOccEnv1_spin);
    coefs_spin.block(nBasisFuncEnv1, nOccEnv1_spin, nBasisFuncEnv2, nBasisFuncEnv2 - nOccEnv2_spin) =
        coefsEnv2_spin.block(0, nOccEnv2_spin, nBasisFuncEnv2, nBasisFuncEnv2 - nOccEnv2_spin);
  };
  for_spin(eigenValues, eigenValuesEnv1, eigenValuesEnv2, nOccEnv1, nOccEnv2) {
    eigenValues_spin.conservativeResize(nBasisFuncEnv1 + nBasisFuncEnv2);
    eigenValues_spin.segment(0, nOccEnv1_spin) = eigenValuesEnv1_spin.segment(0, nOccEnv1_spin);
    eigenValues_spin.segment(nOccEnv1_spin, eigenValuesEnv2_spin.size() - nOccEnv2_spin) =
        eigenValuesEnv2_spin.segment(nOccEnv2_spin, eigenValuesEnv2_spin.size() - nOccEnv2_spin);
  };
  auto coefs_ptr = std::unique_ptr<CoefficientMatrix<SCFMode>>(new CoefficientMatrix<SCFMode>(coefs));
  auto eigenValue_ptr = std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(
      new SpinPolarizedData<SCFMode, Eigen::VectorXd>(eigenValues));
  auto coreOrbitals_ptr = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXi>>(
      _environmentSystems[0]->template getActiveOrbitalController<SCFMode>()->getOrbitalFlags());
  // Build new orbital controller
  _activeSystem->setDiskMode(true);
  auto orbitalSet = std::make_shared<OrbitalController<SCFMode>>(std::move(coefs_ptr), _activeSystem->getBasisController(),
                                                                 std::move(eigenValue_ptr), std::move(coreOrbitals_ptr));
  auto electronStructure =
      std::make_shared<ElectronicStructure<SCFMode>>(orbitalSet, _activeSystem->getOneElectronIntegralController(), nOccEnv1);
  // Add subsystem charge to get final charge for correct number of occupied orbitals
  _activeSystem->setCharge(charge);
  _activeSystem->setElectronicStructure<SCFMode>(electronStructure);
  // Save new electron structure
  _activeSystem->getElectronicStructure<SCFMode>()->toHDF5(_activeSystem->getHDF5BaseName(),
                                                           _activeSystem->getSettings().identifier);
  _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->updateDensityMatrix();
  // Perform Relaxation
  if (relaxation) {
    // Remove occupied contributions to the virtual orbitals via projection
    auto densAct =
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
    auto overlap = _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals();
    for_spin(densAct, nOccEnv1, nOccEnv2, coefs) {
      Eigen::MatrixXd coefsOld = coefs_spin;
      Eigen::MatrixXd diagOnes = Eigen::MatrixXd::Identity(overlap.rows(), overlap.cols());
      unsigned int virtEnv = nBasisFuncEnv2 - nOccEnv2_spin;
      // Set up 1 - Projector
      diagOnes -= (0.5 * densAct_spin * overlap);
      // Project out occupied contributions
      Eigen::MatrixXd temp = diagOnes * coefs_spin.block(0, nOccEnv1_spin, overlap.rows(), virtEnv);
      // Orthogonalize new coefficient against each other:
      Eigen::MatrixXd overlapSqrt = mSqrt_Sym(overlap);
      temp = overlapSqrt * temp;
      Eigen::MatrixXd overlapInverseSqrt = pseudoInversSqrt_Sym(overlap);
      std::vector<Eigen::MatrixXd> test;
      test.push_back(temp);
      Orthogonalization::modifiedGramSchmidtLinDep(test, 0.0);
      temp = test[0];
      temp = overlapInverseSqrt * temp;
      // Transform the virtual orbital coefficients
      coefs_spin.block(0, nOccEnv1_spin, overlap.rows(), virtEnv) = temp;
    };
    // Rediagonalize the fock matrix in supersystem basis
    // 1) For the virtual orbitals
    std::shared_ptr<SystemController> envSysCombined(nullptr);
    auto superSystemGeometry = std::make_shared<Geometry>(*_environmentSystems[1]->getGeometry());
    superSystemGeometry->addAsDummy(*_environmentSystems[0]->getGeometry(), true);
    // Settings and set new name
    auto envSysCombinedSettings = _environmentSystems[0]->getSettings();
    envSysCombinedSettings.name = _environmentSystems[1]->getSettings().name + "_SupersystemBasis";
    envSysCombined = std::make_shared<SystemController>(superSystemGeometry, envSysCombinedSettings);
    unsigned int chargeEnv = 0;
    CoefficientMatrix<SCFMode> coefsEnvNew(envSysCombined->getBasisController());
    for_spin(coefsEnvNew, coefsEnv1, coefsEnv2, nOccEnv1, nOccEnv2) {
      coefsEnvNew_spin.conservativeResize(nBasisFuncEnv1 + nBasisFuncEnv2, nBasisFuncEnv1 + nBasisFuncEnv2);
      coefsEnvNew_spin.block(nBasisFuncEnv1, 0, nBasisFuncEnv2, nOccEnv2_spin) =
          coefsEnv2_spin.block(0, 0, nBasisFuncEnv2, nOccEnv2_spin);
      coefsEnvNew_spin.block(nBasisFuncEnv2, nOccEnv2_spin, nBasisFuncEnv1, nBasisFuncEnv1 - nOccEnv1_spin) =
          coefsEnv1_spin.block(0, nOccEnv1_spin, nBasisFuncEnv1, nBasisFuncEnv1 - nOccEnv1_spin);
    };
    chargeEnv *= (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2 : 1;
    chargeEnv += _environmentSystems[1]->getCharge();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenValuesEnvNew;
    for_spin(eigenValuesEnvNew, eigenValuesEnv1, eigenValuesEnv2, nOccEnv1, nOccEnv2) {
      eigenValuesEnvNew_spin.conservativeResize(nBasisFuncEnv1 + nBasisFuncEnv2);
      eigenValuesEnvNew_spin.segment(0, nOccEnv2_spin) = eigenValuesEnv2_spin.segment(0, nOccEnv2_spin);
      eigenValuesEnvNew_spin.segment(nOccEnv2_spin, eigenValuesEnv1_spin.size() - nOccEnv1_spin) =
          eigenValuesEnv1_spin.segment(nOccEnv1_spin, eigenValuesEnv1_spin.size() - nOccEnv1_spin);
    };
    SpinPolarizedData<SCFMode, Eigen::VectorXi> coreOrbitalsEnvNew;
    for_spin(nOccEnv1, nOccEnv2, coreOrbitalsEnvNew, coreOrbitalsEnv1, coreOrbitalsEnv2) {
      coreOrbitalsEnvNew_spin.conservativeResize(nBasisFuncEnv1 + nBasisFuncEnv2);
      coreOrbitalsEnvNew_spin.segment(0, nOccEnv2_spin) = coreOrbitalsEnv2_spin.segment(0, nOccEnv2_spin);
      coreOrbitalsEnvNew_spin.segment(nOccEnv2_spin, coreOrbitalsEnv1_spin.size() - nOccEnv1_spin) =
          coreOrbitalsEnv1_spin.segment(nOccEnv1_spin, coreOrbitalsEnv1_spin.size() - nOccEnv1_spin);
    };

    auto coefs_ptrEnv = std::unique_ptr<CoefficientMatrix<SCFMode>>(new CoefficientMatrix<SCFMode>(coefsEnvNew));
    auto eigenValue_ptrEnv = std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(
        new SpinPolarizedData<SCFMode, Eigen::VectorXd>(eigenValuesEnvNew));
    auto coreOrbitals_ptrEnv = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXi>>();
    // Build new orbital controller
    auto orbitalSetEnv =
        std::make_shared<OrbitalController<SCFMode>>(std::move(coefs_ptrEnv), envSysCombined->getBasisController(),
                                                     std::move(eigenValue_ptrEnv), std::move(coreOrbitals_ptrEnv));
    auto electronStructureEnv = std::make_shared<ElectronicStructure<SCFMode>>(
        orbitalSetEnv, envSysCombined->getOneElectronIntegralController(), nOccEnv2);
    envSysCombined->setCharge(chargeEnv);
    envSysCombined->setElectronicStructure<SCFMode>(electronStructureEnv);
    envSysCombined->getElectronicStructure<SCFMode>()->getDensityMatrixController()->updateDensityMatrix();
    // calc embedded Fock Matrix
    std::vector<std::shared_ptr<SystemController>> actsys_vec;
    actsys_vec.push_back(_activeSystem);
    auto fock = this->calcEmbeddedFockMatrix(envSysCombined, actsys_vec, settings);
    for_spin(fock, coefs, nOccEnv1, nOccEnv2, eigenValues) {
      unsigned int virtEnv = nBasisFuncEnv2 - nOccEnv2_spin;
      fock_spin = coefs_spin.block(0, nOccEnv1_spin, overlap.rows(), virtEnv).transpose() * fock_spin *
                  coefs_spin.block(0, nOccEnv1_spin, overlap.rows(), virtEnv);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> subspaceSolver(fock_spin);
      Eigen::MatrixXd eigenvectors = subspaceSolver.eigenvectors();
      Eigen::VectorXd eigenvalues = subspaceSolver.eigenvalues();
      coefs_spin.block(0, nOccEnv1_spin, overlap.rows(), virtEnv) =
          coefs_spin.block(0, nOccEnv1_spin, overlap.rows(), virtEnv) * eigenvectors;
      eigenValues_spin.segment(nOccEnv1_spin, virtEnv) = eigenvalues;
    };

    // Rediagonalize the fock matrix in supersystem basis
    // 2) For the occupied orbitals
    std::vector<std::shared_ptr<SystemController>> envSystems;
    envSystems.push_back(_environmentSystems[1]);
    auto fock_new = this->calcEmbeddedFockMatrix(_activeSystem, envSystems, settings);
    for_spin(fock_new, coefs, nOccEnv1, eigenValues) {
      fock_new_spin = coefs_spin.block(0, 0, overlap.rows(), nOccEnv1_spin).transpose() * fock_new_spin *
                      coefs_spin.block(0, 0, overlap.rows(), nOccEnv1_spin);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> subspaceSolver(fock_new_spin);
      Eigen::MatrixXd eigenvectors = subspaceSolver.eigenvectors();
      Eigen::VectorXd eigenvalues = subspaceSolver.eigenvalues();
      coefs_spin.block(0, 0, overlap.rows(), nOccEnv1_spin) =
          coefs_spin.block(0, 0, overlap.rows(), nOccEnv1_spin) * eigenvectors;
      // Set new eigenvalues
      eigenValues_spin.segment(0, nOccEnv1_spin) = eigenvalues;
    };
  } // End of Relaxations
  // The eigenvalues need to be resized after the electron structure is set otherwise the dimensions internally are
  // not matching
  for_spin(eigenValues, nOccEnv1, nOccEnv2) {
    // Change The number of actual eigenvalues
    unsigned int virtEnv = nBasisFuncEnv2 - nOccEnv2_spin;
    eigenValues_spin.segment(nOccEnv1_spin + virtEnv, eigenValues_spin.size() - (nOccEnv1_spin + virtEnv)).array() =
        std::numeric_limits<double>::infinity();
  };
  _activeSystem->getActiveOrbitalController<SCFMode>()->updateOrbitals(coefs, eigenValues);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>
VirtualOrbitalSelectionAlgorithms<SCFMode>::calcEmbeddedFockMatrix(std::shared_ptr<SystemController>& activeSystem,
                                                                   std::vector<std::shared_ptr<SystemController>>& environmentSystem,
                                                                   EmbeddingSettings settings) {
  // list of environment density matrices (their controllers)
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensities;
  for (auto& sys : environmentSystem) {
    if (sys->getSettings().scfMode == SCFMode) {
      envDensities.push_back(sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController());
    }
    else {
      if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        // Build unrestricted DensityMatrixController
        DensityMatrix<SCFMode> uDensMat(sys->getBasisController());
        for_spin(uDensMat) {
          uDensMat_spin = 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
        };
        envDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(uDensMat));
      }
      else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        // Build restricted DensityMatrixController
        DensityMatrix<SCFMode> rDensMat(sys->getBasisController());
        for_spin(rDensMat) {
          rDensMat_spin =
              sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
        };
        envDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(rDensMat));
      }
      else {
        assert(false);
      }
    }
  }
  auto pbePot = FDEPotentialBundleFactory<SCFMode>::produce(
      activeSystem, activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(), environmentSystem,
      envDensities, std::make_shared<EmbeddingSettings>(settings), activeSystem->getGridController(), nullptr, false);
  return pbePot->getFockMatrix(activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix(),
                               std::make_shared<EnergyComponentController>());
}

template class VirtualOrbitalSelectionAlgorithms<Options::SCF_MODES::RESTRICTED>;
template class VirtualOrbitalSelectionAlgorithms<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
