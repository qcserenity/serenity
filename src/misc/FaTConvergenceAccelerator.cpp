/**
 * @file FaTConvergenceAccelerator.cpp
 *
 * @date Mar 29, 2018
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
#include "misc/FaTConvergenceAccelerator.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/diis/DIIS.h"
#include "memory/MemoryManager.h"
#include "misc/SystemSplittingTools.h"
#include "misc/Timing.h"
#include "potentials/BUReconstructionPotential.h"
#include "potentials/bundles/FDEPotentialBundleFactory.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FaTConvergenceAccelerator<SCFMode>::FaTConvergenceAccelerator(unsigned int maxStore, const FreezeAndThawTaskSettings& settings,
                                                              std::vector<std::shared_ptr<SystemController>> activeSystems,
                                                              std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _activeSystems(activeSystems), _environmentSystems(environmentSystems), _settings(settings) {
  // initialize the error vector, diis and density vector
  // calculate cache limit for the vector controller with a buffer of 2 GB
  auto memoManger = MemoryManager::getInstance();
  long long availableMemory =
      (memoManger->getAvailableSystemMemory() > 2e+9) ? memoManger->getAvailableSystemMemory() - 2e+9 : 0;
  double cacheLimit = (double)availableMemory / ((double)maxStore + 3.0);
  _densityVector = std::make_shared<VectorOnDiskStorageController>(cacheLimit, "Density.h5");
  _fpsMinusSPF = std::make_shared<VectorOnDiskStorageController>(cacheLimit, "Error.h5");
  _diis = std::make_shared<DIIS>(maxStore, true);
}
template<Options::SCF_MODES SCFMode>
void FaTConvergenceAccelerator<SCFMode>::accelerateConvergence() {
  takeTime("Freeze and Thaw DIIS");
  _cycle++;
  // map the density matrices to a vector
  for (auto& sys : _activeSystems) {
    unsigned int nBasisFunctions = sys->getBasisController()->getNBasisFunctions();
    DensityMatrix<SCFMode> p =
        sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
    unsigned int vectorSizeFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 1 : 2;
    auto newVectorSegment = std::make_shared<Eigen::VectorXd>(vectorSizeFactor * nBasisFunctions * nBasisFunctions);
    unsigned int lastSegmentIndex = 0;
    for_spin(p) {
      newVectorSegment->block(lastSegmentIndex, 0, nBasisFunctions * nBasisFunctions, 1) =
          Eigen::Map<Eigen::VectorXd>(p_spin.data(), nBasisFunctions * nBasisFunctions);
      lastSegmentIndex += nBasisFunctions * nBasisFunctions;
    };
    _densityVector->storeVectorSegment(newVectorSegment, sys->getSystemName());
  }
  // Calculate RMSD of the density
  double rmsd = calcRMSDofDensity();

  // Start DIIS if the RMSD is below a given threshold.
  // Else damp.
  std::cout << "-------------------------------------" << std::endl;
  std::cout << "Density RMSD (all systems): " << rmsd << std::endl;
  if (rmsd < _settings.diisStart && rmsd > _settings.diisEnd) {
    // Calculate error measure.
    calcFPSminusSPF();
    std::cout << "+++ performing DIIS step +++" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    // perform DIIS step
    _diis->optimize(*_densityVector, *_fpsMinusSPF);

    // Set the optimized density in the system controller
    for (auto& sys : _activeSystems) {
      DensityMatrix<SCFMode> optDensity(sys->getBasisController());
      unsigned int nBasisFunctions = sys->getBasisController()->getNBasisFunctions();
      auto optVectorSegmentD = _densityVector->getVectorSegment(sys->getSystemName());
      unsigned int lastSegmentIndex = 0;
      for_spin(optDensity) {
        optDensity_spin = Eigen::Map<Eigen::MatrixXd>(
            optVectorSegmentD->block(lastSegmentIndex, 0, nBasisFunctions * nBasisFunctions, 1).data(), nBasisFunctions,
            nBasisFunctions);
        lastSegmentIndex += nBasisFunctions * nBasisFunctions;
      };
      sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->setDensityMatrix(optDensity);

    } /* for sys */
  }   /* rmsd < _diisStartError */
  else {
    std::cout << "-------------------------------------" << std::endl;
  }
  timeTaken(3, "Freeze and Thaw DIIS");
}

template<Options::SCF_MODES SCFMode>
void FaTConvergenceAccelerator<SCFMode>::calcFPSminusSPF() {
  /*
   * 1. Loop over active systems.
   * 2. Calculate embedded fock matrix.
   * 3. Calculate FPS-SPF
   * 4. Map the matrix to an vector segment and append the block
   */
  for (unsigned int i = 0; i < _activeSystems.size(); ++i) {
    auto f = calcEmbeddedFockMatrix(i);
    auto p = _activeSystems[i]->template getElectronicStructure<SCFMode>()->getDensityMatrix();
    auto s = _activeSystems[i]->getOneElectronIntegralController()->getOverlapIntegrals();
    unsigned int nBasisFunctions = _activeSystems[i]->getBasisController()->getNBasisFunctions();

    unsigned int vectorSizeFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 1 : 2;
    auto newVectorSegmentFPS_SPF = std::make_shared<Eigen::VectorXd>(vectorSizeFactor * nBasisFunctions * nBasisFunctions);
    auto newVectorSegmentF = std::make_shared<Eigen::VectorXd>(vectorSizeFactor * nBasisFunctions * nBasisFunctions);
    unsigned int lastSegmentIndex = 0;
    for_spin(f, p) {
      FockMatrix<Options::SCF_MODES::RESTRICTED> fps_spf(_activeSystems[i]->getBasisController());
      fps_spf = f_spin * p_spin * s - s * p_spin * f_spin;
      newVectorSegmentFPS_SPF->block(lastSegmentIndex, 0, s.rows() * s.cols(), 1) =
          Eigen::Map<Eigen::VectorXd>(fps_spf.data(), s.rows() * s.cols());
      // Save the fock matrix
      lastSegmentIndex += s.rows() * s.cols();
    };
    _fpsMinusSPF->storeVectorSegment(newVectorSegmentFPS_SPF, _activeSystems[i]->getSystemName());
  }
}
template<Options::SCF_MODES SCFMode>
double FaTConvergenceAccelerator<SCFMode>::calcRMSDofDensity() {
  double rmsd = 0.0;
  // If an old density vector is available calculate the RMSD.
  // Otherwise return inf.
  if (_oldDensityVector) {
    double sumOfSquares = 0.0;
    // Loop over systems
    for (const auto& sys : _activeSystems) {
      std::string label = sys->getSystemName();
      auto differenceSegment = *_densityVector->getVectorSegment(label) - *_oldDensityVector->getVectorSegment(label);
      double tmp = differenceSegment.cwiseProduct(differenceSegment).sum();
      sumOfSquares += tmp;
    }
    rmsd = std::sqrt(sumOfSquares / _densityVector->size());
    // Set new "old" density vector.
    _oldDensityVector = nullptr;
    _oldDensityVector.reset(new VectorOnDiskStorageController(*_densityVector, "OldDensity.h5"));
  }
  else {
    rmsd = std::numeric_limits<double>::infinity();
    _oldDensityVector = std::make_shared<VectorOnDiskStorageController>(*_densityVector, "OldDensity.h5");
  }
  return rmsd;
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> FaTConvergenceAccelerator<SCFMode>::calcEmbeddedFockMatrix(unsigned int activeSystemIndex) {
  // list of environment systems for this active system.
  std::vector<std::shared_ptr<SystemController>> envSystems;
  for (unsigned int i = 0; i < _activeSystems.size(); ++i) {
    if (i != activeSystemIndex)
      envSystems.push_back(_activeSystems[i]);
  }
  for (unsigned int i = 0; i < _environmentSystems.size(); ++i) {
    envSystems.push_back(_environmentSystems[i]);
  }
  // list of environment density matrices (their controllers)
  auto envDensities = SystemSplittingTools<SCFMode>::getEnvironmentDensityControllers(envSystems);
  auto activeSystem = _activeSystems[activeSystemIndex];
  auto grid = activeSystem->getGridController();

  auto fdePot = FDEPotentialBundleFactory<SCFMode>::produce(
      activeSystem, activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(), envSystems,
      envDensities, std::make_shared<EmbeddingSettings>(_settings.embedding), grid);

  return fdePot->getFockMatrix(activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix(),
                               std::make_shared<EnergyComponentController>());
}

template class FaTConvergenceAccelerator<Options::SCF_MODES::RESTRICTED>;
template class FaTConvergenceAccelerator<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
