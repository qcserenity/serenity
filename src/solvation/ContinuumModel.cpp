/**
 * @file ContinuumModel.cpp
 *
 * @author Moritz Bensberg
 * @date Mai 18, 2020
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
#include "solvation/ContinuumModel.h"
/* Include Serenity Internal Headers */
#include "data/grid/CoulombPotentialOnGridCalculator.h"       //calculateElectronNuclei(), calculateElectronElectron()
#include "data/grid/ElectrostaticPotentialOnGridController.h" //getPotential().
#include "data/matrices/DensityMatrixController.h"            //getDensityMatrix()
#include "geometry/MolecularSurfaceController.h"              //Cavity information.
#include "io/FormattedOutputStream.h"                         //Check print level.
#include "io/FormattedOutputStream.h"                         //Filtered output streams.
#include "io/HDF5.h"
#include "math/Matrix.h"          //Needed for getCoordinates of geometry.
#include "misc/Timing.h"          //Timings.
#include "settings/PCMSettings.h" //PCMSettings
#include "solvation/Solvents.h"   //Tabulated solvent data.
/* Include Std and External Headers */
#include <Eigen/Dense> //VectorXd
#include <cmath>       //M_PI

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ContinuumModel<SCFMode>::ContinuumModel(const PCMSettings& pcmSettings,
                                        std::shared_ptr<MolecularSurfaceController> molecularSurface,
                                        std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>> activePotential,
                                        std::vector<std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>> environmentPotentials)
  : _settings(pcmSettings),
    _molecularSurface(molecularSurface),
    _activePotential(activePotential),
    _environmentPotentials(environmentPotentials) {
  _activePotential->getDensityMatrixController()->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  for (auto envPot : _environmentPotentials)
    envPot->getDensityMatrixController()->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  _eps = pcmSettings.eps;
  if (pcmSettings.solvent != Options::PCM_SOLVENTS::EXPLICIT)
    _eps = Solvents::getStaticPermittivity(pcmSettings.solvent);
}

template<Options::SCF_MODES SCFMode>
ContinuumModel<SCFMode>::~ContinuumModel() = default;

template<Options::SCF_MODES SCFMode>
double ContinuumModel<SCFMode>::getActivePCMEnergy() {
  getPCMCharges();
  Timings::takeTime("Implicit Solvation (PCM)       ");
  double energy = calculateEnergy(*_pcmCharges, _activePotential->getPotential());
  Timings::timeTaken("Implicit Solvation (PCM)       ");
  return energy;
}

template<Options::SCF_MODES SCFMode>
double ContinuumModel<SCFMode>::getTotalPCMEnergy() {
  getPCMCharges();
  Timings::takeTime("Implicit Solvation (PCM)       ");
  GridPotential<RESTRICTED> totalElectrostaticPotential = _activePotential->getPotential();
  for (auto envElecPot : _environmentPotentials)
    totalElectrostaticPotential += envElecPot->getPotential();
  double energy = calculateEnergy(*_pcmCharges, totalElectrostaticPotential);
  Timings::timeTaken("Implicit Solvation (PCM)       ");
  return energy;
}
template<Options::SCF_MODES SCFMode>
double ContinuumModel<SCFMode>::calculateEnergy(const GridPotential<RESTRICTED>& pcmCharges,
                                                const GridPotential<RESTRICTED>& potential) {
  double energy = 0.5 * (pcmCharges.array() * potential.array()).sum();
  return energy;
}

template<Options::SCF_MODES SCFMode>
const GridPotential<RESTRICTED>& ContinuumModel<SCFMode>::getPCMCharges() {
  Timings::takeTime("Implicit Solvation (PCM)       ");
  Timings::takeTime(" Tech. -    PCM Surface Charges");
  if (!_pcmCharges) {
    _pcmCharges = std::make_shared<GridPotential<RESTRICTED>>(_molecularSurface);
    if (_molecularSurface->isLoaded()) {
      Eigen::VectorXd loadedPCMCharges;
      std::string filename = _molecularSurface->getChargesPath() + "/PCMChargesData.h5";
      HDF5::H5File file(filename, H5F_ACC_RDONLY);
      HDF5::dataset_exists(file, "pcmcharges");
      HDF5::load(file, "pcmcharges", loadedPCMCharges);
      file.close();
      *_pcmCharges += loadedPCMCharges;
    }
    else {
      if (!_K)
        decomposeCavityMatrix();
      if (_environmentPotentials.size() > 0) {
        GridPotential<RESTRICTED> totalElectrostaticPotential = _activePotential->getPotential();
        for (auto envElecPot : _environmentPotentials)
          totalElectrostaticPotential += envElecPot->getPotential();
        Eigen::setNbThreads(0);
        *_pcmCharges += Eigen::VectorXd(-*_K * totalElectrostaticPotential);
      }
      else {
        Eigen::setNbThreads(0);
        *_pcmCharges += Eigen::VectorXd(-*_K * _activePotential->getPotential());
      }
      // Scaling for CPCM(COSMO)
      if (_settings.solverType == Options::PCM_SOLVER_TYPES::CPCM) {
        double scaling = this->getCPCMScaling();
        *_pcmCharges *= scaling;
      }
    }
    // optionally saving charges to file. This happens every scf cycle for now.
    //  TODO: only output this quantity once the scf is converged
    if (_settings.saveCharges == true) {
      std::string filename = _molecularSurface->getChargesPath() + "/PCMChargesData.h5";
      HDF5::H5File file(filename, H5F_ACC_TRUNC);
      HDF5::save(file, "pcmcharges", *_pcmCharges);
      file.close();
    }
  }
  Timings::timeTaken(" Tech. -    PCM Surface Charges");
  Timings::timeTaken("Implicit Solvation (PCM)       ");
  return *_pcmCharges;
}

template<Options::SCF_MODES SCFMode>
void ContinuumModel<SCFMode>::decomposeCavityMatrix() {
  Timings::takeTime(" Tech. -      PCM Cavity Matrix");
  if (_settings.solverType == Options::PCM_SOLVER_TYPES::CPCM) {
    _K = std::make_shared<Eigen::MatrixXd>(_molecularSurface->getMatrixSinv());
  }
  else if (_settings.solverType == Options::PCM_SOLVER_TYPES::IEFPCM) {
    const Eigen::MatrixXd& S = _molecularSurface->getMatrixS();
    const Eigen::MatrixXd& Ainv = _molecularSurface->getMatrixAinv();
    const Eigen::MatrixXd& D = _molecularSurface->getMatrixD();
    const unsigned int gridSize = S.rows();
    const double epsFrac = 2.0 * M_PI * (_eps + 1) / (_eps - 1);
    Eigen::LLT<Eigen::MatrixXd> llt = ((epsFrac * Ainv - D) * S).llt();
    _K = std::make_shared<Eigen::MatrixXd>(gridSize, gridSize);
    if (llt.info() != Eigen::Success) {
      OutputControl::vOut
          << "Cholesky decomposition failed in PCM calculation! Not positive definite! Trying explicit inversion."
          << std::endl;
      *_K = ((epsFrac * Ainv - D) * S).inverse() * (2.0 * M_PI * Ainv - D);
    }
    else {
      *_K = llt.solve(2.0 * M_PI * Ainv - D);
    }
  }
  else {
    throw SerenityError("The PCM Solver chosen is not implemented yet. Please use CPCM or IEFPCM.");
  }
  Timings::timeTaken(" Tech. -      PCM Cavity Matrix");
}

template<Options::SCF_MODES SCFMode>
const PCMSettings& ContinuumModel<SCFMode>::getPCMSettings() {
  return _settings;
}
template<Options::SCF_MODES SCFMode>
std::shared_ptr<MolecularSurfaceController> ContinuumModel<SCFMode>::getMolecularSurfaceController() {
  return _molecularSurface;
}
template<Options::SCF_MODES SCFMode>
double ContinuumModel<SCFMode>::getCPCMScaling() {
  return (_eps - 1) / (_eps + _settings.correction);
}

template class ContinuumModel<Options::SCF_MODES::RESTRICTED>;
template class ContinuumModel<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
