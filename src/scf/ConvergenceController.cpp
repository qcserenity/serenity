/**
 * @file   ConvergenceController.cpp
 *
 * @date   28. Dezember 2013, 18:13
 * @author Thomas Dresselhaus, M. Boeckers
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
#include "scf/ConvergenceController.h"
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "energies/EnergyComponentController.h"
#include "energies/EnergyContributions.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutput.h"
#include "io/IOOptions.h"
#include "math/diis/ADIIS.h"
#include "math/diis/DIIS.h"
#include "math/linearAlgebra/MatrixFunctions.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "misc/Timing.h"
#include "scf/damper/Damper.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <cmath>
#include <iomanip>
#include <limits>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ConvergenceController<SCFMode>::ConvergenceController(const Settings& settings,
                                                      std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrix,
                                                      std::shared_ptr<OrbitalController<SCFMode>> orbitalController,
                                                      std::shared_ptr<OneElectronIntegralController> oneIntController,
                                                      const std::shared_ptr<EnergyComponentController> energyComponentController)
  : _settings(settings),
    _dmatContr(densityMatrix),
    _orbitalController(orbitalController),
    _oldP(nullptr),
    _oneIntController(oneIntController),
    _energyComponentController(energyComponentController),
    _oldEnergy(std::numeric_limits<double>::infinity()),
    _oldOneElEnergy(std::numeric_limits<double>::infinity()),
    _diisConvMeasure(std::numeric_limits<double>::infinity()),
    _rmsdOfDensity(std::numeric_limits<double>::infinity()),
    _orthoS(nullptr),
    _damping(nullptr),
    _diis(nullptr),
    _adiis(nullptr),
    _diisZoneStart(0),
    _mode("---") {
  _diis = std::make_shared<DIIS>(_settings.scf.diisMaxStore);
  if (settings.scf.useADIIS)
    _adiis = std::make_shared<ADIIS>();
  switch (_settings.scf.damping) {
    case (Options::DAMPING_ALGORITHMS::NONE): {
      _damping = nullptr;
      break;
    }
    case (Options::DAMPING_ALGORITHMS::STATIC): {
      _damping = std::make_shared<Damper<SCFMode>>(_settings.scf.staticDampingFactor);
      break;
    }
    case (Options::DAMPING_ALGORITHMS::SERIES): {
      _damping = std::make_shared<Damper<SCFMode>>(_settings.scf.seriesDampingStart, _settings.scf.seriesDampingStep,
                                                   _settings.scf.seriesDampingEnd, _settings.scf.seriesDampingInitialSteps);
      break;
    }
    case (Options::DAMPING_ALGORITHMS::DYNAMIC): {
      _damping = std::make_shared<Damper<SCFMode>>();
      break;
    }
  }
}

enum class CONVERGENCE_CRITERIA { E_TOT, RMSD_DENSITY, FPSminusSPF };

template<Options::SCF_MODES SCFMode>
bool ConvergenceController<SCFMode>::checkConvergence() {
  bool converged = false;

  ++_cycle;
  if (iOOptions.printSCFCycleInfo) {
    this->printCycleInfo();
  }
  clock_gettime(CLOCK_REALTIME, &_time);

  // check for convergence. Do not allow convergence if a level-shift is stll used.
  if (this->getNConverged() >= _nNeccessaryToConverge && not _levelShiftInLastIteration) {
    converged = true;
  }
  auto currentLevelshift = this->getLevelshift().first;
  _levelShiftInLastIteration = currentLevelshift[0] > 1e-9 || std::abs(currentLevelshift[1] - 1) > 1e-9;
  // Update old energy.
  _oldEnergy = _energyComponentController->getTotalEnergy();
  if (_oldEnergy > 1e+2) {
    throw SerenityError("Positive energy encountered during SCF procedure. The SCF is unlikely to converge!");
  }

  return converged;
}

template<Options::SCF_MODES SCFMode>
void ConvergenceController<SCFMode>::printCycleInfo() {
  const double energy = _energyComponentController->getTotalEnergy();
  double deltaE_abs = std::abs(energy - _oldEnergy);
  if (_cycle == 1) {
    const double inf = std::numeric_limits<double>::infinity();
    OutputControl::m.printf(
        "    Cycle %4s E/a.u. %7s abs(dE)/a.u. %3s rmsd(P)/a.u. %5s [F,P]/a.u. %4s time/min   Mode\n", "", "", "", "", "");
    OutputControl::m.printf("    %4d %16.10f %16.2e %16.2e %16.2e %6i:%02u \n", _cycle, energy, inf, inf,
                            _diisConvMeasure, 0, 0);
    OutputControl::mOut.flush();
  }
  else {
    timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    double sec = (double)(now.tv_sec - _time.tv_sec) + (now.tv_nsec - _time.tv_nsec) * 0.000000001;
    int ms = (int)(sec * 1000) - 1000 * (int)sec;
    std::printf("    %4d %16.10f %16.2e %16.2e %16.2e %6i:%02i:%03i    %3s\n", _cycle, energy, deltaE_abs,
                _rmsdOfDensity, _diisConvMeasure, (int)(sec / 60), (int)(sec) % 60, ms, _mode.c_str());
    OutputControl::mOut.flush();
  }
}

template<Options::SCF_MODES SCFMode>
unsigned int ConvergenceController<SCFMode>::getNConverged() {
  // Prepare convergence map (map is initialized to false by default)
  std::map<CONVERGENCE_CRITERIA, bool> convergenceMap;
  const double energy = _energyComponentController->getTotalEnergy();
  double deltaE_abs = std::abs(energy - _oldEnergy);

  if (deltaE_abs <= _settings.scf.energyThreshold) {
    convergenceMap[CONVERGENCE_CRITERIA::E_TOT] = true;
  }
  else {
    convergenceMap[CONVERGENCE_CRITERIA::E_TOT] = false;
  }
  if (_rmsdOfDensity <= _settings.scf.rmsdThreshold) {
    convergenceMap[CONVERGENCE_CRITERIA::RMSD_DENSITY] = true;
  }
  else {
    convergenceMap[CONVERGENCE_CRITERIA::RMSD_DENSITY] = false;
  }
  if (_diisConvMeasure <= _settings.scf.diisThreshold) {
    convergenceMap[CONVERGENCE_CRITERIA::FPSminusSPF] = true;
  }
  else {
    convergenceMap[CONVERGENCE_CRITERIA::FPSminusSPF] = false;
  }
  unsigned int nConverged = 0;
  for (auto& keyValuePair : convergenceMap) {
    nConverged += keyValuePair.second;
  }
  return nConverged;
}

template<Options::SCF_MODES SCFMode>
double ConvergenceController<SCFMode>::calcRMSDofDensity() {
  double rmsd = 0.0;
  if (_oldP) {
    auto P = _dmatContr->getDensityMatrix();
    DensityMatrix<SCFMode> difference = P - (*_oldP);
    for_spin(difference) {
      double tmp = (difference_spin.array() * difference_spin.array()).sum();
      rmsd = std::max(rmsd, std::sqrt(tmp / (difference_spin.cols() * difference_spin.cols())));
    };
    (*_oldP) = P;
  }
  else {
    rmsd = std::numeric_limits<double>::infinity();
    _oldP = std::make_shared<DensityMatrix<SCFMode>>(_dmatContr->getDensityMatrix());
  }
  return rmsd;
}

/*
 * The error vector to be used in DIIS as suggested by Pulay himself
 * in his second DIIS paper is FPS-SPF in the orthonormal basis.
 * Experience shows that this actually does work pretty well.
 */
template<Options::SCF_MODES SCFMode>
MatrixInBasis<SCFMode> ConvergenceController<SCFMode>::calcFPSminusSPF(FockMatrix<SCFMode>& F) {
  takeTime("Optimizer 1");
  MatrixInBasis<SCFMode> result(F.getBasisController());
  const auto& S = _oneIntController->getOverlapIntegrals();
  if (_orbitalController->getCustomOverlap() && _orbitalController->getCustomTransformMatrix()) {
    // Make sure to not construct density from the inverse of the MO overlap
    _dmatContr->setALMO(nullptr);
    _dmatContr->updateDensityMatrix();
    auto P = _dmatContr->getDensityMatrix();
    _dmatContr->setALMO(_oneIntController);
    auto customS = (*_orbitalController->getCustomOverlap());
    for_spin(result, F, P, customS) {
      result_spin = (F_spin * P_spin * customS_spin);
      result_spin -= (customS_spin * P_spin * F_spin);
    };
  }
  else {
    auto P = _dmatContr->getDensityMatrix();
    for_spin(result, F, P) {
      result_spin = (F_spin * P_spin * S);
      result_spin -= (S * P_spin * F_spin);
    };
  }
  timeTaken(3, "Optimizer 1");
  /*
   * Transform the new error vector to orthogonal basis.
   */
  takeTime("Optimizer 2");
  if (_orbitalController->getCustomOverlap() && _orbitalController->getCustomTransformMatrix()) {
    auto X = (*_orbitalController->getCustomTransformMatrix());
    for_spin(result, X) {
      result_spin = X_spin.transpose() * result_spin * X_spin;
    };
  }
  else {
    if (!_orthoS)
      _orthoS = std::make_unique<Eigen::MatrixXd>(S.householderQr().householderQ());
    for_spin(result) {
      result_spin = (*_orthoS) * result_spin * (*_orthoS).transpose();
    };
  }
  timeTaken(3, "Optimizer 2");
  return result;
}

template<Options::SCF_MODES SCFMode>
std::pair<Eigen::VectorXd, SpinPolarizedData<SCFMode, Eigen::VectorXd>> ConvergenceController<SCFMode>::getLevelshift() {
  /*
   * Generate levelshift data to modify FockMatrix:
   * F_{aa}+=levelshift[0]; levelshift[0] >= 0
   * F_{ia}*=levelshift[1]; 0 <= levelshift[1] <= 1
   * F_{ai}*=levelshift[1]
   *
   * Therefore, levelshift[0] should slowly transition to 0, while
   * levelshift[1] should develop towards 1. This is tried in the following.
   * If you find slow convergence or oscillations during the SCF, try to vary
   * the following damping functions...
   */
  Eigen::VectorXd levelshift(2);
  levelshift[0] = 0.0;
  levelshift[1] = 1.0;
  // Turn off any level-shift if the SCF would be converged otherwise.
  if (this->getNConverged() < _nNeccessaryToConverge) {
    if ((_diisConvMeasure > _settings.scf.diisThreshold * 100.0 and _settings.scf.useLevelshift))
      levelshift[0] = std::max(sqrt(log(_diisConvMeasure + 1)), _settings.scf.minimumLevelshift);
    if (_diisConvMeasure > std::min(_settings.scf.diisThreshold * 100.0, 1e-5) and _settings.scf.useOffDiagLevelshift)
      levelshift[1] = (1.0 / (_diisConvMeasure + 2)) + 0.5;
  }
  return std::pair<Eigen::VectorXd, SpinPolarizedData<SCFMode, Eigen::VectorXd>>(levelshift, _dmatContr->getOccupations());
}

template<>
void ConvergenceController<Options::SCF_MODES::RESTRICTED>::accelerateConvergence(FockMatrix<RESTRICTED>& F,
                                                                                  DensityMatrix<RESTRICTED> D) {
  Timings::takeTime("Tech. -           DIIS/Damping");

  if (_first)
    _orthoS.reset(new Eigen::MatrixXd((*_orbitalController->getTransformMatrix(_oneIntController)).transpose()));

  /*
   * Optimizer part: set up the new error vector
   */
  takeTime("Optimizer error measure");
  symInPlace<RESTRICTED>(F);
  Eigen::MatrixXd errorVector = calcFPSminusSPF(F);
  timeTaken(3, "Optimizer error measure");
  /*
   * Get the convergence measure.
   */
  _diisConvMeasure = errorVector.maxCoeff();
  // calculate root mean square change of density matrix
  _rmsdOfDensity = calcRMSDofDensity();
  /*
   * Switch on DIIS if electronic gradient falls below threshold, else damp
   */
  if (_diisConvMeasure <= _settings.scf.diisStartError && _diis->getNVectorsStored() > 0 && _cycle > 0)
    _diisZoneStart = _cycle;
  _mode = "---";
  if (this->getLevelshift().first[0] > 0.0)
    _mode[2] = 'L';
  if (_diisZoneStart) {
    /*
     * Make an Optimizer step: the Fock matrix is updated.
     */
    takeTime("Optimizer update");
    _adiis = nullptr;
    /* Optimize */
    Eigen::MatrixXd tmp(F);
    _diis->optimize(tmp, errorVector);
    F = tmp;
    _mode[1] = 'D';
    timeTaken(3, "Optimizer update");
  }
  else if (_adiis and !_first) {
    F = _adiis->optimize(F, _dmatContr->getDensityMatrix());
    _mode[1] = 'A';
  }
  // rework: only one instance of _damping, but several functions
  if (_diisConvMeasure >= _settings.scf.endDampErr or _cycle < 2) {
    if (_diisConvMeasure <= 10.0 * _settings.scf.diisStartError && not _diisZoneStart)
      _diis->storeMatrix(F, errorVector);
    if (_settings.scf.damping == Options::DAMPING_ALGORITHMS::SERIES) {
      _damping->arithmeticSeriesDamp(F);
    }
    else if (_settings.scf.damping == Options::DAMPING_ALGORITHMS::DYNAMIC) {
      _damping->dynamicDamp(F, D);
    }
    else if (_settings.scf.damping == Options::DAMPING_ALGORITHMS::STATIC) {
      _damping->staticDamp(F);
    }
    _mode[0] = 'D';
  }
  _first = false;
  if (_cycle % _settings.scf.diisFlush == 0 && _cycle > 0) {
    _diis->reinit();
  }

  /*
   * Done. _F (and thus the SCF procedure) now contains a Fock matrix that
   * should yield way better orbitals.
   */
  Timings::timeTaken("Tech. -           DIIS/Damping");
}

template<>
void ConvergenceController<Options::SCF_MODES::UNRESTRICTED>::accelerateConvergence(FockMatrix<UNRESTRICTED>& F,
                                                                                    DensityMatrix<UNRESTRICTED> D) {
  Timings::takeTime("Tech. -           DIIS/Damping");

  if (_first)
    _orthoS.reset(new Eigen::MatrixXd((*_orbitalController->getTransformMatrix(_oneIntController)).transpose()));

  /*
   * Optimizer part: set up the new error vector
   */
  takeTime("Optimizer error measure");
  symInPlace<UNRESTRICTED>(F);
  auto errorVectors = calcFPSminusSPF(F);
  timeTaken(3, "Optimizer error measure");

  /*
   * Get the convergence measure.
   */
  _diisConvMeasure = errorVectors.alpha.maxCoeff();
  // calculate root mean square change of density matrix
  _rmsdOfDensity = calcRMSDofDensity();

  if (errorVectors.beta.maxCoeff() > _diisConvMeasure)
    _diisConvMeasure = errorVectors.beta.maxCoeff();
  /*
   * Check whether DIIS is or should be switched on.
   */
  if (_diisConvMeasure <= _settings.scf.diisStartError && _diis->getNVectorsStored() > 0 && _cycle > 0)
    _diisZoneStart = _cycle;
  _mode = "---";
  if (this->getLevelshift().first[0] > 0.0)
    _mode[2] = 'L';
  if (_diisZoneStart) {
    takeTime("Optimizer update");
    /*
     * Make an Optimizer step: the Fock matrix is updated.
     */
    _mode[1] = 'D';
    if ((errorVectors.alpha.rows() != errorVectors.beta.rows()) || (errorVectors.alpha.cols() != errorVectors.beta.cols())) {
      throw SerenityError("ConvergenceController: The errorVectors alpha/beta dimensions are not identical!");
    }
    else if ((F.alpha.rows() != F.beta.rows()) || (F.alpha.cols() != F.beta.cols())) {
      throw SerenityError("ConvergenceController: The Fock matrix alpha/beta dimensions are not identical!");
    }

    Matrix<double> eTot(errorVectors.alpha.rows(), 2 * errorVectors.alpha.cols());
    Matrix<double> fTot(F.alpha.rows(), 2 * F.alpha.cols());

    eTot << errorVectors.alpha, errorVectors.beta;
    fTot << F.alpha, F.beta;

    _diis->optimize(fTot, eTot);

    F.alpha = fTot.leftCols(F.alpha.cols());
    F.beta = fTot.rightCols(F.beta.cols());

    timeTaken(3, "Optimizer update");
  }
  // rework: only one instance of _damping, but several functions
  if (_diisConvMeasure >= _settings.scf.endDampErr or _cycle < 2) {
    if (_diisConvMeasure <= 10.0 * _settings.scf.diisStartError && not _diisZoneStart) {
      Matrix<double> eTot(errorVectors.alpha.rows(), 2 * errorVectors.alpha.cols());
      Matrix<double> fTot(F.alpha.rows(), 2 * F.alpha.cols());

      eTot << errorVectors.alpha, errorVectors.beta;
      fTot << F.alpha, F.beta;
      _diis->storeMatrix(fTot, eTot);
    }
    if (_settings.scf.damping == Options::DAMPING_ALGORITHMS::SERIES) {
      _damping->arithmeticSeriesDamp(F);
    }
    else if (_settings.scf.damping == Options::DAMPING_ALGORITHMS::DYNAMIC) {
      _damping->dynamicDamp(F, D);
    }
    else if (_settings.scf.damping == Options::DAMPING_ALGORITHMS::STATIC) {
      _damping->staticDamp(F);
    }
    _mode[0] = 'D';
  }
  _first = false;

  /*
   * Done. _F (and thus the SCF procedure) now contains a Fock matrix that
   * should yield way better orbitals.
   */
  Timings::timeTaken("Tech. -           DIIS/Damping");
}

template class ConvergenceController<Options::SCF_MODES::RESTRICTED>;
template class ConvergenceController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
