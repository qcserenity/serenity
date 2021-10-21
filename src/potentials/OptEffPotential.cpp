/**
 * @file OptEffPotential.cpp
 *
 * @date Dec 12, 2016
 * @author David Schnieders
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
#include "potentials/OptEffPotential.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "data/OrbitalController.h"
#include "data/SpinPolarizedData.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/CoulombPotentialOnGridCalculator.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/MOCalculator.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/FockMatrix.h"
#include "geometry/Geometry.h"
#include "grid/GridController.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/looper/TwoElecThreeCenterIntLooper.h"
#include "integrals/wrappers/Libint.h"
#include "math/optimizer/BFGS.h"
#include "math/optimizer/LBFGS.h"
#include "math/optimizer/NewtonRaphson.h"
#include "math/optimizer/SteepestDescent.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <iostream>
#include <vector>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
OptEffPotential<SCFMode>::OptEffPotential(std::shared_ptr<BasisFunctionOnGridController> basFuncOnGridController,
                                          std::shared_ptr<BasisFunctionOnGridController> potBasFuncOnGridController,
                                          std::shared_ptr<OneElectronIntegralController> oneEIntController,
                                          std::shared_ptr<DensityOnGridCalculator<SCFMode>> densOnGridCalculator,
                                          const SpinPolarizedData<SCFMode, unsigned int>& nOccOrbs,
                                          const double smoothFactor, const double singValThreshold, const double exc)
  : Potential<SCFMode>(basFuncOnGridController->getBasisController()),
    _basisFuncOnGridController(basFuncOnGridController),
    _potBasisFuncOnGridController(potBasFuncOnGridController),
    _oneEIntController(oneEIntController),
    _densOnGridCalculator(densOnGridCalculator),
    _potential(nullptr),
    _potentialOnGrid(nullptr),
    _nOccOrbs(nOccOrbs),
    _smoothFactor(smoothFactor),
    _singValThreshold(singValThreshold),
    _exc(exc) {
}

template<>
void OptEffPotential<Options::SCF_MODES::RESTRICTED>::calculateOEP(
    const DensityOnGrid<Options::SCF_MODES::RESTRICTED>& targetDens,
    std::shared_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED>> resultOrbitals,
    const FockMatrix<Options::SCF_MODES::RESTRICTED>& initialGuess) {
  Timings::takeTime("Tech. -      WY Reconstruction");
  auto gridController = _basisFuncOnGridController->getGridController();
  auto& weights = gridController->getWeights();
  DensityOnGrid<Options::SCF_MODES::RESTRICTED> densMatOnGrid(gridController);
  ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED> scalarOpToMat(_basisFuncOnGridController, 0.0);
  _potential.reset(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basisFuncOnGridController->getBasisController()));
  _potentialOnGrid.reset(new GridPotential<Options::SCF_MODES::RESTRICTED>(_basisFuncOnGridController->getGridController()));

  unsigned int cycle = 0;
  double densDiff;
  double oldDensDiff = std::numeric_limits<double>::infinity();

  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd> OEPCoeffs(
      _potBasisFuncOnGridController->getBasisController()->getNBasisFunctions());
  OEPCoeffs.setZero();

  BFGS optimizer(OEPCoeffs);

  /*
   * BFGS optimization
   */
  auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                  std::shared_ptr<Eigen::MatrixXd> hessian, bool print) {
    (void)hessian;
    (void)print;
    updateDensity(densMatOnGrid, resultOrbitals, initialGuess, parameters);

    gradients = -getGradient(targetDens, densMatOnGrid, OEPCoeffs);

    bool converged = false;
    densDiff = (weights.array() * (densMatOnGrid - targetDens).array().abs()).sum();

    value = densDiff;

    if (fabs(densDiff - oldDensDiff) < 1e-7 or densDiff < 1e-5 or cycle > 50) {
      converged = true;
    }

    cycle++;

    oldDensDiff = densDiff;
    return converged;
  };

  optimizer.optimize(updateFunction);

  /*
   * A few Newton-Raphson steps to reach so stationary point
   */
  cycle = 0;
  oldDensDiff = std::numeric_limits<double>::infinity();
  NewtonRaphson optimizer2(OEPCoeffs, _singValThreshold);
  auto const updateFunction2 = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                   std::shared_ptr<Eigen::MatrixXd> hessian, bool print) {
    (void)print;

    updateDensity(densMatOnGrid, resultOrbitals, initialGuess, parameters);

    gradients = getGradient(targetDens, densMatOnGrid, OEPCoeffs);

    *hessian = *getHessian(resultOrbitals);

    bool converged = false;
    densDiff = (weights.array() * (densMatOnGrid - targetDens).array().abs()).sum();
    value = densDiff;

    if (fabs(densDiff - oldDensDiff) < 1e-5 or cycle > 50) {
      converged = true;
      std::cout << "Remaining density difference is " << densDiff << std::endl;
      std::cout << "After cycle: " << cycle << std::endl;
    }

    oldDensDiff = densDiff;
    cycle++;
    return converged;
  };

  optimizer2.optimize(updateFunction2);

  /*
   * Numerical integration
   */
  updateDensity(densMatOnGrid, resultOrbitals, initialGuess, OEPCoeffs);
  scalarOpToMat.addScalarOperatorToMatrix(*_potential, *_potentialOnGrid);

  Timings::timeTaken("Tech. -      WY Reconstruction");
};

template<>
void OptEffPotential<Options::SCF_MODES::UNRESTRICTED>::calculateOEP(
    const DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>& targetDens,
    std::shared_ptr<OrbitalController<Options::SCF_MODES::UNRESTRICTED>> resultOrbitals,
    const FockMatrix<Options::SCF_MODES::UNRESTRICTED>& initialGuess) {
  Timings::takeTime("Tech. -      WY Reconstruction");
  auto gridController = _basisFuncOnGridController->getGridController();
  auto& weights = gridController->getWeights();
  DensityOnGrid<Options::SCF_MODES::UNRESTRICTED> densMatOnGrid(gridController);
  ScalarOperatorToMatrixAdder<Options::SCF_MODES::UNRESTRICTED> scalarOpToMat(_basisFuncOnGridController, 0.0);
  _potential.reset(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basisFuncOnGridController->getBasisController()));
  _potentialOnGrid.reset(new GridPotential<Options::SCF_MODES::UNRESTRICTED>(_basisFuncOnGridController->getGridController()));
  auto nBasFunc = _potBasisFuncOnGridController->getBasisController()->getNBasisFunctions();

  unsigned int cycle = 0;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> densDiff(std::numeric_limits<double>::infinity());
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, double> oldDensDiff(std::numeric_limits<double>::infinity());

  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd> OEPCoeffs(nBasFunc);
  OEPCoeffs.alpha.setZero();
  OEPCoeffs.beta.setZero();

  Eigen::VectorXd joined(2 * nBasFunc);
  joined << OEPCoeffs.alpha, OEPCoeffs.beta;

  BFGS optimizer(joined);

  /*
   * BFGS optimization
   */

  auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                  std::shared_ptr<Eigen::MatrixXd> hessian, bool print) {
    (void)hessian;
    (void)print;
    OEPCoeffs.alpha = parameters.segment(0, nBasFunc);
    OEPCoeffs.beta = parameters.segment(nBasFunc, nBasFunc);
    updateDensity(densMatOnGrid, resultOrbitals, initialGuess, OEPCoeffs);

    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd> tmpGrads =
        getGradient(targetDens, densMatOnGrid, OEPCoeffs);

    Eigen::VectorXd joinedGradients(2 * nBasFunc);
    joinedGradients << -tmpGrads.alpha, -tmpGrads.beta;
    gradients = joinedGradients;

    bool converged = false;
    for_spin(densDiff, densMatOnGrid, targetDens, oldDensDiff) {
      densDiff_spin = (weights.array() * (densMatOnGrid_spin - targetDens_spin).array().abs()).sum();
      value += densDiff_spin;
      if ((fabs(densDiff_spin - oldDensDiff_spin) < 1e-5 or densDiff_spin < 1e-5) or cycle > 50) {
        converged = true;
      }
    };

    cycle++;
    oldDensDiff = densDiff;

    return converged;
  };

  optimizer.optimize(updateFunction);
  OEPCoeffs.alpha = joined.segment(0, nBasFunc);
  OEPCoeffs.beta = joined.segment(nBasFunc, nBasFunc);

  /*
   * Newton--Raphson optimization
   */

  updateDensity(densMatOnGrid, resultOrbitals, initialGuess, OEPCoeffs);
  cycle = 0;
  oldDensDiff = std::numeric_limits<double>::infinity();
  while (true) {
    auto grads = getGradient(targetDens, densMatOnGrid, OEPCoeffs);
    auto hessian = getHessian(resultOrbitals);
    auto& hess = *hessian;
    for_spin(hess, grads, OEPCoeffs) {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(hess_spin, Eigen::ComputeThinU | Eigen::ComputeThinV);
      svd.setThreshold(_singValThreshold);
      Eigen::MatrixXd step = svd.solve(-grads_spin);
      OEPCoeffs_spin += step;
    };
    updateDensity(densMatOnGrid, resultOrbitals, initialGuess, OEPCoeffs);
    for_spin(densDiff, densMatOnGrid, targetDens) {
      densDiff_spin = (weights.array() * (densMatOnGrid_spin - targetDens_spin).array().abs()).sum();
    };
    bool converged = false;
    for_spin(densDiff, oldDensDiff) {
      if (fabs(densDiff_spin - oldDensDiff_spin) < 1e-5 or cycle > 50) {
        converged = true;
        std::cout << "Remaining density difference is " << densDiff_spin << std::endl;
        std::cout << "After cycle: " << cycle << std::endl;
      }
    };
    oldDensDiff = densDiff;

    oldDensDiff = densDiff;
    if (converged)
      break;

    cycle++;
  }

  /*
   * Numerical integration
   */
  scalarOpToMat.addScalarOperatorToMatrix(*_potential, *_potentialOnGrid);
  Timings::timeTaken("Tech. -      WY Reconstruction");
};

template<Options::SCF_MODES SCFMode>
void OptEffPotential<SCFMode>::calculateOEP(const DensityOnGrid<SCFMode>& targetDens,
                                            std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                                            const GridPotential<SCFMode>& initialGuess) {
  ScalarOperatorToMatrixAdder<SCFMode> scalarOpToMat(_basisFuncOnGridController, 0.0);

  FockMatrix<SCFMode> initialGuessMat(this->_basis);

  scalarOpToMat.addScalarOperatorToMatrix(initialGuessMat, initialGuess);

  calculateOEP(targetDens, resultOrbitals, initialGuessMat);
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd OptEffPotential<SCFMode>::getGeomGradients() {
  Eigen::MatrixXd gradientContr(1, 3);
  gradientContr.setZero();
  // ToDo throw error
  return gradientContr;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, double>
OptEffPotential<SCFMode>::calculateOEPLB(const DensityOnGrid<SCFMode>& targetDens,
                                         std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                                         const FockMatrix<SCFMode>& initialGuess, double damping, unsigned int maxCycles) {
  Timings::takeTime("Tech. -      LB Reconstruction");

  auto gridController(_basisFuncOnGridController->getGridController());
  auto nBlocks(_basisFuncOnGridController->getNBlocks());
  auto nGridPoints(gridController->getNGridPoints());
  auto& weights = gridController->getWeights();

  SpinPolarizedData<SCFMode, double> oldDensDiff(std::numeric_limits<double>::infinity());

  DensityMatrixController<SCFMode> densMatUpdater(resultOrbitals, _nOccOrbs);
  DensityMatrix<SCFMode> densMat = densMatUpdater.getDensityMatrix();
  auto basisController = _basisFuncOnGridController->getBasisController();

  ScalarOperatorToMatrixAdder<SCFMode> scalarOpToMat(_basisFuncOnGridController, 0.0);
  FockMatrix<SCFMode> fockMat(basisController);

  auto& libint = Libint::getInstance();
  auto aoKinInts = libint.compute1eInts(LIBINT_OPERATOR::kinetic, basisController);

  _potential.reset(new FockMatrix<SCFMode>(_basisFuncOnGridController->getBasisController()));
  _potentialOnGrid.reset(new GridPotential<SCFMode>(_basisFuncOnGridController->getGridController()));

  auto& oepGrid = *_potentialOnGrid;

  /*
   * Set grid to optimize
   */

#pragma omp for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    unsigned int blockEnd;
    if (block == (nBlocks - 1)) {
      blockEnd = nGridPoints;
    }
    else {
      blockEnd = _basisFuncOnGridController->getFirstIndexOfBlock(block + 1);
    }
    auto blockStart = _basisFuncOnGridController->getFirstIndexOfBlock(block);
    for (unsigned int gridPt = blockStart; gridPt < blockEnd; gridPt++) {
      for_spin(oepGrid) {
        oepGrid_spin[gridPt] = 1.0;
      };
    }
  }

  for_spin(fockMat) {
    fockMat_spin = aoKinInts;
  };
  fockMat += initialGuess;
  scalarOpToMat.addScalarOperatorToMatrix(fockMat, oepGrid);
  resultOrbitals->updateOrbitals(fockMat, _oneEIntController);
  densMat = densMatUpdater.getDensityMatrix();
  auto densMatOnGrid = _densOnGridCalculator->calcDensityOnGrid(densMat);

  SpinPolarizedData<SCFMode, double> finalDensDiff;

  unsigned int cycle = 0;
  while (true) {
    SpinPolarizedData<SCFMode, Eigen::VectorXd> densDiff(omp_get_max_threads());

    for_spin(densDiff) {
      densDiff_spin.setZero();
    };

      /*
       * Scale the potential
       */
#pragma omp for schedule(dynamic)
    for (unsigned int block = 0; block < nBlocks; block++) {
      unsigned int blockEnd;
      if (block == (nBlocks - 1)) {
        blockEnd = nGridPoints;
      }
      else {
        blockEnd = _basisFuncOnGridController->getFirstIndexOfBlock(block + 1);
      }
      auto blockStart = _basisFuncOnGridController->getFirstIndexOfBlock(block);
      for (unsigned int gridPt = blockStart; gridPt < blockEnd; gridPt++) {
        for_spin(densMatOnGrid, targetDens, oepGrid, densDiff) {
          if (targetDens_spin[gridPt] < 1e-9) {
            double tmp = oepGrid_spin[gridPt];
            oepGrid_spin[gridPt] *= (1 - damping) * densMatOnGrid_spin[gridPt] / 1e-9;
            oepGrid_spin[gridPt] += damping * tmp;
            densDiff_spin[omp_get_thread_num()] +=
                weights[gridPt] * fabs(densMatOnGrid_spin[gridPt] - targetDens_spin[gridPt]);
          }
          else {
            double tmp = oepGrid_spin[gridPt];
            oepGrid_spin[gridPt] *= (1 - damping) * densMatOnGrid_spin[gridPt] / targetDens_spin[gridPt];
            oepGrid_spin[gridPt] += damping * tmp;
            densDiff_spin[omp_get_thread_num()] +=
                weights[gridPt] * fabs(densMatOnGrid_spin[gridPt] - targetDens_spin[gridPt]);
          }
        };
      }
    }
    for_spin(fockMat) {
      fockMat_spin = aoKinInts;
    };
    fockMat += initialGuess;
    scalarOpToMat.addScalarOperatorToMatrix(fockMat, oepGrid);
    resultOrbitals->updateOrbitals(fockMat, _oneEIntController);
    densMat = densMatUpdater.getDensityMatrix();
    densMatOnGrid = _densOnGridCalculator->calcDensityOnGrid(densMat);

    cycle++;
    bool converged = false;
    for_spin(densDiff, oldDensDiff) {
      if (fabs(oldDensDiff_spin - densDiff_spin.sum()) < 1e-7 or densDiff_spin.sum() < 1e-5 or cycle > maxCycles) {
        std::cout << "Remaining density difference is " << densDiff_spin.sum() << std::endl;
        std::cout << "After cycle: " << cycle << std::endl;
        converged = true;
      }
      oldDensDiff_spin = densDiff_spin.sum();
    };
    if (converged) {
      for_spin(finalDensDiff, densDiff) {
        finalDensDiff_spin = densDiff_spin.sum();
      };
      break;
    }
  }

  /*
   * Subtract initial value
   */

#pragma omp for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    unsigned int blockEnd;
    if (block == (nBlocks - 1)) {
      blockEnd = nGridPoints;
    }
    else {
      blockEnd = _basisFuncOnGridController->getFirstIndexOfBlock(block + 1);
    }
    auto blockStart = _basisFuncOnGridController->getFirstIndexOfBlock(block);
    for (unsigned int gridPt = blockStart; gridPt < blockEnd; gridPt++) {
      for_spin(oepGrid) {
        oepGrid_spin[gridPt] -= 1.0;
      };
    }
  }

  /*
   * Numerical integration
   */
  scalarOpToMat.addScalarOperatorToMatrix(*_potential, oepGrid);

  Timings::timeTaken("Tech. -      LB Reconstruction");
  return finalDensDiff;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, double>
OptEffPotential<SCFMode>::calculateOEPLB(const DensityOnGrid<SCFMode>& targetDens,
                                         std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                                         const GridPotential<SCFMode>& initialGuess, double damping, unsigned int maxCycles) {
  ScalarOperatorToMatrixAdder<SCFMode> scalarOpToMat(_basisFuncOnGridController, 0.0);

  FockMatrix<SCFMode> initialGuessMat(this->_basis);

  scalarOpToMat.addScalarOperatorToMatrix(initialGuessMat, initialGuess);

  return calculateOEPLB(targetDens, resultOrbitals, initialGuessMat, damping, maxCycles);
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd>
OptEffPotential<SCFMode>::getGradient(const DensityOnGrid<SCFMode>& targetDens, const DensityOnGrid<SCFMode>& densMatOnGrid,
                                      SpinPolarizedData<SCFMode, Eigen::VectorXd>& OEPCoeffs) {
  const unsigned int nThreads = omp_get_max_threads();
  auto nBasisFunc = _potBasisFuncOnGridController->getBasisController()->getNBasisFunctions();
  auto gridController = _potBasisFuncOnGridController->getGridController();
  auto nGridPoints = gridController->getNGridPoints();
  auto& weights = gridController->getWeights();
  auto nBlocks = _potBasisFuncOnGridController->getNBlocks();
  /*
   * Get Libint
   */
  auto& libint = Libint::getInstance();
  /*
   * Get kinetic integrals for smoothing constraint
   */
  Eigen::MatrixXd aoKineticIntegrals(nBasisFunc, nBasisFunc);
  aoKineticIntegrals = libint.compute1eInts(LIBINT_OPERATOR::kinetic, _potBasisFuncOnGridController->getBasisController());
  /*
   * Prepare gradient
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> gradient(nBasisFunc);
  for_spin(gradient) {
    gradient_spin.setZero();
  };

  /*
   * Thread safety: Generate a gradient vector for every thread
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>> gradVec(nThreads);

#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < nThreads; i++) {
    SpinPolarizedData<SCFMode, Eigen::VectorXd> tmp(nBasisFunc);
    for_spin(tmp) {
      tmp_spin.setZero();
    };
    gradVec[i] = tmp;
  }

  /*
   * Add OEP to the Lagrangian and build gradient
   */
#pragma omp for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    auto& grad = gradVec[omp_get_thread_num()];
    auto basisFuncOnGrid = _potBasisFuncOnGridController->calculateBasisFunctionData(block);
    unsigned int blockEnd;
    if (block == (nBlocks - 1)) {
      blockEnd = nGridPoints;
    }
    else {
      blockEnd = _potBasisFuncOnGridController->getFirstIndexOfBlock(block + 1);
    }
    auto blockStart = _potBasisFuncOnGridController->getFirstIndexOfBlock(block);
    for (unsigned int k = 0; k < nBasisFunc; k++) {
      for_spin(densMatOnGrid, targetDens, grad) {
        grad_spin[k] +=
            (weights.segment(blockStart, blockEnd - blockStart).array() * basisFuncOnGrid->functionValues.col(k).array() *
             (densMatOnGrid_spin.segment(blockStart, blockEnd - blockStart).array() -
              targetDens_spin.segment(blockStart, blockEnd - blockStart).array()))
                .sum();
      };
    }
  }
  /*
   * Sum up threads
   */
  for (unsigned int i = 0; i < nThreads; ++i) {
    gradient += gradVec[i];
  }

  for_spin(gradient, OEPCoeffs) {
    gradient_spin -= (aoKineticIntegrals * OEPCoeffs_spin).eval() * 4 * _smoothFactor;
  };

  return gradient;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>
OptEffPotential<SCFMode>::getHessian(std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals) {
  const unsigned int nThreads = omp_get_max_threads();
  auto nBasisFunc = _potBasisFuncOnGridController->getBasisController()->getNBasisFunctions();
  auto gridController = _potBasisFuncOnGridController->getGridController();
  auto& weights = gridController->getWeights();
  auto nBlocks = _potBasisFuncOnGridController->getNBlocks();

  auto& libint = Libint::getInstance();
  /*
   * Get kinetic integrals for smoothing contraint
   */
  Matrix<double> aoKineticIntegrals(nBasisFunc, nBasisFunc);
  auto aoKinInts = libint.compute1eInts(LIBINT_OPERATOR::kinetic, _potBasisFuncOnGridController->getBasisController());
  /*
   * Prepare Hessian
   */
  auto hessianPtr = std::make_shared<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>(nBasisFunc, nBasisFunc);
  auto& hessian = *hessianPtr;
  for_spin(hessian) {
    hessian_spin.setZero();
  };

  SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>> threeCenterInts(nBasisFunc);

#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < nBasisFunc; i++) {
    for_spin(threeCenterInts, _nOccOrbs) {
      Eigen::MatrixXd tmp(_nOccOrbs_spin,
                          _basisFuncOnGridController->getBasisController()->getNBasisFunctions() - _nOccOrbs_spin);
      tmp.setZero();
      threeCenterInts_spin[i] = tmp;
    };
  }
  SpinPolarizedData<SCFMode, Eigen::VectorXd> orbEnergies = resultOrbitals->getEigenvalues();
  MOCalculator moCalc(_basisFuncOnGridController);
  CoefficientMatrix<SCFMode> coefficients = resultOrbitals->getCoefficients();
  auto MOsOnGrid = moCalc.calcAllMOValuesOnGrid<SCFMode>(coefficients);
  std::vector<SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>>> intsCore(nThreads);

  /*
   * Create a vector of threeCenterInts for every thread
   */
#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < nThreads; ++i) {
    SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>> tmpVec(nBasisFunc);
    for (unsigned int j = 0; j < nBasisFunc; j++) {
      for_spin(tmpVec, _nOccOrbs) {
        Eigen::MatrixXd tmpMat(_nOccOrbs_spin,
                               _basisFuncOnGridController->getBasisController()->getNBasisFunctions() - _nOccOrbs_spin);
        tmpMat.setZero();
        tmpVec_spin[j] = tmpMat;
      };
    }
    intsCore[i] = tmpVec;
  }

  /*
   * Calculation of three center integrals
   */
#pragma omp for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    auto& ints = intsCore[omp_get_thread_num()];
    auto basisFuncOnGrid = _potBasisFuncOnGridController->calculateBasisFunctionData(block);
    auto blockStart = _potBasisFuncOnGridController->getFirstIndexOfBlock(block);
    auto nBlock = basisFuncOnGrid->functionValues.rows();
    for (unsigned int basFuncT = 0; basFuncT < nBasisFunc; basFuncT++) {
      for_spin(_nOccOrbs, MOsOnGrid, ints) {
        for (unsigned int occMo = 0; occMo < _nOccOrbs_spin; occMo++) {
          Eigen::VectorXd intermediateProd(nBlock);
          intermediateProd.setZero();
          intermediateProd.array() = MOsOnGrid_spin.col(occMo).segment(blockStart, nBlock).array() *
                                     basisFuncOnGrid->functionValues.col(basFuncT).array();
          for (unsigned int virtMO = 0;
               virtMO < _basisFuncOnGridController->getBasisController()->getNBasisFunctions() - _nOccOrbs_spin; virtMO++) {
            ints_spin[basFuncT](occMo, virtMO) +=
                (weights.segment(blockStart, nBlock).array() * intermediateProd.array() *
                 MOsOnGrid_spin.col(virtMO + _nOccOrbs_spin).segment(blockStart, nBlock).array())
                    .sum();
          }
        }
      };
    }
  }

  /*
   * Sum up threads
   */
  for (unsigned int i = 0; i < nThreads; ++i) {
    auto& ints = intsCore[i];
    for_spin(ints, threeCenterInts) {
      for (unsigned int j = 0; j < nBasisFunc; j++) {
        threeCenterInts_spin[j] += ints_spin[j];
      }
    };
  }

  /*
   * Calculate Hessian
   */
  for_spin(_nOccOrbs, hessian, orbEnergies, threeCenterInts) {
#pragma omp for schedule(dynamic)
    for (unsigned int basFuncU = 0; basFuncU < nBasisFunc; basFuncU++) {
      for (unsigned int occMo = 0; occMo < _nOccOrbs_spin; occMo++) {
        for (unsigned int virtMo = 0;
             virtMo < _basisFuncOnGridController->getBasisController()->getNBasisFunctions() - _nOccOrbs_spin; virtMo++) {
          double intermediateProd = 4 * threeCenterInts_spin[basFuncU](occMo, virtMo) /
                                    (orbEnergies_spin[occMo] - orbEnergies_spin[virtMo + _nOccOrbs_spin]);
          for (unsigned int basFuncT = 0; basFuncT < nBasisFunc; basFuncT++) {
            hessian_spin(basFuncU, basFuncT) += threeCenterInts_spin[basFuncT](occMo, virtMo) * intermediateProd;
          }
        }
      }
    }
    hessian_spin -= 4 * _smoothFactor * aoKinInts;
  };

  return hessianPtr;
}

template<Options::SCF_MODES SCFMode>
void OptEffPotential<SCFMode>::updateDensity(DensityOnGrid<SCFMode>& densMatOnGrid,
                                             std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                                             const FockMatrix<SCFMode>& initialGuess,
                                             const SpinPolarizedData<SCFMode, Eigen::VectorXd>& OEPCoeffs) {
  /*
   * Get Libint
   */
  auto& libint = Libint::getInstance();

  /*
   * Get stuff related to basisset
   */
  auto basisController(_basisFuncOnGridController->getBasisController());
  auto nBasisFunc(basisController->getNBasisFunctions());
  /*
   * Get stuff related to grid
   */
  auto gridController(_basisFuncOnGridController->getGridController());
  auto nBlocks(_basisFuncOnGridController->getNBlocks());

  /*
   * Get stuff related to the reconstructed orbitalset
   */
  FockMatrix<SCFMode> fockMat(basisController);
  DensityMatrixController<SCFMode> densMatUpdater(resultOrbitals, _nOccOrbs);
  DensityMatrix<SCFMode> densMat = densMatUpdater.getDensityMatrix();
  /*
   * Further tools
   */
  ScalarOperatorToMatrixAdder<SCFMode> scalarOpToMat(_basisFuncOnGridController, 0.0);
  /*
   * Get kinetic integrals for Wu-Yang Lagrangian and the Fock-Matrix generation
   * of the reconstructed system
   */
  Matrix<double> aoKineticIntegrals(nBasisFunc, nBasisFunc);
  aoKineticIntegrals = libint.compute1eInts(LIBINT_OPERATOR::kinetic, basisController);

  GridPotential<SCFMode>& OEPGrid = *_potentialOnGrid;
  for_spin(OEPGrid) {
    OEPGrid_spin.setZero();
  };
    /*
     * Generate new grid representation of the OEP
     */
    // ToDo add prescreening (see MatrixoperatorToGridTransformer.cpp)
#pragma omp for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    auto basisFuncOnGrid = _potBasisFuncOnGridController->calculateBasisFunctionData(block);
    auto blockStart = _potBasisFuncOnGridController->getFirstIndexOfBlock(block);
    auto nBlock = basisFuncOnGrid->functionValues.rows();
    for_spin(OEPGrid, OEPCoeffs) {
      OEPGrid_spin.segment(blockStart, nBlock) = basisFuncOnGrid->functionValues * OEPCoeffs_spin;
    };
  }

  /*
   * Build up new Fock matrix and density of the
   * reconstructed system
   */
  for_spin(fockMat) {
    fockMat_spin = aoKineticIntegrals;
  };
  scalarOpToMat.addScalarOperatorToMatrix(fockMat, OEPGrid);
  fockMat += initialGuess;
  if (_exc > 0.0)
    fockMat += getX(densMat);
  resultOrbitals->updateOrbitals(fockMat, _oneEIntController);
  densMat = densMatUpdater.getDensityMatrix();
  densMatOnGrid = _densOnGridCalculator->calcDensityOnGrid(densMat);
}

// ToDo replace by ExchangePotential.h/.cpp
template<>
FockMatrix<RESTRICTED> OptEffPotential<RESTRICTED>::getX(DensityMatrix<RESTRICTED>& densityMatrix) {
  const unsigned int nBFs = _basisFuncOnGridController->getBasisController()->getNBasisFunctions();

  FockMatrix<RESTRICTED> xMat(_basisFuncOnGridController->getBasisController());

  /*
   * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
   */
  std::vector<FockMatrix<Options::SCF_MODES::RESTRICTED>*> fx;
  const unsigned int nThreads = omp_get_max_threads();
  for (unsigned int i = 0; i < nThreads; ++i) {
    fx.push_back(new FockMatrix<Options::SCF_MODES::RESTRICTED>(_basisFuncOnGridController->getBasisController()));
    fx[i]->setZero();
  }
  /*
   * Function which parses the integrals
   */
  auto distribute = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, const double integral,
                        unsigned int threadId) {
    /*
     * Exchange
     */
    const double exc = integral * 0.5 * _exc;
    const double exc1 = *(densityMatrix.data() + i * nBFs + k) * exc;
    const double exc2 = *(densityMatrix.data() + i * nBFs + l) * exc;
    const double exc3 = *(densityMatrix.data() + j * nBFs + k) * exc;
    const double exc4 = *(densityMatrix.data() + j * nBFs + l) * exc;
    *(fx[threadId]->data() + j * nBFs + l) -= exc1;
    *(fx[threadId]->data() + l * nBFs + j) -= exc1;
    *(fx[threadId]->data() + j * nBFs + k) -= exc2;
    *(fx[threadId]->data() + k * nBFs + j) -= exc2;
    *(fx[threadId]->data() + i * nBFs + l) -= exc3;
    *(fx[threadId]->data() + l * nBFs + i) -= exc3;
    *(fx[threadId]->data() + i * nBFs + k) -= exc4;
    *(fx[threadId]->data() + k * nBFs + i) -= exc4;
  };
  /*
   * Detailed prescreening function
   */
  // Maximum absolute value in densityMatrix
  const auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  const auto maxDens = maxDensMat.maxCoeff();
  auto prescreeningFunc = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, double schwartz) {
    /*
     * Early return for insignificance based on the largest element in the whole density matrix
     */
    if (maxDens * schwartz < 1e-12)
      return true;
    double maxDBlock = maxDensMat(i, k);
    maxDBlock = std::max(maxDBlock, maxDensMat(i, l));
    maxDBlock = std::max(maxDBlock, maxDensMat(j, k));
    maxDBlock = std::max(maxDBlock, maxDensMat(j, l));
    if (maxDBlock * schwartz < 1e-12)
      return true;
    /*
     * Shell quadruple is significant.
     */
    return false;
  };
  /*
   * Construct the looper, which loops over all integrals
   */
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basisFuncOnGridController->getBasisController(), 1e-12);
  /*
   * Run
   */
  looper.loopNoDerivative(distribute, prescreeningFunc, maxDens, nullptr, true);
  for (unsigned int i = 0; i < nThreads; ++i) {
    xMat += *fx[i];
    delete fx[i];
  }
  return xMat;
}

// ToDo replace by ExchangePotential.h/.cpp
template<>
FockMatrix<UNRESTRICTED> OptEffPotential<UNRESTRICTED>::getX(DensityMatrix<UNRESTRICTED>& densityMatrix) {
  FockMatrix<UNRESTRICTED> xMat(_basisFuncOnGridController->getBasisController());

  /*
   * Thread safety issues; create one (partial) Fock matrix for each thread, sum up in the end.
   * ToDo: Don't use raw pointer here.
   */
  std::vector<FockMatrix<Options::SCF_MODES::UNRESTRICTED>*> fx;
  const unsigned int nThreads = omp_get_max_threads();
  for (unsigned int i = 0; i < nThreads; ++i) {
    fx.push_back(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(_basisFuncOnGridController->getBasisController()));
  }
  /*
   * Function which parses the integrals
   */
  auto distribute = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, const double integral,
                        unsigned int threadId) {
    /*
     * Exchange
     */
    const double exc = integral * _exc;

    const double exc1a = densityMatrix.alpha(i, k) * exc;
    const double exc2a = densityMatrix.alpha(i, l) * exc;
    const double exc3a = densityMatrix.alpha(j, k) * exc;
    const double exc4a = densityMatrix.alpha(j, l) * exc;
    fx[threadId]->alpha(j, l) -= exc1a;
    fx[threadId]->alpha(l, j) -= exc1a;
    fx[threadId]->alpha(j, k) -= exc2a;
    fx[threadId]->alpha(k, j) -= exc2a;
    fx[threadId]->alpha(i, l) -= exc3a;
    fx[threadId]->alpha(l, i) -= exc3a;
    fx[threadId]->alpha(i, k) -= exc4a;
    fx[threadId]->alpha(k, i) -= exc4a;
    const double exc1b = densityMatrix.beta(i, k) * exc;
    const double exc2b = densityMatrix.beta(i, l) * exc;
    const double exc3b = densityMatrix.beta(j, k) * exc;
    const double exc4b = densityMatrix.beta(j, l) * exc;
    fx[threadId]->beta(j, l) -= exc1b;
    fx[threadId]->beta(l, j) -= exc1b;
    fx[threadId]->beta(j, k) -= exc2b;
    fx[threadId]->beta(k, j) -= exc2b;
    fx[threadId]->beta(i, l) -= exc3b;
    fx[threadId]->beta(l, i) -= exc3b;
    fx[threadId]->beta(i, k) -= exc4b;
    fx[threadId]->beta(k, i) -= exc4b;
  };
  /*
   * Detailed prescreening function
   */
  // Maximum absolute value in densityMatrix
  const auto maxDensMat = densityMatrix.shellWiseAbsMax().total();
  const double maxDens = maxDensMat.maxCoeff();
  auto prescreeningFunc = [&](const unsigned& i, const unsigned& j, const unsigned& k, const unsigned& l, const double& schwartz) {
    /*
     * Early return for insignificance based on the largest element in the whole density matrix
     */
    if (maxDens * schwartz < 1e-12)
      return true;
    double maxDBlock = maxDensMat(i, j);
    maxDBlock = std::max(maxDBlock, maxDensMat(i, k));
    maxDBlock = std::max(maxDBlock, maxDensMat(i, l));
    maxDBlock = std::max(maxDBlock, maxDensMat(j, k));
    maxDBlock = std::max(maxDBlock, maxDensMat(j, l));
    maxDBlock = std::max(maxDBlock, maxDensMat(k, l));
    if (maxDBlock * schwartz < 1e-12)
      return true;
    /*
     * Shell quadruple is significant.
     */
    return false;
  };
  /*
   * Construct the looper, which loops over all integrals
   */
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basisFuncOnGridController->getBasisController(), 1e-12);
  /*
   * Run
   */
  looper.loopNoDerivative(distribute, prescreeningFunc, maxDens, nullptr, true);
  for (unsigned int i = 0; i < nThreads; ++i) {
    xMat += *fx[i];
    delete fx[i];
  }
  return xMat;
}

template<>
void OptEffPotential<Options::SCF_MODES::RESTRICTED>::calculateOEPCarter(
    const DensityMatrix<Options::SCF_MODES::RESTRICTED>& targetDens,
    std::shared_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED>> resultOrbitals,
    const FockMatrix<Options::SCF_MODES::RESTRICTED>& initialGuess, const unsigned int maxCycles) {
  Timings::takeTime("Tech. -      ZC Reconstruction");

  /*
   * Get Libint
   */
  auto& libint = Libint::getInstance();

  /*
   * Get stuff related to basis set
   */
  auto basisController = _basisFuncOnGridController->getBasisController();
  auto nBasisFunc = basisController->getNBasisFunctions();

  _potential.reset(new FockMatrix<Options::SCF_MODES::RESTRICTED>(basisController));
  auto& pot = *_potential;

  /*
   * Get kinetic integrals for Wu-Yang Lagrangian and the Fock-Matrix generation
   * of the reconstructed system
   */
  Matrix<double> aoKineticIntegrals(nBasisFunc, nBasisFunc);
  aoKineticIntegrals = libint.compute1eInts(LIBINT_OPERATOR::kinetic, basisController);

  /*
   * Get stuff related to the reconstructed orbitalset
   */
  FockMatrix<Options::SCF_MODES::RESTRICTED> fockMat(basisController);
  DensityMatrixController<Options::SCF_MODES::RESTRICTED> densMatUpdater(resultOrbitals, _nOccOrbs);
  DensityMatrix<Options::SCF_MODES::RESTRICTED> densMat(basisController);
  densMat = densMatUpdater.getDensityMatrix();

  Eigen::VectorXd coeffVec = Eigen::Map<Eigen::VectorXd>(pot.data(), nBasisFunc * nBasisFunc);

  LBFGS optimizer(coeffVec);

  unsigned int cycle = 0;
  double densDiffOld = std::numeric_limits<double>::infinity();
  /*
   * BFGS optimization
   */
  auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                  std::shared_ptr<Eigen::MatrixXd> hessian, bool print) {
    (void)hessian;
    (void)print;
    (void)parameters;

    /*
     * Build up new Fock matrix and density of the
     * reconstructed system
     */
    for_spin(fockMat) {
      fockMat_spin = aoKineticIntegrals;
    };
    fockMat += initialGuess;
    pot = Eigen::Map<Eigen::MatrixXd>(coeffVec.data(), nBasisFunc, nBasisFunc);
    fockMat += pot;
    if (_exc > 0.0)
      fockMat += getX(densMat);
    resultOrbitals->updateOrbitals(fockMat, _oneEIntController);
    densMat = densMatUpdater.getDensityMatrix();

    SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd> gradient(nBasisFunc * nBasisFunc);

    DensityMatrix<Options::SCF_MODES::RESTRICTED> targetCopy(targetDens);
    gradients = Eigen::Map<Eigen::VectorXd>(targetCopy.data(), nBasisFunc * nBasisFunc);
    gradients -= Eigen::Map<Eigen::VectorXd>(densMat.data(), nBasisFunc * nBasisFunc);

    double densDiff = 0.0;

    densDiff += gradients.array().abs().sum();

    bool converged = false;
    if (cycle > maxCycles or densDiff < 1e-8 or fabs(densDiffOld - densDiff) < 1e-10) {
      std::cout << "Remaining density difference is " << densDiff << std::endl;
      std::cout << "After cycle: " << cycle << std::endl;
      converged = true;
    }

    densDiffOld = densDiff;

    value = densDiff;

    cycle++;

    return converged;
  };

  optimizer.optimize(updateFunction);

  Timings::timeTaken("Tech. -      ZC Reconstruction");
}

template<>
void OptEffPotential<Options::SCF_MODES::UNRESTRICTED>::calculateOEPCarter(
    const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>& targetDens,
    std::shared_ptr<OrbitalController<Options::SCF_MODES::UNRESTRICTED>> resultOrbitals,
    const FockMatrix<Options::SCF_MODES::UNRESTRICTED>& initialGuess, const unsigned int maxCycles) {
  Timings::takeTime("Tech. -      ZC Reconstruction");

  /*
   * Get Libint
   */
  auto& libint = Libint::getInstance();

  /*
   * Get stuff related to basisset
   */
  auto basisController = _basisFuncOnGridController->getBasisController();
  auto nBasisFunc = basisController->getNBasisFunctions();

  _potential.reset(new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(basisController));
  auto& pot = *_potential;

  /*
   * Get kinetic integrals for Wu-Yang Lagrangian and the Fock-Matrix generation
   * of the reconstructed system
   */
  Matrix<double> aoKineticIntegrals(nBasisFunc, nBasisFunc);
  aoKineticIntegrals = libint.compute1eInts(LIBINT_OPERATOR::kinetic, basisController);

  /*
   * Get stuff related to the reconstructed orbitalset
   */
  FockMatrix<Options::SCF_MODES::UNRESTRICTED> fockMat(basisController);
  DensityMatrixController<Options::SCF_MODES::UNRESTRICTED> densMatUpdater(resultOrbitals, _nOccOrbs);
  DensityMatrix<Options::SCF_MODES::UNRESTRICTED> densMat(basisController);
  densMat = densMatUpdater.getDensityMatrix();

  Eigen::VectorXd coeffs(2 * nBasisFunc * nBasisFunc);
  coeffs.segment(0, nBasisFunc * nBasisFunc) = Eigen::Map<Eigen::VectorXd>(pot.alpha.data(), nBasisFunc * nBasisFunc);
  coeffs.segment(nBasisFunc * nBasisFunc, nBasisFunc * nBasisFunc) =
      Eigen::Map<Eigen::VectorXd>(pot.beta.data(), nBasisFunc * nBasisFunc);

  LBFGS optimizer(coeffs);

  unsigned int cycle = 0;
  double densDiffOld = std::numeric_limits<double>::infinity();
  /*
   * BFGS optimization
   */
  auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                  std::shared_ptr<Eigen::MatrixXd> hessian, bool print) {
    (void)hessian;
    (void)print;
    (void)parameters;

    /*
     * Build up new Fock matrix and density of the
     * reconstructed system
     */
    for_spin(fockMat) {
      fockMat_spin = aoKineticIntegrals;
    };
    fockMat += initialGuess;
    pot.alpha = Eigen::Map<Eigen::MatrixXd>(coeffs.segment(0, nBasisFunc * nBasisFunc).data(), nBasisFunc, nBasisFunc);
    pot.beta = Eigen::Map<Eigen::MatrixXd>(coeffs.segment(nBasisFunc * nBasisFunc, nBasisFunc * nBasisFunc).data(),
                                           nBasisFunc, nBasisFunc);
    fockMat.alpha += pot.alpha;
    fockMat.beta += pot.beta;
    if (_exc > 0.0)
      fockMat += getX(densMat);
    resultOrbitals->updateOrbitals(fockMat, _oneEIntController);
    densMat = densMatUpdater.getDensityMatrix();

    DensityMatrix<Options::SCF_MODES::UNRESTRICTED> targetCopy(targetDens);
    gradients.segment(0, nBasisFunc * nBasisFunc) =
        Eigen::Map<Eigen::VectorXd>(targetCopy.alpha.data(), nBasisFunc * nBasisFunc);
    gradients.segment(nBasisFunc * nBasisFunc, nBasisFunc * nBasisFunc) =
        Eigen::Map<Eigen::VectorXd>(targetCopy.beta.data(), nBasisFunc * nBasisFunc);
    gradients.segment(0, nBasisFunc * nBasisFunc) -=
        Eigen::Map<Eigen::VectorXd>(densMat.alpha.data(), nBasisFunc * nBasisFunc);
    gradients.segment(nBasisFunc * nBasisFunc, nBasisFunc * nBasisFunc) -=
        Eigen::Map<Eigen::VectorXd>(densMat.beta.data(), nBasisFunc * nBasisFunc);

    double densDiff = 0.0;

    densDiff += gradients.array().abs().sum();

    bool converged = false;
    if (cycle > maxCycles or densDiff < 1e-8 or fabs(densDiffOld - densDiff) < 1e-10) {
      std::cout << "Remaining density difference is " << densDiff << std::endl;
      std::cout << "After cycle: " << cycle << std::endl;
      converged = true;
    }
    densDiffOld = densDiff;

    value = densDiff;

    cycle++;

    return converged;
  };

  optimizer.optimize(updateFunction);

  Timings::timeTaken("Tech. -      ZC Reconstruction");
}

template class OptEffPotential<Options::SCF_MODES::RESTRICTED>;
template class OptEffPotential<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
