/**
 * @file FuncPotential.cpp
 *
 * @date Nov 24, 2016
 * @author: Jan Unsleber
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
#include "potentials/FuncPotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "dft/functionals/FunctionalLibrary.h"
#include "dft/functionals/wrappers/PartialDerivatives.h"
#include "misc/Timing.h"
#include "potentials/SAOPPotential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FuncPotential<SCFMode>::FuncPotential(std::shared_ptr<SystemController> system,
                                      std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                                      std::shared_ptr<GridController> grid, Functional functional)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _system(system),
    _dMatController(dMat),
    _grid(grid),
    _functional(functional),
    _potential(nullptr),
    _gridToMatrix(nullptr),
    _densOnGridController(nullptr),
    _basisFunctionOnGridController(nullptr),
    _energy(0.0) {
  this->_grid->addSensitiveObject(ObjectSensitiveClass<Grid>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);

  _basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(system->getSettings(), this->_basis, this->_grid);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
      _basisFunctionOnGridController, system->getSettings().grid.blockAveThreshold);
  _densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densOnGridCalculator, _dMatController);
  _gridToMatrix = std::make_shared<ScalarOperatorToMatrixAdder<SCFMode>>(_basisFunctionOnGridController,
                                                                         system->getSettings().grid.blockAveThreshold);
};

template<Options::SCF_MODES SCFMode>
double FuncPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& /*P*/) {
  if (!_potential)
    getMatrix();
  return _energy;
};

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& FuncPotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System - Functional Pot.");
  if (!_potential) {
    if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::MODELL) {
      auto system = _system.lock();
      SAOPPotential<SCFMode> saopCalc(this->_basis, system, _densOnGridController);
      _potential.reset(new FockMatrix<SCFMode>(this->_basis));
      *_potential = saopCalc.getMatrix();
      _energy = saopCalc.getEnergy(_dMatController->getDensityMatrix());
    }
    else {
      FunctionalLibrary<SCFMode> flib(128);
      auto funcData = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, _functional, _densOnGridController);

      _potential.reset(new FockMatrix<SCFMode>(this->_basis));
      auto& pot = *_potential;
      for_spin(pot) {
        pot_spin.setZero();
      };

      if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::NONE) {
        // Nothing to be done here
      }
      else if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::LDA) {
        _gridToMatrix->addScalarOperatorToMatrix(pot, *funcData.dFdRho);
      }
      else if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
        _gridToMatrix->addScalarOperatorToMatrix(pot, *funcData.dFdRho, *funcData.dFdGradRho);
      }
      else {
        assert(false && "Unsupported functional type as functional!");
      }
      _energy = funcData.energy;
    }
  }
  Timings::timeTaken("Active System - Functional Pot.");
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd FuncPotential<SCFMode>::getGeomGradients() {
  FunctionalLibrary<SCFMode> flib(128);
  auto funcData = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, _functional, _densOnGridController, 1);

  auto system = _system.lock();
  unsigned int nAtoms = system->getAtoms().size();
  auto basisController = system->getAtomCenteredBasisController();
  const unsigned int nBasisFunctions = basisController->getNBasisFunctions();

  unsigned int derivorder = 0;
  if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
    derivorder = _basisFunctionOnGridController->getHighestDerivative();
    if (derivorder < 2)
      _basisFunctionOnGridController->setHighestDerivative(2);
  }

  // Create mapping
  auto orbitalSet = system->template getActiveOrbitalController<SCFMode>();
  std::vector<unsigned int> mapping = basisController->getAtomIndicesOfBasis();

  // Get potential and (in case of GGAs) potential gradients
  const auto& potential = *(funcData.dFdRho);
  auto potentialGradients = funcData.dFdGradRho;

  // Get grid and BFs on it
  unsigned int nBlocks = _basisFunctionOnGridController->getNBlocks();
  const Eigen::VectorXd& weights = _basisFunctionOnGridController->getGridController()->getWeights();

  // Get density Matrix
  auto densityMatrix(system->template getElectronicStructure<SCFMode>()->getDensityMatrix());
  Eigen::MatrixXd gradientContr = Eigen::MatrixXd::Zero(nAtoms, 3);

  Timings::takeTime("Active System - Func. Pot. Grad");
#pragma omp parallel
  {
    Eigen::MatrixXd gradientContrPriv(nAtoms, 3);
    gradientContrPriv.setZero();
    auto gpp = gradientContrPriv.data();
#pragma omp for schedule(dynamic)
    /*
     * Loop over Grid. Data is accessed blockwise.
     */
    for (unsigned int blockIndex = 0; blockIndex < nBlocks; ++blockIndex) {
      const auto& blockOnGridData = _basisFunctionOnGridController->getBlockOnGridData(blockIndex);
      const auto& bfValues = blockOnGridData->functionValues;
      assert(blockOnGridData->derivativeValues);
      const auto& bfDerivatives = *blockOnGridData->derivativeValues;
      const unsigned int blockSize = blockOnGridData->functionValues.rows();
      const unsigned int pointStartIndex = _basisFunctionOnGridController->getFirstIndexOfBlock(blockIndex);
      SpinPolarizedData<SCFMode, Eigen::VectorXd> weightedPot;
      for_spin(weightedPot, potential) {
        weightedPot_spin =
            weights.segment(pointStartIndex, blockSize).cwiseProduct(potential_spin.segment(pointStartIndex, blockSize));
      };

      for (unsigned int mu = 0; mu < nBasisFunctions; ++mu) {
        if (blockOnGridData->negligible[mu])
          continue;
        SpinPolarizedData<SCFMode, Eigen::VectorXd> potMu;
        for_spin(potMu, weightedPot) {
          potMu_spin = weightedPot_spin.cwiseProduct(bfValues.col(mu));
        };
        // auto muAtom = mapping[mu];
        for (unsigned int nu = 0; nu < nBasisFunctions; ++nu) {
          if (blockOnGridData->negligible[nu])
            continue;
          auto nuAtom = mapping[nu];
          // A common prefactor for the calculations below
          for_spin(potMu, densityMatrix) {
            // factor of 2 due to the product rule
            double prefac = 2 * densityMatrix_spin(mu, nu);
            gpp[nuAtom + (0 * nAtoms)] -= prefac * potMu_spin.dot(bfDerivatives.x.col(nu));
            gpp[nuAtom + (1 * nAtoms)] -= prefac * potMu_spin.dot(bfDerivatives.y.col(nu));
            gpp[nuAtom + (2 * nAtoms)] -= prefac * potMu_spin.dot(bfDerivatives.z.col(nu));
          };
          if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
            const auto& bfSecondDer = *blockOnGridData->secondDerivativeValues;

            // variables for the spin macro
            const auto& vdx = potentialGradients->x;
            const auto& vdy = potentialGradients->y;
            const auto& vdz = potentialGradients->z;
            for_spin(vdx, vdy, vdz, densityMatrix) {
              Eigen::VectorXd prefac = 2 * densityMatrix_spin(mu, nu) * weights.segment(pointStartIndex, blockSize);

              Eigen::VectorXd potentialx = vdx_spin.segment(pointStartIndex, blockSize).array() * prefac.array();
              Eigen::VectorXd potentialy = vdy_spin.segment(pointStartIndex, blockSize).array() * prefac.array();
              Eigen::VectorXd potentialz = vdz_spin.segment(pointStartIndex, blockSize).array() * prefac.array();
              Eigen::VectorXd XmuXnu = bfDerivatives.x.col(mu).array() * bfDerivatives.x.col(nu).array();
              Eigen::VectorXd XmuYnu = bfDerivatives.x.col(mu).array() * bfDerivatives.y.col(nu).array();
              Eigen::VectorXd XmuZnu = bfDerivatives.x.col(mu).array() * bfDerivatives.z.col(nu).array();
              Eigen::VectorXd YmuXnu = bfDerivatives.y.col(mu).array() * bfDerivatives.x.col(nu).array();
              Eigen::VectorXd YmuYnu = bfDerivatives.y.col(mu).array() * bfDerivatives.y.col(nu).array();
              Eigen::VectorXd YmuZnu = bfDerivatives.y.col(mu).array() * bfDerivatives.z.col(nu).array();
              Eigen::VectorXd ZmuXnu = bfDerivatives.z.col(mu).array() * bfDerivatives.x.col(nu).array();
              Eigen::VectorXd ZmuYnu = bfDerivatives.z.col(mu).array() * bfDerivatives.y.col(nu).array();
              Eigen::VectorXd ZmuZnu = bfDerivatives.z.col(mu).array() * bfDerivatives.z.col(nu).array();

              gpp[nuAtom + (0 * nAtoms)] -=
                  (potentialx.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xx.col(nu).array())) + XmuXnu)) +
                  (potentialy.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xy.col(nu).array())) + YmuXnu)) +
                  (potentialz.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xz.col(nu).array())) + ZmuXnu));

              gpp[nuAtom + (1 * nAtoms)] -=
                  (potentialx.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xy.col(nu).array())) + XmuYnu)) +
                  (potentialy.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.yy.col(nu).array())) + YmuYnu)) +
                  (potentialz.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.yz.col(nu).array())) + ZmuYnu));

              gpp[nuAtom + (2 * nAtoms)] -=
                  (potentialx.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.xz.col(nu).array())) + XmuZnu)) +
                  (potentialy.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.yz.col(nu).array())) + YmuZnu)) +
                  (potentialz.dot(Eigen::VectorXd((bfValues.col(mu).array() * bfSecondDer.zz.col(nu).array())) + ZmuZnu));
            };
          }
        }
      }
    }
#pragma omp critical
    { gradientContr += gradientContrPriv; }
  } /* END OpenMP parallel */
  Timings::timeTaken("Active System - Func. Pot. Grad");
  if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
    _basisFunctionOnGridController->setHighestDerivative(derivorder);
  }
  return gradientContr;
}

template class FuncPotential<Options::SCF_MODES::RESTRICTED>;
template class FuncPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
