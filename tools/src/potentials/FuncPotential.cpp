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
  std::vector<unsigned int> mapping(nBasisFunctions);
  std::vector<bool> hasElementBeenSet(nBasisFunctions, false);
  const auto& basisIndices = basisController->getBasisIndices();
  // Vector to check whether ALL basis function shells are assigned to an atom index
  for (unsigned int iAtom = 0; iAtom < basisIndices.size(); ++iAtom) {
    const unsigned int firstIndex = basisIndices[iAtom].first;
    const unsigned int endIndex = basisIndices[iAtom].second;
    for (unsigned int mu = firstIndex; mu < endIndex; ++mu) {
      mapping[mu] = iAtom;
      hasElementBeenSet[mu] = true;
    }
  }
  // Check
  for (bool x : hasElementBeenSet) {
    if (not x)
      throw SerenityError("FuncPotential: Missed gradient element in gradient evaluation.");
  }

  // Get potential and (in case of GGAs) potential gradients
  const auto& potential = *(funcData.dFdRho);
  auto potentialGradients = funcData.dFdGradRho;

  // Get grid and BFs on it TODO check whether derivative level needs to be modified
  unsigned int nBlocks = _basisFunctionOnGridController->getNBlocks();
  const auto& weights = _basisFunctionOnGridController->getGridController()->getWeights();

  // Get density Matrix
  auto densityMatrix(system->template getElectronicStructure<SCFMode>()->getDensityMatrix());
  Eigen::MatrixXd gradientContr = Eigen::MatrixXd::Zero(nAtoms, 3);
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
      /*
       * Check ExchangeCorrelationDerivativeCalculator.h for formulas and reference
       */
      for (unsigned int mu = 0; mu < nBasisFunctions; ++mu) {
        auto muAtom = mapping[mu];
        if (blockOnGridData->negligible[mu])
          continue;
        for (unsigned int nu = 0; nu <= mu; ++nu) {
          if (blockOnGridData->negligible[nu])
            continue;
          auto nuAtom = mapping[nu];
          const unsigned int pointStartIndex = _basisFunctionOnGridController->getFirstIndexOfBlock(blockIndex);
          // A common prefactor for the calculations below
          for_spin(potential, densityMatrix) {
            Eigen::VectorXd prefac =
                (mu == nu ? 1.0 : 2.0) * densityMatrix_spin(mu, nu) * weights.segment(pointStartIndex, blockSize);
            Eigen::VectorXd pot = potential_spin.segment(pointStartIndex, blockSize).array() * prefac.array();
            gpp[muAtom + (0 * nAtoms)] -=
                pot.dot(Eigen::VectorXd(bfValues.col(nu).array() * bfDerivatives.x.col(mu).array()));
            gpp[muAtom + (1 * nAtoms)] -=
                pot.dot(Eigen::VectorXd(bfValues.col(nu).array() * bfDerivatives.y.col(mu).array()));
            gpp[muAtom + (2 * nAtoms)] -=
                pot.dot(Eigen::VectorXd(bfValues.col(nu).array() * bfDerivatives.z.col(mu).array()));
            gpp[nuAtom + (0 * nAtoms)] -=
                pot.dot(Eigen::VectorXd(bfValues.col(mu).array() * bfDerivatives.x.col(nu).array()));
            gpp[nuAtom + (1 * nAtoms)] -=
                pot.dot(Eigen::VectorXd(bfValues.col(mu).array() * bfDerivatives.y.col(nu).array()));
            gpp[nuAtom + (2 * nAtoms)] -=
                pot.dot(Eigen::VectorXd(bfValues.col(mu).array() * bfDerivatives.z.col(nu).array()));
          };
          if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
            const auto& bfSecondDer = *blockOnGridData->secondDerivativeValues;

            // variables for the spin macro
            const auto& vdx = potentialGradients->x;
            const auto& vdy = potentialGradients->y;
            const auto& vdz = potentialGradients->z;
            for_spin(vdx, vdy, vdz, densityMatrix) {
              Eigen::VectorXd prefac =
                  (mu == nu ? 1.0 : 2.0) * densityMatrix_spin(mu, nu) * weights.segment(pointStartIndex, blockSize);

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

              gpp[muAtom + (0 * nAtoms)] -=
                  (potentialx.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xx.col(mu).array())) + XmuXnu)) +
                  (potentialy.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xy.col(mu).array())) + XmuYnu)) +
                  (potentialz.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xz.col(mu).array())) + XmuZnu));

              gpp[muAtom + (1 * nAtoms)] -=
                  (potentialx.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xy.col(mu).array())) + YmuXnu)) +
                  (potentialy.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.yy.col(mu).array())) + YmuYnu)) +
                  (potentialz.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.yz.col(mu).array())) + YmuZnu));

              gpp[muAtom + (2 * nAtoms)] -=
                  (potentialx.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.xz.col(mu).array())) + ZmuXnu)) +
                  (potentialy.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.yz.col(mu).array())) + ZmuYnu)) +
                  (potentialz.dot(Eigen::VectorXd((bfValues.col(nu).array() * bfSecondDer.zz.col(mu).array())) + ZmuZnu));

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
  if (_functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
    _basisFunctionOnGridController->setHighestDerivative(derivorder);
  }
  return gradientContr;
}

template class FuncPotential<Options::SCF_MODES::RESTRICTED>;
template class FuncPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
