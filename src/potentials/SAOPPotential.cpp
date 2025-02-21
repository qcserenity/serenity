/**
 * @file   SAOPPotential.cpp
 *
 * @date   Aug 4, 2017
 * @author Moritz Bensberg
 *
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
#include "potentials/SAOPPotential.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridController.h"
#include "data/grid/GridPotential.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/DensityMatrixController.h"
#include "dft/Functional.h"
#include "dft/functionals/BasicFunctionals.h"
#include "dft/functionals/FunctionalLibrary.h"
#include "io/FormattedOutputStream.h"
#include "misc/HelperFunctions.h"
#include "settings/DFTOptions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SAOPPotential<SCFMode>::SAOPPotential(std::shared_ptr<BasisController> basis, std::shared_ptr<SystemController> systemController,
                                      std::shared_ptr<DensityOnGridController<SCFMode>> densityOnGridController)
  : Potential<SCFMode>(basis),
    _densityOnGridController(densityOnGridController),
    _nOCC(systemController->getNOccupiedOrbitals<SCFMode>()),
    _systemController(systemController),
    _potential(nullptr),
    _energy(0.0),
    _grid(densityOnGridController->getGridController()) {
  _occupations = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2 : 1;
  _basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(systemController->getSettings(), this->_basis, this->_grid);
  _gridToMatrix = std::make_shared<ScalarOperatorToMatrixAdder<SCFMode>>(
      _basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& SAOPPotential<SCFMode>::getMatrix() {
  if (!_potential) {
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot) {
      pot_spin.setZero();
    };
    if (_systemController->template hasElectronicStructure<SCFMode>())
      _orbitalController = _systemController->template getActiveOrbitalController<SCFMode>();
    FunctionalLibrary<SCFMode> flib(128);
    if (_orbitalController) {
      // the density
      const auto& density = _densityOnGridController->getDensityOnGrid();
      Eigen::MatrixXi evaluatedGridPoints;
      // the density gradient
      auto& drho_dr = _densityOnGridController->getDensityGradientOnGrid();
      // the potential
      GridPotential<SCFMode> saopGridPotential(_densityOnGridController->getGridController());
      /*
       * 1.a) Calculation of the density matrices P_i
       *   b) Calculate all exp(-2(eps_n-eps_i)
       *   c) Calculate all sqrt(eps_n-eps_i)
       *   d) Calculate LDA potential
       *   e) Calculate BECKEX (Becke 88) and PW91C (Perdew-Wang 91) energy densities
       * 2. Loop over all blocks of the size "blocksize"
       *   Preparing r-specific data
       *   a) Calculate all LB potential fracs -beta x(r)^2 rho^(1/3)(r) / (1+3 beta x(r) ln[(x(r)+sqrt{x(r)^2+1})])
       *   b) Calculate all scaling coefficients q = (psi_i^2/rho) in matrix representation -> blocksize x nOCC
       * 3. Adding up
       */

      /*
       * 1. a) Calculation of the density matrices P_i
       */
      // getting the active orbital set
      const auto& coefficients = _orbitalController->getCoefficients();
      // building a vector of density matrices. Each matrix corresponds to one occupied orbital.
      SpinPolarizedData<SCFMode, std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>> densityMatrices;
      for_spin(_nOCC, coefficients, densityMatrices) {
        MatrixInBasis<Options::SCF_MODES::RESTRICTED> newDensityMatrix(_orbitalController->getBasisController());
        for (unsigned int i = 0; i < _nOCC_spin; ++i) {
          newDensityMatrix.setZero();
          auto currentCoefVector = coefficients_spin.col(i);
          newDensityMatrix += _occupations * currentCoefVector * currentCoefVector.transpose();
          densityMatrices_spin.push_back(newDensityMatrix);
        }
      };
      /*
       * 1. b) Calculate all exp(-2(eps_n-eps_i)^2 and 1 - exp(-2(eps_n-eps_i)^2
       *    c) Calculate all sqrt(eps_n-eps_i)
       */
      // getting the orbital eigenvalues
      const auto& eigenvalues(_orbitalController->getEigenvalues());
      SpinPolarizedData<SCFMode, Eigen::VectorXd> exp_eps(1);
      SpinPolarizedData<SCFMode, Eigen::VectorXd> sqrt_eps(1);
      for_spin(exp_eps, sqrt_eps, eigenvalues, _nOCC) {
        exp_eps_spin.resize(_nOCC_spin);
        exp_eps_spin.setZero();
        sqrt_eps_spin = exp_eps_spin;
        const double homoEigenvalue(eigenvalues_spin[_nOCC_spin - 1]);
        for (unsigned int i = 0; i < _nOCC_spin; ++i) {
          double deltaEps_i = homoEigenvalue - eigenvalues_spin[i];
          exp_eps_spin(i) = exp(-2 * deltaEps_i * deltaEps_i);
          if (deltaEps_i > 1E-6) {
            sqrt_eps_spin(i) = _K * sqrt(deltaEps_i);
          }
        } /* for i */
      };  /* for spin */
      /*
       * 1. d) Calculate LDA potential part of the LB potential
       */
      Functional functional(CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR,
                            {BasicFunctionals::BASIC_FUNCTIONALS::C_VWN}, {1.0});
      auto funcDataVWN = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, _densityOnGridController);
      auto& nu_LDA_c = *funcDataVWN.dFdRho;
      functional = Functional(CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR,
                              {BasicFunctionals::BASIC_FUNCTIONALS::X_SLATER}, {1.0});
      auto funcDataSLATER = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, _densityOnGridController);
      auto& nu_LDA_x = *funcDataSLATER.dFdRho;
      /*
       * 1. e) Calculate BECKEX (Becke 88) and PW91C (Perdew-Wang 91) energy densities
       */
      functional = Functional(CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR,
                              {BasicFunctionals::BASIC_FUNCTIONALS::C_PW91}, {1.0});
      auto funcDataPerdew = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, _densityOnGridController);
      double energyPerdew = funcDataPerdew.energy;
      auto epuvPerdew = funcDataPerdew.epuv;
      functional = Functional(CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR,
                              {BasicFunctionals::BASIC_FUNCTIONALS::X_B88}, {1.0});
      auto funcDataBECKE = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, _densityOnGridController);
      double energyBecke = funcDataBECKE.energy;
      auto epuvBecke = funcDataBECKE.epuv;
      // using the energy provided by the Becke exachange and Perdew correlation functionals
      _energy = energyPerdew + energyBecke;

      /*
       * 2. Loop over all blocks of the size "blocksize"
       */
      // number of blocks
      unsigned int nBlocks = _basisFunctionOnGridController->getNBlocks();
      evaluatedGridPoints.resize(nBlocks, 1);
      evaluatedGridPoints.setZero();
#pragma omp parallel for schedule(static, 1)
      for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
        // getting the basis function values of this block
        auto blockData = _basisFunctionOnGridController->getBlockOnGridData(iBlock);
        const auto& basisFuncBlockVal = blockData->functionValues; // blocksize x nBasisFunc
        const auto& negligible = blockData->negligible;
        // the block size
        unsigned int blocksize = basisFuncBlockVal.rows();
        evaluatedGridPoints(iBlock) = blocksize;
        // the first and last grid index of the current block
        unsigned int firstGridIndex = _basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
        /*
         * 2. a) Calculate all LB potential
         *       alpha nu_x,LDA+ nu_C,LDA - beta x(r)^2 rho^(1/3)(r) / (1+3 beta x(r) ln[(x(r)+sqrt{x(r)^2+1})]),
         *       the reciprocal of rho (needed for the coefficient q)
         *       and calculating nu_hole
         */
        // 1/rho(r) for the current block
        SpinPolarizedData<SCFMode, Eigen::VectorXd> recRho(Eigen::VectorXd::Zero(blocksize));
        Eigen::VectorXd recTotalRho = Eigen::VectorXd::Zero(blocksize);
        for (unsigned int r = 0; r < blocksize; ++r) {
          double totalDens = 0.0;
          for_spin(density, recRho) {
            double rho_r_spin = density_spin(r + firstGridIndex);
            totalDens += rho_r_spin;
            if (rho_r_spin > 1e-9)
              recRho_spin(r) = 1.0 / rho_r_spin;
          };
          if (totalDens > 1E-9)
            recTotalRho(r) = 1.0 / totalDens;
        }
        // the LB potential contribution
        SpinPolarizedData<SCFMode, Eigen::VectorXd> nu_lb(Eigen::VectorXd::Zero(blocksize));
        // nu_hole
        Eigen::VectorXd nu_hole = Eigen::VectorXd::Zero(blocksize);
        // calculating the LDA contribution to the LB potential separately to make use of the for_spin macro.
        for_spin(nu_LDA_x, nu_LDA_c) {
          nu_LDA_x_spin.segment(firstGridIndex, blocksize) = _alpha * nu_LDA_x_spin.segment(firstGridIndex, blocksize) +
                                                             nu_LDA_c_spin.segment(firstGridIndex, blocksize);
        };
        auto& drho_drx = drho_dr.x;
        auto& drho_dry = drho_dr.y;
        auto& drho_drz = drho_dr.z;
        for_spin(drho_drx, drho_dry, drho_drz, density, recRho, nu_lb, nu_LDA_x) {
          const Eigen::VectorXd nablaRhoNorm = (drho_drx_spin.segment(firstGridIndex, blocksize).array() *
                                                    drho_drx_spin.segment(firstGridIndex, blocksize).array() +
                                                drho_dry_spin.segment(firstGridIndex, blocksize).array() *
                                                    drho_dry_spin.segment(firstGridIndex, blocksize).array() +
                                                drho_drz_spin.segment(firstGridIndex, blocksize).array() *
                                                    drho_drz_spin.segment(firstGridIndex, blocksize).array())
                                                   .sqrt();
          const auto one = Eigen::ArrayXd::Constant(blocksize, 1.0);
          if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
            const Eigen::ArrayXd x = nablaRhoNorm.array() * recRho_spin.array().pow(4.0 / 3.0);
            const Eigen::ArrayXd x2 = x * x;
            nu_lb_spin = nu_LDA_x_spin.segment(firstGridIndex, blocksize).array() -
                         (_beta * x2 * density_spin.segment(firstGridIndex, blocksize).array().pow(1.0 / 3.0)) /
                             (one + 3.0 * _beta * x * (x + (x2 + one).sqrt()).log());
          }
          else {
            const Eigen::ArrayXd x = 0.5 * nablaRhoNorm.array() * (2 * recRho_spin).array().pow(4.0 / 3.0);
            const Eigen::ArrayXd x2 = x * x;
            nu_lb_spin = nu_LDA_x_spin.segment(firstGridIndex, blocksize).array() -
                         (_beta * x2 * (0.5 * density_spin.segment(firstGridIndex, blocksize)).array().pow(1.0 / 3.0)) /
                             (one + 3.0 * _beta * x * (x + (x2 + one).sqrt()).log());
          }
        };
        // nu_hole
        nu_hole = (2 * epuvPerdew->segment(firstGridIndex, blocksize).array() +
                   2 * epuvBecke->segment(firstGridIndex, blocksize).array()) *
                  recTotalRho.array();
        /*
         * 2. b) Calculate all scaling coefficients (psi_i^2/rho) in matrix representation -> nOCC x blocksize
         */
        // the qout(r,i) matrix
        SpinPolarizedData<SCFMode, Eigen::MatrixXd> qout;
        const Eigen::MatrixXd pA = constructProjectionMatrix(negligible);
        const Eigen::MatrixXd basisFuncVA = basisFuncBlockVal * pA;
        for_spin(recRho, densityMatrices, qout, _nOCC) {
          qout_spin.resize(blocksize, _nOCC_spin);
          qout_spin.setZero();
          for (unsigned int i = 1; i < _nOCC_spin; ++i) {
            const Eigen::MatrixXd signDensMat_i = pA.transpose() * densityMatrices_spin[i] * pA;
            const unsigned int nSign = signDensMat_i.cols();
            for (unsigned int col = 0; col < nSign; col++) {
              Eigen::VectorXd tmp = Eigen::VectorXd::Zero(blocksize);
              for (unsigned int row = 0; row < nSign; row++) {
                tmp += basisFuncVA.col(row) * signDensMat_i(row, col);
              }
              qout_spin.col(i).array() += tmp.array() * basisFuncVA.col(col).array();
            }
            qout_spin.col(i).array() = qout_spin.col(i).array() * recRho_spin.array();
          }
        };
        /*
         * 3. Adding up
         */
        for_spin(nu_lb, exp_eps, sqrt_eps, qout, saopGridPotential) {
          const Eigen::VectorXd one = Eigen::VectorXd::Constant(exp_eps_spin.size(), 1);
          saopGridPotential_spin.segment(firstGridIndex, blocksize) =
              ((nu_lb_spin * exp_eps_spin.transpose()) +
               ((qout_spin * sqrt_eps_spin) + nu_hole) * (one - exp_eps_spin).transpose())
                  .cwiseProduct(qout_spin)
                  .rowwise()
                  .sum();
        }; /* for spin */
      }    /* for iBlock, parallel for */
      /*
       * 4. Transforming the potential from grid representation to a matrix
       */
      _gridToMatrix->addScalarOperatorToMatrix(pot, saopGridPotential);
    } /* if orbitalsAvailable */
    else {
      OutputControl::vOut << std::endl;
      OutputControl::vOut << "SAOP-Potential: No initial orbital data available."
                          << " Making initial guess using LDA functional." << std::endl
                          << std::endl;
      auto funcData = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS,
                                    CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA),
                                    _densityOnGridController);
      _gridToMatrix->addScalarOperatorToMatrix(pot, *funcData.dFdRho);
      _energy = funcData.energy;
    } /* else if orbitalsAvailable */
  }   /* if !_potential */
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
double SAOPPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& /*P*/) {
  if (!_potential)
    getMatrix();
  return _energy;
};

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd SAOPPotential<SCFMode>::getGeomGradients() {
  assert(false && "No geometry gradients availabe for the SAOP model potential.");
  Eigen::MatrixXd gradientContr(1, 3);
  return gradientContr;
}

template class SAOPPotential<Options::SCF_MODES::RESTRICTED>;
template class SAOPPotential<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
