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
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/GridData.h"
#include "data/grid/MatrixOperatorToGridTransformer.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "dft/functionals/wrappers/XCFun.h"
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
    _orbitalController(systemController->getActiveOrbitalController<SCFMode>()),
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
    XCFun<SCFMode> xcFun(128);
    if (_orbitalController->orbitalsAvailable()) {
      // the density
      const auto& density = _densityOnGridController->getDensityOnGrid();
      Matrix<int> evaluatedGridPoints;
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
          unsigned int nBasisFunc(newDensityMatrix.cols());
          assert(nBasisFunc == currentCoefVector.size());
          for (unsigned int k = 0; k < nBasisFunc; ++k) {
            for (unsigned int l = 0; l < nBasisFunc; ++l) {
              newDensityMatrix(l, k) += _occupations * currentCoefVector(l) * currentCoefVector(k);
            }
          }
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
      auto funcDataVWN = xcFun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, _densityOnGridController);
      auto& nu_LDA_c = *funcDataVWN.dFdRho;
      functional = Functional(CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR,
                              {BasicFunctionals::BASIC_FUNCTIONALS::X_SLATER}, {1.0});
      auto funcDataSLATER = xcFun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, _densityOnGridController);
      auto& nu_LDA_x = *funcDataSLATER.dFdRho;
      /*
       * 1. e) Calculate BECKEX (Becke 88) and PW91C (Perdew-Wang 91) energy densities
       */
      functional = Functional(CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR,
                              {BasicFunctionals::BASIC_FUNCTIONALS::C_PW91}, {1.0});
      auto funcDataPerdew = xcFun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, _densityOnGridController);
      double energyPerdew = funcDataPerdew.energy;
      auto epuvPerdew = funcDataPerdew.epuv;
      functional = Functional(CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR,
                              {BasicFunctionals::BASIC_FUNCTIONALS::X_B86}, {1.0});
      auto funcDataBECKE = xcFun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, _densityOnGridController);
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
        unsigned int lastGridIndex = firstGridIndex + blocksize;
        /*
         * 2. a) Calculate all LB potential
         *       alpha nu_x,LDA+ nu_C,LDA - beta x(r)^2 rho^(1/3)(r) / (1+3 beta x(r) ln[(x(r)+sqrt{x(r)^2+1})]),
         *       the reciprocal of rho (needed for the coefficient q)
         *       and calculating nu_hole
         */
        // 1/rho(r) for the current block
        SpinPolarizedData<SCFMode, Eigen::VectorXd> recRho(blocksize);
        // the LB potential contribution
        SpinPolarizedData<SCFMode, Eigen::VectorXd> nu_lb(blocksize);
        // nu_hole
        Eigen::VectorXd nu_hole(blocksize);
        nu_hole.setZero();
        for_spin(nu_lb) {
          nu_lb_spin.setZero();
        };
        recRho = nu_lb;
        // calculating the LDA contribution to the LB potential separately to make use of the for_spin makro.
        for (unsigned int r = firstGridIndex; r < lastGridIndex; ++r) {
          unsigned int reducedIndex = r - firstGridIndex;
          for_spin(nu_LDA_x, nu_LDA_c) {
            nu_LDA_x_spin[r] = _alpha * nu_LDA_x_spin[r] + nu_LDA_c_spin[r];
          };
          /*
           * Storing the density gradient differently to make further calculations more convenient
           */
          SpinPolarizedData<SCFMode, Eigen::Vector3d> dRho_dr;
          auto& drho_drx = drho_dr.x;
          auto& drho_dry = drho_dr.y;
          auto& drho_drz = drho_dr.z;
          for_spin(drho_drx, drho_dry, drho_drz, dRho_dr) {
            dRho_dr_spin = Eigen::Vector3d(drho_drx_spin(r), drho_dry_spin(r), drho_drz_spin(r));
          };
          for_spin(dRho_dr, density, recRho, nu_lb, nu_LDA_x) {
            if (density_spin(r) > 1E-9) {
              recRho_spin(reducedIndex) = 1 / density_spin(r);
              // lb potential
              if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
                auto x_r = dRho_dr_spin.norm() / (pow(density_spin(r), 4.0 / 3.0));
                nu_lb_spin(reducedIndex) = nu_LDA_x_spin[r] - (_beta * x_r * x_r * pow(density_spin(r), 1.0 / 3.0)) /
                                                                  (1 + 3 * _beta * x_r * log(x_r + sqrt(x_r * x_r + 1)));
              }
              else {
                auto x_r = dRho_dr_spin.norm() / (2 * pow(density_spin(r) / 2, 4.0 / 3.0));
                nu_lb_spin(reducedIndex) = nu_LDA_x_spin[r] - (_beta * x_r * x_r * pow(density_spin(r) / 2, 1.0 / 3.0)) /
                                                                  (1 + 3 * _beta * x_r * log(x_r + sqrt(x_r * x_r + 1)));
              }
            }
          };
          // nu_hole
          double totalDensity_r = 0.0;
          for_spin(density) {
            totalDensity_r += density_spin(r);
          };
          if (totalDensity_r > 1E-9) {
            nu_hole(reducedIndex) = (2 * (*epuvPerdew)[r] + 2 * (*epuvBecke)[r]) / totalDensity_r;
          }
        } /* for r */
        /*
         * 2. b) Calculate all scaling coefficients (psi_i^2/rho) in matrix representation -> nOCC x blocksize
         */
        // the qout(r,i) matrix
        SpinPolarizedData<SCFMode, Matrix<double>> qout;
        for_spin(recRho, densityMatrices, qout, _nOCC) {
          qout_spin.resize(blocksize, _nOCC_spin);
          qout_spin.setZero();
          for (unsigned int i = 1; i < _nOCC_spin; ++i) {
            for (unsigned int col = 0; col < basisFuncBlockVal.cols(); col++) {
              if (negligible[col])
                continue;
              Eigen::VectorXd tmp = Eigen::VectorXd::Zero(blocksize);
              for (unsigned int row = 0; row < basisFuncBlockVal.cols(); row++) {
                if (negligible[row])
                  continue;
                tmp += basisFuncBlockVal.col(row) * densityMatrices_spin[i](row, col);
              }
              qout_spin.col(i).array() += tmp.array() * basisFuncBlockVal.col(col).array();
            }
            qout_spin.col(i).array() = qout_spin.col(i).array() * recRho_spin.array();
          }
        };
        /*
         * 3. Adding up
         */
        for_spin(nu_lb, exp_eps, sqrt_eps, qout, saopGridPotential) {
          Eigen::VectorXd blockPot(blocksize);
          Eigen::VectorXd E = Eigen::VectorXd::Constant(exp_eps_spin.size(), 1);
          //            std::cout << (E-exp_eps_spin)+exp_eps_spin << std::endl;
          blockPot = ((nu_lb_spin * exp_eps_spin.transpose()) +
                      ((qout_spin * sqrt_eps_spin) + nu_hole) * (E - exp_eps_spin).transpose())
                         .cwiseProduct(qout_spin)
                         .rowwise()
                         .sum();
          for (unsigned int r = firstGridIndex; r < lastGridIndex; ++r) {
            unsigned int reducedIndex = r - firstGridIndex;
            saopGridPotential_spin[r] = blockPot(reducedIndex);
          } /* for r */
        };  /* for spin */
      }     /* for iBlock, parallel for */
      // sanity check if all grid points have been evaluated
      assert(evaluatedGridPoints.sum() == static_cast<int>(_grid->getNGridPoints()) &&
             "not all grid points were evaluated! Something is wrong here!");
      /*
       * 4. Transforming the potential from grid representation to a matrix
       */
      _gridToMatrix->addScalarOperatorToMatrix(pot, saopGridPotential);
    } /* if orbtialsAvailable */
    else {
      std::cout << std::endl;
      std::cout << "SAOP-Potential: No initial orbital data available."
                << " Making initial guess using LDA functional." << std::endl;
      std::cout << std::endl;
      auto funcData = xcFun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS,
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
