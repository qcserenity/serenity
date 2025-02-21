/**
 * @file   XCFun.cpp
 *
 * @date   Feb 3, 2017
 * @author M. Boeckers
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
#ifdef SERENITY_USE_XCFUN
/* Include Class Header*/
#include "dft/functionals/wrappers/XCFun.h"
/* Include Serenity Internal Headers */
#include "dft/functionals/BasicFunctionals.h"
#include "grid/GridController.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
XCFun<SCFMode>::XCFun(unsigned int maxBlockSize) : _maxBlockSize(maxBlockSize) {
}

template<Options::SCF_MODES SCFMode>
FunctionalData<SCFMode> XCFun<SCFMode>::calcData(FUNCTIONAL_DATA_TYPE type, const Functional functional,
                                                 const std::shared_ptr<DensityOnGridController<SCFMode>> densityOnGridController,
                                                 unsigned int order) {
  // Check input
  if (order > 3)
    throw SerenityError("XCFun usage is only possible up to 3rd order derivatives w.r.t. electrons/sigma.");
  // Build functional from basic functionals
  auto func = getFunctional(functional);
  // GGA?
  const bool gga = (xcfun_is_gga(func) != 0) ? true : false;
  // Get density
  auto& density = densityOnGridController->getDensityOnGrid();
  // If gga, get gradient of density
  std::shared_ptr<Gradient<DensityOnGrid<SCFMode>>> gradient = nullptr;
  // If gga and potential is requested, get also hessian
  std::shared_ptr<Hessian<DensityOnGrid<SCFMode>>> hessian = nullptr;
  if (gga) {
    if (densityOnGridController->getHighestDerivative() < 1)
      throw SerenityError(
          "XCFun: Density gradient required for GGA functional, but DensityOnGridController does not have it.");
    gradient = std::make_shared<Gradient<DensityOnGrid<SCFMode>>>(densityOnGridController->getDensityGradientOnGrid());
    if (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
      if (densityOnGridController->getHighestDerivative() < 2)
        throw SerenityError(
            "XCFun: Density Hessian required for GGA potential, but DensityOnGridController does not have it.");
      hessian = std::make_shared<Hessian<DensityOnGrid<SCFMode>>>(densityOnGridController->getDensityHessianOnGrid());
    }
  }
  Timings::takeTime("Tech. - XCFun Functional Eval.");
  // Setup XCFun
  int iErr;
  xcfun_vars xcVars = xcfun_vars::XC_N;
  xcfun_mode xcMode = (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) ? XC_POTENTIAL : XC_PARTIAL_DERIVATIVES;
  if (!gga) {
    // If LDA:
    if (SCFMode == Options::SCF_MODES::RESTRICTED) {
      xcVars = xcfun_vars::XC_N;
    }
    else {
      xcVars = xcfun_vars::XC_A_B;
    }
    iErr = xcfun_eval_setup(func, xcVars, XC_PARTIAL_DERIVATIVES, order);
  }
  else {
    // If GGA:
    if (SCFMode == Options::SCF_MODES::RESTRICTED) {
      // If restricted
      if (type == FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS) {
        xcVars = xcfun_vars::XC_N_GNN;
      }
      else if (type == FUNCTIONAL_DATA_TYPE::GRADIENTS) {
        xcVars = xcfun_vars::XC_N_NX_NY_NZ;
      }
      else if (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
        order = 1;
        xcVars = xcfun_vars::XC_N_2ND_TAYLOR;
      }
      else
        throw SerenityError("XCFun (restricted, GGA): Problem with functional data type.");
    }
    else {
      // If unrestricted
      if (type == FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS) {
        xcVars = xcfun_vars::XC_A_B_GAA_GAB_GBB;
      }
      else if (type == FUNCTIONAL_DATA_TYPE::GRADIENTS) {
        xcVars = xcfun_vars::XC_A_B_AX_AY_AZ_BX_BY_BZ;
      }
      else if (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
        order = 1;
        xcVars = xcfun_vars::XC_A_B_2ND_TAYLOR;
      }
      else
        throw SerenityError("XCFun (unrestricted, GGA): Problem with functional data type.");
    }
    iErr = xcfun_eval_setup(func, xcVars, xcMode, order);
  }
  if (iErr != 0)
    throw SerenityError("Failed to set vars, mode and order for XCFun.");
  // In and output dimension
  int nOut = xcfun_output_length(func);
  int nIn = xcfun_input_length(func);

  // Number of grid points and blocks
  const unsigned int nPoints = density.getNGridPoints();
  const unsigned int nBlocks = (unsigned int)ceil((double)nPoints / _maxBlockSize);
  // Build FunctionaData object
  FunctionalData<SCFMode> funcData(order, type, functional, densityOnGridController->getGridController());

  // Loop over blocks of grid points
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    // first index of block
    const unsigned int firstIndex = iBlock * _maxBlockSize;
    // size of this block
    const unsigned int blockSize = determineBlockSize(iBlock, nPoints, nBlocks);
    bool skip = true;
    for_spin(density) {
      skip = skip && (density_spin.segment(firstIndex, blockSize).array().abs().sum() < blockSize * 1e-12);
    };
    if (skip)
      continue;
    // Input for this block
    Eigen::MatrixXd input(nIn, blockSize);
    prepareInput(firstIndex, xcVars, density, gradient, hessian, input);
    // Output object for this block
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(nOut, blockSize);
    double* ptr = input.data();
    double* outpt = output.data();
    for (unsigned int p = 0; p < blockSize; p++) {
      xcfun_eval(func, ptr + p * nIn, outpt + p * nOut);
    }
    // parse data for this block
    parseOutput(order, xcVars, firstIndex, output, funcData);
  } /* Loop over blocks */
  xcfun_delete(func);
  funcData.energy = calcEnergy(funcData.epuv, densityOnGridController->getGridController()->getWeights());
  Timings::timeTaken("Tech. - XCFun Functional Eval.");
  return funcData;
}

template<Options::SCF_MODES SCFMode>
unsigned int XCFun<SCFMode>::determineBlockSize(unsigned int blockIndex, unsigned int nPoints, unsigned int nBlocks) {
  unsigned int blockSize;
  if (blockIndex == nBlocks - 1) {
    blockSize = nPoints % _maxBlockSize;
    if (blockSize == 0)
      blockSize = _maxBlockSize;
  }
  else {
    blockSize = _maxBlockSize;
  }
  return blockSize;
}

template<>
Eigen::MatrixXd
XCFun<Options::SCF_MODES::RESTRICTED>::calculateSigma(const Gradient<DensityOnGrid<Options::SCF_MODES::RESTRICTED>>& gradient,
                                                      const unsigned int& iGridStart, const unsigned int& blocksize) {
  Eigen::MatrixXd sigma(1, blocksize);
  sigma.row(0).array() =
      gradient.x.segment(iGridStart, blocksize).array() * gradient.x.segment(iGridStart, blocksize).array() +
      gradient.y.segment(iGridStart, blocksize).array() * gradient.y.segment(iGridStart, blocksize).array() +
      gradient.z.segment(iGridStart, blocksize).array() * gradient.z.segment(iGridStart, blocksize).array();
  return sigma;
}

template<>
Eigen::MatrixXd
XCFun<Options::SCF_MODES::UNRESTRICTED>::calculateSigma(const Gradient<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>>& gradient,
                                                        const unsigned int& iGridStart, const unsigned int& blocksize) {
  Eigen::MatrixXd sigma(3, blocksize);
  sigma.row(0).array() =
      gradient.x.alpha.segment(iGridStart, blocksize).array() * gradient.x.alpha.segment(iGridStart, blocksize).array() +
      gradient.y.alpha.segment(iGridStart, blocksize).array() * gradient.y.alpha.segment(iGridStart, blocksize).array() +
      gradient.z.alpha.segment(iGridStart, blocksize).array() * gradient.z.alpha.segment(iGridStart, blocksize).array();
  sigma.row(1).array() =
      gradient.x.alpha.segment(iGridStart, blocksize).array() * gradient.x.beta.segment(iGridStart, blocksize).array() +
      gradient.y.alpha.segment(iGridStart, blocksize).array() * gradient.y.beta.segment(iGridStart, blocksize).array() +
      gradient.z.alpha.segment(iGridStart, blocksize).array() * gradient.z.beta.segment(iGridStart, blocksize).array();
  sigma.row(2).array() =
      gradient.x.beta.segment(iGridStart, blocksize).array() * gradient.x.beta.segment(iGridStart, blocksize).array() +
      gradient.y.beta.segment(iGridStart, blocksize).array() * gradient.y.beta.segment(iGridStart, blocksize).array() +
      gradient.z.beta.segment(iGridStart, blocksize).array() * gradient.z.beta.segment(iGridStart, blocksize).array();
  return sigma;
}

template<>
void XCFun<Options::SCF_MODES::RESTRICTED>::prepareInput(
    const unsigned int iGridStart, const xcfun_vars xcVars, const DensityOnGrid<Options::SCF_MODES::RESTRICTED>& density,
    const std::shared_ptr<Gradient<DensityOnGrid<Options::SCF_MODES::RESTRICTED>>> gradient,
    const std::shared_ptr<Hessian<DensityOnGrid<Options::SCF_MODES::RESTRICTED>>> hessian, Eigen::MatrixXd& input) {
  unsigned int n = input.cols();
  if (xcVars == xcfun_vars::XC_N) {
    input.row(0) = density.segment(iGridStart, n);
  }
  else if (xcVars == xcfun_vars::XC_N_GNN) {
    input.row(0) = density.segment(iGridStart, n);
    input.row(1) = calculateSigma(*gradient, iGridStart, n);
  }
  else if (xcVars == xcfun_vars::XC_N_2ND_TAYLOR) {
    input.row(0) = density.segment(iGridStart, n);
    input.row(1) = gradient->x.segment(iGridStart, n);
    input.row(2) = gradient->y.segment(iGridStart, n);
    input.row(3) = gradient->z.segment(iGridStart, n);
    input.row(4) = hessian->xx.segment(iGridStart, n);
    input.row(5) = hessian->xy.segment(iGridStart, n);
    input.row(6) = hessian->xz.segment(iGridStart, n);
    input.row(7) = hessian->yy.segment(iGridStart, n);
    input.row(8) = hessian->yz.segment(iGridStart, n);
    input.row(9) = hessian->zz.segment(iGridStart, n);
  }
  else if (xcVars == xcfun_vars::XC_N_NX_NY_NZ) {
    input.row(0) = density.segment(iGridStart, n);
    input.row(1) = gradient->x.segment(iGridStart, n);
    input.row(2) = gradient->y.segment(iGridStart, n);
    input.row(3) = gradient->z.segment(iGridStart, n);
  }
  else {
    throw SerenityError("XCFun (restricted): Input for requested xcfun_vars not implemented.");
  }
}

template<>
void XCFun<Options::SCF_MODES::RESTRICTED>::parseOutput(const unsigned int order, const xcfun_vars xcVars,
                                                        const unsigned int iGridStart, Eigen::MatrixXd& output,
                                                        FunctionalData<Options::SCF_MODES::RESTRICTED>& funcData) {
  unsigned int n = output.cols();
  // here xcVars can be xcfun_vars::XC_N_2ND_TAYLOR
  if (funcData.getType() == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
    (*funcData.epuv).segment(iGridStart, n) = output.row(0);
    (*funcData.potential).segment(iGridStart, n) = output.row(1);
  }
  else {
    if (xcVars == xcfun_vars::XC_N) {
      (*funcData.epuv).segment(iGridStart, n) = output.row(0);
      if (order >= 1)
        (*funcData.dFdRho).segment(iGridStart, n) = output.row(1);
      if (order >= 2)
        (*funcData.d2FdRho2).segment(iGridStart, n) = output.row(2);
      if (order >= 3)
        (*funcData.d3FdRho3).segment(iGridStart, n) = output.row(3);
    }
    else if (xcVars == xcfun_vars::XC_N_GNN) {
      (*funcData.epuv).segment(iGridStart, n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).segment(iGridStart, n) = output.row(1);
        (*funcData.dFdSigma).segment(iGridStart, n) = output.row(2);
      }
      if (order >= 2) {
        (*funcData.d2FdRho2).segment(iGridStart, n) = output.row(3);
        (*funcData.d2FdRhodSigma).segment(iGridStart, n) = output.row(4);
        (*funcData.d2FdSigma2).segment(iGridStart, n) = output.row(5);
      }
      if (order >= 3) {
        (*funcData.d3FdRho3).segment(iGridStart, n) = output.row(6);
        (*funcData.d3FdRho2dSigma).segment(iGridStart, n) = output.row(7);
        (*funcData.d3FdRhodSigma2).segment(iGridStart, n) = output.row(8);
        (*funcData.d3FdSigma3).segment(iGridStart, n) = output.row(9);
      }
    }
    else if (xcVars == xcfun_vars::XC_N_NX_NY_NZ) {
      (*funcData.epuv).segment(iGridStart, n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).segment(iGridStart, n) = output.row(1);
        (*funcData.dFdGradRho).x.segment(iGridStart, n) = output.row(2);
        (*funcData.dFdGradRho).y.segment(iGridStart, n) = output.row(3);
        (*funcData.dFdGradRho).z.segment(iGridStart, n) = output.row(4);
      }
      if (order >= 2) {
        (*funcData.d2FdRho2).segment(iGridStart, n) = output.row(5);
        (*funcData.d2FdRhodGradRho).x.segment(iGridStart, n) = output.row(6);
        (*funcData.d2FdRhodGradRho).y.segment(iGridStart, n) = output.row(7);
        (*funcData.d2FdRhodGradRho).z.segment(iGridStart, n) = output.row(8);
        (*funcData.d2FdGradRho2).xx.segment(iGridStart, n) = output.row(9);
        (*funcData.d2FdGradRho2).xy.segment(iGridStart, n) = output.row(10);
        (*funcData.d2FdGradRho2).xz.segment(iGridStart, n) = output.row(11);
        (*funcData.d2FdGradRho2).yy.segment(iGridStart, n) = output.row(12);
        (*funcData.d2FdGradRho2).yz.segment(iGridStart, n) = output.row(13);
        (*funcData.d2FdGradRho2).zz.segment(iGridStart, n) = output.row(14);
      }
      if (order >= 3) {
        (*funcData.d3FdRho3).segment(iGridStart, n) = output.row(15);
        (*funcData.d3FdRho2dGradRho).x.segment(iGridStart, n) = output.row(16);
        (*funcData.d3FdRho2dGradRho).y.segment(iGridStart, n) = output.row(17);
        (*funcData.d3FdRho2dGradRho).z.segment(iGridStart, n) = output.row(18);
        (*funcData.d3FdRhodGradRho2).xx.segment(iGridStart, n) = output.row(19);
        (*funcData.d3FdRhodGradRho2).xy.segment(iGridStart, n) = output.row(20);
        (*funcData.d3FdRhodGradRho2).xz.segment(iGridStart, n) = output.row(21);
        (*funcData.d3FdRhodGradRho2).yy.segment(iGridStart, n) = output.row(22);
        (*funcData.d3FdRhodGradRho2).yz.segment(iGridStart, n) = output.row(23);
        (*funcData.d3FdRhodGradRho2).zz.segment(iGridStart, n) = output.row(24);
        (*funcData.d3FdGradRho3).xxx.segment(iGridStart, n) = output.row(25);
        (*funcData.d3FdGradRho3).xxy.segment(iGridStart, n) = output.row(26);
        (*funcData.d3FdGradRho3).xxz.segment(iGridStart, n) = output.row(27);
        (*funcData.d3FdGradRho3).xyy.segment(iGridStart, n) = output.row(28);
        (*funcData.d3FdGradRho3).xyz.segment(iGridStart, n) = output.row(29);
        (*funcData.d3FdGradRho3).xzz.segment(iGridStart, n) = output.row(30);
        (*funcData.d3FdGradRho3).yyy.segment(iGridStart, n) = output.row(31);
        (*funcData.d3FdGradRho3).yyz.segment(iGridStart, n) = output.row(32);
        (*funcData.d3FdGradRho3).yzz.segment(iGridStart, n) = output.row(33);
        (*funcData.d3FdGradRho3).zzz.segment(iGridStart, n) = output.row(34);
      }
    }
    else {
      throw SerenityError("XCFun (restricted): Output for requested xcfun_vars not implemented.");
    }
  }
}

template<>
void XCFun<Options::SCF_MODES::UNRESTRICTED>::prepareInput(
    const unsigned int iGridStart, const xcfun_vars xcVars, const DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>& density,
    const std::shared_ptr<Gradient<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>>> gradient,
    const std::shared_ptr<Hessian<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>>> hessian, Eigen::MatrixXd& input) {
  unsigned int n = input.cols();
  if (xcVars == xcfun_vars::XC_A_B) {
    input.row(0) = density.alpha.segment(iGridStart, n);
    input.row(1) = density.beta.segment(iGridStart, n);
  }
  else if (xcVars == xcfun_vars::XC_A_B_GAA_GAB_GBB) {
    input.row(0) = density.alpha.segment(iGridStart, n);
    input.row(1) = density.beta.segment(iGridStart, n);
    input.bottomRows(3) = calculateSigma(*gradient, iGridStart, n);
  }
  else if (xcVars == xcfun_vars::XC_A_B_2ND_TAYLOR) {
    input.row(0) = density.alpha.segment(iGridStart, n);
    input.row(1) = gradient->x.alpha.segment(iGridStart, n);
    input.row(2) = gradient->y.alpha.segment(iGridStart, n);
    input.row(3) = gradient->z.alpha.segment(iGridStart, n);
    input.row(4) = hessian->xx.alpha.segment(iGridStart, n);
    input.row(5) = hessian->xy.alpha.segment(iGridStart, n);
    input.row(6) = hessian->xz.alpha.segment(iGridStart, n);
    input.row(7) = hessian->yy.alpha.segment(iGridStart, n);
    input.row(8) = hessian->yz.alpha.segment(iGridStart, n);
    input.row(9) = hessian->zz.alpha.segment(iGridStart, n);
    input.row(10) = density.beta.segment(iGridStart, n);
    input.row(11) = gradient->x.beta.segment(iGridStart, n);
    input.row(12) = gradient->y.beta.segment(iGridStart, n);
    input.row(13) = gradient->z.beta.segment(iGridStart, n);
    input.row(14) = hessian->xx.beta.segment(iGridStart, n);
    input.row(15) = hessian->xy.beta.segment(iGridStart, n);
    input.row(16) = hessian->xz.beta.segment(iGridStart, n);
    input.row(17) = hessian->yy.beta.segment(iGridStart, n);
    input.row(18) = hessian->yz.beta.segment(iGridStart, n);
    input.row(19) = hessian->zz.beta.segment(iGridStart, n);
  }
  else if (xcVars == xcfun_vars::XC_A_B_AX_AY_AZ_BX_BY_BZ) {
    input.row(0) = density.alpha.segment(iGridStart, n);
    input.row(1) = density.beta.segment(iGridStart, n);
    input.row(2) = gradient->x.alpha.segment(iGridStart, n);
    input.row(3) = gradient->y.alpha.segment(iGridStart, n);
    input.row(4) = gradient->z.alpha.segment(iGridStart, n);
    input.row(5) = gradient->x.beta.segment(iGridStart, n);
    input.row(6) = gradient->y.beta.segment(iGridStart, n);
    input.row(7) = gradient->z.beta.segment(iGridStart, n);
  }
  else {
    throw SerenityError("XCFun (unrestricted): Input for requested xcfun_vars not implemented.");
  }
}

template<>
void XCFun<Options::SCF_MODES::UNRESTRICTED>::parseOutput(const unsigned int order, const xcfun_vars xcVars,
                                                          const unsigned int iGridStart, Eigen::MatrixXd& output,
                                                          FunctionalData<Options::SCF_MODES::UNRESTRICTED>& funcData) {
  unsigned int n = output.cols();
  if (funcData.getType() == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
    (*funcData.epuv).segment(iGridStart, n) = output.row(0);
    (*funcData.potential).alpha.segment(iGridStart, n) = output.row(1);
    (*funcData.potential).beta.segment(iGridStart, n) = output.row(2);
  }
  else {
    if (xcVars == xcfun_vars::XC_A_B) {
      (*funcData.epuv).segment(iGridStart, n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).alpha.segment(iGridStart, n) = output.row(1);
        (*funcData.dFdRho).beta.segment(iGridStart, n) = output.row(2);
      }
      if (order >= 2) {
        (*funcData.d2FdRho2).aa.segment(iGridStart, n) = output.row(3);
        (*funcData.d2FdRho2).ab.segment(iGridStart, n) = output.row(4);
        (*funcData.d2FdRho2).bb.segment(iGridStart, n) = output.row(5);
      }
      if (order >= 3) {
        (*funcData.d3FdRho3).aaa.segment(iGridStart, n) = output.row(6);
        (*funcData.d3FdRho3).aab.segment(iGridStart, n) = output.row(7);
        (*funcData.d3FdRho3).abb.segment(iGridStart, n) = output.row(8);
        (*funcData.d3FdRho3).bbb.segment(iGridStart, n) = output.row(9);
      }
    }
    else if (xcVars == xcfun_vars::XC_A_B_GAA_GAB_GBB) {
      (*funcData.epuv).segment(iGridStart, n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).alpha.segment(iGridStart, n) = output.row(1);
        (*funcData.dFdRho).beta.segment(iGridStart, n) = output.row(2);
        (*funcData.dFdSigma).aa.segment(iGridStart, n) = output.row(3);
        (*funcData.dFdSigma).ab.segment(iGridStart, n) = output.row(4);
        (*funcData.dFdSigma).bb.segment(iGridStart, n) = output.row(5);
      }
      if (order >= 2) {
        (*funcData.d2FdRho2).aa.segment(iGridStart, n) = output.row(6);
        (*funcData.d2FdRho2).ab.segment(iGridStart, n) = output.row(7);
        (*funcData.d2FdRho2).bb.segment(iGridStart, n) = output.row(11);
        (*funcData.d2FdRhodSigma).aaa.segment(iGridStart, n) = output.row(8);
        (*funcData.d2FdRhodSigma).aab.segment(iGridStart, n) = output.row(9);
        (*funcData.d2FdRhodSigma).abb.segment(iGridStart, n) = output.row(10);
        (*funcData.d2FdRhodSigma).baa.segment(iGridStart, n) = output.row(12);
        (*funcData.d2FdRhodSigma).bab.segment(iGridStart, n) = output.row(13);
        (*funcData.d2FdRhodSigma).bbb.segment(iGridStart, n) = output.row(14);
        (*funcData.d2FdSigma2).aaaa.segment(iGridStart, n) = output.row(15);
        (*funcData.d2FdSigma2).aaab.segment(iGridStart, n) = output.row(16);
        (*funcData.d2FdSigma2).aabb.segment(iGridStart, n) = output.row(17);
        (*funcData.d2FdSigma2).abab.segment(iGridStart, n) = output.row(18);
        (*funcData.d2FdSigma2).abbb.segment(iGridStart, n) = output.row(19);
        (*funcData.d2FdSigma2).bbbb.segment(iGridStart, n) = output.row(20);
      }
      if (order >= 3) {
        (*funcData.d3FdRho3).aaa.segment(iGridStart, n) = output.row(21);
        (*funcData.d3FdRho3).aab.segment(iGridStart, n) = output.row(22);
        (*funcData.d3FdRho3).abb.segment(iGridStart, n) = output.row(26);
        (*funcData.d3FdRho3).bbb.segment(iGridStart, n) = output.row(36);
        (*funcData.d3FdRho2dSigma).aaaa.segment(iGridStart, n) = output.row(23);
        (*funcData.d3FdRho2dSigma).aaab.segment(iGridStart, n) = output.row(24);
        (*funcData.d3FdRho2dSigma).aabb.segment(iGridStart, n) = output.row(25);
        (*funcData.d3FdRho2dSigma).abaa.segment(iGridStart, n) = output.row(27);
        (*funcData.d3FdRho2dSigma).abab.segment(iGridStart, n) = output.row(28);
        (*funcData.d3FdRho2dSigma).abbb.segment(iGridStart, n) = output.row(29);
        (*funcData.d3FdRho2dSigma).bbaa.segment(iGridStart, n) = output.row(37);
        (*funcData.d3FdRho2dSigma).bbab.segment(iGridStart, n) = output.row(38);
        (*funcData.d3FdRho2dSigma).bbbb.segment(iGridStart, n) = output.row(39);
        (*funcData.d3FdRhodSigma2).aaaaa.segment(iGridStart, n) = output.row(30);
        (*funcData.d3FdRhodSigma2).aaaab.segment(iGridStart, n) = output.row(31);
        (*funcData.d3FdRhodSigma2).aaabb.segment(iGridStart, n) = output.row(32);
        (*funcData.d3FdRhodSigma2).aabab.segment(iGridStart, n) = output.row(33);
        (*funcData.d3FdRhodSigma2).aabbb.segment(iGridStart, n) = output.row(34);
        (*funcData.d3FdRhodSigma2).abbbb.segment(iGridStart, n) = output.row(35);
        (*funcData.d3FdRhodSigma2).baaaa.segment(iGridStart, n) = output.row(40);
        (*funcData.d3FdRhodSigma2).baaab.segment(iGridStart, n) = output.row(41);
        (*funcData.d3FdRhodSigma2).baabb.segment(iGridStart, n) = output.row(42);
        (*funcData.d3FdRhodSigma2).babab.segment(iGridStart, n) = output.row(43);
        (*funcData.d3FdRhodSigma2).babbb.segment(iGridStart, n) = output.row(44);
        (*funcData.d3FdRhodSigma2).bbbbb.segment(iGridStart, n) = output.row(45);
        (*funcData.d3FdSigma3).aaaaaa.segment(iGridStart, n) = output.row(46);
        (*funcData.d3FdSigma3).aaaaab.segment(iGridStart, n) = output.row(47);
        (*funcData.d3FdSigma3).aaaabb.segment(iGridStart, n) = output.row(48);
        (*funcData.d3FdSigma3).aaabab.segment(iGridStart, n) = output.row(49);
        (*funcData.d3FdSigma3).aaabbb.segment(iGridStart, n) = output.row(50);
        (*funcData.d3FdSigma3).aabbbb.segment(iGridStart, n) = output.row(51);
        (*funcData.d3FdSigma3).ababab.segment(iGridStart, n) = output.row(52);
        (*funcData.d3FdSigma3).ababbb.segment(iGridStart, n) = output.row(53);
        (*funcData.d3FdSigma3).abbbbb.segment(iGridStart, n) = output.row(54);
        (*funcData.d3FdSigma3).bbbbbb.segment(iGridStart, n) = output.row(55);
      }
    }
    else if (xcVars == xcfun_vars::XC_A_B_AX_AY_AZ_BX_BY_BZ) {
      (*funcData.epuv).segment(iGridStart, n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).alpha.segment(iGridStart, n) = output.row(1);
        (*funcData.dFdRho).beta.segment(iGridStart, n) = output.row(2);
        (*funcData.dFdGradRho).x.alpha.segment(iGridStart, n) = output.row(3);
        (*funcData.dFdGradRho).y.alpha.segment(iGridStart, n) = output.row(4);
        (*funcData.dFdGradRho).z.alpha.segment(iGridStart, n) = output.row(5);
        (*funcData.dFdGradRho).x.beta.segment(iGridStart, n) = output.row(6);
        (*funcData.dFdGradRho).y.beta.segment(iGridStart, n) = output.row(7);
        (*funcData.dFdGradRho).z.beta.segment(iGridStart, n) = output.row(8);
      }
      if (order >= 2) {
        // Now it gets messy...
        (*funcData.d2FdRho2).aa.segment(iGridStart, n) = output.row(9);
        (*funcData.d2FdRho2).ab.segment(iGridStart, n) = output.row(10);
        (*funcData.d2FdRho2).bb.segment(iGridStart, n) = output.row(17);
        (*funcData.d2FdRhodGradRho).x.aa.segment(iGridStart, n) = output.row(11);
        (*funcData.d2FdRhodGradRho).y.aa.segment(iGridStart, n) = output.row(12);
        (*funcData.d2FdRhodGradRho).z.aa.segment(iGridStart, n) = output.row(13);
        (*funcData.d2FdRhodGradRho).x.ab.segment(iGridStart, n) = output.row(14);
        (*funcData.d2FdRhodGradRho).y.ab.segment(iGridStart, n) = output.row(15);
        (*funcData.d2FdRhodGradRho).z.ab.segment(iGridStart, n) = output.row(16);
        (*funcData.d2FdRhodGradRho).x.ba.segment(iGridStart, n) = output.row(18);
        (*funcData.d2FdRhodGradRho).y.ba.segment(iGridStart, n) = output.row(19);
        (*funcData.d2FdRhodGradRho).z.ba.segment(iGridStart, n) = output.row(20);
        (*funcData.d2FdRhodGradRho).x.bb.segment(iGridStart, n) = output.row(21);
        (*funcData.d2FdRhodGradRho).y.bb.segment(iGridStart, n) = output.row(22);
        (*funcData.d2FdRhodGradRho).z.bb.segment(iGridStart, n) = output.row(23);
        (*funcData.d2FdGradRho2).xx.aa.segment(iGridStart, n) = output.row(24);
        (*funcData.d2FdGradRho2).xy.aa.segment(iGridStart, n) = output.row(25);
        (*funcData.d2FdGradRho2).xz.aa.segment(iGridStart, n) = output.row(26);
        (*funcData.d2FdGradRho2).xx.ab.segment(iGridStart, n) = output.row(27);
        (*funcData.d2FdGradRho2).xy.ab.segment(iGridStart, n) = output.row(28);
        (*funcData.d2FdGradRho2).xz.ab.segment(iGridStart, n) = output.row(29);
        (*funcData.d2FdGradRho2).yy.aa.segment(iGridStart, n) = output.row(30);
        (*funcData.d2FdGradRho2).yz.aa.segment(iGridStart, n) = output.row(31);
        (*funcData.d2FdGradRho2).xy.ba.segment(iGridStart, n) = output.row(32);
        (*funcData.d2FdGradRho2).yy.ab.segment(iGridStart, n) = output.row(33);
        (*funcData.d2FdGradRho2).yz.ab.segment(iGridStart, n) = output.row(34);
        (*funcData.d2FdGradRho2).zz.aa.segment(iGridStart, n) = output.row(35);
        (*funcData.d2FdGradRho2).xz.ba.segment(iGridStart, n) = output.row(36);
        (*funcData.d2FdGradRho2).yz.ba.segment(iGridStart, n) = output.row(37);
        (*funcData.d2FdGradRho2).zz.ab.segment(iGridStart, n) = output.row(38);
        (*funcData.d2FdGradRho2).xx.bb.segment(iGridStart, n) = output.row(39);
        (*funcData.d2FdGradRho2).xy.bb.segment(iGridStart, n) = output.row(40);
        (*funcData.d2FdGradRho2).xz.bb.segment(iGridStart, n) = output.row(41);
        (*funcData.d2FdGradRho2).yy.bb.segment(iGridStart, n) = output.row(42);
        (*funcData.d2FdGradRho2).yz.bb.segment(iGridStart, n) = output.row(43);
        (*funcData.d2FdGradRho2).zz.bb.segment(iGridStart, n) = output.row(44);
        // Copy data to fill Hessian...
        (*funcData.d2FdGradRho2).xx.ba.segment(iGridStart, n) = output.row(27);
        (*funcData.d2FdGradRho2).yy.ba.segment(iGridStart, n) = output.row(33);
        (*funcData.d2FdGradRho2).zz.ba.segment(iGridStart, n) = output.row(38);
      }
      if (order >= 3) {
        // there are 8 variables: rho_a, rho_b, xrho_a, yrho_a, zrho_a, xrho_b, yrho_b, zrho_b (with xrho_a being a
        // shorthand notation for the x-component of the gradient of the alpha density) with this order, the three-digit
        // "number" is incremented under the constraint that the following digit is never "smaller" than the preceding
        // e.g. after rho_a xrho_a zrho_b comes rho_a yrho_a yrho_a. still, this is madness... a a a
        (*funcData.d3FdRho3).aaa.segment(iGridStart, n) = output.row(45);
        // a a b
        (*funcData.d3FdRho3).aab.segment(iGridStart, n) = output.row(46);
        // a a xa
        (*funcData.d3FdRho2dGradRho).x.aaa.segment(iGridStart, n) = output.row(47);
        // a a ya
        (*funcData.d3FdRho2dGradRho).y.aaa.segment(iGridStart, n) = output.row(48);
        (*funcData.d3FdRho2dGradRho).z.aaa.segment(iGridStart, n) = output.row(49);
        (*funcData.d3FdRho2dGradRho).x.aab.segment(iGridStart, n) = output.row(50);
        (*funcData.d3FdRho2dGradRho).y.aab.segment(iGridStart, n) = output.row(51);
        (*funcData.d3FdRho2dGradRho).z.aab.segment(iGridStart, n) = output.row(52);
        (*funcData.d3FdRho3).abb.segment(iGridStart, n) = output.row(53);
        (*funcData.d3FdRho2dGradRho).x.aba.segment(iGridStart, n) = output.row(54);
        (*funcData.d3FdRho2dGradRho).y.aba.segment(iGridStart, n) = output.row(55);
        (*funcData.d3FdRho2dGradRho).z.aba.segment(iGridStart, n) = output.row(56);
        (*funcData.d3FdRho2dGradRho).x.abb.segment(iGridStart, n) = output.row(57);
        (*funcData.d3FdRho2dGradRho).y.abb.segment(iGridStart, n) = output.row(58);
        (*funcData.d3FdRho2dGradRho).z.abb.segment(iGridStart, n) = output.row(59);

        (*funcData.d3FdRhodGradRho2).xx.aaa.segment(iGridStart, n) = output.row(60);
        (*funcData.d3FdRhodGradRho2).xy.aaa.segment(iGridStart, n) = output.row(61);
        (*funcData.d3FdRhodGradRho2).xz.aaa.segment(iGridStart, n) = output.row(62);
        (*funcData.d3FdRhodGradRho2).xx.aab.segment(iGridStart, n) = output.row(63);
        (*funcData.d3FdRhodGradRho2).xy.aab.segment(iGridStart, n) = output.row(64);
        (*funcData.d3FdRhodGradRho2).xz.aab.segment(iGridStart, n) = output.row(65);
        (*funcData.d3FdRhodGradRho2).yy.aaa.segment(iGridStart, n) = output.row(66);
        // a ya za
        (*funcData.d3FdRhodGradRho2).yz.aaa.segment(iGridStart, n) = output.row(67);
        // to follow the pattern, this would be yx.aab, but the last two indices are switched because of the Hessian
        // symmetry (the Hessian class only contains an element xy and not yx) a ya xb
        (*funcData.d3FdRhodGradRho2).xy.aba.segment(iGridStart, n) = output.row(68);
        (*funcData.d3FdRhodGradRho2).yy.aab.segment(iGridStart, n) = output.row(69);
        (*funcData.d3FdRhodGradRho2).yz.aab.segment(iGridStart, n) = output.row(70);
        // a za za
        (*funcData.d3FdRhodGradRho2).zz.aaa.segment(iGridStart, n) = output.row(71);
        // next two are permuted again
        // a za xb
        (*funcData.d3FdRhodGradRho2).xz.aba.segment(iGridStart, n) = output.row(72);
        // a za yb
        (*funcData.d3FdRhodGradRho2).yz.aba.segment(iGridStart, n) = output.row(73);
        (*funcData.d3FdRhodGradRho2).zz.aab.segment(iGridStart, n) = output.row(74);
        (*funcData.d3FdRhodGradRho2).xx.abb.segment(iGridStart, n) = output.row(75);
        (*funcData.d3FdRhodGradRho2).xy.abb.segment(iGridStart, n) = output.row(76);
        (*funcData.d3FdRhodGradRho2).xz.abb.segment(iGridStart, n) = output.row(77);
        (*funcData.d3FdRhodGradRho2).yy.abb.segment(iGridStart, n) = output.row(78);
        (*funcData.d3FdRhodGradRho2).yz.abb.segment(iGridStart, n) = output.row(79);
        (*funcData.d3FdRhodGradRho2).zz.abb.segment(iGridStart, n) = output.row(80);

        (*funcData.d3FdRho3).bbb.segment(iGridStart, n) = output.row(81);

        (*funcData.d3FdRho2dGradRho).x.bba.segment(iGridStart, n) = output.row(82);
        (*funcData.d3FdRho2dGradRho).y.bba.segment(iGridStart, n) = output.row(83);
        (*funcData.d3FdRho2dGradRho).z.bba.segment(iGridStart, n) = output.row(84);
        (*funcData.d3FdRho2dGradRho).x.bbb.segment(iGridStart, n) = output.row(85);
        (*funcData.d3FdRho2dGradRho).y.bbb.segment(iGridStart, n) = output.row(86);
        (*funcData.d3FdRho2dGradRho).z.bbb.segment(iGridStart, n) = output.row(87);

        (*funcData.d3FdRhodGradRho2).xx.baa.segment(iGridStart, n) = output.row(88);
        (*funcData.d3FdRhodGradRho2).xy.baa.segment(iGridStart, n) = output.row(89);
        (*funcData.d3FdRhodGradRho2).xz.baa.segment(iGridStart, n) = output.row(90);
        (*funcData.d3FdRhodGradRho2).xx.bab.segment(iGridStart, n) = output.row(91);
        (*funcData.d3FdRhodGradRho2).xy.bab.segment(iGridStart, n) = output.row(92);
        (*funcData.d3FdRhodGradRho2).xz.bab.segment(iGridStart, n) = output.row(93);
        (*funcData.d3FdRhodGradRho2).yy.baa.segment(iGridStart, n) = output.row(94);
        // b ya za
        (*funcData.d3FdRhodGradRho2).yz.baa.segment(iGridStart, n) = output.row(95);
        // permuted
        // b ya xb
        (*funcData.d3FdRhodGradRho2).xy.bba.segment(iGridStart, n) = output.row(96);
        (*funcData.d3FdRhodGradRho2).yy.bab.segment(iGridStart, n) = output.row(97);
        (*funcData.d3FdRhodGradRho2).yz.bab.segment(iGridStart, n) = output.row(98);
        // b za za
        (*funcData.d3FdRhodGradRho2).zz.baa.segment(iGridStart, n) = output.row(99);
        // next two permuted
        // b za xb
        (*funcData.d3FdRhodGradRho2).xz.bba.segment(iGridStart, n) = output.row(100);
        // b za yb
        (*funcData.d3FdRhodGradRho2).yz.bba.segment(iGridStart, n) = output.row(101);
        (*funcData.d3FdRhodGradRho2).zz.bab.segment(iGridStart, n) = output.row(102);
        (*funcData.d3FdRhodGradRho2).xx.bbb.segment(iGridStart, n) = output.row(103);
        (*funcData.d3FdRhodGradRho2).xy.bbb.segment(iGridStart, n) = output.row(104);
        (*funcData.d3FdRhodGradRho2).xz.bbb.segment(iGridStart, n) = output.row(105);
        (*funcData.d3FdRhodGradRho2).yy.bbb.segment(iGridStart, n) = output.row(106);
        (*funcData.d3FdRhodGradRho2).yz.bbb.segment(iGridStart, n) = output.row(107);
        (*funcData.d3FdRhodGradRho2).zz.bbb.segment(iGridStart, n) = output.row(108);

        (*funcData.d3FdGradRho3).xxx.aaa.segment(iGridStart, n) = output.row(109);
        (*funcData.d3FdGradRho3).xxy.aaa.segment(iGridStart, n) = output.row(110);
        (*funcData.d3FdGradRho3).xxz.aaa.segment(iGridStart, n) = output.row(111);
        (*funcData.d3FdGradRho3).xxx.aab.segment(iGridStart, n) = output.row(112);
        (*funcData.d3FdGradRho3).xxy.aab.segment(iGridStart, n) = output.row(113);
        (*funcData.d3FdGradRho3).xxz.aab.segment(iGridStart, n) = output.row(114);
        (*funcData.d3FdGradRho3).xyy.aaa.segment(iGridStart, n) = output.row(115);
        // xa ya za
        (*funcData.d3FdGradRho3).xyz.aaa.segment(iGridStart, n) = output.row(116);
        // permuted
        // xa ya xb
        (*funcData.d3FdGradRho3).xxy.aba.segment(iGridStart, n) = output.row(117);
        (*funcData.d3FdGradRho3).xyy.aab.segment(iGridStart, n) = output.row(118);
        (*funcData.d3FdGradRho3).xyz.aab.segment(iGridStart, n) = output.row(119);
        // xa za za
        (*funcData.d3FdGradRho3).xzz.aaa.segment(iGridStart, n) = output.row(120);
        // // 2 permuted
        // // xa za xb
        (*funcData.d3FdGradRho3).xxz.aba.segment(iGridStart, n) = output.row(121);
        // // xa za yb
        (*funcData.d3FdGradRho3).xyz.aba.segment(iGridStart, n) = output.row(122);
        (*funcData.d3FdGradRho3).xzz.aab.segment(iGridStart, n) = output.row(123);
        (*funcData.d3FdGradRho3).xxx.abb.segment(iGridStart, n) = output.row(124);
        (*funcData.d3FdGradRho3).xxy.abb.segment(iGridStart, n) = output.row(125);
        (*funcData.d3FdGradRho3).xxz.abb.segment(iGridStart, n) = output.row(126);
        (*funcData.d3FdGradRho3).xyy.abb.segment(iGridStart, n) = output.row(127);
        (*funcData.d3FdGradRho3).xyz.abb.segment(iGridStart, n) = output.row(128);
        (*funcData.d3FdGradRho3).xzz.abb.segment(iGridStart, n) = output.row(129);

        (*funcData.d3FdGradRho3).yyy.aaa.segment(iGridStart, n) = output.row(130);
        // ya ya za
        (*funcData.d3FdGradRho3).yyz.aaa.segment(iGridStart, n) = output.row(131);
        // permuted
        // ya ya xb
        (*funcData.d3FdGradRho3).xyy.baa.segment(iGridStart, n) = output.row(132);
        (*funcData.d3FdGradRho3).yyy.aab.segment(iGridStart, n) = output.row(133);
        (*funcData.d3FdGradRho3).yyz.aab.segment(iGridStart, n) = output.row(134);
        // ya za za
        (*funcData.d3FdGradRho3).yzz.aaa.segment(iGridStart, n) = output.row(135);
        // 2 permuted
        // ya za xb
        (*funcData.d3FdGradRho3).xyz.baa.segment(iGridStart, n) = output.row(136);
        // ya za yb
        (*funcData.d3FdGradRho3).yyz.aba.segment(iGridStart, n) = output.row(137);
        // ya za zb
        (*funcData.d3FdGradRho3).yzz.aab.segment(iGridStart, n) = output.row(138);

        // permuted as well (spatial directions only rising instead of spin indices in order)
        // ya xb xb
        (*funcData.d3FdGradRho3).xxy.bba.segment(iGridStart, n) = output.row(139);
        // ya xb yb
        (*funcData.d3FdGradRho3).xyy.bab.segment(iGridStart, n) = output.row(140);
        // ya xb zb
        (*funcData.d3FdGradRho3).xyz.bab.segment(iGridStart, n) = output.row(141);
        // ya yb yb
        (*funcData.d3FdGradRho3).yyy.abb.segment(iGridStart, n) = output.row(142);
        // ya yb zb
        (*funcData.d3FdGradRho3).yyz.abb.segment(iGridStart, n) = output.row(143);
        // ya zb zb
        (*funcData.d3FdGradRho3).yzz.abb.segment(iGridStart, n) = output.row(144);
        // za za za
        (*funcData.d3FdGradRho3).zzz.aaa.segment(iGridStart, n) = output.row(145);
        // many permuted: z can only be in front if it is followed by only z
        // za za xb
        (*funcData.d3FdGradRho3).xzz.baa.segment(iGridStart, n) = output.row(146);
        // za za yb
        (*funcData.d3FdGradRho3).yzz.baa.segment(iGridStart, n) = output.row(147);
        // za za zb
        (*funcData.d3FdGradRho3).zzz.aab.segment(iGridStart, n) = output.row(148);
        // za xb xb
        (*funcData.d3FdGradRho3).xxz.bba.segment(iGridStart, n) = output.row(149);
        // za xb yb
        (*funcData.d3FdGradRho3).xyz.bba.segment(iGridStart, n) = output.row(150);
        // za xb zb
        (*funcData.d3FdGradRho3).xzz.bab.segment(iGridStart, n) = output.row(151);
        // za yb yb
        (*funcData.d3FdGradRho3).yyz.bba.segment(iGridStart, n) = output.row(152);
        // za yb zb
        (*funcData.d3FdGradRho3).yzz.bab.segment(iGridStart, n) = output.row(153);
        // za zb zb
        (*funcData.d3FdGradRho3).zzz.abb.segment(iGridStart, n) = output.row(154);
        // xb xb xb
        (*funcData.d3FdGradRho3).xxx.bbb.segment(iGridStart, n) = output.row(155);
        (*funcData.d3FdGradRho3).xxy.bbb.segment(iGridStart, n) = output.row(156);
        (*funcData.d3FdGradRho3).xxz.bbb.segment(iGridStart, n) = output.row(157);
        (*funcData.d3FdGradRho3).xyy.bbb.segment(iGridStart, n) = output.row(158);
        (*funcData.d3FdGradRho3).xyz.bbb.segment(iGridStart, n) = output.row(159);
        (*funcData.d3FdGradRho3).xzz.bbb.segment(iGridStart, n) = output.row(160);
        (*funcData.d3FdGradRho3).yyy.bbb.segment(iGridStart, n) = output.row(161);
        (*funcData.d3FdGradRho3).yyz.bbb.segment(iGridStart, n) = output.row(162);
        (*funcData.d3FdGradRho3).yzz.bbb.segment(iGridStart, n) = output.row(163);
        (*funcData.d3FdGradRho3).zzz.bbb.segment(iGridStart, n) = output.row(164);
      }
    }
    else {
      throw SerenityError("XCFun (unrestricted): Output for requested xcfun_vars not implemented.");
    }
  }
}

template<Options::SCF_MODES SCFMode>
xcfun_t* XCFun<SCFMode>::getFunctional(Functional functional) {
  // Status integer for XCFun
  int iErr;
  // Build functional from set of basic functionals
  auto functionals = functional.getBasicFunctionals();
  auto mixingFactors = functional.getMixingFactors();
  auto func = xcfun_new();
  for (unsigned int iFunc = 0; iFunc < functionals.size(); ++iFunc) {
    if (functionals[iFunc] == BasicFunctionals::BASIC_FUNCTIONALS::NONE)
      continue;
    iErr = xcfun_set(func, BasicFunctionals::getXCFunAlias(functionals[iFunc]), mixingFactors[iFunc]);
    if (iErr != 0)
      throw SerenityError("Functional " + std::string(BasicFunctionals::getXCFunAlias(functionals[iFunc])) +
                          " unknown to XCFun. Check if alias is set correctly.");
  }
  // Provide additional parameters for RS-Hybrid
  if (functional.isRSHybrid()) {
    iErr = xcfun_set(func, "cam_alpha", functional.getHfExchangeRatio());
    if (iErr != 0)
      throw SerenityError("XCFun: Failed to set cam_alpha.");
    iErr = xcfun_set(func, "cam_beta", functional.getLRExchangeRatio());
    if (iErr != 0)
      throw SerenityError("XCFun: Failed to set cam_beta.");
    iErr = xcfun_set(func, "rangesep_mu", functional.getRangeSeparationParameter());
    if (iErr != 0)
      throw SerenityError("XCFun: Failed to set rangesep_mu.");
  }
  return func;
}

template<Options::SCF_MODES SCFMode>
double XCFun<SCFMode>::calcEnergy(std::shared_ptr<GridData<RESTRICTED>> epuv, const Eigen::VectorXd& weights) {
  unsigned int nPoints = weights.size();
  unsigned int nBlocks = omp_get_max_threads();
  double energy = 0.0;
#pragma omp parallel for schedule(static, 1) reduction(+ : energy)
  for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
    unsigned int n = (unsigned int)(nPoints / nBlocks);
    const unsigned int start = iBlock * n;
    if (iBlock == nBlocks - 1)
      n += nPoints % nBlocks;
    energy += (*epuv).segment(start, n).dot(weights.segment(start, n));
  }
  return energy;
}

template class XCFun<Options::SCF_MODES::RESTRICTED>;
template class XCFun<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
#endif /* SERENITY_USE_XCFUN */
