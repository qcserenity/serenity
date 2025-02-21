/**
 * @file   LibXC.cpp
 *
 * @date   Feb 24, 2020
 * @author Jan P. Unsleber
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
#ifdef SERENITY_USE_LIBXC
/* Include Class Header*/
#include "dft/functionals/wrappers/LibXC.h"
/* Include Serenity Internal Headers */
#include "dft/functionals/BasicFunctionals.h"
#include "grid/GridController.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <xc.h>
#include <xc_funcs.h>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
LibXC<SCFMode>::LibXC(unsigned int maxBlockSize) : _maxBlockSize(maxBlockSize) {
}

template<Options::SCF_MODES SCFMode>
FunctionalData<SCFMode> LibXC<SCFMode>::calcData(FUNCTIONAL_DATA_TYPE type, const Functional functional,
                                                 const std::shared_ptr<DensityOnGridController<SCFMode>> densityOnGridController,
                                                 unsigned int order) {
  // Check input
  if (order > 3) {
    throw SerenityError("LibXC usage is only possible up to 3rd order derivatives w.r.t. electrons/sigma.");
  }
  // Build functional from basic functionals
  // auto func = getFunctional(functional);
  const bool gga = functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA;
  // Get density
  auto& density = densityOnGridController->getDensityOnGrid();
  // If a GGA is requested, get the gradient of the density
  std::shared_ptr<Gradient<DensityOnGrid<SCFMode>>> gradient = nullptr;
  // TODO : implement and enable direct potential evaluation if required
  if (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
    throw SerenityError("Direct Potential evaluation with LibXC not yet implemented.");
  }
  if (gga) {
    if (densityOnGridController->getHighestDerivative() < 1) {
      densityOnGridController->setHighestDerivative(1);
    }
    gradient = std::make_shared<Gradient<DensityOnGrid<SCFMode>>>(densityOnGridController->getDensityGradientOnGrid());
  }
  Timings::takeTime("Tech. - LibXC Functional Eval.");

  // Build FunctionalData object
  auto grid = densityOnGridController->getGridController();
  FunctionalData<SCFMode> funcData(order, type, functional, grid);
  if (type == FUNCTIONAL_DATA_TYPE::GRADIENTS) {
    if (order >= 1) {
      funcData.dFdSigma = std::make_shared<dF_dSigma<SCFMode>>(grid);
    }
    if (order >= 2) {
      funcData.d2FdSigma2 = std::make_shared<d2F_dSigma2<SCFMode>>(grid);
      funcData.d2FdRhodSigma = std::make_shared<d2F_dRhodSigma<SCFMode>>(grid);
    }
    if (order >= 3) {
      funcData.d3FdRho2dSigma = std::make_shared<d3F_dRho2dSigma<SCFMode>>(grid);
      funcData.d3FdRhodSigma2 = std::make_shared<d3F_dRhodSigma2<SCFMode>>(grid);
      funcData.d3FdSigma3 = std::make_shared<d3F_dSigma3<SCFMode>>(grid);
    }
  }
  xc_func_type func;
  // Set density threshold to match xcfun default.
  func.dens_threshold = 1e-14;
  // Translate functional into libxc types
  std::vector<int> functionals;
  for (auto f : functional.getBasicFunctionals()) {
    if (f == BasicFunctionals::BASIC_FUNCTIONALS::NONE)
      continue;
    functionals.push_back(BasicFunctionals::getLibXCAlias(f));
  }
  std::vector<double> mixing = functional.getMixingFactors();
  // Loop over all basic functional
  for (unsigned int f = 0; f < functionals.size(); f++) {
    if (xc_func_init(&func, functionals[f], (SCFMode == RESTRICTED) ? XC_UNPOLARIZED : XC_POLARIZED) != 0)
      throw SerenityError("Error while initializing functional in LibXC.");
    if (functional.isRSHybrid()) {
      //   Currently it looks like only XC_GGA_X_WPBEH has this _omega parameter
      if (functionals[f] == XC_GGA_X_WPBEH) {
        xc_func_set_ext_params_name(&func, "_omega", functional.getRangeSeparationParameter());
      }
    }
    this->eval(funcData, density, gradient, mixing[f], functional.getFunctionalClass(), func, order);
    xc_func_end(&func);
  } /*  Loop over functionals */

  if (gga && (order > 0) && (type == FUNCTIONAL_DATA_TYPE::GRADIENTS)) {
    // Number of grid points and blocks
    const unsigned int nPoints = density.getNGridPoints();
    const unsigned int nBlocks = (unsigned int)ceil((double)nPoints / _maxBlockSize);
    // Loop over blocks of grid points
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
      // first index of block
      const unsigned int firstIndex = iBlock * _maxBlockSize;
      // size of this block
      const unsigned int blockSize = determineBlockSize(iBlock, nPoints, nBlocks);
      this->complete(funcData, *gradient, firstIndex, blockSize);
    }
  }

  funcData.energy = calcEnergy(funcData.epuv, densityOnGridController->getGridController()->getWeights());
  Timings::timeTaken("Tech. - LibXC Functional Eval.");
  return funcData;
}

template<Options::SCF_MODES SCFMode>
std::string LibXC<SCFMode>::version() {
  return xc_version_string();
}

template<Options::SCF_MODES SCFMode>
std::string LibXC<SCFMode>::reference() {
  return xc_reference();
}

template<Options::SCF_MODES SCFMode>
std::string LibXC<SCFMode>::referenceDOI() {
  return xc_reference_doi();
}

template<>
Eigen::MatrixXd LibXC<RESTRICTED>::calculateSigma(const Gradient<DensityOnGrid<RESTRICTED>>& gradient,
                                                  const unsigned int& iGridStart, const unsigned int& blocksize) {
  Eigen::MatrixXd sigma(1, blocksize);
  sigma.array() = gradient.x.segment(iGridStart, blocksize).array() * gradient.x.segment(iGridStart, blocksize).array() +
                  gradient.y.segment(iGridStart, blocksize).array() * gradient.y.segment(iGridStart, blocksize).array() +
                  gradient.z.segment(iGridStart, blocksize).array() * gradient.z.segment(iGridStart, blocksize).array();
  return sigma;
}

template<>
Eigen::MatrixXd LibXC<UNRESTRICTED>::calculateSigma(const Gradient<DensityOnGrid<UNRESTRICTED>>& gradient,
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
void LibXC<RESTRICTED>::complete(const FunctionalData<RESTRICTED>& f, const Gradient<DensityOnGrid<RESTRICTED>>& gradient,
                                 const unsigned int& firstIndex, const unsigned int& blockSize) {
  if (f.dFdSigma) {
    f.dFdGradRho->x.segment(firstIndex, blockSize).array() =
        2.0 * gradient.x.segment(firstIndex, blockSize).array() * f.dFdSigma->segment(firstIndex, blockSize).array();
    f.dFdGradRho->y.segment(firstIndex, blockSize).array() =
        2.0 * gradient.y.segment(firstIndex, blockSize).array() * f.dFdSigma->segment(firstIndex, blockSize).array();
    f.dFdGradRho->z.segment(firstIndex, blockSize).array() =
        2.0 * gradient.z.segment(firstIndex, blockSize).array() * f.dFdSigma->segment(firstIndex, blockSize).array();
  }
  if (f.d2FdRhodSigma) {
    f.d2FdRhodGradRho->x.segment(firstIndex, blockSize).array() =
        2.0 * gradient.x.segment(firstIndex, blockSize).array() * f.d2FdRhodSigma->segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->y.segment(firstIndex, blockSize).array() =
        2.0 * gradient.y.segment(firstIndex, blockSize).array() * f.d2FdRhodSigma->segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->z.segment(firstIndex, blockSize).array() =
        2.0 * gradient.z.segment(firstIndex, blockSize).array() * f.d2FdRhodSigma->segment(firstIndex, blockSize).array();
  }
  if (f.d2FdSigma2) {
    const auto& gx = gradient.x.segment(firstIndex, blockSize).array();
    const auto& gy = gradient.y.segment(firstIndex, blockSize).array();
    const auto& gz = gradient.z.segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->xx.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gx * gx +
        2.0 * f.dFdSigma->segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->yy.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gy * gy +
        2.0 * f.dFdSigma->segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->zz.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gz * gz +
        2.0 * f.dFdSigma->segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->xy.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gx * gy;
    f.d2FdGradRho2->xz.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gx * gz;
    f.d2FdGradRho2->yz.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gy * gz;
  }
  if (f.d3FdRho2dSigma) {
    f.d3FdRho2dGradRho->x.segment(firstIndex, blockSize).array() =
        2. * f.d3FdRho2dSigma->segment(firstIndex, blockSize).array() * gradient.x.segment(firstIndex, blockSize).array();
    f.d3FdRho2dGradRho->y.segment(firstIndex, blockSize).array() =
        2. * f.d3FdRho2dSigma->segment(firstIndex, blockSize).array() * gradient.y.segment(firstIndex, blockSize).array();
    f.d3FdRho2dGradRho->z.segment(firstIndex, blockSize).array() =
        2. * f.d3FdRho2dSigma->segment(firstIndex, blockSize).array() * gradient.z.segment(firstIndex, blockSize).array();
  }
  if (f.d3FdRhodSigma2) {
    const auto& gx = gradient.x.segment(firstIndex, blockSize).array();
    const auto& gy = gradient.y.segment(firstIndex, blockSize).array();
    const auto& gz = gradient.z.segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->xx.segment(firstIndex, blockSize).array() =
        4. * f.d3FdRhodSigma2->segment(firstIndex, blockSize).array() * gx * gx +
        2. * f.d2FdRhodSigma->segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->yy.segment(firstIndex, blockSize).array() =
        4. * f.d3FdRhodSigma2->segment(firstIndex, blockSize).array() * gy * gy +
        2. * f.d2FdRhodSigma->segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->zz.segment(firstIndex, blockSize).array() =
        4. * f.d3FdRhodSigma2->segment(firstIndex, blockSize).array() * gz * gz +
        2. * f.d2FdRhodSigma->segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->xy.segment(firstIndex, blockSize).array() =
        4. * f.d3FdRhodSigma2->segment(firstIndex, blockSize).array() * gx * gy;
    f.d3FdRhodGradRho2->xz.segment(firstIndex, blockSize).array() =
        4. * f.d3FdRhodSigma2->segment(firstIndex, blockSize).array() * gx * gz;
    f.d3FdRhodGradRho2->yz.segment(firstIndex, blockSize).array() =
        4. * f.d3FdRhodSigma2->segment(firstIndex, blockSize).array() * gy * gz;
  }
  if (f.d3FdSigma3) {
    const auto& gx = gradient.x.segment(firstIndex, blockSize).array();
    const auto& gy = gradient.y.segment(firstIndex, blockSize).array();
    const auto& gz = gradient.z.segment(firstIndex, blockSize).array();
    f.d3FdGradRho3->xxx.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gx * gx * gx +
        12 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gx;
    f.d3FdGradRho3->yyy.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gy * gy * gy +
        12 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gy;
    f.d3FdGradRho3->zzz.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gz * gz * gz +
        12 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gz;
    f.d3FdGradRho3->xxy.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gx * gx * gy +
        4 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gy;
    f.d3FdGradRho3->xxz.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gx * gx * gz +
        4 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gz;
    f.d3FdGradRho3->xyy.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gx * gy * gy +
        4 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gx;
    f.d3FdGradRho3->xyz.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gx * gy * gz;
    f.d3FdGradRho3->xzz.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gx * gz * gz +
        4 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gx;
    f.d3FdGradRho3->yyz.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gy * gy * gz +
        4 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gz;
    f.d3FdGradRho3->yzz.segment(firstIndex, blockSize).array() =
        8 * f.d3FdSigma3->segment(firstIndex, blockSize).array() * gy * gz * gz +
        4 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gy;
  }
}

template<>
void LibXC<UNRESTRICTED>::complete(const FunctionalData<UNRESTRICTED>& f, const Gradient<DensityOnGrid<UNRESTRICTED>>& gradient,
                                   const unsigned int& firstIndex, const unsigned int& blockSize) {
  if (f.dFdSigma) {
    f.dFdGradRho->x.alpha.segment(firstIndex, blockSize).array() =
        2.0 * gradient.x.alpha.segment(firstIndex, blockSize).array() * f.dFdSigma->aa.segment(firstIndex, blockSize).array() +
        gradient.x.beta.segment(firstIndex, blockSize).array() * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
    f.dFdGradRho->y.alpha.segment(firstIndex, blockSize).array() =
        2.0 * gradient.y.alpha.segment(firstIndex, blockSize).array() * f.dFdSigma->aa.segment(firstIndex, blockSize).array() +
        gradient.y.beta.segment(firstIndex, blockSize).array() * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
    f.dFdGradRho->z.alpha.segment(firstIndex, blockSize).array() =
        2.0 * gradient.z.alpha.segment(firstIndex, blockSize).array() * f.dFdSigma->aa.segment(firstIndex, blockSize).array() +
        gradient.z.beta.segment(firstIndex, blockSize).array() * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
    f.dFdGradRho->x.beta.segment(firstIndex, blockSize).array() =
        2.0 * gradient.x.beta.segment(firstIndex, blockSize).array() * f.dFdSigma->bb.segment(firstIndex, blockSize).array() +
        gradient.x.alpha.segment(firstIndex, blockSize).array() * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
    f.dFdGradRho->y.beta.segment(firstIndex, blockSize).array() =
        2.0 * gradient.y.beta.segment(firstIndex, blockSize).array() * f.dFdSigma->bb.segment(firstIndex, blockSize).array() +
        gradient.y.alpha.segment(firstIndex, blockSize).array() * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
    f.dFdGradRho->z.beta.segment(firstIndex, blockSize).array() =
        2.0 * gradient.z.beta.segment(firstIndex, blockSize).array() * f.dFdSigma->bb.segment(firstIndex, blockSize).array() +
        gradient.z.alpha.segment(firstIndex, blockSize).array() * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
  }
  if (f.d2FdRhodSigma) {
    f.d2FdRhodGradRho->x.aa.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->aaa.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.x.alpha.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array() * gradient.x.beta.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->y.aa.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->aaa.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.y.alpha.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array() * gradient.y.beta.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->z.aa.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->aaa.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.z.alpha.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array() * gradient.z.beta.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->x.ba.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->baa.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.x.alpha.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array() * gradient.x.beta.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->y.ba.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->baa.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.y.alpha.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array() * gradient.y.beta.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->z.ba.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->baa.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.z.alpha.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array() * gradient.z.beta.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->x.ab.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->abb.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.x.beta.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array() * gradient.x.alpha.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->y.ab.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->abb.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.y.beta.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array() * gradient.y.alpha.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->z.ab.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->abb.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.z.beta.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array() * gradient.z.alpha.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->x.bb.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->bbb.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.x.beta.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array() * gradient.x.alpha.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->y.bb.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->bbb.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.y.beta.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array() * gradient.y.alpha.segment(firstIndex, blockSize).array();
    f.d2FdRhodGradRho->z.bb.segment(firstIndex, blockSize).array() =
        f.d2FdRhodSigma->bbb.segment(firstIndex, blockSize).array() * 2.0 *
            gradient.z.beta.segment(firstIndex, blockSize).array() +
        f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array() * gradient.z.alpha.segment(firstIndex, blockSize).array();
  }
  if (f.d2FdSigma2) {
    const auto& gxa = gradient.x.alpha.segment(firstIndex, blockSize).array();
    const auto& gya = gradient.y.alpha.segment(firstIndex, blockSize).array();
    const auto& gza = gradient.z.alpha.segment(firstIndex, blockSize).array();
    const auto& gxb = gradient.x.beta.segment(firstIndex, blockSize).array();
    const auto& gyb = gradient.y.beta.segment(firstIndex, blockSize).array();
    const auto& gzb = gradient.z.beta.segment(firstIndex, blockSize).array();
    // clang-format off
    f.d2FdGradRho2->xx.aa.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gxa * gxa +
        2.0 * f.dFdSigma->aa.segment(firstIndex, blockSize).array() +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb * gxb +
        4.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa * gxb;
    f.d2FdGradRho2->yy.aa.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gya * gya +
        2.0 * f.dFdSigma->aa.segment(firstIndex, blockSize).array() +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gyb * gyb +
        4.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya * gyb;
    f.d2FdGradRho2->zz.aa.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gza * gza +
        2.0 * f.dFdSigma->aa.segment(firstIndex, blockSize).array() +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gzb * gzb +
        4.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gza * gzb;
    f.d2FdGradRho2->xx.bb.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gxb * gxb +
        2.0 * f.dFdSigma->bb.segment(firstIndex, blockSize).array() +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa * gxa +
        4.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxa * gxb;
    f.d2FdGradRho2->yy.bb.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gyb * gyb +
        2.0 * f.dFdSigma->bb.segment(firstIndex, blockSize).array() +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gya * gya +
        4.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gya * gyb;
    f.d2FdGradRho2->zz.bb.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gzb * gzb +
        2.0 * f.dFdSigma->bb.segment(firstIndex, blockSize).array() +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gza * gza +
        4.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gza * gzb;
    f.d2FdGradRho2->xx.ab.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa * gxa +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb * gxb +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxa * gxb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb * gxa +
        1.0 * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->yy.ab.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya * gya +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gyb * gyb +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gya * gyb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gyb * gya +
        1.0 * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->zz.ab.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gza * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gzb * gzb +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gza * gzb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gzb * gza +
        1.0 * f.dFdSigma->ab.segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->xy.aa.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb * gyb +
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gxa * gya;
    f.d2FdGradRho2->xz.aa.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb * gzb +
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gxa * gza;
    f.d2FdGradRho2->yz.aa.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gyb * gzb +
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gya * gza;
    f.d2FdGradRho2->xy.bb.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa * gya +
        4.0 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gxb * gyb;
    f.d2FdGradRho2->xz.bb.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa * gza +
        4.0 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gxb * gzb;
    f.d2FdGradRho2->yz.bb.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gya * gza +
        4.0 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gyb * gzb;
    f.d2FdGradRho2->xy.ab.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa * gya +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb * gyb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb * gya +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxa * gyb;
    f.d2FdGradRho2->xz.ab.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb * gzb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb * gza +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxa * gzb;
    f.d2FdGradRho2->yz.ab.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gyb * gzb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gyb * gza +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gya * gzb;
    f.d2FdGradRho2->xy.ba.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa * gya +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb * gyb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa * gyb +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxb * gya;
    f.d2FdGradRho2->xz.ba.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb * gzb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa * gzb +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxb * gza;
    f.d2FdGradRho2->yz.ba.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gyb * gzb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gya * gzb +
        4.0 * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gyb * gza;
    // clang-format on

    // Complete duplicates
    f.d2FdGradRho2->xx.ba.segment(firstIndex, blockSize).array() =
        f.d2FdGradRho2->xx.ab.segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->yy.ba.segment(firstIndex, blockSize).array() =
        f.d2FdGradRho2->yy.ab.segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->zz.ba.segment(firstIndex, blockSize).array() =
        f.d2FdGradRho2->zz.ab.segment(firstIndex, blockSize).array();
  }
  if (f.d3FdRho2dSigma) {
    const auto& gxa = gradient.x.alpha.segment(firstIndex, blockSize).array();
    const auto& gya = gradient.y.alpha.segment(firstIndex, blockSize).array();
    const auto& gza = gradient.z.alpha.segment(firstIndex, blockSize).array();
    const auto& gxb = gradient.x.beta.segment(firstIndex, blockSize).array();
    const auto& gyb = gradient.y.beta.segment(firstIndex, blockSize).array();
    const auto& gzb = gradient.z.beta.segment(firstIndex, blockSize).array();

    f.d3FdRho2dGradRho->x.aaa.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->aaaa.segment(firstIndex, blockSize).array() * gxa +
        f.d3FdRho2dSigma->aaab.segment(firstIndex, blockSize).array() * gxb;
    f.d3FdRho2dGradRho->x.aab.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->aaab.segment(firstIndex, blockSize).array() * gxa +
        2 * f.d3FdRho2dSigma->aabb.segment(firstIndex, blockSize).array() * gxb;
    f.d3FdRho2dGradRho->x.aba.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->abaa.segment(firstIndex, blockSize).array() * gxa +
        f.d3FdRho2dSigma->abab.segment(firstIndex, blockSize).array() * gxb;
    f.d3FdRho2dGradRho->x.abb.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->abab.segment(firstIndex, blockSize).array() * gxa +
        2 * f.d3FdRho2dSigma->abbb.segment(firstIndex, blockSize).array() * gxb;
    f.d3FdRho2dGradRho->x.bba.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->bbaa.segment(firstIndex, blockSize).array() * gxa +
        f.d3FdRho2dSigma->bbab.segment(firstIndex, blockSize).array() * gxb;
    f.d3FdRho2dGradRho->x.bbb.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->bbab.segment(firstIndex, blockSize).array() * gxa +
        2 * f.d3FdRho2dSigma->bbbb.segment(firstIndex, blockSize).array() * gxb;

    f.d3FdRho2dGradRho->y.aaa.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->aaaa.segment(firstIndex, blockSize).array() * gya +
        f.d3FdRho2dSigma->aaab.segment(firstIndex, blockSize).array() * gyb;
    f.d3FdRho2dGradRho->y.aab.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->aaab.segment(firstIndex, blockSize).array() * gya +
        2 * f.d3FdRho2dSigma->aabb.segment(firstIndex, blockSize).array() * gyb;
    f.d3FdRho2dGradRho->y.aba.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->abaa.segment(firstIndex, blockSize).array() * gya +
        f.d3FdRho2dSigma->abab.segment(firstIndex, blockSize).array() * gyb;
    f.d3FdRho2dGradRho->y.abb.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->abab.segment(firstIndex, blockSize).array() * gya +
        2 * f.d3FdRho2dSigma->abbb.segment(firstIndex, blockSize).array() * gyb;
    f.d3FdRho2dGradRho->y.bba.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->bbaa.segment(firstIndex, blockSize).array() * gya +
        f.d3FdRho2dSigma->bbab.segment(firstIndex, blockSize).array() * gyb;
    f.d3FdRho2dGradRho->y.bbb.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->bbab.segment(firstIndex, blockSize).array() * gya +
        2 * f.d3FdRho2dSigma->bbbb.segment(firstIndex, blockSize).array() * gyb;

    f.d3FdRho2dGradRho->z.aaa.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->aaaa.segment(firstIndex, blockSize).array() * gza +
        f.d3FdRho2dSigma->aaab.segment(firstIndex, blockSize).array() * gzb;
    f.d3FdRho2dGradRho->z.aab.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->aaab.segment(firstIndex, blockSize).array() * gza +
        2 * f.d3FdRho2dSigma->aabb.segment(firstIndex, blockSize).array() * gzb;
    f.d3FdRho2dGradRho->z.aba.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->abaa.segment(firstIndex, blockSize).array() * gza +
        f.d3FdRho2dSigma->abab.segment(firstIndex, blockSize).array() * gzb;
    f.d3FdRho2dGradRho->z.abb.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->abab.segment(firstIndex, blockSize).array() * gza +
        2 * f.d3FdRho2dSigma->abbb.segment(firstIndex, blockSize).array() * gzb;
    f.d3FdRho2dGradRho->z.bba.segment(firstIndex, blockSize).array() =
        2 * f.d3FdRho2dSigma->bbaa.segment(firstIndex, blockSize).array() * gza +
        f.d3FdRho2dSigma->bbab.segment(firstIndex, blockSize).array() * gzb;
    f.d3FdRho2dGradRho->z.bbb.segment(firstIndex, blockSize).array() =
        f.d3FdRho2dSigma->bbab.segment(firstIndex, blockSize).array() * gza +
        2 * f.d3FdRho2dSigma->bbbb.segment(firstIndex, blockSize).array() * gzb;
  }
  if (f.d3FdRhodSigma2) {
    const auto& gxa = gradient.x.alpha.segment(firstIndex, blockSize).array();
    const auto& gya = gradient.y.alpha.segment(firstIndex, blockSize).array();
    const auto& gza = gradient.z.alpha.segment(firstIndex, blockSize).array();
    const auto& gxb = gradient.x.beta.segment(firstIndex, blockSize).array();
    const auto& gyb = gradient.y.beta.segment(firstIndex, blockSize).array();
    const auto& gzb = gradient.z.beta.segment(firstIndex, blockSize).array();

    // clang-format off
    f.d3FdRhodGradRho2->xx.aaa.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->aaaaa.segment(firstIndex, blockSize).array() * gxa * gxa +
        2.0 * f.d2FdRhodSigma->aaa.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxb * gxb +
        4.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gxa * gxb;
    f.d3FdRhodGradRho2->xx.baa.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->baaaa.segment(firstIndex, blockSize).array() * gxa * gxa +
        2.0 * f.d2FdRhodSigma->baa.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxb * gxb +
        4.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gxa * gxb;
    f.d3FdRhodGradRho2->yy.aaa.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->aaaaa.segment(firstIndex, blockSize).array() * gya * gya +
        2.0 * f.d2FdRhodSigma->aaa.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gyb * gyb +
        4.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gya * gyb;
    f.d3FdRhodGradRho2->yy.baa.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->baaaa.segment(firstIndex, blockSize).array() * gya * gya +
        2.0 * f.d2FdRhodSigma->baa.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gyb * gyb +
        4.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gya * gyb;
    f.d3FdRhodGradRho2->zz.aaa.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->aaaaa.segment(firstIndex, blockSize).array() * gza * gza +
        2.0 * f.d2FdRhodSigma->aaa.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gzb * gzb +
        4.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gza * gzb;
    f.d3FdRhodGradRho2->zz.baa.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->baaaa.segment(firstIndex, blockSize).array() * gza * gza +
        2.0 * f.d2FdRhodSigma->baa.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gzb * gzb +
        4.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gza * gzb;
    f.d3FdRhodGradRho2->xx.abb.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->abbbb.segment(firstIndex, blockSize).array() * gxb * gxb +
        2.0 * f.d2FdRhodSigma->abb.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxa * gxa +
        4.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gxa * gxb;
    f.d3FdRhodGradRho2->xx.bbb.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->bbbbb.segment(firstIndex, blockSize).array() * gxb * gxb +
        2.0 * f.d2FdRhodSigma->bbb.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxa * gxa +
        4.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gxa * gxb;
    f.d3FdRhodGradRho2->yy.abb.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->abbbb.segment(firstIndex, blockSize).array() * gyb * gyb +
        2.0 * f.d2FdRhodSigma->abb.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gya * gya +
        4.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gya * gyb;
    f.d3FdRhodGradRho2->yy.bbb.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->bbbbb.segment(firstIndex, blockSize).array() * gyb * gyb +
        2.0 * f.d2FdRhodSigma->bbb.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gya * gya +
        4.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gya * gyb;
    f.d3FdRhodGradRho2->zz.abb.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->abbbb.segment(firstIndex, blockSize).array() * gzb * gzb +
        2.0 * f.d2FdRhodSigma->abb.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gza * gza +
        4.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gza * gzb;
    f.d3FdRhodGradRho2->zz.bbb.segment(firstIndex, blockSize).array() =
        4.0 * f.d3FdRhodSigma2->bbbbb.segment(firstIndex, blockSize).array() * gzb * gzb +
        2.0 * f.d2FdRhodSigma->bbb.segment(firstIndex, blockSize).array() +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gza * gza +
        4.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gza * gzb;
    f.d3FdRhodGradRho2->xx.aab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gxa * gxa +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gxb * gxb +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gxa * gxb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxb * gxa +
        1.0 * f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->xx.bab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gxa * gxa +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gxb * gxb +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gxa * gxb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxb * gxa +
        1.0 * f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->yy.aab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gya * gya +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gyb * gyb +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gya * gyb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gyb * gya +
        1.0 * f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->yy.bab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gya * gya +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gyb * gyb +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gya * gyb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gyb * gya +
        1.0 * f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->zz.aab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gza * gza +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gzb * gzb +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gza * gzb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gzb * gza +
        1.0 * f.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->zz.bab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gza * gza +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gzb * gzb +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gza * gzb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gzb * gza +
        1.0 * f.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array();
    f.d3FdRhodGradRho2->xy.aaa.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxb * gyb +
        4.0 * f.d3FdRhodSigma2->aaaaa.segment(firstIndex, blockSize).array() * gxa * gya;
    f.d3FdRhodGradRho2->xy.baa.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxb * gyb +
        4.0 * f.d3FdRhodSigma2->baaaa.segment(firstIndex, blockSize).array() * gxa * gya;
    f.d3FdRhodGradRho2->xz.aaa.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxb * gzb +
        4.0 * f.d3FdRhodSigma2->aaaaa.segment(firstIndex, blockSize).array() * gxa * gza;
    f.d3FdRhodGradRho2->xz.baa.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxb * gzb +
        4.0 * f.d3FdRhodSigma2->baaaa.segment(firstIndex, blockSize).array() * gxa * gza;
    f.d3FdRhodGradRho2->yz.aaa.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gyb * gzb +
        4.0 * f.d3FdRhodSigma2->aaaaa.segment(firstIndex, blockSize).array() * gya * gza;
    f.d3FdRhodGradRho2->yz.baa.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gyb * gzb +
        4.0 * f.d3FdRhodSigma2->baaaa.segment(firstIndex, blockSize).array() * gya * gza;
    f.d3FdRhodGradRho2->xy.abb.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxa * gya +
        4.0 * f.d3FdRhodSigma2->abbbb.segment(firstIndex, blockSize).array() * gxb * gyb;
    f.d3FdRhodGradRho2->xy.bbb.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxa * gya +
        4.0 * f.d3FdRhodSigma2->bbbbb.segment(firstIndex, blockSize).array() * gxb * gyb;
    f.d3FdRhodGradRho2->xz.abb.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxa * gza +
        4.0 * f.d3FdRhodSigma2->abbbb.segment(firstIndex, blockSize).array() * gxb * gzb;
    f.d3FdRhodGradRho2->xz.bbb.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxa * gza +
        4.0 * f.d3FdRhodSigma2->bbbbb.segment(firstIndex, blockSize).array() * gxb * gzb;
    f.d3FdRhodGradRho2->yz.abb.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gya * gza +
        4.0 * f.d3FdRhodSigma2->abbbb.segment(firstIndex, blockSize).array() * gyb * gzb;
    f.d3FdRhodGradRho2->yz.bbb.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gya * gza +
        4.0 * f.d3FdRhodSigma2->bbbbb.segment(firstIndex, blockSize).array() * gyb * gzb;
    f.d3FdRhodGradRho2->xy.aab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gxa * gya +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gxb * gyb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxb * gya +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gxa * gyb;
    f.d3FdRhodGradRho2->xy.bab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gxa * gya +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gxb * gyb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxb * gya +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gxa * gyb;
    f.d3FdRhodGradRho2->xz.aab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gxa * gza +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gxb * gzb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxb * gza +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gxa * gzb;
    f.d3FdRhodGradRho2->xz.bab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gxa * gza +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gxb * gzb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxb * gza +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gxa * gzb;
    f.d3FdRhodGradRho2->yz.aab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gya * gza +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gyb * gzb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gyb * gza +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gya * gzb;
    f.d3FdRhodGradRho2->yz.bab.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gya * gza +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gyb * gzb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gyb * gza +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gya * gzb;
    f.d3FdRhodGradRho2->xy.aba.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gxa * gya +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gxb * gyb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxa * gyb +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gxb * gya;
    f.d3FdRhodGradRho2->xy.bba.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gxa * gya +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gxb * gyb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxa * gyb +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gxb * gya;
    f.d3FdRhodGradRho2->xz.aba.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gxa * gza +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gxb * gzb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gxa * gzb +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gxb * gza;
    f.d3FdRhodGradRho2->xz.bba.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gxa * gza +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gxb * gzb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gxa * gzb +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gxb * gza;
    f.d3FdRhodGradRho2->yz.aba.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize).array() * gya * gza +
        2.0 * f.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize).array() * gyb * gzb +
        1.0 * f.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize).array() * gya * gzb +
        4.0 * f.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize).array() * gyb * gza;
    f.d3FdRhodGradRho2->yz.bba.segment(firstIndex, blockSize).array() =
        2.0 * f.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize).array() * gya * gza +
        2.0 * f.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize).array() * gyb * gzb +
        1.0 * f.d3FdRhodSigma2->babab.segment(firstIndex, blockSize).array() * gya * gzb +
        4.0 * f.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize).array() * gyb * gza;
    // clang-format on
  }
  if (f.d3FdSigma3) {
    const auto& gxa = gradient.x.alpha.segment(firstIndex, blockSize).array();
    const auto& gya = gradient.y.alpha.segment(firstIndex, blockSize).array();
    const auto& gza = gradient.z.alpha.segment(firstIndex, blockSize).array();
    const auto& gxb = gradient.x.beta.segment(firstIndex, blockSize).array();
    const auto& gyb = gradient.y.beta.segment(firstIndex, blockSize).array();
    const auto& gzb = gradient.z.beta.segment(firstIndex, blockSize).array();

    // clang-format off
    f.d3FdGradRho3->xxx.aaa.segment(firstIndex, blockSize).array() = 
        12 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gxa +
        6. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxb +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gxa * gxa * gxa +
        12 * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxa * gxb + 
        6. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxb * gxb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxb * gxb;
    f.d3FdGradRho3->xxx.aab.segment(firstIndex, blockSize).array() =
        6. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa +
        2. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxa * gxa +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gxa * gxb +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxb * gxa +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gxb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxb * gxa +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxb * gxb;
    f.d3FdGradRho3->xxx.abb.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxa +
        2. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa +
        6. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxa * gxa +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxa * gxb +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxa * gxb * gxb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxa * gxa +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxb * gxa +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gxb;
    f.d3FdGradRho3->xxx.bbb.segment(firstIndex, blockSize).array() =
        6. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxa +
        12 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gxb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gxa * gxa +
        6. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gxa * gxb +
        12 * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxa * gxb * gxb +
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gxb;

    f.d3FdGradRho3->xxy.aaa.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gya +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gyb +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gxa * gxa * gya +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxa * gyb +
        8. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxb * gya +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gxb * gya +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxb * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxb * gyb;
    f.d3FdGradRho3->xxy.aab.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gyb +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gxa * gyb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxa * gya +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gyb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxb * gyb +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxb * gya +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxb * gya;
    f.d3FdGradRho3->xxy.aba.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxa * gya +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxa * gyb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gxb * gya +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxb * gyb +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gxb * gya +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gyb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gxa * gya +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxa * gyb +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gyb;
    f.d3FdGradRho3->xxy.abb.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxa * gyb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxa * gya +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gyb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxb * gya +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxa * gxb * gyb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gya +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxa * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxa * gya +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gyb +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gya;
    f.d3FdGradRho3->xxy.bba.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gya +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gyb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gya +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gyb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxa * gya +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gxa * gyb +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gya +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gxb * gyb;
    f.d3FdGradRho3->xxy.bbb.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gyb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gya +
        4. * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gyb +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gya +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gxa * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gxa * gya +
        8. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxa * gxb * gyb +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gxb * gya;

    f.d3FdGradRho3->xxz.aaa.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gza +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gzb +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gxa * gxa * gza +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxa * gzb +
        8. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxb * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gxb * gza +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxb * gzb;
    f.d3FdGradRho3->xxz.aab.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gzb +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gza +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gxa * gzb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxa * gza +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxb * gzb +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxb * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxb * gza;
    f.d3FdGradRho3->xxz.aba.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gxa * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxa * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gxb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxb * gzb +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gxb * gza +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gxa * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxa * gzb +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gza +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gzb;
    f.d3FdGradRho3->xxz.abb.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxa * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxa * gza +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxb * gza +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxa * gxb * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gxa * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gxa * gza +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gzb +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gza;
    f.d3FdGradRho3->xxz.bba.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gza +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gzb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gza +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gxa * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gxa * gzb +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gxb * gza +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gxb * gzb;
    f.d3FdGradRho3->xxz.bbb.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gzb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gxb * gza +
        4. * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gzb +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gxa * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gxa * gza +
        8. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxa * gxb * gzb +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gxb * gza;

// derivative of f.d2FdGradRho2->xy.aa w.r.t. ya
    f.d3FdGradRho3->xyy.aaa.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gya +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gyb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gyb * gya +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gyb * gyb +
        4. * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gxa +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gxa * gya * gya +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gya * gyb;
    f.d3FdGradRho3->xyy.aab.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gyb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gya +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gyb * gya +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gya * gyb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gya * gya;
    f.d3FdGradRho3->xyy.abb.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gya * gyb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gya * gya +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gyb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gya +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gya * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gya * gya +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxa +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxa * gyb * gyb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gyb * gya;
    f.d3FdGradRho3->xyy.baa.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gya * gya +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gya * gyb +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gyb * gya +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gyb * gyb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxb +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxb * gya * gya;
    f.d3FdGradRho3->xyy.bab.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gya * gyb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gya * gya +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gyb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gya +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gyb * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gyb * gya +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxb * gya * gyb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gya * gya;
    f.d3FdGradRho3->xyy.bbb.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxa +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gyb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gya +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gya * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gya * gya +
        4. * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gxb +
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gyb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gya;

    f.d3FdGradRho3->xyz.aaa.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gxa * gya * gza +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gza +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gya * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gyb * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gyb * gzb;
    f.d3FdGradRho3->xyz.aab.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gya * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gzb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gya * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gyb * gza;
    f.d3FdGradRho3->xyz.aba.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gyb * gza +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gya * gza +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gyb * gza +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gyb * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gya * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gya * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gya * gzb;
    f.d3FdGradRho3->xyz.abb.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxa * gyb * gzb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gya * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gya * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gya * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gya * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gya * gza;
    f.d3FdGradRho3->xyz.baa.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxb * gya * gza +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gya * gza +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gyb * gza +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gyb * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gya * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gya * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gyb * gzb;
    f.d3FdGradRho3->xyz.bab.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxb * gya * gzb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gya * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gya * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gya * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gyb * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gya * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gyb * gza;
    f.d3FdGradRho3->xyz.bba.segment(firstIndex, blockSize).array() =
        // bbbbaa
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gza +
        // bbbbab
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gzb +
        // bbabaa
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gya * gza +
        // abbbaa
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gyb * gza +
        // bbabab
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gya * gzb +
        // abbbab
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gyb * gzb +
        // ababaa
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gya * gza +
        // ababab
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gya * gzb;
    f.d3FdGradRho3->xyz.bbb.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gzb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gzb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gyb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * (gxb * gya + gxa * gyb) * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gya * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gya * gza;

    f.d3FdGradRho3->xzz.aaa.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxb * gzb * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gzb * gzb +
        4. * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gxa +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gxa * gza * gza +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gza * gzb;
    f.d3FdGradRho3->xzz.aab.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) * gza +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gzb * gza +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxa * gza * gzb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gza * gza;
    f.d3FdGradRho3->xzz.abb.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gza * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gza * gza +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gzb * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gzb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gza * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxb * gza * gza +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxa +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxa * gzb * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gzb * gza;
    f.d3FdGradRho3->xzz.baa.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxa +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gxa * gza * gza +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gza * gzb +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gzb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gzb * gzb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gxb +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gxb * gza * gza;
    f.d3FdGradRho3->xzz.bab.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxa * gza * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gxa * gza * gza +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gzb * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxb * gzb * gza +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gzb * gza +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gxb * gza * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gxb * gza * gza;
    f.d3FdGradRho3->xzz.bbb.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxa +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * (gxb * gza + gxa * gzb) * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gxa * gza * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gxa * gza * gza +
        4. * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gxb +
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gxb * gzb * gzb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gxb * gzb * gza;

    f.d3FdGradRho3->yyy.aaa.segment(firstIndex, blockSize).array() =
        12 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gya +
        6. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gyb +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gya * gya * gya +
        12 * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gya * gyb +
        6. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gyb * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gyb * gyb;
    f.d3FdGradRho3->yyy.aab.segment(firstIndex, blockSize).array() =
        6. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya +
        2. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gyb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gyb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gya * gya +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gya * gya * gyb +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gyb * gya +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gyb * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gyb * gya +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gyb * gyb;
    f.d3FdGradRho3->yyy.abb.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gya +
        2. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gya +
        6. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gyb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gya * gya +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gya * gyb +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gya * gyb * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gya * gya +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gyb * gya +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gyb * gyb * gyb;
    f.d3FdGradRho3->yyy.bbb.segment(firstIndex, blockSize).array() =
        6. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gya +
        12 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gyb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gya * gya * gya +
        6. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gya * gya * gyb +
        12 * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gya * gyb * gyb +
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gyb * gyb * gyb;

    f.d3FdGradRho3->yyz.aaa.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gza +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gzb +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gya * gya * gza +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gya * gzb +
        8. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gyb * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gyb * gyb * gza +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gyb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gyb * gzb;
    f.d3FdGradRho3->yyz.aab.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gzb +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gza +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gya * gya * gzb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gya * gza +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gyb * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gyb * gzb +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gyb * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gyb * gza;
    f.d3FdGradRho3->yyz.aba.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gya * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gya * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gyb * gyb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gyb * gzb +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gya * gyb * gza +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gyb * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gyb * gya * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gya * gzb +
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gza +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gzb;
    f.d3FdGradRho3->yyz.abb.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gya * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gya * gza +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gyb * gyb * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gyb * gza +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gya * gyb * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gyb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gya * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gya * gza +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gzb +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gza;
    f.d3FdGradRho3->yyz.bba.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gyb * gyb * gza +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gyb * gyb * gzb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gza +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gya * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gya * gya * gzb +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gyb * gza +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gya * gyb * gzb;
    f.d3FdGradRho3->yyz.bbb.segment(firstIndex, blockSize).array() =
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gyb * gyb * gzb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gyb * gyb * gza +
        4. * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gzb +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gya * gya * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gya * gya * gza +
        8. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gya * gyb * gzb +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gya * gyb * gza;

    f.d3FdGradRho3->yzz.aaa.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gyb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) * gza +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gyb * gzb * gza +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gzb * gzb +
        4. * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gya +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gya * gza * gza +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gza * gzb;
    f.d3FdGradRho3->yzz.aab.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) * gza +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gyb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gzb * gza +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gya * gza * gzb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gza * gza;
    f.d3FdGradRho3->yzz.abb.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gza * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gza * gza +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gyb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gyb * gzb * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gzb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gza * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gyb * gza * gza +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gya +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gya * gzb * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gzb * gza;
    f.d3FdGradRho3->yzz.baa.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gya +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gya * gza * gza +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gza * gzb +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gyb * gzb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gya * gzb * gzb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gyb +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gyb * gza * gza;
    f.d3FdGradRho3->yzz.bab.segment(firstIndex, blockSize).array() =
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gya * gza * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gya * gza * gza +
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gyb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gyb * gzb * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gyb * gzb * gza +
        1. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gya +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gya * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gya * gzb * gza +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gyb * gza * gzb +
        4. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gyb * gza * gza;
    f.d3FdGradRho3->yzz.bbb.segment(firstIndex, blockSize).array() =
        2. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gya +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) * gzb +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * (gyb * gza + gya * gzb) * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gya * gza * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gya * gza * gza +
        4. * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gyb +
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gyb * gzb * gzb +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gyb * gzb * gza;

    f.d3FdGradRho3->zzz.aaa.segment(firstIndex, blockSize).array() =
        12 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gza +
        6. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gzb +
        8. * f.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize).array() * gza * gza * gza +
        12 * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gza * gza * gzb +
        6. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gza * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gzb * gzb * gzb;
    f.d3FdGradRho3->zzz.aab.segment(firstIndex, blockSize).array() =
        6. * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gza +
        2. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gzb +
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gzb +
        4. * f.d3FdSigma3->aaaaab.segment(firstIndex, blockSize).array() * gza * gza * gza +
        8. * f.d3FdSigma3->aaaabb.segment(firstIndex, blockSize).array() * gza * gza * gzb +
        4. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gza * gzb * gza +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gza * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gzb * gzb * gza +
        2. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gzb * gzb * gzb;
    f.d3FdGradRho3->zzz.abb.segment(firstIndex, blockSize).array() =
        4. * f.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() * gza +
        2. * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gza +
        6. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gzb +
        2. * f.d3FdSigma3->aaabab.segment(firstIndex, blockSize).array() * gza * gza * gza +
        8. * f.d3FdSigma3->aaabbb.segment(firstIndex, blockSize).array() * gza * gza * gzb +
        8. * f.d3FdSigma3->aabbbb.segment(firstIndex, blockSize).array() * gza * gzb * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gzb * gza * gza +
        4. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gzb * gzb * gza +
        4. * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gzb * gzb * gzb;
    f.d3FdGradRho3->zzz.bbb.segment(firstIndex, blockSize).array() =
        6. * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gza +
        12 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gzb +
        1. * f.d3FdSigma3->ababab.segment(firstIndex, blockSize).array() * gza * gza * gza +
        6. * f.d3FdSigma3->ababbb.segment(firstIndex, blockSize).array() * gza * gza * gzb +
        12 * f.d3FdSigma3->abbbbb.segment(firstIndex, blockSize).array() * gza * gzb * gzb +
        8. * f.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize).array() * gzb * gzb * gzb;
    // clang-format on
  }
}

template<>
void LibXC<RESTRICTED>::eval(const FunctionalData<RESTRICTED>& funcData, const DensityOnGrid<RESTRICTED>& density,
                             std::shared_ptr<Gradient<DensityOnGrid<RESTRICTED>>> gradient, const double& prefactor,
                             const CompositeFunctionals::CLASSES ftype, const xc_func_type& func, unsigned int order) {
  // Number of grid points and blocks
  const unsigned int nPoints = density.getNGridPoints();
  const unsigned int nBlocks = (unsigned int)ceil((double)nPoints / _maxBlockSize);
  // Loop over blocks of grid points
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    // first index of block
    const unsigned int firstIndex = iBlock * _maxBlockSize;
    // size of this block
    const unsigned int blockSize = determineBlockSize(iBlock, nPoints, nBlocks);
    bool skip = true;
    skip = skip && (density.segment(firstIndex, blockSize).array().abs().sum() < blockSize * 1e-12);
    if (skip)
      continue;

    if (ftype == CompositeFunctionals::CLASSES::LDA) {
      if (order == 0) {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        xc_lda_exc(&func, blockSize, density.segment(firstIndex, blockSize).data(), epuv.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
      }
      else {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vxc = Eigen::VectorXd::Zero(blockSize);
        xc_lda_exc_vxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), epuv.data(), vxc.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
        funcData.dFdRho->segment(firstIndex, blockSize).array() += prefactor * vxc.array();
        if (order >= 2) {
          Eigen::VectorXd fxc = Eigen::VectorXd::Zero(blockSize);
          xc_lda_fxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), fxc.data());
          funcData.d2FdRho2->segment(firstIndex, blockSize).array() += prefactor * fxc.array();
        }
        if (order >= 3) {
          Eigen::VectorXd kxc = Eigen::VectorXd::Zero(blockSize);
          xc_lda_kxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), kxc.data());
          funcData.d3FdRho3->segment(firstIndex, blockSize).array() += prefactor * kxc.array();
        }
      }
    }
    else if (ftype == CompositeFunctionals::CLASSES::GGA) {
      Eigen::MatrixXd sigma = calculateSigma(*gradient, firstIndex, blockSize);
      if (order == 0) {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        xc_gga_exc(&func, blockSize, density.segment(firstIndex, blockSize).data(), sigma.data(), epuv.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
      }
      else {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vrho = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vsigma = Eigen::VectorXd::Zero(blockSize);
        xc_gga_exc_vxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), sigma.data(), epuv.data(),
                       vrho.data(), vsigma.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
        funcData.dFdRho->segment(firstIndex, blockSize).array() += prefactor * vrho.array();
        funcData.dFdSigma->segment(firstIndex, blockSize).array() += prefactor * vsigma.array();
        if (order >= 2) {
          Eigen::VectorXd v2rho2 = Eigen::VectorXd::Zero(blockSize);
          Eigen::VectorXd v2rhosigma = Eigen::VectorXd::Zero(blockSize);
          Eigen::VectorXd v2sigma2 = Eigen::VectorXd::Zero(blockSize);
          xc_gga_fxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), sigma.data(), v2rho2.data(),
                     v2rhosigma.data(), v2sigma2.data());
          funcData.d2FdRho2->segment(firstIndex, blockSize).array() += prefactor * v2rho2.array();
          funcData.d2FdSigma2->segment(firstIndex, blockSize).array() += prefactor * v2sigma2.array();
          funcData.d2FdRhodSigma->segment(firstIndex, blockSize).array() += prefactor * v2rhosigma.array();
        }
        if (order >= 3) {
          Eigen::VectorXd v3rho3 = Eigen::VectorXd::Zero(blockSize);
          Eigen::VectorXd v3rho2sigma = Eigen::VectorXd::Zero(blockSize);
          Eigen::VectorXd v3rhosigma2 = Eigen::VectorXd::Zero(blockSize);
          Eigen::VectorXd v3sigma3 = Eigen::VectorXd::Zero(blockSize);
          xc_gga_kxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), sigma.data(), v3rho3.data(),
                     v3rho2sigma.data(), v3rhosigma2.data(), v3sigma3.data());
          funcData.d3FdRho3->segment(firstIndex, blockSize) += prefactor * v3rho3;
          funcData.d3FdRho2dSigma->segment(firstIndex, blockSize) += prefactor * v3rho2sigma;
          funcData.d3FdRhodSigma2->segment(firstIndex, blockSize) += prefactor * v3rhosigma2;
          funcData.d3FdSigma3->segment(firstIndex, blockSize) += prefactor * v3sigma3;
        }
      }
    }
  } /* Loop over blocks */
}

template<>
void LibXC<UNRESTRICTED>::eval(const FunctionalData<UNRESTRICTED>& funcData, const DensityOnGrid<UNRESTRICTED>& density,
                               std::shared_ptr<Gradient<DensityOnGrid<UNRESTRICTED>>> gradient, const double& prefactor,
                               const CompositeFunctionals::CLASSES ftype, const xc_func_type& func, unsigned int order) {
  // Number of grid points and blocks
  const unsigned int nPoints = density.getNGridPoints();
  const unsigned int nBlocks = (unsigned int)ceil((double)nPoints / _maxBlockSize);
  auto totalDensity = density.total();
  // Loop over blocks of grid points
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    // first index of block
    const unsigned int firstIndex = iBlock * _maxBlockSize;
    // size of this block
    const unsigned int blockSize = determineBlockSize(iBlock, nPoints, nBlocks);
    bool skip = true;
    skip = skip && (density.alpha.segment(firstIndex, blockSize).array().abs().sum() < blockSize * 1e-12);
    skip = skip && (density.beta.segment(firstIndex, blockSize).array().abs().sum() < blockSize * 1e-12);
    if (skip)
      continue;

    // input density
    Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(2, blockSize);
    rho.row(0) = density.alpha.segment(firstIndex, blockSize);
    rho.row(1) = density.beta.segment(firstIndex, blockSize);

    if (ftype == CompositeFunctionals::CLASSES::LDA) {
      if (order == 0) {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        xc_lda_exc(&func, blockSize, rho.data(), epuv.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
      }
      else {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::MatrixXd vxc = Eigen::MatrixXd::Zero(2, blockSize);
        xc_lda_exc_vxc(&func, blockSize, rho.data(), epuv.data(), vxc.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
        funcData.dFdRho->alpha.segment(firstIndex, blockSize).array() += prefactor * vxc.row(0).array();
        funcData.dFdRho->beta.segment(firstIndex, blockSize).array() += prefactor * vxc.row(1).array();
        if (order >= 2) {
          Eigen::MatrixXd fxc = Eigen::MatrixXd::Zero(3, blockSize);
          xc_lda_fxc(&func, blockSize, rho.data(), fxc.data());
          funcData.d2FdRho2->aa.segment(firstIndex, blockSize).array() += prefactor * fxc.row(0).array();
          funcData.d2FdRho2->ab.segment(firstIndex, blockSize).array() += prefactor * fxc.row(1).array();
          funcData.d2FdRho2->bb.segment(firstIndex, blockSize).array() += prefactor * fxc.row(2).array();
        }
        if (order >= 3) {
          Eigen::MatrixXd kxc = Eigen::MatrixXd::Zero(4, blockSize);
          xc_lda_kxc(&func, blockSize, rho.data(), kxc.data());
          funcData.d3FdRho3->aaa.segment(firstIndex, blockSize) += prefactor * kxc.row(0);
          funcData.d3FdRho3->aab.segment(firstIndex, blockSize) += prefactor * kxc.row(1);
          funcData.d3FdRho3->abb.segment(firstIndex, blockSize) += prefactor * kxc.row(2);
          funcData.d3FdRho3->bbb.segment(firstIndex, blockSize) += prefactor * kxc.row(3);
        }
      }
    }
    else if (ftype == CompositeFunctionals::CLASSES::GGA) {
      // input sigma
      Eigen::MatrixXd sigma = calculateSigma(*gradient, firstIndex, blockSize);
      if (order == 0) {
        // output variables
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        // evaluation
        xc_gga_exc(&func, blockSize, rho.data(), sigma.data(), epuv.data());
        // parsing
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
      }
      else {
        // output variables
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::MatrixXd vrho = Eigen::MatrixXd::Zero(2, blockSize);
        Eigen::MatrixXd vsigma = Eigen::MatrixXd::Zero(3, blockSize);
        // evaluation
        xc_gga_exc_vxc(&func, blockSize, rho.data(), sigma.data(), epuv.data(), vrho.data(), vsigma.data());
        // parsing
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
        funcData.dFdRho->alpha.segment(firstIndex, blockSize).array() += prefactor * vrho.row(0).array();
        funcData.dFdRho->beta.segment(firstIndex, blockSize).array() += prefactor * vrho.row(1).array();
        funcData.dFdSigma->aa.segment(firstIndex, blockSize).array() += prefactor * vsigma.row(0).array();
        funcData.dFdSigma->ab.segment(firstIndex, blockSize).array() += prefactor * vsigma.row(1).array();
        funcData.dFdSigma->bb.segment(firstIndex, blockSize).array() += prefactor * vsigma.row(2).array();
        if (order >= 2) {
          // output variables
          Eigen::MatrixXd v2rho2 = Eigen::MatrixXd::Zero(3, blockSize);
          Eigen::MatrixXd v2rhosigma = Eigen::MatrixXd::Zero(6, blockSize);
          Eigen::MatrixXd v2sigma2 = Eigen::MatrixXd::Zero(6, blockSize);
          // evaluation
          xc_gga_fxc(&func, blockSize, rho.data(), sigma.data(), v2rho2.data(), v2rhosigma.data(), v2sigma2.data());
          // parsing
          funcData.d2FdRho2->aa.segment(firstIndex, blockSize).array() += prefactor * v2rho2.row(0).array();
          funcData.d2FdRho2->ab.segment(firstIndex, blockSize).array() += prefactor * v2rho2.row(1).array();
          funcData.d2FdRho2->bb.segment(firstIndex, blockSize).array() += prefactor * v2rho2.row(2).array();
          funcData.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() += prefactor * v2sigma2.row(0).array();
          funcData.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() += prefactor * v2sigma2.row(1).array();
          funcData.d2FdSigma2->aabb.segment(firstIndex, blockSize).array() += prefactor * v2sigma2.row(2).array();
          funcData.d2FdSigma2->abab.segment(firstIndex, blockSize).array() += prefactor * v2sigma2.row(3).array();
          funcData.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() += prefactor * v2sigma2.row(4).array();
          funcData.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() += prefactor * v2sigma2.row(5).array();
          funcData.d2FdRhodSigma->aaa.segment(firstIndex, blockSize).array() += prefactor * v2rhosigma.row(0).array();
          funcData.d2FdRhodSigma->aab.segment(firstIndex, blockSize).array() += prefactor * v2rhosigma.row(1).array();
          funcData.d2FdRhodSigma->abb.segment(firstIndex, blockSize).array() += prefactor * v2rhosigma.row(2).array();
          funcData.d2FdRhodSigma->baa.segment(firstIndex, blockSize).array() += prefactor * v2rhosigma.row(3).array();
          funcData.d2FdRhodSigma->bab.segment(firstIndex, blockSize).array() += prefactor * v2rhosigma.row(4).array();
          funcData.d2FdRhodSigma->bbb.segment(firstIndex, blockSize).array() += prefactor * v2rhosigma.row(5).array();
        }
        if (order >= 3) {
          // output variables
          Eigen::MatrixXd v3rho3 = Eigen::MatrixXd::Zero(4, blockSize);
          Eigen::MatrixXd v3rho2sigma = Eigen::MatrixXd::Zero(9, blockSize);
          Eigen::MatrixXd v3rhosigma2 = Eigen::MatrixXd::Zero(12, blockSize);
          Eigen::MatrixXd v3sigma3 = Eigen::MatrixXd::Zero(10, blockSize);
          // evaluation
          xc_gga_kxc(&func, blockSize, rho.data(), sigma.data(), v3rho3.data(), v3rho2sigma.data(), v3rhosigma2.data(),
                     v3sigma3.data());
          // parsing
          funcData.d3FdRho3->aaa.segment(firstIndex, blockSize) += prefactor * v3rho3.row(0);
          funcData.d3FdRho3->aab.segment(firstIndex, blockSize) += prefactor * v3rho3.row(1);
          funcData.d3FdRho3->abb.segment(firstIndex, blockSize) += prefactor * v3rho3.row(2);
          funcData.d3FdRho3->bbb.segment(firstIndex, blockSize) += prefactor * v3rho3.row(3);
          funcData.d3FdRho2dSigma->aaaa.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(0);
          funcData.d3FdRho2dSigma->aaab.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(1);
          funcData.d3FdRho2dSigma->aabb.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(2);
          funcData.d3FdRho2dSigma->abaa.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(3);
          funcData.d3FdRho2dSigma->abab.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(4);
          funcData.d3FdRho2dSigma->abbb.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(5);
          funcData.d3FdRho2dSigma->bbaa.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(6);
          funcData.d3FdRho2dSigma->bbab.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(7);
          funcData.d3FdRho2dSigma->bbbb.segment(firstIndex, blockSize) += prefactor * v3rho2sigma.row(8);
          funcData.d3FdRhodSigma2->aaaaa.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(0);
          funcData.d3FdRhodSigma2->aaaab.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(1);
          funcData.d3FdRhodSigma2->aaabb.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(2);
          funcData.d3FdRhodSigma2->aabab.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(3);
          funcData.d3FdRhodSigma2->aabbb.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(4);
          funcData.d3FdRhodSigma2->abbbb.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(5);
          funcData.d3FdRhodSigma2->baaaa.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(6);
          funcData.d3FdRhodSigma2->baaab.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(7);
          funcData.d3FdRhodSigma2->baabb.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(8);
          funcData.d3FdRhodSigma2->babab.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(9);
          funcData.d3FdRhodSigma2->babbb.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(10);
          funcData.d3FdRhodSigma2->bbbbb.segment(firstIndex, blockSize) += prefactor * v3rhosigma2.row(11);
          funcData.d3FdSigma3->aaaaaa.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(0);
          funcData.d3FdSigma3->aaaaab.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(1);
          funcData.d3FdSigma3->aaaabb.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(2);
          funcData.d3FdSigma3->aaabab.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(3);
          funcData.d3FdSigma3->aaabbb.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(4);
          funcData.d3FdSigma3->aabbbb.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(5);
          funcData.d3FdSigma3->ababab.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(6);
          funcData.d3FdSigma3->ababbb.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(7);
          funcData.d3FdSigma3->abbbbb.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(8);
          funcData.d3FdSigma3->bbbbbb.segment(firstIndex, blockSize) += prefactor * v3sigma3.row(9);
        }
      }
    }
  } /* Loop over blocks */
}

template<Options::SCF_MODES SCFMode>
unsigned int LibXC<SCFMode>::determineBlockSize(unsigned int blockIndex, unsigned int nPoints, unsigned int nBlocks) {
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

template<Options::SCF_MODES SCFMode>
double LibXC<SCFMode>::calcEnergy(std::shared_ptr<GridData<RESTRICTED>> epuv, const Eigen::VectorXd& weights) {
  unsigned int nPoints = weights.size();
  unsigned int nBlocks = omp_get_max_threads();
  double energy = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : energy)
  for (unsigned int iBlock = 0; iBlock < nBlocks; iBlock++) {
    unsigned int n = (unsigned int)(nPoints / nBlocks);
    const unsigned int start = iBlock * n;
    if (iBlock == nBlocks - 1)
      n += nPoints % nBlocks;
    energy += (*epuv).segment(start, n).dot(weights.segment(start, n));
  }
  return energy;
}

template class LibXC<RESTRICTED>;
template class LibXC<UNRESTRICTED>;

} /* namespace Serenity */
#endif /* SERENITY_USE_LIBXC */
