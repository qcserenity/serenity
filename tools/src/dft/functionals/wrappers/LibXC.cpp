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
#include "grid/GridController.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <xc.h>
#include <xc_funcs.h>

namespace Serenity {

template<Options::SCF_MODES T>
LibXC<T>::LibXC(unsigned int maxBlockSize) : _maxBlockSize(maxBlockSize) {
}

template<Options::SCF_MODES T>
FunctionalData<T> LibXC<T>::calcData(FUNCTIONAL_DATA_TYPE type, const Functional functional,
                                     const std::shared_ptr<DensityOnGridController<T>> densityOnGridController,
                                     unsigned int order) {
  // Check input
  if (order > 2) {
    throw SerenityError("LibXC usage is only possible up to 2nd order derivatives w.r.t. electrons/sigma.");
  }
  // Build functional from basic functionals
  // auto func = getFunctional(functional);
  const bool gga = functional.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA;
  // Get density
  auto& density = densityOnGridController->getDensityOnGrid();
  // If a GGA is requested, get the gradient of the density
  std::shared_ptr<Gradient<DensityOnGrid<T>>> gradient = nullptr;
  // TODO : implement and enable direct potential evaluation if required
  if (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
    throw SerenityError("Direct Potential evaluation with LibXC not yet implemented.");
  }
  if (gga) {
    if (densityOnGridController->getHighestDerivative() < 1) {
      densityOnGridController->setHighestDerivative(1);
    }
    gradient = std::make_shared<Gradient<DensityOnGrid<T>>>(densityOnGridController->getDensityGradientOnGrid());
  }
  Timings::takeTime("Tech. - LibXC Functional Eval.");

  // Build FunctionalData object
  auto grid = densityOnGridController->getGridController();
  FunctionalData<T> funcData(order, type, functional, grid);
  if (type == FUNCTIONAL_DATA_TYPE::GRADIENTS) {
    if (order >= 1) {
      funcData.dFdSigma = std::make_shared<dF_dSigma<T>>(grid);
    }
    if (order >= 2) {
      funcData.d2FdSigma2 = std::make_shared<d2F_dSigma2<T>>(grid);
      funcData.d2FdRhodSigma = std::make_shared<d2F_dRhodSigma<T>>(grid);
    }
  }
  xc_func_type func;
  // Set density threshold to match xcfun default.
  func.dens_threshold = 1e-14;
  // Translate functional into libxc types
  std::vector<int> functionals;
  for (auto f : functional.getBasicFunctionals())
    functionals.push_back(BasicFunctionals::getLibXCAlias(f));
  std::vector<double> mixing = functional.getMixingFactors();
  // Loop over all basic functional
  for (unsigned int f = 0; f < functionals.size(); f++) {
    if (xc_func_init(&func, functionals[f], (T == RESTRICTED) ? XC_UNPOLARIZED : XC_POLARIZED) != 0)
      throw SerenityError("Error while initializing functional in LibXC.");
    if (functional.isRSHybrid()) {
      // TODO Check how to cover all functionals in LibXC
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
    auto totalDensity = density.total();
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
    f.d2FdGradRho2->xx.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gradient.x.segment(firstIndex, blockSize).array() *
            gradient.x.segment(firstIndex, blockSize).array() +
        2.0 * f.dFdSigma->segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->yy.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gradient.y.segment(firstIndex, blockSize).array() *
            gradient.y.segment(firstIndex, blockSize).array() +
        2.0 * f.dFdSigma->segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->zz.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gradient.z.segment(firstIndex, blockSize).array() *
            gradient.z.segment(firstIndex, blockSize).array() +
        2.0 * f.dFdSigma->segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->xy.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gradient.x.segment(firstIndex, blockSize).array() *
        gradient.y.segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->xz.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gradient.x.segment(firstIndex, blockSize).array() *
        gradient.z.segment(firstIndex, blockSize).array();
    f.d2FdGradRho2->yz.segment(firstIndex, blockSize).array() =
        4.0 * f.d2FdSigma2->segment(firstIndex, blockSize).array() * gradient.y.segment(firstIndex, blockSize).array() *
        gradient.z.segment(firstIndex, blockSize).array();
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
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxb * gya +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxa * gyb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb * gyb +
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gxa * gya;
    f.d2FdGradRho2->xz.aa.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxb * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxa * gzb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxb * gzb +
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gxa * gza;
    f.d2FdGradRho2->yz.aa.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gyb * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gya * gzb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gyb * gzb +
        4.0 * f.d2FdSigma2->aaaa.segment(firstIndex, blockSize).array() * gya * gza;
    f.d2FdGradRho2->xy.bb.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxb * gya +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxa * gyb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa * gya +
        4.0 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gxb * gyb;
    f.d2FdGradRho2->xz.bb.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gxb * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gxa * gzb +
        1.0 * f.d2FdSigma2->abab.segment(firstIndex, blockSize).array() * gxa * gza +
        4.0 * f.d2FdSigma2->bbbb.segment(firstIndex, blockSize).array() * gxb * gzb;
    f.d2FdGradRho2->yz.bb.segment(firstIndex, blockSize).array() =
        2.0 * f.d2FdSigma2->aaab.segment(firstIndex, blockSize).array() * gyb * gza +
        2.0 * f.d2FdSigma2->abbb.segment(firstIndex, blockSize).array() * gya * gzb +
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
      else if (order == 1) {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vxc = Eigen::VectorXd::Zero(blockSize);
        xc_lda_exc_vxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), epuv.data(), vxc.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
        funcData.dFdRho->segment(firstIndex, blockSize).array() += prefactor * vxc.array();
      }
      else if (order == 2) {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vxc = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd fxc = Eigen::VectorXd::Zero(blockSize);
        xc_lda_exc_vxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), epuv.data(), vxc.data());
        xc_lda_fxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), fxc.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
        funcData.dFdRho->segment(firstIndex, blockSize).array() += prefactor * vxc.array();
        funcData.d2FdRho2->segment(firstIndex, blockSize).array() += prefactor * fxc.array();
      }
    }
    else if (ftype == CompositeFunctionals::CLASSES::GGA) {
      if (order == 0) {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::MatrixXd sigma = calculateSigma(*gradient, firstIndex, blockSize);
        xc_gga_exc(&func, blockSize, density.segment(firstIndex, blockSize).data(), sigma.data(), epuv.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
      }
      else if (order == 1) {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vrho = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vsigma = Eigen::VectorXd::Zero(blockSize);
        Eigen::MatrixXd sigma = calculateSigma(*gradient, firstIndex, blockSize);
        xc_gga_exc_vxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), sigma.data(), epuv.data(),
                       vrho.data(), vsigma.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
        funcData.dFdRho->segment(firstIndex, blockSize).array() += prefactor * vrho.array();
        funcData.dFdSigma->segment(firstIndex, blockSize).array() += prefactor * vsigma.array();
      }
      else if (order == 2) {
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vrho = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd vsigma = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd v2rho2 = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd v2rhosigma = Eigen::VectorXd::Zero(blockSize);
        Eigen::VectorXd v2sigma2 = Eigen::VectorXd::Zero(blockSize);
        Eigen::MatrixXd sigma = calculateSigma(*gradient, firstIndex, blockSize);
        xc_gga_exc_vxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), sigma.data(), epuv.data(),
                       vrho.data(), vsigma.data());
        xc_gga_fxc(&func, blockSize, density.segment(firstIndex, blockSize).data(), sigma.data(), v2rho2.data(),
                   v2rhosigma.data(), v2sigma2.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * density.segment(firstIndex, blockSize).array();
        funcData.dFdRho->segment(firstIndex, blockSize).array() += prefactor * vrho.array();
        funcData.dFdSigma->segment(firstIndex, blockSize).array() += prefactor * vsigma.array();
        funcData.d2FdRho2->segment(firstIndex, blockSize).array() += prefactor * v2rho2.array();
        funcData.d2FdSigma2->segment(firstIndex, blockSize).array() += prefactor * v2sigma2.array();
        funcData.d2FdRhodSigma->segment(firstIndex, blockSize).array() += prefactor * v2rhosigma.array();
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

    if (ftype == CompositeFunctionals::CLASSES::LDA) {
      if (order == 0) {
        Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(2, blockSize);
        rho.row(0) = density.alpha.segment(firstIndex, blockSize);
        rho.row(1) = density.beta.segment(firstIndex, blockSize);
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        xc_lda_exc(&func, blockSize, rho.data(), epuv.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
      }
      else if (order == 1) {
        Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(2, blockSize);
        rho.row(0) = density.alpha.segment(firstIndex, blockSize);
        rho.row(1) = density.beta.segment(firstIndex, blockSize);
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::MatrixXd vxc = Eigen::MatrixXd::Zero(2, blockSize);
        xc_lda_exc_vxc(&func, blockSize, rho.data(), epuv.data(), vxc.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
        funcData.dFdRho->alpha.segment(firstIndex, blockSize).array() += prefactor * vxc.row(0).array();
        funcData.dFdRho->beta.segment(firstIndex, blockSize).array() += prefactor * vxc.row(1).array();
      }
      else if (order == 2) {
        Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(2, blockSize);
        rho.row(0) = density.alpha.segment(firstIndex, blockSize);
        rho.row(1) = density.beta.segment(firstIndex, blockSize);
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::MatrixXd vxc = Eigen::MatrixXd::Zero(2, blockSize);
        Eigen::MatrixXd fxc = Eigen::MatrixXd::Zero(3, blockSize);
        xc_lda_exc_vxc(&func, blockSize, rho.data(), epuv.data(), vxc.data());
        xc_lda_fxc(&func, blockSize, rho.data(), fxc.data());
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
        funcData.dFdRho->alpha.segment(firstIndex, blockSize).array() += prefactor * vxc.row(0).array();
        funcData.dFdRho->beta.segment(firstIndex, blockSize).array() += prefactor * vxc.row(1).array();
        funcData.d2FdRho2->aa.segment(firstIndex, blockSize).array() += prefactor * fxc.row(0).array();
        funcData.d2FdRho2->ab.segment(firstIndex, blockSize).array() += prefactor * fxc.row(1).array();
        funcData.d2FdRho2->bb.segment(firstIndex, blockSize).array() += prefactor * fxc.row(2).array();
      }
    }
    else if (ftype == CompositeFunctionals::CLASSES::GGA) {
      if (order == 0) {
        // input variables
        Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(2, blockSize);
        rho.row(0) = density.alpha.segment(firstIndex, blockSize);
        rho.row(1) = density.beta.segment(firstIndex, blockSize);
        Eigen::MatrixXd sigma = calculateSigma(*gradient, firstIndex, blockSize);
        // output variables
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        // evaluation
        xc_gga_exc(&func, blockSize, rho.data(), sigma.data(), epuv.data());
        // parsing
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
      }
      else if (order == 1) {
        // input variables
        Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(2, blockSize);
        rho.row(0) = density.alpha.segment(firstIndex, blockSize);
        rho.row(1) = density.beta.segment(firstIndex, blockSize);
        Eigen::MatrixXd sigma = calculateSigma(*gradient, firstIndex, blockSize);
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
      }
      else if (order == 2) {
        // input variables
        Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(2, blockSize);
        rho.row(0) = density.alpha.segment(firstIndex, blockSize);
        rho.row(1) = density.beta.segment(firstIndex, blockSize);
        Eigen::MatrixXd sigma = calculateSigma(*gradient, firstIndex, blockSize);
        // output variables
        Eigen::VectorXd epuv = Eigen::VectorXd::Zero(blockSize);
        Eigen::MatrixXd vrho = Eigen::MatrixXd::Zero(2, blockSize);
        Eigen::MatrixXd vsigma = Eigen::MatrixXd::Zero(3, blockSize);
        Eigen::MatrixXd v2rho2 = Eigen::MatrixXd::Zero(3, blockSize);
        Eigen::MatrixXd v2rhosigma = Eigen::MatrixXd::Zero(6, blockSize);
        Eigen::MatrixXd v2sigma2 = Eigen::MatrixXd::Zero(6, blockSize);
        // evaluation
        xc_gga_exc_vxc(&func, blockSize, rho.data(), sigma.data(), epuv.data(), vrho.data(), vsigma.data());
        xc_gga_fxc(&func, blockSize, rho.data(), sigma.data(), v2rho2.data(), v2rhosigma.data(), v2sigma2.data());
        // parsing
        funcData.epuv->segment(firstIndex, blockSize).array() +=
            prefactor * epuv.array() * totalDensity.segment(firstIndex, blockSize).array();
        funcData.dFdRho->alpha.segment(firstIndex, blockSize).array() += prefactor * vrho.row(0).array();
        funcData.dFdRho->beta.segment(firstIndex, blockSize).array() += prefactor * vrho.row(1).array();
        funcData.dFdSigma->aa.segment(firstIndex, blockSize).array() += prefactor * vsigma.row(0).array();
        funcData.dFdSigma->ab.segment(firstIndex, blockSize).array() += prefactor * vsigma.row(1).array();
        funcData.dFdSigma->bb.segment(firstIndex, blockSize).array() += prefactor * vsigma.row(2).array();
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
    }
  } /* Loop over blocks */
}

template<Options::SCF_MODES T>
unsigned int LibXC<T>::determineBlockSize(unsigned int blockIndex, unsigned int nPoints, unsigned int nBlocks) {
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

template<Options::SCF_MODES T>
double LibXC<T>::calcEnergy(std::shared_ptr<GridData<RESTRICTED>> epuv, const Eigen::VectorXd& weights) {
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
