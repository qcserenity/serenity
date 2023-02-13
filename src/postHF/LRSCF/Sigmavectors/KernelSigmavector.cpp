/**
 * @file KernelSigmavector.cpp
 *
 * @date Dec 07, 2018
 * @author Johannes Tölle
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
#include "postHF/LRSCF/Sigmavectors/KernelSigmavector.h"

/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "misc/HelperFunctions.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
#include "postHF/LRSCF/Kernel/Kernel.h"
#include "postHF/LRSCF/LRSCFController.h"
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
KernelSigmavector<SCFMode>::KernelSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                              std::vector<Eigen::MatrixXd> b, std::shared_ptr<Kernel<SCFMode>> kernel)
  : Sigmavector<SCFMode>(lrscf, b), _kernel(kernel) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
KernelSigmavector<SCFMode>::calcF(unsigned I, unsigned J,
                                  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Sigmavector: Kernel");

  // This is needed for coupling pattern where the same system is with coupled
  // and uncoupled vectors. The Kernel has only system information while I and J are
  // related to the lrscf Controller, but the systems behind I and J could be indentical
  for (unsigned int ilrscf = 0; ilrscf < this->_lrscf.size(); ilrscf++) {
    if (this->_lrscf[ilrscf]->getSys() == this->_lrscf[I]->getSys()) {
      I = ilrscf;
      break;
    }
  }
  for (unsigned int ilrscf = 0; ilrscf < this->_lrscf.size(); ilrscf++) {
    if (this->_lrscf[ilrscf]->getSys() == this->_lrscf[J]->getSys()) {
      J = ilrscf;
      break;
    }
  }

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  if (!_kernel) {
    WarningTracker::printWarning("A kernel sigma vector was requested with no kernel present.", true);
    return fock;
  }

  // Thread safety.
  std::vector<std::vector<std::vector<MatrixInBasis<SCFMode>>>> Fxc(this->_nThreads);
  for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
    Fxc[iThread] = std::vector<std::vector<MatrixInBasis<SCFMode>>>(this->_nSet);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Fxc[iThread][iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
    }
  }

  bool gga = _kernel->isGGA();
  unsigned derivativeLevel = (gga) ? 1 : 0;

  std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerI;
  std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerJ;

  basisFunctionOnGridControllerI = BasisFunctionOnGridControllerFactory::produce(
      this->_kernel->getBlocksize(I), this->_kernel->getbasFuncRadialThreshold(I), derivativeLevel,
      this->_lrscf[I]->getBasisController(), _kernel->getGridController());

  if (I == J) {
    basisFunctionOnGridControllerJ = basisFunctionOnGridControllerI;
  }
  else {
    basisFunctionOnGridControllerJ = BasisFunctionOnGridControllerFactory::produce(
        this->_kernel->getBlocksize(J), this->_kernel->getbasFuncRadialThreshold(J), derivativeLevel,
        this->_lrscf[J]->getBasisController(), _kernel->getGridController());
  }

  std::vector<std::vector<GridPotential<SCFMode>>> scalarContr(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      scalarContr[iSet].emplace_back(_kernel->getGridController());
    }
  }

  std::vector<std::vector<Gradient<GridPotential<SCFMode>>>> gradientContr(this->_nSet);
  if (gga) {
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        gradientContr[iSet].emplace_back(makeGradient<GridPotential<SCFMode>>(_kernel->getGridController()));
      }
    }
  }

  Timings::takeTime("LRSCF -   Kernel: Contraction");
  this->contractKernel((*densityMatrices), basisFunctionOnGridControllerJ, scalarContr, gradientContr, I, J);
  Timings::timeTaken("LRSCF -   Kernel: Contraction");

  Timings::takeTime("LRSCF -   Kernel: Num. Integ.");
  this->numericalIntegration(Fxc, basisFunctionOnGridControllerI, scalarContr, gradientContr, I);
  Timings::timeTaken("LRSCF -   Kernel: Num. Integ.");

  // Sum over threads.
  for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto& F = (*fock)[iSet][iGuess];
        auto& fxc = Fxc[iThread][iSet][iGuess];
        for_spin(F, fxc) {
          F_spin += fxc_spin;
        };
      }
    }
  }

  Timings::timeTaken("LRSCF -   Sigmavector: Kernel");

  return fock;
}

template<Options::SCF_MODES SCFMode>
void KernelSigmavector<SCFMode>::contractKernel(std::vector<std::vector<MatrixInBasis<SCFMode>>>& dens,
                                                std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                                                std::vector<std::vector<GridPotential<SCFMode>>>& scalarPart,
                                                std::vector<std::vector<Gradient<GridPotential<SCFMode>>>>& gradientPart,
                                                unsigned int I, unsigned int J) {
  unsigned nBlocks = basisFunctionOnGridController->getNBlocks();
  unsigned nb = basisFunctionOnGridController->getNBasisFunctions();
  bool gga = _kernel->isGGA();
  Eigen::setNbThreads(1);

#pragma omp parallel for schedule(dynamic)
  for (unsigned iBlock = 0; iBlock < nBlocks; ++iBlock) {
    auto& blockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
    auto& funcValues = blockData->functionValues;
    auto& gradValues = blockData->derivativeValues;
    unsigned size = blockData->functionValues.rows();
    unsigned start = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
    auto& weights = basisFunctionOnGridController->getGridController()->getWeights();
    auto projector = constructProjectionMatrix(blockData->negligible);

    Eigen::MatrixXd funcProj = funcValues * projector;
    Eigen::MatrixXd gradProjx;
    Eigen::MatrixXd gradProjy;
    Eigen::MatrixXd gradProjz;
    if (gga) {
      gradProjx = gradValues->x * projector;
      gradProjy = gradValues->y * projector;
      gradProjz = gradValues->z * projector;
    }

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        SpinPolarizedData<SCFMode, Eigen::VectorXd> pbb(size);
        Gradient<SpinPolarizedData<SCFMode, Eigen::VectorXd>> pnbb =
            makeGradient<SpinPolarizedData<SCFMode, Eigen::VectorXd>>((gga) ? size : 0);

        auto& p = dens[iSet][iGuess];
        SpinPolarizedData<SCFMode, Eigen::MatrixXd> pb(size, nb);
        for_spin(p, pb, pbb) {
          pb_spin = (p_spin.transpose() + p_spin) * funcProj.transpose();
          pbb_spin = 0.5 * (funcProj * pb_spin).diagonal();
          pbb_spin = pbb_spin.cwiseProduct(weights.segment(start, size));
        };

        if (gga) {
          auto& pnbbx = pnbb.x;
          auto& pnbby = pnbb.y;
          auto& pnbbz = pnbb.z;
          for_spin(pb, pnbbx, pnbby, pnbbz) {
            pnbbx_spin = (gradProjx * pb_spin).diagonal();
            pnbby_spin = (gradProjy * pb_spin).diagonal();
            pnbbz_spin = (gradProjz * pb_spin).diagonal();
            pnbbx_spin = pnbbx_spin.cwiseProduct(weights.segment(start, size));
            pnbby_spin = pnbby_spin.cwiseProduct(weights.segment(start, size));
            pnbbz_spin = pnbbz_spin.cwiseProduct(weights.segment(start, size));
          };
        }
        this->contractBlock(start, size, pbb, pnbb, scalarPart[iSet][iGuess], gradientPart[iSet][iGuess], I, J);
      }
    }
  }
  Eigen::setNbThreads(0);
}

template<Options::SCF_MODES SCFMode>
void KernelSigmavector<SCFMode>::numericalIntegration(std::vector<std::vector<std::vector<MatrixInBasis<SCFMode>>>>& Fxc,
                                                      std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                                                      std::vector<std::vector<GridPotential<SCFMode>>>& scalarPart,
                                                      std::vector<std::vector<Gradient<GridPotential<SCFMode>>>>& gradientPart,
                                                      unsigned int I) {
  const unsigned int nBlocks = basisFunctionOnGridController->getNBlocks();
  double blockAveThreshold = _kernel->getblockAveThreshold(I);
  bool gga = _kernel->isGGA();

  Eigen::setNbThreads(1);
  if (gga) {
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
      unsigned iThread = omp_get_thread_num();
      auto& blockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
      auto& funcValues = blockData->functionValues;
      auto& gradValues = blockData->derivativeValues;
      unsigned start = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
      unsigned size = blockData->functionValues.rows();
      auto projector = constructProjectionMatrix(blockData->negligible);

      Eigen::MatrixXd funcProj = funcValues * projector;
      Eigen::MatrixXd gradProjx = gradValues->x * projector;
      Eigen::MatrixXd gradProjy = gradValues->y * projector;
      Eigen::MatrixXd gradProjz = gradValues->z * projector;

      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& F = Fxc[iThread][iSet][iGuess];

          auto& S = scalarPart[iSet][iGuess];
          auto& X = gradientPart[iSet][iGuess].x;
          auto& Y = gradientPart[iSet][iGuess].y;
          auto& Z = gradientPart[iSet][iGuess].z;

          double average = 0.0;
          for_spin(S, X, Y, Z) {
            average += S_spin.segment(start, size).cwiseAbs().sum();
            average += X_spin.segment(start, size).cwiseAbs().sum();
            average += Y_spin.segment(start, size).cwiseAbs().sum();
            average += Z_spin.segment(start, size).cwiseAbs().sum();
          };

          if (average / size < blockAveThreshold) {
            continue;
          }

          for_spin(F, S, X, Y, Z) {
            Eigen::MatrixXd tmp = 0.5 * S_spin.segment(start, size).asDiagonal() * funcProj;
            tmp.noalias() += X_spin.segment(start, size).asDiagonal() * gradProjx;
            tmp.noalias() += Y_spin.segment(start, size).asDiagonal() * gradProjy;
            tmp.noalias() += Z_spin.segment(start, size).asDiagonal() * gradProjz;
            Eigen::MatrixXd tmp2 = funcProj.transpose() * tmp;
            F_spin.noalias() += projector * (tmp2 + tmp2.transpose()) * projector.transpose();
          };
        }
      }
    }
  }
  else {
#pragma omp parallel for schedule(dynamic)
    for (unsigned iBlock = 0; iBlock < nBlocks; ++iBlock) {
      unsigned iThread = omp_get_thread_num();
      auto& blockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
      auto& funcValues = blockData->functionValues;
      unsigned size = blockData->functionValues.rows();
      unsigned start = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
      auto projector = constructProjectionMatrix(blockData->negligible);

      Eigen::MatrixXd funcProj = funcValues * projector;

      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& S = scalarPart[iSet][iGuess];
          auto& F = Fxc[iThread][iSet][iGuess];

          for_spin(F, S) {
            double average = S_spin.segment(start, size).cwiseAbs().sum() / size;
            if (average < blockAveThreshold) {
              return;
            }
            Eigen::MatrixXd tmp = S_spin.segment(start, size).asDiagonal() * funcProj;
            Eigen::MatrixXd tmp2 = funcProj.transpose() * tmp;
            F_spin.noalias() += projector * tmp2 * projector.transpose();
          };
        }
      }
    }
  }
  Eigen::setNbThreads(0);
}

template<>
void KernelSigmavector<RESTRICTED>::contractBlock(const unsigned int start, const unsigned int size,
                                                  SpinPolarizedData<RESTRICTED, Eigen::VectorXd>& pbb,
                                                  Gradient<SpinPolarizedData<RESTRICTED, Eigen::VectorXd>>& pnbb,
                                                  GridPotential<RESTRICTED>& scalarPart,
                                                  Gradient<GridPotential<RESTRICTED>>& gradientPart, unsigned int I,
                                                  unsigned int J) {
  auto ppptr = _kernel->getPP(I, J, size, start);
  auto& pp = (*ppptr);
  scalarPart.segment(start, size) += pp.cwiseProduct(pbb);
  if (_kernel->isGGA()) {
    auto pgptr = _kernel->getPG(I, J, size, start);
    auto ggptr = _kernel->getGG(I, J, size, start);
    auto& pg = (*pgptr);
    auto& gg = (*ggptr);
    scalarPart.segment(start, size) += pg.x.cwiseProduct(pnbb.x);
    scalarPart.segment(start, size) += pg.y.cwiseProduct(pnbb.y);
    scalarPart.segment(start, size) += pg.z.cwiseProduct(pnbb.z);
    gradientPart.x.segment(start, size) += pg.x.cwiseProduct(pbb);
    gradientPart.y.segment(start, size) += pg.y.cwiseProduct(pbb);
    gradientPart.z.segment(start, size) += pg.z.cwiseProduct(pbb);
    gradientPart.x.segment(start, size) += gg.xx.cwiseProduct(pnbb.x);
    gradientPart.x.segment(start, size) += gg.xy.cwiseProduct(pnbb.y);
    gradientPart.x.segment(start, size) += gg.xz.cwiseProduct(pnbb.z);
    gradientPart.y.segment(start, size) += gg.xy.cwiseProduct(pnbb.x);
    gradientPart.y.segment(start, size) += gg.yy.cwiseProduct(pnbb.y);
    gradientPart.y.segment(start, size) += gg.yz.cwiseProduct(pnbb.z);
    gradientPart.z.segment(start, size) += gg.xz.cwiseProduct(pnbb.x);
    gradientPart.z.segment(start, size) += gg.yz.cwiseProduct(pnbb.y);
    gradientPart.z.segment(start, size) += gg.zz.cwiseProduct(pnbb.z);
  }
}

template<>
void KernelSigmavector<UNRESTRICTED>::contractBlock(const unsigned int start, const unsigned int size,
                                                    SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>& pbb,
                                                    Gradient<SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>>& pnbb,
                                                    GridPotential<UNRESTRICTED>& scalarPart,
                                                    Gradient<GridPotential<UNRESTRICTED>>& gradientPart, unsigned int I,
                                                    unsigned int J) {
  auto ppptr = _kernel->getPP(I, J, size, start);
  auto& pp = (*ppptr);
  scalarPart.alpha.segment(start, size) += pp.aa.cwiseProduct(pbb.alpha);
  scalarPart.alpha.segment(start, size) += pp.ab.cwiseProduct(pbb.beta);
  scalarPart.beta.segment(start, size) += pp.bb.cwiseProduct(pbb.beta);
  scalarPart.beta.segment(start, size) += pp.ab.cwiseProduct(pbb.alpha);
  if (_kernel->isGGA()) {
    auto pgptr = _kernel->getPG(I, J, size, start);
    auto ggptr = _kernel->getGG(I, J, size, start);
    auto& pg = (*pgptr);
    auto& gg = (*ggptr);
    scalarPart.alpha.segment(start, size) += pg.x.aa.cwiseProduct(pnbb.x.alpha);
    scalarPart.alpha.segment(start, size) += pg.y.aa.cwiseProduct(pnbb.y.alpha);
    scalarPart.alpha.segment(start, size) += pg.z.aa.cwiseProduct(pnbb.z.alpha);
    scalarPart.alpha.segment(start, size) += pg.x.ab.cwiseProduct(pnbb.x.beta);
    scalarPart.alpha.segment(start, size) += pg.y.ab.cwiseProduct(pnbb.y.beta);
    scalarPart.alpha.segment(start, size) += pg.z.ab.cwiseProduct(pnbb.z.beta);
    scalarPart.beta.segment(start, size) += pg.x.ba.cwiseProduct(pnbb.x.alpha);
    scalarPart.beta.segment(start, size) += pg.y.ba.cwiseProduct(pnbb.y.alpha);
    scalarPart.beta.segment(start, size) += pg.z.ba.cwiseProduct(pnbb.z.alpha);
    scalarPart.beta.segment(start, size) += pg.x.bb.cwiseProduct(pnbb.x.beta);
    scalarPart.beta.segment(start, size) += pg.y.bb.cwiseProduct(pnbb.y.beta);
    scalarPart.beta.segment(start, size) += pg.z.bb.cwiseProduct(pnbb.z.beta);
    gradientPart.x.alpha.segment(start, size) += pg.x.aa.cwiseProduct(pbb.alpha);
    gradientPart.y.alpha.segment(start, size) += pg.y.aa.cwiseProduct(pbb.alpha);
    gradientPart.z.alpha.segment(start, size) += pg.z.aa.cwiseProduct(pbb.alpha);
    gradientPart.x.alpha.segment(start, size) += pg.x.ba.cwiseProduct(pbb.beta);
    gradientPart.y.alpha.segment(start, size) += pg.y.ba.cwiseProduct(pbb.beta);
    gradientPart.z.alpha.segment(start, size) += pg.z.ba.cwiseProduct(pbb.beta);
    gradientPart.x.beta.segment(start, size) += pg.x.ab.cwiseProduct(pbb.alpha);
    gradientPart.y.beta.segment(start, size) += pg.y.ab.cwiseProduct(pbb.alpha);
    gradientPart.z.beta.segment(start, size) += pg.z.ab.cwiseProduct(pbb.alpha);
    gradientPart.x.beta.segment(start, size) += pg.x.bb.cwiseProduct(pbb.beta);
    gradientPart.y.beta.segment(start, size) += pg.y.bb.cwiseProduct(pbb.beta);
    gradientPart.z.beta.segment(start, size) += pg.z.bb.cwiseProduct(pbb.beta);
    gradientPart.x.alpha.segment(start, size) += gg.xx.aa.cwiseProduct(pnbb.x.alpha);
    gradientPart.x.alpha.segment(start, size) += gg.xy.aa.cwiseProduct(pnbb.y.alpha);
    gradientPart.x.alpha.segment(start, size) += gg.xz.aa.cwiseProduct(pnbb.z.alpha);
    gradientPart.x.alpha.segment(start, size) += gg.xx.ab.cwiseProduct(pnbb.x.beta);
    gradientPart.x.alpha.segment(start, size) += gg.xy.ab.cwiseProduct(pnbb.y.beta);
    gradientPart.x.alpha.segment(start, size) += gg.xz.ab.cwiseProduct(pnbb.z.beta);
    gradientPart.x.beta.segment(start, size) += gg.xx.bb.cwiseProduct(pnbb.x.beta);
    gradientPart.x.beta.segment(start, size) += gg.xy.bb.cwiseProduct(pnbb.y.beta);
    gradientPart.x.beta.segment(start, size) += gg.xz.bb.cwiseProduct(pnbb.z.beta);
    gradientPart.x.beta.segment(start, size) += gg.xx.ba.cwiseProduct(pnbb.x.alpha);
    gradientPart.x.beta.segment(start, size) += gg.xy.ba.cwiseProduct(pnbb.y.alpha);
    gradientPart.x.beta.segment(start, size) += gg.xz.ba.cwiseProduct(pnbb.z.alpha);
    gradientPart.y.alpha.segment(start, size) += gg.xy.aa.cwiseProduct(pnbb.x.alpha);
    gradientPart.y.alpha.segment(start, size) += gg.yy.aa.cwiseProduct(pnbb.y.alpha);
    gradientPart.y.alpha.segment(start, size) += gg.yz.aa.cwiseProduct(pnbb.z.alpha);
    gradientPart.y.alpha.segment(start, size) += gg.xy.ab.cwiseProduct(pnbb.x.beta);
    gradientPart.y.alpha.segment(start, size) += gg.yy.ab.cwiseProduct(pnbb.y.beta);
    gradientPart.y.alpha.segment(start, size) += gg.yz.ab.cwiseProduct(pnbb.z.beta);
    gradientPart.y.beta.segment(start, size) += gg.xy.bb.cwiseProduct(pnbb.x.beta);
    gradientPart.y.beta.segment(start, size) += gg.yy.bb.cwiseProduct(pnbb.y.beta);
    gradientPart.y.beta.segment(start, size) += gg.yz.bb.cwiseProduct(pnbb.z.beta);
    gradientPart.y.beta.segment(start, size) += gg.xy.ba.cwiseProduct(pnbb.x.alpha);
    gradientPart.y.beta.segment(start, size) += gg.yy.ba.cwiseProduct(pnbb.y.alpha);
    gradientPart.y.beta.segment(start, size) += gg.yz.ba.cwiseProduct(pnbb.z.alpha);
    gradientPart.z.alpha.segment(start, size) += gg.xz.aa.cwiseProduct(pnbb.x.alpha);
    gradientPart.z.alpha.segment(start, size) += gg.yz.aa.cwiseProduct(pnbb.y.alpha);
    gradientPart.z.alpha.segment(start, size) += gg.zz.aa.cwiseProduct(pnbb.z.alpha);
    gradientPart.z.alpha.segment(start, size) += gg.xz.ab.cwiseProduct(pnbb.x.beta);
    gradientPart.z.alpha.segment(start, size) += gg.yz.ab.cwiseProduct(pnbb.y.beta);
    gradientPart.z.alpha.segment(start, size) += gg.zz.ab.cwiseProduct(pnbb.z.beta);
    gradientPart.z.beta.segment(start, size) += gg.xz.bb.cwiseProduct(pnbb.x.beta);
    gradientPart.z.beta.segment(start, size) += gg.yz.bb.cwiseProduct(pnbb.y.beta);
    gradientPart.z.beta.segment(start, size) += gg.zz.bb.cwiseProduct(pnbb.z.beta);
    gradientPart.z.beta.segment(start, size) += gg.xz.ba.cwiseProduct(pnbb.x.alpha);
    gradientPart.z.beta.segment(start, size) += gg.yz.ba.cwiseProduct(pnbb.y.alpha);
    gradientPart.z.beta.segment(start, size) += gg.zz.ba.cwiseProduct(pnbb.z.alpha);
  }
}

template class KernelSigmavector<Options::SCF_MODES::RESTRICTED>;
template class KernelSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
