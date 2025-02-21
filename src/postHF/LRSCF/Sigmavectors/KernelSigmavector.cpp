/**
 * @file KernelSigmavector.cpp
 *
 * @date Dec 07, 2018
 * @author Niklas Niemeyer, Johannes Toelle
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
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
KernelSigmavector<SCFMode>::KernelSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                              std::vector<Eigen::MatrixXd> b, std::shared_ptr<Kernel<SCFMode>> kernel,
                                              std::shared_ptr<Kernel<UNRESTRICTED>> ukernel)
  : Sigmavector<SCFMode>(lrscf, b),
    _kernel(kernel),
    _ukernel(ukernel),
    _blockAveThreshold(lrscf[0]->getLRSCFSettings().grid.blockAveThreshold),
    _isGGA(_kernel->isGGA()) {
  this->contractSupersystemDensity(lrscf);
}

template<Options::SCF_MODES SCFMode>
KernelSigmavector<SCFMode>::KernelSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                              std::shared_ptr<Kernel<SCFMode>> kernel)
  : Sigmavector<SCFMode>(lrscf),
    _kernel(kernel),
    _blockAveThreshold(lrscf[0]->getLRSCFSettings().grid.blockAveThreshold),
    _isGGA(_kernel->isGGA()) {
}

template<Options::SCF_MODES SCFMode>
void KernelSigmavector<SCFMode>::contractSupersystemDensity(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf) {
  Timings::takeTime("LRSCF -   Sigmavector: Kernel");

  // The inter-subsystem kernel is the same for all inter-subsystem blocks of the response matrix.
  // For this reason, the perturbed densities for the supersystem kernel, i.e. the off diagonal
  // blocks need only be evaluated once and not again for each single block of the response matrix.
  // This way, a lot of computation time can be saved for supersystems composed of many subsystems.
  // Note that this simplification cannot be applied in case of a partial response matrix construction
  // (where every single block needs to be evaluated due to the symmetry exploit) or in the case of
  // mixed exact--approximate embedding. Niklas

  bool doIntegrateSupersystem = this->_lrscf.size() > 1;
  doIntegrateSupersystem = doIntegrateSupersystem && !this->_lrscf[0]->getLRSCFSettings().partialResponseConstruction;
  doIntegrateSupersystem = doIntegrateSupersystem && !_kernel->usesMixedEmbedding();
  doIntegrateSupersystem = doIntegrateSupersystem && !this->_lrscf[0]->getLRSCFSettings().noCoupling;
  bool isNotCC2 = this->_lrscf[0]->getLRSCFSettings().method == Options::LR_METHOD::TDDFT ||
                  this->_lrscf[0]->getLRSCFSettings().method == Options::LR_METHOD::TDA;
  doIntegrateSupersystem = doIntegrateSupersystem && isNotCC2;

  if (doIntegrateSupersystem) {
    if (!_kernel) {
      WarningTracker::printWarning("A kernel sigma vector was requested with no kernel present.", true);
    }

    _supersystem_scalar.resize(this->_nSet);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        _supersystem_scalar[iSet].emplace_back(_kernel->getGridController());
      }
    }
    _supersystem_gradient.resize(this->_nSet);
    if (_isGGA) {
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          _supersystem_gradient[iSet].emplace_back(makeGradient<GridPotential<SCFMode>>(_kernel->getGridController()));
        }
      }
    }

    unsigned derivativeLevel = _isGGA ? 1 : 0;

    for (unsigned I = 0; I < lrscf.size(); ++I) {
      auto basisFunctionOnGridControllerI = BasisFunctionOnGridControllerFactory::produce(
          this->_kernel->getBlocksize(I), this->_kernel->getbasFuncRadialThreshold(I), derivativeLevel,
          this->_lrscf[I]->getBasisController(), _kernel->getGridController());

      auto densities = *this->calcP(I);
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& D = densities[iSet][iGuess];
          for_spin(D) {
            D_spin += D_spin.transpose().eval();
          };
        }
      }

      Timings::takeTime("LRSCF -   Kernel: Contraction");
      // Insert 0 and 1 here to make sure that only the supersystem density derivatives are evaluated.
      // The actual numbers do not matter as long as they are different (see Kernel.cpp).
      this->contractKernel(densities, basisFunctionOnGridControllerI, _supersystem_scalar, _supersystem_gradient, 0, 1);
      Timings::timeTaken("LRSCF -   Kernel: Contraction");
    }
  }
  Timings::timeTaken("LRSCF -   Sigmavector: Kernel");
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
KernelSigmavector<SCFMode>::calcF(unsigned I, unsigned J,
                                  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Sigmavector: Kernel");

  // This is needed for coupling pattern where the same system is with coupled
  // and uncoupled vectors. The Kernel has only system information while I and J are
  // related to the lrscf Controller, but the systems behind I and J could be identical
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

  if (!_kernel || (I != J && _supersystem_scalar.size() > 0)) {
    Timings::timeTaken("LRSCF -   Sigmavector: Kernel");
    return nullptr;
  }

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
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

  unsigned derivativeLevel = _isGGA ? 1 : 0;

  auto basisFunctionOnGridControllerI = BasisFunctionOnGridControllerFactory::produce(
      this->_kernel->getBlocksize(I), this->_kernel->getbasFuncRadialThreshold(I), derivativeLevel,
      this->_lrscf[I]->getBasisController(), _kernel->getGridController());
  auto basisFunctionOnGridControllerJ = basisFunctionOnGridControllerI;

  if (I != J) {
    basisFunctionOnGridControllerJ = BasisFunctionOnGridControllerFactory::produce(
        this->_kernel->getBlocksize(J), this->_kernel->getbasFuncRadialThreshold(J), derivativeLevel,
        this->_lrscf[J]->getBasisController(), _kernel->getGridController());
  }

  std::vector<std::vector<GridPotential<SCFMode>>> scalar(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      scalar[iSet].emplace_back(_kernel->getGridController());
    }
  }

  std::vector<std::vector<Gradient<GridPotential<SCFMode>>>> gradient(this->_nSet);
  if (_isGGA) {
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        gradient[iSet].emplace_back(makeGradient<GridPotential<SCFMode>>(_kernel->getGridController()));
      }
    }
  }

  auto& densities = (*densityMatrices);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      auto& D = densities[iSet][iGuess];
      for_spin(D) {
        D_spin += D_spin.transpose().eval();
      };
    }
  }

  Timings::takeTime("LRSCF -   Kernel: Contraction");
  this->contractKernel(densities, basisFunctionOnGridControllerJ, scalar, gradient, I, J);
  Timings::timeTaken("LRSCF -   Kernel: Contraction");

  // Add supersystem contraction to this intra-subsystem contraction.
  if (I == J && _supersystem_scalar.size() > 0) {
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        scalar[iSet][iGuess] += _supersystem_scalar[iSet][iGuess];
        if (_isGGA) {
          gradient[iSet][iGuess].x += _supersystem_gradient[iSet][iGuess].x;
          gradient[iSet][iGuess].y += _supersystem_gradient[iSet][iGuess].y;
          gradient[iSet][iGuess].z += _supersystem_gradient[iSet][iGuess].z;
        }
      }
    }
  }

  Timings::takeTime("LRSCF -   Kernel: Num. Integ.");
  this->numericalIntegration(Fxc, basisFunctionOnGridControllerI, scalar, gradient);
  Timings::timeTaken("LRSCF -   Kernel: Num. Integ.");

  Timings::timeTaken("LRSCF -   Sigmavector: Kernel");

  // Sum over threads.
  Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      auto& F = (*fock)[iSet][iGuess];
      for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
        auto& fxc = Fxc[iThread][iSet][iGuess];
        for_spin(F, fxc) {
          F_spin += fxc_spin;
        };
      }
      for_spin(F) {
        F_spin += F_spin.transpose().eval();
      };
    }
  }
  Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");

  return fock;
}

template<Options::SCF_MODES SCFMode>
void KernelSigmavector<SCFMode>::contractKernel(std::vector<std::vector<MatrixInBasis<SCFMode>>>& dens,
                                                std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                                                std::vector<std::vector<GridPotential<SCFMode>>>& scalar,
                                                std::vector<std::vector<Gradient<GridPotential<SCFMode>>>>& gradient,
                                                unsigned int I, unsigned int J) {
  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
  for (unsigned iBlock = 0; iBlock < basisFunctionOnGridController->getNBlocks(); ++iBlock) {
    auto& blockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
    auto& funcValues = blockData->functionValues;
    auto& gradValues = blockData->derivativeValues;
    auto& average = blockData->averageFunctionValues;

    const unsigned start = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
    const unsigned nb = funcValues.cols();
    const unsigned np = funcValues.rows();
    Eigen::VectorXd contr(np);

    const auto weights = basisFunctionOnGridController->getGridController()->getWeights().segment(start, np);

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        SpinPolarizedData<SCFMode, Eigen::VectorXd> scal(Eigen::VectorXd::Zero(np));
        auto grad = makeGradient<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(Eigen::VectorXd::Zero(np));
        auto& gradx = grad.x;
        auto& grady = grad.y;
        auto& gradz = grad.z;

        auto& D = dens[iSet][iGuess];
        for_spin(D, scal, gradx, grady, gradz) {
          const auto dens_ptr = D_spin.data();
          for (unsigned i = 0; i < nb; ++i) {
            for (unsigned j = 0; j < nb; ++j) {
              double density = dens_ptr[i * nb + j];
              if (std::abs(density * average[i] * average[j]) > _blockAveThreshold) {
                contr = funcValues.col(j) * density;
                scal_spin += funcValues.col(i).cwiseProduct(contr);
                if (_isGGA) {
                  gradx_spin += gradValues->x.col(i).cwiseProduct(contr);
                  grady_spin += gradValues->y.col(i).cwiseProduct(contr);
                  gradz_spin += gradValues->z.col(i).cwiseProduct(contr);
                }
              }
            }
          }
          scal_spin = 0.5 * weights.cwiseProduct(scal_spin);
          if (_isGGA) {
            gradx_spin = weights.cwiseProduct(gradx_spin);
            grady_spin = weights.cwiseProduct(grady_spin);
            gradz_spin = weights.cwiseProduct(gradz_spin);
          }
        };
        this->contractBlock(start, np, scal, grad, scalar[iSet][iGuess], gradient[iSet][iGuess], I, J);
      }
    }
  }
  Eigen::setNbThreads(0);
}

template<Options::SCF_MODES SCFMode>
void KernelSigmavector<SCFMode>::numericalIntegration(std::vector<std::vector<std::vector<MatrixInBasis<SCFMode>>>>& Fxc,
                                                      std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                                                      std::vector<std::vector<GridPotential<SCFMode>>>& scalar,
                                                      std::vector<std::vector<Gradient<GridPotential<SCFMode>>>>& gradient) {
  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < basisFunctionOnGridController->getNBlocks(); ++iBlock) {
    unsigned iThread = omp_get_thread_num();
    auto& blockData = basisFunctionOnGridController->getBlockOnGridData(iBlock);
    auto& funcValues = blockData->functionValues;
    auto& gradValues = blockData->derivativeValues;
    auto& average = blockData->averageFunctionValues;

    const unsigned start = basisFunctionOnGridController->getFirstIndexOfBlock(iBlock);
    const unsigned nb = funcValues.cols();
    const unsigned np = funcValues.rows();
    Eigen::VectorXd contr(np);

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto& fock = Fxc[iThread][iSet][iGuess];
        auto& scal = scalar[iSet][iGuess];
        auto& gradx = gradient[iSet][iGuess].x;
        auto& grady = gradient[iSet][iGuess].y;
        auto& gradz = gradient[iSet][iGuess].z;

        for_spin(fock, scal, gradx, grady, gradz) {
          const auto F = fock_spin.data();
          double scal_sum = scal_spin.segment(start, np).cwiseAbs().sum();
          for (unsigned i = 0; i < nb; ++i) {
            for (unsigned j = 0; j < nb; ++j) {
              if (scal_sum * average[i] * average[j] > _blockAveThreshold) {
                contr = 0.5 * scal_spin.segment(start, np).cwiseProduct(funcValues.col(j));
                if (_isGGA) {
                  contr += gradx_spin.segment(start, np).cwiseProduct(gradValues->x.col(j));
                  contr += grady_spin.segment(start, np).cwiseProduct(gradValues->y.col(j));
                  contr += gradz_spin.segment(start, np).cwiseProduct(gradValues->z.col(j));
                }
                F[i * nb + j] += contr.dot(funcValues.col(i));
              }
            }
          }
        };
      }
    }
  }
  Eigen::setNbThreads(0);
}

template<>
void KernelSigmavector<RESTRICTED>::contractBlock(const unsigned start, const unsigned size,
                                                  SpinPolarizedData<RESTRICTED, Eigen::VectorXd>& p,
                                                  Gradient<SpinPolarizedData<RESTRICTED, Eigen::VectorXd>>& g,
                                                  GridPotential<RESTRICTED>& scalar,
                                                  Gradient<GridPotential<RESTRICTED>>& gradient, unsigned I, unsigned J) {
  // Singlet case.
  if (!_ukernel) {
    auto pp = _kernel->getPP(I, J, size, start);
    scalar.segment(start, size) += pp->cwiseProduct(p);
    if (_isGGA) {
      auto pg = _kernel->getPG(I, J, size, start);
      auto gg = _kernel->getGG(I, J, size, start);
      scalar.segment(start, size) += pg->x.cwiseProduct(g.x);
      scalar.segment(start, size) += pg->y.cwiseProduct(g.y);
      scalar.segment(start, size) += pg->z.cwiseProduct(g.z);
      gradient.x.segment(start, size) += pg->x.cwiseProduct(p);
      gradient.y.segment(start, size) += pg->y.cwiseProduct(p);
      gradient.z.segment(start, size) += pg->z.cwiseProduct(p);
      gradient.x.segment(start, size) += gg->xx.cwiseProduct(g.x);
      gradient.x.segment(start, size) += gg->xy.cwiseProduct(g.y);
      gradient.x.segment(start, size) += gg->xz.cwiseProduct(g.z);
      gradient.y.segment(start, size) += gg->xy.cwiseProduct(g.x);
      gradient.y.segment(start, size) += gg->yy.cwiseProduct(g.y);
      gradient.y.segment(start, size) += gg->yz.cwiseProduct(g.z);
      gradient.z.segment(start, size) += gg->xz.cwiseProduct(g.x);
      gradient.z.segment(start, size) += gg->yz.cwiseProduct(g.y);
      gradient.z.segment(start, size) += gg->zz.cwiseProduct(g.z);
    }
  }
  // Triplet case.
  else {
    auto pp = _ukernel->getPP(I, J, size, start);
    scalar.segment(start, size) += (pp->aa - pp->ab).cwiseProduct(p);
    if (_isGGA) {
      auto pg = _ukernel->getPG(I, J, size, start);
      auto gg = _ukernel->getGG(I, J, size, start);
      scalar.segment(start, size) += (pg->x.aa - pg->x.ab).cwiseProduct(g.x);
      scalar.segment(start, size) += (pg->y.aa - pg->y.ab).cwiseProduct(g.y);
      scalar.segment(start, size) += (pg->z.aa - pg->z.ab).cwiseProduct(g.z);
      gradient.x.segment(start, size) += (pg->x.aa - pg->x.ab).cwiseProduct(p);
      gradient.y.segment(start, size) += (pg->y.aa - pg->y.ab).cwiseProduct(p);
      gradient.z.segment(start, size) += (pg->z.aa - pg->z.ab).cwiseProduct(p);
      gradient.x.segment(start, size) += (gg->xx.aa - gg->xx.ab).cwiseProduct(g.x);
      gradient.x.segment(start, size) += (gg->xy.aa - gg->xy.ab).cwiseProduct(g.y);
      gradient.x.segment(start, size) += (gg->xz.aa - gg->xz.ab).cwiseProduct(g.z);
      gradient.y.segment(start, size) += (gg->xy.aa - gg->xy.ab).cwiseProduct(g.x);
      gradient.y.segment(start, size) += (gg->yy.aa - gg->yy.ab).cwiseProduct(g.y);
      gradient.y.segment(start, size) += (gg->yz.aa - gg->yz.ab).cwiseProduct(g.z);
      gradient.z.segment(start, size) += (gg->xz.aa - gg->xz.ab).cwiseProduct(g.x);
      gradient.z.segment(start, size) += (gg->yz.aa - gg->yz.ab).cwiseProduct(g.y);
      gradient.z.segment(start, size) += (gg->zz.aa - gg->zz.ab).cwiseProduct(g.z);
    }
  }
}

template<>
void KernelSigmavector<UNRESTRICTED>::contractBlock(const unsigned start, const unsigned size,
                                                    SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>& p,
                                                    Gradient<SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>>& g,
                                                    GridPotential<UNRESTRICTED>& scalar,
                                                    Gradient<GridPotential<UNRESTRICTED>>& gradient, unsigned I, unsigned J) {
  auto pp = _kernel->getPP(I, J, size, start);
  scalar.alpha.segment(start, size) += pp->aa.cwiseProduct(p.alpha);
  scalar.alpha.segment(start, size) += pp->ab.cwiseProduct(p.beta);
  scalar.beta.segment(start, size) += pp->bb.cwiseProduct(p.beta);
  scalar.beta.segment(start, size) += pp->ab.cwiseProduct(p.alpha);
  if (_isGGA) {
    auto pg = _kernel->getPG(I, J, size, start);
    auto gg = _kernel->getGG(I, J, size, start);
    scalar.alpha.segment(start, size) += pg->x.aa.cwiseProduct(g.x.alpha);
    scalar.alpha.segment(start, size) += pg->y.aa.cwiseProduct(g.y.alpha);
    scalar.alpha.segment(start, size) += pg->z.aa.cwiseProduct(g.z.alpha);
    scalar.alpha.segment(start, size) += pg->x.ab.cwiseProduct(g.x.beta);
    scalar.alpha.segment(start, size) += pg->y.ab.cwiseProduct(g.y.beta);
    scalar.alpha.segment(start, size) += pg->z.ab.cwiseProduct(g.z.beta);
    scalar.beta.segment(start, size) += pg->x.ba.cwiseProduct(g.x.alpha);
    scalar.beta.segment(start, size) += pg->y.ba.cwiseProduct(g.y.alpha);
    scalar.beta.segment(start, size) += pg->z.ba.cwiseProduct(g.z.alpha);
    scalar.beta.segment(start, size) += pg->x.bb.cwiseProduct(g.x.beta);
    scalar.beta.segment(start, size) += pg->y.bb.cwiseProduct(g.y.beta);
    scalar.beta.segment(start, size) += pg->z.bb.cwiseProduct(g.z.beta);
    gradient.x.alpha.segment(start, size) += pg->x.aa.cwiseProduct(p.alpha);
    gradient.y.alpha.segment(start, size) += pg->y.aa.cwiseProduct(p.alpha);
    gradient.z.alpha.segment(start, size) += pg->z.aa.cwiseProduct(p.alpha);
    gradient.x.alpha.segment(start, size) += pg->x.ba.cwiseProduct(p.beta);
    gradient.y.alpha.segment(start, size) += pg->y.ba.cwiseProduct(p.beta);
    gradient.z.alpha.segment(start, size) += pg->z.ba.cwiseProduct(p.beta);
    gradient.x.beta.segment(start, size) += pg->x.ab.cwiseProduct(p.alpha);
    gradient.y.beta.segment(start, size) += pg->y.ab.cwiseProduct(p.alpha);
    gradient.z.beta.segment(start, size) += pg->z.ab.cwiseProduct(p.alpha);
    gradient.x.beta.segment(start, size) += pg->x.bb.cwiseProduct(p.beta);
    gradient.y.beta.segment(start, size) += pg->y.bb.cwiseProduct(p.beta);
    gradient.z.beta.segment(start, size) += pg->z.bb.cwiseProduct(p.beta);
    gradient.x.alpha.segment(start, size) += gg->xx.aa.cwiseProduct(g.x.alpha);
    gradient.x.alpha.segment(start, size) += gg->xy.aa.cwiseProduct(g.y.alpha);
    gradient.x.alpha.segment(start, size) += gg->xz.aa.cwiseProduct(g.z.alpha);
    gradient.x.alpha.segment(start, size) += gg->xx.ab.cwiseProduct(g.x.beta);
    gradient.x.alpha.segment(start, size) += gg->xy.ab.cwiseProduct(g.y.beta);
    gradient.x.alpha.segment(start, size) += gg->xz.ab.cwiseProduct(g.z.beta);
    gradient.x.beta.segment(start, size) += gg->xx.bb.cwiseProduct(g.x.beta);
    gradient.x.beta.segment(start, size) += gg->xy.bb.cwiseProduct(g.y.beta);
    gradient.x.beta.segment(start, size) += gg->xz.bb.cwiseProduct(g.z.beta);
    gradient.x.beta.segment(start, size) += gg->xx.ba.cwiseProduct(g.x.alpha);
    gradient.x.beta.segment(start, size) += gg->xy.ba.cwiseProduct(g.y.alpha);
    gradient.x.beta.segment(start, size) += gg->xz.ba.cwiseProduct(g.z.alpha);
    gradient.y.alpha.segment(start, size) += gg->xy.aa.cwiseProduct(g.x.alpha);
    gradient.y.alpha.segment(start, size) += gg->yy.aa.cwiseProduct(g.y.alpha);
    gradient.y.alpha.segment(start, size) += gg->yz.aa.cwiseProduct(g.z.alpha);
    gradient.y.alpha.segment(start, size) += gg->xy.ab.cwiseProduct(g.x.beta);
    gradient.y.alpha.segment(start, size) += gg->yy.ab.cwiseProduct(g.y.beta);
    gradient.y.alpha.segment(start, size) += gg->yz.ab.cwiseProduct(g.z.beta);
    gradient.y.beta.segment(start, size) += gg->xy.bb.cwiseProduct(g.x.beta);
    gradient.y.beta.segment(start, size) += gg->yy.bb.cwiseProduct(g.y.beta);
    gradient.y.beta.segment(start, size) += gg->yz.bb.cwiseProduct(g.z.beta);
    gradient.y.beta.segment(start, size) += gg->xy.ba.cwiseProduct(g.x.alpha);
    gradient.y.beta.segment(start, size) += gg->yy.ba.cwiseProduct(g.y.alpha);
    gradient.y.beta.segment(start, size) += gg->yz.ba.cwiseProduct(g.z.alpha);
    gradient.z.alpha.segment(start, size) += gg->xz.aa.cwiseProduct(g.x.alpha);
    gradient.z.alpha.segment(start, size) += gg->yz.aa.cwiseProduct(g.y.alpha);
    gradient.z.alpha.segment(start, size) += gg->zz.aa.cwiseProduct(g.z.alpha);
    gradient.z.alpha.segment(start, size) += gg->xz.ab.cwiseProduct(g.x.beta);
    gradient.z.alpha.segment(start, size) += gg->yz.ab.cwiseProduct(g.y.beta);
    gradient.z.alpha.segment(start, size) += gg->zz.ab.cwiseProduct(g.z.beta);
    gradient.z.beta.segment(start, size) += gg->xz.bb.cwiseProduct(g.x.beta);
    gradient.z.beta.segment(start, size) += gg->yz.bb.cwiseProduct(g.y.beta);
    gradient.z.beta.segment(start, size) += gg->zz.bb.cwiseProduct(g.z.beta);
    gradient.z.beta.segment(start, size) += gg->xz.ba.cwiseProduct(g.x.alpha);
    gradient.z.beta.segment(start, size) += gg->yz.ba.cwiseProduct(g.y.alpha);
    gradient.z.beta.segment(start, size) += gg->zz.ba.cwiseProduct(g.z.alpha);
  }
}

template class KernelSigmavector<Options::SCF_MODES::RESTRICTED>;
template class KernelSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
