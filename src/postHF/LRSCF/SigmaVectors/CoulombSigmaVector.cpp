/**
 * @file CoulombSigmaVector.cpp
 *
 * @date Dec 07, 2018
 * @author Michael Boeckers, Johannes Toelle
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
#include "postHF/LRSCF/SigmaVectors/CoulombSigmaVector.h"

/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/CustomBasisController.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "integrals/looper/CoulombInteractionIntLooper.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "misc/Timing.h"
#include "potentials/CoulombPotential.h"
#include "potentials/ERIPotential.h"
#include "settings/Settings.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CoulombSigmaVector<SCFMode>::CoulombSigmaVector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                std::vector<Eigen::MatrixXd> b, const double densityScreeningThreshold)
  : SigmaVector<SCFMode>(lrscf, b, densityScreeningThreshold) {
}

template<Options::SCF_MODES SCFMode>
CoulombSigmaVector<SCFMode>::CoulombSigmaVector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                const double densityScreeningThreshold)
  : SigmaVector<SCFMode>(lrscf, densityScreeningThreshold) {
}

/*
Separates the Coulomb evaluation for RESTRICTED and UNRESTRICTED because of performance reasons.
The following calculation is performed: \f$ \tilde{F}_{ij} = \sum_{k l} \tilde{P}_{kl} (ij|kl) \f$
For LRSCF problems, the (pseudo) density matrix is non symmetric. Every integral (ij|kl) is needed for two
Fock-matrix elements, F_{ij} and F_{kl}. Since the resulting Coulomb Fock-like matrix
is symmetric, we can add P_{kl} to P_{lk} (and P_{ij} to P_{ji} to save some time.
For the diagonal elements of P and F, we then need to divide by a factor of two, which is
done by the the factor perm.*/

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>
CoulombSigmaVector<Options::SCF_MODES::RESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Fock-like matrix: J");
  // Number of AOs in subsystem I and J
  const unsigned int nBFs_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
  const unsigned int nBFs_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();
  // Set dimensions for Fock like matrices
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> fock(
      new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>(this->_nSet));
  for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }
  // Evaluation of the coulomb contribution with the RI approximation (based on CoulombPotential.h and
  // CoulombInteractionPotential.h)
  if (this->_lrscf[I]->getSysSettings().dft.densityFitting == Options::DENS_FITS::RI) {
    // Evaluates the intra-subsystem coulomb contribution
    if (I == J) {
      auto sys = this->_lrscf[I]->getSys();
      auto basisController = sys->getAtomCenteredBasisController();
      auto auxBasisController = sys->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
      auto& factory = RI_J_IntegralControllerFactory::getInstance();
      auto ri_j_IntController = factory.produce(basisController, auxBasisController);
      // Number of aux functions
      const unsigned int nAuxFunctions = ri_j_IntController->getAuxBasisController()->getNBasisFunctions();
      // Simple Prescreening based on max density
      double maxDens = 0.0;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          maxDens = std::max(maxDens, (*densityMatrices)[iSet][iGuess].array().abs().maxCoeff());
        }
      }
      // Use smallest prescreening threshold of subsystem I and J
      double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                              this->_lrscf[J]->getSysSettings().basis.integralThreshold);
      // Prescreening function
      auto prescreeningFunc = [&](unsigned int i, unsigned int j, unsigned int nI, unsigned int nJ, double schwartz) {
        (void)i;
        (void)j;
        (void)nI;
        (void)nJ;
        if (maxDens * schwartz < prescreeningThreshold)
          return true;
        return false;
      };

#ifdef _OPENMP
      std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMat[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMat[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctions, omp_get_max_threads());
          sumMat[iSet][iGuess].setZero();
        }
      }
#else
      std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMat[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMat[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctions, 1);
          sumMat[iSet][iGuess].setZero();
        }
      }
#endif

      auto const loopEvalFunction = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                        const double& integral, const unsigned int threadId) {
        const unsigned int ij = i * nBFs_I + j;
        const unsigned int ji = j * nBFs_I + i;
        // Perm
        double perm = 1.0;
        perm *= (i == j) ? 0.5 : 1.0;
        const double coul = perm * integral;
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pp = (*densityMatrices)[iSet][iGuess];
            double coul1 = 0;
            coul1 += *(pp.data() + ij);
            coul1 += *(pp.data() + ji);
            coul1 *= coul;
            sumMat[iSet][iGuess](K, threadId) += coul1;
          }
        }
      };

      ri_j_IntController->loopOver3CInts(loopEvalFunction, prescreeningFunc);

      std::vector<std::vector<Eigen::VectorXd>> coefficients(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        coefficients[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          coefficients[iSet][iGuess] = ri_j_IntController->getInverseM() * sumMat[iSet][iGuess].rowwise().sum();
        }
      }

      // Thread safety
      std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>*> f;
#ifdef _OPENMP
      const unsigned int nThreads = omp_get_max_threads();
      f.push_back(&(*fock));
      for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
        f.push_back(new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>);
        (*f[iThread]).resize(this->_nSet);
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
          }
        }
      }
#else
      f.push_back(&(*fock));
#endif

      auto const loopEvalFunction2 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         const double& integral, const unsigned int threadId) {
        (void)threadId; // no warning

        const unsigned int ij = i * nBFs_I + j;
        const unsigned int ji = j * nBFs_I + i;

        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pf = (*f[threadId])[iSet][iGuess];
            double coul = 0.0;
            if (i == j) {
              coul = integral * coefficients[iSet][iGuess](K);
              *(pf.data() + ij) += coul;
            }
            else {
              coul = integral * coefficients[iSet][iGuess](K);
              *(pf.data() + ij) += coul;
              *(pf.data() + ji) += coul;
            }
          }
        }
      };

      ri_j_IntController->loopOver3CInts(loopEvalFunction2, prescreeningFunc);

#ifdef _OPENMP
      for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
          }
        }
        delete f[iThread];
      }
#endif

      // Evaluates the inter-subsystem coulomb contribution (I != J)
    }
    else {
      auto& factory = RI_J_IntegralControllerFactory::getInstance();
      // System I
      auto sysI = this->_lrscf[I]->getSys();
      auto basisControllerI = sysI->getAtomCenteredBasisController();
      auto auxBasisControllerI = sysI->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
      auto ri_j_IntControllerI = factory.produce(basisControllerI, auxBasisControllerI);
      const unsigned int nAuxFunctionsI = ri_j_IntControllerI->getAuxBasisController()->getNBasisFunctions();
      // System J
      auto sysJ = this->_lrscf[J]->getSys();
      auto basisControllerJ = sysJ->getAtomCenteredBasisController();
      auto auxBasisControllerJ = sysJ->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
      auto ri_j_IntControllerJ = factory.produce(basisControllerJ, auxBasisControllerJ);
      const unsigned int nAuxFunctionsJ = ri_j_IntControllerJ->getAuxBasisController()->getNBasisFunctions();
      // super system
      Basis superBas;
      Basis superAuxBas;
      superBas.reserve(basisControllerI->getBasis().size() + basisControllerJ->getBasis().size());
      superAuxBas.reserve(auxBasisControllerI->getBasis().size() + auxBasisControllerJ->getBasis().size());
      superBas.insert(superBas.end(), basisControllerI->getBasis().begin(), basisControllerI->getBasis().end());
      superBas.insert(superBas.end(), basisControllerJ->getBasis().begin(), basisControllerJ->getBasis().end());
      superAuxBas.insert(superAuxBas.end(), auxBasisControllerI->getBasis().begin(), auxBasisControllerI->getBasis().end());
      superAuxBas.insert(superAuxBas.end(), auxBasisControllerJ->getBasis().begin(), auxBasisControllerJ->getBasis().end());
      // SupersystemBasisController
      auto superSystemBasisController =
          std::shared_ptr<BasisController>(new CustomBasisController(superBas, "tmpSuperBas"));
      auto superSystemAuxBasisController =
          std::shared_ptr<BasisController>(new CustomBasisController(superAuxBas, "tmpSuperAuxBas"));

#ifdef _OPENMP
      std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMat[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMat[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctionsI, omp_get_max_threads());
          sumMat[iSet][iGuess].setZero();
        }
      }
#else
      std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMat[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMat[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctionsI, 1);
          sumMat[iSet][iGuess].setZero();
        }
      }
#endif

#ifdef _OPENMP
      std::vector<std::vector<Eigen::MatrixXd>> sumMatEnv(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMatEnv[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMatEnv[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctionsJ, omp_get_max_threads());
          sumMatEnv[iSet][iGuess].setZero();
        }
      }
#else
      std::vector<std::vector<Eigen::MatrixXd>> sumMatEnv(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMatEnv[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMatEnv[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctionsJ, 1);
          sumMatEnv[iSet][iGuess].setZero();
        }
      }
#endif
      // Supersystem ri_j controller
      auto ri_j_IntController =
          RI_J_IntegralControllerFactory::getInstance().produce(superSystemBasisController, superSystemAuxBasisController);
      const auto inverseM = ri_j_IntController->getInverseM();

      TwoElecThreeCenterIntLooper looper1(libint2::Operator::coulomb, 0, basisControllerJ, auxBasisControllerI, 1E-10);

      auto const loopEvalFunction1 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadId) {
        const unsigned int ij = i * nBFs_J + j;
        const unsigned int ji = j * nBFs_J + i;
        // Perm
        double perm = 1.0;
        perm *= (i == j) ? 0.5 : 1.0;
        const double coul = perm * integral[0];
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pp = (*densityMatrices)[iSet][iGuess];
            double coul1 = 0;
            coul1 += *(pp.data() + ij);
            coul1 += *(pp.data() + ji);
            coul1 *= coul;
            sumMat[iSet][iGuess](K, threadId) += coul1;
          }
        }
      };

      looper1.loop(loopEvalFunction1);

      TwoElecThreeCenterIntLooper looper2(libint2::Operator::coulomb, 0, basisControllerJ, auxBasisControllerJ, 1E-10);

      auto const loopEvalFunction2 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadId) {
        const unsigned int ij = i * nBFs_J + j;
        const unsigned int ji = j * nBFs_J + i;
        // Perm
        double perm = 1.0;
        perm *= (i == j) ? 0.5 : 1.0;
        const double coul = perm * integral[0];
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pp = (*densityMatrices)[iSet][iGuess];
            // auto sumptr = sumMatEnv[iSet][iGuess].data();
            double coul1 = 0;
            coul1 += *(pp.data() + ij);
            coul1 += *(pp.data() + ji);
            coul1 *= coul;
            sumMatEnv[iSet][iGuess](K, threadId) += coul1;
          }
        }
      };

      looper2.loop(loopEvalFunction2);

      std::vector<std::vector<Eigen::VectorXd>> coefficients(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        coefficients[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          Eigen::VectorXd sumPMuNu_DMuNuAct = sumMat[iSet][iGuess].rowwise().sum();
          sumMat[iSet][iGuess].resize(1, 1);
          Eigen::VectorXd sumPMuNu_DMuNuEnv = sumMatEnv[iSet][iGuess].rowwise().sum();
          sumMatEnv[iSet][iGuess].resize(1, 1);
          Eigen::VectorXd sumPMuNu_DMuNu(sumPMuNu_DMuNuAct.rows() + sumPMuNu_DMuNuEnv.rows());
          sumPMuNu_DMuNu << sumPMuNu_DMuNuAct, sumPMuNu_DMuNuEnv;
          coefficients[iSet][iGuess] = inverseM * sumPMuNu_DMuNu.eval();
        }
      }

      // Thread safety
      std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>*> f;
#ifdef _OPENMP
      const unsigned int nThreads = omp_get_max_threads();
      f.push_back(&(*fock));
      for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
        f.push_back(new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>);
        (*f[iThread]).resize(this->_nSet);
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
          }
        }
      }
#else
      f.push_back(&(*fock));
#endif

      TwoElecThreeCenterIntLooper looper3(libint2::Operator::coulomb, 0, basisControllerI, auxBasisControllerI, 1E-10);

      auto const loopEvalFunction3 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadID) {
        const unsigned int ij = i * nBFs_I + j;
        const unsigned int ji = j * nBFs_I + i;

        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pf = (*f[threadID])[iSet][iGuess];
            double coul = 0.0;
            coul = integral[0] * coefficients[iSet][iGuess](K);
            *(pf.data() + ij) += coul;
            if (i != j) {
              coul = integral[0] * coefficients[iSet][iGuess](K);
              *(pf.data() + ji) += coul;
            }
          }
        }
      };

      looper3.loop(loopEvalFunction3);

      TwoElecThreeCenterIntLooper looper4(libint2::Operator::coulomb, 0, basisControllerI, auxBasisControllerJ, 1E-10);
      auto const loopEvalFunction4 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadID) {
        const unsigned int ij = i * nBFs_I + j;
        const unsigned int ji = j * nBFs_I + i;

        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pf = (*f[threadID])[iSet][iGuess];
            double coul = 0.0;
            coul = integral[0] * coefficients[iSet][iGuess](K + nAuxFunctionsI);
            *(pf.data() + ij) += coul;
            if (i != j) {
              coul = integral[0] * coefficients[iSet][iGuess](K + nAuxFunctionsI);
              *(pf.data() + ji) += coul;
            }
          }
        }
      };
      looper4.loop(loopEvalFunction4);

#ifdef _OPENMP
      // const unsigned int nThreads = omp_get_max_threads();
      for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
          }
        }
        delete f[iThread];
      }
#endif
    }

    // If no RI is used:
  }
  else {
    // Thread safety
    std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>*> f;
#ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    f.push_back(&(*fock));
    for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
      f.push_back(new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>);
      (*f[iThread]).resize(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
        }
      }
    }
#else
    f.push_back(&(*fock));
#endif

    // Function to calculate Coulomb pseudo-Fock matrix (for I=J)
    //
    // F_{(ij} = \sum_{kl} (ij|kl) P_{kl}
    //
    // Here, each integral is multiplied with a density matrix element. For LRSCF problems,
    // the (pseudo) density matrix is non symmetric. Every integral (ij|kl) is needed for two
    // Fock-matrix elements, F_{ij} and F_{kl}. Since the resulting Coulomb Fock-like matrix
    // is symmetric, we can add P_{kl} to P_{lk} (and P_{ij} to P_{ji} to save some time.
    // For the diagonal elements of P and F, we then need to divide by a factor of two, which is
    // done by the the factor perm.
    auto distribute_II = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, const double integral,
                             unsigned int threadId) {
      const unsigned int kl = k * nBFs_I + l;
      const unsigned int lk = l * nBFs_I + k;
      const unsigned int ij = i * nBFs_I + j;
      const unsigned int ji = j * nBFs_I + i;
      double perm = 1.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;
      perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
      const double coul = perm * integral;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& pp = (*densityMatrices)[iSet][iGuess];
          auto& pf = (*f[threadId])[iSet][iGuess];
          double coul1 = 0;
          double coul2 = 0;
          coul1 += *(pp.data() + kl);
          coul1 += *(pp.data() + lk);
          coul2 += *(pp.data() + ij);
          coul2 += *(pp.data() + ji);
          coul1 *= coul;
          coul2 *= coul;
          *(pf.data() + ij) += coul1;
          *(pf.data() + ji) += coul1;
          *(pf.data() + kl) += coul2;
          *(pf.data() + lk) += coul2;
        }
      }
    };

    // Function to calculate Coulomb pseudo-Fock matrix (for I!=J)
    auto distribute_IJ = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, const double integral,
                             unsigned int threadId) {
      const unsigned int kl = k * nBFs_J + l;
      const unsigned int lk = l * nBFs_J + k;
      const unsigned int ij = i * nBFs_I + j;
      const unsigned int ji = j * nBFs_I + i;
      // Permutations
      double perm = 1.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;
      double coul = perm * integral;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& pp = (*densityMatrices)[iSet][iGuess];
          auto& pf = (*f[threadId])[iSet][iGuess];
          double coulJ = 0;
          coulJ += *(pp.data() + kl);
          coulJ += *(pp.data() + lk);
          coulJ *= coul;
          *(pf.data() + ij) += coulJ;
          *(pf.data() + ji) += coulJ;
        }
      }
    };

    double maxDens = 0.0;
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        maxDens = std::max(maxDens, (*densityMatrices)[iSet][iGuess].array().abs().maxCoeff());
      }
    }
    // Use smalles prescreening threshold of subsystem I and J
    double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                            this->_lrscf[J]->getSysSettings().basis.integralThreshold);

    auto prescreeningFunc = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, unsigned int nI,
                                unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
      (void)i;
      (void)j;
      (void)k;
      (void)l;
      (void)nI;
      (void)nJ;
      (void)nK;
      (void)nL;
      if (maxDens * schwartz < prescreeningThreshold)
        return true;
      return false;
    };

    // Calculate pseudo Fock matrices
    if (I == J) {
      TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                        prescreeningThreshold);
      looper.loopNoDerivative(distribute_II, prescreeningFunc);
    }
    else if (I != J) {
      CoulombInteractionIntLooper looper(libint2::Operator::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                         this->_lrscf[J]->getBasisController(), prescreeningThreshold);
      looper.loopNoDerivative(distribute_IJ, prescreeningFunc);
    }
    else {
      assert(false);
    }
#ifdef _OPENMP
    // const unsigned int nThreads = omp_get_max_threads();
    for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
        }
      }
      delete f[iThread];
    }
#endif
  }
  Timings::timeTaken("LRSCF -   Fock-like matrix: J");
  return fock;
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>
CoulombSigmaVector<Options::SCF_MODES::UNRESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Fock-like matrix: J");

  // Number of AOs in subsystem I and J
  const unsigned int nBFs_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
  const unsigned int nBFs_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();
  // Set dimensions for Fock like matrices
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> fock(
      new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>(this->_nSet));
  for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // TODO: Add the spin contributions in the begining with .total() -> should be FASTER!
  if (this->_lrscf[I]->getSysSettings().dft.densityFitting == Options::DENS_FITS::RI) {
    if (I == J) {
      auto sys = this->_lrscf[I]->getSys();
      auto basisController = sys->getAtomCenteredBasisController();
      auto auxBasisController = sys->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
      auto& factory = RI_J_IntegralControllerFactory::getInstance();
      auto ri_j_IntController = factory.produce(basisController, auxBasisController);

      // Function to calculate Coulomb pseudo-Fock matrix (for I=J)
      //
      // F_{(ij} = \sum_{kl} (ij|kl) P_{kl}
      //
      // Here, each integral is multiplied with a density matrix element. For LRSCF problems,
      // the (pseudo) density matrix is non symmetric. Every integral (ij|kl) is needed for two
      // Fock-matrix elements, F_{ij} and F_{kl}. Since the resulting Coulomb Fock-like matrix
      // is symmetric, we can add P_{kl} to P_{lk} (and P_{ij} to P_{ji} to save some time.
      // For the diagonal elements of P and F, we then need to divide by a factor of two, which is
      // done by the the factor perm.
      const unsigned int nAuxFunctions = ri_j_IntController->getAuxBasisController()->getNBasisFunctions();

      double maxDens = 0.0;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& P = (*densityMatrices)[iSet][iGuess];
          for_spin(P) {
            maxDens = std::max(maxDens, P_spin.array().abs().maxCoeff());
          };
        }
      }
      // Use smalles prescreening threshold of subsystem I and J
      double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                              this->_lrscf[J]->getSysSettings().basis.integralThreshold);

      auto prescreeningFunc = [&](unsigned int i, unsigned int j, unsigned int nI, unsigned int nJ, double schwartz) {
        (void)i;
        (void)j;
        (void)nI;
        (void)nJ;
        if (maxDens * schwartz < prescreeningThreshold)
          return true;
        return false;
      };

#ifdef _OPENMP
      std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMat[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMat[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctions, omp_get_max_threads());
          sumMat[iSet][iGuess].setZero();
        }
      }
#else
      std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMat[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMat[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctions, 1);
          sumMat[iSet][iGuess].setZero();
        }
      }
#endif

      auto const loopEvalFunction = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                        const double& integral, const unsigned int threadId) {
        const unsigned int ij = i * nBFs_I + j;
        const unsigned int ji = j * nBFs_I + i;
        // Perm
        double perm = 1.0;
        perm *= (i == j) ? 0.5 : 1.0;
        const double coul = perm * integral;
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pp = (*densityMatrices)[iSet][iGuess];
            double coul1 = 0;
            coul1 += *(pp.alpha.data() + ij);
            coul1 += *(pp.alpha.data() + ji);

            coul1 += *(pp.beta.data() + ij);
            coul1 += *(pp.beta.data() + ji);

            coul1 *= coul;
            sumMat[iSet][iGuess](K, threadId) += coul1;
          }
        }
      };

      ri_j_IntController->loopOver3CInts(loopEvalFunction, prescreeningFunc);

      std::vector<std::vector<Eigen::VectorXd>> coefficients(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        coefficients[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          coefficients[iSet][iGuess] = ri_j_IntController->getInverseM() * sumMat[iSet][iGuess].rowwise().sum();
        }
      }

      // Thread safety
      std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>*> f;
#ifdef _OPENMP
      const unsigned int nThreads = omp_get_max_threads();
      f.push_back(&(*fock));
      for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
        f.push_back(new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>);
        (*f[iThread]).resize(this->_nSet);
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
          }
        }
      }
#else
      f.push_back(&(*fock));
#endif

      auto const loopEvalFunction2 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         const double& integral, const unsigned int threadId) {
        (void)threadId; // no warning

        const unsigned int ij = i * nBFs_I + j;
        const unsigned int ji = j * nBFs_I + i;

        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pf = (*f[threadId])[iSet][iGuess];
            double coul = 0.0;
            if (i == j) {
              coul = integral * coefficients[iSet][iGuess](K);
              *(pf.alpha.data() + ij) += coul;
              *(pf.beta.data() + ij) += coul;
            }
            else {
              coul = integral * coefficients[iSet][iGuess](K);
              *(pf.alpha.data() + ij) += coul;
              *(pf.alpha.data() + ji) += coul;

              *(pf.beta.data() + ij) += coul;
              *(pf.beta.data() + ji) += coul;
            }
          }
        }
      };

      ri_j_IntController->loopOver3CInts(loopEvalFunction2, prescreeningFunc);

#ifdef _OPENMP
      // const unsigned int nThreads = omp_get_max_threads();
      for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
          }
        }
        delete f[iThread];
      }
#endif

      // if (I!=J) in the RI case
    }
    else {
      auto& factory = RI_J_IntegralControllerFactory::getInstance();
      // System I
      auto sysI = this->_lrscf[I]->getSys();
      auto basisControllerI = sysI->getAtomCenteredBasisController();
      auto auxBasisControllerI = sysI->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
      auto ri_j_IntControllerI = factory.produce(basisControllerI, auxBasisControllerI);
      const unsigned int nAuxFunctionsI = ri_j_IntControllerI->getAuxBasisController()->getNBasisFunctions();
      // System J
      auto sysJ = this->_lrscf[J]->getSys();
      auto basisControllerJ = sysJ->getAtomCenteredBasisController();
      auto auxBasisControllerJ = sysJ->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
      auto ri_j_IntControllerJ = factory.produce(basisControllerJ, auxBasisControllerJ);
      const unsigned int nAuxFunctionsJ = ri_j_IntControllerJ->getAuxBasisController()->getNBasisFunctions();
      // super system
      Basis superBas;
      Basis superAuxBas;
      superBas.reserve(basisControllerI->getBasis().size() + basisControllerJ->getBasis().size());
      superAuxBas.reserve(auxBasisControllerI->getBasis().size() + auxBasisControllerJ->getBasis().size());
      superBas.insert(superBas.end(), basisControllerI->getBasis().begin(), basisControllerI->getBasis().end());
      superBas.insert(superBas.end(), basisControllerJ->getBasis().begin(), basisControllerJ->getBasis().end());
      superAuxBas.insert(superAuxBas.end(), auxBasisControllerI->getBasis().begin(), auxBasisControllerI->getBasis().end());
      superAuxBas.insert(superAuxBas.end(), auxBasisControllerJ->getBasis().begin(), auxBasisControllerJ->getBasis().end());
      // SupersystemBasisController
      auto superSystemBasisController =
          std::shared_ptr<BasisController>(new CustomBasisController(superBas, "tmpSuperBas"));
      auto superSystemAuxBasisController =
          std::shared_ptr<BasisController>(new CustomBasisController(superAuxBas, "tmpSuperAuxBas"));

#ifdef _OPENMP
      std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMat[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMat[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctionsI, omp_get_max_threads());
          sumMat[iSet][iGuess].setZero();
        }
      }
#else
      std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMat[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMat[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctionsI, 1);
          sumMat[iSet][iGuess].setZero();
        }
      }
#endif

#ifdef _OPENMP
      std::vector<std::vector<Eigen::MatrixXd>> sumMatEnv(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMatEnv[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMatEnv[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctionsJ, omp_get_max_threads());
          sumMatEnv[iSet][iGuess].setZero();
        }
      }
#else
      std::vector<std::vector<Eigen::MatrixXd>> sumMatEnv(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        sumMatEnv[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          sumMatEnv[iSet][iGuess] = Eigen::MatrixXd(nAuxFunctionsJ, 1);
          sumMatEnv[iSet][iGuess].setZero();
        }
      }
#endif
      // Supersystem ri_j controller
      auto ri_j_IntController =
          RI_J_IntegralControllerFactory::getInstance().produce(superSystemBasisController, superSystemAuxBasisController);
      const auto inverseM = ri_j_IntController->getInverseM();

      TwoElecThreeCenterIntLooper looper1(libint2::Operator::coulomb, 0, basisControllerJ, auxBasisControllerI, 1E-10);

      auto const loopEvalFunction1 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadId) {
        const unsigned int ij = i * nBFs_J + j;
        const unsigned int ji = j * nBFs_J + i;
        // Perm
        double perm = 1.0;
        perm *= (i == j) ? 0.5 : 1.0;
        const double coul = perm * integral[0];
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pp = (*densityMatrices)[iSet][iGuess];
            double coul1 = 0;
            coul1 += *(pp.alpha.data() + ij);
            coul1 += *(pp.alpha.data() + ji);
            coul1 += *(pp.beta.data() + ij);
            coul1 += *(pp.beta.data() + ji);
            coul1 *= coul;
            sumMat[iSet][iGuess](K, threadId) += coul1;
          }
        }
      };

      looper1.loop(loopEvalFunction1);

      TwoElecThreeCenterIntLooper looper2(libint2::Operator::coulomb, 0, basisControllerJ, auxBasisControllerJ, 1E-10);

      auto const loopEvalFunction2 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadId) {
        const unsigned int ij = i * nBFs_J + j;
        const unsigned int ji = j * nBFs_J + i;
        // Perm
        double perm = 1.0;
        perm *= (i == j) ? 0.5 : 1.0;
        const double coul = perm * integral[0];
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pp = (*densityMatrices)[iSet][iGuess];
            // auto sumptr = sumMatEnv[iSet][iGuess].data();
            double coul1 = 0;
            coul1 += *(pp.alpha.data() + ij);
            coul1 += *(pp.alpha.data() + ji);
            coul1 += *(pp.beta.data() + ij);
            coul1 += *(pp.beta.data() + ji);
            coul1 *= coul;
            sumMatEnv[iSet][iGuess](K, threadId) += coul1;
          }
        }
      };

      looper2.loop(loopEvalFunction2);

      std::vector<std::vector<Eigen::VectorXd>> coefficients(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        coefficients[iSet].resize(this->_nGuess);
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          Eigen::VectorXd sumPMuNu_DMuNuAct = sumMat[iSet][iGuess].rowwise().sum();
          sumMat[iSet][iGuess].resize(1, 1);
          Eigen::VectorXd sumPMuNu_DMuNuEnv = sumMatEnv[iSet][iGuess].rowwise().sum();
          sumMatEnv[iSet][iGuess].resize(1, 1);
          Eigen::VectorXd sumPMuNu_DMuNu(sumPMuNu_DMuNuAct.rows() + sumPMuNu_DMuNuEnv.rows());
          sumPMuNu_DMuNu << sumPMuNu_DMuNuAct, sumPMuNu_DMuNuEnv;
          coefficients[iSet][iGuess] = inverseM * sumPMuNu_DMuNu.eval();
        }
      }

      // Thread safety
      std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>*> f;
#ifdef _OPENMP
      const unsigned int nThreads = omp_get_max_threads();
      f.push_back(&(*fock));
      for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
        f.push_back(new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>);
        (*f[iThread]).resize(this->_nSet);
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
          }
        }
      }
#else
      f.push_back(&(*fock));
#endif

      TwoElecThreeCenterIntLooper looper3(libint2::Operator::coulomb, 0, basisControllerI, auxBasisControllerI, 1E-10);

      auto const loopEvalFunction3 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadID) {
        const unsigned int ij = i * nBFs_I + j;
        const unsigned int ji = j * nBFs_I + i;

        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            // auto& pp = (*densityMatrices)[iSet][iGuess];
            auto& pf = (*f[threadID])[iSet][iGuess];
            // auto coeff = coefficients[iSet][iGuess].data();
            double coul = 0.0;
            coul = integral[0] * coefficients[iSet][iGuess](K);
            *(pf.alpha.data() + ij) += coul;
            *(pf.beta.data() + ij) += coul;
            if (i != j) {
              // coul = integral[0] * coeff[K];
              coul = integral[0] * coefficients[iSet][iGuess](K);
              *(pf.alpha.data() + ji) += coul;
              *(pf.beta.data() + ji) += coul;
            }
          }
        }
      };

      looper3.loop(loopEvalFunction3);

      TwoElecThreeCenterIntLooper looper4(libint2::Operator::coulomb, 0, basisControllerI, auxBasisControllerJ, 1E-10);
      auto const loopEvalFunction4 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadID) {
        const unsigned int ij = i * nBFs_I + j;
        const unsigned int ji = j * nBFs_I + i;

        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            auto& pf = (*f[threadID])[iSet][iGuess];
            double coul = 0.0;
            coul = integral[0] * coefficients[iSet][iGuess](K + nAuxFunctionsI);
            *(pf.alpha.data() + ij) += coul;
            *(pf.beta.data() + ij) += coul;

            if (i != j) {
              coul = integral[0] * coefficients[iSet][iGuess](K + nAuxFunctionsI);
              *(pf.alpha.data() + ji) += coul;
              *(pf.beta.data() + ji) += coul;
            }
          }
        }
      };
      looper4.loop(loopEvalFunction4);

#ifdef _OPENMP
      // const unsigned int nThreads = omp_get_max_threads();
      for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
        for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
          for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
            (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
          }
        }
        delete f[iThread];
      }
#endif
    }

    // If no RI is used:
  }
  else {
    // Thread safety
    std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>*> f;
#ifdef _OPENMP
    const unsigned int nThreads = omp_get_max_threads();
    f.push_back(&(*fock));
    for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
      f.push_back(new std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>);
      (*f[iThread]).resize(this->_nSet);
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          (*f[iThread])[iSet].emplace_back(this->_lrscf[I]->getBasisController());
        }
      }
    }
#else
    f.push_back(&(*fock));
#endif

    // Function to calculate Coulomb pseudo-Fock matrix (for I=J)
    //
    // F_{(ij} = \sum_{kl} (ij|kl) P_{kl}
    //
    // Here, each integral is multiplied with a density matrix element. For LRSCF problems,
    // the (pseudo) density matrix is non symmetric. Every integral (ij|kl) is needed for two
    // Fock-matrix elements, F_{ij} and F_{kl}. Since the resulting Coulomb Fock-like matrix
    // is symmetric, we can add P_{kl} to P_{lk} (and P_{ij} to P_{ji} to save some time.
    // For the diagonal elements of P and F, we then need to divide by a factor of two, which is
    // done by the the factor perm.
    auto distribute_II = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, const double integral,
                             unsigned int threadId) {
      const unsigned int kl = k * nBFs_I + l;
      const unsigned int lk = l * nBFs_I + k;
      const unsigned int ij = i * nBFs_I + j;
      const unsigned int ji = j * nBFs_I + i;
      double perm = 1.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;
      perm *= (i == k) ? (j == l ? 0.5 : 1.0) : 1.0;
      const double coul = perm * integral;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& pp = (*densityMatrices)[iSet][iGuess];
          auto& pf = (*f[threadId])[iSet][iGuess];

          double coul1 = 0.0;
          double coul2 = 0.0;

          coul1 += *(pp.alpha.data() + kl);
          coul1 += *(pp.alpha.data() + lk);
          coul2 += *(pp.alpha.data() + ij);
          coul2 += *(pp.alpha.data() + ji);
          coul1 += *(pp.beta.data() + kl);
          coul1 += *(pp.beta.data() + lk);
          coul2 += *(pp.beta.data() + ij);
          coul2 += *(pp.beta.data() + ji);

          coul1 *= coul;
          coul2 *= coul;

          *(pf.alpha.data() + ij) += coul1;
          *(pf.alpha.data() + ji) += coul1;
          *(pf.alpha.data() + kl) += coul2;
          *(pf.alpha.data() + lk) += coul2;
          *(pf.beta.data() + ij) += coul1;
          *(pf.beta.data() + ji) += coul1;
          *(pf.beta.data() + kl) += coul2;
          *(pf.beta.data() + lk) += coul2;
        }
      }
    };

    // Function to calculate Coulomb pseudo-Fock matrix (for I!=J)
    auto distribute_IJ = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, const double integral,
                             unsigned int threadId) {
      const unsigned int kl = k * nBFs_J + l;
      const unsigned int lk = l * nBFs_J + k;
      const unsigned int ij = i * nBFs_I + j;
      const unsigned int ji = j * nBFs_I + i;
      // Permutations
      double perm = 1.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;
      double coul = perm * integral;
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& pp = (*densityMatrices)[iSet][iGuess];
          auto& pf = (*f[threadId])[iSet][iGuess];
          double coulJ = 0;

          coulJ += *(pp.alpha.data() + kl);
          coulJ += *(pp.alpha.data() + lk);
          coulJ += *(pp.beta.data() + kl);
          coulJ += *(pp.beta.data() + lk);

          coulJ *= coul;

          *(pf.alpha.data() + ij) += coulJ;
          *(pf.alpha.data() + ji) += coulJ;
          *(pf.beta.data() + ij) += coulJ;
          *(pf.beta.data() + ji) += coulJ;
        }
      }
    };

    double maxDens = 0.0;
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto& P = (*densityMatrices)[iSet][iGuess];
        for_spin(P) {
          maxDens = std::max(maxDens, P_spin.array().abs().maxCoeff());
        };
      }
    }
    // Use smalles prescreening threshold of subsystem I and J
    double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                            this->_lrscf[J]->getSysSettings().basis.integralThreshold);

    auto prescreeningFunc = [&](unsigned int i, unsigned int j, unsigned int k, unsigned int l, unsigned int nI,
                                unsigned int nJ, unsigned int nK, unsigned int nL, double schwartz) {
      (void)i;
      (void)j;
      (void)k;
      (void)l;
      (void)nI;
      (void)nJ;
      (void)nK;
      (void)nL;
      if (maxDens * schwartz < prescreeningThreshold)
        return true;
      return false;
    };

    // Calculate pseudo Fock matrices
    if (I == J) {
      TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                        prescreeningThreshold);
      looper.loopNoDerivative(distribute_II, prescreeningFunc);
    }
    else if (I != J) {
      CoulombInteractionIntLooper looper(libint2::Operator::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                         this->_lrscf[J]->getBasisController(), prescreeningThreshold);
      looper.loopNoDerivative(distribute_IJ, prescreeningFunc);
    }
    else {
      assert(false);
    }
#ifdef _OPENMP
    // const unsigned int nThreads = omp_get_max_threads();
    for (unsigned int iThread = 1; iThread < nThreads; ++iThread) {
      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          (*f[0])[iSet][iGuess] += (*f[iThread])[iSet][iGuess];
        }
      }
      delete f[iThread];
    }
#endif
  }
  Timings::timeTaken("LRSCF -   Fock-like matrix: J");
  return fock;
}

template class CoulombSigmaVector<Options::SCF_MODES::RESTRICTED>;
template class CoulombSigmaVector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
