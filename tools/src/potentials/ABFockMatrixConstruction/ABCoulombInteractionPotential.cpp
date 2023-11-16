/**
 * @file ABCoulombInteractionPotential.cpp
 *
 * @date May 8, 2018
 * @author Moritz Bensberg
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
#include "potentials/ABFockMatrixConstruction/ABCoulombInteractionPotential.h"
/* Include Serenity Internal Headers */
#include "basis/ABShellPairCalculator.h"
#include "basis/BasisFunctionMapper.h"
#include "basis/CustomBasisController.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "integrals/looper/ABTwoElecThreeCenterIntLooper.h"
#include "integrals/wrappers/Libint.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ABCoulombInteractionPotential<SCFMode>::ABCoulombInteractionPotential(
    std::shared_ptr<SystemController> actSystem, std::shared_ptr<BasisController> basisA,
    std::shared_ptr<BasisController> basisB,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensityMatrixController, bool topDown,
    Options::DENS_FITS densFitJ, std::shared_ptr<BasisController> auxBasisAB,
    std::vector<std::shared_ptr<BasisController>> envAuxBasisController)
  : ABPotential<SCFMode>(basisA, basisB),
    _actSystem(actSystem),
    _envDMatController(envDensityMatrixController),
    _libint(Libint::getSharedPtr()),
    _topDown(topDown),
    _mode(densFitJ),
    _auxBasisControllerAB(auxBasisAB),
    _envAuxBasisController(envAuxBasisController) {
  // Setting the notifying system up and sanity checks.
  // A,B basis
  this->_basisA->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_basisB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  // Densities
  for (const auto& mat : envDensityMatrixController) {
    assert(mat);
    mat->getDensityMatrix().getBasisController()->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
    mat->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  }
  // RI specials
  if (_mode == Options::DENS_FITS::RI) {
    // AB aux.
    assert(_auxBasisControllerAB);
    // env aux.
    for (const auto& envAuxBas : _envAuxBasisController) {
      assert(envAuxBas);
      envAuxBas->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
    }
    assert(_envAuxBasisController.size() == _envDMatController.size());
    _auxBasisControllerAB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
    for (const auto& envAux : _envAuxBasisController) {
      envAux->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
    }
    // Special treatment for top-down calculations. This will effectively reduce the basis set in _superABBasis.
    _superABBasis = this->_basisB;
    if (!_topDown) {
      Basis superABBasis;
      superABBasis.reserve(this->_basisA->getBasis().size() + this->_basisB->getBasis().size());
      superABBasis.insert(superABBasis.end(), this->_basisA->getBasis().begin(), this->_basisA->getBasis().end());
      superABBasis.insert(superABBasis.end(), this->_basisB->getBasis().begin(), this->_basisB->getBasis().end());
      _superABBasis = std::shared_ptr<BasisController>(new CustomBasisController(superABBasis, "tmpSuperABBasis"));
    }
  }
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode>& ABCoulombInteractionPotential<SCFMode>::getMatrix() {
  if (!_abPotential) {
    _abPotential.reset(new SPMatrix<SCFMode>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
    SPMatrix<SCFMode>& f_AB = *_abPotential;
    for_spin(f_AB) {
      f_AB_spin.setZero();
    };
    if (_mode == Options::DENS_FITS::NONE) {
      calculateFockMatrixNoFitting();
    }
    else if (_mode == Options::DENS_FITS::RI) {
      calculateFockMatrixRI();
    }
    else {
      assert(false && "The required type of density fitting is not implemented!");
    }
  } /* if !abPotential */
  return *_abPotential;
}

template<Options::SCF_MODES SCFMode>
void ABCoulombInteractionPotential<SCFMode>::calculateFockMatrixNoFitting() {
  // intialize libint
  double maxD = 1;
  for (const auto& densityC : _envDMatController)
    maxD = std::max(maxD, densityC->getDensityMatrix().total().array().abs().maxCoeff());
  // intialize libint
  LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb;
  _libint->keepEngines(op, 0, 4);
  _libint->initialize_plain(op, 4, std::numeric_limits<double>::epsilon(), maxD,
                            std::max(this->_basisA->getMaxNumberOfPrimitives(), this->_basisB->getMaxNumberOfPrimitives()));
  for (const auto& densityC : _envDMatController) {
    auto totDensC = densityC->getDensityMatrix().total();
    auto actSystem = _actSystem.lock();
    auto intThreshold = actSystem->getSettings().basis.integralThreshold;
    SPMatrix<SCFMode>& f_AB = *_abPotential;
    auto const basisContC = densityC->getDensityMatrix().getBasisController();

    if (!_shellPairsAB)
      _shellPairsAB = ABShellPairCalculator::calculateShellPairData_AB(this->_basisA, this->_basisB);
    // Calculate AB ShellPairData
    const auto& shellPairsC = basisContC->getShellPairData();
    const auto& basisA = this->_basisA->getBasis();
    const auto& basisB = this->_basisB->getBasis();
    const auto& basisC = basisContC->getBasis();

    /*
     * Thread safety: one buffer for each thread
     */
#ifdef _OPENMP
    std::vector<Eigen::MatrixXd> ints(omp_get_max_threads());
    std::vector<SPMatrix<Options::SCF_MODES::RESTRICTED>> eriContr(
        omp_get_max_threads(),
        SPMatrix<Options::SCF_MODES::RESTRICTED>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
#else
    std::vector<Eigen::MatrixXd> ints(1);
    std::vector<SPMatrix<Options::SCF_MODES::RESTRICTED>> eriContr(
        1, SPMatrix<Options::SCF_MODES::RESTRICTED>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
    eriContr[0].setZero();
#endif
    // loops over shells
#pragma omp parallel for schedule(static, 1)
    for (int pIndex = _shellPairsAB->size() - 1; pIndex >= 0; --pIndex) {
      const auto& p = (*_shellPairsAB)[pIndex];

#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      const unsigned int a = p.bf1;
      const unsigned int b = p.bf2;
      const auto& basA = *basisA[a];
      const auto& basB = *basisB[b];

      for (const auto& q : *shellPairsC) {
        /*
         * Simple Prescreening, break out of loops early
         */
        if (p.factor * q.factor < intThreshold)
          break;
        const unsigned int c1 = q.bf1;
        const unsigned int c2 = q.bf2;
        const auto& basC1 = *basisC[c1];
        const auto& basC2 = *basisC[c2];
        const unsigned int firstC1 = basisContC->extendedIndex(c1);
        const unsigned int firstC2 = basisContC->extendedIndex(c2);

        if (firstC2 > firstC1)
          continue;

        // calculate integrals
        bool significant = _libint->compute(op, 0, basA, basB, basC1, basC2, ints[threadId]);
        const Eigen::VectorXd integral = ints[threadId].col(0);
        if (significant) {
          // unpack and run
          unsigned int counter = 0;
          for (unsigned int A = 0; A < basisA[a]->getNContracted(); ++A) {
            const unsigned int aa = this->_basisA->extendedIndex(a) + A;
            for (unsigned int B = 0; B < basisB[b]->getNContracted(); ++B) {
              const unsigned int bb = this->_basisB->extendedIndex(b) + B;
              for (unsigned int C1 = 0; C1 < basisC[c1]->getNContracted(); ++C1) {
                const unsigned int cc1 = basisContC->extendedIndex(c1) + C1;
                for (unsigned int C2 = 0; C2 < basisC[c2]->getNContracted(); ++C2, ++counter) {
                  const unsigned int cc2 = basisContC->extendedIndex(c2) + C2;
                  if (cc2 > cc1)
                    continue;

                  double pre = (cc1 != cc2) ? 2.0 : 1.0;
                  const double set(integral(counter));
                  eriContr[threadId](aa, bb) += pre * totDensC(cc1, cc2) * set;

                } /* primitives of c2 -> C2  */
              }   /* primitives of c1 -> C1  */
            }     /* primitives of a -> A  */
          }       /* primitives of b -> B  */
        }         /* if (prescreen) */
      }           /* q/shellPairsC */
    }             /* p/shellPairsAB */

#ifdef _OPENMP
    // sum over all threads
    for (int t = 0; t < omp_get_max_threads(); ++t) {
      for_spin(f_AB) {
        f_AB_spin += eriContr[t];
      };
    }
#else
    for_spin(f_AB) {
      f_AB_spin += eriContr[0];
    };
#endif
  } /* for densityC */
  // finalize libint
  _libint->freeEngines(op, 0, 4);
}

template<Options::SCF_MODES SCFMode>
void ABCoulombInteractionPotential<SCFMode>::calculateFockMatrixRI() {
  /*
   * Calculate the Coulomb interaction using the RI-approximation
   * What happens here:
   *   1. Loop over all interacting densities C.
   *   2. If C != A
   *      --> Three different basis sets and two aux. basis sets.
   *      --> We need to transform the density of C using the AB-aux. basis and
   *          the C-aux. basis
   *   3. If C == A
   *      --> This is effectively a outer diagonal block of the normal Coulomb contribution
   *          of a one system with itself.
   *      --> Two different basis sets and one aux. basis set.
   *      --> Uses caching for the integrals and is much cheaper than the case above.
   */
  SPMatrix<SCFMode>& f_AB = *_abPotential;

  auto nABAuxBasFunc = _auxBasisControllerAB->getNBasisFunctions();
  BasisFunctionMapper abBasisMapper(_superABBasis);
  BasisFunctionMapper abAuxBasisMapper(_auxBasisControllerAB);

  /*
   * 1. Loop over all interacting densities C.
   */
  for (unsigned int iC = 0; iC < _envDMatController.size(); ++iC) {
    const auto& densityC = _envDMatController[iC]->getDensityMatrix().total();
    const auto& basisControllerC = densityC.getBasisController();
    const auto& auxBasisControllerC = _envAuxBasisController[iC];

    // Simple prescreenfunction
    // Maximum absolute value in densityMatrix
    auto actSystem = _actSystem.lock();
    const double threshold = actSystem->getSettings().basis.integralThreshold;
    const auto maxDens = densityC.shellWiseAbsMax();
    auto prescreeningFunc = [&](const unsigned i, const unsigned j, const unsigned, const double schwartz) {
      if (std::abs(maxDens(i, j) * schwartz) < threshold)
        return true;
      return false;
    };
#ifdef _OPENMP
    // create a vector of matrices for each thread
    Eigen::MatrixXd sumMat(nABAuxBasFunc, omp_get_max_threads());
    std::vector<SPMatrix<RESTRICTED>> eriContr(
        omp_get_max_threads(), SPMatrix<RESTRICTED>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
#else
    // or just one
    Eigen::MatrixXd sumMat(nABAuxBasFunc, 1);
    std::vector<SPMatrix<RESTRICTED>> eriContr(
        1, SPMatrix<RESTRICTED>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
    eriContr[0].setZero();
#endif
    sumMat.setZero();

    if (this->_basisA != basisControllerC) {
      // Build A+B+C aux-basis in order to calculate the inverse two center matrix.
      auto superABCAuxBasisController = abAuxBasisMapper.getCombinedBasis(auxBasisControllerC);
      auto differentialAuxBasisController = abAuxBasisMapper.getDifferentialBasis(auxBasisControllerC);
      const unsigned int nDiffAuxBasisFunc =
          (differentialAuxBasisController) ? differentialAuxBasisController->getNBasisFunctions() : 0;
#ifdef _OPENMP
      Eigen::MatrixXd sumMatC(nDiffAuxBasisFunc, omp_get_max_threads());
#else
      Eigen::MatrixXd sumMatC(nDiffAuxBasisFunc, 1);
#endif
      sumMatC.setZero();
      TwoElecThreeCenterIntLooper looper1(LIBINT_OPERATOR::coulomb, 0, basisControllerC, _auxBasisControllerAB,
                                          basisControllerC->getPrescreeningThreshold());
      auto const loopEvalFunction1 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadId) {
        sumMat(K, threadId) += (i == j ? 1.0 : 2.0) * integral[0] * densityC(i, j);
      };
      looper1.loop(loopEvalFunction1, prescreeningFunc, maxDens.maxCoeff());
      Eigen::VectorXd sumPMuNu_DMuNuAB = sumMat.rowwise().sum();
      sumMat.resize(1, 1);
      Eigen::VectorXd sumPMuNu_DMuNu(sumPMuNu_DMuNuAB.rows() + nDiffAuxBasisFunc);
      if (differentialAuxBasisController) {
        TwoElecThreeCenterIntLooper looper2(LIBINT_OPERATOR::coulomb, 0, basisControllerC,
                                            differentialAuxBasisController, basisControllerC->getPrescreeningThreshold());
        auto const loopEvalFunction2 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                           Eigen::VectorXd& integral, const unsigned int threadId) {
          sumMatC(K, threadId) += (i == j ? 1.0 : 2.0) * integral[0] * densityC(i, j);
        };
        looper2.loop(loopEvalFunction2, prescreeningFunc, maxDens.maxCoeff());
        Eigen::VectorXd sumPMuNu_DMuNuC = sumMatC.rowwise().sum();
        sumMatC.resize(1, 1);
        sumPMuNu_DMuNu << sumPMuNu_DMuNuAB, sumPMuNu_DMuNuC;
      }
      else {
        sumPMuNu_DMuNu << sumPMuNu_DMuNuAB;
      } // else differentialAuxBasisController
      auto ri_j_IntController =
          RI_J_IntegralControllerFactory::getInstance().produce(_superABBasis, superABCAuxBasisController);
      Eigen::VectorXd coefficients = ri_j_IntController->getLLTMetric().solve(sumPMuNu_DMuNu).eval();
      auto coeffptr = coefficients.data();
      double prescreeningThresholdA = this->_basisA->getPrescreeningThreshold();
      double prescreeningThresholdB = this->_basisB->getPrescreeningThreshold();
      double prescreeningThreshold = std::min(prescreeningThresholdA, prescreeningThresholdB);
      ABTwoElecThreeCenterIntLooper looper3(LIBINT_OPERATOR::coulomb, 0, this->_basisA, this->_basisB,
                                            _auxBasisControllerAB, prescreeningThreshold);
      auto const loopEvalFunction3 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         Eigen::VectorXd& integral, const unsigned int threadID) {
        eriContr[threadID](i, j) += integral[0] * coeffptr[K];
      };
      looper3.loop(loopEvalFunction3, coefficients.array().abs().maxCoeff());
      if (differentialAuxBasisController) {
        ABTwoElecThreeCenterIntLooper looper4(LIBINT_OPERATOR::coulomb, 0, this->_basisA, this->_basisB,
                                              differentialAuxBasisController, prescreeningThreshold);
        auto const loopEvalFunction4 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                           Eigen::VectorXd& integral, const unsigned int threadID) {
          eriContr[threadID](i, j) += integral[0] * coeffptr[K + nABAuxBasFunc];
        };
        looper4.loop(loopEvalFunction4, coefficients.array().abs().maxCoeff());
      } // if differentialAuxBasisController
    }   /* A != C */
    else {
      if (!_ri_j_IntegralController_A_AB) {
        _ri_j_IntegralController_A_AB =
            RI_J_IntegralControllerFactory::getInstance().produce(this->_basisA, _auxBasisControllerAB);
      }
      auto const loopEvalFunction1 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         const double& integral, const unsigned int threadId) {
        sumMat(K, threadId) += (i == j ? 1.0 : 2.0) * integral * densityC(i, j);
      };
      _ri_j_IntegralController_A_AB->loopOver3CInts(loopEvalFunction1, prescreeningFunc);
      Eigen::VectorXd sumPMuNu_DMuNuAB = sumMat.rowwise().sum();
      sumMat.resize(1, 1);
      const auto inverseM = _ri_j_IntegralController_A_AB->getInverseM();
      Eigen::VectorXd coefficients = inverseM * sumPMuNu_DMuNuAB.eval();
      if (!_ri_j_IntegralController_A_B_AB) {
        _ri_j_IntegralController_A_B_AB =
            RI_J_IntegralControllerFactory::getInstance().produce(this->_basisA, _auxBasisControllerAB, this->_basisB);
      }
      auto const loopEvalFunction3 = [&](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                         const double& integral, const unsigned int threadID) {
        eriContr[threadID](i, j) += integral * coefficients(K);
      };
      auto coeffMax = _auxBasisControllerAB->shellWiseAbsMax(coefficients);
      auto prescreeningCoeff = [&](const unsigned, const unsigned, const unsigned K, const double schwartz) {
        if (std::abs(coeffMax(K) * schwartz) < threshold)
          return true;
        return false;
      };
      _ri_j_IntegralController_A_B_AB->loopOver3CInts(loopEvalFunction3, prescreeningCoeff);
    } /* A == C */
#ifdef _OPENMP
    // sum over all threads
    for (int t = 0; t < omp_get_max_threads(); ++t) {
      for_spin(f_AB) {
        f_AB_spin += eriContr[t];
      };
    }
#else
    for_spin(f_AB) {
      f_AB_spin += eriContr[0];
    };
#endif
  } /* for densityC */
}

template class ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED>;
template class ABCoulombInteractionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
