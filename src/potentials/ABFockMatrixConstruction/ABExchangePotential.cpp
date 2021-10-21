/**
 * @file ABExchangePotential.cpp
 *
 * @date Jun 20, 2018
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
#include "potentials/ABFockMatrixConstruction/ABExchangePotential.h"
/* Include Serenity Internal Headers */
#include "basis/ABShellPairCalculator.h"
#include "basis/Basis.h"
#include "integrals/wrappers/Libint.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ABExchangePotential<SCFMode>::ABExchangePotential(std::shared_ptr<SystemController> actSystem,
                                                  std::shared_ptr<BasisController> basisA,
                                                  std::shared_ptr<BasisController> basisB,
                                                  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> dMats,
                                                  double exchangeRatio)
  : ABPotential<SCFMode>(basisA, basisB),
    _actSystem(actSystem),
    _libint(Libint::getSharedPtr()),
    _densityMatrices(dMats),
    _exchangeRatio(exchangeRatio) {
  // Basis
  this->_basisA->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_basisB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  // density matrices
  for (const auto& dMat : _densityMatrices) {
    dMat->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  }
}

template<>
SPMatrix<Options::SCF_MODES::RESTRICTED>& ABExchangePotential<Options::SCF_MODES::RESTRICTED>::getMatrix() {
  if (!_abPotential) {
    const auto scfMode = Options::SCF_MODES::RESTRICTED;
    auto actSystem = _actSystem.lock();
    auto intThreshold = actSystem->getSettings().basis.integralThreshold;
    unsigned int nBasisA = this->_basisA->getNBasisFunctions();
    unsigned int nBasisB = this->_basisB->getNBasisFunctions();
    _abPotential.reset(new SPMatrix<scfMode>(nBasisA, nBasisB));
    SPMatrix<scfMode>& f_AB = *_abPotential;
    for_spin(f_AB) {
      f_AB_spin.setZero();
    };
    double pre = 0.5 * _exchangeRatio;
    double maxD = 1;
    for (const auto& densityC : _densityMatrices)
      maxD = std::max(maxD, densityC->getDensityMatrix().total().array().abs().maxCoeff());
    // intialize libint
    LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb;
    _libint->keepEngines(op, 0, 4);
    _libint->initialize_plain(op, 4, std::numeric_limits<double>::epsilon(), maxD,
                              std::max(_basisA->getMaxNumberOfPrimitives(), _basisB->getMaxNumberOfPrimitives()));

    for (const auto& densityC : _densityMatrices) {
      auto densityMatrixC = densityC->getDensityMatrix();
      const auto& basisContC = densityC->getDensityMatrix().getBasisController();
      const Eigen::MatrixXd densMatrixCTotal = densityC->getDensityMatrix().shellWiseAbsMax();

      // Calculate shellPairData
      const auto shellPairsAC = ABShellPairCalculator::calculateShellPairData_AB(this->_basisA, basisContC);
      const auto shellPairsBC = ABShellPairCalculator::calculateShellPairData_AB(this->_basisB, basisContC);

      const auto& basisA = this->_basisA->getBasis();
      const auto& basisB = this->_basisB->getBasis();
      const auto& basisC = basisContC->getBasis();

      /*
       * Thread safety: one buffer for each thread
       */
#ifdef _OPENMP
      std::vector<Eigen::MatrixXd> ints(omp_get_max_threads());
      std::vector<SPMatrix<scfMode>> eriContr(
          omp_get_max_threads(), SPMatrix<scfMode>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
#else
      std::vector<Eigen::MatrixXd> ints(1);
      std::vector<SPMatrix<scfMode>> eriContr(
          1, SPMatrix<scfMode>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
#endif

#pragma omp parallel for schedule(dynamic)
      // loops over shells
      for (int pIndex = shellPairsAC->size() - 1; pIndex >= 0; --pIndex) {
        auto& p = (*shellPairsAC)[pIndex];

#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        const unsigned int a = p.bf1;
        const unsigned int c1 = p.bf2;
        const auto& basA = *basisA[a];
        const auto& basC1 = *basisC[c1];

        for (const auto& q : *shellPairsBC) {
          /*
           * Simple Prescreening, break out of loops early
           */
          if (p.factor * q.factor < intThreshold)
            break;
          const unsigned int c2 = q.bf2;
          double densMax = densMatrixCTotal(c1, c2);
          if (densMax * p.factor * q.factor < intThreshold)
            continue;
          const auto& basC2 = *basisC[c2];
          const unsigned int b = q.bf1;
          const auto& basB = *basisB[b];
          bool significant = _libint->compute(op, 0, basA, basC1, basB, basC2, ints[threadId]);
          const Eigen::VectorXd integral = ints[threadId].col(0);
          if (significant) {
            // unpack and run
            unsigned int counter = 0;
            for (unsigned int A = 0; A < basisA[a]->getNContracted(); ++A) {
              const unsigned int aa = this->_basisA->extendedIndex(a) + A;
              for (unsigned int C1 = 0; C1 < basisC[c1]->getNContracted(); ++C1) {
                const unsigned int cc1 = basisContC->extendedIndex(c1) + C1;
                for (unsigned int B = 0; B < basisB[b]->getNContracted(); ++B) {
                  const unsigned int bb = this->_basisB->extendedIndex(b) + B;
                  for (unsigned int C2 = 0; C2 < basisC[c2]->getNContracted(); ++C2, ++counter) {
                    const unsigned int cc2 = basisContC->extendedIndex(c2) + C2;
                    //                    if (cc2>cc1) continue;
                    const double set(integral(counter));
                    auto& eris = eriContr[threadId];
                    eris(aa, bb) -= pre * densityMatrixC(cc1, cc2) * set;

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
        const auto& eris = eriContr[t];
        f_AB += eris;
      }
#else
      const auto& eris = eriContr[0];
      f_AB += eris;
#endif
    } /* for densityC */
    // finalize libint
    _libint->freeEngines(op, 0, 4);
  } /* if !_abPotential */
  return *_abPotential;
}
template<>
SPMatrix<Options::SCF_MODES::UNRESTRICTED>& ABExchangePotential<Options::SCF_MODES::UNRESTRICTED>::getMatrix() {
  if (!_abPotential) {
    const auto scfMode = Options::SCF_MODES::UNRESTRICTED;
    auto actSystem = _actSystem.lock();
    auto intThreshold = actSystem->getSettings().basis.integralThreshold;
    unsigned int nBasisA = this->_basisA->getNBasisFunctions();
    unsigned int nBasisB = this->_basisB->getNBasisFunctions();
    _abPotential.reset(new SPMatrix<scfMode>(nBasisA, nBasisB));
    SPMatrix<scfMode>& f_AB = *_abPotential;
    for_spin(f_AB) {
      f_AB_spin.setZero();
    };
    double pre = 1.0 * _exchangeRatio;
    double maxD = 1;
    for (const auto& densityC : _densityMatrices)
      maxD = std::max(maxD, densityC->getDensityMatrix().total().array().abs().maxCoeff());
    // intialize libint
    LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb;
    _libint->keepEngines(op, 0, 4);
    _libint->initialize_plain(op, 4, std::numeric_limits<double>::epsilon(), maxD,
                              std::max(_basisA->getMaxNumberOfPrimitives(), _basisB->getMaxNumberOfPrimitives()));

    for (const auto& densityC : _densityMatrices) {
      auto densityMatrixC = densityC->getDensityMatrix();
      const auto& basisContC = densityC->getDensityMatrix().getBasisController();
      const Eigen::MatrixXd densMatrixCTotal = densityC->getDensityMatrix().shellWiseAbsMax().total();

      // Calculate shellPairData
      const auto shellPairsAC = ABShellPairCalculator::calculateShellPairData_AB(this->_basisA, basisContC);
      const auto shellPairsBC = ABShellPairCalculator::calculateShellPairData_AB(this->_basisB, basisContC);

      const auto& basisA = this->_basisA->getBasis();
      const auto& basisB = this->_basisB->getBasis();
      const auto& basisC = basisContC->getBasis();

      /*
       * Thread safety: one buffer for each thread
       */
#ifdef _OPENMP
      std::vector<Eigen::MatrixXd> ints(omp_get_max_threads());
      std::vector<SPMatrix<scfMode>> eriContr(
          omp_get_max_threads(), SPMatrix<scfMode>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
#else
      std::vector<Eigen::MatrixXd> ints(1);
      std::vector<SPMatrix<scfMode>> eriContr(
          1, SPMatrix<scfMode>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
#endif

#pragma omp parallel for schedule(static, 1)
      // loops over shells
      for (int pIndex = shellPairsAC->size() - 1; pIndex >= 0; --pIndex) {
        auto& p = (*shellPairsAC)[pIndex];

#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        const unsigned int a = p.bf1;
        const unsigned int c1 = p.bf2;
        const auto& basA = *basisA[a];
        const auto& basC1 = *basisC[c1];

        for (const auto& q : *shellPairsBC) {
          /*
           * Simple Prescreening, break out of loops early
           */
          if (p.factor * q.factor < intThreshold)
            break;
          const unsigned int c2 = q.bf2;
          double densMax = densMatrixCTotal(c1, c2);
          if (densMax * p.factor * q.factor < intThreshold)
            continue;
          const auto& basC2 = *basisC[c2];
          const unsigned int b = q.bf1;
          const auto& basB = *basisB[b];
          bool significant = _libint->compute(op, 0, basA, basC1, basB, basC2, ints[threadId]);
          const Eigen::VectorXd integral = ints[threadId].col(0);
          if (significant) {
            // unpack and run
            unsigned int counter = 0;
            for (unsigned int A = 0; A < basisA[a]->getNContracted(); ++A) {
              const unsigned int aa = this->_basisA->extendedIndex(a) + A;
              for (unsigned int C1 = 0; C1 < basisC[c1]->getNContracted(); ++C1) {
                const unsigned int cc1 = basisContC->extendedIndex(c1) + C1;
                for (unsigned int B = 0; B < basisB[b]->getNContracted(); ++B) {
                  const unsigned int bb = this->_basisB->extendedIndex(b) + B;
                  for (unsigned int C2 = 0; C2 < basisC[c2]->getNContracted(); ++C2, ++counter) {
                    const unsigned int cc2 = basisContC->extendedIndex(c2) + C2;
                    //                    if (cc2>cc1) continue;

                    const double set(integral(counter));
                    auto& eris = eriContr[threadId];
                    eris.alpha(aa, bb) -= pre * densityMatrixC.alpha(cc1, cc2) * set;
                    eris.beta(aa, bb) -= pre * densityMatrixC.beta(cc1, cc2) * set;

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
        const auto& eris = eriContr[t];
        f_AB.alpha += eris.alpha;
        f_AB.beta += eris.beta;
      }
#else
      const auto& eris = eriContr[0];
      f_AB.alpha += eris.alpha;
      f_AB.beta += eris.beta;
#endif
    } /* for densityC */
    // finalize libint
    _libint->freeEngines(op, 0, 4);
  } /* if !_abPotential */
  return *_abPotential;
}
template class ABExchangePotential<Options::SCF_MODES::RESTRICTED>;
template class ABExchangePotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
