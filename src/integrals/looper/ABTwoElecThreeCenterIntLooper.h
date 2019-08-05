/**
 * @file ABTwoElecThreeCenterIntLooper.h
 *
 * @date Jun 20, 2018
 * @author Moritz Bensberg
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef INTEGRALS_LOOPER_ABTWOELECTHREECENTERINTLOOPER_H_
#define INTEGRALS_LOOPER_ABTWOELECTHREECENTERINTLOOPER_H_

/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "integrals/wrappers/Libint.h"
#include "basis/ABShellPairCalculator.h"
/* Include Std and External Headers */
#include <memory>


namespace Serenity {
/**
 * @class ABTwoElecThreeCenterIntLooper ABTwoElecThreeCenterIntLooper.h
 * @brief Looper class for the evaluation of AB-three-center-integrals.
 */
class ABTwoElecThreeCenterIntLooper  {
public:
  /**
   * @brief Constructor.
   * @param op The libint-operator.
   * @param deriv The order of the derivative.
   * @param basisA The basis controller A.
   * @param basisB The basis controller B.
   * @param auxbasis The auxillary basis controller.
   * @param prescreeningThreshold The prescreening threshold.
   * @param kRange The evaluation range for the auxiallary functions.
   */
  ABTwoElecThreeCenterIntLooper(
      libint2::Operator op,
      const unsigned int deriv,
      std::shared_ptr<BasisController> basisA,
      std::shared_ptr<BasisController> basisB,
      std::shared_ptr<BasisController> auxbasis,
      double prescreeningThreshold,
      std::pair<unsigned int, unsigned int> kRange = {0,0}) :
      _op(op),
      _deriv(deriv),
      _basisControllerA(basisA),
      _basisControllerB(basisB),
      _auxbasis(auxbasis),
      _prescreeningThreshold(prescreeningThreshold),
      _kRange(kRange){
    assert(_basisControllerA);
    assert(_basisControllerB);
    assert(_auxbasis);
  }
  /**
   * @brief Destructor.
   */
  virtual ~ABTwoElecThreeCenterIntLooper()=default;
  /**
   * @brief Performs the loop.
   * @param loopEvalFunction The evaluation function for the looper.
   *        For examples and more information see TwoElecThreeCenterIntLooper.h
   */
  template<class Func>
  __attribute__((always_inline)) inline void loop(Func loopEvalFunction){
    const auto& shellPairsAB = ABShellPairCalculator::calculateShellPairData_AB(
        _basisControllerA,
        _basisControllerB);
    // intialize libint
    auto& libint = Libint::getInstance();
    if(_op==libint2::Operator::delta){
      libint.initialize(_op,_deriv,4);
    } else{
      libint.initialize(_op,_deriv,3);
    }
    auto& basisA = _basisControllerA->getBasis();
    auto& basisB = _basisControllerB->getBasis();
    auto& auxbasis = _auxbasis->getBasis();
    const auto& riPrescreeningFactors = _auxbasis->getRIPrescreeningFactors();
    if (_kRange.first ==  _auxbasis->getNBasisFunctions()) return;
    if (_kRange.second == 0) _kRange.second =  _auxbasis->getNBasisFunctions();



#ifdef _OPENMP
   std::vector<Eigen::MatrixXd> integrals(omp_get_max_threads());
#else
   std::vector<Eigen::MatrixXd> integrals(1);
#endif

#pragma omp parallel for schedule(static, 1)
     for (int kIndex = _auxbasis->reducedIndex(_kRange.second-1); kIndex >=(int)_auxbasis->reducedIndex(_kRange.first); --kIndex){
      auto& k = (*riPrescreeningFactors)[kIndex];

#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      const unsigned int shellK = k.bf1;
      const auto& auxbasK= *auxbasis[shellK];
      const unsigned int nK = auxbasis[shellK]->getNContracted();
      // loops over shells
      for (auto& p : *shellPairsAB) {

        if(p.factor*k.factor < _prescreeningThreshold) break;

//        const bool swap = (basisA[p.bf1]->getAngularMomentum()<basisB[p.bf2]->getAngularMomentum());

        const unsigned int iA = p.bf1;
        const unsigned int jB = p.bf2;

        const auto& basIA= *basisA[iA];
        const auto& basJB= *basisB[jB];

        const unsigned int nIA = basisA[iA]->getNContracted();
        const unsigned int nJB = basisB[jB]->getNContracted();
          // calculate integrals
          bool take;
          if(_op==libint2::Operator::delta){
            take=libint.compute(_op,_deriv,
                auxbasK,
                basIA,
                basJB,
                libint2::Shell::unit(),
                integrals[threadId]);
          }else{
            take=libint.compute(_op,_deriv,
                auxbasK,
                basIA,
                basJB,
                integrals[threadId]);
          }
          if (take){
            // unpack and run
            for (unsigned int K=0;K<nK;++K){
              const unsigned int kk = _auxbasis->extendedIndex(shellK)+K;
              if(kk >= _kRange.second) continue;
              if(kk < _kRange.first) continue;
              for (unsigned int A=0;A<nIA;++A){
                const unsigned int aa = _basisControllerA->extendedIndex(iA)+A;
                for (unsigned int B=0;B<nJB;++B){
                  const unsigned int bb = _basisControllerB->extendedIndex(jB)+B;
//                  if (ii>jj) continue;
                  Eigen::VectorXd set(integrals[threadId].row(K*nJB*nIA+A*nJB+B));
                  loopEvalFunction(aa,bb,kk,set,threadId);
                } /* primitives of bb -> B */
              } /* primitives of aa -> A */
            } /* primitives of k -> K */
          } /* if (prescreen) */
        } /* p/shellpair */
      } /* K/shellK */
    // finalize libint
    libint.finalize(_op,_deriv,3);
  }


private:
  /// @brief The kernel/operator as libint enum.
  libint2::Operator _op;
  /// @brief The derivative level.
  const unsigned int _deriv;
  /// @brief The basis.
  std::shared_ptr<BasisController> _basisControllerA;
  /// @brief The basis.
  std::shared_ptr<BasisController> _basisControllerB;
  /// @brief The auxiliary basis.
  std::shared_ptr<BasisController> _auxbasis;
  /// @brief A threshold for the Schwartz conditions.
  double _prescreeningThreshold;
  /// @brief The range for the K shells.
  std::pair<unsigned int, unsigned int> _kRange;
};


} /* namespace Serenity */
#endif /* INTEGRALS_LOOPER_ABTWOELECTHREECENTERINTLOOPER_H_ */
