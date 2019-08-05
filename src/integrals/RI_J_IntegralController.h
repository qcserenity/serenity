/**
 * @file RI_J_IntegralController.h
 *
 * @date Mar 8, 2016
 * @author Kevin Klahr
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

#ifndef BASICS_INTEGRALS_RI_J_INTEGRALCONTROLLER_H_
#define BASICS_INTEGRALS_RI_J_INTEGRALCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "integrals/wrappers/Libint.h"
#include "math/Matrix.h"
#include "memory/MemoryManager.h"
#include "notification/ObjectSensitiveClass.h"
#include "integrals/looper/TwoElecThreeCenterIntLooper.h"
#include "integrals/looper/ABTwoElecThreeCenterIntLooper.h"
/* Include Std and External Headers */
#include <memory>
#include <mutex>
#include <vector>


namespace Serenity {

/* Forward declarations */
class Basis;
class BasisController;
class ElectronRepulsionIntegralController;
class MatrixInverter;
class MemoryManager;
class SpinFreePotential;
class TotalDensityMatrix;
/**
 * @class RI_J_IntegralController RI_J_IntegralController.h
 * @brief A controller that handles the integrals for the RI-J approximation.
 * All of the 2 center ERIs as well as all or parts of the three center ERIs
 * are calculated once when needed and stored here for further use.
 *
 * Implementation according to:
 * [1] Weigend, F.; Kattannek, M.; Ahlrichs, R.; J.Chem.Phys (2009), 130, 164106 (eq. 4)
 *
 * See also:
 * [2] Neese, F.; J.Comput.Chem. (2003), 24, 1740 (Scheme 3)
 *
 */

class RI_J_IntegralController : public ObjectSensitiveClass<Basis>{

public:
  /**
   * @brief Constructor.
   * @param basisController of the Basis in which the densityMatrix and Fock matrix are/will be
   *                        expressed
   * @param auxBasisController for the Basis, to which the density is fitted
   * @param densityMatrix see IncrementalPotentialAdder for why we want a reference to the density
   *                      matrix here
   * @param eriController used for prescreening
   */
  RI_J_IntegralController(
      std::shared_ptr<BasisController> basisControllerA,
      std::shared_ptr<BasisController> auxBasisController,
      std::shared_ptr<BasisController> basisControllerB = nullptr);
  /**
   * @brief Default destructor.
   */
  virtual ~RI_J_IntegralController() = default;



  /**
   * @brief Loops over 2e- 3-center integrals.
   *
   * This function expects a lambda function as input.
   * The function will be called for each integral
   * (be aware of symmetry!).
   * The function has to expect four arguments:
   * @code
   *   void std::function<void (const unsigned int& i,
   *   const unsigned int& j,
   *   const unsigned int& K,
   *   const double& intValue,
   *   const unsigned int threadId )>
   * @endcode
   * Here i and j are the indices of the 'normal' basis while
   * K is the index of the auxiliary basis.
   * All indices are unfolded (running over all contracted functions
   * in all shells.)
   * The routine uses the symmetry within i and j thus all function calls
   * should expect i<=j and account for the 'missing' calls.
   *
   * The function can than do anything with the given integrals,
   * e.g. sum up:
   * @code
   *  double sum = 0.0;
   *  auto const sumInts = [&sum]
   *                       (const unsigned int&  i,
   *                        const unsigned int&  j,
   *                        const unsigned int&  K,
   *                        const double& intValue,
   *                        const unsigned int threadId) {
   *  sum += intValues(0);
   *  if (i!=j) sum += intValues[0];
   *  };
   *  rijController.loopOver3CInts(sumInts);
   * @endcode
   *
   * The difference between this function and the looper is that this version
   * allows the RI_J_Controller to cache some of the values.
   * The loop then runs over cached values and afterwards all the other values
   * will be recalculated.
   *
   *
   * @param loopEvalFunction The function to use each integral, see above for extended description.
   */
  template<class Func>
  __attribute__((always_inline)) inline void loopOver3CInts(Func distributionFunction) {
    loopOver3CInts(distributionFunction, [](
      unsigned int,unsigned int,unsigned int,unsigned int,double ){return false;});
  }

  template<class Func, class PrescreenFunc>
  __attribute__((always_inline)) inline void loopOver3CInts(Func& loopEvalFunction, PrescreenFunc& prescreenFunc){
    // Determine the mode the method is run in.
    const bool twoBasisMode = _basisControllerB != nullptr;
    auto const interm = [&loopEvalFunction](const unsigned int&  i,
                                const unsigned int&  j,
                                const unsigned int&  K,
                                Eigen::VectorXd&  integral,
                                const unsigned int threadId) {
      loopEvalFunction(i,j,K,integral.data()[0],threadId);
    };


    // cache integrals if none are cached
    if (!_cache) cache3CInts();
    if (!_cache) {
      // if there is nothing cached now there was no space
      // then run the entire loop
      if(!twoBasisMode) {
        TwoElecThreeCenterIntLooper looper(libint2::Operator::coulomb,0,_basisControllerA,_auxBasisController,1E-10);
        looper.loop(interm,prescreenFunc);
      } else {
        ABTwoElecThreeCenterIntLooper looper(libint2::Operator::coulomb,0,_basisControllerA,_basisControllerB,_auxBasisController,1E-10);
        looper.loop(interm);
      }
    } else {
      // if there was a cache, use it,
      auto data = _cache->data();
      const unsigned int colsize = _cache->rows();
      const unsigned int nBFs_A = _basisControllerA->getNBasisFunctions();
      for (unsigned int i=0;i<nBFs_A;++i){
        const unsigned int jEnd = (twoBasisMode)? _basisControllerB->getNBasisFunctions()-1 : i;
#pragma omp parallel for schedule(static, 1)
        for (unsigned int j=0;j<=jEnd;++j){
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
          const unsigned long long col = (!twoBasisMode)? (unsigned long long)((i*(i+1)/2)+j)*colsize :
              (unsigned long long)((i * (_basisControllerB->getNBasisFunctions())+j)*colsize);
          for (unsigned int K=0;K<colsize;++K){
            loopEvalFunction(i,j,K,data[(unsigned long long)col+K],threadId);
          }
        }
      }
      // then run the rest.
      if (_cache->rows()!=_auxBasisController->getNBasisFunctions()){
        if (!twoBasisMode) {
          TwoElecThreeCenterIntLooper looper(libint2::Operator::coulomb,0,
              _basisControllerA,_auxBasisController,1E-10,
              std::pair<unsigned int,unsigned int>(_cache->rows(),
                  _auxBasisController->getNBasisFunctions()));
          looper.loop(interm,prescreenFunc);
        } else {
          ABTwoElecThreeCenterIntLooper looper(libint2::Operator::coulomb,0,
            _basisControllerA,_basisControllerB,_auxBasisController,1E-10,
            std::pair<unsigned int,unsigned int>(_cache->rows(),
                _auxBasisController->getNBasisFunctions()));
        looper.loop(interm);
        }
      }
    }
  };

  /**
   * @brief Empties the cache.
   */
  void clearCache(){
    _cache.reset(nullptr);
  }

private:

  /**
   * @brief Checks if 3 center ints can be cached.
   *        Note: expects the cache to be a nullptr.
   */
  inline void cache3CInts(){
    const bool twoBasisMode = _basisControllerB != nullptr;
    assert(!_cache);
    // check how many ij sets can be stored
    const unsigned int nBFs_A = _basisControllerA->getNBasisFunctions();
    auto memManager = MemoryManager::getInstance();

    long long memPerBlock = (!twoBasisMode)? nBFs_A*(nBFs_A+1)/2 * sizeof(double) :
        nBFs_A * _basisControllerB->getNBasisFunctions() * sizeof(double);
    long long freeMem = memManager->getAvailableSystemMemory();
    unsigned long long nBlocks = 0.85*freeMem/memPerBlock; // keep 15% of memory here for other things that will be cached
    // return if there was nothing to store
    if (freeMem<0) nBlocks=0;
    if (nBlocks==0) return;
    if (nBlocks > _auxBasisController->getNBasisFunctions()){
      nBlocks = _auxBasisController->getNBasisFunctions();
    }
    // store everything that there was space for
    const unsigned int blockSize = (!twoBasisMode)? nBFs_A*(nBFs_A+1)/2 : nBFs_A * _basisControllerB->getNBasisFunctions();
    _cache.reset(new Eigen::MatrixXd( Eigen::MatrixXd::Zero(nBlocks,blockSize)));
    auto cache = _cache->data();
    const auto basisControllerB = _basisControllerB;
    auto const loopEvalFunction = [&cache,&nBlocks,&twoBasisMode,&basisControllerB]
                               (const unsigned int&  i,
                                const unsigned int&  j,
                                const unsigned int&  K,
                                Eigen::VectorXd&  integral,
                                const unsigned int threadId){
      (void)threadId; //no warning, please
      const unsigned long long ijIndex = (!twoBasisMode)? (unsigned long long)(nBlocks*((i*(i+1)/2)+j)) :
          (unsigned long long)(nBlocks*(i * (basisControllerB->getNBasisFunctions()) + j));
      cache[(unsigned long long)K+ijIndex] = integral.data()[0];
    };
    if (!twoBasisMode) {
      TwoElecThreeCenterIntLooper looper(libint2::Operator::coulomb,0,
          _basisControllerA,_auxBasisController,1E-10,
          std::pair<unsigned int,unsigned int>(0,_cache->rows()));
      looper.loop(loopEvalFunction);
    } else {
      ABTwoElecThreeCenterIntLooper looper(libint2::Operator::coulomb,0,
          _basisControllerA,_basisControllerB,_auxBasisController,1E-10,
          std::pair<unsigned int,unsigned int>(0,_cache->rows()));
      looper.loop(loopEvalFunction);
    }
  }

  /// @brief The 2e3e integral cache.
  std::unique_ptr<Eigen::MatrixXd> _cache;

public:


  const Matrix<double>& getInverseM();

  const std::shared_ptr<BasisController> getBasisController();

  const std::shared_ptr<BasisController> getBasisControllerB();

  const std::shared_ptr<BasisController> getAuxBasisController();

  void calculate2CenterIntegrals();

  void initialize();

  void notify() override final {
    _2cIntsAvailable = false;
    _inverseM.resize(0,0);
    _cache=nullptr;
  }

private:

  bool _2cIntsAvailable;

  const std::shared_ptr<BasisController> _basisControllerA;
  std::shared_ptr<BasisController> _basisControllerB;
  const std::shared_ptr<BasisController> _auxBasisController;


  /*
   * The number of basis functions inside the 'normal' basis, i.e. the basis in which e.g.
   * the density and fock matrix are expressed.
   */
  const unsigned int _nBasisFunctions;
  /*
   * The number of basis functions forming the auxiliary basis
   */
  const unsigned int _nAuxFunctions;
  const unsigned int _nAuxFunctionsRed;

  std::vector<double> normFactors;
  /*
   * The inverted Coulomb metric of the auxiliary basis
   */
  Matrix<double> _inverseM;


  const std::shared_ptr<Libint> _libint = Libint::getSharedPtr();

  std::shared_ptr<MemoryManager> _memManager;

};

} /* namespace Serenity */
#endif /* BASICS_INTEGRALS_RI_J_INTEGRALCONTROLLER_H_ */
