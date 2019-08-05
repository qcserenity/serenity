/**
 * @file   Libint.h
 *
 * @date   28. Juli 2013, 14:12
 * @author Thomas Dresselhaus, Jan Unsleber
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
#ifndef LIBINT_H
#define	LIBINT_H

/* Include Serenity Internal Headers */
#include "math/Matrix.h"
#include "geometry/Atom.h"
#include "basis/BasisController.h"

/* Include Std and External Headers */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <../../../ext/libint/include/libint2.hpp>
#pragma GCC diagnostic pop
#include <mutex>

namespace Serenity {
using namespace Serenity;
/**
 * @class IntegralType
 * @brief A small class to wrap all the possible integral types.
 */
struct IntegralType
{
  /**
   * @brief Default constructor.
   */
  IntegralType() : op() , deriv() ,nCenter(){};
  /**
   * @brief Constructor.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   */
  IntegralType(libint2::Operator op,
      const unsigned int deriv,
      const unsigned int nCenter) : op(op), deriv(deriv),nCenter(nCenter-2) {}
  /// @brief The kernel/operator as libint enum.
  libint2::Operator op;
  /// @brief The derivative level.
  unsigned int deriv;
  /// @brief The number of centers (2-4).
 unsigned int nCenter;
};

/**
 * @class Libint Libint.h
 * @brief Wrapper for the libint2 library.
 *
 * This class follows the singleton design pattern, it is constructed only once
 * for the entire program run, no copies should ever be made
 *
 * For more information on singletons:
 * http://stackoverflow.com/questions/1008019/c-singleton-design-pattern
 *
 * (Yes, this singleton is not correctly destroyed however
 *  we will live with this for now)
 *
 * Currently interfaced: libint-2.2.0-beta1
 * (https://github.com/evaleev/libint)
 *
 * Usage:
 *  The three step usage of libint always looks like this:
 *  1. Initialize libint
 *     In order to compute integrals libint has to be initialized.
 *     (We need to tell it that we will need integrals soon, and which
 *      type of integrals)
 *     The function initialize(...) does just that. Go figure!
 *     The function needs three arguments, first the operator inside the
 *     integral, e.g. the overlap operator.
 *     Secondly the number of geometric derivatives (0=integrals,1=first derivative, ...).
 *     And last but not least the number of basis functions the integral has
 *     thus two for (i|j) or four for (ij|kl) or three for (ij|K)
 *  2. Computation of integrals
 *     The compute() function computes the integrals, using the shells
 *     and again the information about the operator and derivative.
 *     The integrals will be returned as matrix where the rows are one set.
 *     A set contains one value per basis function combination.
 *     The first set is always integral values (thus int(i,0) will give i'th-integral)
 *     The number of sets depends on the number of derivatives and is ordered as follows:
 *     0->[ints],1->[dx,dz,dz],2->[dxdx,dxdy,dxdz,dydy,dydz,dzdz],... .
 *     Here e.g. derivative level 2 also contains all the previous values.
 *     Inside one set the integrals are ordered by their shells, with the first shell being
 *     the outer index (loop).
 *     A mocked example:
 *     @code
 *      compute(libint2::Operator,0,shellI,shellJ,shellK,shellL,int);
 *      for (shellI){
 *       for (shellJ){
 *        for (shellK){
 *         for (shellL ++counter){
 *          (void)integrals(counter,0);
 *         }
 *        }
 *       }
 *      }
 *     @endcode
 *  3. Finalization
 *     Each initialization has to be met with the fitting finalization in order
 *     to tell libint that we are done with those integrals for now.
 *
 *  Note: The construction of the engines actual has a decent overhead,
 *        many subsequent calls to initialize() and finalize() will thus be slow.
 *        here the keepEngines() and freeEngines() functions can be used to keep
 *        engines alive even though finalize is called.
 *        A prime example for this is the SCF where the Coulomb engines are
 *        kept alive until convergence in order not to initialize them in every
 *        iteration.
 */
class Libint {

public:
  /**
   * @brief One of two singleton 'Constructors'
   *
   * @return A reference to the only existing version of the Libint object.
   */
  static Libint& getInstance(){
    return *Libint::getSharedPtr();
  }
  /**
   * @brief One of two Singleton 'Constructors'
   * @return A shared pointer to the only existing version of the Libint object.
   */
  static std::shared_ptr<Libint> getSharedPtr() {
    static std::shared_ptr<Libint> instance(new Libint);
    return instance;
  }
  /**
   * @brief Deleted -> singleton.
   */
  Libint(Libint const&) = delete;
  /**
   * @brief Deleted -> singleton.
   */
  void operator=(Libint const&) = delete;

  /**
   * @brief Destructor
   */
  virtual ~Libint();

  /**
   * @brief Initializes the engines for a specific type of integral.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   * @param pointCharges The charges to be added to the engines
   *                     (needed for Operator::nuclear integrals).
   */
  void initialize(libint2::Operator op,
      const unsigned int deriv,
      const unsigned int nCenter,
      const std::vector<std::pair<double,std::array<double,3>>> pointCharges,
      double mu);

  /**
   * TODO: What is mu?
   *
   * @brief Initializes the engines for a specific type of integral.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   * @param atoms The atoms to be added to the engines
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   */
  void initialize(libint2::Operator op,
      const unsigned int deriv,
      const unsigned int nCenter,
      const std::vector<std::shared_ptr<Atom> >& atoms = std::vector<std::shared_ptr<Atom> >(0),
      double mu = 0.0);

  /**
   * @brief Cleans the engines for a specific type of integral.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   */
  void finalize(libint2::Operator op,
      unsigned int deriv,
      unsigned int nCenter);

  /**
   * @brief Keeps a set of engines form being destroyed by finalize calls.
   *
   * Every Libint::keepEngines has to be matched by a Libint::freeEngines .
   *
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   */
  void keepEngines(libint2::Operator op,
      const unsigned int deriv,
      const unsigned int nCenter);

  /**
   * @brief Allows the destruction and destroys a set of engines if they are present.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   */
  void freeEngines(libint2::Operator op,
      const unsigned int deriv,
      const unsigned int nCenter);

  /**
   * @brief Destroys all engines, call with caution!
   *
   * This is supposed to be called before error exist of the program in order
   * to shorten core dumps and allow for a minimum of cleanliness upon crash.
   */
  void clearAllEngines();

  /**
   * @brief Computes a set (operator/derivative) of integrals for the given shells.
   * @param op     The kernel/operator as libint enum.
   * @param deriv  The derivative level.
   * @param a      The first shell.
   * @param b      The second shell.
   * @param ints   The integrals stored per set (see class description for more information).
   * @return Boolean. True if integrals were calculated, false if they were screened out.
   */
  bool compute(libint2::Operator op,
        unsigned int deriv,
        const libint2::Shell& a,
        const libint2::Shell& b,
        Eigen::MatrixXd& ints);
  /**
   * @brief Computes a set (operator/derivative) of integrals for the given shells.
   * @param op     The kernel/operator as libint enum.
   * @param deriv  The derivative level.
   * @param a      The first shell.
   * @param b      The second shell.
   * @param c      The third shell.
   * @param ints   The integrals stored per set (see class description for more information).
   * @return Boolean. True if integrals were calculated, false if they were screened out.
   */
  bool compute(libint2::Operator op,
        unsigned int deriv,
        const libint2::Shell& a,
        const libint2::Shell& b,
        const libint2::Shell& c,
        Eigen::MatrixXd& ints);
  /**
   * @brief Computes a set (operator/derivative) of integrals for the given shells.
   * @param op     The kernel/operator as libint enum.
   * @param deriv  The derivative level.
   * @param a      The first shell.
   * @param b      The second shell.
   * @param c      The third shell.
   * @param d      The fourth shell.
   * @param ints   The integrals stored per set (see class description for more information).
   * @return Boolean. True if integrals were calculated, false if they were screened out.
   */
  bool compute(libint2::Operator op,
        unsigned int deriv,
        const libint2::Shell& a,
        const libint2::Shell& b,
        const libint2::Shell& c,
        const libint2::Shell& d,
        Eigen::MatrixXd& ints);

  /**
   * @brief Shorthand for an entire set of 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis  The basis.
   * @param pointCharges  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @return The entire set if integrals.
   */
  Eigen::MatrixXd compute1eInts(libint2::Operator op,
                                std::shared_ptr<BasisController> basis,
                                const std::vector<std::pair<double,std::array<double,3>>> pointCharges);

  /**
   * @brief Shorthand for an entire set of 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis  The basis.
   * @param pointCharges  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @return The entire set if integrals.
   */
  Eigen::MatrixXd compute1eInts(libint2::Operator op,
                                std::shared_ptr<BasisController> basis,
                                const std::vector<std::shared_ptr<Atom> >& atoms =
                                     std::vector<std::shared_ptr<Atom> >(0));
  /**
   * @brief Shorthand for an entire set of 'interaction' 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis1 The first basis.
   * @param basis2 The sacond basis.
   * @param atoms  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @return The entire set of 'interaction' integrals.
   */
  Eigen::MatrixXd compute1eInts(libint2::Operator op,
                                std::shared_ptr<BasisController> basis1,
                                std::shared_ptr<BasisController> basis2,
                                const std::vector<std::pair<double,std::array<double,3>>> pointCharges);

  /**
   * @brief Shorthand for an entire set of 'interaction' 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis1 The first basis.
   * @param basis2 The second basis.
   * @param atoms  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @return The entire set of 'interaction' integrals.
   */
  Eigen::MatrixXd compute1eInts(libint2::Operator op,
                                std::shared_ptr<BasisController> basis1,
                                std::shared_ptr<BasisController> basis2,
                                const std::vector<std::shared_ptr<Atom> >& atoms =
                                     std::vector<std::shared_ptr<Atom> >(0));

private:
  /**
   * @brief Private destructor -> singleton.
   */
  Libint();
  /**
   * @brief Map for all the engines.
   */
  std::map<IntegralType ,std::vector<std::unique_ptr<libint2::Engine> > > _engines;
  /**
   * @brief Map for all the engines.
   */
  std::map<IntegralType ,unsigned int > _keep;
  /**
   * @brief Help variable for parallelization critical steps.
   */
  std::mutex _lock;
  /**
   * Assuming all used basis sets have <= 20 primitives per basis function.
   * If more are used just increase this number; Libint objects will then
   * only take slightly more memory.
   */
  const unsigned int N_PRIM_MAX = 20;

};

} /* namespace Serenity */
#endif	/* LIBINT_H */
