/**
 * @file   Libint.h
 *
 * @date   28. Juli 2013, 14:12
 * @author Thomas Dresselhaus, Jan Unsleber
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
#ifndef LIBINT_H
#define LIBINT_H

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace libint2 {
enum class Operator;
class Engine;
struct Shell;
} // namespace libint2

namespace Serenity {

/* Forward Declarations */
class BasisController;
class Atom;
class Point;
class ShellPairData;

enum class LIBINT_OPERATOR {
  overlap,
  kinetic,
  nuclear,
  emultipole1,
  emultipole2,
  delta,
  coulomb,
  cgtg,
  cgtg_x_coulomb,
  delcgtg2,
  erf_coulomb,
  erfc_coulomb
};

/**
 * @class IntegralType
 * @brief A small class to wrap all the possible integral types.
 */
struct IntegralType {
  /**
   * @brief Default constructor.
   */
  IntegralType() : op(), deriv(), nCenter(){};
  /**
   * @brief Constructor.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   */
  IntegralType(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter)
    : op(op), deriv(deriv), nCenter(nCenter - 2) {
  }
  /// @brief The kernel/operator as libint enum.
  LIBINT_OPERATOR op;
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
 * Currently interfaced: libint-2.7.0-beta6
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
 *     0->[ints],1->[dx,dy,dz],2->[dxdx,dxdy,dxdz,dydy,dydz,dzdz],... .
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
  static Libint& getInstance() {
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
   * @param mu           Range separation parameter if operator is erf_coulomb.
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   */
  void initialize(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter,
                  const std::vector<std::pair<double, std::array<double, 3>>>& pointCharges, double mu,
                  double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                  unsigned int maxNPrim = N_PRIM_MAX);

  /**
   * @brief Initializes the engines for a specific type of integral.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   * @param pointCharges The charges to be added to the engines
   *                     (needed for Operator::nuclear integrals).
   * @param mu           Range separation parameter if operator is erf_coulomb.
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   */
  void initialize(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter,
                  const std::vector<std::pair<double, Point>>& pointCharges, double mu = 0.0,
                  double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                  unsigned int maxNPrim = N_PRIM_MAX);

  /**
   * @brief Initializes the engines for a specific type of integral.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   * @param atoms The atoms to be added to the engines
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @param mu           Range separation parameter if operator is erf_coulomb.
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   */
  void initialize(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter,
                  const std::vector<std::shared_ptr<Atom>>& atoms = std::vector<std::shared_ptr<Atom>>(0),
                  double mu = 0.0, double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                  unsigned int maxNPrim = N_PRIM_MAX);

  /**
   * @brief Initializes the engines for a specific type of integral.
   * @param op The kernel/operator as libint enum.
   * @param nCenter The number of centers (2-4).
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   */
  void initialize_plain(LIBINT_OPERATOR op, const unsigned int nCenter,
                        double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                        unsigned int maxNPrim = N_PRIM_MAX);
  /**
   * @brief Initializes the engines for a specific type of integral.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   * @param multipoleOrigin The gauge-origin for multipole integrals.
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   */
  void initialize(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter, const Point multipoleOrigi,
                  double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                  unsigned int maxNPrim = N_PRIM_MAX);

  /**
   * @brief Cleans the engines for a specific type of integral.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   */
  void finalize(LIBINT_OPERATOR op, unsigned int deriv, unsigned int nCenter);

  /**
   * @brief Keeps a set of engines form being destroyed by finalize calls.
   *
   * Every Libint::keepEngines has to be matched by a Libint::freeEngines .
   *
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   */
  void keepEngines(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter);

  /**
   * @brief Allows the destruction and destroys a set of engines if they are present.
   * @param op The kernel/operator as libint enum.
   * @param deriv The derivative level.
   * @param nCenter The number of centers (2-4).
   */
  void freeEngines(LIBINT_OPERATOR op, const unsigned int deriv, const unsigned int nCenter);

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
   * @param normAux Flag for normalization of integrals corresponding to auxiliary basis functions.
   * @return Boolean. True if integrals were calculated, false if they were screened out.
   */
  bool compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
               Eigen::MatrixXd& ints, bool normAux = true);
  /**
   * @brief Computes a set (operator/derivative) of integrals for the given shells.
   * @param op     The kernel/operator as libint enum.
   * @param deriv  The derivative level.
   * @param a      The first shell.
   * @param b      The second shell.
   * @param ints   The integrals stored per set (see class description for more information).
   * @param normAux Flag for normalization of integrals corresponding to auxiliary basis functions.
   * @return Boolean. True if integrals were calculated, false if they were screened out.
   */
  bool compute(LIBINT_OPERATOR op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
               Eigen::MatrixXd& ints, bool normAux = true);
  /**
   * @brief Computes a set (operator/derivative) of integrals for the given shells.
   * @param op     The kernel/operator as libint enum.
   * @param deriv  The derivative level.
   * @param a      The first shell.
   * @param b      The second shell.
   * @param c      The third shell.
   * @param ints   The integrals stored per set (see class description for more information).
   * @param normAux Flag for normalization of integrals corresponding to auxiliary basis functions.
   * @return Boolean. True if integrals were calculated, false if they were screened out.
   */
  bool compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
               const libint2::Shell& c, Eigen::MatrixXd& ints, bool normAux = true);
  /**
   * @brief Computes a set (operator/derivative) of integrals for the given shells.
   * @param op     The kernel/operator as libint enum.
   * @param deriv  The derivative level.
   * @param a      The first shell.
   * @param b      The second shell.
   * @param c      The third shell.
   * @param ints   The integrals stored per set (see class description for more information).
   * @param normAux Flag for normalization of integrals corresponding to auxiliary basis functions.
   * @return Boolean. True if integrals were calculated, false if they were screened out.
   */
  bool compute(LIBINT_OPERATOR op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
               const libint2::Shell& c, Eigen::MatrixXd& ints, bool normAux = true);
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
  bool compute(libint2::Operator op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
               const libint2::Shell& c, const libint2::Shell& d, Eigen::MatrixXd& ints);
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
  bool compute(LIBINT_OPERATOR op, unsigned int deriv, const libint2::Shell& a, const libint2::Shell& b,
               const libint2::Shell& c, const libint2::Shell& d, Eigen::MatrixXd& ints);

  /**
   * @brief Shorthand for an entire set of 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis  The basis.
   * @param pointCharges  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   * @param shellPairData The shell pairs to calculate the integrals for. If not provided, all integrals are calculated.
   * @return The integrals.
   */
  Eigen::MatrixXd compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis,
                                const std::vector<std::pair<double, std::array<double, 3>>> pointCharges,
                                double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                                unsigned int maxNPrim = N_PRIM_MAX,
                                std::shared_ptr<std::vector<ShellPairData>> shellPairData = nullptr);

  /**
   * @brief Shorthand for an entire set of 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis  The basis.
   * @param pointCharges  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   * @return The integrals.
   */
  Eigen::MatrixXd compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis,
                                const std::vector<std::pair<double, Point>>& pointCharges,
                                double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                                unsigned int maxNPrim = N_PRIM_MAX,
                                std::shared_ptr<std::vector<ShellPairData>> shellPairData = nullptr);

  /**
   * @brief Shorthand for an entire set of 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis  The basis.
   * @param pointCharges  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   * @return The integrals.
   */
  Eigen::MatrixXd compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis,
                                const std::vector<std::shared_ptr<Atom>>& atoms = std::vector<std::shared_ptr<Atom>>(0),
                                double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                                unsigned int maxNPrim = N_PRIM_MAX);
  /**
   * @brief Shorthand for an entire set of 'interaction' 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis1 The first basis.
   * @param basis2 The second basis.
   * @param atoms  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   * @return The entire set of 'interaction' integrals.
   */
  Eigen::MatrixXd compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis1,
                                std::shared_ptr<BasisController> basis2,
                                const std::vector<std::pair<double, std::array<double, 3>>> pointCharges,
                                double precision = std::numeric_limits<double>::epsilon(), double maxD = 10,
                                unsigned int maxNPrim = N_PRIM_MAX);

  /**
   * @brief Shorthand for an entire set of 'interaction' 1e integrals.
   * @param op     The kernel/operator as libint enum.
   * @param basis1 The first basis.
   * @param basis2 The second basis.
   * @param atoms  The atoms to be added to the engine
   *                     (needed for Operator::nuclear integrals).
   *                     [default=None]
   * @param mu Range separation parameter if operator is erf_coulomb.
   * @param precision    The precision requested for the fully contracted integrals
   *                     (as linear combination of primitives).
   * @param maxD         Maximum coefficient to be contracted with the integrals.
   * @param maxNPrim     Maximum number of primitives in any basis function shell.
   * @return The entire set of 'interaction' integrals.
   */
  Eigen::MatrixXd compute1eInts(LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis1,
                                std::shared_ptr<BasisController> basis2,
                                const std::vector<std::shared_ptr<Atom>>& atoms = std::vector<std::shared_ptr<Atom>>(0),
                                double mu = 0.0, double precision = std::numeric_limits<double>::epsilon(),
                                double maxD = 10, unsigned int maxNPrim = N_PRIM_MAX);

  /**
   * @brief Engines for 2e-4c integrals.
   * @param The kernel/operator as libint enum.
   * @return A reference to the libint engines for 2e-4c integrals.
   */
  std::vector<std::unique_ptr<libint2::Engine>>& getFourCenterEngines(LIBINT_OPERATOR op);

  /**
   * @brief  Resolve the internal operator to its libint2 representation.
   * @param op The operator.
   */
  static libint2::Operator resolveLibintOperator(LIBINT_OPERATOR op);
  /**
   * @brief  Resolve the internal operator to its Serenity representation.
   * @param op The operator.
   */
  static LIBINT_OPERATOR resolveLibintOperator(libint2::Operator op);
  /**
   * @brief Getter for the maximum number of primitive functions in a shell.
   * @return The maximum number of primitive functions in a shell.
   */
  static unsigned int getNPrimMax();

 private:
  /**
   * @brief Private destructor -> singleton.
   */
  Libint();
  /**
   * @brief Map for all the engines.
   */
  std::map<IntegralType, std::vector<std::unique_ptr<libint2::Engine>>> _engines;
  /**
   * @brief Map for all the engines.
   */
  std::map<IntegralType, unsigned int> _keep;
  /**
   * @brief Help variable for parallelization critical steps.
   */
  std::mutex _lock;
  /**
   * Assuming all used basis sets have <= 23 primitives per basis function.
   * If more are used just increase this number; Libint objects will then
   * only take slightly more memory.
   */
  static const unsigned int N_PRIM_MAX = 23;
  /**
   * Libint will truncate the primitive functions at some point. Thus, we have to make sure that
   * the error introduced like this in the final contracted integral set does not exceed our required
   * precision.
   * This function returns the precision needed for the primitive function evaluation.
   * The current implementation is rather conservative.
   */
  long double getFinalPrecision(double precision, double maxD, double nPrim, unsigned int nCenters);
};

} /* namespace Serenity */
#endif /* LIBINT_H */
