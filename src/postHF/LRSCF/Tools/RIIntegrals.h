/**
 * @file RIIntegrals.h
 *
 * @date Feb 14, 2020
 * @author Niklas Niemeyer
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

#ifndef LRSCF_RIINTEGRALS
#define LRSCF_RIINTEGRALS

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "integrals/wrappers/Libint.h"
#include "settings/BasisOptions.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/LRSCFOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

enum class SCF_MODES;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

class SystemController;
class TwoElecThreeCenterIntLooper;
class TwoElecThreeCenterCalculator;
class Geometry;
class BasisController;

/**
 * @class RIIntegrals RIIntegrals.h
 * @brief Performs the calculation of the delta energy sigma vectors
 */
template<Options::SCF_MODES SCFMode>
class RIIntegrals {
 public:
  /**
   * @brief Constructor
   * @param lrscf The LRSCFController.
   * @param op The libint operator.
   * @param mu The range-separation parameter.
   * @param calcJia Are (ia|Q) integrals also requested?
   * @param pStart Start MO index for (pq|Q) integrals.
   * @param pEnd End MO index for (pq|Q) integrals.
   * @param geo Custom geometry for auxiliary basis controller.
   */
  RIIntegrals(std::shared_ptr<LRSCFController<SCFMode>> lrscf, LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb, double mu = 0.0,
              bool calcJia = true, unsigned pStart = 0, unsigned pEnd = 0, std::shared_ptr<Geometry> geo = nullptr);

  /**
   * @brief Constructor
   * @param sys The SystemController.
   * @param op The libint operator.
   * @param mu The range-separation parameter.
   * @param calcJia Are (ia|Q) integrals also requested?
   * @param pStart Start MO index for (pq|Q) integrals.
   * @param pEnd End MO index for (pq|Q) integrals.
   * @param nafThresh NAF Threshold.
   * @param geo Custom geometry for auxiliary basis controller.
   */
  RIIntegrals(std::shared_ptr<SystemController> sys, LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb, double mu = 0.0,
              bool calcJia = true, unsigned pStart = 0, unsigned pEnd = 0, double nafThresh = 0.0,
              std::shared_ptr<Geometry> geo = nullptr, Options::DENS_FITS densFitCache = Options::DENS_FITS::RI);

  /**
   * @brief Destructor.
   */
  virtual ~RIIntegrals() = default;

  /**
   * @brief Triggers AO integral caching.
   */
  void cacheAOIntegrals();

  /**
   * @brief Returns the number of occupied orbitals.
   * @return The number of occupied orbitals.
   */
  SpinPolarizedData<SCFMode, unsigned> getNOccupied();

  /**
   * @brief Returns the number of virtual orbitals.
   * @return The number of virtual orbitals.
   */
  SpinPolarizedData<SCFMode, unsigned> getNVirtual();

  /**
   * @brief Returns the number of basis functions.
   * @return The number of basis functions.
   */
  unsigned getNBasisFunctions();

  /**
   * @brief Returns the number of aux-transformed auxiliary basis functions.
   * @return The number of aux-transformed auxiliary basis functions.
   */
  unsigned getNTransformedAuxBasisFunctions();

  /**
   * @brief Returns the number of actual auxiliary basis functions.
   * @return The number of actual auxiliary basis functions.
   */
  unsigned getNAuxBasisFunctions();

  /**
   * @brief Returns (ij|Q) integral pointer. Already contracted with the auxiliary transformation matrix.
   * @return The (ij|Q) integral pointer.
   */
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> getJijPtr();

  /**
   * @brief Returns (ia|Q) integral pointer. Already contracted with the auxiliary transformation matrix.
   * @return The (ia|Q) integral pointer.
   */
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> getJiaPtr();

  /**
   * @brief Returns (pq|Q) integral pointer. Already contracted with the auxiliary transformation matrix.
   * @return The (pq|Q) integral pointer.
   */
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> getJpqPtr();

  /**
   * @brief Returns the auxiliary transformation matrix [sqrt(V^-1/2) * NAF].
   * @return The auxiliary transformation matrix [sqrt(V^-1/2) * NAF].
   */
  std::shared_ptr<Eigen::MatrixXd> getAuxTrafoPtr();

  /**
   * @brief Returns the auxiliary metric matrix [V].
   * @return The auxiliary transformation matrix [V].
   */
  std::shared_ptr<Eigen::MatrixXd> getAuxMetricPtr();

  /**
   * @brief Returns the 3c-looper for those integrals which are not stored.
   * @return The 3c-looper for those integrals which are not stored.
   */
  std::shared_ptr<TwoElecThreeCenterCalculator> getIntegralPtr();

  /**
   * @brief Returns true if all 3c-MO integrals are stored.
   * @return True if all 3c-MO integrals are stored.
   */
  bool isFullyMOCached();

  /**
   * @brief Returns true if all 3c-AO integrals are stored.
   * @return True if all 3c-AO integrals are stored.
   */
  bool isFullyAOCached();

  /**
   * @brief Removes all stored MO integrals from memory.
   */
  void clearMOCache();

  /**
   * @brief Removes all stored AO integrals from memory.
   */
  void clearAOCache();

  /**
   * @brief Returns start MO index for (pq|Q) integrals.
   * @return The start MO index for (pq|Q) integrals.
   */
  int getPStart();

  /**
   * @brief Returns end MO index for (pq|Q) integrals.
   * @return The end MO index for (pq|Q) integrals.
   */
  int getPEnd();

  /**
   * @brief Returns a pointer to the custom geometry.
   * @return Pointer to the custom geometry.
   */
  std::shared_ptr<Geometry> getGeo();

  /**
   * @brief Sets a pointer to the custom geometry.
   */
  void setGeo(std::shared_ptr<Geometry> geo);

 private:
  ///@brief Initialize function.
  void printInfo();

  ///@brief Calculates integrals.
  void calculateIntegrals();

  ///@brief Applies NAF approximation.
  void applyNAFApproximation();

  ///@brief The LRSCF Controller.
  std::weak_ptr<LRSCFController<SCFMode>> _lrscf;

  ///@brief The kernel/operator as libint enum.
  LIBINT_OPERATOR _op;

  ///@brief range separation parameter mu
  double _mu;

  ///@brief Number of occupied orbitals.
  SpinPolarizedData<SCFMode, unsigned> _no;

  ///@brief Number of virtual orbitals.
  SpinPolarizedData<SCFMode, unsigned> _nv;

  ///@brief Number of basis functions.
  unsigned long _nb;

  ///@brief Number of auxiliary functions of the unmodified auxiliary basis.
  unsigned long _nxb;

  ///@brief Number of actual auxiliary functions (e.g. after NAF approx).
  unsigned long _nx;

  ///@brief Start custom MO index.
  unsigned long _pStart;

  ///@brief End custom MO index.
  unsigned long _pEnd;

  ///@brief Indicate whether (ia|Q) is requested.
  bool _calcJia;

  ///@brief Indicate whether (ia|Q) is requested.
  bool _calcJpq;

  ///@brief Integral lists (ij|Q).
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jij = nullptr;

  ///@brief Integral lists (ia|Q).
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jia = nullptr;

  ///@brief Integral lists (pq|Q).
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jpq = nullptr;

  ///@brief [V^{-1/2}] [* naf-matrix].
  std::shared_ptr<Eigen::MatrixXd> _auxTrafo = nullptr;

  ///@brief [V^{1}] [* naf-matrix].
  std::shared_ptr<Eigen::MatrixXd> _metric = nullptr;

  ///@brief Looper for RI integrals.
  std::shared_ptr<TwoElecThreeCenterCalculator> _integrals = nullptr;

  ///@brief Basis controller.
  std::shared_ptr<BasisController> _basContr = nullptr;

  ///@brief Auxiliary basis controller.
  std::shared_ptr<BasisController> _auxBasContr = nullptr;

  ///@brief Custom geometry for auxiliary basis.
  std::shared_ptr<Geometry> _geo;

  ///@brief Does this thing hold all 3-center MO integrals?
  bool _fullyMOCached;

}; // class RIIntegrals
} // namespace Serenity

#endif /* LRSCF_RIINTEGRALS */
