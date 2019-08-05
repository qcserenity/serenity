/**
 * @file   OneElectronIntegralController.h
 *
 * @date   30. Dezember 2013, 19:51
 * @author Thomas Dresselhaus
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
#ifndef ONEELECTRONINTEGRALCONTROLLER_H
#define	ONEELECTRONINTEGRALCONTROLLER_H
/* Include Serenity Internal Headers */
#include "notification/ObjectSensitiveClass.h"
#include "settings/Options.h"
#include "data/matrices/MatrixInBasis.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class Basis;
class BasisController;
class Geometry;
template <Options::SCF_MODES> class MatrixInBasis;
/**
 * @class OneElectronIntegralController OneElectronIntegralController.h
 *
 * Objects of this class provide and store integrals over two basis functions,
 * namely the overlap and one-electron integrals. Integrals are calculated lazy,
 * i.e. only on the first request.
 */
class OneElectronIntegralController : public ObjectSensitiveClass<Basis> {
public:
  /**
   * @param oneIntCalculator
   * @param basis
   * @param geometry         A Geometry fitting to the Basis must be specified for the nuclei-
   *                         electron attraction integrals.
   */
  OneElectronIntegralController(
      std::shared_ptr<BasisController> basisController,
      std::shared_ptr<const Geometry> geometry);
  // Not defaulted to avoid includes
  virtual ~OneElectronIntegralController();
  /**
   * Forget about the (probably outdated) held integrals
   */
  void clearOneInts();

  void notify() override final {
    this->clearOneInts();
  }
  /**
   * @brief Direct access to an overlap integral (slow)
   *
   * @param i,j (Extended) Indices of the basis functions
   */
  double S(unsigned int i, unsigned int j);
  /**
   * @brief Direct access to a one-electron integral (slow)
   *
   * , i.e. kinetic plus potential (nuclei-electron attraction) energy integrals.
   *
   * @param i,j (Extended) Indices of the basis functions
   */
  double h(unsigned int i, unsigned int j);
  /* Getters */
  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   *
   * @returns the matrix typically denoted with h. Sum of kinetic and potential (nuclei-electron
   *          attraction) energy integrals.
   */
  const MatrixInBasis<RESTRICTED>& getOneElectronIntegrals() {
    if (!_oneElectronIntegralsTotal) {
      calcHCoreIntegrals();
    }
    return *_oneElectronIntegralsTotal;
  }
  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   *
   * @returns the overlap matrix; typically denoted with S.
   */
  const MatrixInBasis<RESTRICTED>& getOverlapIntegrals() {
    if (!_overlapIntegrals) {
      calcOverlapIntegrals();
    }
    return *_overlapIntegrals;
  }
  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   *
   * @return Returns the ECP integrals.
   */
  const MatrixInBasis<RESTRICTED>& getECPIntegrals() {
    if (!_ecpIntegrals) {
      calcHCoreIntegrals();
    }
    return *_ecpIntegrals;
  }
  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   *
   * @return Returns the sum of kinetic and nuclear integrals.
   */
  const MatrixInBasis<RESTRICTED>& getNucKinIntegrals() {
    if (!_oneElectronIntegrals) {
      calcHCoreIntegrals();
    }
    return *_oneElectronIntegrals;
  }
  /**
   * @returns the Controller of the Basis for which integrals are calculated.
   */
  std::shared_ptr<BasisController> getBasisController() const {
    return _basisController;
  }
  /* Print methods */
  /**
   * @brief Prints to stdout
   */
  void printOneElectronIntegrals();
  /**
   * @brief Prints to stdout
   */
  void printOverlapIntegrals();

private:
  const std::shared_ptr<BasisController> _basisController;
  const std::shared_ptr<const Geometry> _geometry;

  std::unique_ptr<MatrixInBasis<RESTRICTED> > _overlapIntegrals;
  std::unique_ptr<MatrixInBasis<RESTRICTED> > _oneElectronIntegrals;
  std::unique_ptr<MatrixInBasis<RESTRICTED> > _oneElectronIntegralsTotal;
  std::unique_ptr<MatrixInBasis<RESTRICTED> > _ecpIntegrals;

  bool _calcECPs;

  void calcHCoreIntegrals();
  void calcOverlapIntegrals();
  void cleanMem();
};

} /* namespace Serenity */
#endif	/* ONEELECTRONINTEGRALCONTROLLER_H */
