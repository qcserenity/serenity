/**
 * @file   OneElectronIntegralController.h
 *
 * @date   30. Dezember 2013, 19:51
 * @author Thomas Dresselhaus
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
#ifndef ONEELECTRONINTEGRALCONTROLLER_H
#define ONEELECTRONINTEGRALCONTROLLER_H
/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h"
#include "geometry/Point.h"
#include "notification/ObjectSensitiveClass.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class Basis;
class BasisController;
class Geometry;
class ExternalChargeController;
template<Options::SCF_MODES>
class MatrixInBasis;
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
   * @brief Constructor
   * @param basisController            The reference basis.
   * @param geometry                   The reference geometry, needed for the nuclei-electron attraction integrals.
   * @param externalChargeController   The external charges.
   */
  OneElectronIntegralController(std::shared_ptr<BasisController> basisController, std::shared_ptr<const Geometry> geometry,
                                std::shared_ptr<ExternalChargeController> externalChargeController);

  // Not defaulted to avoid includes
  virtual ~OneElectronIntegralController();
  /**
   * Forget about the (probably outdated) held integrals
   */
  void clearOneInts();

  void notify() override final {
    this->clearOneInts();
  }

  /* Getters */
  /**
   *
   * @brief If not yet present, the necessary integrals are calculated on call.
   *
   * @returns the matrix typically denoted with h. Sum of kinetic and potential (nuclei-electron
   *          attraction) energy integrals.
   */
  const MatrixInBasis<RESTRICTED>& getOneElectronIntegrals() {
    if (!_oneElectronIntegralsTotal) {
      this->calcHCoreIntegrals();
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
      this->calcOverlapIntegrals();
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
      this->calcECPIntegrals();
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
   * @brief If not yet present, the necessary integrals are calculated on call.
   *
   * @return Returns the kinetic integrals.
   */
  const MatrixInBasis<RESTRICTED>& getKinIntegrals() {
    if (!_kinIntegrals) {
      calcKinIntegrals();
    }
    return *_kinIntegrals;
  }

  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   *
   * @return Returns the nuclear-electron interaction integrals.
   */
  const MatrixInBasis<RESTRICTED>& getNucIntegrals() {
    if (!_nucIntegrals) {
      calcNucIntegrals();
    }
    return *_nucIntegrals;
  }

  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   *
   * @return Returns the external charge integrals.
   */
  const MatrixInBasis<RESTRICTED>& getExtChargeIntegrals();

  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   * @return Returns the electric dipole integrals in length representation of the form <mu|-r|nu>
   *         for a given set of basis functions mu and nu.\n
   *         Note that the negative sign is included to be consistent with response theory.
   */
  const std::vector<MatrixInBasis<RESTRICTED>>& getDipoleLengths(Point gaugeOrigin = Point(0, 0, 0)) {
    if (!_diplen || !gaugeOrigin.isSamePoint(_gaugeOrigin, 1e-6)) {
      this->calcDipoleLengths(gaugeOrigin);
      this->calcDipoleVelocities(gaugeOrigin);
      this->calcDipoleMagnetics(gaugeOrigin);
      _gaugeOrigin = gaugeOrigin;
    }
    return *_diplen;
  }

  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   * @return Returns the dipole integrals in velocity representation of the form <mu|-p|nu>
   *         for a given set of basis functions mu and nu.\n
   */
  const std::vector<MatrixInBasis<RESTRICTED>>& getDipoleVelocities(Point gaugeOrigin = Point(0, 0, 0)) {
    if (!_dipvel || !gaugeOrigin.isSamePoint(_gaugeOrigin, 1e-6)) {
      this->calcDipoleLengths(gaugeOrigin);
      this->calcDipoleVelocities(gaugeOrigin);
      this->calcDipoleMagnetics(gaugeOrigin);
      _gaugeOrigin = gaugeOrigin;
    }
    return *_dipvel;
  }

  /**
   * @brief If not yet present, the necessary integrals are calculated on call.
   * @return Returns the angular momentum integrals of the form 0.5 <mu|L|nu>
   *         for a given set of basis functions mu and nu.\n
   *         Note that the factor of 0.5 is included to be consistent with response theory.
   */
  const std::vector<MatrixInBasis<RESTRICTED>>& getDipoleMagnetics(Point gaugeOrigin = Point(0, 0, 0)) {
    if (!_angmom || !gaugeOrigin.isSamePoint(_gaugeOrigin, 1e-6)) {
      this->calcDipoleLengths(gaugeOrigin);
      this->calcDipoleVelocities(gaugeOrigin);
      this->calcDipoleMagnetics(gaugeOrigin);
      _gaugeOrigin = gaugeOrigin;
    }
    return *_angmom;
  }

  /**
   * @returns the Controller of the Basis for which integrals are calculated.
   */
  std::shared_ptr<BasisController> getBasisController() const {
    return _basisController;
  }
  /**
   * @brief Getter for the external charges if any. If none are available, an empty vector is returned.
   * @return The vector of external charges.
   */
  const std::vector<std::pair<double, Point>>& getExternalCharges();
  /**
   * @brief Setter to cache integrals over external charges.
   * @param integrals The integrals.
   */
  void cacheExtChargeIntegrals(const MatrixInBasis<RESTRICTED>& integrals);

 private:
  const std::shared_ptr<BasisController> _basisController;
  const std::shared_ptr<const Geometry> _geometry;
  Point _gaugeOrigin = Point(0, 0, 0);

  std::unique_ptr<MatrixInBasis<RESTRICTED>> _overlapIntegrals;
  std::unique_ptr<MatrixInBasis<RESTRICTED>> _oneElectronIntegrals;
  std::unique_ptr<MatrixInBasis<RESTRICTED>> _oneElectronIntegralsTotal;
  std::unique_ptr<MatrixInBasis<RESTRICTED>> _kinIntegrals;
  std::unique_ptr<MatrixInBasis<RESTRICTED>> _nucIntegrals;
  std::unique_ptr<MatrixInBasis<RESTRICTED>> _ecpIntegrals;
  std::unique_ptr<MatrixInBasis<RESTRICTED>> _extChargeIntegrals;

  // Dipole integrals, vector has 3 elements for x, y, z components
  std::unique_ptr<std::vector<MatrixInBasis<RESTRICTED>>> _diplen;
  std::unique_ptr<std::vector<MatrixInBasis<RESTRICTED>>> _dipvel;
  std::unique_ptr<std::vector<MatrixInBasis<RESTRICTED>>> _angmom;

  std::shared_ptr<ExternalChargeController> _externalChargeController;

  bool _calcECPs;

  void calcHCoreIntegrals();
  void calcKinIntegrals();
  void calcNucIntegrals();
  void calcECPIntegrals();
  void calcOverlapIntegrals();
  void calcExternalChargeIntegrals();

  void calcDipoleLengths(Point gaugeOrigin);
  void calcDipoleVelocities(Point gaugeOrigin);
  void calcDipoleMagnetics(Point gaugeOrigin);
};

} /* namespace Serenity */
#endif /* ONEELECTRONINTEGRALCONTROLLER_H */
