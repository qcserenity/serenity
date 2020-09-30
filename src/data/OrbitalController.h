/**
 * @file   OrbitalController.h
 *
 * @date   last rework Jan 16. 2017
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
#ifndef ORBITALCONTROLLER_H_
#define ORBITALCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "notification/NotifyingClass.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <iostream>
#include <memory>
#include <vector>

namespace Serenity {
class OneElectronIntegralController;
/**
 * @class OrbitalController OrbitalController.h
 * @brief holds a completely defined set of orbitals in analytical form.
 */
template<Options::SCF_MODES T>
class OrbitalController : public NotifyingClass<OrbitalController<T>>, public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @param coefficientMatrix with data defined for the basis
   * @param basis             for which the orbitals in coefficientMatrix are defined.
   * @param eigenvalues       the orbital energies
   */
  OrbitalController(std::unique_ptr<CoefficientMatrix<T>> coefficients, std::shared_ptr<BasisController> basisController,
                    std::unique_ptr<SpinPolarizedData<T, Eigen::VectorXd>> eigenvalues);
  /**
   * @brief provides an empty set of orbitals waiting to be filled
   * @param basis for which the orbitals in coefficientMatrix are defined.
   */
  explicit OrbitalController(std::shared_ptr<BasisController> basisController);
  /**
   * @param orig Explicit copy constructor.
   *     (Is explicit to avoid unintentional copies.)
   */
  explicit OrbitalController(const OrbitalController<T>& orig);
  /**
   * @brief Constructor from HDF5 file
   * @param filePath The HDF5 file containing the data
   * @param basisController The basisController of the running system.
   */
  OrbitalController(std::string filePath, std::shared_ptr<BasisController> basisController, std::string id);
  ///@brief Move constructor
  OrbitalController(OrbitalController<T>&&) = default;
  virtual ~OrbitalController();

  /// @returns the controller for the underlying basis
  std::shared_ptr<BasisController> getBasisController() const {
    return _basisController;
  }

  /**
   * @brief Sets the mode for the orbital data storage.
   *
   * @param diskmode     Iff false keeps data in memory, else writes and reads to disk.
   * @param fBaseName    The base name for HDF5 files. Usually the systems path and systems name.
   *                     (This can be a dummy string is diskmode is set to be off.)
   * @param id           The system id. (This can be a dummy string is diskmode is set to be off.)
   */
  void setDiskMode(bool diskmode, std::string fBaseName, std::string id);

  /// @returns the coefficients determining the orbitals in connection with the basis.
  CoefficientMatrix<T> getCoefficients();

  /**
   * @brief Updates the orbitals using custom coefficients and eigenvalues.
   * @param updatedCoefficients The new coefficients.
   * @param updatedEigenvalues The new orbital energies.
   */
  void updateOrbitals(const CoefficientMatrix<T>& updatedCoefficients,
                      const SpinPolarizedData<T, Eigen::VectorXd>& updatedEigenvalues);

  /// @returns the orbital energies.
  SpinPolarizedData<T, Eigen::VectorXd> getEigenvalues();
  /**
   * @returns the number of orbitals (occupied AND virtual) in this OrbitalController
   */
  unsigned int getNOrbitals() const;

  /**
   * @brief Erases and recreates the owned molecular orbital coefficients and orbital energies.
   *
   * Does this based on the owned Fock-like matrix and overlap integrals.
   *
   * @param fockMatrix The Fock matrix to update the orbitals with.
   * @param oneIntController One electron integrals.
   *                        (The overlap is needed to transform the Fock matrix into an orthogonal
   *                         basis.)
   */
  void updateOrbitals(const FockMatrix<T>& fockMatrix, std::shared_ptr<OneElectronIntegralController> oneIntController);
  /**
   * @brief Erases and recreates the owned molecular orbital coefficients and orbital energies.
   *        Adds a levelshift:
   *        The level shift increases the orbital energies of the virtual orbitals. The effect of
   *        this can be understood by perturbation theory. If we write the new coefficients as sum
   *        of the old coefficients plus a contribution from the virtual orbitals,
   *        \f[
   *          c_{vi}^{(n+1)} = c_{vi}^{(n)} + \sum_a d_{ia}^{(n+1)} c_{va}^{(n)} \; ,
   *        \f]
   *        insert this into the Fock equations and multiply by \f$ c_{\mu a}^* \f$ we obtain
   *        an equation for the expansion coefficients
   *        \f[
   *          d_{ia}^{(n+1)} = - \frac{F_{ai}^{(n)}}{\epsilon_a - \epsilon_i} \; .
   *        \f]
   *        The orbital energy \f$\epsilon_a \f$ is replaced by \f$\epsilon_a = \epsilon_a + b \f$
   *        with \f$ b > 0\f$. This will reduce the changes in the MO coefficients and the tendency
   *        for the MO coefficients to oscillate is reduced. While damping reduces the step of the
   *        SCF iteration for all orbital equally, the level shift reduces the occupied-virtual mixing
   *        more for pairs with small energy differences.
   *        \n
   *        The level shift can be fixed or determined dynamically. Here, we determine the level shift
   *        based on the owned Fock-like matrix and overlap integrals
   *        \f[
   *          b = 0.1 \sqrt{log(d_{conv} + 2)}
   *        \f]
   *        where \f$ d_{conv} \f$ is a scalar measure for convergence. See function getLevelShift()
   *        in the ConvergenceController.cpp.
   *        \n
   *        Ref.: H. B. Schlegel and J. J. W. McDouall, Do you have SCF Stability and Convergence
   *              Problems? in Computational Advances in Organic Chemistry: Molecular Structure and Reactivity,
   *              167-185, 1991
   *
   * @param levelshift std::pair of the levelshift energy in Hartree and the number of electrons.
   * @param fockMatrix The Fock matrix to update the orbitals with.
   * @param oneIntController One electron integrals.
   *                        (The overlap is needed to transform the Fock matrix into an orthogonal
   *                         basis.)
   */
  void updateOrbitals(const std::pair<Eigen::VectorXd, SpinPolarizedData<T, Eigen::VectorXd>> levelshift,
                      const FockMatrix<T>& fockMatrix, std::shared_ptr<OneElectronIntegralController> oneIntController);

  ///@brief Notification
  void notify() {
    _X.resize(0, 0);
    _firstIteration = true;
  }

  /**
   * @brief Triggers the notification.
   */
  void externalNotifyObjects() {
    this->notifyObjects();
  }
  /**
   * @brief Saves an OrbitalController to file.
   * @param fBaseName The basename of the HDF5 files.
   */
  void toHDF5(std::string fBaseName, std::string id);
  /**
   * @brief Reads an OrbitalController from file
   * @param fBaseName The basename of the HDF5 files.
   */
  void fromHDF5(std::string fBaseName, std::string id);

  /**
   * @brief Reads only coefficients from file
   * @param fBaseName The basename of the HDF5 files.
   */
  void coefficientsfromHDF5(std::string fBaseName, std::string id);

  /**
   * @brief Reads only eigenvalues from file
   * @param fBaseName The basename of the HDF5 files.
   */
  void eigenvaluesfromHDF5(std::string fBaseName, std::string id);

  /**
   * @brief Sets the threshold for the canonical orthogonalization.
   * @param threshold The new threshold.
   */
  void setCanOrthTh(const double& threshold) {
    _canOrthTh = threshold;
    _X.resize(0, 0);
  }

  /**
   * @brief Sets iff F is expected to be in an orthogonal basis [default: false].
   * @param value The new value.
   */
  void setExpectedOthogonality(bool value) {
    _fIsInOthoBasis = value;
  }

  /**
   * @brief Sets a custom overlap to be used.
   *
   * This custom overlap can is allowed to be for alpha and beta orbitals,
   * for this reason a MatrixInBasis is expected.
   *
   * @param S The overlap.
   */
  void useCustomOverlap(const MatrixInBasis<T> S) {
    _customS.reset(new MatrixInBasis<T>(S));
  }

  /**
   * @brief Check if orbitals are already calculated.
   * @return Returns true iff orbitals are available.
   */
  bool orbitalsAvailable() {
    return !_firstIteration;
  }

  /**
   * @brief Getter for the transformation matrix to the basis in which the eigenvalue problem is solved.
   * @param oneIntController The OneElectronIntegralController of the system.
   * @return The transformation matrix.
   */
  std::shared_ptr<Eigen::MatrixXd> getTransformMatrix(std::shared_ptr<OneElectronIntegralController> oneIntController) {
    if (!(_X.cols() > 0))
      calculateTransformationX(oneIntController);
    assert(_X.cols() > 0);
    return std::make_shared<Eigen::MatrixXd>(_X);
  }
  /**
   * @brief Getter for the inverse transformation matrix to the canonical
   *        basis.
   */
  std::shared_ptr<Eigen::MatrixXd> getTransformMatrixInverse(std::shared_ptr<OneElectronIntegralController> oneIntController) {
    if (!(_Xinv.cols() > 0))
      calculateTransformationX(oneIntController);
    assert(_X.cols() > 0);
    return std::make_shared<Eigen::MatrixXd>(_Xinv);
  }

 private:
  std::unique_ptr<CoefficientMatrix<T>> _coefficients;
  const std::shared_ptr<BasisController> _basisController;
  std::unique_ptr<SpinPolarizedData<T, Eigen::VectorXd>> _eigenvalues;
  double _canOrthTh;
  std::shared_ptr<OneElectronIntegralController> _oneIntController;
  bool _firstIteration = true;
  Eigen::MatrixXd _X;
  Eigen::MatrixXd _Xinv;
  bool _linearDependent;
  unsigned int _nZero;
  bool _keepInMemory = true;
  bool _fIsInOthoBasis = false;
  std::unique_ptr<MatrixInBasis<T>> _customS = nullptr;
  std::string _fBaseName;
  std::string _id;
  void calculateTransformationX(std::shared_ptr<OneElectronIntegralController> oneIntController);
};

} /* namespace Serenity */
#endif /* ORBITALCONTROLLER_H_ */
