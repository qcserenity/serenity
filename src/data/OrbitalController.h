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
template<Options::SCF_MODES SCFMode>
class OrbitalController : public NotifyingClass<OrbitalController<SCFMode>>, public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @param coefficientMatrix with data defined for the basis.
   * @param basis             for which the orbitals in coefficientMatrix are defined.
   * @param eigenvalues       the orbital energies.
   * @param isCoreOrbital     Flag for core orbitals.
   */
  OrbitalController(std::unique_ptr<CoefficientMatrix<SCFMode>> coefficients, std::shared_ptr<BasisController> basisController,
                    std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>> eigenvalues,
                    std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXi>> isCoreOrbital);
  /**
   * @param isCoreOrbital     Flag for core orbitals.
   * @param coefficients      with data defined for the basis.
   * @param basisController   for which the orbitals in coefficients are defined.
   * @param eigenvalues       the orbital energies.
   */
  OrbitalController(std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXi>> isCoreOrbital,
                    std::unique_ptr<CoefficientMatrix<SCFMode>> coefficients, std::shared_ptr<BasisController> basisController,
                    std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>> eigenvalues);
  /**
   * @param coefficientMatrix with data defined for the basis
   * @param basis             for which the orbitals in coefficientMatrix are defined.
   * @param eigenvalues       the orbital energies
   * @param nCoreElectrons    The number of core electrons (Assigns the core orbitals by eigenvalue).
   */
  OrbitalController(std::unique_ptr<CoefficientMatrix<SCFMode>> coefficients, std::shared_ptr<BasisController> basisController,
                    const SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues, unsigned int nCoreElectrons);
  /**
   * @brief provides an empty set of orbitals waiting to be filled
   * @param basis for which the orbitals in coefficientMatrix are defined.
   * @param nCoreOrbitals The number of core orbitals (Assigns the first n/2 orbitals as core).
   */
  explicit OrbitalController(std::shared_ptr<BasisController> basisController,
                             const SpinPolarizedData<SCFMode, unsigned int> nCoreElectrons);
  /**
   * @param orig Explicit copy constructor.
   *     (Is explicit to avoid unintentional copies.)
   */
  explicit OrbitalController(const OrbitalController<SCFMode>& orig);
  /**
   * @brief Constructor from HDF5 file
   * @param filePath The HDF5 file containing the data
   * @param basisController The basisController of the running system.
   */
  OrbitalController(std::string filePath, std::shared_ptr<BasisController> basisController, std::string id);
  ///@brief Move constructor
  OrbitalController(OrbitalController<SCFMode>&&) = default;
  ///@brief Default destructor.
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
  CoefficientMatrix<SCFMode> getCoefficients();

  /**
   * @brief Updates the orbitals using custom coefficients and eigenvalues.
   * @param updatedCoefficients The new coefficients.
   * @param updatedEigenvalues  The new orbital energies.
   * @param coreOrbitals        The core orbital flags.
   */
  void updateOrbitals(const CoefficientMatrix<SCFMode>& updatedCoefficients,
                      const SpinPolarizedData<SCFMode, Eigen::VectorXd>& updatedEigenvalues,
                      SpinPolarizedData<SCFMode, Eigen::VectorXi> coreOrbitals);

  /**
   * @brief Updates the orbitals using custom coefficients and eigenvalues.
   * @param updatedCoefficients The new coefficients.
   * @param updatedEigenvalues  The new orbital energies.
   */
  void updateOrbitals(const CoefficientMatrix<SCFMode>& updatedCoefficients,
                      const SpinPolarizedData<SCFMode, Eigen::VectorXd>& updatedEigenvalues);

  /**
   * @brief Resorts the coefficients and eigenvalues acording to the MOM/IMOM procedures.
   * @param c  The new coefficients.
   * @param eps   The new orbital energies.
   * @param momMatrix     Represents occupied orbitals according to MOM/IMOM procedure.
   * @param overlapMatrix The overlap integrals.
   */
  void applyMOMProcedure(CoefficientMatrix<SCFMode>& c, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eps,
                         const SPMatrix<SCFMode> momMatrix, const MatrixInBasis<RESTRICTED> overlapMatrix);

  /// @returns the orbital energies.
  SpinPolarizedData<SCFMode, Eigen::VectorXd> getEigenvalues();
  /// @returns Returns the core/Rydberg orbital flags: 1: Core orbital, 2: Rydberg orbital, 0: valence orbital
  SpinPolarizedData<SCFMode, Eigen::VectorXi> getOrbitalFlags();
  /// @brief Returns the number of core orbitals.
  SpinPolarizedData<SCFMode, unsigned int> getNCoreOrbitals();
  /**
   * @brief Set all orbitals with an orbital eigenvalue lower than the cut-off as core orbitals.
   * @param energyCutOff The energy cut-off
   */
  void setCoreOrbitalsByEnergyCutOff(double energyCutOff);
  /**
   * @brief Set virtual orbitals to valence and Rydberg-like orbitals by energy cut-off.
   * @param energyCutOff The energy cut-off. Orbitals with eigenvalues larger than the
   *                     given threshold are considered as Rydberg orbitals.
   */
  void setRydbergOrbitalsByEnergyCutOff(double energyCutOff);
  /**
   * @brief Set the orbitals with lowest eigenvalues as core orbitals.
   * @param nCoreOrbitals The number of core orbitals.
   */
  void setCoreOrbitalsByNumber(unsigned int nCoreOrbitals);
  /**
   * @brief Set virtual orbitals to valence and Rydberg-like orbitals.
   * @param nRydbergOrbitals The number of non-valence orbitals.
   */
  void setRydbergOrbitalsByNumber(unsigned int nRydbergOrbitals);
  /**
   * @brief Set the first N orbitals as core orbitals irrespective of their eigenvalue.
   * @param nCoreOrbitals The number of core orbitals.
   */
  void setCoreOrbitalsFirstN(const SpinPolarizedData<SCFMode, unsigned int>& nCoreOrbitals);
  /**
   * @brief Get the ranges for valence and core orbitals among the set of occupied orbitals.
   * @param nOcc The number of occupied orbitals.
   * @return The ranges for core and valence orbitals.
   */
  std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>>
  getValenceOrbitalIndices(SpinPolarizedData<SCFMode, unsigned int> nOcc);
  /**
   * @brief Get the ranges for virtual valence and Rydberg orbitals.
   * @param nOcc The number of occupied orbitals.
   * @return The ranges for Rydberg and valence orbitals.
   */
  std::pair<SpinPolarizedData<SCFMode, std::vector<unsigned int>>, SpinPolarizedData<SCFMode, std::vector<unsigned int>>>
  getVirtualValenceOrbitalIndices(SpinPolarizedData<SCFMode, unsigned int> nOcc);
  /**
   * @brief Get the range for all valence orbitals.
   * @return The range of all valence orbitals.
   */
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> getAllValenceOrbitalIndices();
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
   * @param momMatrix Represents occupied orbitals according to MOM/IMOM procedure.
   */
  void updateOrbitals(const FockMatrix<SCFMode>& fockMatrix, std::shared_ptr<OneElectronIntegralController> oneIntController,
                      std::shared_ptr<SPMatrix<SCFMode>> momMatrix = nullptr);
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
   *        based on the owned Fock-like matrix and overlap integrals. See function getLevelShift()
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
   * @param momMatrix Represents occupied orbitals according to MOM/IMOM procedure.
   */
  void updateOrbitals(const std::pair<Eigen::VectorXd, SpinPolarizedData<SCFMode, Eigen::VectorXd>> levelshift,
                      const FockMatrix<SCFMode>& fockMatrix, std::shared_ptr<OneElectronIntegralController> oneIntController,
                      std::shared_ptr<SPMatrix<SCFMode>> momMatrix = nullptr);

  ///@brief Notification
  void notify();
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
   * @brief Reads only the core orbital flags from file
   * @param fBaseName The basename of the HDF5 files.
   */
  void coreOrbitalsfromHDF5(std::string fBaseName, std::string id);

  /**
   * @brief Sets the threshold for the canonical orthogonalization.
   * @param threshold The new threshold.
   */
  void setCanOrthTh(const double& threshold) {
    _canOrthTh = threshold;
    _X.resize(0, 0);
  }
  /**
   * @brief Getter for the custom overlap.
   * @return The overlap matrix.
   */
  std::shared_ptr<MatrixInBasis<SCFMode>> getCustomOverlap() {
    return _customS;
  }
  /**
   * @brief Sets a custom overlap to be used.
   *
   * This custom overlap can is allowed to be for alpha and beta orbitals,
   * for this reason a MatrixInBasis is expected.
   *
   * @param S The overlap.
   */
  void useCustomOverlap(const MatrixInBasis<SCFMode> S, bool calcX = false) {
    _customS = std::make_shared<MatrixInBasis<SCFMode>>(S);
    if (calcX)
      calculateCustomTransformationX();
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
   * @brief Getter for the transformation matrix to the basis in which the eigenvalue problem is solved.
   * @return The transformation matrix.
   */
  std::shared_ptr<MatrixInBasis<SCFMode>> getCustomTransformMatrix() {
    return _customX;
  }

  /**
   * @brief Constructs a core orbital vector according to the eigenvalue up to the given number of core electrons.
   * @param nCoreElectrons The number of core electrons.
   * @param eigenvalues    The eigenvalues.
   * @return The core orbital vector.
   */
  static SpinPolarizedData<SCFMode, Eigen::VectorXi>
  getCoreOrbitalsByEigenvalue(unsigned int nCoreElectrons, const SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues);

 private:
  std::unique_ptr<CoefficientMatrix<SCFMode>> _coefficients;
  const std::shared_ptr<BasisController> _basisController;
  std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>> _eigenvalues;
  std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXi>> _orbitalFlags;
  double _canOrthTh;
  bool _firstIteration = true;
  Eigen::MatrixXd _X;
  Eigen::MatrixXd _Xinv;
  bool _linearDependent;
  unsigned int _nZero;
  bool _keepInMemory = true;
  bool _fIsInOthoBasis = false;
  std::shared_ptr<MatrixInBasis<SCFMode>> _customS = nullptr;
  std::shared_ptr<MatrixInBasis<SCFMode>> _customX = nullptr;
  std::shared_ptr<MatrixInBasis<SCFMode>> _customXinv = nullptr;
  std::string _fBaseName;
  std::string _id;
  void calculateTransformationX(std::shared_ptr<OneElectronIntegralController> oneIntController);
  void calculateCustomTransformationX();
};

} /* namespace Serenity */
#endif /* ORBITALCONTROLLER_H_ */
