/**
 * @file LRSCFController.h
 *
 * @date Dec 06, 2018
 * @author Michael Boeckers, Niklas Niemeyer, Johannes Toelle
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

#ifndef LRSCF_LRSCFCONTROLLER
#define LRSCF_LRSCFCONTROLLER

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "settings/BasisOptions.h"

namespace Serenity {

class SystemController;
class GridController;
namespace Options {
enum class LRSCF_TYPE;
enum class LR_METHOD;
} // namespace Options

template<Options::SCF_MODES SCFMode>
class CC2Controller;

template<Options::SCF_MODES SCFMode>
class RIIntegrals;

struct Settings;

struct LRSCFTaskSettings;

/**
 * @class LRSCFController LRSCFController.h
 * @brief The LRSCFController holds the information about a particular\n
 *        active subsystem in the LRSCFTask (e.g.: coefficients, occupied-, virtual orbitals ...).\n
 *        Here, orbital spaces are truncated based on the input in the LRSCFTask. This changes
 *        CoefficientMatrix, occupied orbitals, virtual orbitals and eigenvalues etc.., which are
 *        controlled by this controller. \n
 */
template<Options::SCF_MODES SCFMode>
class LRSCFController : public std::enable_shared_from_this<LRSCFController<SCFMode>> {
 public:
  /**
   * @brief Constructor
   *
   * @param system The system associated with the LRSCFController.
   * @param settings The LRSCFSettings.
   */
  LRSCFController(std::shared_ptr<SystemController> system, const LRSCFTaskSettings& settings);

  /**
   * @brief Default destructor.
   */
  virtual ~LRSCFController() = default;

  /**
   * @brief Returns the number of occupied orbitals stored and modified in the LRSCFController.
   * @return Number of occupied orbitals stored and modified in the LRSCFController.
   */
  SpinPolarizedData<SCFMode, unsigned int> getNOccupied();

  /**
   * @brief Returns the number of virtual orbitals stored and modified in the LRSCFController.
   * @return Number of virtual orbitals stored and modified in the LRSCFController.
   */
  SpinPolarizedData<SCFMode, unsigned int> getNVirtual();

  /**
   * @brief Sets the number of occupied orbitals stored
   */
  void setNOccupied(SpinPolarizedData<SCFMode, unsigned int> nOcc);

  /**
   * @brief Sets the number of virtual orbitals stored
   */
  void setNVirtual(SpinPolarizedData<SCFMode, unsigned int> nVirt);

  /**
   * @brief Returns the reference CoefficientMatrix stored in the LRSCFController.
   * @return Reference CoefficientMatrix stored in the LRSCFController.
   */
  CoefficientMatrix<SCFMode>& getCoefficients();

  /**
   * @brief Returns the reference singles-transformed particle CoefficientMatrix stored in the LRSCFController.
   * @return Reference singles-transformed particle CoefficientMatrix stored in the LRSCFController.
   */
  CoefficientMatrix<SCFMode>& getParticleCoefficients();

  /**
   * @brief Returns the reference singles-transformed hole CoefficientMatrix stored in the LRSCFController.
   * @return Reference singles-transformed hole CoefficientMatrix stored in the LRSCFController.
   */
  CoefficientMatrix<SCFMode>& getHoleCoefficients();

  /**
   * @brief Sets the reference CoefficientMatrix stored in the LRSCFController.
   */
  void setCoefficients(CoefficientMatrix<SCFMode> coeff);

  /**
   * @brief Returns the grid controller (with the system's grid settings, not the LRSCFTask grid settings!).
   * @return Grid Controller.
   */
  std::shared_ptr<GridController> getGridController();

  /**
   * @brief Returns the corresponding eigenvalues to the reference orbitals (can be modified in the LRSCFController).
   * anton: what if the reference orbitals are not canonical? then they don't have an energy
   * @return Corresponding eigenvalues to the reference orbitals.
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> getEigenvalues();

  /**
   * @brief Set the corresponding eigenvalues to the reference orbitals
   */
  void setEigenvalues(SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenvalues);

  /**
   * @brief Returns the controller of the AO basis which is used to express the MOs.
   * @return Controller of the AO basis which is used to express the MOs.
   */
  std::shared_ptr<BasisController> getBasisController(Options::BASIS_PURPOSES basisPurpose = Options::BASIS_PURPOSES::DEFAULT);

  /**
   * @brief Returns the MO Fock matrix of the reference orbitals.
   * @return MO Fock matrix of the reference orbitals.
   */
  std::shared_ptr<SPMatrix<SCFMode>> getMOFockMatrix();

  /**
   * @brief Returns true, if the MO Fock matrix is diagonal or not.
   * @return MO Fock matrix diagonal or not.
   */
  bool isMOFockMatrixDiagonal();

  /**
   * @brief Getter for system settings.
   * @return Settings of the reference system.
   */
  const Settings& getSysSettings();

  /**
   * @brief Returns underlying system controller.
   *
   * Note: Controls objects, which might not work in combination with LRSCF routines
   *       (e.g. when restricting the MO space). Use with care!
   *
   * @return Reference system.
   */
  std::shared_ptr<SystemController> getSys();

  /**
   * @brief Returns the LRSCFTaskSettings.
   * @return Settings of the LRSCFTask this controller was created in.
   */
  const LRSCFTaskSettings& getLRSCFSettings();

  /**
   * @brief Returns the underlying theoretical method for this LRSCFController (i.e. TDDFT or ADC2).
   * @return Method enum for this controller.
   */
  Options::LR_METHOD getResponseMethod();

  /**
   * @brief Returns the eigenvectors corresponding to the excitation energies.
   * @param type The type of the underlying response calculation (iso, uncoupled or coupled).
   * @return Eigenvectors corresponding to the excitation energies.
   */
  std::shared_ptr<std::vector<Eigen::MatrixXd>> getExcitationVectors(Options::LRSCF_TYPE type);

  /**
   * @brief Returns the excitation energies corresponding to the eigenvectors.
   * @param type The type of the underlying response calculation (iso, uncoupled or coupled).
   * @return Excitation energies corresponding to the eigenvectors.
   */
  std::shared_ptr<Eigen::VectorXd> getExcitationEnergies(Options::LRSCF_TYPE type);

  /**
   * @brief Sets the solution of the response problem.
   * @param eigenvectors The eigenvectors ..
   * @param eigenvalues .. and corresponding eigenvalues.
   * @param type The type of the underlying response calculation (iso, uncoupled or coupled).
   */
  void setSolution(std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors,
                   std::shared_ptr<Eigen::VectorXd> eigenvalues, Options::LRSCF_TYPE type);

  /**
   * @brief Sets the environment systems; Knowledge about the environment systems, Important for EOSigmaVector
   * @param envSystems The environment systems.
   */
  void setEnvSystems(std::vector<std::shared_ptr<SystemController>> envSystems);

  /**
   * @brief Returns the environment systems.
   * @return The environment systems.
   */
  std::vector<std::shared_ptr<SystemController>> getEnvSystems();

  /**
   * @brief Edit reference for exclude projection where the individual calculation is carried out in LRSCFTask.
   * @param indexWhiteList New reference orbital indices.
   */
  void editReference(SpinPolarizedData<SCFMode, std::vector<unsigned int>> indexWhiteList);

  /**
   * @brief Initializes an CC2Controller object to this LRSCFController.
   */
  void initializeCC2Controller();

  /**
   * @brief Returns the underlying CC2Controller calculator.
   * @return A pointer to the underlyng CC2Controller calculator.
   * @brief Finalizes the CC2Controller attached to this LRSCFController.
   */
  void finalizeCC2Controller();

  /**
   * @brief Returns the underlying CC2Controller calculator.
   * @return A pointer to the underlyng CC2Controller calculator.
   */
  std::shared_ptr<CC2Controller<SCFMode>> getCC2Controller();

  /**
   * @brief Setup RI integral cache.
   */
  void initializeRIIntegrals(LIBINT_OPERATOR op, double mu, bool calcJia);

  /**
   * @brief Finalizes the RIIntegrals attached to this LRSCFController.
   */
  void finalizeRIIntegrals(LIBINT_OPERATOR op);

  /**
   * @brief The RI integral cache.
   * @return RI integral cache.
   */
  std::shared_ptr<RIIntegrals<SCFMode>> getRIIntegrals(LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb);

  /**
   * @brief The RPA screening used for BSE calculations.
   * @param eia The virtual-occupied orbital energy differences.
   */
  void calculateScreening(const Eigen::VectorXd& eia);

  /**
   * @brief Returns the RPA screening matrix in auxiliary basis.
   * @return The RPA screening matrix in auxiliary basis.
   */
  std::shared_ptr<Eigen::MatrixXd> getScreeningAuxMatrix();

  /**
   * @brief Returns the environment screening.
   * @return The environment approximate RPA screening.
   */
  std::shared_ptr<Eigen::MatrixXd> getEnvTrafo();

  /**
   * @brief Sets the cached inverse aux metric.
   */
  void setInverseMetric(std::shared_ptr<Eigen::MatrixXd> metric);

  /**
   * @brief Returns the cached inverse aux metric.
   * @return The cached inverse aux metric.
   */
  std::shared_ptr<Eigen::MatrixXd> getInverseMetric();

  /**
   * @brief Sets the cached inverse erf aux metric.
   */
  void setInverseErfMetric(std::shared_ptr<Eigen::MatrixXd> metric);

  /**
   * @brief Returns the cached inverse erf aux metric.
   * @return The cached inverse aux metric.
   */
  std::shared_ptr<Eigen::MatrixXd> getInverseErfMetric();

  /**
   * @brief Applies frozen core approximation.
   */
  void applyFrozenCore();

  /**
   * @brief Applies frozen virtual approximation (based on energy difference with HOMO).
   */
  void applyFrozenVirtual();

  /**
   * @brief Applies core only approximation.
   */
  void applyCoreOnly();

  /**
   * @brief Sets up this LRSCF Controller to be the reference for a Spin-Flip calculation.
   */
  void setupSpinFlipReference();

  /**
   * @brief Rotates the reference orbitals according to the found instabilites.
   *
   */
  void rotateOrbitalsSCFInstability();

  std::shared_ptr<std::vector<MatrixInBasis<SCFMode>>> getUnrelaxedDiffDensities() const {
    return _unrelaxedDiffDensities;
  }

  std::shared_ptr<std::vector<MatrixInBasis<SCFMode>>> getRelaxedDiffDensities() {
    return _relaxedDiffDensities;
  }

 private:
  // The system controller
  std::shared_ptr<SystemController> _system;
  // User defined LRSCF settings
  const LRSCFTaskSettings& _settings;
  // The environment systems
  std::vector<std::shared_ptr<SystemController>> _envSystems;
  // Number of occupied orbitals
  SpinPolarizedData<SCFMode, unsigned int> _nOcc;
  // Number of virtual orbitals
  SpinPolarizedData<SCFMode, unsigned int> _nVirt;
  // Fock matrix in MO basis
  std::shared_ptr<SPMatrix<SCFMode>> _fock;
  // List of orbital indices included in LRSCF calculation
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> _indexWhiteList;
  // The grid used for numerical integration
  std::shared_ptr<GridController> _grid;
  // A set of reference orbitals
  CoefficientMatrix<SCFMode> _coefficients;
  // A set of reference orbitals
  CoefficientMatrix<SCFMode> _particleCoefficients;
  // A set of reference orbitals
  CoefficientMatrix<SCFMode> _holeCoefficients;
  // Corresponding orbital energies
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _orbitalEnergies;
  // X (_excitationVectors[0]) and Y (_excitationVectors[1]) excitation vectors.
  std::shared_ptr<std::vector<Eigen::MatrixXd>> _excitationVectors;

  std::shared_ptr<std::vector<MatrixInBasis<SCFMode>>> _unrelaxedDiffDensities =
      std::make_shared<std::vector<MatrixInBasis<SCFMode>>>();
  std::shared_ptr<std::vector<MatrixInBasis<SCFMode>>> _relaxedDiffDensities =
      std::make_shared<std::vector<MatrixInBasis<SCFMode>>>();
  // Type of _excitationVectors (isolated,uncoupled,coupled)
  Options::LRSCF_TYPE _type;
  // Corresponding excitation energies
  std::shared_ptr<Eigen::VectorXd> _excitationEnergies;
  // load excitation vectors and excitation energies from HDF5
  void loadFromH5(Options::LRSCF_TYPE type);
  // the underlying CC2Controller
  std::shared_ptr<CC2Controller<SCFMode>> _CC2Controller;
  // the conventional ri integrals associated with this controller
  std::shared_ptr<RIIntegrals<SCFMode>> _riints;
  // the erf ri integrals associated with this controller
  std::shared_ptr<RIIntegrals<SCFMode>> _riErfints;
  // screening matrix for BSE
  std::shared_ptr<Eigen::MatrixXd> _screening;
  // environment screening contribution BSE
  std::shared_ptr<Eigen::MatrixXd> _envTransformation;

  std::shared_ptr<Eigen::MatrixXd> _inverseM;
  std::shared_ptr<Eigen::MatrixXd> _inverseErfM;

}; // class LRSCFController
} // namespace Serenity

#endif /* LRSCF_LRSCFCONTROLLER */
