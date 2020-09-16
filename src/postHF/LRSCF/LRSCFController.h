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
#include "settings/Options.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

namespace Options {
enum class LRSCF_TYPE;
}

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
class LRSCFController {
 public:
  /**
   * @brief Constructor
   *
   * @param system The system associated with the LRSCFController.
   * @param settings The LRSCFSettings.
   */
  LRSCFController(std::shared_ptr<SystemController> system, LRSCFTaskSettings& settings);

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
   * @brief Returns an occupation vector stored and modified in the LRSCFController.
   * @return Occupation vector stored and modified in the LRSCFController.
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd>& getOccupation();

  /**
   * @brief Returns the reference CoefficientMatrix stored in the LRSCFController.
   * @return Reference CoefficientMatrix stored in the LRSCFController.
   */
  CoefficientMatrix<SCFMode>& getCoefficients();

  /**
   * @brief Returns the grid controller.
   * @return Grid Controller.
   */
  std::shared_ptr<GridController> getGridController();

  /**
   * @brief Returns the corresponding eigenvalues to the reference orbitals (can be modified in the LRSCFController).
   * @return Corresponding eigenvalues to the reference orbitals.
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> getEigenvalues();

  /**
   * @brief Returns the controller of the AO basis which is used to express the MOs.
   * @return Controller of the AO basis which is used to express the MOs.
   */
  std::shared_ptr<BasisController> getBasisController();

  /**
   * @brief Returns the Fock matrix of the non-canonical reference orbitals.
   * @return Fock matrix of the non-canonical reference orbitals.
   */
  std::shared_ptr<MatrixInBasis<SCFMode>> getFockNonCanon();

  /**
   * @brief Set the Fock matrix for non-canonical orbitals.
   * @param fockNonCanon Fock matrix of the non-canonical reference orbitals.
   */
  void setFockNonCanon(std::shared_ptr<MatrixInBasis<SCFMode>> fockNonCanon);

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
   * @return Settings of the LRSCFTask this Controller was created in.
   */
  const LRSCFTaskSettings& getLRSCFSettings();

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
   * @brief Restricts orbital space according to the user's input.
   */
  void editReference();

  /**
   * @brief Edit reference for exclude projection where the individual calculation is carried out in LRSCFTask.
   * @param indexWhiteList New reference orbital indices.
   * @param system The reference system.
   * @param type The type of the underlying response calculation (iso, uncoupled or coupled).
   */
  void editReference(SpinPolarizedData<SCFMode, std::vector<unsigned int>> indexWhiteList, unsigned int system,
                     Options::LRSCF_TYPE type);

  /**
   * @brief Returns the (spin-polarized) list of orbials indices included in the calculation.
   * @return The (spin-polarized) list of orbials indices included in the calculation.
   */
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> getReferenceOrbitals();

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

 private:
  // The system controller
  std::shared_ptr<SystemController> _system;
  // User defined LRSCF settings
  LRSCFTaskSettings& _settings;
  // The environment systems
  std::vector<std::shared_ptr<SystemController>> _envSystems;
  // A occupation vector for the set of reference orbitals
  SpinPolarizedData<SCFMode, Eigen::VectorXi> _occupation;
  // Non canocical fock matrix
  std::shared_ptr<MatrixInBasis<SCFMode>> _fockNonCanon;
  // List of orbital indices included in LRSCF calculation
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> _indexWhiteList;
  // The grid used for numerical integration
  std::shared_ptr<GridController> _grid;
  // A set of reference orbitals
  CoefficientMatrix<SCFMode> _coefficients;
  // Corresponding orbital energies
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _orbitalEnergies;
  // X (_excitationVectors[0]) and Y (_excitationVectors[1]) excitation vectors.
  std::shared_ptr<std::vector<Eigen::MatrixXd>> _excitationVectors;
  // Type of _excitationVectors (isolated,uncoupled,coupled)
  Options::LRSCF_TYPE _type;
  // Corresponding excitation energies
  std::shared_ptr<Eigen::VectorXd> _excitationEnergies;
  // Cast user input into SpinPolarizedData to be used in a for_spin loop
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> getSet();
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> getExclude();

  void loadFromH5(Options::LRSCF_TYPE type);

  // Write indexwhitelist to H5 file
  void indexToH5(Options::LRSCF_TYPE type);

}; // class LRSCFController
} // namespace Serenity

#endif /* LRSCF_LRSCFCONTROLLER */
