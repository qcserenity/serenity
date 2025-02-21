/**
 * @file   ElectronicStructure.h
 *
 * @date   last rework Nov 29. 2016
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
#ifndef ELECTRONICSTRUCTURE_H_
#define ELECTRONICSTRUCTURE_H_
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/FockMatrix.h"
#include "energies/EnergyComponentController.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory> //smart ptr.

namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class OrbitalController;
template<Options::SCF_MODES SCFMode>
class NAddFuncPotential;
class OneElectronIntegralController;
template<Options::SCF_MODES SCFMode>
class PotentialBundle;
class Geometry;
class ExternalChargeController;

enum class ES_STATE { INITIAL = 1, GUESS = 2, CONVERGED = 3, FAILED = 4, OTHER = 5 };

/**
 * @class ElectronicStructure ElectronicStructure.h
 * @brief The data structure for results of quantum chemical calculations
 *
 * An electronic structure (in this program) consists of the following parts.
 * The electron density matrix, and the (molecular) orbitals.
 * Implicitly it thus also contains a basis that it is defined in.
 * Furthermore it contains the Fock matrix, and the energies associated with
 * the potentials held in the Fock matrix.
 *
 * For technical reasons the electronic structure also consists of a controller
 * for the one electron integrals which implies that knowledge of the geometry
 * and again the basis is present.
 *
 */
template<Options::SCF_MODES SCFMode>
class ElectronicStructure {
 public:
  /* ===============================
   *   Constructors and Destructor
   * =============================== */

  /**
   * @brief Constructor. Initialization without any orbital/density information.
   * @param oneEIntController The one electron integrals.
   * @param nOccupiedOrbitals The number of occupied orbitals.
   * @param nCoreElectrons    The number of core orbitals.
   */
  ElectronicStructure(std::shared_ptr<OneElectronIntegralController> oneEIntController,
                      const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals,
                      const SpinPolarizedData<SCFMode, unsigned int> nCoreElectrons);
  /**
   * @brief Constructor. Initialization without any orbital/density information. The one electron integral controller
   *        is constructed on the fly.
   * @param basisController   The basis controller.
   * @param geometry          The geometry.
   * @param nOccupiedOrbitals The number of occupied orbitals.
   * @param nCoreElectrons    The number of core orbitals.
   */
  ElectronicStructure(std::shared_ptr<BasisController> basisController, std::shared_ptr<const Geometry> geometry,
                      std::shared_ptr<ExternalChargeController> externalCharges,
                      const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals,
                      const SpinPolarizedData<SCFMode, unsigned int> nCoreElectrons);
  /**
   * @brief Constructor. Initialization using existing orbitals. Orbitals and density will be directly available.
   * @param molecularOrbitals  The molecular orbitals.
   * @param oneEIntController  The one electron integrals.
   * @param nOccupiedOrbitals  The number of occupied orbitals.
   */
  ElectronicStructure(std::shared_ptr<OrbitalController<SCFMode>> molecularOrbitals,
                      std::shared_ptr<OneElectronIntegralController> oneEIntController,
                      const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals);
  /**
   * @brief Copy constructor/conversion constructor between RESTRICTED and UNRESTRICTED.
   * @param other The other electronic structure.
   */
  ElectronicStructure(std::shared_ptr<ElectronicStructure<RESTRICTED>> other);
  /**
   * @brief Copy constructor/conversion constructor between RESTRICTED and UNRESTRICTED.
   * @param other The other electronic structure.
   */
  ElectronicStructure(std::shared_ptr<ElectronicStructure<UNRESTRICTED>> other);
  /**
   * @brief Constructor from HDF5 file.
   * @param fBaseName The basename of the HDF5 files.
   * @param basis The BasisController of the current System.
   * @param geometry The geometry.
   * @param externalCharges The external charges.
   * @param id The string id identifier.
   */
  ElectronicStructure(std::string fBaseName, std::shared_ptr<BasisController> basis, std::shared_ptr<const Geometry> geometry,
                      std::shared_ptr<ExternalChargeController> externalCharges, std::string id);

  /// @brief Default destructor.
  virtual ~ElectronicStructure() = default;

  /* ===============================
   *      Density and Orbitals
   * =============================== */

  /// @returns Returns the density matrix
  DensityMatrix<SCFMode> getDensityMatrix() const {
    return std::move(_densityMatrixController->getDensityMatrix());
  }
  /// @returns the densityMatrixController
  const std::shared_ptr<DensityMatrixController<SCFMode>> getDensityMatrixController() const {
    return _densityMatrixController;
  }
  /// @returns the orbitals
  inline std::shared_ptr<OrbitalController<SCFMode>> getMolecularOrbitals() const {
    return _molecularOrbitals;
  }

  /// @returns the orbitals
  inline void setMolecularOrbitals(std::shared_ptr<OrbitalController<SCFMode>> newOrbitals) {
    _molecularOrbitals = newOrbitals;
    _densityMatrixController->attachOrbitals(_molecularOrbitals, _densityMatrixController->getOccupations());
  }

  /* ===============================
   *             Energies
   * =============================== */

  /// @returns the energy
  inline double getEnergy() const {
    return _energyComponentController->getTotalEnergy();
  }
  /// @returns the energy
  inline double getEnergy(ENERGY_CONTRIBUTIONS energyType) const {
    return _energyComponentController->getEnergyComponent(energyType);
  }
  /// @returns true if all the components for the energy exist
  inline bool checkEnergy(ENERGY_CONTRIBUTIONS energyType) const {
    return _energyComponentController->checkEnergyComponentFromChildren(energyType);
  }
  /// @returns the energyController
  inline std::shared_ptr<EnergyComponentController> getEnergyComponentController() const {
    return _energyComponentController;
  }

  /// @returns The non-additive kinetic potential (used for top-down potential reconstruction apporoach)
  inline std::shared_ptr<NAddFuncPotential<SCFMode>> getNaddKinPotential() const {
    return _naddKinPotential;
  }

  /// @brief Sets the non-additive kinetic potential (used for top-down potential reconstruction apporoach)
  inline void setNaddKinPotential(std::shared_ptr<NAddFuncPotential<SCFMode>> naddKinPotential) {
    _naddKinPotential = naddKinPotential;
  }

  /* ==============================
   *   Potentials and Fock Matrix
   * ============================== */

  /**
   * @brief Attaches a set of potentials to the electronic structure.
   * @param potentialBundle The set of potentials this electronic structure was generated with.
   */
  inline void attachPotentials(std::shared_ptr<PotentialBundle<SCFMode>> potentialBundle) {
    this->_potentials = potentialBundle;
  }
  ///@return Returns true if potentials are available.
  inline bool potentialsAvailable() {
    return (bool)_potentials;
  }
  /**
   * @brief Getter for the potentials used to generate this electronic structure.
   *
   * ElectronicStructures read from disk do not contain this information, thus
   *   potentials are not available.
   * Check for the existence using potentialsAvailable();
   *
   * @return Returns the potentials used to generate the FockMatrix.
   */
  inline std::shared_ptr<PotentialBundle<SCFMode>> getPotentialBundle() {
    if (!_potentials)
      throw SerenityError("No potentials available in the electronic structure.");
    return _potentials;
  }
  ///@brief whether a fock matrix exists
  bool checkFock() {
    bool exists = false;
    if (_fockMatrix != nullptr)
      exists = true;
    return exists;
  }
  ///@brief sets a Fock matrix
  void setFockMatrix(FockMatrix<SCFMode>& fock);
  ///@brief reads the Fock matrix from disk if it exists
  void fockFromHDF5(std::string fBaseName, std::string id);
  /**
   * @brief Getter for the Fock matrix.
   *
   * The getter will attempt to use attached potentials to generate
   *   a new Fock matrix. Please use potentialsAvailable() to check if this is possible
   *   before calling this function.
   *
   * @return The Fock matrix requested.
   */
  FockMatrix<SCFMode> getFockMatrix();
  /// @returns Returns the controller of the overlap and H_core integrals.
  inline std::shared_ptr<OneElectronIntegralController> getOneElectronIntegralController() const {
    return _oneEIntController;
  }

  /* ==============================
   *      I/O and Storage Mode
   * ============================== */

  /**
   * @brief Diskmode getter.
   * @return Returns the diskmode, which is true the data of this electronic structure is manly kept on disk.
   */
  inline bool getDiskMode() const {
    return _diskmode;
  }
  /**
   * @brief Sets the mode for the electronic structure data storage.
   *
   * Before setting this to disk mode (value: false) the file ID and file path have
   * to be set via setFilePathInfo().
   *
   * @param diskmode     Iff false keeps data in memory, else writes and reads to disk.
   * @param fBaseName    The base name for HDF5 files. Usually the systems path and systems name.
   *                     (This can be a dummy string is diskmode is set to be off.)
   * @param id           The system id. (This can be a dummy string is diskmode is set to be off.)
   */
  void setDiskMode(bool diskmode, std::string fBaseName, std::string id);
  /**
   * @brief Prints MO energies and occupations.
   */
  void printMOEnergies() const;
  /**
   * @brief Saves the ElectronicStructure to file.
   * @param fBaseName The basename of the HDF5 files.
   */
  void toHDF5(std::string fBaseName, std::string id);
  /**
   * @brief Getter for the number of occupied orbitals.
   * @return The number of occupied orbitals.
   */
  const SpinPolarizedData<SCFMode, unsigned int>& getNOccupiedOrbitals() {
    return _nOccupiedOrbitals;
  }

  ES_STATE state;

 private:
  bool _diskmode;
  std::shared_ptr<OneElectronIntegralController> _oneEIntController;
  SpinPolarizedData<SCFMode, unsigned int> _nOccupiedOrbitals;
  std::shared_ptr<OrbitalController<SCFMode>> _molecularOrbitals;
  std::shared_ptr<DensityMatrixController<SCFMode>> _densityMatrixController;
  std::shared_ptr<EnergyComponentController> _energyComponentController;
  std::shared_ptr<PotentialBundle<SCFMode>> _potentials;
  std::shared_ptr<NAddFuncPotential<SCFMode>> _naddKinPotential;
  std::shared_ptr<FockMatrix<SCFMode>> _fockMatrix;
  std::string _fBaseName = "";
  std::string _id = "";
};

} /* namespace Serenity */
#endif /* ELECTRONICSTRUCTURE_H_ */
