/**
 * @file   SystemController.h
 * @author Thomas Dresselhaus, Jan Unsleber
 *
 * @date   20. Juli 2015, 16:51
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
#ifndef SYSTEMCONTROLLER_H
#define SYSTEMCONTROLLER_H
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "dft/functionals/CompositeFunctionals.h" //XC Functionals
#include "settings/BasisOptions.h"                //Default arguments.
#include "settings/ElectronicStructureOptions.h"  //RESTRICTED/UNRESTRICTED
#include "settings/GridOptions.h"                 //Default arguments.
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory> //smart ptr.
#include <vector>

namespace Serenity {
/* Forward Declarations */
class Atom;
class System;
class Geometry;
class Point;
class ExternalChargeController;
class BasisController;
class AtomCenteredBasisController;
class OneElectronIntegralController;
class GridController;
class AtomCenteredGridController;
class CDIntegralController;
template<Options::SCF_MODES>
class OrbitalController;
template<Options::SCF_MODES>
class PotentialBundle;
struct Settings;
template<Options::SCF_MODES>
class ElectronicStructure;
enum class MOLECULAR_SURFACE_TYPES;
template<Options::SCF_MODES>
class ScfTask;
class MolecularSurfaceController;
template<Options::SCF_MODES>
class ElectrostaticPotentialOnGridController;
class IntegralCachingController;
struct CUSTOMFUNCTIONAL;
/**
 * @class SystemController SystemController.h
 * @brief A quite complex class managing all data associated with a System
 *
 * where the System is basically defined by a Geometry and charge and spin.
 * This controller provides (or takes in) objects which will be associated with the System, but it
 * does not fully control those. E.g. if the geometry changes it is not this controller that
 * takes care of adapting the basis, integrals and so on. Check the Geometry and Basis class and
 * suitable controllers to find these functionalities.
 *
 * Note: do not construct an additional shared_ptr/unique_ptr of this class, instead use the getSharedPtr() method.
 *       Furthermore, always construct this class via make_shared<SystemController>(...).
 */
class SystemController : public std::enable_shared_from_this<SystemController> {
 public:
  /**
   * @brief Constructor.
   * @param settings
   */
  SystemController(Settings settings);
  /**
   * @brief Constructor.
   * @param geometry
   * @param settings
   */
  SystemController(std::shared_ptr<Geometry> geometry, Settings settings);
  /**
   * @brief Destructor. Set as default in the .cpp.
   */
  ~SystemController();

  /**
   * @brief Adds up two Systems. The new system will have a
   *        simply added geometry, charge and spin. Settings
   *        will be taken from the system on the left hand side
   *        of the operator.
   *        Furthermore, the orbitals and orbital energies are combined,
   *        ordered as follows: lhs_occ,rhs_occ,lhs_virt,rhs_virt.\n
   *
   *        Example:\n
   *        @code
   *        std::shared_ptr<SystemController> sys1;
   *        std::shared_ptr<SystemController> sys2;
   *        std::shared_ptr<SystemController> sys3= *sys1 + *sys2;
   *        @endcode
   * @param rhs The second system.
   * @return A shared pointer on the new system
   */
  std::shared_ptr<SystemController> operator+(SystemController& rhs);

  /**
   * A SystemController holds a reference to the only shared_ptr
   * that should ever be used.
   * Warning: do not create additional shared/unique_ptr
   * @return Returns a shared pointer to the current instance of this class
   */
  inline std::shared_ptr<SystemController> getSharedPtr() {
    return shared_from_this();
  }
  /**
   * @param The new SCF mode [RESTRICTED,UNRESTRICTED].
   */
  void setSCFMode(Options::SCF_MODES mode);
  /**
   * @brief Getter for the system SCFMode.
   * @return The SCFMode.
   */
  Options::SCF_MODES getSCFMode();
  /**
   * @returns the name of the controlled molecular system. It should be unique.
   */
  std::string getSystemName();
  /**
   * @param name the new name.
   */
  void setSystemName(std::string name);
  /**
   * @returns the id of the controlled molecular system.
   */
  std::string getSystemIdentifier();
  /**
   * @returns Returns the path to saved files.
   */
  std::string getSystemPath();
  /**
   * @returns Returns the path to saved files.
   */
  std::string getHDF5BaseName();
  /**
   * @returns the underlying configuration. All data received from this SystemController is
   *          constructed by using this configuration.
   */
  const Settings& getSettings() const;
  /**
   * @brief Based on the charges of the nuclei the number of electrons is adjusted accordingly.
   *
   * @param charge the charge of the whole system. Cannot be greater than the sum of charges of
   *               the nuclei
   */
  void setCharge(int charge);
  /**
   * @brief Sets the systems spin.
   *
   * @param spin the spin of the whole system.
   */
  void setSpin(int spin);
  /**
   * @brief Sets the systems electric field.
   * @param position the field vector
   * @param fStrength the fieldstrength
   * @param analytical the decider if field shall be determined analytically oder numerically
   */
  void setElectricField(std::vector<double> position, double fStrength, bool analytical, bool use);
  /**
   * @brief Sets the DFT XC functional of the system.
   *
   * @param XCfunc the XC functional.
   */
  void setXCfunctional(CompositeFunctionals::XCFUNCTIONALS XCfunc);
  /**
   * @brief Sets the DFT XC functional of the system.
   *
   * @param customfunc the XC functional.
   */
  void setXCfunctional(CUSTOMFUNCTIONAL customfunc);

  /**
   * @brief Reduces memory storage, clearing temporary data, dumping results to disk.
   *
   * This function switched the existing electronic structure to be dumped to disk, and
   *   read/written on each request that is made, also all cached integrals for this system
   *   will be flushed.
   *   This reduces the memory usage significantly for big systems, or a case with many
   *   used systems, however it means reevaluation of many properties/integrals, be sure
   *   that this is what you want when activating this mode.
   *
   * @param diskmode The new disk mode, true equals storage on disk,
   *                  false results in mostly memory storage.
   */
  void setDiskMode(bool diskmode);

  /*
   * Getters
   */
  /// @returns the charge of the molecular system
  int getCharge() const;
  /// @return the number of core electrons/non-valence electrons from the Geometry.
  ///         This is NOT the information stored in the orbital controller!
  unsigned int getNCoreElectrons() const;
  /// @returns the SCF type of the last SCF Calculation, i.e. restricted or unrestricted
  enum Options::SCF_MODES getLastSCFMode() const;
  /// @returns The SCF type assigned in the system settings.
  enum Options::SCF_MODES getSCFMode() const;
  /**
   * @returns The expectation value of the S_z-Operator, i.e. the excess of alpha electrons over
   *          beta electrons. May also be negative, i.e. there may also be more beta electrons than
   *          alpha electrons.
   */
  int getSpin() const;
  /**
   * @brief returns the number of alpha and beta electrons
   * @returns the number of alpha and beta electrons
   */
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int> getNAlphaAndBetaElectrons() const;

  /// @returns true iff the system is open shell (i.e. \< S_z \> != 0).
  bool isOpenShell() const;
  /**
   *  @returns the number of electrons. Note that this number may not be the same as the sum of
   *           the nuclear charges (modified by the total charge of the system) if effective core
   *           potentials are used (the ECPs reduce the effective nuclear charges and thus also
   *           the number of electrons).
   *
   *  See also getNECPElectrons() in Atom.h
   */
  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode, unsigned int> getNElectrons() const;
  /// @returns the number of occupied molecular orbitals
  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode, unsigned int> getNOccupiedOrbitals() const;
  /// @returns the number of virtual molecular orbitals
  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode, unsigned int> getNVirtualOrbitals();
  /// @returns the number of truncated virtual molecular orbitals (only inclusion of specific virtuals)
  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode, unsigned int> getNVirtualOrbitalsTruncated();

  /**
   * @param   basisPurpose what you want to use the basis for. May have different specifications,
   *                       e.g. it may be much larger if you want an auxiliary basis to evaluate ´
   *                       Coulombic interactions or something.
   * @returns the associated basis (for the specified basisPurpose)
   */
  std::shared_ptr<BasisController> getBasisController(Options::BASIS_PURPOSES basisPurpose = Options::BASIS_PURPOSES::DEFAULT) const;
  /**
   * @param   auxBasisPurpose The type of integrals the auxiliary basis is tailored towards.
   * @param   dfMode          The mode of density-fitting used (i.e. RI, ACD, ACCD)
   * @returns the associated auxiliary basis controller (for the specified basisPurpose)
   */
  std::shared_ptr<BasisController> getAuxBasisController(Options::AUX_BASIS_PURPOSES auxBasisPurpose,
                                                         Options::DENS_FITS dfMode) const;
  /**
   * @param   auxBasisPurpose The type of integrals the auxiliary basis is tailored towards.
   * @param   dfMode          The mode of density-fitting used (i.e. RI, ACD, ACCD)
   * @returns the associated basis purpose option (for the specified auxiliary basis purpose)
   */
  Options::BASIS_PURPOSES resolveAuxBasisPurpose(Options::AUX_BASIS_PURPOSES auxBasisPurpose, Options::DENS_FITS dfMode) const;
  /**
   * Currently we only use atom-centered basis sets, but that may change in the future, e.g. if you
   * think about embedded calculations where you have basis functions on atoms this system does not
   * know anything about. Because of that we have this additional function call and interface.
   *
   * @param   basisPurpose what you want to use the basis for. May have different specifications,
   *                       e.g. it may be much larger if you want an auxiliary basis to evaluate ´
   *                       Coulombic interactions or something.
   * @returns the associated basis (for the specified basisPurpose) with the atom-related
   *          information.
   */
  std::shared_ptr<AtomCenteredBasisController>
  getAtomCenteredBasisController(Options::BASIS_PURPOSES basisPurpose = Options::BASIS_PURPOSES::DEFAULT) const;
  /**
   * @param   basisPurpose what you want to use the basis for. May have different specifications,
   *                       e.g. it may be much larger if you want an auxiliary basis to evaluate ´
   *                       Coulombic interactions or something.
   * @returns the used controller for one-electron integrals (for the basis with the specified
   *          basisPurpose)
   */
  std::shared_ptr<OneElectronIntegralController>
  getOneElectronIntegralController(Options::BASIS_PURPOSES basisPurpose = Options::BASIS_PURPOSES::DEFAULT) const;
  /// @returns the underlying geometry
  std::shared_ptr<Geometry> getGeometry() const;
  /**
   * @param gridPurpose In some cases it is useful to have a smaller integration grid than the final
   *                    one. E.g. during the SCF-cycles in DFT calculations (or during geometry
   *                    optimizations) a less accurate integration grid is often sufficient, because
   *                    the total energy is much more sensitive to the integration grid than the
   *                    final orbitals. Thus one can save a lot of computational time when
   *                    intermediately working with a smaller grid. However, a final solution should
   *                    always be created with the 'normal' (i.e. DEFAULT) integration grid.
   * @returns the attached gridController with the specifications for gridPurpose
   */
  std::shared_ptr<GridController> getGridController(Options::GRID_PURPOSES gridPurpose = Options::GRID_PURPOSES::DEFAULT) const;
  /**
   * @param gridController The new grid controller.
   * @param gridPurpose In some cases it is useful to have a smaller integration grid than the final
   *                    one. E.g. during the SCF-cycles in DFT calculations (or during geometry
   *                    optimizations) a less accurate integration grid is often sufficient, because
   *                    the total energy is much more sensitive to the integration grid than the
   *                    final orbitals. Thus one can save a lot of computational time when
   *                    intermediately working with a smaller grid. However, a final solution should
   *                    always be created with the 'normal' (i.e. DEFAULT) integration grid.
   */
  void setGridController(std::shared_ptr<GridController> gridController,
                         Options::GRID_PURPOSES gridPurpose = Options::GRID_PURPOSES::DEFAULT) const;
  /**
   * @brief see the method getGridController(); this method returns an atom-centered version
   */
  std::shared_ptr<AtomCenteredGridController>
  getAtomCenteredGridController(Options::GRID_PURPOSES gridPurpose = Options::GRID_PURPOSES::DEFAULT) const;
  /// @returns the currently active OrbitalController
  template<Options::SCF_MODES SCFMode>
  std::shared_ptr<OrbitalController<SCFMode>> getActiveOrbitalController();
  ;
  /**
   *  @brief Runs a SCF if no active ElectronicStructure is available.
   *  @returns the currently active ElectronicStructure
   */
  template<Options::SCF_MODES SCFMode>
  std::shared_ptr<ElectronicStructure<SCFMode>> getElectronicStructure();
  /**
   *  @brief Checks if an ElectronicStructure is available.
   *  @returns Returns true if an ElectronicStructure is present.
   */
  template<Options::SCF_MODES SCFMode>
  bool hasElectronicStructure();
  /*
   * Forwarded getters
   */
  /// @returns the geometry's atoms
  const std::vector<std::shared_ptr<Atom>>& getAtoms() const;
  /// @returns the geometry's number of atoms
  unsigned int getNAtoms() const;
  /**
   * @brief This function generates a set of potentials for basic HF/DFT SCF runs.
   *
   * If there is no electronic structure present when this function is called then
   * there will be an initial guess generated. The storage of the PotentialBundle inside ElectronicStructure is done in
   * the SCF. This is due to the fact that all of the potentials need reference objects.
   *
   * @return Returns a set of potentials.
   */
  template<Options::SCF_MODES SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES Theory>
  std::shared_ptr<PotentialBundle<SCFMode>> getPotentials(Options::GRID_PURPOSES grid = Options::GRID_PURPOSES::DEFAULT);
  /**
   *  @brief Getter for the Cholesky integral Controller of the system.
   *  @returns The Cholesky integral controller.
   */
  std::shared_ptr<CDIntegralController> getCDIntegralController();
  /*
   * Setters
   */
  /// @param electronicStructure will be used in this System from now on
  template<Options::SCF_MODES SCFMode>
  void setElectronicStructure(std::shared_ptr<ElectronicStructure<SCFMode>>);

  /**
   * @brief Set a new BasisController
   * It is not possible to override an existing basis controller
   * assign it before it is ever needed, or don't assign it!
   * @param basisController A new basis Controller
   * @param basisPurpose    The purpose of the new BasisController
   */
  void setBasisController(std::shared_ptr<AtomCenteredBasisController> basisController,
                          Options::BASIS_PURPOSES basisPurpose = Options::BASIS_PURPOSES::DEFAULT);

  /**
   * @brief Loads an ElectronicStructure from file. System objects
   *        like Geometries, spin and charge are handled in the
   *        Constructor of this class, since they should not change
   *        in the lifetime of this class.
   * @param system A reference to the system.
   */
  void fromHDF5(std::string loadPath);

  /**
   * @brief Prints system information to screen.
   */
  void print();
  /**
   * @brief Setter for a molecular surface/cavity.
   * @param surface       The new surface.
   * @param surfaceType   The surface type (Active or FDE)
   */
  void setMolecularSurface(std::shared_ptr<MolecularSurfaceController> surface, MOLECULAR_SURFACE_TYPES surfaceType);
  /**
   * @brief Getter for the molecular surface.
   * @param surfaceType The surface type.
   * @return The molecular surface.
   */
  std::shared_ptr<MolecularSurfaceController> getMolecularSurface(MOLECULAR_SURFACE_TYPES surfaceType);
  /**
   * @brief Getter for the electrostatic potential on molecular surface controller.
   * @param surfaceType The surface type.
   * @return The electrostatic potential controller.
   */
  template<Options::SCF_MODES SCFMode>
  std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>
  getElectrostaticPotentialOnMolecularSurfaceController(MOLECULAR_SURFACE_TYPES surfaceType);
  /**
   * @brief Setter for the electrostatic potential on molecular surface controller.
   * @param potential    The new electrostatic potential on molecular surface controller.
   * @param surfaceType  The surface type.
   */
  template<Options::SCF_MODES SCFMode>
  void setElectrostaticPotentialOnMolecularSurfaceController(std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>> potential,
                                                             MOLECULAR_SURFACE_TYPES surfaceType);
  /**
   * @brief Set the "use" variable in the system specific settings to newMode
   */
  void setSystemContinuumModelMode(bool newMode);
  /**
   * @brief Forwarded getter for _settings.pcm.use.
   */
  bool getSystemContinuumModelMode();
  /**
   * @brief Setter for the electronic structure method.
   * @param method The method.
   */
  void setElectronicStructureMethod(Options::ELECTRONIC_STRUCTURE_THEORIES method);
  /**
   * @brief Getter for the integral caching controller, i.e. the 4-center integral cache.
   * @return The inegral caching controller.
   */
  std::shared_ptr<IntegralCachingController> getIntegralCachingController();
  /**
   * @brief Clear the integral caching controller.
   */
  void clear4CenterCache();
  /**
   * @return Returns true if there are external point charges to be considered for the system.
   */
  bool hasExternalCharges() const;

  /**
   * @brief Set the gradients of the point charges.
   * @param pointChargeGradients The point charge gradients.
   */
  void setPointChargeGradients(const Eigen::MatrixXd& pointChargeGradients);

  /**
   * @brief Getter for the point charge gradients. May throw an error if none are available.
   * @return The derivative of the energy with respect to the point charge positions.
   */
  const Eigen::MatrixXd& getPointChargeGradients();
  /**
   * @brief Getter for any external charges to the system. These can be set through the settings.
   * @return The external charges.
   */
  std::shared_ptr<ExternalChargeController> getExternalChargeController() const;

 private:
  void produceBasisController(Options::BASIS_PURPOSES basisPurpose) const;
  void produceGridController(Options::GRID_PURPOSES gridPurpose) const;
  void produceMolecularSurface();
  void produceMolecularVanDerWaalsSurface();
  template<Options::SCF_MODES SCFMode>
  void produceElectrostaticPotentialOnMolecularSurfaceController(MOLECULAR_SURFACE_TYPES surfaceType);

  template<Options::SCF_MODES SCFMode>
  void produceScfTask();

  std::unique_ptr<System> _system;

  std::unique_ptr<ScfTask<Options::SCF_MODES::RESTRICTED>> _restrictedScfTask;
  std::unique_ptr<ScfTask<Options::SCF_MODES::UNRESTRICTED>> _unrestrictedScfTask;

  std::shared_ptr<CDIntegralController> _cdIntController;
  std::shared_ptr<IntegralCachingController> _integralCachingController;
};

} /* namespace Serenity */
#endif /* SYSTEMCONTROLLER_H */
