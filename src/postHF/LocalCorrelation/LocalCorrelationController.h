/**
 * @file LocalCorrelationController.h
 *
 * @date May 13, 2019
 * @author Moritz Bensberg
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
#ifndef POSTHF_LOCALCORRELATION_LOCALCORRELATIONCONTROLLER_H_
#define POSTHF_LOCALCORRELATION_LOCALCORRELATIONCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "data/OrbitalPair.h"                  //Orbital pair definition.
#include "data/SingleSubstitution.h"           //Singles definition.
#include "data/matrices/FockMatrix.h"          //Fock matrix definition.
#include "settings/CorrelatedMethodsOptions.h" //PNO_SETTINGS definition.
#include "settings/EmbeddingSettings.h"        //EmbeddingSettings.
/* Include Std and External Headers */
#include <limits> //numeric limits for default settings.
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {
/* Forward Declarations */
class SystemController;
class PAOController;
class PAOSelector;
class PNOConstructor;
class MO3CenterIntegralController;
class QuasiCanonicalPAODomainConstructor;
class TNOConstructor;
class OrbitalTriple;
class DomainOverlapMatrixController;
class SparseMapsController;
class OrbitalPairSet;
class OrbitalTripleSet;

// clang-format off
/**
 * @class LocalCorrelationSettings LocalCorrelationController.h
 * @brief The settings for the local correlation calculation.\n
 *
 * Settings used in ORCA:\n
 *   see JCTC 11, 1525-1539 (2015) table 1. And from Orca manual in brackets if deviating.\n
 *                T_CutPairs      T_CutPNO       T_CutMKN\n
 *   LoosePNO        10-3           10-6            10-3\n
 *   NormalPNO       10-4           3.33 10-7       10-3\n
 *   TightPNO        10-5           10-7            10-4(10-3)\n\n\n
 *
 *   Default is Normal-PNO
 *   Some of the settings refer directly to JCTC 11, 1525-1539 (2015) and J. Chem. Phys. 143, 034108 (2015).
 *
 * Settings:\n
 *   projectedEnvironment                    --- Levelshift occ. environment orbitals.\n
 *   useBPAlgorithm                          --- Use Boughton--Pulay algorithm for PAO selection.\n
 *   completenessThreshold                   --- Boughton--Pulay completeness threshold\n
 *   doiPairThreshold                        --- T_CutDOij\n
 *   doiPAOThreshold                         --- T_CutDO\n
 *   collinearDipoleScaling                  --- T_CutPre\n
 *   ccsdPairThreshold                       --- T_CutPair\n
 *   triplesSCMP2Scaling                     --- Scaling (for ccsdPairThreshold) for including weak pairs into the
 * triples. pnoThreshold                            --- T_CutPNO\n tnoThreshold                            ---
 * T_CutTNO\n singlesPNOFactor                        --- 0.03T_CutPNO see JCP 138, 034106 (2013) page 7\n
 *   orbitalToShellThreshold                 --- T_CutC\n
 *   mullikenThreshold                       --- T_CutMKN\n
 *   fockMatrixPrescreeningThreshold        --- F_Cut (LMP2 only)\n
 *   crudeDomainFactor                       --- Factor for approximate aux.-domain construction in prescreening.
 *   crudeStrongTripFactor                   --- Factor for the aux.-domain construction of strong triples.
 *   crudeWeakTripFactor                     --- Factor for the aux.-domain construction of weak triples.
 *   doiNetThreshold                         --- MnP truncation of DOIs.\n
 *   paoOrthogonalizationThreshold           --- Threshold for the canonical PAO orthogonalization.\n
 *   paoNormalizationThreshold               --- Renormalization threshold for PAOs.\n
 *   maximumMemoryRatio                      --- Maximum ratio of memory used for integral storage.\n
 *   dumpIntegrals                           --- Dump integrals to files and do not delete them.\n
 *   diisStartResidual                       --- start the DIIS after the residual is below this threshold.\n
 *   dampingFactor                           --- Initial damping for the amplitude optimization.
 *   dampingChange                           --- Damping change during amplitude optimization.
 *   finalDamping                            --- Final damping during amplitude optimization.
 *   diisMaxStore                            --- Maximum number of DIIS vector stored during amplitude optimization.
 *   setFaiZero                              --- Force the F_ai block to be zero.
 *   pnoCoreScaling                          --- Scaling factor for pairs/singles that include core-like orbitals.
 *   useFrozenCore                           --- Use the frozen-core approximation.
 *   energyCutOff                            --- Orbital energy cut-off to determine core-like orbitals.
 *   useTriplesCoreScaling                   --- Use the pnoCoreScaling for triples that contain core-like orbitals.
 *   pnoSettings                             --- PNO macro setting (LOOSE,NORMAL,TIGHT)
 *   method                                  --- Local-correlation method used.
 *   topDownReconstruction                   --- Enforce top-down ansatz for potential reconstruction.
 *   linearScalingSigmaVector                --- Build the sigma vector in DLPNO-CCSD directly from PNO-based integrals.
 *   extendedDomainScaling                   --- include additional pairs as close pairs in the sparse map / extended
 *                                               domain construction.
 *   enforceHFFockian                        ---  Enforce the use of the HF Fock operator.
 *   reuseFockMatrix                         --- If true, we will try to read the Fock matrix form disk.
 *   lowMemory                               --- Limit the number of 3-center integrals stored in memory and recalculate
 *                                               integrals more often.
 *   useProjectedOccupiedOrbitals            --- Add a projection operator to the Fock matrix to remove
 *                                               environment orbitals. This is only recommended if the environment
 *                                               orbitals are not orthogonal to the active system orbitals.
 *   ignoreMemoryConstraints                 --- If true, Serenity will assume that there is always enough memory for
 *                                               everything it tries to do.
 */
// clang-format on
struct LocalCorrelationSettings {
  LocalCorrelationSettings()
    : projectedEnvironment(false),
      useBPAlgorithm(false),
      completenessThreshold(0.02),
      doiPairThreshold(1e-5),
      doiPAOThreshold(std::numeric_limits<double>::infinity()),
      collinearDipoleScaling(0.01),
      ccsdPairThreshold(std::numeric_limits<double>::infinity()),
      triplesSCMP2Scaling(0.1),
      pnoThreshold(std::numeric_limits<double>::infinity()),
      tnoThreshold(1e-9),
      singlesPNOFactor(0.03),
      orbitalToShellThreshold(std::numeric_limits<double>::infinity()),
      mullikenThreshold(std::numeric_limits<double>::infinity()),
      crudeDomainFactor(10.0),
      crudeStrongTripFactor(10),
      crudeWeakTripFactor(100),
      fockMatrixPrescreeningThresholdd(1e-5),
      doiNetThreshold(1e-7),
      paoOrthogonalizationThreshold(1e-6),
      paoNormalizationThreshold(1e-6),
      maximumMemoryRatio(0.8),
      dumpIntegrals(false),
      diisStartResidual(1.0),
      dampingFactor(0.4),
      dampingChange(0.1),
      finalDamping(0.0),
      diisMaxStore(10),
      setFaiZero(true),
      pnoCoreScaling(0.01),
      useFrozenCore(false),
      useTriplesCoreScaling(false),
      pnoSettings(Options::PNO_SETTINGS::NORMAL),
      method(Options::PNO_METHOD::DLPNO_CCSD),
      topDownReconstruction(true),
      linearScalingSigmaVector(true),
      extendedDomainScaling(1),
      enforceHFFockian(false),
      reuseFockMatrix(true),
      lowMemory(false),
      useProjectedOccupiedOrbitals(false),
      ignoreMemoryConstraints(false) {
  }

 public:
  REFLECTABLE((bool)projectedEnvironment, // Levelshift occ. environment orbitals.
              /** ==== Bougthon--Pulay TYPE PAO Selection ==== **/
              (bool)useBPAlgorithm,          // Use Boughton--Pulay algorithm for PAO selection.
              (double)completenessThreshold, // Boughton--Pulay completeness threshold
              /** ====   DLPNO PARAMETERS   ==== **/
              (double)doiPairThreshold,        // T_CutDOij
              (double)doiPAOThreshold,         // T_CutDO
              (double)collinearDipoleScaling,  // T_CutPre
              (double)ccsdPairThreshold,       // T_CutPair
              (double)triplesSCMP2Scaling,     // T_CutMP2 factor
              (double)pnoThreshold,            // T_CutPNO
              (double)tnoThreshold,            // T_CutTNO
              (double)singlesPNOFactor,        // 0.03T_CutPNO
              (double)orbitalToShellThreshold, // T_CutC
              (double)mullikenThreshold,       // T_CutMKN
              (double)crudeDomainFactor,       // 100T_CutMKN see JCP 144, 024109 (2016) page 5
              (double)crudeStrongTripFactor, (double)crudeWeakTripFactor,
              (double)fockMatrixPrescreeningThresholdd, // F_Cut (LMP2 only)
              /** ====    GENERAL  PARAMETERS    ==== **/
              (double)doiNetThreshold,               // MnP truncation of DOIs
              (double)paoOrthogonalizationThreshold, // Threshold for the canonical PAO orthogonalization.
              (double)paoNormalizationThreshold,     // Renormalization threshold for PAOs
              (double)maximumMemoryRatio,            // Maximum ratio of memory used for integral storage.
              (bool)dumpIntegrals,                   // Dump integrals to files and do not delete them.
              (double)diisStartResidual,             // start the DIIS after the residual is below this threshold.
              (double)dampingFactor,                 // Initial damping for LMP2 or CCSD residual evaluation.
              (double)dampingChange,                 // Damping change.
              (double)finalDamping,                  // Final damping ratio.
              (unsigned int)diisMaxStore,            // Number of vectors stored for the DIIS.
              (bool)setFaiZero,                      // Force the F_ai block to be zero.
              (double)pnoCoreScaling,      // PNO threshold scaling for pairs/singles containing core-like orbitals.
              (bool)useFrozenCore,         // Use the frozen core approximation.
              (bool)useTriplesCoreScaling, // Scale the TNO-threshold for triples containing core-like orbitals with
                                           // pnoCoreScaling.
              (Options::PNO_SETTINGS)pnoSettings, // The PNO-macro setting.
              (Options::PNO_METHOD)method, // Flag for the local-correlation method (LMP2, DLPNO-CCSD, DLPNO-CCSD(T0) ...)
              (bool)topDownReconstruction,    // Enforce top-down ansatz for potential reconstruction.
              (bool)linearScalingSigmaVector, // Build the sigma vector in DLPNO-CCSD directly from PNO-based integrals.
              (double)extendedDomainScaling,  // include additional pairs as close pairs in the sparse map / extended
                                              // domain construction.
              (bool)enforceHFFockian,         // Enforce the use of the HF Fock operator.
              (bool)reuseFockMatrix,          // If true, we will try to read the Fock matrix form disk.
              (bool)lowMemory, // Limit the number of 3-center integrals stored in memory and recalculate integrals more
                               // often.
              (bool)useProjectedOccupiedOrbitals, // Use virtual orbitals that are occupied environment orbitals.
              (bool)ignoreMemoryConstraints // If true, Serenity assumes that it has enough memory to handle all orbital
                                            // pairs.
  )
  /** ==== FOCK MATRIX  CONSTRUCTION ==== **/
  EmbeddingSettings embeddingSettings;
  /**
   * @brief Parse the settings from the input to an instance of this class.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("LC")) {
      visit_each(*this, v);
      return true;
    }
    else if (this->embeddingSettings.visitAsBlockSettings(v, blockname)) {
      return true;
    }
    return false;
  }
  // Resolve PNO-macro settings.
  void resolvePNOSettings();

 private:
  // Check if a PNO setting was changed manually in the input.
  void checkSetting(double& setting, double defaultValue);
  // Check all PNO settings for manual changes.
  void resolvePNOSettingsSet(const Eigen::VectorXd values);
};

/**
 * @class LocalCorrelationController LocalCorrelationController.h
 * @brief A class that functions as a central hub for all local correlation calculation.
 *        This class constructs PAOs, does the prescreening, initializes orbital pairs
 *        and single substitutions etc.
 */
class LocalCorrelationController {
 public:
  /**
   * @brief Constructor.
   * @param system The associated system controller.
   * @param settings The local correlation settings.
   * @param environmentSystems Possible environment systems.
   * @param fockMatrix A possibly already calculated Fock matrix.
   * @param initialPairs The set of orbitals which should be correlated.
   *                     If non are given, all orbitals will be correlated.
   */
  LocalCorrelationController(std::shared_ptr<SystemController> system, LocalCorrelationSettings settings,
                             std::vector<std::shared_ptr<SystemController>> environmentSystems = {},
                             std::shared_ptr<FockMatrix<Options::SCF_MODES::RESTRICTED>> fockMatrix = nullptr,
                             std::vector<std::shared_ptr<OrbitalPair>> initialPairs = {},
                             Eigen::VectorXd orbitalWiseMullikenThresholds = Eigen::VectorXd::Zero(0),
                             Eigen::VectorXd orbitalToShellThresholds = Eigen::VectorXd::Zero(0),
                             Eigen::VectorXd orbitalWiseDOIPAOThresholds = Eigen::VectorXd::Zero(0));

  /**
   * @brief Getter for the active system controller.
   * @return The active system controller.
   */
  std::shared_ptr<SystemController> getActiveSystemController() {
    return _activeSystem;
  }
  /**
   * @brief Getter for the PAO controller.
   * @return The PAO controller.
   */
  std::shared_ptr<PAOController> getPAOController() {
    assert(_paoController);
    return _paoController;
  }
  /**
   * @brief Getter for the associated fock matrix. If environment systems
   *        are given, this will be an embedded fock matrix.
   * @return The fock matrix.
   */
  const FockMatrix<Options::SCF_MODES::RESTRICTED>& getFockMatrix() {
    assert(_fock);
    return *_fock;
  }
  /**
   * @brief Getter for the sparse maps controller.
   * @return The sparse maps controller.
   */
  std::shared_ptr<SparseMapsController> getSparseMapController(bool closeOnly = true);
  /**
   * @brief Getter for the triples sparse maps controller.
   * @return The sparse map controller.
   */
  std::shared_ptr<SparseMapsController> getTriplesSparseMapController();
  /**
   * @brief Getter for the linear scaling MO 3 center integral controller.
   * @param closeOnly If true, only close orbital pairs are concidered in the sparse map construction for prescreening.
   * @return The MO 3 center integral controller.
   */
  std::shared_ptr<MO3CenterIntegralController> getMO3CenterIntegralController(bool closeOnly = true);
  /**
   * @brief Getter for the linear scaling MO 3 center integral controller for the triples.
   * @return The MO 3 center integral controller.
   */
  std::shared_ptr<MO3CenterIntegralController> getTriplesMO3CenterIntegralController();
  /**
   * @brief Delete MO3 center integral controller to free some memory or initiate reconstruction upon
   *        call of the getter function.
   */
  void removeMO3CenterIntegralController();
  /**
   * @brief Delete approximate MO3 center integral controller to free some memory or initiate reconstruction upon
   *        call of the getter function.
   */
  void removeApproximateMO3CenterIntegralController();
  /**
   * @brief Getter for the approximate linear scaling MO 3 center integral controller.
   * @return The MO 3 center integral controller.
   */
  std::shared_ptr<MO3CenterIntegralController> getApproximateMO3CenterIntegralController();
  /**
   * @brief Getter for all orbital pairs.
   * @return All orbital pairs close+distant+very-distant.
   */
  std::vector<std::shared_ptr<OrbitalPair>> getAllOrbitalPairs() {
    return _allOrbitalPairs;
  }
  /**
   * @brief Getter for a subset of the orbital pairs.
   * @param type The orbital pair type.
   * @return The orbital pairs with the given type.
   */
  std::vector<std::shared_ptr<OrbitalPair>> getOrbitalPairs(OrbitalPairTypes type);

  /**
   * @brief Default destructor.
   */
  ~LocalCorrelationController();
  /**
   * @brief Construct k-sets for the close pairs.
   */
  void buildOrbitalPairCouplingMap();
  /**
   * @brief Construct the kl-pair lists.
   */
  void buildKLOrbitalPairs();
  /**
   * @brief Produces the PNO constructor for this correlation set up.
   * @param ssScaling Same spin scaling factor.
   * @param osScaling Opposite spin scaling factor.
   * @return The PNO constructor.
   */
  std::shared_ptr<PNOConstructor> producePNOConstructor(double ssScaling = 1.0, double osScaling = 1.0);
  /**
   * @brief Produces the quasi canonical PAO constructor.
   * @param ssScaling Same spin scaling factor.
   * @param osScaling Opposite spin scaling factor.
   * @param clear     Clear integrals and amplitudes after pair energy calculation.
   * @return The QCPAO constructor.
   */
  std::shared_ptr<QuasiCanonicalPAODomainConstructor> produceQCPAOConstructor(double ssScaling = 1.0,
                                                                              double osScaling = 1.0, bool clear = true);
  /**
   * @brief Produces the triple-natural orbital constructor.
   * @return The TNO constructor.
   */
  std::shared_ptr<TNOConstructor> produceTNOConstructor();
  /**
   * @brief Getter for the local correlation settings.
   * @return The local correlation settings.
   */
  const LocalCorrelationSettings& getSettings() {
    return _settings;
  }
  /**
   * @brief Construct singles for this correlation set up.
   * @param orbitalPairTypes The orbital pair types.
   */
  void initializeSingles(std::vector<OrbitalPairTypes> orbitalPairTypes = {OrbitalPairTypes::CLOSE});
  /**
   * @brief Getter for the singles.
   * @return The singles.
   */
  std::vector<std::shared_ptr<SingleSubstitution>> getSingles() {
    if (_singles.size() == 0)
      buildSingles();
    return _singles;
  }
  /**
   * @brief Divide the close orbital pair list into distant and close orbital pairs.
   *        This usually done for CCSD calculations.
   */
  void selectDistantOrbitalPairs();
  /**
   * @brief Getter for the sets of close orbital pairs. It can be expected that
   *        all integrals of one set can be kept in memory simultaneously.
   * @return The close orbital pair sets.
   */
  std::vector<std::shared_ptr<OrbitalPairSet>> getCloseOrbitalPairSets();
  /**
   * @brief Getter for the orbital triples sets (set-wise representation of the orbital triples).
   *        This will trigger the construction of the fitting domains for the orbital triples!
   * @return The orbital triple sets.
   */
  std::vector<std::shared_ptr<OrbitalTripleSet>> getOrbitalTripleSets();

  /**
   * @brief Build the singles for all close pairs.
   */
  void buildSingles();

  /**
   * @brief Getter for the orbital pair index map.
   * @return The orbital pair indices as a matrix occ x occ. Negative
   *         index denotes that the pair is not considered as "close".
   */
  const Eigen::MatrixXi& getOrbitalPairIndices();
  /**
   * @brief Getter for the orbital pair index maps.
   * @param orbitalPairType The orbital pair type.
   * @return The associated index map.
   */
  const Eigen::MatrixXi& getOrbitalPairIndices(OrbitalPairTypes orbitalPairType);

  /**
   * @brief Getter for a specific orbital pari ij.
   * @param i
   * @param j
   * @param types The types of orbital pairs which are considered.
   * @return The orbital pair or a nullptr if the orbital pair does not belong to one of the given types.
   */
  std::shared_ptr<OrbitalPair> getOrbitalPair(unsigned int i, unsigned int j,
                                              std::vector<OrbitalPairTypes> types = {OrbitalPairTypes::CLOSE,
                                                                                     OrbitalPairTypes::DISTANT,
                                                                                     OrbitalPairTypes::VERY_DISTANT});

  /**
   * @brief Getter for the triples. Note that this list will depend on the prescreening done
   *        by calling functions like selectDistantOrbitalPairs() and will be updated accordingly.
   *        Thus a copy of the pointer list is parsed and not only a reference! Note that the triples may
   *        not be fully initialized yet!
   * @return The triples.
   */
  std::vector<std::shared_ptr<OrbitalTriple>> getOrbitalTriples() {
    if (_orbitalTriples.size() == 0 && _allOrbitalPairs.size() > 1)
      _orbitalTriples = constructOrbitalTriplets();
    return _orbitalTriples;
  }
  /**
   * @brief Getter for the domain overlap matrix controller.
   * @return The domain overlap matrix controller. May be constructed upon call.
   *         This will lead to the calculation of all overlap matrix integrals!
   */
  std::shared_ptr<DomainOverlapMatrixController> getDomainOverlapMatrixController();
  /**
   * @brief Set the orbital triples list to the given list.
   * @param triples The triples.
   */
  void setOrbitalTriples(std::vector<std::shared_ptr<OrbitalTriple>> triples);

  const Eigen::MatrixXd& getPairEnergyMatrix();
  void writePairEnergies(std::string postfix = "");

  std::string getPairIntegralFileName();

 private:
  // The active system controller.
  std::shared_ptr<SystemController> _activeSystem;
  // The settings.
  LocalCorrelationSettings _settings;
  // The environment systems.
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  // The fock matrix.
  std::shared_ptr<FockMatrix<Options::SCF_MODES::RESTRICTED>> _fock;
  // The PAO controller.
  std::shared_ptr<PAOController> _paoController;
  // For every occupied orbital i (a/b) the indices of the pairs in which it is a part of.
  // Negative if not paired.
  std::shared_ptr<Eigen::MatrixXi> _orbitalPairIndices;
  std::shared_ptr<Eigen::MatrixXi> _distantTriplesOrbitalPairIndices;
  std::shared_ptr<Eigen::MatrixXi> _distantOrbitalPairIndices;
  std::shared_ptr<Eigen::MatrixXi> _veryDistantOrbitalPairIndices;
  std::shared_ptr<Eigen::MatrixXi> _sparseMapClosePairIndices;
  // The sparse maps controller.
  std::shared_ptr<SparseMapsController> _sparseMapsController;

  std::shared_ptr<SparseMapsController> _triplesSparseMapsController;
  // The MO 3 center integral controller.
  std::shared_ptr<MO3CenterIntegralController> _accurateMo3CenterIntegralController;
  // An approximated version of the MO 3 center integral controller with smaller fitting domains.
  std::shared_ptr<MO3CenterIntegralController> _approximateMo3CenterIntegralController;
  // The MO3 center integral controller for the triplet integrals.
  std::shared_ptr<MO3CenterIntegralController> _tripletMo3CenterIntegralController;
  // List of all orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _allOrbitalPairs;
  // List of the close orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _closeOrbitalPairs;
  std::vector<std::shared_ptr<OrbitalPair>> _sparseMapConstructionPairs;
  // List of the distant orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _distantOrbitalPairs = {};
  std::vector<std::shared_ptr<OrbitalPair>> _distantOrbitalPairsTriples = {};
  // List of the very-distant orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _veryDistantOrbitalPairs;
  // List of the singles.
  std::vector<std::shared_ptr<SingleSubstitution>> _singles = {};
  // The indices of the singles.
  std::shared_ptr<Eigen::VectorXi> _singlesIndices;
  // The orbital triples.
  std::vector<std::shared_ptr<OrbitalTriple>> _orbitalTriples;
  // The DomainOverlapMatrixController given access to all overlap integrals.
  std::shared_ptr<DomainOverlapMatrixController> _domainOverlapMatrixController;
  // Construct the fock matrix.
  void constructFockMatrix();
  // Produce the PAO selecter specified in the settings.
  std::shared_ptr<PAOSelector> producePAOSelector();
  // Prints information about the number of PAOs per orbital.
  void printPAOInfo(const std::shared_ptr<Eigen::SparseMatrix<int>> occupiedToPAOOrbitalMap);
  // Build a list of all orbital pairs within the system.
  std::vector<std::shared_ptr<OrbitalPair>> buildInitialOrbitalPairs();
  // Build the PAO pair domains for all close pairs.
  inline void buildPAOPairDomains(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs,
                                  const SparseMap& occupiedToPAOOrbitalMap);
  // Constructes the orbital triples list.
  std::vector<std::shared_ptr<OrbitalTriple>> constructOrbitalTriplets();
  // Constructs the orbita pair index map for the given orbital pair type.
  Eigen::MatrixXi buildOrbitalPairIndices(OrbitalPairTypes orbitalPairType);
  // The sets of close orbital pairs.
  std::vector<std::shared_ptr<OrbitalPairSet>> _closeOrbitalPairSets;
  // The sets of orbital triples.
  std::vector<std::shared_ptr<OrbitalTripleSet>> _orbitalTripleSets;
  // Map between occupied orbitals and PAOs.
  std::shared_ptr<SparseMap> _occupiedToPAOOrbitalMap;
  // Print PNO-settings info to the output (depends on print level).
  void printSettings();
  // Domain construction thresholds:
  // Fitting domains.
  Eigen::VectorXd _orbitalWiseMullikenThresholds;
  // Coefficient truncation.
  Eigen::VectorXd _orbitalToShellThresholds;
  // PAO truncation.
  Eigen::VectorXd _orbitalWiseDOIPAOThresholds;
  // The matrix (occ x occ) of the pair energies.
  std::shared_ptr<Eigen::MatrixXd> _pairEnergyMatrix;
  // Name of the pair integrals file for CCSD.
  std::string _pairIntegralFileName = "PairIntegrals.h5";
  /**
   * @brief Sort the orbital pairs/triples according to their fitting domain. Grouping orbital pairs/triples into
   *        groups that have a similar fitting domain.
   * @tparam OrbitalTuple OrbitalPair or OrbitalTriple.
   * @param orbitalPairs The orbital tuples to be sorted.
   * @return The sorted list of orbital tuples.
   */
  template<class OrbitalTuple>
  std::vector<std::shared_ptr<OrbitalTuple>> sortOrbitalPairs(std::vector<std::shared_ptr<OrbitalTuple>> orbitalPairs);
  /**
   * @brief Set the extended fitting domains for the given list of orbital pairs.
   * @param orbitalPairs The orbital pairs.
   * @param occToK The fitting domain map (occ to aux. function).
   */
  void setExtendedAuxDomain(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs, const SparseMap& occToK);
  /**
   * @brief Extract (get and remove from list) the orbital pair from the given list with the maximum overlap in the
   *        fitting domain as the given reference fitting domain.
   * @tparam OrbitalTuple OrbitalPair or OrbitalTriple.
   * @param orbitalPairs The list of orbital triples to extract from.
   * @param referenceDomain The reference domain.
   * @return The orbital tuple with the maximum overlap.
   */
  template<class OrbitalTuple>
  std::shared_ptr<OrbitalTuple> extractOrbitalPairWithMaxOverlapDomain(std::vector<std::shared_ptr<OrbitalTuple>>& orbitalPairs,
                                                                       const Eigen::VectorXi& referenceDomain);
  /**
   * @brief Extract (get and remove from list) all orbital tuples from the given list for which the fitting domain is
   *        fully contained in the given reference fitting domain.
   * @tparam OrbitalTuple OrbitalPair or OrbitalTriple.
   * @param orbitalPairs The orbital truple list to extract from.
   * @param referenceDomain The reference fitting domain.
   * @return The list of orbital tuples fully contained.
   */
  template<class OrbitalTuple>
  std::vector<std::shared_ptr<OrbitalTuple>>
  extractOrbitalPairsFullyContained(std::vector<std::shared_ptr<OrbitalTuple>>& orbitalPairs,
                                    const Eigen::VectorXi& referenceDomain);
  /**
   * @brief Get the overlap between two fitting domains.
   * @param domainI The first domain.
   * @param domainJ the second domain.
   * @return The number aux. shells in the intersection of both domains.
   */
  int getDomainOverlap(const Eigen::VectorXi& domainI, const Eigen::SparseVector<int>& domainJ);
};

} /* namespace Serenity */

#endif /* POSTHF_LOCALCORRELATION_LOCALCORRELATIONCONTROLLER_H_ */
