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

typedef std::vector<std::shared_ptr<OrbitalPair>> OrbitalPairSet;

/* Forward Declarations */
class SystemController;
class PAOController;
class PAOSelecter;
class PNOConstructor;
class MO3CenterIntegralController;
class QuasiCanonicalPAODomainConstructor;
class TNOConstructor;
class OrbitalTriple;
class DomainOverlapMatrixController;
class SparseMapsController;

/*
 * Flags for the PNO method that is used.
 */
enum class PNO_METHOD { DLPNO_MP2, DLPNO_CCSD, DLPNO_CCSD_T0 };

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
 *
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
 *   fockMatrixPreescreeningThreshold        --- F_Cut (LMP2 only)\n
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
 *   useTriplesCoreScaling                   --- Use the pnoCoreScaling for triples that contain core-like orbtials.
 *   pnoSettings                             --- PNO macro setting (LOOSE,NORMAL,TIGHT)
 *   method                                  --- Local-correlation method used.
 * Some of the settings refer directly to JCTC 11, 1525-1539 (2015) and J. Chem. Phys. 143, 034108 (2015).
 */
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
      fockMatrixPreescreeningThreshold(1e-5),
      doiNetThreshold(1e-7),
      paoOrthogonalizationThreshold(1e-6),
      paoNormalizationThreshold(1e-6),
      maximumMemoryRatio(0.5),
      dumpIntegrals(false),
      diisStartResidual(5e-2),
      dampingFactor(0.4),
      dampingChange(0.1),
      finalDamping(0.0),
      diisMaxStore(10),
      setFaiZero(true),
      pnoCoreScaling(0.01),
      useFrozenCore(false),
      energyCutOff(-5.0),
      useTriplesCoreScaling(false),
      pnoSettings(Options::PNO_SETTINGS::NORMAL),
      method(PNO_METHOD::DLPNO_CCSD),
      topDownReconstruction(true) {
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
              (double)fockMatrixPreescreeningThreshold, // F_Cut (LMP2 only)
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
              (double)energyCutOff,        // Energy cut of for assigning core-like orbitals.
              (bool)useTriplesCoreScaling, // Scale the TNO-threshold for triples containing core-like orbitals with
                                           // pnoCoreScaling.
              (Options::PNO_SETTINGS)pnoSettings, // The PNO-macro setting.
              (PNO_METHOD)method, // Flag for the local-correlation method (LMP2, DLPNO-CCSD, DLPNO-CCSD(T0) ...)
              (bool)topDownReconstruction // Enforce top-down ansatz for potential reconstruction.
  )
  /** ==== FOCK MATRIX  CONSTRUCTION ==== **/
  EmbeddingSettings embeddingSettings;
  /**
   * @brief Parse the settings from the input an instance of this class.
   * @param c The settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  bool visitSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("LC")) {
      visit_each(*this, v);
      return true;
    }
    else if (!blockname.compare("EMB")) {
      visit_each(this->embeddingSettings, v);
      return true;
    }
    return false;
  }
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
                             std::vector<std::shared_ptr<OrbitalPair>> initialPairs = {});

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
  std::shared_ptr<SparseMapsController> getSparseMapController(bool closeOnly = false);
  /**
   * @brief Getter for the linear scaling MO 3 center integral controller.
   * @return The MO 3 center integral controller.
   */
  std::shared_ptr<MO3CenterIntegralController> getMO3CenterIntegralController(bool closeOnly = false);
  /**
   * @brief Delete MO3 center integral controller to free some memory or initiate reconstruction upon
   *        call of the getter function.
   */
  void removeMO3CenterIntegralControlle() {
    _accurateMo3CenterIntegralController = nullptr;
  }
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
  std::vector<std::shared_ptr<OrbitalPair>> getOrbitalPairs(OrbitalPairTypes type) {
    switch (type) {
      case OrbitalPairTypes::CLOSE:
        return _closeOrbitalPairs;
      case OrbitalPairTypes::DISTANT_TRIPLES:
        return _distantOrbitalPairsTriples;
      case OrbitalPairTypes::DISTANT:
        return _distantOrbitalPairs;
      case OrbitalPairTypes::VERY_DISTANT:
        return _veryDistantOrbitalPairs;
    }
    assert(false && "Case not handled in switch statement");
    return _allOrbitalPairs;
  }

  /**
   * @brief Default destructor.
   */
  ~LocalCorrelationController() = default;
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
   * @return The QCPAO constructor.
   */
  std::shared_ptr<QuasiCanonicalPAODomainConstructor> produceQCPAOConstructor(double ssScaling = 1.0, double osScaling = 1.0);
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
   */
  void initializeSingles();
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
  std::vector<OrbitalPairSet> getCloseOrbitalPairSets();

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
   *        Thus a copy of the pointer list is parsed and not only a reference!
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
  // The sparse maps controller.
  std::shared_ptr<SparseMapsController> _sparseMapsController;
  // The MO 3 center integral controller.
  std::shared_ptr<MO3CenterIntegralController> _accurateMo3CenterIntegralController;
  // An approximated version of the MO 3 center integral controller with smaller fitting domains.
  std::shared_ptr<MO3CenterIntegralController> _approximateMo3CenterIntegralController;
  // List of all orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _allOrbitalPairs;
  // List of the close orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _closeOrbitalPairs;
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
  std::shared_ptr<PAOSelecter> producePAOSelecter();
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
  std::vector<OrbitalPairSet> _closeOrbitalPairSets;
  // Map between occupied orbitals and PAOs.
  std::shared_ptr<SparseMap> _occupiedToPAOOrbitalMap;
  // Resolve PNO-macro settings.
  void resolvePNOSettings();
  // Check if a PNO setting was changed manually in the input.
  void checkSetting(double& setting, double defaultValue);
  // Check all PNO settings for manual changes.
  void resolvePNOSettingsSet(const Eigen::VectorXd values);
  // Print PNO-settings info to the output (depends on print level).
  void printSettings();
};

} /* namespace Serenity */

#endif /* POSTHF_LOCALCORRELATION_LOCALCORRELATIONCONTROLLER_H_ */
