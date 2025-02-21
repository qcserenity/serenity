/**
 * @file SparseMapsController.h
 *
 * @date May 9, 2019
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

#ifndef DATA_SPARSEMAPSCONTROLLER_H_
#define DATA_SPARSEMAPSCONTROLLER_H_
/* Include Std and External Headers */
#include <Eigen/SparseCore> //Sparse maps
#include <memory>           //smrt_ptr
#include <vector>           //std::vector

namespace Serenity {
/* Type definitions */
typedef Eigen::SparseMatrix<int> SparseMap;

/* Forward Declarations */
class PAOController;
class OrbitalPair;
class PAOSelector;
class SystemController;
class OrbitalTriple;

/**
 * @class SparseMapsController SparseMapsController.h
 * @brief A class that generates and holds all sparse maps for a given
 *        system and set of orbitals.\n
 *        For the general concept of sparse maps used in prescreening and
 *        the definitions of the maps controlled by this class see:\n
 *          J.Chem.Phys. 143, 034108 (2015).
 */
class SparseMapsController {
 public:
  /**
   * @brief Constructor.
   * @param system The associated system.
   * @param paoController The PAO controller.
   * @param occupiedToPAOOrbitalMap The PAO selection for the occupied orbitals.
   * @param closeOrbitalPairs The set of close orbital pairs.
   * @param distantOrbitalPairs The set of distant orbital pairs.
   * @param orbitalWiseMullikenThresholdss The Mulliken population thresholds which are used to truncate
   *                                        the auxiliary basis.
   * @param orbitalToShellThresholds The Mulliken population threshold which are used to truncate
   *                                 the AO basis for each occupied orbital.
   * @param strongTripletMullikenScaling Mulliken population scaling for strong triplets.
   * @param weakTripletMullikenScaling   Mulliken population scaling for weak triplets.
   * @param klListExtension Extend the extended pair-lists using the kl-sets.
   * @param triples The triples list for the extended map construction.
   */
  SparseMapsController(std::shared_ptr<SystemController> system, std::shared_ptr<PAOController> paoController,
                       std::shared_ptr<SparseMap> occupiedToPAOOrbitalMap,
                       std::vector<std::shared_ptr<OrbitalPair>> closeOrbitalPairs,
                       std::vector<std::shared_ptr<OrbitalPair>> distantOrbitalPairs,
                       Eigen::VectorXd orbitalWiseMullikenThresholds, Eigen::VectorXd orbitalToShellThresholds,
                       double strongTripletMullikenScaling = 10, double weakTripletMullikenScaling = 100,
                       bool klListExtension = false, std::vector<std::shared_ptr<OrbitalTriple>> triples = {});
  /**
   * @brief Constructor.
   * @param system The system controller.
   * @param occupiedCoefficients The orbital coefficients of the occupied orbitals.
   * @param virtualCoefficients The virtual orbital coefficients.
   * @param occupiedToPAOOrbitalMap The PAO selection for the occupied orbitals.
   * @param closeOrbitalPairs The set of close orbital pairs.
   * @param distantOrbitalPairs The set of distant orbital pairs.
   * @param orbitalWiseMullikenThresholdss The Mulliken population thresholds which are used to truncate
   *                                        the auxiliary basis.
   * @param orbitalToShellThresholds The Mulliken population threshold which are used to truncate
   *                                 the AO basis for each occupied orbital.
   * @param strongTripletMullikenScaling Mulliken population scaling for strong triplets.
   * @param weakTripletMullikenScaling   Mulliken population scaling for weak triplets.
   * @param klListExtension Extend the extended pair-lists using the kl-sets.
   * @param triples The triples list for the extended map construction.
   */
  SparseMapsController(std::shared_ptr<SystemController> system, std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients,
                       std::shared_ptr<Eigen::MatrixXd> virtualCoefficients, std::shared_ptr<SparseMap> occupiedToPAOOrbitalMap,
                       std::vector<std::shared_ptr<OrbitalPair>> closeOrbitalPairs,
                       std::vector<std::shared_ptr<OrbitalPair>> distantOrbitalPairs,
                       Eigen::VectorXd orbitalWiseMullikenThresholds, Eigen::VectorXd orbitalToShellThresholds,
                       double strongTripletMullikenScaling = 10, double weakTripletMullikenScaling = 100,
                       bool klListExtension = false, std::vector<std::shared_ptr<OrbitalTriple>> triples = {});
  /**
   * @brief Default destructor.
   */
  ~SparseMapsController() = default;
  /**
   * @brief Getter for the map occ->atom.
   * @return The map occ->atom.
   */
  const SparseMap& getOccToAtomMap();
  /**
   * @brief Getter for the map atom->PAO.
   * @return The map atom->PAO.
   */
  const SparseMap& getAtomToPAOMap();
  /**
   * @brief Getter for the map occ->PAO.
   * @return The map occ->PAO.
   */
  const SparseMap& getOccToPAOMap();
  /**
   * @brief Getter for the map shell->occ.
   * @return The map shell->occ.
   */
  const SparseMap& getShellToOccMap();
  /**
   * @brief Getter for the map occ->aux. shell.
   * @return The map occ->aux. shell.
   */
  const SparseMap& getOccToAuxShellMap();
  /**
   * @brief Getter for the map occ->aux. shell constructed with
   *        the triples fitting domain thresholds.
   * @param weak Defines the threshold used for the map construction.
   *        Weak and strong thresholds differ.
   * @return The map occ->aux. shell.
   */
  const SparseMap& getTripletOccToAuxShellMap(bool weak);
  /**
   * @brief Getter for the extended map occ->aux. shell.
   *        Map extension is done with close and distant pairs.
   * @return The extended map occ->aux. shell.
   */
  const SparseMap& getExtendedOccToAuxShellMap();
  /**
   * @brief Getter for the extended map occ->PAO.
   * @return The extended map occ->PAO.
   */
  const SparseMap& getExtendedOccToPAOMap();
  /**
   * @brief Getter for the extended map aux. shell->AO shell
   *        for the AO functions used for the occupied orbitals.
   * @return The extended map aux. shell->AO shell.
   */
  const SparseMap& getExtendedKtoRhoMap();
  /**
   * @brief Getter for the extended map aux. shell->AO shell
   *        for the AO functions used for the virtual orbitals.
   * @return The extended map aux. shell->AO shell.
   */
  const SparseMap& getExtendedKtoSigmaMap();

  /**
   * @brief Getter for the extended map aux. shell->PAO
   * @return The extended map aux. shell->PAO.
   */
  const SparseMap& getExtendedKtoPAOMap();

  /**
   * @brief Getter for the map AO shell->PAO.
   * @return The map AO shell->PAO.
   */
  const SparseMap& getShellToPAOMap();
  /**
   * @brief Getter for the map atom->aux shell.
   * @return The map atom->aux shell.
   */
  const SparseMap& getAtomToAuxShellMap();
  /**
   * @brief Getter for the (triples) extended map between aux. shells and PAOs.
   * @return The map.
   */
  const SparseMap& getExtendedKtoPAOMap_triples();
  /**
   * @brief Getter for the (triples) extended map between aux. shells and virt. shells.
   * @return The map.
   */
  const SparseMap& getExtendedKtoSigmaMap_triples();
  /**
   * @brief Getter for the (triples) extended map between aux. shells and occupied shells.
   * @return The map.
   */
  const SparseMap& getExtendedKtoRhoMap_triples();
  /**
   * @brief Getter for the (triples) extended map between occupied orbitals and PAOs.
   * @return The map.
   */
  const SparseMap& getExtendedOccToPAOMap_triples();
  /**
   * @brief Getter for the (triples) extended map between occupied orbitals and aux. shells.
   * @return The map.
   */
  const SparseMap& getExtendedOccToAuxShellMap_triples();

 private:
  // The system controller.
  std::shared_ptr<SystemController> _system;
  // The orbital coefficients of the occupied orbitals.
  std::shared_ptr<Eigen::MatrixXd> _occupiedCoefficients;
  // The virtual orbital coefficients, e.g., PAO coefficients.
  std::shared_ptr<Eigen::MatrixXd> _virtualCoefficients;
  /*
   * Example: xToYMap --> x-->columns, y-->rows
   */
  // Maps occ. orbitals --> atoms (localization)
  std::shared_ptr<SparseMap> _occupiedOrbitalToAtomMap;
  // Maps occ. orbitals --> PAOs
  std::shared_ptr<SparseMap> _occupiedToPAOOrbitalMap;
  // Maps atoms --> PAOs
  std::shared_ptr<SparseMap> _atomToPAOMap;
  // Maps AO-shells --> occ. orbitals
  std::shared_ptr<SparseMap> _shellToOccMap;
  // Maps AO-shells --> PAOs
  std::shared_ptr<SparseMap> _shellToPAOMap;
  // Maps atoms --> aux. shells
  std::shared_ptr<SparseMap> _atomToAuxShellMap;
  // Map occ --> aux. shell. (pairs.)
  std::shared_ptr<SparseMap> _occToKMap;
  // Map occ --> aux. shell. (strong triples)
  std::shared_ptr<SparseMap> _strongTripletOccToKMap;
  // Map occ --> aux. shell. (weak triples)
  std::shared_ptr<SparseMap> _weakTripletOccToKMap;
  // Extended map occ --> aux. shell (all pairs)
  std::shared_ptr<SparseMap> _extendedOccToK;
  // Extended map occ --> aux. shell (close pairs)
  std::shared_ptr<SparseMap> _closeExtendedOccToK;
  // Extended map occ --> PAOs (all pairs)
  std::shared_ptr<SparseMap> _extendedOccToPAO;
  // Extended map aux. shell --> basis function shell (all pairs, occupied orbitals)
  std::shared_ptr<SparseMap> _extendedAuxShellToRho;
  // Extended map aux. shell --> basis function shell (all pairs, PAOs)
  std::shared_ptr<SparseMap> _extendedAuxShellToSigma;
  // Extended map aux. shell --> PAOs (all pairs)
  std::shared_ptr<SparseMap> _extendedAuxShellToPAO;
  // Triples extended map occ --> aux. shell
  std::shared_ptr<SparseMap> _extendedOccToK_triples;
  // Triples extended map occ --> PAOs (all pairs)
  std::shared_ptr<SparseMap> _extendedOccToPAO_triples;
  // Triples extended map aux. shell --> basis function shell (all pairs, occupied orbitals)
  std::shared_ptr<SparseMap> _extendedAuxShellToRho_triples;
  // Triples extended map aux. shell --> basis function shell (all pairs, PAOs)
  std::shared_ptr<SparseMap> _extendedAuxShellToSigma_triples;
  // Triples extended map aux. shell --> PAOs (all pairs)
  std::shared_ptr<SparseMap> _extendedAuxShellToPAO_triples;

  // Construct atom-wise representation of orbital-wise thresholds from atom to orbital map.
  Eigen::VectorXd convertToAtomWiseThresholds(Eigen::VectorXd orbitalWiseThresholds);
  // Constructs a shell-->orbital map, based on a Mulliken net-population analysis.
  std::shared_ptr<SparseMap> constructShellToOrbitalMap(const Eigen::MatrixXd& coefficients, const Eigen::VectorXd& thresholds);
  /*
   * Constructs an extended map of type occ->x by defining the new occ_i->x as
   * the union of all occ_j->x for which the pair ij exists.
   *
   * Use all orbital pairs for the map extension.
   */
  std::shared_ptr<SparseMap> buildExtendedMap(const Eigen::SparseMatrix<int>& initialMap);
  // Use only close orbital pairs for the map extension.
  std::shared_ptr<SparseMap> buildCloseExtendedMap(const Eigen::SparseMatrix<int>& initialMap);
  // Use the triples in order to extend the map.
  std::shared_ptr<SparseMap> buildTriplesExtendedMap(const Eigen::SparseMatrix<int>& initialMap);
  std::shared_ptr<SparseMap> buildTriplesExtendedMap(const Eigen::SparseMatrix<int>& strongMap,
                                                     const Eigen::SparseMatrix<int>& weakMap);
  // Construct a prescreening map for occupied orbital --> atom assignment.
  SparseMap buildOccToAtomMap(Eigen::VectorXd thresholds);
  // The various orbital pair lists.
  // All orbital pairs: _orbitalPairs = _closeOrbitalPairs + _distantOrbitalPairs
  std::vector<std::shared_ptr<OrbitalPair>> _orbitalPairs;
  // Close orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _closeOrbitalPairs;
  // Distant orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _distantOrbitalPairs;
  // Prescreening thresholds for pair aux. domain construction (orbital-wise)
  const Eigen::VectorXd _orbitalWiseMullikenThresholds;
  // Prescreening threshold for orbital coefficients.
  const Eigen::VectorXd _orbitalToShellThresholds;
  // Prescreening threshold for strong triplet aux. domain construction.
  double _strongTripletMullikenScaling;
  // Prescreening threshold for weak triplet aux. domain construction.
  double _weakTripletMullikenScaling;
  // Getter for the orbital wise Mulliken charges.
  const Eigen::MatrixXd& getOrbitalWiseMullikenPopulations();
  // The orbital-wise mulliken populations.
  std::shared_ptr<Eigen::MatrixXd> _orbitalWiseMullikenPopulations;
  // Flag for the construction of pair-extended maps by the kl-pair lists.
  bool _klListExtension;
  // Triples
  std::vector<std::shared_ptr<OrbitalTriple>> _triples;
};

} /* namespace Serenity */

#endif /* DATA_SPARSEMAPSCONTROLLER_H_ */
