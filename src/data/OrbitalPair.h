/**
 * @file OrbitalPair.h
 *
 * @date May 3, 2019
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

#ifndef DATA_ORBITALPAIR_H_
#define DATA_ORBITALPAIR_H_
/* Include Serenity Internal Headers */
#include "math/Matrix.h" //Four index objects.
/* Include Std and External Headers */
#include <Eigen/Dense>      // Dense matrices.
#include <Eigen/SparseCore> // Sparse matrices.
#include <limits>
#include <memory> // smrt_ptr
#include <string> // Integral file names.
#include <vector> // std::vector

namespace H5 {
class H5File;
} // namespace H5

namespace Serenity {

/* Forward Declarations */
namespace HDF5 {
using H5File = H5::H5File;
} // namespace HDF5
class SingleSubstitution;
class CouplingOrbitalSet;
class KLOrbitalSet;
class DomainOverlapMatrixController;

/**
 * @brief Enum class for different orbital pair types.
 *   CLOSE:           Strong pairs which are iterated in the CCSD equations.
 *   DISTANT_TRIPLES: Pairs for which the semi-canonical MP2 pair energy is below
 *                    the scaled CCSD-pair truncation threshold. These pairs may
 *                    be used for the determination of orbital pairs for the triples
 *                    correction.
 *   DISTANT:         Pairs for which the semi-canonical MP2 pari energy is below
 *                    the CCSD-pair truncation threshold. These pairs may be treated
 *                    using SC-MP2 or Local MP2
 *   VERY_DISTANT:    Pairs for which only a dipole approximation is calculated.
 *   SPARSE_MAP_CLOSE:The close orbital pairs used in the sparse-map construction.
 */
enum class OrbitalPairTypes { CLOSE, DISTANT_TRIPLES, DISTANT, VERY_DISTANT, SPARSE_MAP_CLOSE };
/**
 * @brief Enum class with all integral sets used in DLPNO-CCSD.
 */
enum class FOUR_CENTER_INTEGRAL_TYPE {
  ia_jb,               // Pair
  ac_bd,               // Pair
  ka_bc,               // kSet
  ij_ab,               // Pair
  ja_ik,               // kSet
  ia_jk,               // kSet
  ij_ak,               // kSet
  ik_jl,               // Pair
  ia_kc,               // kSet
  ik_ca,               // kSet
  ja_kc,               // kSet
  jk_ca,               // kSet
  ki_la,               // Pair
  kj_la,               // Pair
  ia_bc,               // Pair
  ja_bc,               // Pair
  jc_ab,               // Pair
  ic_ab,               // Pair
  iaS_jbSX2_M_ij_aSbS, // Pair
  s_ij_kj,             // kSet
  s_ij_ik,             // kSet
  s_ij_kl,             // Pair
  s_ij_i,              // Pair
  s_ij_j,              // Pair
  s_ij_k,              // kSet
  s_i_j,               // Pair
  ab_kcX2_M_ak_bc
}; // kSet

/**
 * @class OrbitalPair OrbitalPair.h
 * @brief A storage class that is concerned with the double amplitudes etc. in local-coupled cluster.
 *        Holds integrals, amplitudes, pair energies etc. and manages the storage of these on disk
 *        if desired.
 */
class OrbitalPair {
 public:
  /**
   * @brief Constructor.
   * @param orbital_i Index of the occupied orbital i.
   * @param orbital_j Index of the occupied orbital j.
   * @param pnoThreshold The PNO threshold used for this pair.
   */
  OrbitalPair(unsigned int orbital_i, unsigned int orbital_j, double pnoThreshold, double mainPairEnergyThreshold,
              double collinearScaling, double sparseMapsConstructionPairScaling = std::numeric_limits<double>::infinity());
  /**
   * @brief Destructor.
   */
  ~OrbitalPair();
  ///@brief The indices of the pair.
  const unsigned int i;
  ///@brief The indices of the pair.
  const unsigned int j;
  ///@brief The amplitudes in PAO pair basis.
  Eigen::MatrixXd t_ij;
  ///@brief The DLPNO-CCSD residual.
  Eigen::MatrixXd residual;
  ///@brief If false, this pair will never be used for orbital triples construction.
  ///       This is useful for embedding.
  bool eligibleForTriples = true;

  /* ========INTEGRAL SETS ======== */
  // Matrices in "pure" ij-Pair PNO basis.
  ///@brief The exchange integrals in the PAO pair basis.
  Eigen::MatrixXd k_ij; // (ia|jb) --> ia_jb
  ///@brief (ac|bd) integrals, a,b,c,d in [ij].
  std::unique_ptr<Matrix<Eigen::MatrixXd>> ac_bd; //(ac|bd) ~ 30x30x30x30 ~ 3.5 MB per pair, only the upper triangle is
                                                  // saved!
  ///@brief (ij|ab) integrals, a,b in [ij].
  Eigen::MatrixXd ij_ab; //(ij|ab)
  ///@brief (ik|jl) integrals.
  Eigen::MatrixXd ik_jl; //(ik|jl)
  // Integrals in mixed singles/doubles PNO basis
  ///@brief (ki|la) integrals, a in [j]
  std::vector<Eigen::MatrixXd> ki_la; // (ki|la) matrices-> a * nK x nK, a in [j]
  ///@brief (kj|la) integrals, a in [i]
  std::vector<Eigen::MatrixXd> kj_la; // (kj|la) matrices-> a * nK x nK, a in [i]
  ///@brief (ia|bc) integrals, a,b in [ij] and c in [j]
  std::vector<Eigen::MatrixXd> ia_bc; // a,b in [ij] and c in [j]
  ///@brief (ia|bc) integrals, a,b in [ij] and c in [i]
  std::vector<Eigen::MatrixXd> ja_bc; // a,b in [ij] and c in [i]
  ///@brief (jc|ab) integrals, a,b in [ij] and c in [i]
  std::vector<Eigen::MatrixXd> jc_ab; // a,b in [ij] and c in [i]
  ///@brief (ic|ab) integrals, a,b in [ij] and c in [j]
  std::vector<Eigen::MatrixXd> ic_ab; // a,b in [ij] and c in [j]
  ///@brief 2(ia|jb)-(ij|ab), b in [j] and a in [i]
  Eigen::MatrixXd iaS_jbSX2_M_ij_aSbS; // b in [j] and a in [i]

  /* ================================ */
  ///@brief The transformation to the linear independent quasi canonical PAO/PNO basis.
  Eigen::MatrixXd toPAODomain;
  ///@brief The indices of the atoms selected for the orbital pair.
  Eigen::SparseVector<int> paoDomain;
  ///@brief Sparse projection to the PAO domain used by this pair.
  Eigen::SparseMatrix<double> domainProjection;
  ///@brief The uncoupled linear (in terms of LMP2) term in T_ij.
  Eigen::MatrixXd uncoupledTerm;
  ///@brief The dipole approximation to the pair energy.
  double dipolePairEnergy = 0.0;
  ///@brief The upper bound to the dipole approximation assuming collinear dipoles.
  double dipoleCollinearPairEnergy = 0.0;
  ///@brief The (fully converged) local MP2 pair energy.
  double lMP2PairEnergy = 0.0;
  ///@brief Semi-canonical MP2 pair energy
  double scMP2PairEnergy = 0.0;
  ///@brief The DLPNO-CCSD pair energy.
  double dlpnoCCSDPairEnergy = 0.0;
  ///@brief The differential overlap of the pair.
  double doi = 0.0;
  ///@brief The PNO-truncation error approximation.
  double deltaPNO = 0.0;
  ///@brief The recovered density matrix trace after PNO truncation.
  double pnoNormError = 0.0;
  ///@brief Number of local fitting functions.
  unsigned int nAuxFunctions = 0;
  ///@brief Orbital pair type. For options, see above.
  OrbitalPairTypes type = OrbitalPairTypes::CLOSE;
  ///@brief Sets representing the orbitals k that couple with this pair.
  std::vector<std::shared_ptr<CouplingOrbitalSet>> coupledPairs;
  ///@brief Sets that represent the pairs k and l that couple with this pair.
  std::vector<std::shared_ptr<KLOrbitalSet>> klPairSets;
  ///@brief Singles corresponding to i.
  std::shared_ptr<SingleSubstitution> singles_i = nullptr;
  ///@brief Singles corresponding to j.
  std::shared_ptr<SingleSubstitution> singles_j = nullptr;

  // Dressed quantities:
  ///@brief Effective double amplitudes.
  Eigen::MatrixXd tau_ij;
  ///@brief Dressed pair--pair interaction part.
  Eigen::MatrixXd y_ij;
  ///@brief Dressed pair--pair interaction part.
  Eigen::MatrixXd y_ji;
  ///@brief Dressed pair--pair interaction part.
  Eigen::MatrixXd z_ij;
  ///@brief Dressed pair--pair interaction part.
  Eigen::MatrixXd z_ji;

  // Fock matrix and dressed fock matrix blocks.
  ///@brief ij fock matrix block.
  double f_ij = 0.0;
  ///@brief Dressed ij fock matrix block.
  double tf_ij = 0.0;
  ///@brief Dressed ji fock matrix block.
  double tf_ji = 0.0;
  ///@brief Dressed ij fock matrix block with singles.
  double ttf_ij = 0.0;
  ///@brief Dressed ji fock matrix block with singles.
  double ttf_ji = 0.0;
  ///@brief ab fock matrix block.
  Eigen::VectorXd f_ab;
  ///@brief Dressed ab fock matrix block.
  Eigen::MatrixXd tf_ab;

  /**
   * @brief Write all integrals to a file.
   */
  void writeIntegralsToFile(HDF5::H5File& file);
  /**
   * @brief Load all integrals from a file.
   */
  void loadIntegralsFromFile(HDF5::H5File& file);
  /**
   * @brief Remove all (ac|bd) and k-Set integrals from main memory.
   */
  void flushIntegrals();
  /**
   * @brief Remoce all integrals and intermediates except K_ij from main memory.
   *
   * Only integrals K_ij and amplitudes t_ij remain.
   */
  void cleanUp();
  /**
   * @brief Returns an estimate of the memory required for the orbital pair.
   * @param mp2Memory Flag if the pair is only used for MP2
   * @param sigmaInts If true, the memory for the integrals used for the sigma vector is included.
   * @return The memory estimate.
   */
  double getMemoryRequirement(bool mp2Memory, bool sigmaInts = false);
  /**
   * @brief Getter for the overlap matrix between virtual domains of this pair and the i single.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_ij_i();
  /**
   * @brief Getter for the overlap matrix between virtual domains of this pair and the j single.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_ij_j();
  /**
   * @brief Getter for the overlap matrix between virtual domains of the i and j singles.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_i_j();
  /**
   * @brief Setter for the domain overlap matrix controller.
   * @param domainSController The domain overlap matrix controller.
   */
  void setOverlapMatrixController(std::shared_ptr<DomainOverlapMatrixController> domainSController);
  /**
   * @brief Getter for the local MP2 pair energy.
   * @return The pair energy.
   */
  double getLMP2PairEnergy();
  /**
   * @brief Getter for the local CCSD pair energy.
   * @return The pair energy.
   */
  double getCCSDPairEnergy();
  /**
   * @brief Getter for the pair energy.
   * @return Returns the CCSD pair energy if available. Otherwise the MP2 pair energy is returned.
   */
  double getPairEnergy();
  /**
   * @brief Getter for the orbital pair specific PNO threshold.
   * @return The PNO threshold.
   */
  double getPNOThreshold();
  /**
   * @brief Getter for the orbital pair specific pair energy truncation threshold.
   * @return The truncation threshold.
   */
  double getPairEnergyThreshold();
  /**
   * @brief Getter for the orbital pair specific truncation threshold for the dipole approximation to the pair energy.
   * @return The truncation threshold.
   */
  double getCollinearDipolePairThreshold();
  /**
   * @brief Getter for the truncation threshold for the "close" pairs used in the sparse map construction.
   * @return The truncation threshold.
   */
  double getSparseMapConstructionPairThreshold();
  /**
   * @brief Getter for the extended fitting domain.
   * @return The fitting domain.
   */
  const Eigen::SparseVector<int>& getFittingDomain();
  /**
   * @brief Setter for the extended fitting domain.
   * @param fittingDomain The extended fitting domain.
   */
  void setFittingDomain(const Eigen::SparseVector<int>& fittingDomain);

 private:
  ///@brief The orbital pair specific PNO threshold.
  double _pnoThreshold;
  ///@brief The CCSD pair truncation threshold.
  double _mainPairEnergyThreshold;
  ///@brief The truncation threshold for the dipole approximation.
  double _collinearDipoleApproxThreshold;
  ///@brief The pair energy threshold used for the sparse map construction.
  double _sparseMapsConstructionPairThreshold;
  /**
   * @brief The file names.
   */
  std::map<FOUR_CENTER_INTEGRAL_TYPE, std::string> _groupNames = {
      {FOUR_CENTER_INTEGRAL_TYPE::ia_jb, "ia_jb"},
      {FOUR_CENTER_INTEGRAL_TYPE::ac_bd, "ac_bd"},
      {FOUR_CENTER_INTEGRAL_TYPE::ka_bc, "ka_bc"},
      {FOUR_CENTER_INTEGRAL_TYPE::ij_ab, "ij_ab"},
      {FOUR_CENTER_INTEGRAL_TYPE::ja_ik, "ja_ik"},
      {FOUR_CENTER_INTEGRAL_TYPE::ia_jk, "ia_jk"},
      {FOUR_CENTER_INTEGRAL_TYPE::ij_ak, "ij_ak"},
      {FOUR_CENTER_INTEGRAL_TYPE::ik_jl, "ik_jl"},
      {FOUR_CENTER_INTEGRAL_TYPE::ia_kc, "ia_kc"},
      {FOUR_CENTER_INTEGRAL_TYPE::ik_ca, "ik_ca"},
      {FOUR_CENTER_INTEGRAL_TYPE::ja_kc, "ja_kc"},
      {FOUR_CENTER_INTEGRAL_TYPE::jk_ca, "jk_ca"},
      {FOUR_CENTER_INTEGRAL_TYPE::ki_la, "ki_la"},
      {FOUR_CENTER_INTEGRAL_TYPE::kj_la, "kj_la"},
      {FOUR_CENTER_INTEGRAL_TYPE::ia_bc, "ia_bc"},
      {FOUR_CENTER_INTEGRAL_TYPE::ja_bc, "ja_bc"},
      {FOUR_CENTER_INTEGRAL_TYPE::jc_ab, "jc_ab"},
      {FOUR_CENTER_INTEGRAL_TYPE::ic_ab, "ic_ab"},
      {FOUR_CENTER_INTEGRAL_TYPE::iaS_jbSX2_M_ij_aSbS, "iaS_jbSX2_M_ij_aSbS"},
      {FOUR_CENTER_INTEGRAL_TYPE::s_ij_kj, "s_ij_kj"},
      {FOUR_CENTER_INTEGRAL_TYPE::s_ij_ik, "s_ij_ik"},
      {FOUR_CENTER_INTEGRAL_TYPE::s_ij_kl, "s_ij_kl"},
      {FOUR_CENTER_INTEGRAL_TYPE::s_ij_i, "s_ij_i"},
      {FOUR_CENTER_INTEGRAL_TYPE::s_ij_j, "s_ij_j"},
      {FOUR_CENTER_INTEGRAL_TYPE::s_ij_k, "s_ij_k"},
      {FOUR_CENTER_INTEGRAL_TYPE::s_i_j, "s_i_j"},
      {FOUR_CENTER_INTEGRAL_TYPE::ab_kcX2_M_ak_bc, "ab_kcX2_M_ak_bc"}};
  /**
   * @brief Flag whether the integrals have been written to disk.
   */
  bool _integralsOnDisk = false;
  /**
   * @brief Write the (ac|bd) integrals to file.
   * @param file The file.
   */
  void write_acbd(HDF5::H5File& file, std::string id);
  /**
   * @brief Load the (ac|bd) integrals from file.
   * @param file The file.
   */
  void load_acbd(HDF5::H5File& file, std::string id);

  ///@brief The PNO overlap between singles and this pair
  std::shared_ptr<Eigen::MatrixXd> _s_ij_i;
  std::shared_ptr<Eigen::MatrixXd> _s_ij_j;
  ///@brief The PNO overlap between the singles
  std::shared_ptr<Eigen::MatrixXd> _s_i_j;
  ///@brief The domain overlap matrix controller.
  std::weak_ptr<DomainOverlapMatrixController> _domainSController;
  ///@brief The pair identification string within the HDF5 file.
  std::string _id = "";
  ///@brief The extended fitting domain.
  Eigen::SparseVector<int> _extendedAuxDomain;
};

} /* namespace Serenity */

#endif /* DATA_ORBITALPAIR_H_ */
