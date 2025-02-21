/**
 * @file   Ao2MoExchangeIntegralTransformer.h
 *
 * @date   Dec. 23, 2018
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

#ifndef AO2MOEXCHANGEINTEGRALTRANSFORMER_H_
#define AO2MOEXCHANGEINTEGRALTRANSFORMER_H_
/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h"           //Definition of MatrixInBasis
#include "integrals/MO3CenterIntegralController.h" //Three center integrals
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dense matrices
#include <Eigen/SparseCore> //Sparse matrices
#include <memory>           //smrt_ptr

namespace Serenity {
/* Forward Declarations */
class BasisController;
class AtomCenteredBasisController;
class OrbitalPair;
class QuasiCanonicalPAODomainConstructor;
class MO3CenterIntegralController;
class OrbitalPairSet;

/**
 * @class Ao2MoExchangeIntegralTransformer Ao2MoExchangeIntegralTransformer.h
 * @brief Performs PAO/AO to PNO/MO transformations for integrals needed in
 *        LMP2 and CCSD. The RI approximation is used though out.
 */
class Ao2MoExchangeIntegralTransformer {
 private:
  /*
   * Purely static. Never instantiated.
   */
  Ao2MoExchangeIntegralTransformer();
  ~Ao2MoExchangeIntegralTransformer();

 public:
  /**
   * @brief Calculates the exchange integrals (ai|bj) via sparse map prescreening.
   * @param auxBasisController The axiliary basis.
   * @param mo3CenterIntegralController The MO 3 center integral controller.
   * @param orbitalPairs The orbital pairs.
   * @param pnoConstructor The PNO/QCPAO constructor
   */
  static void transformExchangeIntegrals(std::shared_ptr<BasisController> auxBasisController,
                                         std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                                         std::vector<std::shared_ptr<OrbitalPair>>& orbitalPairs,
                                         std::shared_ptr<QuasiCanonicalPAODomainConstructor> pnoConstructor);

  /**
   * @brief Calculates all integrals needed for DLPNO-CCSD [except (ia|jb)].
   *        It is assumed that the PNO basis was already constructed.
   * @param auxBasisController The auxiliary basis controller.
   * @param mo3CenterIntegralController The MO 3 center integral controller.
   * @param orbitalPairs The orbital pairs.
   * @param calcSigmaVectorInts Calculate the integrals necessary for the sigma vector construction.
   * @param lowMemory A flag for low memory usage. 3-Center integrals are recalculated more often.
   * @param memoryDemandForSet The memory demand for storing all four center integrals of the orbital pair set.
   * @param ignoreMemoryHandling If true, Serenity assumes that it always has enough memory.
   */
  static void transformAllIntegrals(std::shared_ptr<BasisController> auxBasisController,
                                    std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                                    std::vector<std::shared_ptr<OrbitalPair>>& orbitalPairs, bool calcSigmaVectorInts,
                                    bool lowMemory = false, double memoryDemandForSet = 0.0,
                                    bool ignoreMemoryHandling = false);

  /**
   * @brief Calculates all integrals needed for DLPNO-CCSD [except (ia|jb)].
   *        It is assumed that the PNO basis was already constructed.
   * @param auxBasisController The auxiliary basis controller.
   * @param mo3CenterIntegralController The MO 3 center integral controller.
   * @param orbitalPairSets The orbital pairs as sets.
   * @param dumpIntegrals If true the integrals for each pair are written to disk.
   * @param calcSigmaVectorInts Calculate the integrals necessary for the sigma vector construction.
   * @param lowMemory A flag for low memory usage. 3-Center integrals are recalculated more often.
   * @param ignoreMemoryHandling If true, Serenity assumes that it always has enough memory.
   */
  static void transformAllIntegrals(std::shared_ptr<BasisController> auxBasisController,
                                    std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                                    std::vector<std::shared_ptr<OrbitalPairSet>> orbitalPairSets, bool dumpIntegrals,
                                    std::string pairIntegralFileName, bool calcSigmaVectorInts, bool lowMemory = false,
                                    bool ignoreMemoryHandling = false);

  /**
   * @brief Calculates the two center integrals (K|Q).
   * @param metric An empty MatrixInBasis which defines the
   *               basis set and will be filled.
   */
  static void calculateTwoCenterIntegrals(MatrixInBasis<Options::SCF_MODES::RESTRICTED>& metric);
  /**
   * @brief Calculates the exchange integrals (ai|bj) using full four center integrals.
   *        This is only sensible for testing purposes.
   * @param basisController The basis controller.
   * @param aoCoefficients The coefficients of the occupied orbitals.
   * @param paoController The PAO controller.
   * @param orbitalPairs The orbital pairs.
   * @param pnoConstructor The PNO/QCPAO constructor.
   */
  static void transformExchangeIntegrals(std::shared_ptr<BasisController> basisController,
                                         const Eigen::MatrixXd& aoCoefficients, std::shared_ptr<PAOController> paoController,
                                         std::vector<std::shared_ptr<OrbitalPair>>& orbitalPairs,
                                         std::shared_ptr<QuasiCanonicalPAODomainConstructor> pnoConstructor);

 private:
  /*
   *           ==== A Technical Note on all Integral Evaluations ====
   *
   * All integral functions that calculate four index integrals of type (rs|pq)
   * are structured in the same way:
   *   1. The three center integrals (rs|K) and (pq|K) in PAO/MO basis are
   *      extracted from the sparse lists.
   *   2. The Cholesky decomposition for one of the integral sets is calculated.
   *   3. The three center integrals and decomposed integrals are recombined to
   *      give the four center integral
   *         (rs|pq) = llt[(rs|K)] * (pq|K) = (rs|K)(K|Q)^(-1)(Q|pq) (implicit summations)
   *      All this is done employing fast matrix-matrix multiplications.
   *
   *   Step 1 is handled by the "getter" functions get_abK, get_iaK and get_ikK
   *   which extract integrals over two, one and no virtual index, respectively.
   */

  /**
   * @brief Calculates the integrals with at least two virtual indices:
   *  (ac|bd), a,b,c,d in [ij].
   *  (ka|bc) c in [ij],
   *  (ia|bc) c in [j],
   *  (ja|bc) c in [i],
   *  (ic|ab) c in [j],
   *  (jc|ab) c in [i],
   *  (ij|ab) a,b in [ij]
   *  with a,b in [ij]
   */
  static void
  calculate_2Virt_integrals(std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
                            const Eigen::LLT<Eigen::MatrixXd>& llt_metric, const Eigen::SparseVector<int>& pairDomainToK,
                            const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
                            const MO3CenterIntegrals& m_abKs, const Eigen::SparseMatrix<double>& extendedPAODomainProjection_T,
                            const Eigen::SparseMatrix<double>& extendedAuxProjectionMatrix_T,
                            const Eigen::VectorXd& invV_ijK, const MO3CenterIntegrals& iaK,
                            const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
                            bool calcSigmaVectorInts, const MO3CenterIntegrals& abK);

  /**
   * @brief Calculates the integrals
   *  (ja|ik), (ia|jk), (ki|la) and (kj|la)
   *  with a in [ij].
   */
  static void
  calculate_3Occ_Integrals(std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
                           const Eigen::LLT<Eigen::MatrixXd>& llt_metric, const Eigen::SparseVector<int>& pairDomainToK,
                           const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
                           const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
                           const MO3CenterIntegrals& iaK, const MO3CenterIntegrals& klK);
  /**
   * @brief Calculates the integrals over mixed PNO domains.
   *        Note: This is the most demanding step in the DLPNO-CCSD integral
   *              transformation.
   *  (ik|ca) c in [kj],
   *  (jk|ca) c in [ik],
   *  (ia|kc) c in [kj],
   *  (ja|kc) c in [ik],
   *  for all a in [ij]
   */
  static void
  calculateMixedIntegrals(std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
                          const MatrixInBasis<RESTRICTED>& metric,
                          const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
                          const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
                          const SparseMap& occToK, const MO3CenterIntegrals& iaK, const MO3CenterIntegrals& klK,
                          const Eigen::SparseMatrix<double>& extendedPAODomainProjection_T,
                          const Eigen::SparseVector<int>& extendedAuxDomain,
                          const Eigen::SparseMatrix<double>& extendedAuxProjectionMatrix_T,
                          const std::vector<Eigen::MatrixXd>& mEx_amuQ, bool calcSigmaVectorInts);

  /**
   * @brief Calculates (ij|K)(K|Q)^(-1).
   * @return (ij|K)(K|Q)^(-1).
   */
  static Eigen::VectorXd
  calculate_invV_ijK(std::shared_ptr<OrbitalPair> pair, const Eigen::SparseVector<int>& pairDomainToK,
                     std::shared_ptr<BasisController> auxBasisController, const Eigen::LLT<Eigen::MatrixXd>& llt_metric,
                     const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
                     const MO3CenterIntegrals& klK);
  /**
   * @brief Build a sparse projection matrix from the selected aux. domain to the full domain.
   * @param auxBasisController The auxiliary basis controller.
   * @param auxDomain The auxiliary basis function domain.
   * @return The sparse projection matrix.
   */
  static Eigen::SparseMatrix<double> buildSparseAuxProjection(std::shared_ptr<BasisController> auxBasisController,
                                                              const Eigen::SparseVector<int>& auxDomain);

  /**
   * @brief The integral list (ac|K) with a as the vector index is resorted to c as the vector index and c is
   *        transformed to the given PNO basis.
   * @param ints The initial integral list.
   * @param toPNO The transformation for the second index.
   * @param sparseProjection The PAO projection for the second index.
   * @return The transformed and sorted integrals.
   */
  static MO3CenterIntegrals flipAndTransform(const MO3CenterIntegrals& ints, const Eigen::MatrixXd& toPNO,
                                             const Eigen::SparseMatrix<double>& sparseProjection);
  /**
   * @brief The integral list (ac|K) with a as the vector index is resorted to c as the vector index and c is
   *        transformed to the given PNO basis.
   * @param ints The initial integral list.
   * @param toPNO The transformation for the second index.
   * @param sparseProjection The PAO projection for the second index.
   * @param auxProjection The projection matrix to further reduce the auxiliary domain.
   * @return The transformed and sorted integrals.
   */
  static MO3CenterIntegrals flipAndTransform(const MO3CenterIntegrals& ints, const Eigen::MatrixXd& toPNO,
                                             const Eigen::SparseMatrix<double>& sparseProjection,
                                             const Eigen::SparseMatrix<double>& auxProjection);

 public:
  /**
   * @brief Extracts the integrals (ab|K) from the linear scaling
   *        integral list abK.
   * @param auxBasisController The auxiliary basis controller.
   * @param nLocalAux The number of auxiliary basis functions.
   * @param auxDomain The fitting domain as a shell-wise list.
   * @param a_paoDomain The local PAO domain for the outer dimension.
   * @param b_paoDomain The local PAO domain for the inner dimension.
   * @param toPNO_a The transformation to the PNO domain for the outer dimension (vector length).
   * @param toPNO_b The transformation to the PNO domain for the inner dimension (matrix rows).
   * @param abK The linear scaling integral list (ab|K)
   * @return The (ab|K) integrals over PNOs.
   */
  static std::vector<Eigen::MatrixXd>
  get_abK(std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
          const Eigen::SparseMatrix<double>& projectionMatrix_a, const Eigen::SparseMatrix<double>& projectionMatrix_b,
          const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
          const Eigen::MatrixXd& toPNO_a, const Eigen::MatrixXd& toPNO_b, const MO3CenterIntegrals& abK);
  /**
   * @brief Extracts the integrals (a mu|K) from the linear scaling integral list abK.
   *        Only the index a (vector index) is transformed to the PNO basis.
   * @param auxBasisController The auxiliary basis controller.
   * @param nLocalAux The number of auxiliary basis functions.
   * @param auxDomain The fitting domain as a shell-wise list.
   * @param projectionMatrix_a The local PAO projection matrix for the outer dimension.
   * @param projectionMatrix_b The local PAO projection matrix for the inner dimension.
   * @param k_PAOToFullPAOMaps PAO projection matrices for the integral lists.
   * @param toPNO_a The transformation to the PNO domain for the outer dimension (vector length).
   * @param abK The linear scaling integral list (ab|K)
   * @return The (ab|K) integrals over PNOs.
   */
  static std::vector<Eigen::MatrixXd>
  get_abK(std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
          const Eigen::SparseMatrix<double>& projectionMatrix_a, const Eigen::SparseMatrix<double>& projectionMatrix_b,
          const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
          const Eigen::MatrixXd& toPNO_a, const MO3CenterIntegrals& abK);
  /**
   * @brief Extracts the integrals (a mu|K) from the linear scaling integral list abK
   *        and the integrals (cd|K). For (a mu|K) only the index a (vector index) is
   *        transformed to the PNO basis. For the list (cd|K), the set of Ks have to be
   *        a subset of the aux. functions used for (a mu|K). The same applies for the
   *        PAO sets used to expand the integral sets.
   * @param auxBasisController The auxiliary basis controller.
   * @param nLocalAux The number of auxiliary basis functions.
   * @param auxDomain The fitting domain as a shell-wise list.
   * @param projectionMatrix_a The local PAO projection matrix for the outer dimension.
   * @param projectionMatrix_b The local PAO projection matrix for the inner dimension.
   * @param projectionMatrix_c The local PAO projection matrix for the outer dimension of the second list.
   * @param projectionMatrix_d The local PAO projection matrix for the inner dimension of the second list.
   * @param k_PAOToFullPAOMaps PAO projection matrices for the integral lists.
   * @param toPNO_a The transformation to the PNO domain for the outer dimension (vector length).
   * @param toPNO_c The transformation to the PNO domain for the outer dimension (vector length) second list.
   * @param toPNO_d The transformation to the PNO domain for the inner dimension second list.
   * @param nLocalAuxQ The number of auxiliary basis functions of the second list.
   * @param auxDomainQ The fitting domain as a shell-wise list for the second list.
   * @param abK The linear scaling integral list (ab|K)
   * @param cdK The second integral list to be filled.
   * @return
   */
  static std::vector<Eigen::MatrixXd>
  get_abK(std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux, const Eigen::SparseVector<int>& auxDomain,
          const Eigen::SparseMatrix<double>& projectionMatrix_a, const Eigen::SparseMatrix<double>& projectionMatrix_b,
          const Eigen::SparseMatrix<double>& projectionMatrix_c, const Eigen::SparseMatrix<double>& projectionMatrix_d,
          const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps, const Eigen::MatrixXd& toPNO_a,
          const Eigen::MatrixXd& toPNO_c, const Eigen::MatrixXd& toPNO_d, unsigned int nLocalAuxQ,
          const Eigen::SparseVector<int>& auxDomainQ, const MO3CenterIntegrals& abK, MO3CenterIntegrals& cdK);

  /**
   * @brief Analogous to the function above. Extracts integrals (ia|K).
   * @param auxBasisController The auxiliary basis controller.
   * @param nLocalAux The number of auxiliary basis functions.
   * @param auxDomain The fitting domain as a shell-wise list.
   * @param projectionMatrix_a Projection to the important PAOs for PNOs a.
   * @param k_PAOToFullPAOMaps Projection to the integrals calculated for a given aux. shell.
   * @param reducedOccIndices The indices of the occupied orbitals considered for each aux. shell.
   * @param toPNO_a The transformation to the PNO domain.
   * @param i The occupied orbital index i.
   * @param iaK The (ia|K) integrals over MOs and PAOs.
   * @return The integrals (ia|K) in PNO basis over the local fitting domain.
   */
  static Eigen::MatrixXd get_iaK(std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux,
                                 const Eigen::SparseVector<int>& auxDomain, const Eigen::SparseMatrix<double>& projectionMatrix_a,
                                 const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
                                 const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
                                 const Eigen::MatrixXd& toPNO_a, const unsigned int i, const MO3CenterIntegrals& iaK);
  /**
   * @brief Analogous to the function above. Extracts integrals (ia|K) for a set of i.
   * @param auxBasisController The auxiliary basis controller.
   * @param nLocalAux The number of auxiliary basis functions.
   * @param auxDomain The fitting domain as a shell-wise list.
   * @param projectionMatrix_a Projection to the important PAOs for PNOs a.
   * @param k_PAOToFullPAOMaps Projection to the integrals calculated for a given aux. shell.
   * @param reducedOccIndices The indices of the occupied orbitals considered for each aux. shell.
   * @param toPNO_a The transformation to the PNO domain.
   * @param iIndices The occupied orbital indices i.
   * @param iaK The (ia|K) integrals over MOs and PAOs.
   * @return The integrals (ia|K) in PNO basis over the local fitting domain.
   */
  static std::vector<Eigen::MatrixXd>
  get_iaK(std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux,
          const Eigen::SparseVector<int>& auxDomain, const Eigen::SparseMatrix<double>& projectionMatrix_a,
          const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_PAOToFullPAOMaps,
          const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
          const Eigen::MatrixXd& toPNO_a, const std::map<unsigned int, unsigned int>& iIndices, const MO3CenterIntegrals& iaK);
  /**
   * @brief Analogous to the function above. Extracts integrals (ik|K). Duplicates in k are not possible!
   * @param auxBasisController The auxiliary basis controller.
   * @param nLocalAux The number of auxiliary basis functions.
   * @param auxDomain The fitting domain as a shell-wise list.
   * @param reducedOccIndices The indices of the occupied orbitals considered for each aux. shell.
   * @param kIndices The ks to be extracted as a non-redundant #-map.
   * @param i The occupied orbital index i.
   * @param klK All integrals (ik|K) integrals over MOs.
   * @return The integral list (ik|K) over the local fitting domain.
   */
  static Eigen::MatrixXd get_ikK(std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux,
                                 const Eigen::SparseVector<int>& auxDomain,
                                 const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
                                 const std::map<unsigned int, unsigned int>& kIndices, const unsigned int i,
                                 const MO3CenterIntegrals& klK);
  /**
   * @brief Analogous to the function above. Extracts integrals (ik|K)
   * @param auxBasisController The auxiliary basis controller.
   * @param nLocalAux The number of auxiliary basis functions.
   * @param auxDomain The fitting domain as a shell-wise list.
   * @param reducedOccIndices The indices of the occupied orbitals considered for each aux. shell.
   * @param kIndices The ks to be extracted.
   * @param i The occupied orbital index i.
   * @param klK All integrals (ik|K) integrals over MOs.
   * @return The integral list (ik|K) over the local fitting domain.
   */
  static Eigen::MatrixXd get_ikK(std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux,
                                 const Eigen::SparseVector<int>& auxDomain,
                                 const std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>& reducedOccIndices,
                                 const std::vector<unsigned int> kIndices, const unsigned int i,
                                 const MO3CenterIntegrals& klK);
};
} /* namespace Serenity */
#endif /* AO2MOEXCHANGEINTEGRALTRANSFORMER_H_ */
