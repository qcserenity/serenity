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
typedef std::vector<std::shared_ptr<OrbitalPair>> OrbitalPairSet;

/* Forward Declarations */
class BasisController;
class AtomCenteredBasisController;
class OrbitalPair;
class QuasiCanonicalPAODomainConstructor;
class MO3CenterIntegralController;

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
   */
  static void transformAllIntegrals(std::shared_ptr<BasisController> auxBasisController,
                                    std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                                    std::vector<std::shared_ptr<OrbitalPair>>& orbitalPairs);

  /**
   * @brief Calculates all integrals needed for DLPNO-CCSD [except (ia|jb)].
   *        It is assumed that the PNO basis was already constructed.
   * @param auxBasisController The auxiliary basis controller.
   * @param mo3CenterIntegralController The MO 3 center integral controller.
   * @param orbitalPairSets The orbital pairs as sets.
   * @param dumpIntegrals If true the integrals for each pair are written to disk.
   */
  static void transformAllIntegrals(std::shared_ptr<BasisController> auxBasisController,
                                    std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                                    std::vector<OrbitalPairSet>& orbitalPairSets, bool dumpIntegrals);

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
   * @brief Calculates the integrals (ac|bd), a,b,c,d in [ij].
   */
  static void calculate_acbd_integrals(std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
                                       const Eigen::LLT<Eigen::MatrixXd>& llt_metric,
                                       const Eigen::SparseVector<int>& pairDomainToK,
                                       const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
                                       const MO3CenterIntegrals& abK);
  /**
   * @brief Calculates the integrals:
   *  (ka|bc) c in [ij],
   *  (ia|bc) c in [j],
   *  (ja|bc) c in [i],
   *  (ic|ab) c in [j],
   *  (jc|ab) c in [i]
   *  with a,b in [ij]
   */
  static void calculate_kabcANDib_acANDjb_ac_integrals(std::shared_ptr<OrbitalPair> pair,
                                                       std::shared_ptr<BasisController> auxBasisController,
                                                       const Eigen::LLT<Eigen::MatrixXd>& llt_metric,
                                                       const Eigen::SparseVector<int>& pairDomainToK,
                                                       const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
                                                       const std::vector<Eigen::VectorXi>& reducedOccIndices,
                                                       const MO3CenterIntegrals& abK, const MO3CenterIntegrals& iaK);

  /**
   * @brief Calculates the integrals
   *  (ij|ab) a,b in [ij]
   */
  static void calculate_ijab_integrals(std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
                                       const Eigen::SparseVector<int>& pairDomainToK,
                                       const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
                                       const MO3CenterIntegrals& abK, const Eigen::VectorXd& invV_ijK);

  /**
   * @brief Calculates the integrals
   *  (ja|ik), (ia|jk), (ki|la) and (kj|la)
   *  with a in [ij].
   */
  static void calculate_ikjaANDikjlANDkilaANDkjla_integrals(std::shared_ptr<OrbitalPair> pair,
                                                            std::shared_ptr<BasisController> auxBasisController,
                                                            const Eigen::LLT<Eigen::MatrixXd>& llt_metric,
                                                            const Eigen::SparseVector<int>& pairDomainToK,
                                                            const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
                                                            const std::vector<Eigen::VectorXi>& reducedOccIndices,
                                                            const MO3CenterIntegrals& iaK, const MO3CenterIntegrals& klK);
  /**
   * @brief Calculates the integrals
   *  (ij|ka) with a in [ij]
   */
  static void calculate_ijka_integrals(std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
                                       const Eigen::SparseVector<int>& pairDomainToK,
                                       const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
                                       const std::vector<Eigen::VectorXi>& reducedOccIndices,
                                       const MO3CenterIntegrals& iaK, const Eigen::VectorXd& invV_ijK);
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
  static void calculateMixedIntegrals(std::shared_ptr<OrbitalPair> pair, std::shared_ptr<BasisController> auxBasisController,
                                      const MatrixInBasis<RESTRICTED>& metric,
                                      const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
                                      const std::vector<Eigen::VectorXi>& reducedOccIndices, const SparseMap& occToK,
                                      const MO3CenterIntegrals& iaK, const MO3CenterIntegrals& abK,
                                      const MO3CenterIntegrals& klK);

  /**
   * @brief Calculates (ij|K)(K|Q)^(-1).
   * @return (ij|K)(K|Q)^(-1).
   */
  static Eigen::VectorXd
  calculate_invV_ijK(std::shared_ptr<OrbitalPair> pair, const Eigen::SparseVector<int>& pairDomainToK,
                     std::shared_ptr<BasisController> auxBasisController, const Eigen::LLT<Eigen::MatrixXd>& llt_metric,
                     const std::vector<Eigen::VectorXi>& reducedOccIndices, const MO3CenterIntegrals& klK);
  /**
   * @brief Build a sparse projection matrix from the selected aux. domain to the full domain.
   * @param auxBasisController The auxiliary basis controller.
   * @param auxDomain The auxiliary basis function domain.
   * @return The sparse projection matrix.
   */
  static Eigen::SparseMatrix<double> buildSparseAuxProjection(std::shared_ptr<BasisController> auxBasisController,
                                                              const Eigen::SparseVector<int>& auxDomain);

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
          const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps, const Eigen::MatrixXd& toPNO_a,
          const Eigen::MatrixXd& toPNO_b, const MO3CenterIntegrals& abK);
  /**
   * @brief Extracts the integrals (ab|K) and (ac|K) from the linear scaling
   *        integral list abK, where both integral sets are represented over
   *        the same fitting domain.
   *
   *  Note: It is slightly faster to extract the (ab|K) and (ac|K) integrals
   *        with this functions if both sets are needed instead of calling the
   *        get_abK(...) function above twice.
   *
   * @param auxBasisController The auxiliary basis controller.
   * @param nLocalAux The number of auxiliary basis functions.
   * @param auxDomain The fitting domain as a shell-wise list.
   * @param projectionMatrix_a Projection to the important PAOs for PNOs a.
   * @param projectionMatrix_b Projection to the important PAOs for PNOs b.
   * @param projectionMatrix_c Projection to the important PAOs for PNOs c.
   * @param k_PAOToFullPAOMaps Projection to the integrals calculated for a given aux. shell.
   * @param toPNO_a The transformation to the PNO domain for the outer dimension (vector length).
   * @param toPNO_b The transformation to the PNO domain for the inner dimension (matrix rows) of index b.
   * @param toPNO_c The transformation to the PNO domain for the inner dimension (matrix rows) of index c.
   * @param abK The linear scaling integral list (ab|K)
   * @return The (ab|K) and (ac|K) integrals over PNOs.
   *         pair.first-->(ab|K) and pair.second-->(ac|K)
   */
  static std::pair<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>>
  get_abK(std::shared_ptr<BasisController> auxBasisController, unsigned int nLocalAux,
          const Eigen::SparseVector<int>& auxDomain, const Eigen::SparseMatrix<double>& projectionMatrix_a,
          const Eigen::SparseMatrix<double>& projectionMatrix_b, const Eigen::SparseMatrix<double>& projectionMatrix_c,
          const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps, const Eigen::MatrixXd& toPNO_a,
          const Eigen::MatrixXd& toPNO_b, const Eigen::MatrixXd& toPNO_c, const MO3CenterIntegrals& abK);
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
                                 const std::vector<Eigen::SparseMatrix<double>>& k_PAOToFullPAOMaps,
                                 const std::vector<Eigen::VectorXi>& reducedOccIndices, const Eigen::MatrixXd& toPNO_a,
                                 const unsigned int i, const MO3CenterIntegrals& iaK);
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
                                 const Eigen::SparseVector<int>& auxDomain, const std::vector<Eigen::VectorXi>& reducedOccIndices,
                                 std::vector<unsigned int> kIndices, const unsigned int i, const MO3CenterIntegrals& klK);
};
} /* namespace Serenity */
#endif /* AO2MOEXCHANGEINTEGRALTRANSFORMER_H_ */
