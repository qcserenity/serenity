/**
 * @file OrbitalTriple.h
 *
 * @author Moritz Bensberg
 * @date Oct 31, 2019
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

#ifndef DATA_ORBITALTRIPLE_H_
#define DATA_ORBITALTRIPLE_H_
/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h" // Parsing of the coulomb metric for the integral calculation.
/* Include Std and External Headers */
#include <Eigen/Dense>      // Dense matrices.
#include <Eigen/SparseCore> // Sparse matrices.
#include <memory>           // smrt_ptr
#include <vector>           // std::vector

namespace Serenity {
typedef std::vector<Eigen::MatrixXd> MO3CenterIntegrals;
/* Forward Declartaions */
class OrbitalPair;
class SingleSubstitution;
class MO3CenterIntegralController;

/**
 * @class OrbitalTriple OrbitalTriple.h
 * @brief A class that is concerned with the triple-wise triples correction for DLPNO-CCSD(T0).
 *        It calculates the necessary integrals for the semi-canonical energy correction
 *        increment as well as calculates the latter.\n\n
 *
 *    References for the semi-canonical triples correction:\n
 *      J. Chem. Phys. 139, 134101 (2013)\n
 *      Chem. Phys. Lett. 178, 462 (1991)\n
 *      or the pdf document concerned with DLPNO-CCSD(T0)
 */
class OrbitalTriple {
 public:
  /**
   * @brief Constructor. An orbital triple is given by the three pairs ij, ik and jk.
   * @param ikPair   The ik-pair.
   * @param jkPair   The jk-pair.
   * @param ijPair   The ij-pair.
   * @param ilPairs  List of il-pairs.
   * @param jlPairs  List of jl-pairs.
   * @param klPairs  List of kl-pairs.
   * @param i        Triples index i.
   * @param j        Triples index j.
   * @param k        Triples index k.
   */
  OrbitalTriple(std::shared_ptr<OrbitalPair> ikPair, std::shared_ptr<OrbitalPair> jkPair, std::shared_ptr<OrbitalPair> ijPair,
                std::vector<std::shared_ptr<OrbitalPair>> ilPairs, std::vector<std::shared_ptr<OrbitalPair>> jlPairs,
                std::vector<std::shared_ptr<OrbitalPair>> klPairs, unsigned int i, unsigned int j, unsigned int k);

  /**
   * @brief Getter for semi-canonical triples correction increment.
   * @return The semi-canonical energy correction increment.
   */
  double getTripleEnergy() {
    if (!_triplesCorrection)
      _triplesCorrection = std::make_shared<double>(calculateTripleEnergy());
    return *_triplesCorrection;
  }
  /**
   * @brief Getter for the PAO domain.
   * @return The PAO domain as a sparse list.
   *         1 entry if the PAO is selected.
   *         0 entry else.
   */
  const Eigen::SparseVector<int>& getPAODomain() {
    assert(_paoDomain.nonZeros() != 0);
    return _paoDomain;
  }

  /**
   * @brief Getter for the ij-pair.
   * @return The ij-pair.
   */
  std::shared_ptr<OrbitalPair> getIJPair() {
    return _ijPair;
  }
  /**
   * @brief Getter for the ik-pair.
   * @return The ik-pair.
   */
  std::shared_ptr<OrbitalPair> getIKPair() {
    return _ikPair;
  }
  /**
   * @brief Getter for the jk-pair.
   * @return The jk-pair.
   */
  std::shared_ptr<OrbitalPair> getJKPair() {
    return _jkPair;
  }
  /**
   * @brief Getter for the kl-pairs.
   * @return The kl-pairs.
   */
  std::vector<std::shared_ptr<OrbitalPair>> getKLPairs() {
    return _klPairs;
  }
  /**
   * @brief Getter for the il-pairs.
   * @return The il-pairs.
   */
  std::vector<std::shared_ptr<OrbitalPair>> getILPairs() {
    return _ilPairs;
  }
  /**
   * @brief Getter for the jl-pairs.
   * @return The jl-pairs.
   */
  std::vector<std::shared_ptr<OrbitalPair>> getJLPairs() {
    return _jlPairs;
  }
  /**
   * @brief Getter for the i-singles.
   * @return The i-singles.
   */
  std::shared_ptr<SingleSubstitution> getISingle() {
    return _iSingle;
  }
  /**
   * @brief Getter for the j-singles.
   * @return The j-singles.
   */
  std::shared_ptr<SingleSubstitution> getJSingle() {
    return _jSingle;
  }
  /**
   * @brief Getter for the k-singles.
   * @return The k-singles.
   */
  std::shared_ptr<SingleSubstitution> getKSingle() {
    return _kSingle;
  }
  /**
   * @brief Getter for the i-index.
   * @return The i-index.
   */
  unsigned int getI() {
    return _i;
  }
  /**
   * @brief Getter for the j-index.
   * @return The j-index.
   */
  unsigned int getJ() {
    return _j;
  }
  /**
   * @brief Getter for the k-index.
   * @return The k-index.
   */
  unsigned int getK() {
    return _k;
  }
  /**
   * @brief Default destructor.
   */
  virtual ~OrbitalTriple() = default;
  /**
   * @brief Setter for the TNO basis, pseudo eigenvalues and projections.
   * @param trafo The TNO coefficients in the redundant PAO basis.
   * @param eps The pseudo-eigenvalues.
   * @param s_ij_ijk The overlap matrix between ij-PNOs and ijk-TNOs.
   * @param s_ik_ijk The overlap matrix between ik-PNOs and ijk-TNOs.
   * @param s_jk_ijk The overlap matrix between jk-PNOs and ijk-TNOs.
   * @param s_il_ijk The overlap matrices between il-PNOs and ijk-TNOs.
   * @param s_kl_ijk The overlap matrices between kl-PNOs and ijk-TNOs.
   * @param s_jl_ijk The overlap matrices between jl-PNOs and ijk-TNOs.
   * @param s_i_ijk  The overlap matrix between i-PNOs and ijk-TNOs.
   * @param s_j_ijk  The overlap matrix between j-PNOs and ijk-TNOs.
   * @param s_k_ijk  The overlap matrix between k-PNOs and ijk-TNOs.
   * @param paoProjection The PAO projection matrix.
   */
  void setTNOTransformation(Eigen::MatrixXd trafo, Eigen::VectorXd eps, Eigen::MatrixXd s_ij_ijk,
                            Eigen::MatrixXd s_ik_ijk, Eigen::MatrixXd s_jk_ijk, std::vector<Eigen::MatrixXd> s_il_ijk,
                            std::vector<Eigen::MatrixXd> s_kl_ijk, std::vector<Eigen::MatrixXd> s_jl_ijk,
                            Eigen::MatrixXd s_i_ijk, Eigen::MatrixXd s_j_ijk, Eigen::MatrixXd s_k_ijk,
                            Eigen::SparseMatrix<double> paoProjection) {
    assert(!_tnoInit);
    _toPAODomain = trafo;
    _tnoEigenvalues = eps;
    _tnoInit = true;
    _s_ij_ijk = s_ij_ijk;
    _s_ik_ijk = s_ik_ijk;
    _s_jk_ijk = s_jk_ijk;
    _s_il_ijk = s_il_ijk;
    _s_kl_ijk = s_kl_ijk;
    _s_jl_ijk = s_jl_ijk;
    _s_i_ijk = s_i_ijk;
    _s_j_ijk = s_j_ijk;
    _s_k_ijk = s_k_ijk;
    _domainProjection = paoProjection;
  }
  /**
   * @brief Getter for the TNO coefficients in the redundant PAO basis.
   * @return The TNO coefficients.
   */
  const Eigen::MatrixXd& getTNOCoefficients() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _toPAODomain;
  }
  /**
   * @brief Getter for the TNO eigenvalues.
   * @return The TNO pseudo eigenvalues.
   */
  const Eigen::VectorXd& getTNOEigenvalues() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _tnoEigenvalues;
  }
  /**
   * @brief Calculates the integrals for the given triple.
   * @param auxBasisController The aux. basis controller.
   * @param mo3CenterIntegralController The MO 3 center integral controller.
   * @param coulombMetric The aux. coulomb metric.
   */
  void calculateIntegrals(std::shared_ptr<BasisController> auxBasisController,
                          std::shared_ptr<MO3CenterIntegralController> mo3CenterIntegralController,
                          const MO3CenterIntegrals& iaK, const MO3CenterIntegrals& abK, const MO3CenterIntegrals& klK,
                          const MatrixInBasis<RESTRICTED>& coulombMetric);
  /**
   * @brief Getter for the integrals (ia|bc).
   *        Note: a,b,c,d in [ijk] over TNOs.
   * @return The integrals (ia|bc)
   */
  const std::vector<Eigen::MatrixXd>& get_ia_bc() {
    assert(_integralsCalc);
    return _ia_bc;
  }
  /**
   * @brief Getter for the integrals (ja|bc).
   *        Note: a,b,c,d in [ijk] over TNOs.
   * @return The integrals.
   */
  const std::vector<Eigen::MatrixXd>& get_ja_bc() {
    assert(_integralsCalc);
    return _ja_bc;
  }
  /**
   * @brief Getter for the integrals (ka|bc).
   *        Note: a,b,c,d in [ijk] over TNOs.
   * @return The integrals.
   */
  const std::vector<Eigen::MatrixXd>& get_ka_bc() {
    assert(_integralsCalc);
    return _ka_bc;
  }
  /**
   * @brief Getter for the integrals (ja|il).
   *        Note: a,b,c,d in [ijk] over TNOs,
   *        l in weak/strong pair kl exists.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ja_il() {
    assert(_integralsCalc);
    return _ja_il;
  }
  /**
   * @brief Getter for the integrals (ia|jl).
   *        Note: a,b,c,d in [ijk] over TNOs,
   *        l in weak/strong pair kl exists.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ia_jl() {
    assert(_integralsCalc);
    return _ia_jl;
  }
  /**
   * @brief Getter for the integrals (ka|il).
   *        Note: a,b,c,d in [ijk] over TNOs,
   *        l in weak/strong pair jl exists.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ka_il() {
    assert(_integralsCalc);
    return _ka_il;
  }
  /**
   * @brief Getter for the integrals (ia|kl).
   *        Note: a,b,c,d in [ijk] over TNOs,
   *        l in weak/strong pair jl exists.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ia_kl() {
    assert(_integralsCalc);
    return _ia_kl;
  }
  /**
   * @brief Getter for the integrals (ka|jl).
   *        Note: a,b,c,d in [ijk] over TNOs,
   *        l in weak/strong pair il exists.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ka_jl() {
    assert(_integralsCalc);
    return _ka_jl;
  }
  /**
   * @brief Getter for the integrals (ja|kl).
   *        Note: a,b,c,d in [ijk] over TNOs,
   *        l in weak/strong pair kl exists.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ja_kl() {
    assert(_integralsCalc);
    return _ja_kl;
  }
  /**
   * @brief Getter for the integrals (ia|jb).
   *        Note: a,b,c,d in [ijk] over TNOs.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ia_jb() {
    assert(_integralsCalc);
    return _ia_jb;
  }
  /**
   * @brief Getter for the integrals (ia|kb).
   *        Note: a,b,c,d in [ijk] over TNOs.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ia_kb() {
    assert(_integralsCalc);
    return _ia_kb;
  }
  /**
   * @brief Getter for the integrals (ja|kb).
   *        Note: a,b,c,d in [ijk] over TNOs.
   * @return The integrals.
   */
  const Eigen::MatrixXd& get_ja_kb() {
    assert(_integralsCalc);
    return _ja_kb;
  }
  /**
   * @brief Getter for strong/weak pair indices il.
   * @return The indices.
   */
  const std::vector<unsigned int>& getILIndices() {
    return _ilIndices;
  }
  /**
   * @brief Getter for strong/weak pair indices jl.
   * @return The indices.
   */
  const std::vector<unsigned int>& getJLIndices() {
    return _jlIndices;
  }
  /**
   * @brief Getter for strong/weak pair indices kl.
   * @return The indices.
   */
  const std::vector<unsigned int>& getKLIndices() {
    return _klIndices;
  }
  /**
   * @brief Getter for the overlap matrix ij-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_ij_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_ij_ijk;
  }
  /**
   * @brief Getter for the overlap matrix ik-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_ik_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_ik_ijk;
  }
  /**
   * @brief Getter for the overlap matrix jk-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_jk_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_jk_ijk;
  }
  /**
   * @brief Getter for the overlap matrix i-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_i_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_i_ijk;
  }
  /**
   * @brief Getter for the overlap matrix j-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_j_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_j_ijk;
  }
  /**
   * @brief Getter for the overlap matrix k-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_k_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_k_ijk;
  }
  /**
   * @brief Getter for the overlap matrix il-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const std::vector<Eigen::MatrixXd>& getS_il_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_il_ijk;
  }
  /**
   * @brief Getter for the overlap matrix kl-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const std::vector<Eigen::MatrixXd>& getS_kl_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_kl_ijk;
  }
  /**
   * @brief Getter for the overlap matrix jl-PNOs to ijk-TNOs.
   * @return The overlap matrix.
   */
  const std::vector<Eigen::MatrixXd>& getS_jl_ijk() {
    if (!_tnoInit)
      throw SerenityError("Call to TNO based method before initializing TNOs");
    return _s_jl_ijk;
  }
  /**
   * @brief Getter for the triples type.
   * @return Returns true if the triple is a weak triple. False, otherwise.
   */
  bool isWeak();
  /**
   * @brief Getter for the number of auxiliary functions used in the triple.
   *        This is initialized during the integral calculation.
   * @return The number of auxiliary functions.
   */
  unsigned int getNLocalAuxiliaryFunctions();
  /**
   * @brief Getter for the number of TNOs used.
   * @return The number of TNOs.
   */
  unsigned int getNTNOs();
  /**
   * @brief Setter for the triples specific fitting domain.
   * @param fittingDomain The fittinge domain.
   */
  void setFittingDomain(const Eigen::SparseVector<int>& fittingDomain);
  /**
   * @brief Getter for the triple specific fitting domain.
   * @return The fitting domain.
   */
  const Eigen::SparseVector<int>& getFittingDomain() const;

 private:
  // The overlap matrices. See the associated getter functions.
  Eigen::MatrixXd _s_ij_ijk;
  std::vector<Eigen::MatrixXd> _s_il_ijk;
  Eigen::MatrixXd _s_ik_ijk;
  std::vector<Eigen::MatrixXd> _s_kl_ijk;
  Eigen::MatrixXd _s_jk_ijk;
  std::vector<Eigen::MatrixXd> _s_jl_ijk;
  Eigen::MatrixXd _s_i_ijk;
  Eigen::MatrixXd _s_j_ijk;
  Eigen::MatrixXd _s_k_ijk;
  /**
   * @brief Calculates the semi-canonical triple correction.
   * @return
   */
  double calculateTripleEnergy();
  /**
   * @brief Removes all integrals from memory.
   */
  void cleanUp();
  /*
   * All references to equations correspond to J. Chem. Phys. 139, 134101 (2013).
   */
  // Calculates the main intermediate from Eq. (8).
  std::vector<Eigen::MatrixXd> calculateW();
  // Calculates Eq. (5)
  std::vector<Eigen::MatrixXd> calculateV(const std::vector<Eigen::MatrixXd>& W);
  // Calculates Eq. (4)
  std::vector<Eigen::MatrixXd> calculateX(const std::vector<Eigen::MatrixXd>& W, const std::vector<Eigen::MatrixXd>& V);
  // Calculates Eq. (2)
  std::vector<Eigen::MatrixXd> calculateY(const std::vector<Eigen::MatrixXd>& V);
  // Calculates Eq. (3)
  std::vector<Eigen::MatrixXd> calculateZ(const std::vector<Eigen::MatrixXd>& V);
  // Construct the map occupied orbital to list index in indices map.
  std::map<unsigned int, unsigned int> getIndexMaps(std::vector<unsigned int>& indices);

  // The PAO domain.
  Eigen::SparseVector<int> _paoDomain;
  // The auxiliary fitting domain.
  Eigen::SparseVector<int> _fittingDomain;
  // The projection matrix to the PAO domain.
  Eigen::SparseMatrix<double> _domainProjection;
  // TNO coefficients in the redundant PAO basis.
  Eigen::MatrixXd _toPAODomain;
  // TNO pseudo eigenvalues.
  Eigen::VectorXd _tnoEigenvalues;

  /*
   * Integrals.
   *   a,b,c,d in [ijk] over TNOs.
   * For the integrals (ia|jl) all l are included for which a close or distant pair
   * kl exists.
   */
  // Three virtual indices.
  std::vector<Eigen::MatrixXd> _ia_bc;
  std::vector<Eigen::MatrixXd> _ja_bc;
  std::vector<Eigen::MatrixXd> _ka_bc;
  // One virtual index.
  Eigen::MatrixXd _ja_il;
  Eigen::MatrixXd _ia_jl;
  Eigen::MatrixXd _ka_il;
  Eigen::MatrixXd _ia_kl;
  Eigen::MatrixXd _ka_jl;
  Eigen::MatrixXd _ja_kl;
  // Two virtual indices.
  Eigen::MatrixXd _ia_jb;
  Eigen::MatrixXd _ia_kb;
  Eigen::MatrixXd _ja_kb;
  // Shared_ptr. The orbital pair will never hold a
  // pointer to a triplet.
  // The ik-pair.
  std::shared_ptr<OrbitalPair> _ikPair;
  // The jk-pair.
  std::shared_ptr<OrbitalPair> _jkPair;
  // The ij-pair.
  std::shared_ptr<OrbitalPair> _ijPair;
  // The il-pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _ilPairs;
  // The jl-pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _jlPairs;
  // The kl-pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _klPairs;
  // The i-single.
  std::shared_ptr<SingleSubstitution> _iSingle;
  // The j-single.
  std::shared_ptr<SingleSubstitution> _jSingle;
  // The k-single.
  std::shared_ptr<SingleSubstitution> _kSingle;
  // The occ--occ Fock matrix elements.
  double _f_ii;
  double _f_jj;
  double _f_kk;
  // The l indices for il/jl and kl occupied orbital indices.
  std::vector<unsigned int> _ilIndices;
  std::vector<unsigned int> _jlIndices;
  std::vector<unsigned int> _klIndices;
  // Flag for a weak triple (at least one pair ij,ik or jk is not close).
  bool _weakTriple = false;
  // The occupied orbital indices.
  unsigned int _i;
  unsigned int _j;
  unsigned int _k;
  // The triples energy increment.
  std::shared_ptr<double> _triplesCorrection = nullptr;
  // Debug flags to check for TNO initialization and integral transformation.
  bool _tnoInit = false;
  bool _integralsCalc = false;
  // Number of local auxiliary basis functions.
  unsigned int _nLocalAux;
};

} /* namespace Serenity */

#endif /* DATA_ORBITALTRIPLE_H_ */
