/**
 * @file CouplingOrbitalSet.h
 *
 * @date 1 Jun 2019
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

#ifndef POSTHF_LOCALCORRELATION_COUPLINGORBITALSET_H_
#define POSTHF_LOCALCORRELATION_COUPLINGORBITALSET_H_

/* Include Serenity Internal Headers */
#include "data/OrbitalPair.h" //Definition of OrbitalPair.
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices
#include <memory>      //smrt_ptr

namespace H5 {
class H5File;
} // namespace H5

namespace Serenity {
namespace HDF5 {
using H5File = H5::H5File;
} // namespace HDF5
class DomainOverlapMatrixController;
class SingleSubstitution;

/**
 * @class CouplingOrbitalSet CouplingOrbitalSet.h
 * @brief A class that represents a ik/kj set for a given
 *        orbital pair ij and occupied orbital with index k.\n
 *        This is constantly needed in both DLPNO-MP2 and DLPNO-CCSD.\n
 *        The class holds the associated integrals and provides access
 *        to the underlying pairs/singles.
 */
class CouplingOrbitalSet {
 public:
  /**
   * @brief Constructor.
   * @param parentPair The ij-pair.
   * @param ikPair     The orbital pair ik or ki.
   * @param kjPair     The orbital pair kj or jk.
   * @param k          The occupied orbital index k.
   */
  CouplingOrbitalSet(std::shared_ptr<OrbitalPair> parentPair, std::shared_ptr<OrbitalPair> ikPair,
                     std::shared_ptr<OrbitalPair> kjPair, unsigned int k);
  /**
   * @brief Default destructor.
   */
  ~CouplingOrbitalSet() = default;
  /* ===== INTEGRALS ===== */
  Eigen::MatrixXd ia_kc; // K^ik_a_ij c_kj
  Eigen::MatrixXd ja_kc; // K^jk_a_ij c_ik
  Eigen::MatrixXd ik_ca; // J^ik_ c_kj a_ij
  Eigen::MatrixXd jk_ca; // J^jk_ c_ik a_ij

  Eigen::VectorXd ij_ak; //(ij|ka) a in [ij]
  Eigen::VectorXd ja_ik; //(ja|ik) a in [ij]
  Eigen::VectorXd ia_jk; //(jk|ia) a-->rows

  Eigen::VectorXd ij_akX2_M_ia_jk; // 2(ij|ak) - (ia|jk), a in [k]
  Eigen::VectorXd ij_akX2_M_ja_ik; // 2(ij|ak) - (ja|ik), a in [k]

  std::vector<Eigen::MatrixXd> ka_bc; //(ka|bc) struc: c * a x b || a,b,c in [ij]
  // 2(ab|kc)-(ak|bc) This is currently only used for tests.
  std::vector<Eigen::MatrixXd> ab_kcX2_M_ak_bc; // 2(ab|kc) - (ak|bc), a,b in [ij], c in [k] struc: c * a x b
  /* ===================== */
  /**
   * @brief Getter for the k-singles.
   * @return The k-singles.
   */
  inline std::shared_ptr<SingleSubstitution> getKSingles() {
    auto ikPair = _ikPair.lock();
    return (ikPair->i == _k) ? ikPair->singles_i : ikPair->singles_j;
  }
  /**
   * @brief Getter for the ik-pair. Returns nullptr if ik == ij.
   * @return The ik-pair.
   */
  inline std::shared_ptr<OrbitalPair> getIKPairScreened() {
    return (!_ik_eq_ij) ? _ikPair.lock() : nullptr;
  }
  /**
   * @brief Getter for the ik-pair.
   * @return The ik-pair.
   */
  inline std::shared_ptr<OrbitalPair> getIKPair() {
    return _ikPair.lock();
  }
  /**
   * @brief Getter for the kj-pair.
   * @return The kj-pair.
   */
  inline std::shared_ptr<OrbitalPair> getKJPair() {
    return _kjPair.lock();
  }
  /**
   * @brief Getter for the kj-pair. Return nullptr if kj == ij.
   * @return The kj-pair.
   */
  inline std::shared_ptr<OrbitalPair> getKJPairScreened() {
    return (!_kj_eq_ij) ? _kjPair.lock() : nullptr;
  }
  /**
   * @brief Getter for the index of the orbital k
   * @return The index k.
   */
  inline unsigned int getK() {
    return _k;
  }

  /**
   * @brief Getter for the overlap matrix between ij pair and ik pair.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_ij_ik();
  /**
   * @brief Getter for the overlap matrix between ij pair and kj pair.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_ij_kj();
  /**
   * @brief Getter for the overlap matrix between ij pair and k singles.
   * @return The overlap matrix.
   */
  const Eigen::MatrixXd& getS_ij_k();

  /**
   * @brief Getter for the domain overlap between ik and kj.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_ik_kj();
  /**
   * @brief Getter for the domain overlap between kj and ik.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_kj_ik();
  /**
   * @brief Getter for the domain overlap between ik and j.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_ik_j();
  /**
   * @brief Getter for the domain overlap between kj and j.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_kj_j();
  /**
   * @brief Getter for the domain overlap between kj and i.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_kj_i();
  /**
   * @brief Getter for the domain overlap between kj and j.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_ik_i();
  /**
   * @brief Getter for the domain overlap between j and k.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_j_k();
  /**
   * @brief Getter for the domain overlap between k and j.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_i_k();
  /**
   * @brief Getter for the domain overlap between kj and k.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_kj_k();
  /**
   * @brief Getter for the domain overlap between ik and k.
   * @return The domain overlap.
   */
  const Eigen::MatrixXd& getS_ik_k();

  /**
   * @brief Write the integrals to disk.
   * @param file The target file.
   * @param groupNames The names of the integrals.
   * @param id The integral set id.
   */
  void writeIntegralsToFile(HDF5::H5File& file, std::map<FOUR_CENTER_INTEGRAL_TYPE, std::string>& groupNames, std::string id);
  /**
   * @brief Load the integrals from disk.
   * @param file The source file.
   * @param groupNames The HDF5-names within the file.
   * @param id The integral set id.
   */
  void loadIntegralsFromFile(HDF5::H5File& file, std::map<FOUR_CENTER_INTEGRAL_TYPE, std::string>& groupNames, std::string id);
  /**
   * @brief Removes all integrals from main memory.
   */
  void flushIntegrals();
  /**
   * @brief Get the memory requirement for holding this object.
   * @return The memeory requirement.
   */
  double getMemoryRequirement(bool sigmaInts);
  /**
   * @brief Setter for the domain overlap matrix controller.
   * @param domainSController
   */
  void setOverlapMatrixController(std::shared_ptr<DomainOverlapMatrixController> domainSController);

 private:
  // Overlap matrices between PNO domains.
  std::shared_ptr<Eigen::MatrixXd> _s_ij_kj;
  std::shared_ptr<Eigen::MatrixXd> _s_ij_ik;
  std::shared_ptr<Eigen::MatrixXd> _s_ij_k;
  std::shared_ptr<Eigen::MatrixXd> _s_ik_kj;
  std::shared_ptr<Eigen::MatrixXd> _s_kj_ik;
  std::shared_ptr<Eigen::MatrixXd> _s_ik_j;
  std::shared_ptr<Eigen::MatrixXd> _s_kj_i;
  std::shared_ptr<Eigen::MatrixXd> _s_kj_j;
  std::shared_ptr<Eigen::MatrixXd> _s_ik_i;
  std::shared_ptr<Eigen::MatrixXd> _s_i_k;
  std::shared_ptr<Eigen::MatrixXd> _s_j_k;

  ///@brief The k index.
  unsigned int _k;
  ///@brief The ij-pair.
  const std::weak_ptr<OrbitalPair> _ijPair;
  ///@brief The ik-pair.
  const std::weak_ptr<OrbitalPair> _ikPair;
  ///@brief The kj-pair.
  const std::weak_ptr<OrbitalPair> _kjPair;
  ///@brief The domain overlap matrix controller.
  std::weak_ptr<DomainOverlapMatrixController> _domainSController;
  ///@brief A flag whether ik==ij.
  bool _ik_eq_ij;
  ///@brief A flag whether kj==ij.
  bool _kj_eq_ij;
  /**
   * @brief Write a vector of matrices to file.
   * @param file The file.
   * @param ints The integral vector.
   * @param name The integral name.
   */
  void write_vectorSet(HDF5::H5File& file, std::vector<Eigen::MatrixXd>& ints, std::string name);
  /**
   * @brief Load a vector of matrices from file.
   * @param file The file.
   * @param ints The integral vector to be filled.
   * @param name The interal name.
   * @param nRows The number of rows of the matrices.
   * @param nCols The number of columns of the matrices.
   */
  void load_vectorSet(HDF5::H5File& file, std::vector<Eigen::MatrixXd>& ints, std::string name, unsigned int nRows,
                      unsigned int nCols);
  ///@brief a flag to load the sigma vector integrals if they were written to disk.
  bool _loadSinglesSigmaInts = false;
};

} /* namespace Serenity */

#endif /* POSTHF_LOCALCORRELATION_COUPLINGORBITALSET_H_ */
