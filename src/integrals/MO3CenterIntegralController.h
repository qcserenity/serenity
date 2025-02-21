/**
 * @file MO3CenterIntegralController.h
 *
 * @date May 7, 2019
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
#ifndef INTEGRALS_MO3CENTERINTEGRALCONTROLLER_H_
#define INTEGRALS_MO3CENTERINTEGRALCONTROLLER_H_

/* Include Serenity Internal Headers */
#include "data/SparseMapsController.h" //Getter for the sparse maps controller.
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dense matrices
#include <Eigen/SparseCore> //Sparse matrices
#include <memory>           //smrt_ptr
#include <string>           //std::string

namespace Serenity {

typedef std::vector<Eigen::MatrixXd> MO3CenterIntegrals;

/* Forward Declarations */
class BasisController;
class PAOController;

/**
 * @brief The various types of MO 3 center integrals.
 *
 * In this notation i,k and l denote occupied MOs, a and b virtual MOs
 * and K an auxiliary function. The integrals can be understood as\n
 *     \f$ \int \int \mathrm{d}r_1 \mathrm{d}r_2 i(r_1) a(r_1) \frac{1}{r_{12}} K(r_2). \f$
 */

enum class MO3CENTER_INTS { ia_K, kl_K, ab_K };
/**
 * @brief Classification of the two non-auxiliary functions: i,k--> Occupied, a,b-->virtual.
 */
enum class ORBITAL_TYPE { OCCUPIED, VIRTUAL };

/**
 * @class MO3CenterIntegralController MO3CenterIntegralController.h
 * @brief Calculates and stores 3 center MO integrals.
 */
class MO3CenterIntegralController {
 public:
  /**
   * @brief Constructor.
   * @param auxiliaryBasisController The auxiliary basis controller.
   * @param basisController The AO basis controller.
   * @param sparseMaps Sparse maps controller that holds all prescreening information.
   * @param paoController The PAO controller.
   * @param occupiedCoefficients The occupied coefficients.
   * @param fBaseName The base name of the associated system.
   * @param id The id of the system.
   * @param triplesMode If true, the controller will use functions dedicated to calculating integrals for the triples.
   */
  MO3CenterIntegralController(std::shared_ptr<BasisController> auxiliaryBasisController,
                              std::shared_ptr<BasisController> basisController,
                              const std::shared_ptr<SparseMapsController> sparseMaps,
                              std::shared_ptr<PAOController> paoController, std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients,
                              std::string fBaseName, std::string id, bool triplesMode = false);
  /**
   * @brief Constructor.
   * @param auxiliaryBasisController The auxiliary basis controller.
   * @param basisController The AO basis controller.
   * @param sparseMaps Sparse maps controller that holds all prescreening information.
   * @param virtualCoefficients The virtual orbitals coefficients.
   * @param occupiedCoefficients The occupied coefficients.
   * @param fBaseName The base name of the associated system.
   * @param id The id of the system.
   * @param triplesMode If true, the controller will use functions dedicated to calculating integrals for the triples.
   */
  MO3CenterIntegralController(std::shared_ptr<BasisController> auxiliaryBasisController,
                              std::shared_ptr<BasisController> basisController,
                              const std::shared_ptr<SparseMapsController> sparseMaps,
                              std::shared_ptr<Eigen::MatrixXd> virtualCoefficients,
                              std::shared_ptr<Eigen::MatrixXd> occupiedCoefficients, std::string fBaseName,
                              std::string id, bool triplesMode = false);
  /**
   * @brief Non-default destructor.
   *
   * Deletes integral files which may have been written to disk.
   */
  ~MO3CenterIntegralController();
  /**
   * @brief Getter for the sparse maps controller.
   * @return The sparse maps controller.
   */
  const std::shared_ptr<SparseMapsController> getSparseMapsController() {
    assert(_sparseMaps);
    return _sparseMaps;
  }

  /**
   * @brief Getter for the 3 center integrals by type and domain.
   * TODO: Add functionality to only load a part of the integral sets.
   *
   * @param mo3CenterType The integral type.
   * @param kDomain The domain of auxiliary functions which should at least be available (over shells of aux. functions).
   * @return The MO 3 center integrals.
   */
  const MO3CenterIntegrals& getMO3CenterInts(MO3CENTER_INTS mo3CenterType, const Eigen::SparseVector<int>& kDomain,
                                             const Eigen::SparseVector<int> paoDomain = Eigen::SparseVector<int>(0)) {
    // Check what is missing.
    Eigen::SparseVector<int> missing = getMissingDomain(mo3CenterType, kDomain);
    // Remove unused integrals.
    Eigen::SparseVector<int> unusedDomain = getUnusedDomain(kDomain);
    if (paoDomain.size() != 0) {
      unusedDomain = Eigen::VectorXi::Constant(unusedDomain.size(), 1).sparseView();
      missing = kDomain;
    }
    this->removeIntegralsByDomain(unusedDomain, mo3CenterType);
    _integrals[mo3CenterType].second = kDomain;
    // Calculate missing integrals.
    if (missing.nonZeros() != 0)
      calculateIntegrals(mo3CenterType, missing, paoDomain);
    return *_integrals[mo3CenterType].first;
  }
  /**
   * @brief Removes the integrals from memory and write them to disk
   * @param mo3CenterType The integral type to be removed.
   */
  void removeFromMemory(MO3CENTER_INTS mo3CenterType);

  /**
   * @brief Set the disk mode. This writes everything to disk if diskMode=true.
   * @param diskMode The diskMode.
   */
  void setDiskMode(bool diskMode);
  /**
   * @brief Removes all integrals from memory.
   */
  void flushAllIntegrals() {
    flushIntegrals(MO3CENTER_INTS::ab_K);
    flushIntegrals(MO3CENTER_INTS::ia_K);
    flushIntegrals(MO3CENTER_INTS::kl_K);
  }
  /**
   * @brief Removes a specific set of integrals from memory.
   * @param mo3CenterType The integral set type.
   */
  void flushIntegrals(MO3CENTER_INTS mo3CenterType);

  /**
   * @brief Getter for the projection matrix to the significant occupied/virtual orbitals
   *        for each auxiliary basis function K.
   *
   *  For each K only a subset of integrals with a selection of virtual and occupied functions
   *  is calculated. The projection matrix is the map between the small set and full set of
   *  functions for the given orbital type.
   *
   * @param orbitalType The orbital type: Occupied or virtual.
   * @return The projection matrix and the associated indices.
   */
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> getProjection(ORBITAL_TYPE orbitalType);
  /**
   * @brief Getter for the indices of the significant occupied/virtual orbitals
   *        for each auxiliary basis function K.
   *
   *  For each K only a subset of integrals with a selection of virtual and occupied functions
   *  is calculated. The indices which are calculated for each K are given for the specified orbital
   *  type by the map. If the orbital is not in the map, it was never calculated.
   *
   * @param orbitalType The orbital type: Occupied or virtual.
   * @return The projection matrix and the associated indices.
   */
  std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>> getIndices(ORBITAL_TYPE orbitalType);
  /**
   * @brief Remove the integrals domain-wise.
   * @param kDomain The domain to be removed.
   * @param mo3CenterType The integral set for which integrals should be removed.
   */
  void removeIntegralsByDomain(const Eigen::SparseVector<int>& kDomain, MO3CENTER_INTS mo3CenterType);
  /**
   * @brief Getter for the total memory requirement for the given integral type. This function is static so that
   *        it can be used without an instance of this class.
   * @param type The integral type.
   * @param sparseMaps The sparse map controller encoding which integrals need to be calculated.
   * @param auxBasisController The auxiliary basis controller.
   * @param kDomain The auxiliary domain to calculate integrals on.
   * @param triplesMode Boolean if the integrals are needed for the triples calculation.
   * @return The memory requirement.
   */
  static double getMemoryRequirement(MO3CENTER_INTS type, std::shared_ptr<SparseMapsController> sparseMaps,
                                     std::shared_ptr<BasisController> auxBasisController,
                                     const Eigen::SparseVector<int>& kDomain, bool triplesMode = false);
  /**
   * @brief Getter for the memory requirement of all three integral types. Assuming triplesMode = false.
   * @param sparseMaps The sparse map controller encoding which integrals need to be calculated.
   * @param auxBasisController The auxiliary basis controller.
   * @param kDomain The auxiliary domain to calculate integrals on.
   * @return The total memory requirement.
   */
  static double getTotalMemoryRequirement(std::shared_ptr<SparseMapsController> sparseMaps,
                                          std::shared_ptr<BasisController> auxBasisController,
                                          const Eigen::SparseVector<int>& kDomain);

 private:
  /**
   * @brief Calculate projection matrix to the subset which corresponds to the given map.
   * @param map The map.
   * @return The projection matrix.
   */
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> getProjectionMatrices(const SparseMap& map);

  /**
   * @brief Calculate the mapping indices of a given projection matrix and stores them
   *        easily accessible within a map.
   * @param k_redToFullMaps The set of projection matrices.
   * @return The indices.
   */
  std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>>
  getReducedIndices(const std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>>& k_redToFullMaps);
  /**
   * @brief Load the integrals.
   * @param mo3CenterType The integral type.
   */
  void loadIntegrals(MO3CENTER_INTS mo3CenterType, const Eigen::SparseVector<int>& kDomain);

  /**
   * @brief Calculate the integrals.
   * @param mo3CenterType The integral type.
   */
  void calculateIntegrals(MO3CENTER_INTS mo3CenterType, const Eigen::SparseVector<int>& kDomain,
                          const Eigen::SparseVector<int>& paoDomain);

  /**
   * @brief Write the integrals to disk.
   * @param mo3CenterType The integral type.
   */
  void writeToDisk(MO3CENTER_INTS mo3CenterType);

  /**
   * @brief Checks whether an integral file is already present and accessible.
   * @param mo3CenterType The integral type.
   * @return True if the file is present, false if not.
   */
  bool checkDisk(MO3CENTER_INTS mo3CenterType);

  /**
   * @brief Selects the maps used for the integral calculation.
   * @param mo3CenterType The integral type.
   * @return The prescreening maps.
   */
  std::pair<const std::shared_ptr<SparseMap>, const std::shared_ptr<SparseMap>> selectPrescreeningMaps(MO3CENTER_INTS mo3CenterType);
  ///@brief Map from significant set of PAOs to total set.
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> _projection_virt;
  ///@brief Map from significant set of occupied orbitals to total set.
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> _projection_occ;
  ///@brief Indices of PAOs within the significant set of PAOs.
  std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>> _indices_virt;
  ///@brief Indices of occupied orbitals within the significant set of orbitals.
  std::vector<std::shared_ptr<std::map<unsigned int, unsigned int>>> _indices_occ;

  /**
   * @brief The integrals and their currently available K-domain.
   */
  std::map<MO3CENTER_INTS, std::pair<std::shared_ptr<MO3CenterIntegrals>, Eigen::SparseVector<int>>> _integrals = {
      {MO3CENTER_INTS::ia_K, std::make_pair(nullptr, Eigen::SparseVector<int>(0))},
      {MO3CENTER_INTS::kl_K, std::make_pair(nullptr, Eigen::SparseVector<int>(0))},
      {MO3CENTER_INTS::ab_K, std::make_pair(nullptr, Eigen::SparseVector<int>(0))}};
  /**
   * @brief The file names.
   */
  std::map<MO3CENTER_INTS, std::string> _fileNames = {
      {MO3CENTER_INTS::ia_K, "iaK"}, {MO3CENTER_INTS::kl_K, "klK"}, {MO3CENTER_INTS::ab_K, "abK"}};

  /**
   * @brief Transform the AO 3-center integrals to MO basis.
   * @param mo3CenterType The integral type.
   * @param result The transformed integrals are written into this object.
   * @param c_K The occupied coefficients.
   * @param pao_K The virtual coefficients.
   * @param i_K The AO integrals.
   */
  inline void transformCoefficients(MO3CENTER_INTS mo3CenterType, Eigen::MatrixXd& result, const Eigen::MatrixXd& c_K,
                                    const Eigen::MatrixXd& pao_K, const Eigen::MatrixXd& i_K);
  /**
   * @brief Getter for the aux. domains over which no integrals have been calculated yet.
   * @param mo3CenterType The integral type.
   * @param kDomain The requested aux. domain.
   * @return The missing aux. domain.
   */
  Eigen::SparseVector<int> getMissingDomain(MO3CENTER_INTS mo3CenterType, const Eigen::SparseVector<int>& kDomain);

  Eigen::SparseVector<int> getUnusedDomain(const Eigen::SparseVector<int>& kDomain);
  ///@brief Auxiliary basis controller.
  std::shared_ptr<BasisController> _auxiliaryBasisController;
  ///@brief Basis controller.
  std::shared_ptr<BasisController> _basisController;
  ///@brief Sparse maps controller.
  std::shared_ptr<SparseMapsController> _sparseMaps;
  ///@brief Sparse map between aux. function and basis function which originates from
  ///       an occupied orbital.
  std::shared_ptr<SparseMap> _kToRhoMap = nullptr;
  ///@brief Sparse map between aux. function and basis function which originates from
  ///       a virtual orbital.
  std::shared_ptr<SparseMap> _kToSigmaMap = nullptr;
  ///@brief The PAO controller.
  std::shared_ptr<Eigen::MatrixXd> _virtualCoefficients;
  ///@brief The coefficients of the occupied orbitals.
  std::shared_ptr<Eigen::MatrixXd> _occupiedCoefficients;
  ///@brief The base name of the system (Needed to write files).
  std::string _fBaseName = "";
  ///@brief The system ID.
  std::string _id = "";
  ///@brief The disk mode.
  bool _diskMode = false;
  ///@brief Use triples prescreening maps indstead of pair-based maps.
  bool _triplesMode = false;
  ///@brief Print general information about the number of MO integrals to be stored and the AO integrals to be calculated.
  void printInfo(const SparseMap& kToRhoMap, const SparseMap& kToSigmaMap, const SparseMap& kToOccMap,
                 const SparseMap& kToPAOMap, MO3CENTER_INTS type, const Eigen::SparseVector<int>& kDomain);
};

} /* namespace Serenity */

#endif /* INTEGRALS_MO3CENTERINTEGRALCONTROLLER_H_ */
