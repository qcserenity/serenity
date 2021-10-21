/**
 * @file   CDIntegralController.h
 *
 * @date   Jun 28, 2018
 * @author Lars Hellmann
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef CDINTEGRALCONTROLLER_H_
#define CDINTEGRALCONTROLLER_H_

/* Include Serenity Internal Headers */
#include "integrals/wrappers/Libint.h"
#include "settings/Settings.h"

namespace Serenity {

/* Forward declarations */
class CDStorageController;
class Geometry;
class BasisController;
class CustomBasisController;
class CombinedShellPair;

/**
 * @class CDIntegralController CDIntegralController.h
 *
 * @brief A class that manages everything regarding the Cholesky decomposition of matrices that
 * correspond to one system (for example the two electron four center integrals).
 * It holds a storage controller for every matrix decomposed and information regarding the
 * decomposition and disk mode.
 */
class CDIntegralController : public std::enable_shared_from_this<CDIntegralController> {
 public:
  /**
   * @brief Constructor
   * @param settings The settings of the corresponding system.
   */
  CDIntegralController(const Settings& settings);

  /**
   * @brief Getter for the Cholesky decomposition threshold.
   * @return The decomposition threshold.
   */
  double getDecompositionThreshold();

  /**
   * @brief Getter for the disk mode. (If true all Cholesky vectors are stored on disk.)
   * @return The disk mode.
   */
  std::shared_ptr<bool> getDiskMode();

  /**
   * @brief Sets the disk mode to true.
   */
  void setDiskMode();

  /**
   * @brief Sets the disk mode to false.
   */
  void unsetDiskMode();

  /**
   * @brief Getter for the storage controller for a certain matrix. If the requested storage controller
   *        does not exist, it is initialized.
   * @param label The identifier of the corresponding matrix
   * @return The storage controller.
   */
  std::shared_ptr<CDStorageController> getStorageController(std::string label);

  /**
   * @brief Getter for the integrals produced by generateACDVectors()
   * @param basisController The basis controller for the default basis.
   * @param auxBasisController The basis controller for the auxiliary basis.
   */
  bool getACDVectors(std::shared_ptr<BasisController> basisController, std::shared_ptr<BasisController> auxBasisController);

  /**
   * @brief Generates a set of vectors similar to those obtained by Cholesky decomposition using an
   *        auxiliary basis (normally the acd/accd basis): L^J_{ij} = \sum_K (ij|K)(K|J)^{-1/2}
   * @param basisController The basis controller for the default basis.
   * @param auxBasisController The basis controller for the auxiliary basis.
   */
  void generateACDVectors(std::shared_ptr<BasisController> basisController, std::shared_ptr<BasisController> auxBasisController);

  /**
   * @brief Generates the atomic Cholesky decomposition basis and writes it to file
   * @param geom The geometry providing atom types and basis sets to generate the acd basis from
   * @param op_label The label for the operator used.
   * @param op The operator used to calculate the integrals used to generate the basis.
   */
  void generateACDBasis(std::shared_ptr<Geometry> geom, std::string op_label = "", LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb);

  /**
   * @brief Generates the atomic Cholesky decomposition basis and writes it to file
   * @param geom The geometry providing atom types and basis sets to generate the acd basis from
   * @param label The label of the basis the acd basis should be generate from.
   * @param op_label The label for the operator used.
   * @param op The operator used to calculate the integrals used to generate the basis.
   *
   */
  void generateACDBasis(std::shared_ptr<Geometry> geom, std::string label, std::string op_label = "",
                        LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb);

  /**
   * @brief Generates the atomic-compact Cholesky decomposition basis and writes it to file
   * @param geom The geometry providing atom types and basis sets to generate the accd basis from
   * @param op_label The label for the operator used.
   * @param op The operator used to calculate the integrals used to generate the basis.
   */
  void generateACCDBasis(std::shared_ptr<Geometry> geom, std::string op_label = "",
                         LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb);

  /**
   * @brief Generates the atomic-compact Cholesky decomposition basis and writes it to file
   * @param geom The geometry providing atom types and basis sets to generate the accd basis from.
   * @param label The label of the basis the accd basis should be generate from.
   * @param op_label The label for the operator used.
   * @param op The operator used to calculate the integrals used to generate the basis.
   */
  void generateACCDBasis(std::shared_ptr<Geometry> geom, std::string label, std::string op_label = "",
                         LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb);

  /**
   * @brief Gets all unique atom types and a corresponding index from the given geometry.
   * @param geom The geometry providing atom types and basis sets to generate the accd basis from.
   * @param label The label of the basis the atomtypes should be collected for.
   * @return A vector of atomtypes and an index where that type is found in the provided geometry.
   */
  static std::vector<std::pair<std::string, unsigned int>> getAtomTypes(std::shared_ptr<Geometry> geom, std::string label);

  /**
   * @brief Checks if the file provided in path contains basis functions for all atomtypes provided at the correct
   *        Cholesky decomposition threshold.
   * @param path The path of the basis file.
   * @param atomTypes A list of all needed atomtypes.
   * @return True if basis functions for all atomtypes are found at the correct threshold.
   */
  bool checkBasisFile(std::string path, std::vector<std::pair<std::string, unsigned int>> atomTypes);

  /**
   * @brief Reads the shell split from the basis file.
   * @param path The path to the basis file.
   * @param element The element the shell split should be read.
   * @return The shell split. An equal index means the basis functions belong to the same original shell.
   *                   This includes the added basis functions with different angular momentum.
   */
  std::vector<unsigned int> getShellSplit(std::string path, std::string element);

  /**
   * @brief Generates pseudo orbital coefficients from a matrix
   * @param p The (density) matrix the coefficients should match
   * @return It returns the std::pair of a matrix (containing the pseudo coefficients) and a vector (containing the
   * signs that might arise from negative eigenvalus of p).
   */
  std::pair<Eigen::VectorXd, Eigen::MatrixXd> generatePseudoCoefficients(Eigen::MatrixXd p);

  /**
   * @brief Generates additional shells to obtain a complete shell product in spherical basis.
   * @param tmpShellPair The incomplete shell pair product.
   * @return A vector of shells that represent the complete product of the original shell pair.
   */
  std::vector<std::shared_ptr<CombinedShellPair>>
  generateCompleteCombinedSphericalShell(std::shared_ptr<CombinedShellPair> tmpShellPair);

  /**
   * @brief Splits a shell into different subshells if the number of primitive functions is to large to be handled by
   * libint as one.
   * @param tmpShellPair The shell pair containing the primitive functions to be split.
   * @return A vector of the split subshells.
   */
  std::vector<std::shared_ptr<const CombinedShellPair>> splitPrimitives(std::shared_ptr<CombinedShellPair> tmpShellPair);

  /**
   * @brief Deletes all data held by this controller, both in memory and on disk. (Including all CDStorageController and
   * the Cholesky vectors they hold)
   */
  void cleanup();

 private:
  // The settings of the system
  const Settings& _settings;
  // The threshold of the Cholesky decomposition
  double _decompositionThreshold;
  // A map of all storage controllers for this system
  std::map<std::string, std::shared_ptr<CDStorageController>> _cdStorageController;
  // Switch for the disk Mode. A shared pointer of this is given to each storage controller
  std::shared_ptr<bool> _diskMode;
  // Initializes a new storage controller
  void addMatrix(std::string label);
  // returns the character corresponding to an angular momentum
  char getAngularMomentumChar(unsigned int angularMom);
};

} /* namespace Serenity */

#endif /* CDINTEGRALCONTROLLER_H_ */
