/**
 * @file OrbitalPairSet.h
 *
 * @date Feb. 18, 2021
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

#ifndef SERENITY_ORBITALPAIRSET_H
#define SERENITY_ORBITALPAIRSET_H
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dense vectors.
#include <Eigen/SparseCore> //Sparse matrices.
#include <memory>           //smrt_ptr.
#include <vector>           //std::vector.

namespace H5 {
class H5File;
} // namespace H5
namespace Serenity {

/* Forward Declarations */
class OrbitalPair;
namespace HDF5 {
using H5File = H5::H5File;
} // namespace HDF5
/**
 * @class
 * @brief Handles a set of (close) orbital pairs for writing and reading from disk.
 *
 * This simply provides a convenient interface to handle a list of orbital pairs simultaneously.
 */
class OrbitalPairSet : public std::vector<std::shared_ptr<OrbitalPair>> {
 public:
  /**
   * @brief Constructor.
   */
  OrbitalPairSet();
  /**
   * @brief Default destructor.
   */
  ~OrbitalPairSet();
  /**
   * @brief Write all k-Set and (ac|bd) integrals to file.
   * @param file The file to write to.
   */
  void toHDF5(HDF5::H5File& file);
  /**
   * @brief Read all k-Set and (ac|bd) integrals from file.
   * @param file The file to read from.
   */
  void fromHDF5(HDF5::H5File& file);
  /**
   * @brief Check if integrals are available for all pairs.
   * @return True if integrals are kept in main memory, false otherwise.
   */
  bool integralsReady();
  /**
   * @brief Tell the orbital-pair set that the pair integrals have been modified externally.
   * @param inMemory New value for inMemory flag (see integralsReady())
   */
  void setInMemory(bool inMemory);
  /**
   * @brief Get the memory demand for the total orbital set.
   * @param If true, the memory for the integrals used for the sigma vector is included.
   * @return The memory demand.
   */
  double memoryDemand(bool sigmaInts);
  /**
   * @brief Calls "flushIntegrals" for each pair.
   */
  void removeInteralsFromMemory();
  /**
   * @brief Getter for the union of the pairs extended fitting domain.
   * @return The union of the extended fitting domains.
   */
  Eigen::SparseVector<int> getTotalFittingDomain();

 private:
  // Flag for the availability of integrals on disk
  bool _onDisk = false;
  // True if integrals are in main memory, false otherwise.
  bool _inMemory = false;
  // Memory requirement of the total orbital-pair set.
  std::shared_ptr<double> _memory;
};

} /* namespace Serenity */

#endif // SERENITY_ORBITALPAIRSET_H
