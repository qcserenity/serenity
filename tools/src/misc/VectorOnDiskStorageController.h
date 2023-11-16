/**
 * @file VectorOnDiskStorageController.h
 *
 * @date Apr 19, 2018
 * @author Moritz Bensberg
 *
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

#ifndef MISC_VECTORONDISKSTORAGECONTROLLER_H_
#define MISC_VECTORONDISKSTORAGECONTROLLER_H_

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <map>
#include <memory>
#include <stack>
#include <vector>

namespace H5 {
class H5File;
} // namespace H5

namespace Serenity {
namespace HDF5 {
using H5File = H5::H5File;
} // namespace HDF5

/**
 * @class VectorOnDiskStorageController VectorOnDiskStorageController.h
 * @brief A controller which controls a vector which can partially be stored
 *        on disk or memory. The vector is divided into fragments which can be
 *        identified via their label. The fragments can be altered in length of
 *        content if one desires. Additional segments can be added ad any time.
 */
class VectorOnDiskStorageController {
 public:
  /**
   * @brief Constructor.
   * @param cacheLimit The cache limit.
   * @param fileName The appendix for the file names. Must end with ".h5".
   */
  VectorOnDiskStorageController(double cacheLimit, std::string fileName);

  /**
   * @brief Default copy constructor is invalid!
   * @param copy
   */
  VectorOnDiskStorageController(const VectorOnDiskStorageController& copy) = delete;
  /**
   * @brief Copy constructor.
   * @param old
   * @param fileName
   */
  VectorOnDiskStorageController(VectorOnDiskStorageController& old, std::string fileName);
  /**
   * @brief Destructor. Removes remaining files from the disk.
   */
  virtual ~VectorOnDiskStorageController();

  /**
   * @brief Returns the part of the vector with the label <label>
   * @param label The label of the vector part.
   * @return The vector part.
   */
  std::shared_ptr<Eigen::VectorXd> getVectorSegment(std::string label);

  /**
   * @brief Stores a vector part.
   * @param vectorPart The vector part.
   * @param label The label of the vector part.
   */
  void storeVectorSegment(std::shared_ptr<Eigen::VectorXd> vectorPart, std::string label);

  /**
   * @brief Calculates the scalar product of the underlying vectors.
   * @param rhs The second vector.
   * @return Returns the scalar product.
   */
  double operator*(VectorOnDiskStorageController& rhs);

  /**
   * @brief Returns the total number of stored doubles.
   * @return
   */
  inline long int size() {
    long int totalSize = 0;
    for (const auto& it : _segmentSizes) {
      totalSize += it.second;
    }
    return totalSize;
  }

  /**
   * @brief Returns the cache limit.
   *        The cache limit is the maximum amount of memory which will be allocated for storing vector
   *        segments in memory.
   * @return
   */
  inline double getCacheLimit() {
    return _cacheLimit;
  }

  /**
   * @brief Returns a list of all segment labels.
   * @return
   */
  std::vector<std::string> getLabelList();

  /**
   * @brief Returns the appendix for the HDF5 files.
   * @return
   */
  inline std::string getHDF5FileName() {
    return _fileName;
  }

 private:
  /**
   * @brief Loads the vector segment with label <label>.
   * @param label The vector segment-label.
   * @return The vector segment.
   */
  std::shared_ptr<Eigen::VectorXd> loadVectorSegment(std::string label);

  /**
   * @brief Saves a segment to tmp/<label>_<_fileName>
   * @param vectorSegment The vector segment which is saved.
   * @param label The label of the vector segment.
   */
  void saveVectorSegment(Eigen::VectorXd& vectorSegment, std::string label);

  /**
   * @brief Removes a vector segment with the label <label> from the cache.
   *        The vector segment will remain in a file on the disk and can be loaded
   *        again.
   * @param label The label of the vector segment.
   */
  void removeVectorSegmentFromCache(std::string label);

  /**
   * @brief Deletes a vector segment from the controller.
   * @param label The label of the vector segment.
   */
  void deleteVectorSegmentFromController(std::string label);

  /// @brief The memory limit for the cache. If the limit is reached all additional data will
  /// be written to scratch.
  double _cacheLimit;
  /// @brief Memory counter.
  long int _nDoubles;
  /// @brief The appendix for the file names.
  std::string _fileName;
  /// @brief The cache.
  std::map<std::string, std::shared_ptr<Eigen::VectorXd>> _vectorPool;
  /// @brief The sizes of the vector segments.
  std::map<std::string, int> _segmentSizes;
  /// @brief The files for the vector segments.
  std::map<std::string, std::shared_ptr<H5::H5File>> _files;
};

} /* namespace Serenity */

#endif /* MISC_VECTORONDISKSTORAGECONTROLLER_H_ */
