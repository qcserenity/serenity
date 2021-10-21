/**
 * @file   CDStorageController.h
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

#ifndef CDSTORAGECONTROLLER_H_
#define CDSTORAGECONTROLLER_H_

/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "io/HDF5.h"
/* Include Std and External Headers */
#include <Eigen/Eigen>
#include <limits>
#include <memory>

namespace Serenity {

/* Forward declarations */
class BasisController;
class CDIntegralController;

/**
 * @class CDStorageController CDStorageController.h
 *
 * @brief A class that manages (Cholesky-)Vectors that belong to one matrix.
 * In this stage it is possible to store vectors and load them again. However,
 * editing already stored vectors is not supported.
 *
 * This class primarily uses row major vectors, since HDF5 stores all
 * information in a row major fashion. Therefore, column major objects are
 * converted to row major, creating additional overhead.
 */
class CDStorageController : public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @brief Constructor. Calls addMatrix() and creates files on disk for data storage.
   *
   * @param systemName The name of the system the vectors (and matrix) belong to
   * @param label The name of the matrix
   * @param cdIntController The Cholesky integral controller of the system
   */
  CDStorageController(std::string systemName, std::string label, std::shared_ptr<CDIntegralController> cdIntController);

  /**
   * @brief Destructor. Removes all files created by this controller from disk in addition to the default destructor.
   */
  virtual ~CDStorageController();

  /**
   * @brief Getter for a vector held in memory. Use loadBatch() to load vectors from disk to memory.
   *
   * @param index The (column) index of the vector
   * @return The requested vector
   */
  std::shared_ptr<Eigen::RowVectorXd> loadVector(const int index);

  /**
   * @brief Stores a vector in memory or on disk if needed
   *
   * @param index The (column) index of the vector
   * @param choleskyVector The vector to be stored
   */
  void storeVector(const int index, std::shared_ptr<Eigen::RowVectorXd> choleskyVector);

  /**
   * @brief Loads the diagonal of the matrix from disk
   *
   * @return The diagonal of the matrix
   */
  std::shared_ptr<Eigen::VectorXd> loadDiag();

  /**
   * @brief Stores the diagonal of the matrix on disk
   *
   * @param vector The diagonal as a vector to be stored
   */
  void storeDiag(std::shared_ptr<Eigen::VectorXd> vector);

  /**
   * @brief Getter for the number of Cholesky vectors stored by this controller
   *
   * @return The number of already stored vectors
   */
  unsigned int getNVectors();

  /**
   * @brief Dumps vectors from memory to disk and frees the used memory
   */
  void freeMemory();

  /**
   * @brief Setter for the disk mode.
   */
  void setDiskMode();

  /**
   * @brief A function that makes sure that all vector in the requested batch are held in memory
   *
   * @param startIndex Start index for the requested batch
   * @param batchSize Size of the requested batch
   * @return The actual size of the batch loaded to memory starting from the given startIndex. This can be smaller than
   * the requested batch size!
   */
  unsigned int loadBatch(unsigned int startIndex, unsigned int batchSize = 99999999);

  /**
   * @brief deletes the last batch that was loaded from memory. _lastBatchStart and _lastBatchSize are used to identify
   * that batch.
   */
  void freeLastBatch();

  /**
   * @brief Sets the dimension of the vectors stored. This has to be called if the diagonal is not stored first.
   *
   * @param dim The dimension of the vectors that are stored.
   */
  void setVectorDimension(unsigned int dim);

  /**
   * @brief If the dimension of the Cholesky vectors is known previous to their calculation this function is used to
   * allocate all vectors. This also enables index blockwise storage using storeSegment()
   * @param nVec Number of Cholesky vectors
   * @param dim Length of the Cholesky vectors
   */
  void allocateVectors(unsigned int nVec, unsigned int dim);

  /**
   * @brief Stores Segments of Cholesky vectors on disk.
   * @param index The index of the Cholesky vector to add the segment to.
   * @param firstIndex First index that will be written to.
   * @param val Vector holding the values added to the segment.
   * @param stride Stride between the entries.
   */
  void storeSegment(unsigned int index, unsigned int firstIndex, std::shared_ptr<std::vector<double>> val,
                    unsigned int stride = 1);

  /**
   * @brief Flushes the data to disk if HDF5 kept it in memory to this point.
   */
  void flushFile();
  /**
   * @brief Notify function
   */
  void notify();
  /**
   * @brief Adds this Storage Controller to the list of objects notified when the BasisController changes.
   * @param basCont The BasisController to which the CDStorageController has to be sensitive.
   */
  void addSensitiveBasis(std::shared_ptr<BasisController> basCont);

  /**
   * @brief Sets the flag that the stored vectors are up to date.
   */
  void setUpToDate();
  /**
   * @brief Gets the flag if the stored vectors are up to date.
   * @return Returns true if the stored vectors are up to date.
   */
  bool getUpToDate();

 private:
  // return the number of all vectors kept in memory
  unsigned int loadBatchMem(unsigned int startIndex, unsigned int batchSize);

  // load vectors in batch from disk
  unsigned int loadBatchDisk(unsigned int startIndex, unsigned int batchSize);

  // updates the array of existing and pre-opened datasets
  void updateDataSets();

  // initializes the matrix and prepares the corresponding HDF5 file
  void addMatrix();

  // The name of the system
  std::string _systemName;
  // Switch for storing vectors on disk or in memory
  std::shared_ptr<bool> _diskMode;
  // The Cholesky controller for the system
  std::weak_ptr<CDIntegralController> _cdIntegralController;
  // The pool of vectors that are stored in memory
  std::shared_ptr<std::vector<std::shared_ptr<Eigen::RowVectorXd>>> _vectorPool;
  // The label of the matrix
  std::string _label;
  // Definition of the doubletype for easier HDF5 handling
  hid_t _doubletype;
  // Dimension of the vectors
  unsigned int _vectorDim;
  // The HDF5 file the vectors are stored in
  std::shared_ptr<H5::H5File> _file;
  // List of pre-opened HDF5 datasets
  std::vector<H5::DataSet> _dataSets;
  // start index of the last batch loaded to memory
  unsigned int _lastBatchStart;
  // size of the last batch loaded to memory
  unsigned int _lastBatchSize;
  // Flag if diskspace for all vectors is allocated
  bool _vectorsAllocated;
  // Flag if the stored vectors are up to date
  bool _vectorsUpToDate;
};

} /* namespace Serenity */

#endif /* CDSTORAGECONTROLLER_H_ */
